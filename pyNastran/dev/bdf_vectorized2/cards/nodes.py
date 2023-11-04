from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.femutils.utils import duplicates
from pyNastran.bdf.bdf import GRID

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank,
    components_or_blank)
from pyNastran.bdf.field_writer_8 import print_float_8, set_string8_blank_if_default
from pyNastran.bdf.field_writer_16 import print_float_16, set_string16_blank_if_default
from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.nptyping_interface import NDArrayN3int, NDArrayN3float, NDArrayNint


class Nodes:
    def __init__(self, model):
        self.model = model
        self.grid = model.grid
        self.spoints = model.spoints
        self.epoints = model.epoints
        self.xyz_cid0 = None

    #@property
    #def nids(self):
        #self.grid.make_current()
        #return self.grid.nid

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.grid):
            self.grid.write_card(size, is_double, bdf_file)
        #if len(self.spoints):
            #self.spoints.write_card(size, is_double, bdf_file)
        #if len(self.epoints):
            #self.epoints.write_card(size, is_double, bdf_file)

    def __len__(self):
        """returns the number of nodes"""
        return len(self.model.grid) + len(self.spoints) + len(self.epoints)

    def __repr__(self):
        return self.repr_indent('')

    def repr_indent(self, indent=''):
        msg = '%s<Nodes>:\n' % indent
        msg += '%s  GRID: %s\n' % len(indent, self.grid)
        msg += '%s  SPOINT: %s\n' % len(indent, self.spoints)
        msg += '%s  EPOINT: %s\n' % len(indent, self.epoints)

    def get_by_nid(self, nid):
        #self.grid.make_current()
        inid = np.searchsorted(self.nid, nid)
        return self[inid]

    def get_displacement_index_xyz_cp_cd(self,
                                         fdtype: str='float64',
                                         idtype: str='int32') -> tuple[dict[int, NDArrayNint],
                                                                       dict[int, NDArrayNint],
                                                                       NDArrayN3float, NDArrayN3int]:
        """
        Get index and transformation matricies for nodes with
        their output in coordinate systems other than the global.
        Used in combination with ``OP2.transform_displacements_to_global``

        Parameters
        ----------
        fdtype : str
            the type of xyz_cp
        idtype : str
            the type of nid_cp_cd

        Returns
        -------
        icd_transform : dict{int cd : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their output (`CD`) in that
            coordinate system.
        icp_transform : dict{int cp : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their input (`CP`) in that
            coordinate system.
        xyz_cp : (n, 3) float ndarray
            points in the CP coordinate system
        nid_cp_cd : (n, 3) int ndarray
            node id, CP, CD for each node

        Examples
        --------
        # assume GRID 1 has a CD=10, CP=0
        # assume GRID 2 has a CD=10, CP=0
        # assume GRID 5 has a CD=50, CP=0
        >>> model.point_ids
        [1, 2, 5]
        >>> out = model.get_displacement_index_xyz_cp_cd()
        >>> icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
        >>> nid_cp_cd
        [
           [1, 0, 10],
           [2, 0, 10],
           [5, 0, 50],
        ]
        >>> icd_transform[10]
        [0, 1]

        >>> icd_transform[50]
        [2]
        """
        self.grid.make_current()
        nids_cd_transform = defaultdict(list)  # type: dict[int, np.ndarray]
        nids_cp_transform = defaultdict(list)  # type: dict[int, np.ndarray]

        nnodes = len(self.model.grid)
        nspoints = 0
        nepoints = 0
        spoints = None
        epoints = None
        nrings = len(self.model.ringaxs)
        if self.model.spoints:
            spoints = list(self.model.spoints)
            nspoints = len(spoints)
        if self.model.epoints:
            epoints = list(self.model.epoints)
            nepoints = len(epoints)

        if nnodes + nspoints + nepoints + nrings == 0:
            msg = 'nnodes=%s nspoints=%s nepoints=%s nrings=%s' % (nnodes, nspoints, nepoints, nrings)
            raise ValueError(msg)

        xyz_cp = np.zeros((nnodes + nspoints + nepoints, 3), dtype=fdtype)
        nid_cp_cd = np.zeros((nnodes + nspoints + nepoints, 3), dtype=idtype)

        xyz_cp[:nnodes, :] = self.model.grid.xyz
        nid_cp_cd[:nnodes, 0] = self.model.grid.nid
        nid_cp_cd[:nnodes, 1] = self.model.grid.cp
        nid_cp_cd[:nnodes, 2] = self.model.grid.cd

        for nid, cp, cd in nid_cp_cd:
            nids_cp_transform[cp].append(nid)
            nids_cd_transform[cd].append(nid)

        i = nnodes

        if nspoints:
            for nid in sorted(spoints):
                nid_cp_cd[i, 0] = nid
                i += 1
        if nepoints:
            for nid in sorted(epoints):
                nid_cp_cd[i, 0] = nid
                i += 1
        assert nid_cp_cd[:, 0].min() > 0, nid_cp_cd[:, 0].tolist()
        #assert nid_cp_cd[:, 0].min() > 0, nid_cp_cd[:, 0].min()

        # sorting
        nids = nid_cp_cd[:, 0]
        isort = nids.argsort()
        nid_cp_cd = nid_cp_cd[isort, :]
        xyz_cp = xyz_cp[isort, :]

        icp_transform = {}
        icd_transform = {}
        nids_all = nid_cp_cd[:, 0]

        # get the indicies of the xyz array where the nodes that
        # need to be transformed are
        for cd, nids in sorted(nids_cd_transform.items()):
            if cd in [0, -1]:
                continue
            nids = np.array(nids)
            icd_transform[cd] = np.where(np.in1d(nids_all, nids))[0]

        for cp, nids in sorted(nids_cp_transform.items()):
            if cp in [-1]:
                continue
            nids = np.array(nids)
            icp_transform[cp] = np.where(np.in1d(nids_all, nids))[0]
        return icd_transform, icp_transform, xyz_cp, nid_cp_cd

    def get_node_index(self, nids, allow0=False):
        """maps the requested nodes to their desired index in the array"""
        self.grid.make_current()
        nids = np.asarray(nids)
        nids_ravel = nids.ravel()

        sorted_nodes = self.nids
        #print('sorted_nodes = %s' % sorted_nodes.tolist())
        #print(nids)
        #print(nids_ravel)
        i = np.searchsorted(sorted_nodes, nids_ravel)

        if allow0:
            izero = np.where(nids_ravel == 0)[0]
            i[izero] = -1
        else:
            if not np.array_equal(sorted_nodes[i], nids_ravel):
                msg = (
                    '  nids:\n%s\n'
                    '  self.nid = %s\n'
                    '  i = %s\n'
                    '  nid[i]    = %s\n'
                    '  nid_ravel = %s\n'
                    '  nids_new:\n%s\n' % (
                        nids, sorted_nodes.tolist(),
                        i,
                        sorted_nodes[i].tolist(),
                        nids_ravel.tolist(),
                        nids.reshape(nids.shape)))
                raise RuntimeError(msg)
        i = i.reshape(nids.shape)
        return i


class GRIDv:
    """
    +------+-----+----+----+----+----+----+----+------+
    |   1  |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
    +======+=====+====+====+====+====+====+====+======+
    | GRID | NID | CP | X1 | X2 | X3 | CD | PS | SEID |
    +------+-----+----+----+----+----+----+----+------+
    """
    card_name = 'GRID'
    def __init__(self, model):
        """initializes the vectorized GRID class"""
        self.model = model
        self.is_current = True
        self.nid = np.array([], dtype='int32')
        self.xyz = np.array([], dtype='float64')
        self.cp = np.array([], dtype='int32')
        self.cd = np.array([], dtype='int32')
        self.ps = np.array([], dtype='|U8')
        self.seid = np.array([], dtype='int32')

        self._nid = []
        self._xyz = []
        self._cp = []
        self._cd = []
        self._ps = []
        self._seid = []
        self.comment = defaultdict(str)

    def add(self, nid: int, xyz: list[float], cp: int=0, cd: int=0, ps: str='', seid: int=0, comment: str=''):
        """
        Creates the GRID card

        Parameters
        ----------
        nid : int
            node id
        cp : int; default=0
            the xyz coordinate frame
        xyz : (3, ) float ndarray; default=None -> [0., 0., 0.]
            the xyz/r-theta-z/rho-theta-phi values
        cd : int; default=0
            the analysis coordinate frame
        ps : str; default=''
            Additional SPCs in the analysis coordinate frame (e.g. '123').
            This corresponds to DOF set ``SG``.
        seid : int; default=0
            superelement id
            TODO: how is this used by Nastran???
        comment : str; default=''
            a comment for the card
        """
        self.is_current = False
        self._nid.append(nid)
        self._xyz.append(xyz)
        self._cp.append(cp)
        self._cd.append(cd)
        self._ps.append(ps)
        self._seid.append(seid)
        if comment:
            self.comment[nid] = _format_comment(comment)

    def add_card(self, card, comment: str=''):
        """
        Adds a GRID card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        nfields = len(card)
        #print('card = %s' % card)
        #: Node ID
        nid = integer(card, 1, 'nid')

        #: Grid point coordinate system
        cp = integer_or_blank(card, 2, 'cp', 0)

        #: node location in local frame
        xyz = [
            double_or_blank(card, 3, 'x1', 0.),
            double_or_blank(card, 4, 'x2', 0.),
            double_or_blank(card, 5, 'x3', 0.)]

        if nfields > 6:
            #: Analysis coordinate system
            cd = integer_or_blank(card, 6, 'cd', 0)

            #: SPC constraint
            ps = components_or_blank(card, 7, 'ps', '')
            #u(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            seid = integer_or_blank(card, 8, 'seid', 0)
            assert len(card) <= 9, 'len(GRID card) = %i\ncard=%s' % (len(card), card)
        else:
            cd = 0
            ps = ''
            seid = 0
        self.add(nid, xyz, cp, cd, ps, seid, comment=comment)

    def get_by_nid(self, nid):
        self.make_current()
        inid = np.searchsorted(self.nid, nid)
        return self[inid]

    def check_if_current(self, nid, nids):
        """we split this up to reason about it easier"""
        if self.is_current:
            if nid in nids:
                # card exists, so we use that slot
                add_card = False
            else:
                add_card = True
        else:
            add_card = True
        return add_card

    def update(self, grid):
        """functions like a dictionary"""
        nid = grid.nid

        add_grid = self.check_if_current(nid, self.nid)
        if add_grid:
            self.add(nid, grid.xyz, cp=grid.cp, cd=grid.cd,
                     ps=grid.ps, seid=grid.seid, comment=grid.comment)
            self.is_current = False
        else:
            inid = np.where(nid == self.nid)[0]
            self.nid[inid] = grid.nid
            self.xyz[inid] = grid.xyz
            self.cp[inid] = grid.cp
            self.cd[inid] = grid.cd
            self.ps[inid] = grid.ps
            self.seid[inid] = grid.seid
            #self.comment[nid] = comment
            #self.is_current = True  # implicit

    def make_current(self):
        """creates an array of the GRID points"""
        if not self.is_current:
            nnid = len(self.nid)
            if nnid > 0: # there are already nodes in self.nid
                self.nid = np.hstack([self.nid, self._nid])
                self.xyz = np.vstack([self.xyz, self._xyz])
                self.cp = np.hstack([self.cp, self._cp])
                self.cd = np.hstack([self.cd, self._cd])
                self.ps = np.hstack([self.ps, self._ps])
                self.seid = np.hstack([self.seid, self._seid])
                # don't need to handle comments
            else:
                self.nid = np.array(self._nid, dtype='int32')
                self.xyz = np.array(self._xyz, dtype='float64')
                self.cp = np.array(self._cp, dtype='int32')
                self.cd = np.array(self._cd, dtype='int32')
                self.ps = np.array(self._ps, dtype='|U8')
                self.seid = np.array(self._seid, dtype='int32')

            unid = np.unique(self.nid)
            if len(self.nid) != len(unid):
                duplicate_nodes = duplicates(self.nid)
                msg = ('there are duplicate nodes\n'
                       'nid =%s; n=%s\n'
                       'unid=%s; n=%s\n'
                       'duplicates=%s' % (
                    self.nid, len(self.nid),
                    unid, len(unid),
                    duplicate_nodes))
                raise RuntimeError(msg)

            isort = np.argsort(self.nid)
            self.nid = self.nid[isort]
            self.xyz = self.xyz[isort, :]
            self.cp = self.cp[isort]
            self.cd = self.cd[isort]
            self.ps = self.ps[isort]
            self.seid = self.seid[isort]

            self._nid = []
            self._xyz = []
            self._cp = []
            self._cd = []
            self._ps = []
            self._seid = []
            self.is_current = True

    def cross_reference(self, model: BDF) -> None:
        """does this do anything?"""
        self.make_current()

    def write_card(self, size: int=8, is_double: bool=False, bdf_file=None):
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int; default=8
            the size of the card (8/16)
        is_double : bool; default=False
            should this card be written with double precision

        Returns
        -------
        msg : str
            the card as a string
        """
        assert bdf_file is not None
        if size == 8:
            return self.write_card_8(bdf_file=bdf_file)
        else:
            return self.write_card_16(is_double, bdf_file=bdf_file)

    def write_card_8(self, bdf_file=None) -> str:
        """
        Writes a GRID card in 8-field format
        """
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for nid, xyz, cp, cd, ps, seid in zip(
            self.nid, self.xyz, self.cp, self.cd, self.ps, self.seid):
            cps = set_string8_blank_if_default(cp, 0)
            if [cd, ps, seid] == [0, '', 0]:
                # default
                print_float_8(xyz[0])
                print_float_8(xyz[1])
                print_float_8(xyz[2])
                msgi = 'GRID    %8i%8s%s%s%s\n' % (
                    nid, cps,
                    print_float_8(xyz[0]),
                    print_float_8(xyz[1]),
                    print_float_8(xyz[2]),
                )
            else:
                cds = set_string8_blank_if_default(cd, 0)
                seids = set_string8_blank_if_default(seid, 0)
                msgi = 'GRID    %8i%8s%s%s%s%s%8s%s\n' % (
                    nid, cps,
                    print_float_8(xyz[0]),
                    print_float_8(xyz[1]),
                    print_float_8(xyz[2]),
                    cds, ps, seids)
            msg += self.comment[nid] + msgi
        bdf_file.write(msg)
        return msg

    def write_card_16(self, is_double: bool=False, bdf_file=None) -> str:
        """
        Writes a GRID card in 16-field format
        """
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for nid, xyz, cp, cd, ps, seid in zip(
            self.nid, self.xyz, self.cp, self.cd, self.ps, self.seid):
            cps = set_string8_blank_if_default(cp, 0)
            cp = set_string16_blank_if_default(cp, 0)
            cd = set_string16_blank_if_default(cd, 0)
            seid = set_string16_blank_if_default(seid, 0)

            if is_double:
                if [cd, ps, seid] == [0, '', 0]:
                    msgi = ('GRID*   %16i%16s%16s%16s\n'
                            '*       %16s\n' % (
                                nid,
                                cp,
                                print_scientific_double(xyz[0]),
                                print_scientific_double(xyz[1]),
                                print_scientific_double(xyz[2])))
                else:
                    msgi = ('GRID*   %16i%16s%16s%16s\n'
                            '*       %16s%16s%16s%16s\n' % (
                                nid,
                                cp,
                                print_scientific_double(xyz[0]),
                                print_scientific_double(xyz[1]),
                                print_scientific_double(xyz[2]),
                                cd, ps, seid))
            else:
                if [cd, self.ps, self.seid] == [0, '', 0]:
                    msgi = ('GRID*   %16i%16s%16s%16s\n'
                            '*       %16s\n' % (
                                nid,
                                cp,
                                print_float_16(xyz[0]),
                                print_float_16(xyz[1]),
                                print_float_16(xyz[2])))
                else:
                    msgi = ('GRID*   %16i%16s%16s%16s\n'
                            '*       %16s%16s%16s%16s\n' % (
                                nid,
                                cp,
                                print_float_16(xyz[0]),
                                print_float_16(xyz[1]),
                                print_float_16(xyz[2]),
                                cd, ps, seid))
                msg += self.comment[nid] + msgi
        bdf_file.write(msg)
        return msg

    def __len__(self):
        self.make_current()
        return len(self.nid)

    def __repr__(self):
        self.make_current()
        msg = 'GRIDv; ngrids=%s:\n' % len(self.nid)
        msg += '  nid = %s\n' % self.nid

        ucp = np.unique(self.cp)
        if len(ucp) == 1 and ucp[0] == 0:
            msg += '  ucp = %s\n' % ucp
        else:
            msg += '  cp = %s\n' % self.cp

        ucd = np.unique(self.cd)
        if len(ucd) == 1 and ucd[0] == 0:
            msg += '  ucd = %s\n' % ucd
        else:
            msg += '  cd = %s\n' % self.cd

        ups = np.unique(self.ps)
        if len(ups) == 1 and ups[0] == '':
            msg += '  ups = %s\n' % ups
        else:
            msg += '  ps = %s\n' % self.ps

        useid = np.unique(self.seid)
        if len(useid) == 1 and useid[0] == 0:
            msg += '  useid = %s\n' % useid
        else:
            msg += '  seid = %s\n' % self.seid
        #msg += '  xyz =\n%s' % self.xyz
        return msg

    #def __iter__(self):
        #pass
    #def __next__(self):
        #pass
    #def __items__(self):
        #pass
    #def __keys__(self):
        #pass
    #def __values__(self):
        #pass

    def __getitem__(self, i):
        """this works on index"""
        self.make_current()
        nid = self.nid[i]
        return GRID(nid, self.xyz[i], cp=self.cp[i], cd=self.cd[i],
                    ps=self.ps[i], seid=self.seid[i], comment=self.comment[nid])

    def __setitem__(self, i, value):
        pass
    #def __delitem__(self, i):
        #pass
