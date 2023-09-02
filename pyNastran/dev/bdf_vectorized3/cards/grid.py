from __future__ import annotations
from itertools import zip_longest
from collections import defaultdict
from typing import Callable, Any, TYPE_CHECKING
import numpy as np
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8, print_float_8 # , print_field_8
from pyNastran.bdf.field_writer_16 import print_scientific_16 # , print_field_16 # print_card_16,
from pyNastran.bdf.field_writer_double import print_scientific_double
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.nodes import compress_xpoints
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, components_or_blank,
    integer_or_string, blank)

from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int
from pyNastran.dev.bdf_vectorized3.cards.base_card import VectorizedBaseCard

from pyNastran.femutils.coord_transforms import (
    xyz_to_rtz_array,
    xyz_to_rtp_array, # xyz to xxx transforms
    rtz_to_xyz_array as cylindrical_to_rectangular,
    rtp_to_xyz_array as spherical_to_rectangular, # xxx to xyz transforms
    #rtz_to_rtp_array, rtp_to_rtz_array, # rtp/rtz and rtz/rtp transforms
)

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class XPOINT(VectorizedBaseCard):
    """
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |   1    |  2  |  3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+======+=====+=====+=====+=====+=====+=====+
    | SPOINT | ID1 | THRU | ID2 |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    | SPOINT | ID1 | ID1  | ID3 | ID4 | ID5 | ID6 | ID7 | ID8 |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |        | ID8 | etc. |     |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.ids = np.array([], dtype='int32')

    def add(self, ids: np.ndarray, comment: str=''):
        if isinstance(ids, integer_types):
            ids = [ids]
        self.cards.append((ids, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        ids = []
        for i in range(1, len(card)):
            field = integer_or_string(card, i, 'ID%d' % i)
            ids.append(field)
        self.cards.append((ids, comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        ids = []
        for i, (idsi, comment) in enumerate(self.cards):
            idsi2 = expand_thru(idsi)
            ids.extend(idsi2)
            self.comment[i] = comment

        try:
            self.ids = np.array(ids, dtype='int32')
        except OverflowError:
            self.ids = np.array(ids, dtype='int64')
        self._sort()
        self.cards = []

    def _sort(self) -> None:
        iarg = np.argsort(self.ids)
        uarg = np.unique(iarg)
        #nvalues = len(self.node_id)
        if not np.array_equal(uarg, iarg):
            self.ids = self.ids[iarg]

    def write(self, size: int=8, is_double: bool=False) -> str:
        lines = []
        if size == 8:
            print_card = print_card_8
        #node_id = array_str(self.node_id, size=8)
        lists_fields = compress_xpoints(self.type, self.ids)
        for fields in lists_fields:
            msg = print_card(fields)
            lines.append(msg)
        return ''.join(lines)

class SPOINT(XPOINT):
    _id_name = 'spoint_id'
    @property
    def spoint_id(self):
        return self.ids
    @spoint_id.setter
    def spoint_id(self, spoint_id: np.ndarray):
        self.ids = spoint_id

    #def slice_card(self, i: np.ndarray) -> SPOINT:
        #"""uses a node_index to extract SPOINTs"""
        #assert len(self.ids) > 0, self.ids
        #i = np.atleast_1d(np.asarray(i, dtype=self.ids.dtype))
        #i.sort()
        #spoint = SPOINT(self.model)
        #spoint.n = len(i)
        #spoint.ids = self.ids[i]
        #return spoint

class GRDSET(BaseCard):
    """
    Defines default options for fields 3, 7, 8, and 9 of all GRID entries.

    +--------+-----+----+----+----+----+----+----+------+
    |    1   |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
    +========+=====+====+====+====+====+====+====+======+
    | GRDSET |     | CP |    |    |    | CD | PS | SEID |
    +--------+-----+----+----+----+----+----+----+------+

    """
    type = 'GRDSET'

    #: allows the get_field method and update_field methods to be used
    _field_map = {2:'cp', 6:'cd', 7:'ps', 8:'seid'}

    @classmethod
    def _init_from_empty(cls):
        cp = 0
        cd = 1
        ps = 0
        seid = 0
        return GRDSET(cp, cd, ps, seid, comment='')

    def __init__(self, cp: int=0, cd: int=0, ps: int=0, seid: int=0, comment: str=''):
        """
        Creates the GRDSET card

        Parameters
        ----------
        cp : int; default=0
            the xyz coordinate frame
        cd : int; default=0
            the analysis coordinate frame
        ps : int; default=0
            Additional SPCs in the analysis coordinate frame (e.g. 123).
            This corresponds to DOF set ``SG``.
        seid : int; default=0
            superelement id
            TODO: how is this used by Nastran???
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: Output Coordinate System
        self.cp = cp

        #: Analysis coordinate system
        self.cd = cd

        #: Default SPC constraint on undefined nodes
        self.ps = ps

        #: Superelement ID
        self.seid = seid

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a GRDSET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Grid point coordinate system
        blank(card, 1, 'blank')

        cp = integer_or_blank(card, 2, 'cp', default=0)
        blank(card, 3, 'blank')
        blank(card, 4, 'blank')
        blank(card, 5, 'blank')
        cd = integer_or_blank(card, 6, 'cd', default=0)

        ps = integer_or_blank(card, 7, 'ps', default=0)
        seid = integer_or_blank(card, 8, 'seid', default=0)
        assert len(card) <= 9, f'len(GRDSET card) = {len(card):d}\ncard={card}'
        return GRDSET(cp, cd, ps, seid, comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref: bool
            has this model been cross referenced

        """
        assert isinstance(self.cp, integer_types), 'cp=%r' % self.cp
        assert isinstance(self.cd, integer_types), 'cd=%r' % self.cd
        assert isinstance(self.ps, integer_types), 'ps=%r' % self.ps
        assert isinstance(self.seid, integer_types), 'seid=%r' % self.seid

    def write(self, size=8):
        cp = set_blank_if_default(self.cp, 0)
        cd = set_blank_if_default(self.cd, 0)
        ps = set_blank_if_default(self.ps, 0)
        seid = set_blank_if_default(self.SEid(), 0)
        list_fields = ['GRDSET', None, cp, None, None, None, cd, ps, seid]
        return print_card_8(list_fields)


class GRID(VectorizedBaseCard):
    _id_name = 'node_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self._is_sorted = False
        self.node_id = np.array([], dtype='int32')
        self.cp = np.array([], dtype='int32')
        self.cd = np.array([], dtype='int32')
        self.ps = np.array([], dtype='int32')
        self.seid = np.array([], dtype='int32')
        self.xyz = np.zeros((0, 3), dtype='float64')
        self._xyz_cid0 = np.zeros((0, 3), dtype='float64')

    def add(self, nid: int, xyz: np.ndarray,
            cp: int=0, cd: int=0,
            ps: int=0, seid: int=0, comment: str=''):
        ps = 0 if ps == '' else ps
        assert isinstance(ps, int), ps
        self.cards.append((nid, xyz, cp, cd, ps, seid, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')
        nfields = len(card)
        #: Node ID
        nid = integer(card, 1, 'nid')

        GRID_CP_DEFAULT = GRID_CD_DEFAULT = GRID_PS_DEFAULT = GRID_SEID_DEFAULT = -100

        #: Grid point coordinate system
        cp = integer_or_blank(card, 2, 'cp', default=GRID_CP_DEFAULT)

        #: node location in local frame
        xyz = [
            double_or_blank(card, 3, 'x1', default=0.),
            double_or_blank(card, 4, 'x2', default=0.),
            double_or_blank(card, 5, 'x3', default=0.),
        ]

        if nfields > 6:
            #: Analysis coordinate system
            cd = integer_or_blank(card, 6, 'cd', default=GRID_CD_DEFAULT)

            #: SPC constraint
            ps = int(components_or_blank(card, 7, 'ps', default=GRID_PS_DEFAULT))
            #u(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            seid = integer_or_blank(card, 8, 'seid', default=GRID_SEID_DEFAULT)
            assert len(card) <= 9, f'len(GRID card) = {len(card):d}\ncard={card}'
        else:
            cd = GRID_CD_DEFAULT
            ps = GRID_PS_DEFAULT
            seid = GRID_SEID_DEFAULT
            assert len(card) >= 2, f'len(GRID card) = {len(card):d}\ncard={card}'
        self.cards.append((nid, xyz, cp, cd, ps, seid, comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        if self.debug:
            self.model.log.debug('parse GRID')

        try:
            node_id, cp, cd, xyz, ps, seid = self._setup(
                ncards, self.cards, 'int32', 'float64')
        except OverflowError:
            node_id, cp, cd, xyz, ps, seid = self._setup(
                ncards, self.cards, 'int64', 'float64')

        self._save(node_id, cp, cd, xyz, ps, seid)

        grdset = self.model.grdset
        if grdset is not None:
            cp_ = grdset.cp
            cd_ = grdset.cd
            ps_ = grdset.ps
            seid_ = grdset.seid
        else:
            cp_ = 0
            cd_ = 0
            ps_ = 0
            seid_ = 0

        GRID_CP_DEFAULT = GRID_CD_DEFAULT = GRID_PS_DEFAULT = GRID_SEID_DEFAULT = -100
        data_temp_default = [
            (self.cp, GRID_CP_DEFAULT, cp_),
            (self.cd, GRID_CD_DEFAULT, cd_),
            (self.ps, GRID_PS_DEFAULT, ps_),
            (self.seid, GRID_SEID_DEFAULT, seid_),
        ]
        for data, temp_value, default_value in data_temp_default:
            ibad = np.where(data == temp_value)[0]
            if len(ibad):
                data[ibad] = default_value

        assert self.cp.min() >= 0, self.cp.min()

        # -1 is for fluid
        # -3 is for ADAPT/PVAL/PINTS/SINT/GMBNDS
        assert self.cd.min() >= -3, self.cd.min()
        self.sort()
        icp0 = np.where(self.cp == 0)[0]
        self._xyz_cid0[icp0, :] = self.xyz[icp0, :]
        assert self.xyz.shape == self._xyz_cid0.shape
        self.cards = []

    def _setup(self, ncards: int, cards: list[Any],
               idtype: str, fdtype: str):
        node_id = np.zeros(ncards, dtype=idtype)
        cp = np.zeros(ncards, dtype=idtype)
        cd = np.zeros(ncards, dtype=idtype)
        xyz = np.zeros((ncards, 3), dtype=fdtype)
        ps = np.zeros(ncards, dtype=idtype)
        seid = np.zeros(ncards, dtype=idtype)

        for i, card in enumerate(cards):
            (nidi, xyzi, cpi, cdi, psi, seidi, commenti) = card
            node_id[i] = nidi
            cp[i] = cpi
            cd[i] = cdi
            xyz[i, :] = xyzi
            ps[i] = psi
            seid[i] = seidi
            #comment[i] = commenti
        return node_id, cp, cd, xyz, ps, seid

    def _save(self,
              node_id: np.ndarray, cp: np.ndarray, cd: np.ndarray,
              xyz: np.ndarray, ps: np.ndarray, seid: np.ndarray,
              ) -> None:
        ncards = len(node_id)
        ncards_existing = len(self.node_id)

        _xyz_cid0 = np.full((ncards, 3), np.nan, dtype='float64')
        if ncards_existing != 0:
            node_id = np.hstack([self.node_id, node_id])
            cp = np.hstack([self.cp, cp])
            cd = np.hstack([self.cd, cd])
            xyz = np.vstack([self.xyz, xyz])
            ps = np.hstack([self.ps, ps])
            seid = np.hstack([self.seid, seid])
            _xyz_cid0 = np.vstack([self._xyz_cid0, _xyz_cid0])

        self.node_id = node_id
        self.cp = cp
        self.cd = cd
        self.xyz = xyz
        self.ps = ps
        self.seid = seid
        self._xyz_cid0 = _xyz_cid0
        self.n = len(node_id)

        #self.sort()
        #self.cards = []

    def has_ps(self) -> bool:
        return self.n > 0 and self.ps.max() > 0

    #def slice_by_node_id(self, node_id: np.ndarray) -> GRID:
        #inid = self._node_index(node_id)
        #return self.slice_card(inid)

    def slice_card_by_node_id(self, node_id: np.ndarray) -> GRID:
        """uses a node_ids to extract GRIDs"""
        inid = self.index(node_id)
        #assert len(self.node_id) > 0, self.node_id
        #i = np.searchsorted(self.node_id, node_id)
        grid = self.slice_card_by_index(inid)
        return grid

    def slice_card_by_index(self, i: np.ndarray) -> GRID:
        """uses a node_index to extract GRIDs"""
        assert self.xyz.shape == self._xyz_cid0.shape
        assert len(self.node_id) > 0, self.node_id
        i = np.atleast_1d(np.asarray(i, dtype=self.node_id.dtype))
        i.sort()
        grid = GRID(self.model)
        self.__apply_slice__(grid, i)
        return grid

    def __apply_slice__(self, grid: GRID, i: np.ndarray) -> None:
        grid.n = len(i)
        grid._is_sorted = self._is_sorted
        grid.node_id = self.node_id[i]
        grid.cp = self.cp[i]
        grid.cd = self.cd[i]
        grid.ps = self.ps[i]
        grid.seid = self.seid[i]

        grid.xyz = self.xyz[i, :]
        grid._xyz_cid0 = self._xyz_cid0[i, :]
        assert grid.xyz.shape == grid._xyz_cid0.shape

    def get_cp_transform(self) -> dict[int, list[int]]:
        nid_cp_transform = defaultdict(list)
        icp = np.where(self.cp != 0)[0]
        for nid, cp in zip(self.node_id[icp], self.cp[icp]):
            nid_cp_transform[cp].append(nid)
        return dict(nid_cp_transform)

    def get_cd_transform(self) -> dict[int, list[int]]:
        nid_cd_transform = defaultdict(list)
        icd = np.where(self.cd != 0)[0]
        for nid, cd in zip(self.node_id[icd], self.cd[icd]):
            nid_cd_transform[cd].append(nid)
        return dict(nid_cd_transform)

    def get_nid_cp_cd_xyz_transforms(self):
        nnodes = len(self.node_id)
        nid_cp_cd = np.zeros((nnodes, 3), dtype=self.node_id.dtype)
        nid_cp_cd[:, 0] = self.node_id
        nid_cp_cd[:, 1] = self.cp
        nid_cp_cd[:, 2] = self.cd

        #for nid, node in sorted(model.nodes.items()):
            #cd = node.Cd()
            #cp = node.Cp()
            #nids_cp_transform[cp].append(nid)
            #nids_cd_transform[cd].append(nid)
            #nid_cp_cd[i, :] = [nid, cp, cd]
            #xyz_cp[i, :] = node.xyz
            #i += 1
        nids_cp_transform = self.get_cp_transform()
        nids_cd_transform = self.get_cd_transform()
        return nids_cp_transform, nids_cd_transform, nid_cp_cd, self.xyz

    @property
    def _coord_basic(self):
        coord_basic = self.cp.max() == 0 and self.cp.min() == 0 and self.cd.max() == 0 and self.cd.min() == 0
        return coord_basic

    def write(self, size: int=8, is_double: bool=False) -> str:
        if len(self.node_id) == 0:
            return ''
        if size == 8:
            msg = self.write_8()
        elif is_double:
            msg = self.write_double()
        else:
            msg = self.write_16()
        return msg

    def write_8(self):
        coord_basic, extra_basic, is_basic = self._write_flags()
        lines = []
        node_id = array_str(self.node_id, size=8)
        if is_basic:
            cps = ''
            for nid, xyz in zip_longest(node_id, self.xyz):
                #x, y, z = xyz
                msg = 'GRID    %8s%8s%8s%8s%8s\n' % (
                    nid, cps,
                    print_float_8(xyz[0]),
                    print_float_8(xyz[1]),
                    print_float_8(xyz[2]))
                lines.append(msg)
        elif coord_basic: # cp=cd=0
            cps = array_default_int(self.cp, size=8, default=0)
            cds = array_default_int(self.cd, size=8, default=0)
            for nid, cp, cd, xyz in zip_longest(node_id, cps, cds, self.xyz):
                #x, y, z = xyz
                #msg = print_card_8(['GRID', nid, cp, x, y, z, cd])
                msg = 'GRID    %8s%8s%8s%8s%8s%8s\n' % (
                    nid, cp,
                    print_float_8(xyz[0]),
                    print_float_8(xyz[1]),
                    print_float_8(xyz[2]),
                    cd)
                lines.append(msg)
        else:
            cps = array_default_int(self.cp, size=8, default=0)
            cds = array_default_int(self.cd, size=8, default=0)
            pss = array_default_int(self.ps, size=8, default=0)
            seids = array_default_int(self.seid, size=8, default=0)
            for nid, cp, cd, xyz, ps, seid in zip_longest(node_id, cps, cds, self.xyz, pss, seids):
                #msg = print_card_8(['GRID', nid, cp, x, y, z, cd])
                msg = 'GRID    %8s%8s%8s%8s%8s%8s%8s%8s\n' % (
                    nid, cp,
                    print_float_8(xyz[0]),
                    print_float_8(xyz[1]),
                    print_float_8(xyz[2]),
                    cd, ps, seid)
                lines.append(msg)
        return ''.join(lines)

    def _write_flags(self):
        coord_basic = self._coord_basic
        extra_basic = self.ps.max() == 0 and self.ps.min() == 0 and self.seid.max() == self.seid.min()
        is_basic = coord_basic and extra_basic
        return coord_basic, extra_basic, is_basic

    def write_16(self) -> str:
        return _write_grid_large(self, print_scientific_16)

    def write_double(self) -> str:
        return _write_grid_large(self, print_scientific_double)

    def is_equal_by_index(self, inode1, inode2) -> bool:
        is_equal = True
        for key in ['node_id', 'cp', 'cd', 'ps', 'seid', 'xyz']:
            values = getattr(self, key)
            val1 = values[inode1]
            val2 = values[inode2]
            if not np.array_equal(val1, val2):
                print(f'key={key!r} are not equal')
                is_equal = False
                break
        return is_equal
        #self.node_id = np.array([], dtype='int32')
        #self.cp = np.array([], dtype='int32')
        #self.cd = np.array([], dtype='int32')
        #self.ps = np.array([], dtype='int32')
        #self.seid = np.array([], dtype='int32')
        #self.xyz = np.zeros((0, 3), dtype='float64')

    def sort(self) -> None:
        # get the sort index
        ids = self.node_id
        iarg = np.argsort(ids)

        # sort the sort index
        uarg = np.unique(iarg)

        ids_sorted = ids[iarg]
        uids_sorted = np.unique(ids_sorted)

        if not np.array_equal(uarg, iarg):
            print(iarg.tolist())
            print('->  ', ids_sorted.tolist())
            print('--> ', uids_sorted.tolist())
            assert (iarg - uarg).sum() == 0, (iarg - uarg).sum()
            #if len(iarg) == len(uarg):
                ## if the lengths are the same, we can use dumb sorting

                ## check on dumb sort; no missing entries and no duplicates
                #assert (iarg - uarg).sum() == 0, (iarg - uarg).sum()

            if len(ids_sorted) != len(uids_sorted):
                # Now that we've sorted the ids, we'll take the delta id with the next id.
                # We check the two neighoring values when there is a duplicate,
                # so if they're they same in xyz, cp, cd, etc., the results are the same
                # we don't need to check unique ids
                dnode = ids_sorted[1:] - ids_sorted[:-1] # high - low
                #dnode = np.diff(ids_sorted)
                #inode = (dnode == 0)
                inode1 = np.hstack([False, dnode == 0])
                inode2 = np.hstack([dnode == 0, False])
                #ids_unique1 = ids_sorted[inode1] # high node
                #ids_unique2 = ids_sorted[inode2] # low  node

                if not self.is_equal_by_index(iarg[inode1], iarg[inode2]):
                    raise RuntimeError('bad duplicates')

                #ids_sorted[inode]

                #dnode2 = ids_sorted[:-1] - ids_sorted[1:] # low-high
                #inode2 = np.hstack([dnode2 == 0, False])

                self.__apply_slice__(self, iarg[inode1])
            else:
                # check if the sort index is sorted
                # True: no sorting is required
                # False: we need to sort
                self.node_id = self.node_id[iarg]
                self.cp = self.cp[iarg]
                self.xyz = self.xyz[iarg, :]
                self.ps = self.ps[iarg]
                self.cd = self.cd[iarg]
                self.seid = self.seid[iarg]
        unid = np.unique(self.node_id)

        if len(self.node_id) != len(unid):
            # we need to filter nodes
            msg = f'len(self.node_id)={len(self.node_id)} != len(unid)={len(unid)}'
            print(msg)
            #x = 1

        #self.remove_duplicates(inplace=False)
        self._is_sorted = True

    #def get_global_position_by_node_id(self, nids: np.ndarray) -> np.ndarray:
        #xyz_cid0 = self.xyz_cid0()
        #inid = np.searchsorted(self.node_id, nids)
        #return xyz_cid0[inid, :]

    def get_position_by_node_id(self, node_id: np.ndarray) -> np.ndarray:
        """returns the global position of the node"""
        xyz_cid0 = self._xyz_cid0
        inid = self.index(node_id)
        return xyz_cid0[inid, :]

    #def index(self, node_id: np.ndarray, safe: bool=False) -> np.ndarray:
        #assert len(self.node_id) > 0, self.node_id
        #node_id = np.atleast_1d(np.asarray(node_id, dtype=self.node_id.dtype))
        #inid = np.searchsorted(self.node_id, node_id)
        #if safe:
            #ibad = inid >= len(self.node_id)
            #if sum(ibad):
                ##self.model.log.error(f'bad nids; node_id={node_id[ibad]}')
                #raise RuntimeError(f'bad nids; node_id={node_id[ibad]}')
            #inids_leftover = inid[~ibad]
            #if len(inids_leftover):
                #actual_nids = self.node_id[inids_leftover]
                #assert np.array_equal(actual_nids, node_id)
        #return inid

    def xyz_cid0(self):
        """not validated"""
        if not self._is_sorted:
            self.sort()
        xyz_cid0 = self._xyz_cid0


        xyz = self.xyz
        #self.model.log.info(f'xyz = {xyz_cid0.shape}')
        assert xyz.shape == xyz_cid0.shape, f'xyz.shape={xyz.shape} xyz_cid0.shape={xyz_cid0.shape}'
        #icp = (self.cp > 0)
        #cps_to_transform = self.cp[icp]
        # TODO: filter?
        ucp = np.unique(self.cp)

        coord = self.model.coord
        coord_types = coord.coord_type
        coord_ids = coord.coord_id
        origins = coord.origin
        iaxes = coord.i
        jaxes = coord.j
        kaxes = coord.k
        for cp in ucp:
            if cp == 0:
                continue
            icp = np.where(self.cp == cp)[0]
            #nnodes = len(icp)

            icid = np.where(coord_ids == cp)[0]
            assert len(icid) == 1, f'CORDx={cp} not found; referenced by CP field of GRID'
            #coord_type = np.asscalar(coord_types[icid])  #deprecated in numpy 1.16
            coord_type = coord_types[icid].item()

            #assert coord_type in {'R'}, coord_type
            assert coord_type in {'R', 'C', 'S'}, coord_type  # output coord
            origin = origins[icid, :]
            iaxis = iaxes[icid, :]
            jaxis = jaxes[icid, :]
            kaxis = kaxes[icid, :]
            beta = np.vstack([iaxis, jaxis, kaxis])
            assert origin.shape == (1, 3), origin.shape
            #print('--------------------------')
            #print(f'cp = {cp}')
            #print(f'beta = {beta}')
            #print(f'origin = {origin}')
            if 0:
                xyz_cid = xyz[icp, :] @ beta + origin
                if coord_type == 'R':
                    xyz_cid0i = xyz_cid
                elif coord_type == 'C':
                    xyz_cid0i = cylindrical_to_rectangular(xyz_cid)
                elif coord_type == 'S':
                    xyz_cid0i = spherical_to_rectangular(xyz_cid)
                else:
                    raise NotImplementedError(coord_type)
            else:
                # transform node from CORD2C into a rectangular frame
                # then transform it to the global
                if coord_type == 'R':
                    xyz_cidr = xyz[icp, :]
                elif coord_type == 'C':
                    xyz_cidr = cylindrical_to_rectangular(xyz[icp, :])
                elif coord_type == 'S':
                    xyz_cidr = spherical_to_rectangular(xyz[icp, :])
                else:
                    raise NotImplementedError(coord_type)
                xyz_cid0i = xyz_cidr @ beta + origin
            xyz_cid0[icp, :] = xyz_cid0i
        assert not np.any(np.isnan(xyz_cid0))
        assert xyz_cid0.shape[0] > 0, xyz_cid0.shape
        return xyz_cid0


def _write_grid_large(grid: GRID, function: Callable[float]) -> str:
    coord_basic, extra_basic, is_basic = grid._write_flags()
    lines = []
    node_id = array_str(grid.node_id, size=16)
    if is_basic:
        for nid, xyz in zip_longest(node_id, grid.xyz):
            msg = ('GRID*   %16s                %16s%16s\n'
                   '*       %16s\n' % (
                       nid,
                       print_scientific_16(xyz[0]),
                       print_scientific_16(xyz[1]),
                       print_scientific_16(xyz[2])))
            lines.append(msg)
    elif coord_basic: # cp=cd=0
        cps = array_default_int(grid.cp, size=16, default=0)
        cds = array_default_int(grid.cd, size=16, default=0)
        for nid, cp, cd, xyz in zip_longest(node_id, cps, cds, grid.xyz):
            msg = ('GRID*   %16s%16s%16s%16s\n'
                   '*       %16s%16s\n' % (
                       nid,
                       cp,
                       print_scientific_16(xyz[0]),
                       print_scientific_16(xyz[1]),
                       print_scientific_16(xyz[2]),
                       cd))
            lines.append(msg)
    else:
        cps = array_default_int(grid.cp, size=16, default=0)
        cds = array_default_int(grid.cd, size=16, default=0)
        pss = array_default_int(grid.ps, size=16, default=0)
        seids = array_default_int(grid.seid, size=16, default=0)
        for nid, cp, cd, xyz, ps, seid in zip_longest(node_id, cps, cds, grid.xyz, pss, seids):
            msg = ('GRID*   %16s%16s%16s%16s\n'
                   '*       %16s%16s%16s%16s\n' % (
                       nid,
                       cp,
                       print_scientific_16(xyz[0]),
                       print_scientific_16(xyz[1]),
                       print_scientific_16(xyz[2]),
                       cd, ps, seid))
            lines.append(msg)
    return ''.join(lines)
