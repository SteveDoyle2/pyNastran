from __future__ import print_function
from collections import defaultdict
import numpy as np
from pyNastran.bdf.bdf import GRID

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank,
    components_or_blank)
from pyNastran.bdf.field_writer_8 import print_float_8, set_string8_blank_if_default
from pyNastran.bdf.field_writer_16 import print_float_16, set_string16_blank_if_default
from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.base_card import _format_comment


class GRIDv(object):
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
        self.is_current = False
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

    def add(self, nid, xyz, cp=0, cd=0, ps='', seid=0, comment=''):
        # type: (int, Union[None, List[float], np.ndarray], int, int, str, int, str) -> None
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

    def add_card(self, card, comment=''):
        # type: (Any, str) -> GRID
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

    def _make_current(self):
        """creates an array of the GRID points"""
        if not self.is_current:
            if len(self.nid) > 0: # there are already nodes in self.nid
                self.nid = np.hstack([self.nid, self._nid])
                self.xyz = np.vstack([self.xyz, self._xyz])
                self.cp = np.hstack([self.cp, self._cp])
                self.cd = np.hstack([self.cd, self._cd])
                self.ps = np.hstack([self.ps, self._ps])
                self.seid = np.hstack([self.seid, self._seid])
                # don't need to handle comments
            else:
                self.nid = np.array(self._nid)
                self.xyz = np.array(self._xyz)
                self.cp = np.array(self._cp)
                self.cd = np.array(self._cd)
                self.ps = np.array(self._ps)
                self.seid = np.array(self._seid)
            assert len(self.nid) == len(np.unique(self.nid))

            isort = np.argsort(self.nid)
            self.nid = self.nid[isort]
            self.xyz = self.xyz[isort, :]
            self.cp = self.cp[isort]
            self.cd = self.cd[isort]
            self.ps = self.ps[isort]
            self.seid = self.seid[isort]

            #print(self.nid)
            self._nid = []
            self._xyz = []
            self._cp = []
            self._cd = []
            self._ps = []
            self._seid = []
            self.is_current = True
        #else:
            #print('no GRIDs')

    def cross_reference(self, model):
        """does this do anything?"""
        self._make_current()

    def write_card(self, size=8, is_double=False):
        # type: (int, bool) -> str
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
        if size == 8:
            return self.write_card_8()
        else:
            return self.write_card_16(is_double)

    def write_card_8(self):
        # type: () -> str
        """
        Writes a GRID card in 8-field format
        """
        self._make_current()
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
        return msg

    def write_card_16(self, is_double=False):
        # type: (bool) -> str
        """
        Writes a GRID card in 16-field format
        """
        self._make_current()
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
        return msg

    def __len__(self):
        self._make_current()
        return len(self.nid)

    def __repr__(self):
        self._make_current()
        msg = 'GRID_Vector; ngrids=%s:\n' % len(self.nid)
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
        self._make_current()
        nid = self.nid[i]
        return GRID(nid, self.xyz[i], cp=self.cp[i], cd=self.cd[i],
                    ps=self.ps[i], seid=self.seid[i], comment=self.comment[nid])

    def get_grid_by_nid(self, nid):
        self._make_current()
        inid = np.searchsorted(self.nid, nid)
        return self[inid]

    def __setitem__(self, i, value):
        pass
    #def __delitem__(self, i):
        #pass


