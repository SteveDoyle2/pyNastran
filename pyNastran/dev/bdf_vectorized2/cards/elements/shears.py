from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank)
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class ShearElement:
    """base class for CSHEAR"""
    card_name = ''
    def __init__(self, model):
        """initializes the ShearElement"""
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='float64')

        self._eid = []
        self._pid = []
        self._nids = []
        self.comment = defaultdict(str)

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

    #def get_element_by_eid(self, eid):
        #self.make_current()
        #ieid = np.searchsorted(eid, self.eid)
        #return self[ieid]

    def make_current(self):
        """creates an array of the GRID points"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.pid = np.vstack([self.pid, self._pid])
                self.nids = np.hstack([self.nids, self._nids])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.pid = np.array(self._pid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
            assert len(self.eid) == len(np.unique(self.eid))

            isort = np.argsort(self.eid)
            self.eid = self.eid[isort]
            self.pid = self.pid[isort]
            self.nids = self.nids[isort, :]

            self._eid = []
            self._pid = []
            self._nids = []
            self.is_current = True

    def cross_reference(self, model: BDF) -> None:
        """does this do anything?"""
        self.make_current()

    def __len__(self):
        """returns the number of elements"""
        return len(self.eid) + len(self._eid)

    def repr_indent(self, indent=''):
        self.make_current()
        neids = len(self.eid)
        if neids == 0:
            return '%s%sv; nelements=%s' % (indent, self.card_name, neids)
        msg = '%s%sv; nelements=%s:\n' % (indent, self.card_name, neids)
        msg += '%s  eid = %s\n' % (indent, self.eid)

        upid = np.unique(self.pid)
        if len(upid) == 1 and upid[0] == 0:
            msg += '%s  upid = %s\n' % (indent, upid)
        else:
            msg += '%s  pid = %s\n' % (indent, self.pid)

        #msg += '  nid =\n%s' % self.nid
        return msg

    def __repr__(self):
        return self.repr_indent('')


class CSHEARv(ShearElement):
    """
    +--------+-------+-------+----+----+----+----+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |
    +========+=======+=======+=====+===+====+====+
    | CSHEAR |  EID  |  PID  | N1 | N2 | N3 | N4 |
    +--------+-------+-------+----+----+----+----+
    """
    card_name = 'CSHEAR'
    nnodes = 4
    nthickness = 4

    def add(self, eid, pid, nids, comment=''):
        """
        Creates a CSHEAR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHEAR)
        nids : list[int, int, int, int]
            node ids
        comment : str; default=''
            a comment for the card
        """
        self.model.shears.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CSHEAR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer_or_blank(card, 3, 'n1'),
                integer_or_blank(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4')]
        assert len(card) <= 7, 'len(CSHEAR card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids)

    #def update(self, grid):
        #"""functions like a dictionary"""
        #nid = grid.nid
        #add_card = self.check_if_current(eid, self.eid)
        #if add_card:
            #self.add(nid, grid.xyz, cp=grid.cp, cd=grid.cd,  # add_cquad4
                     #ps=grid.ps, seid=grid.seid, comment=grid.comment)
            #self.is_current = False
        #else:
            #inid = np.where(nid == self.nid)[0]
            #self.nid[inid] = grid.nid
            #self.xyz[inid] = grid.xyz
            #self.cp[inid] = grid.cp
            #self.cd[inid] = grid.cd
            #self.ps[inid] = grid.ps
            #self.seid[inid] = grid.seid
            #self.comment[nid] = comment
            #self.is_current = True  # implicit

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
    #def __getitem__(self, i):
        #"""this works on index"""
        #self.make_current()
        #eid = self.eid[i]
        #return GRID(nid, self.xyz[i], cp=self.cp[i], cd=self.cd[i],
                    #ps=self.ps[i], seid=self.seid[i], comment=self.comment[nid])

    #def __setitem__(self, i, value):
        #pass
    #def __delitem__(self, i):
        #pass
    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            data = [eid, pid] + nodes.tolist()
            msgi = ('CSHEAR  %8i%8i%8i%8i%8i%8i\n' % tuple(data))
            msg += self.comment[eid] + msgi
        bdf_file.write(msg)
        return msg


class Shears:
    """
    Stores CSHEAR elements that exist in 3D space
    """
    def __init__(self, model):
        self.model = model
        self.cshear = model.cshear
        self._eids = set()

    def add(self, eid):
        if eid not in self._eids:
            self._eids.add(eid)
        else:
            raise RuntimeError('eid=%s is duplicated' % eid)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.cshear):
            self.cshear.write_card(size, is_double, bdf_file)

    def make_current(self):
        self.cshear.make_current()

    def __len__(self):
        return len(self.cshear)

    def repr_indent(self, indent=''):
        msg = '%s<Shears> : nelements=%s\n' % (indent, len(self))
        msg += '%s  CSHEAR: %s\n' % (indent, len(self.cshear))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')

