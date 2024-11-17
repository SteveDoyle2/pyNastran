from __future__ import annotations
from typing import TYPE_CHECKING
from collections import defaultdict
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class Rods:
    """initializes the Rods"""
    def __init__(self, model):
        self.model = model
        self.conrod = model.conrod
        self.crod = model.crod
        self.ctube = model.ctube
        self._eids = set()

    def add(self, eid):
        if eid not in self._eids:
            self._eids.add(eid)
        else:
            raise RuntimeError('eid=%s is duplicated' % eid)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.conrod):
            self.conrod.write_card(size, is_double, bdf_file)
        if len(self.crod):
            self.crod.write_card(size, is_double, bdf_file)
        if len(self.ctube):
            self.ctube.write_card(size, is_double, bdf_file)

    def make_current(self):
        self.conrod.make_current()
        self.crod.make_current()
        self.ctube.make_current()

    def __len__(self):
        return len(self.conrod) + len(self.crod) + len(self.ctube)

    def repr_indent(self, indent='  '):
        msg = '%s<Rods> : nelements=%s\n' % (indent, len(self))
        msg += '%s  CONROD:  %s\n' % (indent, len(self.conrod))
        msg += '%s  CROD  :  %s\n' % (indent, len(self.crod))
        msg += '%s  CTUBE :  %s\n' % (indent, len(self.ctube))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')


class RodElement:
    """base class for CONROD, CROD, and CTUBE"""
    card_name = ''

    def check_if_current(self, eid, eids):
        """we split this up to reason about it easier"""
        if self.is_current:
            if eid in eids:
                # card exists, so we use that slot
                add_card = False
            else:
                add_card = True
        else:
            add_card = True
        return add_card

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
        msg = '%s%sv; nelements=%s\n' % (indent, self.card_name, neids)
        msg += '%s  eid = %s\n' % (indent, self.eid)

        if hasattr(self, 'pid'):
            upid = np.unique(self.pid)
            if len(upid) == 1:
                msg += '%s  upid = %s\n' % (indent, upid)
            else:
                msg += '%s  pid = %s\n' % (indent, self.pid)
        else:
            msg += '%s  A = %s\n' % (indent, self.A)
            msg += '%s  j = %s\n' % (indent, self.j)
            msg += '%s  c = %s\n' % (indent, self.c)
            msg += '%s  nsm = %s\n' % (indent, self.nsm)
        return msg

        #umcid = np.unique(self.mcid)
        #if len(umcid) == 1 and umcid[0] == 0:
            #msg += '  umcid = %s\n' % umcid
        #else:
            #msg += '  umcid = %s\n' % umcid
            #msg += '  mcid = %s\n' % self.mcid

        #utheta = np.unique(self.theta)
        #if len(utheta) == 1 and umcid[0] == 0:
            #msg += '  utheta = %s\n' % utheta
        #else:
            #msg += '  theta = %s\n' % self.theta
        #msg += '  is_theta = %s\n' % self.is_theta
        #msg += '  nid =\n%s' % self.nid
        #return msg

    def __repr__(self):
        return self.repr_indent(indent='')


class CONRODv(RodElement):
    """
    +--------+-----+-----+----+-----+---+---+---+-----+
    |   1    |  2  |  3  |  4 |  5  | 6 | 7 | 8 |  9  |
    +========+=====+=====+====+=====+===+===+===+=====+
    | CONROD | EID | N1  | N2 | MID | A | J | C | NSM |
    +--------+-----+-----+----+-----+---+---+---+-----+
    """
    card_name = 'CONROD'

    def __init__(self, model):
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')
        self.mid = np.array([], dtype='int32')
        self.A = np.array([], dtype='float64')
        self.j = np.array([], dtype='float64')
        self.c = np.array([], dtype='float64')
        self.nsm = np.array([], dtype='float64')

        self._eid = []
        self._nids = []
        self._mid = []
        self._A = []
        self._j = []
        self._c = []
        self._nsm = []
        self.comment = defaultdict(str)

    def add(self, eid, mid, nids, A=0., j=0., c=0., nsm=0., comment=''):
        """
        Creates a CONROD card

        Parameters
        ----------
        eid : int
            element id
        mid : int
            material id
        nids : list[int, int]
            node ids
        A : float
            area
        j : float; default=0.
            polar moment of inertia
        c : float; default=0.
            stress factor
        nsm : float; default=0.
            non-structural mass per unit length
        comment : str; default=''
            a comment for the card
        """
        self.model.rods.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._nids.append(nids)
        self._mid.append(mid)
        self._A.append(A)
        self._j.append(j)
        self._c.append(c)
        self._nsm.append(nsm)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a CONROD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        nids = [integer(card, 2, 'n1'),
                integer(card, 3, 'n2')]
        mid = integer(card, 4, 'mid')
        A = double_or_blank(card, 5, 'A', 0.0)
        j = double_or_blank(card, 6, 'j', 0.0)
        c = double_or_blank(card, 7, 'c', 0.0)
        nsm = double_or_blank(card, 8, 'nsm', 0.0)
        self.add(eid, mid, nids, A, j, c, nsm, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, mid, nodes, A, j, c, nsm in zip(self.eid, self.mid,
                                                 self.nids, self.A, self.j, self.c, self.nsm):
            j = set_blank_if_default(j, 0.0)
            c = set_blank_if_default(c, 0.0)
            nsm = set_blank_if_default(nsm, 0.0)
            list_fields = ['CONROD', eid] + nodes.tolist() + [mid, A, j, c, nsm]
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.mid = np.vstack([self.mid, self._mid])
                self.nids = np.hstack([self.nids, self._nids])

                self.A = np.hstack([self.A, self._A])
                self.j = np.hstack([self.j, self._j])
                self.c = np.hstack([self.c, self._c])
                self.nsm = np.hstack([self.nsm, self._nsm])

                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.mid = np.array(self._mid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
                self.A = np.array(self._A, dtype='float64')
                self.j = np.array(self._j, dtype='float64')
                self.c = np.array(self._c, dtype='float64')
                self.nsm = np.array(self._nsm, dtype='float64')
            assert len(self.eid) == len(np.unique(self.eid))
            self._eid = []
            self._mid = []
            self._nids = []
            self._A = []
            self._j = []
            self._c = []
            self._nsm = []
            self.is_current = True


class CRODv(RodElement):
    """
    +------+-----+-----+----+----+
    |   1  |  2  |  3  |  4 |  5 |
    +======+=====+=====+====+====+
    | CROD | EID | PID | N1 | N2 |
    +------+-----+-----+----+----+
    """
    card_name = 'CROD'

    def __init__(self, model):
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')

        self._eid = []
        self._pid = []
        self._nids = []
        self._dofs = []
        self.comment = defaultdict(str)

    def add(self, eid, pid, nids, comment=''):
        """
        Creates a CROD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PROD)
        nids : list[int, int]
            node ids
        comment : str; default=''
            a comment for the card
        """
        self.model.solids.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a CROD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2')]
        assert len(card) == 5, 'len(CROD card) = %i\ncard=%s' % (len(card), str(card))
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            list_fields = ['CROD', eid, pid] + nodes.tolist()
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
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
            self._eid = []
            self._pid = []
            self._nids = []
            self.is_current = True


class CTUBEv(RodElement):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CELAS3 | EID | PID | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    card_name = 'CTUBE'

    def __init__(self, model):
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')

        self._eid = []
        self._pid = []
        self._nids = []
        self.comment = defaultdict(str)

    def add(self, eid, pid, nids, comment=''):
        """
        Creates a CTUBE card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id
        nids : list[int, int]
            node ids
        comment : str; default=''
            a comment for the card
        """
        self.model.rods.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a CTUBE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2')]
        assert len(card) == 5, 'len(CTUBE card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            list_fields = ['CTUBE', eid, pid, nodes[0], nodes[1]]
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
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
            self._eid = []
            self._pid = []
            self._nids = []
            self.is_current = True
