from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class Springs:
    """intializes the Springs"""
    def __init__(self, model):
        self.model = model
        self.celas1 = model.celas1
        self.celas2 = model.celas2
        self.celas3 = model.celas3
        self.celas4 = model.celas4
        self._eids = set()

    def add(self, eid):
        if eid not in self._eids:
            self._eids.add(eid)
        else:
            raise RuntimeError('eid=%s is duplicated' % eid)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.celas1):
            self.celas1.write_card(size, is_double, bdf_file)
        if len(self.celas2):
            self.celas2.write_card(size, is_double, bdf_file)
        if len(self.celas3):
            self.celas3.write_card(size, is_double, bdf_file)
        if len(self.celas4):
            self.celas4.write_card(size, is_double, bdf_file)

    def make_current(self):
        self.celas1.make_current()
        self.celas2.make_current()
        self.celas3.make_current()
        self.celas4.make_current()

    def __len__(self):
        return(len(self.celas1) + len(self.celas2) +
               len(self.celas3) + len(self.celas4))

    def repr_indent(self, indent='  '):
        msg = '%s<Springs> : nelements=%s\n' % (indent, len(self))
        msg += '%s  CELAS1:  %s\n' % (indent, len(self.celas1))
        msg += '%s  CELAS2:  %s\n' % (indent, len(self.celas2))
        msg += '%s  CELAS3:  %s\n' % (indent, len(self.celas3))
        msg += '%s  CELAS4:  %s\n' % (indent, len(self.celas4))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')


class SpringElement:
    """base class for CELAS1, CELAS2, CELAS3, and CELAS4"""
    card_name = ''

    def __init__(self, model):
        self.model = model
        self.is_current = True

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
            msg += '%s  k = %s\n' % (indent, self.k)

        if hasattr(self, 's'):
            msg += '%s  s = %s\n' % (indent, self.s)
            msg += '%s  ge = %s\n' % (indent, self.ge)
        return msg.rstrip()

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


class CELAS1(SpringElement):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CELAS1 | EID | PID | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+
    """
    card_name = 'CELAS1'

    def __init__(self, model):
        SpringElement.__init__(self, model)
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')
        self.dofs = np.array([], dtype='int32')

        self._eid = []
        self._pid = []
        self._nids = []
        self._dofs = []
        self.comment = defaultdict(str)

    def add(self, eid, pid, nids, dofs, comment=''):
        """
        Creates a CELAS1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : List[int]
            GRID/SPOINT ids; n=2
        dofs : List[int]
            node ids; n=2
        comment : str; default=''
            a comment for the card
        """
        self.model.solids.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        self._dofs.append(dofs)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CELAS1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'g1'), integer_or_blank(card, 5, 'g2', 0)]

        #: component number
        c1 = integer_or_blank(card, 4, 'c1', 0)
        c2 = integer_or_blank(card, 6, 'c2', 0)
        assert len(card) <= 7, 'len(CELAS1 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, [c1, c2], comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes, dofs in zip(self.eid, self.pid, self.nids, self.dofs):
            list_fields = ['CELAS1', eid, pid, nodes[0], dofs[0], nodes[1], dofs[1]]
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
                self.dofs = np.hstack([self.dofs, self._dofs])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.pid = np.array(self._pid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
                self.dofs = np.array(self._dofs, dtype='int32')
            assert len(self.eid) == len(np.unique(self.eid))
            assert self.nids[:, 0].min() >= 0, self.nids
            self._eid = []
            self._pid = []
            self._nids = []
            self._dofs = []
            self.is_current = True


class CELAS2(SpringElement):
    """
    +--------+-----+-----+----+----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |  8 |  9 |
    +========+=====+=====+====+====+====+====+====+====+
    | CELAS2 | EID |  K  | G1 | C1 | G2 | C2 | GE | S  |
    +--------+-----+-----+----+----+----+----+----+----+
    """
    card_name = 'CELAS2'

    def __init__(self, model):
        SpringElement.__init__(self, model)
        self.eid = np.array([], dtype='int32')
        self.k = np.array([], dtype='float64')
        self.nids = np.array([], dtype='int32')
        self.dofs = np.array([], dtype='int32')
        self.ge = np.array([], dtype='float64')
        self.s = np.array([], dtype='float64')

        self._eid = []
        self._k = []
        self._nids = []
        self._dofs = []
        self._ge = []
        self._s = []
        self.comment = defaultdict(str)

    def add(self, eid, k, nids, dofs, ge, s, comment=''):
        """
        Creates a CELAS1 card

        Parameters
        ----------
        eid : int
            element id
        k : float
            spring stiffness
        nids : List[int]
            GRID/SPOINT ids; n=2
        dofs : List[int]
            node ids; n=2
        ge : float
            ???
        s : float
            ???
        comment : str; default=''
            a comment for the card
        """
        self.model.solids.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._k.append(k)
        self._nids.append(nids)
        self._dofs.append(dofs)
        self._ge.append(ge)
        self._s.append(s)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CELAS2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        k = double(card, 2, 'k')
        nids = [integer_or_blank(card, 3, 'g1', 0),
                integer_or_blank(card, 5, 'g2', 0)]
        c1 = integer_or_blank(card, 4, 'c1', 0)
        c2 = integer_or_blank(card, 6, 'c2', 0)
        ge = double_or_blank(card, 7, 'ge', 0.)
        s = double_or_blank(card, 8, 's', 0.)
        assert len(card) <= 9, 'len(CELAS2 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, k, nids, [c1, c2], ge, s, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, k, nodes, dofs, ge, s in zip(self.eid, self.k, self.nids, self.dofs, self.ge, self.s):
            list_fields = ['CELAS2', eid, k, nodes[0], dofs[0], nodes[1], dofs[1], ge, s]
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.k = np.vstack([self.k, self._k])
                self.nids = np.hstack([self.nids, self._nids])
                self.dofs = np.hstack([self.dofs, self._dofs])
                self.ge = np.hstack([self.ge, self._ge])
                self.s = np.hstack([self.s, self._s])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.k = np.array(self._k, dtype='float64')
                self.nids = np.array(self._nids, dtype='int32')
                self.dofs = np.array(self._dofs, dtype='int32')
                self.ge = np.array(self._ge, dtype='float64')
                self.s = np.array(self._s, dtype='float64')
            assert len(self.eid) == len(np.unique(self.eid))
            assert self.nids[:, 0].min() >= 0, self.nids
            self._eid = []
            self._k = []
            self._nids = []
            self._dofs = []
            self._ge = []
            self._s = []
            self.is_current = True


class CELAS3(SpringElement):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CELAS3 | EID | PID | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    card_name = 'CELAS3'

    def __init__(self, model):
        SpringElement.__init__(self, model)
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')

        self._eid = []
        self._pid = []
        self._nids = []
        self.comment = defaultdict(str)

    def add(self, eid, pid, nids, comment=''):
        """
        Creates a CELAS3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : List[int]
            SPOINT ids; n=2
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

    def add_card(self, card, comment=''):
        """
        Adds a CELAS3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)

        s1 = integer_or_blank(card, 3, 's1', 0)
        s2 = integer_or_blank(card, 4, 's2', 0)
        assert len(card) <= 5, 'len(CELAS3 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, [s1, s2], comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            list_fields = ['CELAS3', eid, pid, nodes[0], nodes[1]]
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
            assert self.nids[:, 0].min() >= 0, self.nids
            self._eid = []
            self._pid = []
            self._nids = []
            self.is_current = True


class CELAS4(SpringElement):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CELAS4 | EID |  K  | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    card_name = 'CELAS4'

    def __init__(self, model):
        SpringElement.__init__(self, model)
        self.eid = np.array([], dtype='int32')
        self.k = np.array([], dtype='float64')
        self.nids = np.array([], dtype='int32')

        self._eid = []
        self._k = []
        self._nids = []
        self.comment = defaultdict(str)

    def add(self, eid, k, nids, comment=''):
        """
        Creates a CELAS4 card

        Parameters
        ----------
        eid : int
            element id
        k : float
            spring stiffness
        nids : List[int]
            SPOINT ids; n=2
        comment : str; default=''
            a comment for the card
        """
        self.model.springs.add(eid)
        self.is_current = False
        self._eid.append(eid)
        assert isinstance(k, float), k
        self._k.append(k)
        self._nids.append(nids)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CELAS4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        k = double(card, 2, 'k')
        s1 = integer_or_blank(card, 3, 's1', 0)
        s2 = integer_or_blank(card, 4, 's2', 0)
        assert len(card) <= 5, 'len(CELAS4 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, k, [s1, s2], comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, k, nodes in zip(self.eid, self.k, self.nids):
            list_fields = ['CELAS4', eid, k, nodes[0], nodes[1]]
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.k = np.vstack([self.k, self._k])
                self.nids = np.hstack([self.nids, self._nids])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.k = np.array(self._k, dtype='float64')
                self.nids = np.array(self._nids, dtype='int32')
            assert len(self.eid) == len(np.unique(self.eid))
            assert self.nids[:, 0].min() >= 0, self.nids
            self._eid = []
            self._k = []
            self._nids = []
            self.is_current = True
