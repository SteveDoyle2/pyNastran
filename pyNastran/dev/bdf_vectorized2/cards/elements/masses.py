from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class Masses:
    """initializes the Masses"""
    def __init__(self, model):
        self.model = model
        #self.conm1 = model.conm1
        self.conm2 = model.conm2
        #self.cmass1 = model.cmass1
        #self.cmass2 = model.cmass2
        #self.cmass3 = model.cmass3
        #self.cmass4 = model.cmass4
        self._eids = set()

    def add(self, eid):
        if eid not in self._eids:
            self._eids.add(eid)
        else:
            raise RuntimeError('eid=%s is duplicated' % eid)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        #if len(self.conm1):
            #self.conm1.write_card(size, is_double, bdf_file)
        if len(self.conm2):
            self.conm2.write_card(size, is_double, bdf_file)
        #if len(self.cmass1):
            #self.cmass1.write_card(size, is_double, bdf_file)
        #if len(self.cmass2):
            #self.cmass2.write_card(size, is_double, bdf_file)
        #if len(self.cmass3):
            #self.cmass3.write_card(size, is_double, bdf_file)
        #if len(self.cmass4):
            #self.cmass4.write_card(size, is_double, bdf_file)

    def make_current(self):
        self.conm1.make_current()
        self.conm2.make_current()
        self.cmass1.make_current()
        self.cmass2.make_current()
        self.cmass3.make_current()
        self.cmass4.make_current()

    @property
    def groups(self):
        """gets the sub-element groups"""
        groups = [
            #self.conm1,
            self.conm2,
            #self.cmass1,
            #self.cmass2,
            #self.cmass3,
            #self.cmass4,
        ]
        return groups

    @property
    def elements(self):
        """gets all the elements in the expected sorted order"""
        elems = []
        for group in self.groups:
            elems += group.elements
        return elems

    def __len__(self):
        """gets the number of elements"""
        return sum([len(group) for group in self.groups])

    def repr_indent(self, indent='  '):
        msg = '%s<Masses> : nelements=%s\n' % (indent, len(self))
        #msg += '%s  CONM1 :  %s\n' % (indent, len(self.conm1))
        msg += '%s  CONM2 :  %s\n' % (indent, len(self.conm2))
        #msg += '%s  CMASS1:  %s\n' % (indent, len(self.cmass1))
        #msg += '%s  CMASS2:  %s\n' % (indent, len(self.cmass2))
        #msg += '%s  CMASS3:  %s\n' % (indent, len(self.cmass3))
        #msg += '%s  CMASS4:  %s\n' % (indent, len(self.cmass4))
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


class CONM2v(RodElement):
    """
    +-------+--------+-------+-------+---------+------+------+------+
    |   1   |    2   |    3  |   4   |    5    |  6   |  7   |   8  |
    +=======+========+=======+=======+=========+======+======+======+
    | CONM2 |   EID  |  NID  |  CID  |  MASS   |  X1  |  X2  |  X3  |
    +-------+--------+-------+-------+---------+------+------+------+
    |       |   I11  |  I21  |  I22  |   I31   |  I32 |  I33 |      |
    +-------+--------+-------+-------+---------+------+------+------+
    | CONM2 | 501274 | 11064 |       | 132.274 |      |      |      |
    +-------+--------+-------+-------+---------+------+------+------+
    """
    card_name = 'CONM2'

    def __init__(self, model):
        self.model = model
        self.is_current = False
        self.eid = np.array([], dtype='int32')
        self.nid = np.array([], dtype='int32')
        self.cid = np.array([], dtype='int32')
        self.x = np.array([], dtype='float64')
        self.I = np.array([], dtype='float64')
        self.mass = np.array([], dtype='float64')

        self._eid = []
        self._nid = []
        self._cid = []
        self._x = []
        self._I = []
        self._mass = []
        self.comment = defaultdict(str)

    def add(self, eid, nid, mass, cid=0, X=None, I=None, comment=''):
        """
        Creates a CONM2 card

        Parameters
        ----------
        eid : int
           element id
        nid : int
           node id
        mass : float
           the mass of the CONM2
        cid : int; default=0
           coordinate frame of the offset (-1=absolute coordinates)
        X : (3, ) list[float]; default=None -> [0., 0., 0.]
            xyz offset vector relative to nid
        I : (6, ) list[float]; default=None -> [0., 0., 0., 0., 0., 0.]
            mass moment of inertia matrix about the CG
            I11, I21, I22, I31, I32, I33 = I
        comment : str; default=''
            a comment for the card
        """
        self.model.rods.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._nid.append(nid)
        self._cid.append(cid)
        self._x.append(X)
        self._I.append(I)
        self._mass.append(mass)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a CONM2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        nid = integer(card, 2, 'nid')
        cid = integer_or_blank(card, 3, 'cid', 0)
        mass = double_or_blank(card, 4, 'mass', 0.)

        X = [
            double_or_blank(card, 5, 'x1', 0.0),
            double_or_blank(card, 6, 'x2', 0.0),
            double_or_blank(card, 7, 'x3', 0.0),
        ]

        I = [
            double_or_blank(card, 9, 'I11', 0.0),
            double_or_blank(card, 10, 'I21', 0.0),
            double_or_blank(card, 11, 'I22', 0.0),
            double_or_blank(card, 12, 'I31', 0.0),
            double_or_blank(card, 13, 'I32', 0.0),
            double_or_blank(card, 14, 'I33', 0.0),
        ]
        assert len(card) <= 15, f'len(CONM2 card) = {len(card):d}\ncard={card}'
        self.add(eid, nid, mass, cid=cid, X=X, I=I, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, nid, cid, mass, X, I in zip(self.eid, self.nid, self.cid, self.mass, self.x, self.I):
            Inew = []
            for i in I:
                if i == 0.:
                    Inew.append(None)
                else:
                    Inew.append(i)
            Xnew = []
            for x in X:
                if x == 0.:
                    Xnew.append(None)
                else:
                    Xnew.append(x)

            cid = set_blank_if_default(cid, 0)
            list_fields = (['CONM2', eid, nid, cid, mass] + Xnew +
                           [None] + Inew)

            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.cid = np.vstack([self.cid, self._cid])
                self.nid = np.hstack([self.nid, self._nid])
                self.x = np.hstack([self.x, self._x])
                self.I = np.hstack([self.I, self._I])
                self.mass = np.hstack([self.mass, self._mass])

                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.nid = np.array(self._nid, dtype='int32')
                self.cid = np.array(self._cid, dtype='int32')
                self.x = np.array(self._x, dtype='float64')
                self.I = np.array(self._I, dtype='float64')
                self.mass = np.array(self._mass, dtype='float64')
            assert len(self.eid) == len(np.unique(self.eid))
            self._eid = []
            self._nid = []
            self._cid = []
            self._x = []
            self._I = []
            self._mass = []
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
        self.is_current = False
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

    def add_card(self, card, comment=''):
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
