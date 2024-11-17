from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class Bushes:
    """initializes the Bushes"""
    def __init__(self, model):
        self.model = model
        self.cbush = model.cbush
        #self.cbush1d = model.cbush1d
        #self.cbush2d = model.cbush2d
        self._eids = set()

    def add(self, eid):
        if eid not in self._eids:
            self._eids.add(eid)
        else:
            raise RuntimeError('eid=%s is duplicated' % eid)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.cbush):
            self.cbush.write_card(size, is_double, bdf_file)
        #if len(self.cbush1d):
            #self.cbush1d.write_card(size, is_double, bdf_file)
        #if len(self.cbush2d):
            #self.cbush2d.write_card(size, is_double, bdf_file)

    def make_current(self):
        self.cbush.make_current()
        #self.cbush1d.make_current()
        #self.cbush2d.make_current()

    def __len__(self):
        return len(self.cbush) #+ len(self.cbush1d) + len(self.cbush2d)

    def repr_indent(self, indent='  '):
        msg = '%s<Bushes> : nelements=%s\n' % (indent, len(self))
        msg += '%s  CBUSH:  %s\n' % (indent, len(self.cbush))
        #msg += '%s  CBUSH1D:  %s\n' % (indent, len(self.cbush1d))
        #msg += '%s  CBUSH2D:  %s\n' % (indent, len(self.cbush2d))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')


class BushElement:
    """base class for CBUSH, CBUSH1D, CBUSH2D"""
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


class CBUSHv(BushElement):
    """
    Generalized Spring-and-Damper Connection

    Defines a generalized spring-and-damper structural element that
    may be nonlinear or frequency dependent.
    +-------+-----+------+----+----+-------+----+----+-----+
    |   1   |  2  |  3   |  4 |  5 |   6   |  7 |  8 |  9  |
    +=======+=====+======+====+====+=======+====+====+=====+
    | CBUSH | EID | PID  | GA | GB | GO/X1 | X2 | X3 | CID |
    +-------+-----+------+----+----+-------+----+----+-----+
    |       |  S  | OCID | S1 | S2 |  S3   |    |    |     |
    +-------+-----+------+----+----+-------+----+----+-----+
    """
    card_name = 'CBUSH'

    def __init__(self, model):
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')
        self.cid = np.array([], dtype='int32')
        self.ocid = np.array([], dtype='int32')
        self.s = np.array([], dtype='float64')
        self.si = np.array([], dtype='float64')
        #self.nsm = np.array([], dtype='float64')

        self._eid = []
        self._pid = []
        self._nids = []
        self._cid = []
        self._ocid = []
        self._s = []
        self._si = []
        #self._nsm = []
        self.comment = defaultdict(str)

    def add(self, eid, pid, nids, x, g0, cid=None, s=0.5, ocid=-1, si=None, comment=''):
        """
        Creates a CBUSH card

        Parameters
        ----------
        eid : int
            Element id
        pid : int
            Property id (PBUSH)
        nids : list[int, int]
            node ids; connected grid points at ends A and B
            The nodes may be coincident, but then cid is required.
        x : list[float, float, float]; None
            List : the directional vector used to define the stiffnesses
                   or damping from the PBUSH card
            None : use g0
        g0 : int/None
            int : the directional vector used to define the stiffnesses
                  or damping from the PBUSH card
            None : use x
        cid : int; default=None
            Element coordinate system identification. A 0 means the basic
            coordinate system. If CID is blank, then the element coordinate
            system is determined from GO or Xi.
        s: float; default=0.5
            Location of spring damper (0 <= s <= 1.0)
        ocid : int; default=-1
            Coordinate system identification of spring-damper offset.
            (Integer > -1; Default = -1, which means the offset
            point lies on the line between GA and GB)
        si : list[float, float, float]; default=None
            Components of spring-damper offset in the OCID coordinate system
            if OCID > 0.
            None : [None, None, None]
        comment : str; default=''
            a comment for the card
        """
        if cid is None:
            cid = -1
        self.model.bushes.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        self._cid.append(cid)
        self._ocid.append(ocid)
        self._s.append(s)
        self._si.append(si)
        #self._nsm.append(nsm)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a CBUSH card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer_or_blank(card, 4, 'gb')

        #: Element coordinate system identification. A 0 means the basic
        #: coordinate system. If CID is blank, then the element coordinate
        #: system is determined from GO or Xi.
        #: (default=blank=element-based)
        cid = integer_or_blank(card, 8, 'cid')

        x1_g0 = integer_double_or_blank(card, 5, 'x1_g0')
        if isinstance(x1_g0, integer_types):
            g0 = x1_g0
            x = None
        elif isinstance(x1_g0, float):
            g0 = None
            x1 = x1_g0
            x2 = double_or_blank(card, 6, 'x2', 0.0)
            x3 = double_or_blank(card, 7, 'x3', 0.0)
            x = [x1, x2, x3]
            if not isinstance(cid, integer_types):
                assert max(x) != min(x), 'x=%s' % x
        else:
            g0 = None
            x = [None, None, None]

        #: Location of spring damper (0 <= s <= 1.0)
        s = double_or_blank(card, 9, 's', 0.5)
        #: Coordinate system identification of spring-damper offset. See
        #: Remark 9. (Integer > -1; Default = -1, which means the offset
        #: point lies on the line between GA and GB
        ocid = integer_or_blank(card, 10, 'ocid', -1)
        #: Components of spring-damper offset in the OCID coordinate system
        #: if OCID > 0.
        si = [double_or_blank(card, 11, 's1'),
              double_or_blank(card, 12, 's2'),
              double_or_blank(card, 13, 's3')]
        assert len(card) <= 14, 'len(CBUSH card) = %i\ncard=%s' % (len(card), card)
        #return CBUSH(eid, pid, [ga, gb], x, g0, cid=cid, s=s, ocid=ocid, si=si, comment=comment)
        return self.add(eid, pid, [ga, gb], x, g0, cid=cid, s=s, ocid=ocid, si=si, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        assert len(self.eid) == len(self.pid)
        assert len(self.eid) == len(self.cid)
        assert len(self.eid) == len(self.ocid)
        assert len(self.eid) == len(self.s)
        assert len(self.eid) == len(self.si)

        for eid, pid, nodes, cid, ocid, s, si in zip(self.eid, self.pid, self.nids, self.cid, self.ocid, self.s, self.si):
            ga, gb = nodes
            #if self.g0 is not None:
                #x = [self.g0, None, None]
            #else:
                #x = self.x

            x = [None, None, None]
            ocid = set_blank_if_default(ocid, -1)
            s = set_blank_if_default(s, 0.5)
            list_fields = (['CBUSH', eid, pid, ga, gb] +
                           x + [cid, s, ocid] + si.tolist())

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
                self.cid = np.hstack([self.cid, self._cid])
                self.ocid = np.hstack([self.ocid, self._ocid])
                self.s = np.hstack([self.s, self._s])
                self.si = np.hstack([self.si, self._si])
                #self.nsm = np.hstack([self.nsm, self._nsm])

                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.pid = np.array(self._pid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
                self.cid = np.array(self._cid, dtype='int32')
                self.ocid = np.array(self._ocid, dtype='int32')
                self.s = np.array(self._s, dtype='float64')
                self.si = np.array(self._si, dtype='float64')
                #self.nsm = np.array(self._nsm, dtype='float64')
            assert len(self.eid) == len(np.unique(self.eid))
            #print(self.nid)
            self._eid = []
            self._pid = []
            self._nids = []
            self._cid = []
            self._ocid = []
            self._s = []
            self._si = []
            #self._A = []
            #self._j = []
            #self._c = []
            #self._nsm = []
            self.is_current = True
