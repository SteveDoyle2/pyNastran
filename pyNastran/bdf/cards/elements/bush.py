# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All bush elements are defined in this file.  This includes:

 * CBUSH
 * CBUSH1D
 * CBUSH2D

All bush elements are BushElement and Element objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Element
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_double_or_blank, double_or_blank,
    string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8


class BushElement(Element):
    def __init__(self):
        self.cid = None
        Element.__init__(self)

    def Cid(self):
        if self.cid is None:
            return None
        elif isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def Mass(self):
        return 0.

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        return [tuple(sorted(self.node_ids))]

class CBUSH(BushElement):
    type = 'CBUSH'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', 8:'cid', 9:'s', 10:'ocid'
    }

    def _update_field_helper(self, n, value):
        if n == 11:
            self.si[0] = value
        elif n == 12:
            self.si[1] = value
        elif n == 13:
            self.si[2] = value
        else:
            if self.g0 is not None:
                if n == 5:
                    self.g0 = value
                else:
                    raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))
            else:
                if n == 5:
                    self.x[0] = value
                elif n == 6:
                    self.x[1] = value
                elif n == 7:
                    self.x[2] = value
                else:
                    raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, pid, ga, gb, x, g0, cid, s, ocid, si, comment=''):
        BushElement.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.ga = ga
        self.gb = gb
        self.x = x
        self.g0 = g0
        self.cid = cid
        self.s = s
        self.ocid = ocid
        self.si = si

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer_or_blank(card, 4, 'gb')

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
            assert max(x) != min(x), 'x=%s' % x
        else:
            g0 = None
            x = [None, None, None]

        #: Element coordinate system identification. A 0 means the basic
        #: coordinate system. If CID is blank, then the element coordinate
        #: system is determined from GO or Xi.
        #: (default=blank=element-based)
        cid = integer_or_blank(card, 8, 'cid')
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
        return CBUSH(eid, pid, ga, gb, x, g0, cid, s, ocid, si, comment=comment)

    @classmethod
    def add_op2_data(cls, data, f, comment=''):
        ((eid, pid, ga, gb, cid, s, ocid, si), x, g0) = data
        return CBUSH(eid, pid, ga, gb, x, g0, cid, s, ocid, si, comment=comment)

    def Eid(self):
        return self.eid

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def nodes(self):
        return [self.ga, self.gb]

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def _verify(self, xref=False):
        ga = self.Ga()
        gb = self.Gb()
        cid = self.Cid()
        ocid = self.OCid()
        pid = self.Pid()
        #si = self.si
        assert isinstance(ga, integer_types), 'ga=%r' % ga
        assert isinstance(gb, integer_types) or gb is None, 'gb=%r' % gb
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(cid, integer_types) or cid is None, 'cid=%r' % cid
        assert isinstance(ocid, integer_types), 'ocid=%r' % ocid

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, integer_types) or self.gb is None:
            return self.gb
        return self.gb_ref.nid

    def OCid(self):
        if self.ocid is None:
            return None
        elif isinstance(self.ocid, integer_types):
            return self.ocid
        return self.ocid_ref.cid

    def Cid(self):
        if self.cid is None:
            return None
        elif isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CBUSH eid=%s' % self.eid
        self.ga = model.Node(self.ga, msg=msg)
        self.ga_ref = self.ga
        if self.gb is not None:
            self.gb = model.Node(self.gb, msg=msg)
            self.gb_ref = self.gb
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid
        if self.cid is not None:
            self.cid = model.Coord(self.cid, msg=msg)
            self.cid_ref = self.cid

    def uncross_reference(self):
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.pid = self.Pid()
        self.cid = self.Cid()
        if self.cid is not None:
            del self.cid_ref
        del self.ga_ref, self.gb_ref, self.pid_ref

    def raw_fields(self):
        if self.g0 is not None:
            x = [self.g0, None, None]
        else:
            x = self.x
        list_fields = (['CBUSH', self.eid, self.Pid(), self.Ga(), self.Gb()] + x +
                       [self.Cid(), self.s, self.ocid] + self.si)
        return list_fields

    def repr_fields(self):
        if self.g0 is not None:
            x = [self.g0, None, None]
        else:
            x = self.x

        ocid = set_blank_if_default(self.OCid(), -1)
        s = set_blank_if_default(self.s, 0.5)
        list_fields = (['CBUSH', self.eid, self.Pid(), self.Ga(), self.Gb()] +
                       x + [self.Cid(), s, ocid] + self.si)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CBUSH1D(BushElement):
    type = 'CBUSH1D'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', 5:'cid',
    }

    def __init__(self, eid, pid, ga, gb, cid, comment=''):
        if comment:
            self._comment = comment
        BushElement.__init__(self)
        self.eid = eid
        self.pid = pid
        self.ga = ga
        self.gb = gb
        self.cid = cid

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer_or_blank(card, 4, 'gb')
        cid = integer_or_blank(card, 5, 'cid')
        assert len(card) <= 6, 'len(CBUSH1D card) = %i\ncard=%s' % (len(card), card)
        return CBUSH1D(eid, pid, ga, gb, cid, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #eid = data[0]
        #pid = data[1]
        #ga = data[2]
        #gb = data[3]
        #raise NotImplementedError(data)
        #return CBUSH1D(eid, pid, ga, gb, cid, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CBUSH1D eid=%s' % self.eid
        self.ga = model.Node(self.ga, msg=msg)
        self.ga_ref = self.ga
        if self.gb:
            self.gb = model.Node(self.gb, msg=msg)
            self.gb_ref = self.gb
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid
        if self.cid is not None:
            self.cid = model.Coord(self.cid)
            self.cid_ref = self.cid

    def uncross_reference(self):
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.cid = self.Cid()
        self.pid = self.Pid()
        del self.ga_ref, self.gb_ref, self.cid_ref, self.pid_ref

    def _verify(self, xref=False):
        ga = self.Ga()
        gb = self.Gb()
        cid = self.Cid()
        pid = self.Pid()
        assert isinstance(ga, integer_types), 'ga=%r' % ga
        assert isinstance(gb, integer_types) or gb is None, 'gb=%r' % gb
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(cid, integer_types) or cid is None, 'cid=%r' % cid

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        #elif self.ga is None:
            #return None
        return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, integer_types):
            return self.gb
        elif self.gb is None:
            return None
        return self.gb_ref.nid

    @property
    def nodes(self):
        return [self.ga, self.gb]

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def raw_fields(self):
        list_fields = ['CBUSH1D', self.eid, self.Pid(), self.Ga(), self.Gb(),
                       self.Cid()]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CBUSH2D(BushElement):
    """
    2-D Linear-Nonlinear Connection
    Defines the connectivity of a two-dimensional Linear-Nonlinear element.
    """
    type = 'CBUSH2D'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', 5:'cid', 6:'plane', 7:'sptid',
    }

    def __init__(self, eid, pid, ga, gb, cid, plane, sptid, comment=''):
        BushElement.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.ga = ga
        self.gb = gb
        self.cid = cid
        self.plane = plane
        self.sptid = sptid
        if self.plane not in ['XY', 'YZ', 'ZX']:
            msg = ("plane not in required list, plane=|%s|\n"
                   "expected planes = ['XY','YZ','ZX']" % self.plane)
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid')
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        cid = integer_or_blank(card, 5, 'cid', 0)
        plane = string_or_blank(card, 6, 'plane', 'XY')
        sptid = integer_or_blank(card, 7, 'sptid')
        assert len(card) <= 8, 'len(CBUSH2D card) = %i\ncard=%s' % (len(card), card)
        return CBUSH2D(eid, pid, ga, gb, cid, plane, sptid, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #eid = data[0]
        #pid = data[1]
        #ga = data[2]
        #gb = data[3]
        #raise NotImplementedError(data)
        #return CBUSH2D(eid, pid, ga, gb, cid, plane, sptid, comment=comment)

    def _verify(self, xref=False):
        ga = self.Ga()
        gb = self.Gb()
        cid = self.Cid()
        pid = self.Pid()
        plane = self.plane
        assert isinstance(ga, integer_types), 'ga=%r' % ga
        assert isinstance(gb, integer_types), 'gb=%r' % gb
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(cid, integer_types), 'cid=%r' % cid
        assert self.plane in ['XY', 'YZ', 'ZX'], 'plane=%r' % plane

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, integer_types):
            return self.gb
        return self.gb_ref.nid

    @property
    def nodes(self):
        return [self.ga, self.gb]

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CBUSH2D eid=%s' % self.eid
        self.ga = model.Node(self.ga, msg=msg)
        self.ga_ref = self.ga
        self.gb = model.Node(self.gb, msg=msg)
        self.gb_ref = self.gb
        self.pid = model.Property(self.pid)
        self.pid_ref = self.pid
        if self.cid is not None:
            self.cid = model.Coord(self.cid, msg=msg)
            self.cid_ref = self.cid
        if self.sptid is not None:
            pass

    def uncross_reference(self):
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.cid = self.Cid()
        self.pid = self.Pid()
        del self.ga_ref, self.gb_ref, self.cid_ref, self.pid_ref

    def raw_fields(self):
        list_fields = ['CBUSH2D', self.eid, self.Pid(), self.Ga(), self.Gb(),
                       self.Cid(), self.plane, self.sptid]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)
