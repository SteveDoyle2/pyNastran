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
#import sys
#from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.bdf.assign_type import (integer, integer_or_blank,
    integer_double_or_blank, double, double_or_blank, string_or_blank)


class BushElement(Element):
    def __init__(self, card, data):
        self.cid = None
        Element.__init__(self, card, data)

    def Cid(self):
        if self.cid is None:
            return None
        elif isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    #def Ga(self):
        #print dir(self)
        #if isinstance(self.ga, int):
            #return self.ga
        #return self.ga.nid

    #def Gb(self):
        #if isinstance(self.gb, int):
            #return self.gb
        #return self.gb.nid

    #def NodeIDs(self):
        #print self.nodeIDs()
        #return [self.Ga(), self.Gb()]


class CBUSH(BushElement):
    type = 'CBUSH'

    def __init__(self, card=None, data=None, comment=''):
        BushElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.ga = integer(card, 3, 'ga')
            self.gb = integer_or_blank(card, 4, 'gb')

            x1G0 = integer_double_or_blank(card, 5, 'x1_g0')
            if isinstance(x1G0, int):
                self.g0 = x1G0
                self.x = None
            elif isinstance(x1G0, float):
                self.g0 = None
                x1 = x1G0
                x2 = double(card, 6, 'x2')
                x3 = double(card, 7, 'x3')
                self.x = [x1, x2, x3]
            else:
                self.g0 = None
                self.x = [None, None, None]

            ## Element coordinate system identification. A 0 means the basic
            ## coordinate system. If CID is blank, then the element coordinate
            ## system is determined from GO or Xi.
            ## (default=blank=element-based)
            self.cid = integer_or_blank(card, 8, 'cid')
            ## Location of spring damper (0 <= s <= 1.0)
            self.s = double_or_blank(card, 9, 's', 0.5)
            ## Coordinate system identification of spring-damper offset. See
            ## Remark 9. (Integer > -1; Default = -1, which means the offset
            ## point lies on the line between GA and GB
            self.ocid = integer_or_blank(card, 10, 'ocid', -1)
            ## Components of spring-damper offset in the OCID coordinate system
            ## if OCID > 0.
            self.si = [double_or_blank(card, 11, 's1'),
                       double_or_blank(card, 12, 's2'),
                       double_or_blank(card, 13, 's3')]
            assert len(card) <= 13, 'len(CBUSH card) = %i' % len(card)
        else:
            self.eid = data[0]
            raise NotImplementedError('CBUSH data...')

    def Ga(self):
        if isinstance(self.ga, int):
            return self.ga
        return self.ga.nid

    def Gb(self):
        if isinstance(self.gb, int):
            return self.gb
        return self.gb.nid

    def OCid(self):
        if self.ocid is None:
            return None
        elif isinstance(self.ocid, int):
            return self.ocid
        return self.ocid.cid

    def Cid(self):
        if self.cid is None:
            return None
        elif isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def cross_reference(self, model):
        msg = ' which is required by CBUSH eid=%s' % self.eid
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        if self.cid is not None:
            self.cid = model.Coord(self.cid, msg=msg)

    def rawFields(self):
        if self.g0 is not None:
            x = [self.g0, None, None]
        else:
            x = self.x
        list_fields = (['CBUSH', self.eid, self.Pid(), self.Ga(), self.Gb()] + x +
                  [self.Cid(), self.s, self.ocid] + self.si)
        return list_fields

    def reprFields(self):
        if self.g0 is not None:
            x = [self.g0, None, None]
        else:
            x = self.x

        ocid = set_blank_if_default(self.OCid(), -1)
        s = set_blank_if_default(self.s, 0.5)
        list_fields = (['CBUSH', self.eid, self.Pid(), self.Ga(), self.Gb()] +
                  x + [self.Cid(), s, ocid] + self.si)
        return list_fields


class CBUSH1D(BushElement):
    type = 'CBUSH1D'

    def __init__(self, card=None, data=None, comment=''):
        BushElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.ga = integer(card, 3, 'ga')
            self.gb = integer_or_blank(card, 4, 'gb')
            self.cid = integer_or_blank(card, 5, 'cid')
            assert len(card) <= 5, 'len(CBUSH1D card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.ga = data[2]
            self.gb = data[3]

    def cross_reference(self, model):
        msg = ' which is required by CBUSH1D eid=%s' % self.eid
        self.ga = model.Node(self.ga, msg=msg)
        if self.gb:
            self.gb = model.Node(self.gb, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        if self.cid is not None:
            self.cid = model.Coord(self.cid)

    def Ga(self):
        if isinstance(self.ga, int):
            return self.ga
        #elif self.ga is None:
            #return None
        return self.ga.nid

    def Gb(self):
        if isinstance(self.gb, int):
            return self.gb
        elif self.gb is None:
            return None
        return self.gb.nid

    def rawFields(self):
        list_fields = ['CBUSH1D', self.eid, self.Pid(), self.Ga(), self.Gb(),
                  self.Cid()]
        return list_fields

    #def reprFields(self):
        #return self.rawFields()


class CBUSH2D(BushElement):
    """
    2-D Linear-Nonlinear Connection
    Defines the connectivity of a two-dimensional Linear-Nonlinear element.
    """
    type = 'CBUSH2D'

    def __init__(self, card=None, data=None, comment=''):
        BushElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid')
            self.ga = integer(card, 3, 'ga')
            self.gb = integer(card, 4, 'gb')

            self.cid = integer_or_blank(card, 5, 'cid', 0)
            self.plane = string_or_blank(card, 6, 'plane', 'XY')
            if self.plane not in ['XY', 'YZ', 'ZX']:
                msg = ("plane not in required list, plane=|%s|\n"
                       "expected planes = ['XY','YZ','ZX']" % self.plane)
                raise RuntimeError(msg)
            assert len(card) <= 6, 'len(CBUSH2D card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.ga = data[2]
            self.gb = data[3]

    def Ga(self):
        if isinstance(self.ga, int):
            return self.ga
        return self.ga.nid

    def Gb(self):
        if isinstance(self.gb, int):
            return self.gb
        return self.gb.nid

    def rawFields(self):
        list_fields = ['CBUSH2D', self.eid, self.Pid(), self.Ga(), self.Gb(),
                  self.Cid(), self.plane, self.sptid]
        return list_fields

    def cross_reference(self, model):
        msg = ' which is required by CBUSH2D eid=%s' % self.eid
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)
        #self.pid = model.Property(self.pid)
        if self.cid is not None:
            self.cid = model.Coord(self.cid, msg=msg)

    #def reprFields(self):
        #return self.rawFields()