## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
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
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    integer_double_or_blank, double_or_blank, string_or_blank) # double


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
                x2 = double_or_blank(card, 6, 'x2', 0.0)
                x3 = double_or_blank(card, 7, 'x3', 0.0)
                self.x = [x1, x2, x3]
                assert max(self.x) != min(self.x), 'x=%s' % self.x
            else:
                self.g0 = None
                self.x = [None, None, None]

            #: Element coordinate system identification. A 0 means the basic
            #: coordinate system. If CID is blank, then the element coordinate
            #: system is determined from GO or Xi.
            #: (default=blank=element-based)
            self.cid = integer_or_blank(card, 8, 'cid')
            #: Location of spring damper (0 <= s <= 1.0)
            self.s = double_or_blank(card, 9, 's', 0.5)
            #: Coordinate system identification of spring-damper offset. See
            #: Remark 9. (Integer > -1; Default = -1, which means the offset
            #: point lies on the line between GA and GB
            self.ocid = integer_or_blank(card, 10, 'ocid', -1)
            #: Components of spring-damper offset in the OCID coordinate system
            #: if OCID > 0.
            self.si = [double_or_blank(card, 11, 's1'),
                       double_or_blank(card, 12, 's2'),
                       double_or_blank(card, 13, 's3')]
            assert len(card) <= 14, 'len(CBUSH card) = %i' % len(card)
        else:
            self.eid = data[0]
            raise NotImplementedError('CBUSH data...')

    def nodeIDs(self):
        return [self.Ga(), self.Gb()]

    def _verify(self, xref=False):
        ga = self.Ga()
        gb = self.Gb()
        cid = self.Cid()
        ocid = self.OCid()
        pid = self.Pid()
        #si = self.si
        assert isinstance(ga, int), 'ga=%r' % ga
        assert isinstance(gb, int), 'gb=%r' % gb
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(cid, int) or cid is None, 'cid=%r' % cid
        assert isinstance(ocid, int), 'ocid=%r' % ocid

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
            assert len(card) <= 6, 'len(CBUSH1D card) = %i' % len(card)
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

    def _verify(self, xref=False):
        ga = self.Ga()
        gb = self.Gb()
        cid = self.Cid()
        pid = self.Pid()
        assert isinstance(ga, int), 'ga=%r' % ga
        assert isinstance(gb, int), 'gb=%r' % gb
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(cid, int) or cid is None, 'cid=%r' % cid

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

    def nodeIDs(self):
        return [self.Ga(), self.Gb()]

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
            assert len(card) <= 7, 'len(CBUSH2D card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.ga = data[2]
            self.gb = data[3]

    def _verify(self, xref=False):
        ga = self.Ga()
        gb = self.Gb()
        cid = self.Cid()
        pid = self.Pid()
        plane = self.plane
        assert isinstance(ga, int), 'ga=%r' % ga
        assert isinstance(gb, int), 'gb=%r' % gb
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(cid, int), 'cid=%r' % cid
        assert self.plane in ['XY', 'YZ', 'ZX'], 'plane=%r' % plane

    def Ga(self):
        if isinstance(self.ga, int):
            return self.ga
        return self.ga.nid

    def Gb(self):
        if isinstance(self.gb, int):
            return self.gb
        return self.gb.nid

    def nodeIDs(self):
        return [self.Ga(), self.Gb()]

    def cross_reference(self, model):
        msg = ' which is required by CBUSH2D eid=%s' % self.eid
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)
        #self.pid = model.Property(self.pid)
        if self.cid is not None:
            self.cid = model.Coord(self.cid, msg=msg)

    def rawFields(self):
        list_fields = ['CBUSH2D', self.eid, self.Pid(), self.Ga(), self.Gb(),
                  self.Cid(), self.plane, self.sptid]
        return list_fields

    #def reprFields(self):
        #return self.rawFields()