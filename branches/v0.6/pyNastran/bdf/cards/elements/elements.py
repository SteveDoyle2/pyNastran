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
All ungrouped elements are defined in this file.  This includes:

 * CFAST
 * CGAP
 * CRAC2D
 * CRAC3D

All ungrouped elements are Element objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from pyNastran.bdf.cards.baseCard import Element
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    integer_double_or_blank, double_or_blank, string)  # double
 
class CFAST(Element):
    type = 'CFAST'

    def __init__(self, card=None, data=None, comment=''):
        Element.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.Type = string(card, 3, 'Type')
            self.ida = integer(card, 4, 'ida')
            self.idb = integer(card, 5, 'idb')
            self.gs = integer_or_blank(card, 6, 'gs')
            self.ga = integer_or_blank(card, 7, 'ga')
            self.gb = integer_or_blank(card, 8, 'gb')
            self.xs = double_or_blank(card, 9, 'xs')
            self.ys = double_or_blank(card, 10, 'ys')
            self.zs = double_or_blank(card, 11, 'zs')
            assert len(card) <= 12, 'len(CFAST card) = %i' % len(card)
        else:
            raise NotImplementedError(data)
        #if self.Type=='PROP': # PSHELL/PCOMP  ida & idb

    def cross_reference(self, model):
        msg = ' which is required by CFAST eid=%s' % self.eid
        self.pid = model.Property(self.pid, msg=msg)
        self.gs = model.Node(self.gs, msg=msg)
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)

    def rawFields(self):
        list_fields = ['CFAST', self.eid, self.Pid(), self.Type, self.ida, self.idb,
                  self.gs, self.ga, self.gb, self.xs, self.ys, self.zs]
        return list_fields

    def reprFields(self):
        return self.rawFields()


class CGAP(Element):
    type = 'CGAP'

    def __init__(self, card=None, data=None, comment=''):
        Element.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.ga = integer_or_blank(card, 3, 'ga')
            self.gb = integer_or_blank(card, 4, 'gb')
            x1G0 = integer_double_or_blank(card, 5, 'x1_g0')
            if isinstance(x1G0, int):
                self.g0 = x1G0
                self.x = None
                self.cid = None
            elif isinstance(x1G0, float):
                self.g0 = None
                x1 = x1G0
                x2 = double_or_blank(card, 6, 'x2', 0.0)
                x3 = double_or_blank(card, 7, 'x3', 0.0)
                self.x = [x1, x2, x3]
                self.cid = integer_or_blank(card, 8, 'cid', 0)
            else:
                #raise RuntimeError('invalid CGAP...x1/g0 = |%s|' %(x1G0))
                self.g0 = None
                self.x = [None, None, None]
                self.cid = None
            assert len(card) <= 9, 'len(CGAP card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.ga = data[2]
            self.gb = data[3]
            self.g0 = data[4]
            x1 = data[5]
            x2 = data[6]
            x3 = data[7]
            self.x = [x1, x2, x3]
            self.cid = data[8]

    def cross_reference(self, model):
        msg = ' which is required by CGAP eid=%s' % self.eid
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)
        if self.g0:
            self.g0 = model.Node(self.g0, msg=msg)
            self.x = self.g0.Position()
        self.pid = model.Property(self.pid, msg=msg)
        if self.cid:
            self.cid = model.Coord(self.cid, msg=msg)

    def Cid(self):
        if isinstance(self.cid, int) or self.cid is None:
            return self.cid
        return self.cid.cid

    def Ga(self):
        if isinstance(self.ga, int):
            return self.ga
        return self.ga.nid

    def Gb(self):
        if isinstance(self.gb, int):
            return self.gb
        return self.gb.nid

    def rawFields(self):
        if self.g0 is not None:
            x = [self.g0, None, None]
        else:
            x = self.x
        list_fields = (['CGAP', self.eid, self.Pid(), self.Ga(), self.Gb()] + x +
                  [self.Cid()])
        return list_fields


class CrackElement(Element):
    type = 'Crack'

    def __init__(self, card, data):
        pass

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self. type, self.eid)
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)


class CRAC2D(CrackElement):
    type = 'CRAC2D'

    def __init__(self, card=None, data=None, comment=''):
        CrackElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer(card, 2, 'pid')
            nids = [integer(card, 3, 'n1'), integer(card, 4, 'n2'),
                    integer(card, 5, 'n3'), integer(card, 6, 'n4'),
                    integer(card, 7, 'n5'), integer(card, 8, 'n6'),
                    integer(card, 9, 'n7'), integer(card, 10, 'n8'),
                    integer(card, 11, 'n9'), integer(card, 12, 'n10'),
                    integer_or_blank(card, 13, 'n11'),
                    integer_or_blank(card, 14, 'n12'),
                    integer_or_blank(card, 15, 'n13'),
                    integer_or_blank(card, 16, 'n14'),
                    integer_or_blank(card, 17, 'n15'),
                    integer_or_blank(card, 18, 'n16'),
                    integer_or_blank(card, 19, 'n17'),
                    integer_or_blank(card, 20, 'n18')]
            assert len(card) <= 21, 'len(CRAC2D card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 18

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def rawFields(self):
        list_fields = ['CRAC2D', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields


class CRAC3D(CrackElement):
    type = 'CRAC3D'

    def __init__(self, card=None, data=None, comment=''):
        CrackElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer(card, 2, 'pid')
            # required 1-10, 19-28
            # optional 11-18, 29-36, 37-64
            # all/none 37-46 
            nids = card.fields(3, 67)  # cap at +3 = 67
            assert len(card) <= 67, 'len(CRAC3D card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 64

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def rawFields(self):
        list_fields = ['CRAC3D', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields