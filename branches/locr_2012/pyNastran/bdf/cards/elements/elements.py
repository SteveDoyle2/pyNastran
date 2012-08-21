# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from pyNastran.bdf.cards.baseCard import Element


class CFAST(Element):
    type = 'CFAST'

    def __init__(self, card=None, data=None):
        Element.__init__(self, card, data)
        self.eid = card.field(1)
        self.pid = card.field(2)
        self.Type = card.field(3)
        self.ida = card.field(4)
        self.idb = card.field(5)
        self.gs = card.field(6)
        self.ga = card.field(7)
        self.gb = card.field(8)
        self.xs = card.field(9)
        self.ys = card.field(10)
        self.zs = card.field(11)

        #if self.Type=='PROP': # PSHELL/PCOMP  ida & idb

    def cross_reference(self, model):
        self.pid = model.Property(self.pid)
        self.gs = model.Node(self.gs)
        self.ga = model.Node(self.ga)
        self.gb = model.Node(self.gb)

    def rawFields(self):
        fields = ['CFAST', self.eid, self.Pid(), self.Type, self.ida, self.idb, self.gs, self.ga, self.gb,
                  self.xs, self.ys, self.zs]
        return fields

    def reprFields(self):
        return self.rawFields()


class CGAP(Element):
    type = 'CGAP'

    def __init__(self, card=None, data=None):
        Element.__init__(self, card, data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            self.ga = card.field(3)
            self.gb = card.field(4)
            x1G0 = card.field(5)
            if isinstance(x1G0, int):
                self.g0 = x1G0
                self.x = None
                self.cid = None
            elif isinstance(x1G0, float):
                self.g0 = None
                x1 = x1G0
                x2 = card.field(6)
                x3 = card.field(7)
                self.x = [x1, x2, x3]
                self.cid = card.field(8, 0)
            else:
                #raise RuntimeError('invalid CGAP...x1/g0 = |%s|' %(x1G0))
                self.g0 = None
                self.x = [None, None, None]
                self.cid = None
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
        ###

    def cross_reference(self, model):
        self.ga = model.Node(self.ga)
        self.gb = model.Node(self.gb)
        if self.g0:
            self.g0 = model.Node(self.g0)
            self.x = self.g0.Position()
        self.pid = model.Property(self.pid)
        if self.cid:
            self.cid = model.Coord(self.cid)
        ###
    def Cid(self):
        if isinstance(self.cid, int) or self.cid is None:
            return self.cid
        return self.cid.cid

    def rawFields(self):
        if self.g0 is not None:
            x = [self.g0, None, None]
        else:
            x = self.x
        fields = ['CGAP', self.eid, self.Pid(), self.ga, self.gb] + x + [self.Cid()]
        return fields


class CrackElement(Element):
    def __init__(self, card, data):
        pass

    def cross_reference(self, model):
        self.nodes = model.Nodes(self.nodes)
        self.pid = model.Property(self.pid)

    def rawFields(self):
        fields = [self.type, self.eid, self.Pid()] + self.nodeIDs()
        return fields


class CRAC2D(CrackElement):
    type = 'CRAC2D'

    def __init__(self, card=None, data=None):
        CrackElement.__init__(self, card, data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3, 21)  # caps at 18
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 18


class CRAC3D(CrackElement):
    type = 'CRAC3D'

    def __init__(self, card=None, data=None):
        CrackElement.__init__(self, card, data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3, 67)  # cap at +3 = 67
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 64
