# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
#from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element


class BushElement(Element):
    def __init__(self, card, data):
        self.cid = None
        Element.__init__(self, card, data)

    def Cid(self):
        if isinstance(self.cid, int):
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

    def __init__(self, card=None, data=None):
        BushElement.__init__(self, card, data)

        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            self.ga = card.field(3)
            self.gb = card.field(4)
            x1G0 = card.field(5)
            if isinstance(x1G0, int):
                self.g0 = x1G0
                self.x = None
            elif isinstance(x1G0, float):
                self.g0 = None
                x1 = x1G0
                x2 = card.field(6)
                x3 = card.field(7)
                self.x = [x1, x2, x3]
            else:
                #raise RuntimeError('invalid CBUSH...x1/g0 = |%s|' %(x1G0))
                self.g0 = None
                self.x = [None, None, None]
            ###
            ## Element coordinate system identification. A 0 means the basic
            ## coordinate system. If CID is blank, then the element coordinate
            ## system is determined from GO or Xi.
            self.cid = card.field(8, 0)
            ## Location of spring damper (0 <= s <= 1.0)
            self.s = card.field(9, 0.5)
            ## Coordinate system identification of spring-damper offset. See
            ## Remark 9. (Integer > -1; Default = -1, which means the offset
            ## point lies on the line between GA and GB
            self.ocid = card.field(10, -1)
            ## Components of spring-damper offset in the OCID coordinate system
            ## if OCID > 0.
            self.si = card.fields(i=11, j=13)
        else:
            self.eid = data[0]
            raise NotImplementedError('CBUSH data...')
        #self.prepareNodeIDs(nids,allowEmptyNodes=True)
        #assert len(self.nodes)==2

    #def OCid(self):
        #if isinstance(self.ocid,int):
            #return self.ocid
        #return self.ocid.cid

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def cross_reference(self, model):
        self.nodes = model.Nodes(self.nodes)
        #self.pid = model.Property(self.pid)
        self.cid = model.Coord(self.cid)

    def rawFields(self):
        if self.g0 is not None:
            x = [self.g0, None, None]
        else:
            x = self.x
        fields = ['CBUSH', self.eid, self.Pid(), self.ga, self.gb] + x + [self.Cid(),
                                                                          self.s, self.ocid] + self.si
        return fields

    def reprFields(self):
        if self.g0 is not None:
            x = [self.g0, None, None]
        else:
            x = self.x

        cid = set_blank_if_default(self.Cid(), 0)
        fields = ['CBUSH', self.eid, self.Pid(), self.ga, self.gb] + x + [cid,
                                                                          self.s, self.ocid] + self.si
        return fields


class CBUSH1D(BushElement):
    type = 'CBUSH1D'

    def __init__(self, card=None, data=None):
        BushElement.__init__(self, card, data)
        if card:
            self.eid = int(card.field(1))
            self.pid = int(card.field(2, self.eid))
            nids = card.fields(3, 5)
            self.cid = card.field(5, 0)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        ###
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 2

    def cross_reference(self, model):
        self.nodes = model.Nodes(self.nodes)
        self.pid = model.Property(self.pid)
        self.cid = model.Coord(self.cid)

    def rawFields(self):
        nodeIDs = self.nodeIDs()
        fields = ['CBUSH1D', self.eid, self.Pid(), nodeIDs[0],
                  nodeIDs[1], self.Cid()]
        return fields

    def reprFields(self):
        nodeIDs = self.nodeIDs()
        cid = set_blank_if_default(self.Cid(), 0)
        fields = ['CBUSH1D', self.eid, self.Pid(), nodeIDs[0], nodeIDs[1], cid]
        return fields


class CBUSH2D(BushElement):
    """
    2-D Linear-Nonlinear Connection
    Defines the connectivity of a two-dimensional Linear-Nonlinear element.
    """
    type = 'CBUSH2D'

    def __init__(self, card=None, data=None):
        BushElement.__init__(self, card, data)
        if card:
            self.eid = int(card.field(1))
            self.pid = int(card.field(2))
            nids = card.fields(3, 5)
            self.cid = int(card.field(5))
            self.plane = card.field(6, 'XY')
            self.sptid = card.field(7)
            if self.plane not in ['XY', 'YZ', 'ZX']:
                msg = 'plane not in required list, plane=|%s|\n' % (self.plane)
                msg += "expected planes = ['XY','YZ','ZX']"
                raise RuntimeError(msg)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        ###
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 2

    def rawFields(self):
        nodeIDs = self.nodeIDs()
        fields = ['CBUSH1D', self.eid, self.Pid(), nodeIDs[0],
                  nodeIDs[0], self.Cid(), self.plane, self.sptid]
        return fields

    def cross_reference(self, model):
        self.nodes = model.Nodes(self.nodes)
        #self.pid = model.Property(self.pid)
        self.cid = model.Coord(self.cid)

    #def reprFields(self):
        #return self.rawFields()
###
