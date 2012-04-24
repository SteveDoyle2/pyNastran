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
#import sys

# pyNastran
from baseCard import Element

from elementsRigid  import *
from elementsSolid  import *
from    bars.elementsBars    import *
from  plates.elementsShell   import *
from springs.elementsSprings import *
from    mass.elementsMass    import *


class BushingElement(Element):
    def __init__(self,card,data):
        Element.__init__(self,card,data)

class CBUSH(BushingElement):
    type = 'CBUSH'
    def __init__(self,card=None,data=None):
        BushingElement.__init__(self,card,data)
        
        if card:
            fields = card.fields()
            #print "fields = ",fields
            self.eid = card.field(1)
            self.pid = card.field(2)
            self.ga = card.field(3)
            self.gb = card.field(4)
            x1G0    = card.field(5)
            if isinstance(x1G0,int):
                self.g0 = x1G0
                self.x = None
            elif isinstance(x1G0,float):
                self.g0  = None
                x1  = x1G0
                x2  = card.field(6)
                x3  = card.field(7)
                self.x   = [x1,x2,x3]
            else:
                #raise Exception('invalid CBUSH...x1/g0 = |%s|' %(x1G0))
                self.g0 = None
                self.x  = [None,None,None]
            ###
            ## Element coordinate system identification. A 0 means the basic
            ## coordinate system. If CID is blank, then the element coordinate system
            ## is determined from GO or Xi.
            self.cid  = card.field(8,0)
            ## Location of spring damper (0 <= s <= 1.0)
            self.s    = card.field(9,0.5)
            ## Coordinate system identification of spring-damper offset. See
            ## Remark 9. (Integer > -1; Default = -1, which means the offset point lies
            ## on the line between GA and GB
            self.ocid = card.field(10,-1)
            ## Components of spring-damper offset in the OCID coordinate system if OCID > 0.
            self.si   = card.fields(i=11,j=13)
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
        if isinstance(self.cid,int):
            return self.cid
        return self.cid.cid

    def crossReference(self,model):
        #self.nodes = model.Nodes(self.nodes)
        self.pid = model.Property(self.pid)
        self.cid = model.Coord(self.cid)
        
    def rawFields(self):
        if self.g0 is not None:
            x = [self.g0,None,None]
        else:
            x = self.x
        fields = ['CBUSH',self.eid,self.Pid(),self.ga,self.gb]+x+[self.Cid(),
                          self.s,self.ocid]+self.si
        return fields

    def reprFields(self):
        return self.rawFields()

class CFAST(Element):
    type = 'CFAST'
    def __init__(self,card=None,data=None):
        Element.__init__(self,card,data)
        self.eid  = card.field(1)
        self.pid  = card.field(2)
        self.Type = card.field(3)
        self.ida = card.field(4)
        self.idb = card.field(5)
        self.gs  = card.field(6)
        self.ga  = card.field(7)
        self.gb  = card.field(8)
        self.xs  = card.field(9)
        self.ys  = card.field(10)
        self.zs  = card.field(11)
        
        #if self.Type=='PROP': # PSHELL/PCOMP  ida & idb

    def crossReference(self,model):
        self.pid = model.Property(self.pid)
        #self.gs = model.Node(self.gs)
        #self.ga = model.Node(self.ga)
        #self.gb = model.Node(self.gb)
    
    def rawFields(self):
        fields = ['CFAST',self.eid,self.Pid(),self.Type,self.ida,self.idb,self.gs,self.ga,self.gb,
                          self.xs,self.ys,self.zs]
        return fields
    
    def reprFields(self):
        return self.rawFields()

class DamperElement(Element):
    def __init__(self,card,data):
        Element.__init__(self,card,data)

class CDAMP1(DamperElement):
    type = 'CDAMP1'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)

            nids = [card.field(3),card.field(5)]

            ## component number
            self.c1 = card.field(4)
            self.c2 = card.field(6)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[4]]
            self.c1  = data[3]
            self.c2  = data[5]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def isSameCard(self,elem):
        if self.type!=elem.type:  return False
        fields1 = [self.eid,self.Pid()]+self.nodeIDs()+[self.c1,self.c2]
        fields2 = [elem.eid,elem.Pid()]+elem.nodeIDs()+[elem.c1,elem.c2]
        return self.isSameFields(fields1,fields2)

    def B(self):
        return self.pid.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['CDAMP1',self.eid,self.Pid(),nodes[0],self.c1,nodes[1],self.c2]
        return fields

class CDAMP2(DamperElement):
    type = 'CDAMP2'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)

            ## Value of the scalar damper (Real)
            self.b   = card.field(2)

            nids = [card.field(3),card.field(5)]

            ## component number
            self.c1 = card.field(4)
            self.c2 = card.field(6)
        else:
            self.eid = data[0]
            self.b   = data[1]
            nids     = [data[2],data[4]]
            self.c1  = data[3]
            self.c2  = data[5]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def B(self):
        return self.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        
    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['CDAMP2',self.eid,self.b,nodes[0],self.c1,nodes[1],self.c2]
        return fields

class CDAMP3(DamperElement):
    type = 'CDAMP3'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = [card.field(3),card.field(4)]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[3]]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def B(self):
        return self.pid.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['CDAMP3',self.eid,self.pid,nodes[0],nodes[1]]
        return fields

class CDAMP4(DamperElement):
    type = 'CDAMP4'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)
            ## Value of the scalar damper (Real)
            self.b   = card.field(2)
            nids = [card.field(3),card.field(4)]
        else:
            self.eid = data[0]
            self.b   = data[1]
            nids     = [data[2],data[3]]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def B(self):
        return self.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        
    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['CDAMP4',self.eid,self.b,nodes[0],nodes[1]]
        return fields

class CDAMP5(DamperElement):
    type = 'CDAMP5'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)
            ## Property ID
            self.pid = card.field(2)
            nids = [card.field(3),card.field(4)]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[3]]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['CDAMP5',self.eid,self.Pid(),nodes[0],nodes[1]]
        return fields

class CGAP(Element):
    type = 'CGAP'
    def __init__(self,card=None,data=None):
        Element.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            self.ga = card.field(3)
            self.gb = card.field(4)
            x1G0    = card.field(5)
            if isinstance(x1G0,int):
                self.g0 = x1G0
                self.x = None
                self.cid = None
            elif isinstance(x1G0,float):
                self.g0  = None
                x1  = x1G0
                x2  = card.field(6)
                x3  = card.field(7)
                self.x   = [x1,x2,x3]
                self.cid = card.field(8,0)
            else:
                #raise Exception('invalid CGAP...x1/g0 = |%s|' %(x1G0))
                self.g0 = None
                self.x  = [None,None,None]
                self.cid = None
            ###
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.ga  = data[2]
            self.gb  = data[3]
            self.g0  = data[4]
            x1  = data[5]
            x2  = data[6]
            x3  = data[7]
            self.x   = [x1,x2,x3]
            self.cid = data[8]
        ###

    def crossReference(self,mesh):
        self.ga = mesh.Node(self.ga)
        self.gb = mesh.Node(self.gb)
        if self.g0:
            self.g0 = mesh.Node(self.g0)
            self.x  = self.g0.Position()
        self.pid = mesh.Property(self.pid)
        if self.cid:
            self.cid = mesh.Coord(self.cid)
        ###

    def Cid(self):
        if isinstance(self.cid,int) or self.cid is None:
            return self.cid
        return self.cid.cid

    def rawFields(self):
        if self.g0 is not None:
            x = [self.g0,None,None]
        else:
            x = self.x
        fields = ['CGAP',self.eid,self.Pid(),self.ga,self.gb]+x+[self.Cid()]
        return fields

class CRAC2D(Element):
    type = 'CRAC2D'
    def __init__(self,card=None,data=None):
        Element.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,21) # caps at 18
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==18

    def rawFields(self):
        fields = ['CRAC2D',self.eid,self.Pid()]+self.nodes
        return fields

class CRAC3D(Element):
    type = 'CRAC3D'
    def __init__(self,card=None,data=None):
        Element.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,67) # cap at +3 = 67
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==64

    def rawFields(self):
        fields = ['CRAC3D',self.eid,self.Pid()]+self.nodes
        return fields
