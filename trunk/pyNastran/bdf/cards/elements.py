#import sys

# pyNastran
from pyNastran.bdf.cards.baseCard import Element

from pyNastran.bdf.cards.elementsRigid import *
from pyNastran.bdf.cards.elementsSolid import *
from pyNastran.bdf.cards.bars.elementsBars import *
from pyNastran.bdf.cards.plates.elementsShell import *
from pyNastran.bdf.cards.springs.elementsSprings import *
from pyNastran.bdf.cards.mass.elementsMass import *
from pyNastran.bdf.cards.bush.elementsBush import *

#------------------------------------------------------------------------------
class CFAST(Element):  ## @todo xref
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

#------------------------------------------------------------------------------
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

            nids = [card.field(3,0),card.field(5,0)]

            ## component number
            self.c1 = card.field(4,0)
            self.c2 = card.field(6,0)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[4]]
            self.c1  = data[3]
            self.c2  = data[5]
        ###
        assert self.c1 in [0,1,2,3,4,5,6],'c1=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' %(str(self),self.c1)
        assert self.c2 in [0,1,2,3,4,5,6],'c2=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' %(str(self),self.c2)
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
        nodes = self.nodeIDs(allowEmptyNodes=True)
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

            nids = [card.field(3,0),card.field(5,0)]

            ## component number
            self.c1 = card.field(4,0)
            self.c2 = card.field(6,0)
        else:
            self.eid = data[0]
            self.b   = data[1]
            nids     = [data[2],data[4]]
            self.c1  = data[3]
            self.c2  = data[5]
        ###
        assert self.c1 in [0,1,2,3,4,5,6],'c1=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' %(str(self),self.c1)
        assert self.c2 in [0,1,2,3,4,5,6],'c2=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' %(str(self),self.c2)
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def B(self):
        return self.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes,allowEmptyNodes=True)
        
    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        fields = ['CDAMP2',self.eid,self.b,nodes[0],self.c1,nodes[1],self.c2]
        return fields

class CDAMP3(DamperElement):
    type = 'CDAMP3'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = [card.field(3,0),card.field(4,0)]
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
        self.nodes = model.Nodes(self.nodes,allowEmptyNodes=True)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
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
            nids = [card.field(3,0),card.field(4,0)]
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
        self.nodes = model.Nodes(self.nodes,allowEmptyNodes=True)
        
    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
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
            nids = [card.field(3,0),card.field(4,0)]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[3]]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes,allowEmptyNodes=True)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        fields = ['CDAMP5',self.eid,self.Pid(),nodes[0],nodes[1]]
        return fields

#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
class CrackElement(Element):
    def __init__(self,card,data):
        pass

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.pid = self.Property(self.pid)

    def rawFields(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()
        return fields

class CRAC2D(CrackElement):
    type = 'CRAC2D'
    def __init__(self,card=None,data=None):
        CrackElement.__init__(self,card,data)
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

class CRAC3D(CrackElement):
    type = 'CRAC3D'
    def __init__(self,card=None,data=None):
        CrackElement.__init__(self,card,data)
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
