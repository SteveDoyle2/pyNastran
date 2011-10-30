#import sys
from numpy import array
#from numpy import cross,dot
#from numpy.linalg import norm

# my code
from AirQuick.general.generalMath import Area,Triangle_AreaCentroidNormal,Normal
from baseCard import Element

from elementsBars  import *
from elementsRigid import *
from elementsShell import *
from elementsSolid import *


class SpringElement(Element):
    def __init__(self,card):
        Element.__init__(self,card)
        self.eid = card.field(1)

class CELAS1(SpringElement):
    type = 'CELAS1'
    def __init__(self,card):
        SpringElement.__init__(self,card)
        self.pid = card.field(2)

        nids = [card.field(3),card.field(5)]
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==2

        ## property ID
        self.pid = card.field(2,self.eid)

        ## component number
        self.c1 = card.field(4)
        self.c2 = card.field(5)

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid(),self.nodes[0],self.c1,self.nodes[1],self.c2]
        return self.printCard(fields)

class CELAS2(SpringElement):
    type = 'CELAS2'
    def __init__(self,card):
        SpringElement.__init__(self,card)
        nids = [card.field(3),card.field(5)]
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==2

        ## stiffness of the scalar spring
        self.k   = card.field(2)
        
        ## component number
        self.c1 = card.field(4)
        self.c2 = card.field(5)
        
        ## damping coefficient
        self.ge = card.field(6)
        
        ## stress coefficient
        self.s  = card.field(7)

    def __repr__(self):
        fields = [self.type,self.eid,self.k,self.nodes[0],self.c1,self.nodes[1],self.c2,self.ge,self.s]
        return self.printCard(fields)

class CSHEAR(Element):
    type = 'CSHEAR'
    def __init__(self,card):
        Element.__init__(self,card)
        nids = card.fields(3,7)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==4

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodes
        return self.printCard(fields)

class CRAC2D(Element):
    type = 'CRAC2D'
    def __init__(self,card):
        Element.__init__(self,card)

        nids = card.fields(3,21) # caps at 18
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==18

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodes
        return self.printCard(fields)

class CRAC3D(Element):
    type = 'CRAC3D'
    def __init__(self,card):
        Element.__init__(self,card)

        nids = card.fields(3,67) # cap at +3 = 67
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==64

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodes
        return self.printCard(fields)
        
class CVISC(CROD):
    type = 'CVISC'
    def __init__(self,card):
        CROD.__init__(self,card)
    ###
    def __repr__(self):  # not done...
        fields = ['CVISC',self.eid]
###

class CONM2(Element): # v0.1 not done
    type = 'CONM2'
    # 'CONM2    501274  11064          132.274'
    def __init__(self,card):
        Element.__init__(self,card)
        #self.nids  = [ card[1] ]
        #del self.nids
        #self.pid = None
        self.eid  = card.field(1)
        self.gid  = card.field(2)
        self.cid  = card.field(3,0)
        self.mass = card.field(4)
        self.X    = array(card.fields(5,8,[0.,0.,0.]))
        self.I    = card.fields(9,15,[0.]*6)

    def crossReference(self,mesh):
        """
        @warning only supports cid=0
        """
        if self.cid==0:
            pass
        #else:
        #    raise Exception('CONM2 cid !=0 is not coded...')
        ###
    def __repr__(self):
        #I = []
        #for i in self.I:
        #    if i==0.:
        fields = [self.type,self.eid,self.gid,self.cid,self.mass]+list(self.X)+self.I
        return self.printCard(fields)

   
