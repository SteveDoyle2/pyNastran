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

    def length_noXref(self,n1=None,n2=None):
        """
        Returns the length of a bar/rod/beam element
        \f[ \large \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  } \f]
        @param self the object pointer
        @param n1 a Node object (default=None)
        @param n2 a Node object (default=None)
        @note
            if n1 AND n2 are both none (the default), then the model must
            be cross-referenced already
        """
        #print self.type
        L = norm(n1.Position()-n2.Position())
        return L

    def length(self):
        """
        Returns the length of a bar/rod/beam element
        \f[ \large \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  } \f]
        @param self the object pointer
        @note
            the model must be cross-referenced already
        """
        #print self.type
        return self.length_noXref(self.nodes[1],self.nodes[0])

    def mass(self):
        return 0.0

class PointElement(Element):
    def __init__(self,card):
        Element.__init__(self,card)

class CELAS1(SpringElement):
    type = 'CELAS1'
    def __init__(self,card):
        SpringElement.__init__(self,card)
        self.eid = card.field(2)

        nids = [card.field(3),card.field(5)]
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==2

        ## property ID
        self.pid = card.field(2,self.eid)

        ## component number
        self.c1 = card.field(4)
        self.c2 = card.field(5)

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.pid   = model.Property(self.pid)
        
    def __repr__(self):
        nodes = self.nodeIDs()
        fields = [self.type,self.eid,self.Pid(),nodes[0],self.c1,nodes[1],self.c2]
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

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        
    def __repr__(self):
        nodes = self.nodeIDs()
        fields = [self.type,self.eid,self.k,nodes[0],self.c1,nodes[1],self.c2,self.ge,self.s]
        return self.printCard(fields)

class CSHEAR(Element): # not integrated
    type = 'CSHEAR'
    def __init__(self,card):
        Element.__init__(self,card)
        self.pid = card.field(2)
        nids = card.fields(3,7)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==4

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodes
        return self.printCard(fields)

class CRAC2D(Element): # not integrated
    type = 'CRAC2D'
    def __init__(self,card):
        Element.__init__(self,card)
        self.pid = card.field(2)

        nids = card.fields(3,21) # caps at 18
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==18

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodes
        return self.printCard(fields)

class CRAC3D(Element): # not integrated
    type = 'CRAC3D'
    def __init__(self,card):
        Element.__init__(self,card)
        self.pid = card.field(2)

        nids = card.fields(3,67) # cap at +3 = 67
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==64

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodes
        return self.printCard(fields)
        
class CVISC(CROD): # not integrated
    type = 'CVISC'
    def __init__(self,card):
        CROD.__init__(self,card)
    ###
    def __repr__(self):  # not done...
        fields = ['CVISC',self.eid]
###

class CONM2(PointElement): # v0.1 not done
    """
    @todo support cid != 0
    """
    type = 'CONM2'
    # 'CONM2    501274  11064          132.274'
    def __init__(self,card):
        PointElement.__init__(self,card)
        #self.nids  = [ card[1] ]
        #del self.nids
        #self.pid = None
        self.eid  = card.field(1)
        self.gid  = card.field(2)
        self.cid  = card.field(3,0)
        self.Mass = card.field(4,0.0)
        self.X    = array(card.fields(5,8,[0.,0.,0.]))
        self.I    = card.fields(9,15,[0.]*6)

    def mass(self):
        return self.mass

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
        fields = [self.type,self.eid,self.gid,self.cid,self.Mass]+list(self.X)+self.I
        return self.printCard(fields)

   
