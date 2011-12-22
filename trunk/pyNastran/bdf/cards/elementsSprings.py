import sys
from numpy.linalg import norm

from baseCard import Element

class SpringElement(Element):
    def __init__(self,card,data):
        Element.__init__(self,card,data)

    def Length_noXref(self,n1=None,n2=None):
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

    def Centroid(self):
        p = (self.nodes[1].Position()-self.nodes[0].Position())/2.
        return p

    def K(self):
        raise Exception('K not implemented in the %s class' %(self.type))

    def Length(self):
        """
        Returns the length of a bar/rod/beam element
        \f[ \large \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  } \f]
        @param self the object pointer
        @note
            the model must be cross-referenced already
        """
        #print self.type
        return self.Length_noXref(self.nodes[1],self.nodes[0])

    def Mass(self):
        return 0.0

class CELAS1(SpringElement):
    type = 'CELAS1'
    def __init__(self,card=None,data=None):
        SpringElement.__init__(self,card,data)
        if card:
            self.eid = card.field(1)

            ## property ID
            self.pid = card.field(2,self.eid)

            nids = [card.field(3),card.field(5)]
            ## component number
            self.c1 = card.field(4)
            self.c2 = card.field(6)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[3]]
            self.c1 = data[4]
            self.c2 = data[5]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def K(self):
        return self.pid.k

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.pid   = model.Property(self.pid)
        
    def __repr__(self):
        nodes = self.nodeIDs()
        fields = ['CELAS1',self.eid,self.Pid(),nodes[0],self.c1,nodes[1],self.c2]
        return self.printCard(fields)

class CELAS2(SpringElement):
    type = 'CELAS2'
    def __init__(self,card=None,data=None):
        SpringElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)

            ## stiffness of the scalar spring
            self.k   = card.field(2)

            nids = [card.field(3),card.field(5)]

            ## component number
            self.c1 = card.field(4)
            self.c2 = card.field(6)

            ## damping coefficient
            self.ge = card.field(7)

            ## stress coefficient
            self.s  = card.field(8)
        else:
            self.eid = data[0]
            self.k   = data[1]
            nids     = [data[2],data[4]]
            self.c1  = data[3]
            self.c2  = data[5]
            self.ge  = data[6]
            self.s   = data[7]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def K(self):
        return self.k

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        
    def __repr__(self):
        nodes = self.nodeIDs()
        fields = ['CELAS2',self.eid,self.k,nodes[0],self.c1,nodes[1],self.c2,self.ge,self.s]
        return self.printCard(fields)

class CELAS3(SpringElement):
    type = 'CELAS3'
    def __init__(self,card=None,data=None):
        SpringElement.__init__(self,card,data)

        if card:
            #nids = [card.field(3),card.field(5)]
            #self.prepareNodeIDs(nids)
            #assert len(self.nodes)==2

            self.eid = card.field(1)
            ## property ID
            self.pid = card.field(2,self.eid)

            ## Scalar point identification numbers
            self.s1 = card.field(3)
            self.s2 = card.field(4)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.s1  = data[2]
            self.s2  = data[3]
        ###

    def K(self):
        return self.pid.k

    def crossReference(self,model):
        pass
        #self.nodes = model.Nodes(self.nodes)
        self.pid   = model.Property(self.pid)
        
    def __repr__(self):
        #nodes = self.nodeIDs()
        fields = ['CELAS3',self.eid,self.Pid(),self.s1,self.s2]
        return self.printCard(fields)

class CELAS4(SpringElement):
    type = 'CELAS4'
    def __init__(self,card=None,data=None):
        SpringElement.__init__(self,card,data)
        
        if card:
            #nids = [card.field(3),card.field(5)]
            #self.prepareNodeIDs(nids)
            #assert len(self.nodes)==2

            self.eid = card.field(1)

            ## stiffness of the scalar spring
            self.k   = card.field(2)

            ## Scalar point identification numbers
            self.s1 = card.field(3)
            self.s2 = card.field(4)
        else:
            self.eid = data[0]
            self.k   = data[1]
            self.s1  = data[2]
            self.s2  = data[3]
        ###

    def K(self):
        return self.k

    def crossReference(self,model):
        pass
        #self.nodes = model.Nodes(self.nodes)

    def __repr__(self):
        fields = ['CELAS4',self.eid,self.k,self.s1,self.s2]
        return self.printCard(fields)

