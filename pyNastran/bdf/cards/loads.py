#import sys
from numpy import array
from numpy.linalg import norm

from baseCard import BaseCard

class Load(BaseCard):
    type = 'DefaultLoad'
    def __init__(self,card):
        #self.type = card[0]
        self.lid  = card.field(1)

    #def normalize(self,v):
    #    #print "v = ",v
    #    return v/norm(v)

    def __repr__(self):
        fields = [self.type,self.lid]
        return self.printCard(fields)


class LOAD(Load):
    def __init__(self,card):
        fields = card.fields()
        self.lid = fields[1]
        self.id = self.lid

        self.s   = float(fields[2]) # overall scale factor
        nFields = len(fields)-2
        nLoads  = nFields/2
        print "nFields = ",nFields
        print "nLoads  = ",nLoads
        self.loads = fields[3:] # alternating of scale factor & load set ID
        assert len(self.loads)%2==0
    
    def __repr__(self):
        fields = ['LOAD',self.lid,self.s]+self.loads
        return self.printCard(fields)


class OneDeeLoad(Load): # FORCE/MOMENT
    def __init__(self,card):
        Load.__init__(self,card)

    def normalize(self):
        """
        adjust the vector to a unit length
        scale up the magnitude of the vector
        """
        if self.mag != 0.0:  # enforced displacement
            normXYZ = norm(self.xyz)
            mag = self.mag*normXYZ
            self.mag *= normXYZ
            self.xyz = self.xyz/normXYZ

class FORCE(OneDeeLoad):
    def __init__(self,card):
        """
        FORCE          3       1            100.      0.      0.      1.
        """
        OneDeeLoad.__init__(self,card)
        self.node = card.field(2)
        self.cid  = card.field(3,0)
        self.mag  = card.field(4)
        xyz = card.fields(5,8,[0.,0.,0.])
        assert len(xyz)==3,'xyz=%s' %(xyz)

        self.xyz = array(xyz)
        
        #print "node = ",self.node
        #print "mag  = ",self.mag
        #print "xyz  = ",self.xyz
        
        #self.type = card[0]
        #self.lid  = card[1]
        #print "card = ",card
        #print printCard(card)
        #print "self..."
        #print self
        #sys.exit()

    def __repr__(self):
        fields = ['FORCE',self.lid,self.node,self.cid,self.mag] + list(self.xyz)
        #print printCard(fields)
        return self.printCard(fields)

class PLOAD(Load):
    def __init__(self,card):
        self.lid = card.field(1)
        self.p   = card.field(2)
        self.nodes = card.fields(3,8)
        assert len(self.nodes)==4
    
    def __repr__(self):
        fields = ['PLOAD',self.lid,self.p]+self.nodes
        return self.printCard(fields)

class PLOAD1(Load):
    validTypes = ['FX','FY','FZ','FXE','FYE','FZE',
                  'MX','MY','MZ','MXE','MYE','MZE']
    validScales = ['LE','FR','LEPR','FRPR']
    def __init__(self,card):
        self.lid   = card.field(1)
        self.eid   = card.field(2)
        self.type  = card.field(3)
        self.scale = card.field(4)
        assert self.type in validTypes,  '%s is an invalid type on the PLOAD1 card' %(self.type)
        assert self.scale in validScales,'%s is an invalid scale on the PLOAD1 card' %(self.scale)
        self.x1 = card.field(5)
        self.p1 = card.field(6)
        self.x2 = card.field(7)
        self.p2 = card.field(8)
    
    def __repr__(self):
        fields = ['PLOAD1',self.lid,self.eid,self.type,self.scale,self.x1,self.p1,self.x2,self.p2]
        return self.printCard(fields)

class PLOAD2(Load):  # todo:  support THRU
    def __init__(self,card):
        self.lid = card.field(1)
        self.p   = card.field(2)
        self.nodes = card.fields(3,9)
        assert len(self.nodes)==6
        
        if card.field(4)=='THRU':
            print "found a THRU on PLOAD2"
    
    def __repr__(self):
        fields = ['PLOAD2',self.lid,self.p]+self.nodes
        return self.printCard(fields)

class PLOAD4(Load):  # todo:  support THRU, not done...
    def __init__(self,card):
        self.lid = card.field(1)
        self.eid = card.field(2)
        p1 = card.field(3)
        p  = card.fields(4,7,[p1,p1,p1])
        p = [p1]+p
        self.nodes = card.fields(3,9)
        assert len(self.nodes)==6

        if card.field(7)=='THRU':
            print "found a THRU on PLOAD4"
    
    def __repr__(self):
        fields = ['PLOAD4',self.lid,self.p]+self.nodes
        return self.printCard(fields)
