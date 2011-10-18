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
        self.loads = fields[3:] # alternating of scale factor & load set ID #
        assert len(self.loads)%2==0
    
    def __repr__(self):
        fields = ['LOAD',self.id,self.s]+self.loads
        return self.printCard(fields)
        

class FORCE(Load):
    def __init__(self,card):
#FORCE          3       1            100.      0.      0.      1.
        Load.__init__(self,card)
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

    def normalize(self):
        normXYZ = norm(self.xyz)
        mag = self.mag*normXYZ
        self.mag *= normXYZ
        self.xyz = self.xyz/normXYZ

    def __repr__(self):
        fields = ['FORCE',self.lid,self.node,self.cid,self.mag] + list(self.xyz)
        #print printCard(fields)
        return self.printCard(fields)

