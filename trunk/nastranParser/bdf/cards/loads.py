#import sys
from numpy import array
from numpy.linalg import norm

class Load(object):
    def __init__(self,card):
        #self.type = card[0]
        self.lid  = card[1]

    #def normalize(self,v):
    #    #print "v = ",v
    #    return v/norm(v)

    def __repr__(self):
        fields = [self.type,self.lid]
        return printCard(fields)


class FORCE(Load):
    def __init__(self,card):
#FORCE          3       1            100.      0.      0.      1.
        Load.__init__(self,card)
        self.node = card[2]
        self.cid  = card[3]
        self.mag  = card[4]
        self.xyz = array(card[5:8])
        
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
        return printCard(fields)

