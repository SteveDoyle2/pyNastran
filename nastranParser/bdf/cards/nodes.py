#import sys
#from numpy import array,cross,dot
from numpy import array
#from numpy.linalg import norm

# my code
from baseCard import BaseCard

class GRID(BaseCard):
    def __init__(self,card):
        self.nid = int(card.field(1))
        self.id  = self.nid

        self.cid = card.field(2,0)
        xyz = card.fields(3,6,[0.,0.,0.])  # TODO:  is standard nastran???
        #displayCard(card)
        #print "xyz = ",xyz

        self.xyz = array(xyz)
        self.cd = card.field(6,0)
        self.ps = card.field(7,0)
        self.seid = card.field(8,0)

    def Position(self):
        return self.xyzGlobal

    def crossReference(self,coord):
        #print str(self)
        self.xyzGlobal = coord.transformToGlobal(self.xyz)
        #return self.

    def __repr__(self):
        cid  = self.setBlankIfDefault(self.cid, 0)
        cd   = self.setBlankIfDefault(self.cd,  0)
        ps   = self.setBlankIfDefault(self.ps,  0)
        seid = self.setBlankIfDefault(self.seid,0)
        fields = ['GRID',self.nid,cid]+list(self.xyz)+[cd,ps,seid]
        #print "fields = ",fields
        return self.printCard(fields)

#class SPOINT
#class RINGAX
