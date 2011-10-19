#import sys
from numpy import array,cross,dot
from numpy.linalg import norm

# my code
#from pyNastran.bdf.bdf_helper import BDF_Card
from BDF_Card import BDF_Card
from baseCard import BaseCard

class Coord(BaseCard):
    def __init__(self,card):

        #self.type = card[0]
        self.cid  = card.field(1)
        self.dunno = card.field(2)
        #print "cid = ",self.cid

    def normalize(self,v):
        #print "v = ",v
        return v/norm(v)

    def __repr__(self):
        fields = [self.type,self.cid]
        return self.printCard(fields)

class CORD2R(Coord):  # working...
    type = 'CORD2R'
    def __init__(self,card=['CORD2R',0,0,  0.,0.,0.,  0.,0.,1., 1.,0.,0.]):
        if isinstance(card,list):
            card = BDF_Card(card)

        Coord.__init__(self,card)
        
        #print card
        self.eo = array( card.fields(3,6 ,[0.,0.,0.]) )
        self.ez = array( card.fields(6,9 ,[0.,0.,0.]) )
        self.ex = array( card.fields(9,12,[0.,0.,0.]) )
        assert len(self.eo)==3,self.eo
        assert len(self.eo)==3,self.ez
        assert len(self.eo)==3,self.ex
        
        print "eo = ",self.eo
        print "ez = ",self.ez
        print "ex = ",self.ex

        self.ez0 = self.normalize(self.ez-self.eo)
        self.ex0 = self.normalize(self.ex-self.eo)
        self.ey0 = cross(self.ez0,self.ex0)
        
        print norm(self.ex)
        #print norm(self.ey)
        print norm(self.ez)
        #print card
        #print str(self)

        
    def transformToGlobal(self,p):
        if self.cid==0:
            return p
        
        p2 = p-self.eo
        
        #Bij = Bip*j
        ex = self.ex0
        ey = self.ey0
        ez = self.ez0
        gx = array([1.,0.,0.])
        gy = array([0.,1.,0.])
        gz = array([0.,0.,1.])
        
        matrix = array([[dot(gx,ex),dot(gx,ey),dot(gx,ez)],
                        [dot(gy,ex),dot(gy,ey),dot(gy,ez)],
                        [dot(gz,ex),dot(gz,ey),dot(gz,ez)]])
        print "p = ",p
        print "matrix = ",matrix
        p2 = dot(p,matrix)
        p3 = p2+self.eo
        print "p2 = ",p2
        
        #print str(self)
        #sys.exit('stop')
        return p

    def __repr__(self):
        eo = self.eo
        ez = self.ez
        ex = self.ex
        fields = [self.type,self.cid,self.dunno] +list(eo)+list(ez) +list(ex)

        #print "z1=%s z2=%s z3=%s" %(self.z1,self.z2,self.z3)
        #print "fields = ",fields
        return self.printCard(fields)
