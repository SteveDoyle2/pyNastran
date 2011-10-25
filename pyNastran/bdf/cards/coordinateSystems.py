#import sys
from numpy import array,cross,dot
from numpy.linalg import norm

# my code
#from pyNastran.bdf.bdf_helper import BDF_Card
from pyNastran.bdf.BDF_Card import BDF_Card
from baseCard import BaseCard

class Coord(BaseCard):
    def __init__(self,card):
        self.cid  = card.field(1)

        #self.type = card[0]
        #print "cid = ",self.cid

    def normalize(self,v):
        #print "v = ",v
        return v/norm(v)

    def __repr__(self):
        fields = [self.type,self.cid]
        return self.printCard(fields)

class CORD2R(Coord):  # working...needs more testing
    type = 'CORD2R'
    def __init__(self,card=['CORD2R',0,0,  0.,0.,0.,  0.,0.,1., 1.,0.,0.]):
        if isinstance(card,list):
            card = BDF_Card(card)

        Coord.__init__(self,card)
        self.rid = card.field(2,0)
        
        #print card
        self.eo = array( card.fields(3,6 ,[0.,0.,0.]) )
        self.ez = array( card.fields(6,9 ,[0.,0.,0.]) )
        self.ex = array( card.fields(9,12,[0.,0.,0.]) )
        assert len(self.eo)==3,self.eo
        assert len(self.ez)==3,self.ez
        assert len(self.ex)==3,self.ex
        
        #print "eo = ",self.eo
        #print "ez = ",self.ez
        #print "ex = ",self.ex

        self.ez0 = self.normalize(self.ez-self.eo)
        self.ex0 = self.normalize(self.ex-self.eo)
        self.ey0 = cross(self.ez0,self.ex0)
        
        #print norm(self.ex)
        #print norm(self.ey)
        #print norm(self.ez)
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
        return p

    def __repr__(self):
        eo = self.eo
        ez = self.ez
        ex = self.ex
        fields = ['CORD2R',self.cid,self.rid] +list(eo)+list(ez) +list(ex)

        #print "z1=%s z2=%s z3=%s" %(self.z1,self.z2,self.z3)
        #print "fields = ",fields
        return self.printCard(fields)


class CORD2C(Coord):  # not done...
    type = 'CORD2C'
    def __init__(self,card=['CORD2C',0,0,  0.,0.,0.,  0.,0.,1., 1.,0.,0.]):
        if isinstance(card,list):
            card = BDF_Card(card)

        Coord.__init__(self,card)
        self.rid = card.field(2,0)
        
        #print card
        self.eo = array( card.fields(3,6 ,[0.,0.,0.]) )
        self.ez = array( card.fields(6,9 ,[0.,0.,0.]) )
        self.ex = array( card.fields(9,12,[0.,0.,0.]) )
        assert len(self.eo)==3,self.eo
        assert len(self.ez)==3,self.ez
        assert len(self.ex)==3,self.ex
        
        #print "eo = ",self.eo
        #print "ez = ",self.ez
        #print "ex = ",self.ex

        self.ez0 = self.normalize(self.ez-self.eo)
        self.ex0 = self.normalize(self.ex-self.eo)
        self.ey0 = cross(self.ez0,self.ex0)

        #print card
        #print str(self)

        
    def transformToGlobal(self,p):
        (R,theta,z) = p
        raise Exception('cylindrical coordinate system...point R=%s theta=%s z=%s' %(R,theta,z))
        #assert R!=0.0
        x2 = self.eo[0]*self.ex+R*cos(theta)
        y2 = self.eo[1]*self.ex+R*sin(theta)
        z2 = self.eo[2]*self.ez+z*self.ez
        return array([x2,y2,z2])
        
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
        return p

    def __repr__(self):
        eo = self.eo
        ez = self.ez
        ex = self.ex
        fields = ['CORD2C',self.cid,self.rid] +list(eo)+list(ez) +list(ex)
        print "fields = ",fields
        raise Exception('not coded...')

        #print "z1=%s z2=%s z3=%s" %(self.z1,self.z2,self.z3)
        return self.printCard(fields)
