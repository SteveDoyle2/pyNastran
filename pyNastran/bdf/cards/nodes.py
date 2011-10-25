#import sys
#from numpy import array,cross,dot
from numpy import array
#from numpy.linalg import norm

# my code
from baseCard import BaseCard

class Node(BaseCard): # base class
    def __init__(self,card):
        pass
    def crossReference(self,mesh):
        raise Exception('%s hasnt implemented a crossReference method' %(self.type))
    def __repr__(self):
        raise Exception('%s hasnt implemented a __repr__ method' %(self.type))

class SPOINT(Node):
    """
    SPOINT ID1 ID2 ID3 ID4 ID5 ID6 ID7 ID8
    or
    SPOINT ID1 THRU ID2
    SPOINT 5   THRU 649
    """
    type = 'SPOINT'
    def __init__(self,card):
        Node.__init__(self,card)
        fields  = card.fields(1)
        nFields = card.nFields()
        
        self.spoints = []
        i = 0
        while i<nFields: # =1 ???
            if fields[i]=='THRU':
                self.spoints += [fields[i-1],fields[i]+1]
                i+=1
            else:
                self.spoints.append(fields[i])
            i+=1
        ###
        
    def __repr__(self):
        ## @todo support THRU in output
        fields = ['SPOINT']+self.spoints
        return self.printCard(fields)
        

class GRDSET(Node):
    type = 'GRDSET'
    def __init__(self,card):
        self.cid  = card.field(2,0)
        self.ps   = card.field(7,0)
        self.seid = card.field(8,0)

    def __repr__(self):
        cid  = self.setBlankIfDefault(self.cid, 0)
        cd   = self.setBlankIfDefault(self.cd,  0)
        ps   = self.setBlankIfDefault(self.ps,  0)
        seid = self.setBlankIfDefault(self.seid,0)
        fields = ['GRDSET',None,cid,None,None,None,cd,ps,seid]

class GRID(Node):
    type = 'GRID'
    def __init__(self,card):
        Node.__init__(self,card)
        self.nid = int(card.field(1))
        self.cid = card.field(2,0)
        xyz = card.fields(3,6,[0.,0.,0.])  # TODO:  is standard nastran???
        #displayCard(card)
        #print "xyz = ",xyz

        self.xyz = array(xyz)
        self.cd = card.field(6,0)
        self.ps = card.field(7,0)
        self.seid = card.field(8,0)

        #print "xyz = ",self.xyz
        #print "cd = ",self.cd
        #print "ps = ",self.ps

    def Position(self):
        return self.xyzGlobal

    def crossReference(self,mesh,grdset=None):
        #print str(self)
        if grdset: # update using a gridset object
            if not self.cid:  self.cid  = grdset.cid
            if not self.cd:   self.cd   = grdset.cd
            if not self.ps:   self.ps   = grdset.ps
            if not self.seid: self.seid = grdset.seid
        
        coord = mesh.Coord(self.cid)
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

#class RINGAX
