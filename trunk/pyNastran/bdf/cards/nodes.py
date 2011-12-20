#import sys
#from numpy import array,cross,dot
from numpy import array
#from numpy.linalg import norm

# my code
from baseCard import BaseCard

class Ring(BaseCard): # base class
    def __init__(self,card,data):
        assert card is None or data is None

class Node(BaseCard): # base class
    def __init__(self,card,data):
        assert card is None or data is None

    def crossReference(self,mesh):
        raise Exception('%s hasnt implemented a crossReference method' %(self.type))

    def Cp(self):
        if isinstance(self.cp,int):
            return self.cp
        else:
            return self.cp.cid
        ###

    def Cd(self):
        if isinstance(self.cd,int):
            return self.cd
        else:
            return self.cd.cid
        ###
        
    def __repr__(self):
        raise Exception('%s hasnt implemented a __repr__ method' %(self.type))

class RINGAX(Ring):
    """
    Defines a ring for conical shell problems
    RINGAX ID R    Z    PS
    RINGAX 3  2.0 -10.0 162
    """
    type = 'RINGAX'
    def __init__(self,card):
        Node.__init__(self,card)
        self.nid = card.field(1)
        #2
        self.R   = card.field(3)
        self.z   = card.field(4)
        #5
        #6
        self.ps  = card.field(7)

    def Position(self):
        return array(0.,0.,0.)

    def __repr__(self):
        fields = ['RINGAX',self.nid,None,self.R,self.z,None,None,self.ps]
        return self.printCard(fields)
        
    
class SPOINT(Node):
    """
    SPOINT ID1 ID2 ID3 ID4 ID5 ID6 ID7 ID8
    or
    SPOINT ID1 THRU ID2
    SPOINT 5   THRU 649
    """
    type = 'SPOINT'
    def __init__(self,card=None,data=None):
        Node.__init__(self,card,data)
        fields  = card.fields(1)
        nFields = card.nFields()
        
        self.spoints = []
        self.dof = '23456'  # the constrained DOF
        i = 0
        while i<nFields: # =1 ???
            if fields[i]=='THRU':
                self.spoints += [fields[i-1],fields[i]+1]
                i+=1
            else:
                self.spoints.append(fields[i])
            i+=1
        ###

    def Position(self):
        return array(0.,0.,0.)

    def __repr__(self):
        ## @todo support THRU in output
        fields = ['SPOINT']+self.spoints
        return self.printCard(fields)
        

class GRDSET(Node):
    """
    Defines default options for fields 3, 7, 8, and 9 of all GRID entries.
    """
    type = 'GRDSET'
    def __init__(self,card=None,data=None):
        ## Grid point coordinate system
        self.cp  = card.field(2,0)
        
        ## Analysis coordinate system
        self.cd   = card.field(6,0)
        
        ## Default SPC constraint on undefined nodes
        self.ps   = card.field(7,0)
        
        ## Superelement ID
        self.seid = card.field(8,0)

    def crossReference(self):
        cp = mesh.Coord(self.cp)
        cd = mesh.Coord(self.cd)
        #seid = mesh.Super(self.seid)

    def __repr__(self):
        cp   = self.setBlankIfDefault(self.cp,  0)
        cd   = self.setBlankIfDefault(self.cd,  0)
        ps   = self.setBlankIfDefault(self.ps,  0)
        seid = self.setBlankIfDefault(self.seid,0)
        fields = ['GRDSET',None,cp,None,None,None,cd,ps,seid]
        return self.printCard(fields)

class GRIDB(Node):
    type = 'GRIDB'
    def __init__(self,card=None,data=None):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        Node.__init__(self,card,data)

        if card:
            raise Exception('not implemented...')
        else:
            print data
            self.nid  = data[0]
            self.phi  = data[1]
            self.cd   = data[2]
            self.ps   = data[3]
            self.idf  = data[4]
        ###
        assert self.nid > 0,  'nid=%s'  %(self.nid)
        assert self.phi >= 0, 'phi=%s'  %(self.phi)
        assert self.cd  >= 0, 'cd=%s'   %(self.cd)
        assert self.ps  >= 0, 'ps=%s'   %(self.ps)
        assert self.idf >= 0, 'idf=%s'  %(self.idf)

    def __repr__(self):
        #phi = self.setBlankIfDefault(self.phi,0.0)
        cd  = self.setBlankIfDefault(self.cd, 0)
        ps  = self.setBlankIfDefault(self.ps, 0)
        idf = self.setBlankIfDefault(self.idf,0)
        fields = ['GRIDB',self.nid,None,None,self.phi,None,cd,ps,idf]
        #print "fields = ",fields
        return self.printCard(fields)

class GRID(Node):
    type = 'GRID'
    def __init__(self,card=None,data=None):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        Node.__init__(self,card,data)

        if card:
            ## Node ID
            self.nid = int(card.field(1))

            ## Grid point coordinate system
            self.cp = card.field(2,0)

            xyz = card.fields(3,6,[0.,0.,0.])  ## @todo is standard nastran to set <0,0,0>as the defaults???
            #displayCard(card)
            #print "xyz = ",xyz
            self.xyz = array(xyz)

            ## Analysis coordinate system
            self.cd = card.field(6,0)

            ## SPC constraint
            self.ps = card.field(7,0)

            ## Superelement ID
            self.seid = card.field(8,0)

            #print "xyz = ",self.xyz
            #print "cd = ",self.cd
            #print "ps = ",self.ps
        else:
            #print data
            self.nid  = data[0]
            self.cp   = data[1]
            self.xyz  = array(data[2:5])
            self.cd   = data[5]
            self.ps   = data[6]
            self.seid = data[7]
            assert len(self.xyz)==3
        ###
        assert self.nid  >   0,  'nid=%s' %(self.nid)
        assert self.cp   >=  0,  'cp=%s'  %(self.cp)
        assert self.cd   >= -1,  'cd=%s'  %(self.cd)
        assert self.ps   >=  0,  'ps=%s'  %(self.ps)
        assert self.seid >=  0,'seid=%s'  %(self.seid)

    def Position(self,debug=False):
        #print type(self.cp)
        return self.cp.transformToGlobal(self.xyz,debug=debug)

    def PositionWRT(self,mesh,cid,debug=False):
        #print type(self.cp)
        coord = mesh.Coord(cid)
        return coord.transformToGlobal(self.xyz,debug=debug)

    def crossReference(self,mesh,grdset=None):
        #print str(self)
        if grdset: # update using a gridset object
            if not self.cp:   self.cp   = grdset.cp
            if not self.cd:   self.cd   = grdset.cd
            if not self.ps:   self.ps   = grdset.ps
            if not self.seid: self.seid = grdset.seid
        self.cp = mesh.Coord(self.cp)
        self.cd = mesh.Coord(self.cd)
        #self.xyzGlobal = coord.transformToGlobal(self.xyz)
        #return self.

    def __repr__(self):
        cp   = self.setBlankIfDefault(self.Cp(), 0)
        cd   = self.setBlankIfDefault(self.Cd(), 0)
        ps   = self.setBlankIfDefault(self.ps,  0)
        seid = self.setBlankIfDefault(self.seid,0)
        fields = ['GRID',self.nid,cp]+list(self.xyz)+[cd,ps,seid]
        #print "fields = ",fields
        return self.printCard(fields)

#class RINGAX
