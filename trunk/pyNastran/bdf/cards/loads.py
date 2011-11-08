#import sys
from numpy import array,cross
from numpy.linalg import norm

from baseCard import BaseCard

class Load(BaseCard):
    """defines the DefaultLoad class"""
    type = 'DefLoad'
    def __init__(self,card):
        #self.type = card[0]
        self.lid  = card.field(1)

    #def normalize(self,v):
    #    #print "v = ",v
    #    return v/norm(v)

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        else:
            return self.cid.cid
        ###

    def nodeIDs(self,nodes=None):
        """returns nodeIDs for repr functions"""
        if not nodes:
           nodes = self.nodes
        if isinstance(nodes[0],int):
            #print 'if'
            return [node     for node in nodes]
        else:
            #print 'else'
            return [node.nid for node in nodes]
        ###

    def __repr__(self):
        fields = [self.type,self.lid]
        return self.printCard(fields)


class LSEQ(BaseCard):  # how does this work...
    """
    Defines a sequence of static load sets
    """
    type = 'LSEQ'
    def __init__(self,card):
        self.sid  = card.field(1)
        self.exciteID = card.field(2)
        self.lid = card.field(3)
        self.tid = card.field(4)


    def nodeIDs(self,nodes=None):
        """returns nodeIDs for repr functions"""
        if not nodes:
           nodes = self.nodes
        if isinstance(nodes[0],int):
            #print 'if'
            return [node     for node in nodes]
        else:
            #print 'else'
            return [node.nid for node in nodes]
        ###

    def __repr__(self):
        fields = ['LSEQ',self.lid]
        return self.printCard(fields)

class LOAD(Load):
    type = 'LOAD'
    def __init__(self,card):
        """
        @todo parse the loads data to have scale factor and load
        """
        fields   = card.fields()
        
        ## load ID
        self.lid = card.field(1)
        #self.id  = self.lid

        ## overall scale factor
        self.s  = card.field(2)
        
        nFields = len(fields)-2
        nLoads  = nFields/2
        #print "nFields = ",nFields
        #print "nLoads  = ",nLoads
        
        loads = card.fields(3) # temp list
        assert len(loads)%2==0
        
        ## individual scale factors (corresponds to loadIDs)
        self.scaleFactors = []

        ## individual loadIDs (corresponds to scaleFactors)
        self.loadIDs = []

        for i in range(0,nLoads,2):   # alternating of scale factor & load set ID
            self.scaleFactors.append(loads[i  ])
            self.loadIDs.append(     loads[i+1])
    
    def __repr__(self):
        fields = ['LOAD',self.lid,self.s]
        for scaleFactor,loadID in zip(self.scaleFactors,self.loadIDs):
            fields += [scaleFactor,loadID]
        return self.printCard(fields)


class OneDeeLoad(Load): # FORCE/MOMENT
    type = '1D_Load'
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

class Force(OneDeeLoad):
    type = '1D_Load'
    def __init__(self,card):
        OneDeeLoad.__init__(self,card)

    def F(self):
        return self.xyz*self.mag

class Moment(OneDeeLoad):
    type = 'Moment'
    def __init__(self,card):
        OneDeeLoad.__init__(self,card)

    def M(self):
        return self.xyz*self.mag


class FORCE(Force):
    type = 'FORCE'
    def __init__(self,card):
        """
        FORCE          3       1            100.      0.      0.      1.
        """
        Force.__init__(self,card)
        self.node = card.field(2)
        self.cid  = card.field(3,0)
        self.mag  = card.field(4)

        xyz = card.fields(5,8,[0.,0.,0.])
        assert len(xyz)==3,'xyz=%s' %(xyz)
        self.xyz = array(xyz)

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        pass

    def __repr__(self):
        cid = self.setBlankIfDefault(self.cid,0)
        fields = ['FORCE',self.lid,self.node,cid,self.mag] + list(self.xyz)
        return self.printCard(fields)

class FORCE1(Force):
    """
    Defines a static concentrated force at a grid point by specification of a magnitude and
    two grid points that determine the direction.
    """
    type = 'FORCE1'
    def __init__(self,card):
        Force.__init__(self,card)
        self.node = card.field(2)
        self.mag  = card.field(3)
        self.g1   = card.field(4)
        self.g2   = card.field(5)

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        self.node = mdodel.Node(self.node)

        self.xyz = model.Node(self.g2).Position() - model.Node(self.g1).Position()
        self.Normalize()

    def NodeID(self):
        if isinstance(self.node,int):
            return self.node
        return self.node.nid

    def __repr__(self):
        (node,g1,g2) = self.nodeIDs([self.node,self.g1,self.g2])
        fields = ['FORCE1',self.lid,node,self.mag,g1,g2]
        return self.printCard(fields)

class FORCE2(Force):
    """
    Defines a static concentrated force at a grid point by specification of a magnitude and
    four grid points that determine the direction.
    """
    type = 'FORCE2'
    def __init__(self,card):
        """
        FORCE2 SID G F G1 G2 G3 G4
        """
        Force.__init__(self,card)
        self.node = card.field(2)
        self.mag  = card.field(3)
        self.g1   = card.field(4)
        self.g2   = card.field(5)
        self.g3   = card.field(5)
        self.g4   = card.field(5)

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        self.node = mdodel.Node(self.node)

        v12 = model.Node(self.g2).Position() - model.Node(self.g1).Position()
        v34 = model.Node(self.g4).Position() - model.Node(self.g3).Position()
        v12 = v12/norm(v12)
        v34 = v34/norm(v34)
        self.xyz = cross(v12,v34)
        self.Normalize()

    def NodeID(self):
        if isinstance(self.node,int):
            return self.node
        return self.node.nid

    def __repr__(self):
        (node,g1,g2,g3,g4) = self.nodeIDs([self.node,self.g1,self.g2,self.g3,self.g4])
        fields = ['FORCE2',self.lid,node,self.mag,g1,g2,g3,g4]
        return self.printCard(fields)

class MOMENT(Moment):    # can i copy the force init without making the MOMENT a FORCE ???
    type = 'MOMENT'
    def __init__(self,card):
        """
        Defines a static concentrated moment at a grid point by specifying a scale factor and
        a vector that determines the direction.

        MOMENT SID G CID M    N1  N2  N3
        MOMENT 2   5   6 2.9 0.0 1.0 0.0
        """
        Moment.__init__(self,card)
        self.node = card.field(2)
        self.cid  = card.field(3,0)
        self.mag  = card.field(4)

        xyz = card.fields(5,8,[0.,0.,0.])
        assert len(xyz)==3,'xyz=%s' %(xyz)
        self.xyz = array(xyz)

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        pass

    def __repr__(self):
        cid = self.setBlankIfDefault(self.cid,0)
        fields = ['MOMENT',self.lid,self.node,cid,self.mag] + list(self.xyz)
        return self.printCard(fields)

class MOMENT1(Moment):
    type = 'MOMENT1'
    def __init__(self,card):
        """
        Defines a static concentrated moment at a grid point by specifying a magnitude and
        two grid points that determine the direction

        MOMENT1 SID G M G1 G2
        """
        Moment.__init__(self,card)
        self.node = card.field(2)
        self.mag  = card.field(3)
        self.g1   = card.field(4)
        self.g2   = card.field(5)
        self.g3   = card.field(6)
        self.g4   = card.field(7)

        xyz = card.fields(5,8,[0.,0.,0.])
        assert len(xyz)==3,'xyz=%s' %(xyz)
        self.xyz = array(xyz)

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        self.node = mdodel.Node(self.node)
        self.xyz = model.Node(self.g2).Position() - model.Node(self.g1).Position()
        self.Normalize()

    def __repr__(self):
        (node,g1,g2) = self.nodeIDs([self.node,self.g1,self.g2])
        fields = ['MOMENT1',self.lid,node,self.mag,g1,g2]
        return self.printCard(fields)


class MOMENT2(Moment):
    type = 'MOMENT2'
    def __init__(self,card):
        """
        Defines a static concentrated moment at a grid point by specification of a magnitude
        and four grid points that determine the direction.

        MOMENT2 SID G M G1 G2 G3 G4
        """
        Moment.__init__(self,card)
        self.node = card.field(2)
        self.mag  = card.field(3)
        self.g1   = card.field(4)
        self.g2   = card.field(5)
        self.g3   = card.field(6)
        self.g4   = card.field(7)

        xyz = card.fields(5,8,[0.,0.,0.])
        assert len(xyz)==3,'xyz=%s' %(xyz)
        self.xyz = array(xyz)

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        (self.g1,self.g2,self.g3,self.g4) = model.Nodes(self.g1,self.g2,self.g3,self.g4)
        v12 = g2.Position()-g1.Position()
        v34 = g4.Position()-g3.Position()
        v12 = v12/norm(v12)
        v34 = v34/norm(v34)
        self.xyz = cross(v12,v34)

    def __repr__(self):
        (node,g1,g2,g3,g4) = self.nodeIDs([self.node,self.g1,self.g2,self.g3,self.g4])
        fields = ['MOMENT2',self.lid,node,self.mag,g1,g2,g3,g4]
        return self.printCard(fields)

class PLOAD(Load):
    type = 'PLOAD'
    def __init__(self,card):
        self.lid = card.field(1)
        self.p   = card.field(2)
        self.nodes = card.fields(3,8)
        assert len(self.nodes)==4
    
    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        pass

    def __repr__(self):
        fields = ['PLOAD',self.lid,self.p]+self.nodeIDs()
        return self.printCard(fields)

class PLOAD1(Load):
    type = 'PLOAD'
    validTypes = ['FX','FY','FZ','FXE','FYE','FZE',
                  'MX','MY','MZ','MXE','MYE','MZE']
    validScales = ['LE','FR','LEPR','FRPR']
    def __init__(self,card):
        self.lid   = card.field(1)
        self.eid   = card.field(2)
        self.Type  = card.field(3)
        self.scale = card.field(4)
        assert self.Type  in self.validTypes, '%s is an invalid type on the PLOAD1 card' %(self.type)
        assert self.scale in self.validScales,'%s is an invalid scale on the PLOAD1 card' %(self.scale)
        self.x1 = card.field(5)
        self.p1 = card.field(6)
        self.x2 = card.field(7)
        self.p2 = card.field(8)
    
    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        pass

    def __repr__(self):
        fields = ['PLOAD1',self.lid,self.eid,self.Type,self.scale,self.x1,self.p1,self.x2,self.p2]
        return self.printCard(fields)

class PLOAD2(Load):  # todo:  support THRU
    type = 'PLOAD2'
    def __init__(self,card):
        self.lid = card.field(1)
        self.p   = card.field(2)
        self.nodes = card.fields(3,9)
        assert len(self.nodes)==6
        
        if card.field(4)=='THRU':
            print "found a THRU on PLOAD2"
            pass
    
    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        pass

    def __repr__(self):
        fields = ['PLOAD2',self.lid,self.p]+self.nodeIDs()
        return self.printCard(fields)

class PLOAD4(Load):
    type = 'PLOAD4'
    def __init__(self,card):
        self.lid = card.field(1)
        self.eid = card.field(2)
        p1 = card.field(3)
        p  = card.fields(4,7,[p1,p1,p1])
        self.p = [p1]+p
        
        if card.field(7)=='THRU':
            #print "found a THRU on PLOAD4"
            pass
            eid2 = card.field(8)
            self.eids= self.expandThru([self.eid,'THRU',eid2])
            self.g3 = None
            self.g4 = None
        else:   # used for CPENTA, CHEXA
            self.eids = None
            self.g3 = card.field(7)
            self.g4 = card.field(8)
        ###
        
        ## Coordinate system identification number. See Remark 2. (Integer >= 0;Default=0)
        self.cid     = card.field(9,0)
        self.NVector = card.fields(10,13,[0.,0.,0.])
        self.sorl    = card.field(13,'SURF')
    
    def crossReference(self,model):
        self.cid = model.Coord(self.cid)
        if self.g1: self.g1 = model.Node(self.g1)
        if self.g3: self.g3 = model.Node(self.g3)
        if self.g4: self.g4 = model.Node(self.g4)
        if self.eids:
            self.eids = model.Elements(self.eids)

    def __repr__(self):
        cid  = self.setBlankIfDefault(self.cid,0)
        sorl = self.setBlankIfDefault(self.sorl,'SURF')
        p2 = self.setBlankIfDefault(self.p[1],self.p[0])
        p3 = self.setBlankIfDefault(self.p[2],self.p[0])
        p4 = self.setBlankIfDefault(self.p[3],self.p[0])
        fields = ['PLOAD4',self.lid,self.p[0],p2,p3,p4]

        if self.g3:
            (g3,g4) = self.nodeIDs([self.g3,self.g4])
            fields.append(g3)
            fields.append(g4)
        else:
            fields.append('THRU')
            eid = self.eids[-1]
            fields.append(self.getElementIDs([eid]) )
        fields.append(cid)
        fields += list(self.NVector)
        fields.append(sorl)
        #print "fields = ",fields
        return self.printCard(fields)
