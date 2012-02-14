#import sys
from numpy import array
#from numpy import cross,dot
#from numpy.linalg import norm

# my code
from baseCard import Element

from elementsRigid  import *
from elementsSolid  import *
from bars.elementsBars       import *
from plates.elementsShell    import *
from springs.elementsSprings import *

class BushingElement(Element):
    def __init__(self,card,data):
        Element.__init__(self,card,data)

class CBUSH(BushingElement):
    type = 'CBUSH'
    def __init__(self,card=None,data=None):
        BushingElement.__init__(self,card,data)
        
        if card:
            fields = card.fields()
            #print "fields = ",fields
            self.eid = card.field(1)
            self.pid = card.field(2)
            self.ga = card.field(3)
            self.gb = card.field(4)
            x1G0    = card.field(5)
            if isinstance(x1G0,int):
                self.g0 = x1G0
                self.x = None
            elif isinstance(x1G0,float):
                self.g0  = None
                x1  = x1G0
                x2  = card.field(6)
                x3  = card.field(7)
                self.x   = [x1,x2,x3]
            else:
                #raise Exception('invalid CGAP...x1/g0 = |%s|' %(x1G0))
                self.g0 = None
                self.x  = [None,None,None]
            ###
            ## Element coordinate system identification. A 0 means the basic
            ## coordinate system. If CID is blank, then the element coordinate system
            ## is determined from GO or Xi.
            self.cid  = card.field(8,0)
            ## Location of spring damper (0 <= s <= 1.0)
            self.s    = card.field(10,0.5)
            ## Coordinate system identification of spring-damper offset. See
            ## Remark 9. (Integer > -1; Default = -1, which means the offset point lies
            ## on the line between GA and GB
            self.ocid = card.field(11,-1)
            ## Components of spring-damper offset in the OCID coordinate system if OCID > 0.
            self.si   = card.fields(i=12,j=14)
        else:
            self.eid = data[0]
            raise NotImplementedError('CBUSH data...')
        #self.prepareNodeIDs(nids,allowEmptyNodes=True)
        #assert len(self.nodes)==2

    #def OCid(self):
        #if isinstance(self.ocid,int):
            #return self.ocid
        #return self.ocid.cid

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        return self.cid.cid

    def crossReference(self,model):
        #self.nodes = model.Nodes(self.nodes)
        self.pid = model.Property(self.pid)
        self.cid = model.Coord(self.cid)
        
    def rawFields(self):
        if self.g0 is not None:
            x = [self.g0,None,None]
        else:
            x = self.x
        fields = ['CBUSH',self.eid,self.Pid(),self.ga,self.gb]+x+[self.Cid(),
                          self.s,self.ocid]+self.si
        return fields

    def reprFields(self):
        return self.rawFields()

class DamperElement(Element):
    def __init__(self,card,data):
        Element.__init__(self,card,data)

class CDAMP1(DamperElement):
    type = 'CDAMP1'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)

            nids = [card.field(3),card.field(5)]

            ## component number
            self.c1 = card.field(4)
            self.c2 = card.field(6)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[4]]
            self.c1  = data[3]
            self.c2  = data[5]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def isSameCard(self,elem):
        if self.type!=elem.type:  return False
        fields1 = [self.eid,self.Pid()]+self.nodeIDs()+[self.c1,self.c2]
        fields2 = [elem.eid,elem.Pid()]+elem.nodeIDs()+[elem.c1,elem.c2]
        return self.isSameFields(fields1,fields2)

    def B(self):
        return self.pid.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['CDAMP1',self.eid,self.Pid(),nodes[0],self.c1,nodes[1],self.c2]
        return fields

class CDAMP2(DamperElement):
    type = 'CDAMP2'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)

            ## Value of the scalar damper (Real)
            self.b   = card.field(2)

            nids = [card.field(3),card.field(5)]

            ## component number
            self.c1 = card.field(4)
            self.c2 = card.field(6)
        else:
            self.eid = data[0]
            self.b   = data[1]
            nids     = [data[2],data[4]]
            self.c1  = data[3]
            self.c2  = data[5]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def B(self):
        return self.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        
    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['CDAMP2',self.eid,self.b,nodes[0],self.c1,nodes[1],self.c2]
        return fields

class CDAMP3(DamperElement):
    type = 'CDAMP3'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = [card.field(3),card.field(4)]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[3]]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def B(self):
        return self.pid.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['CDAMP3',self.eid,self.pid,nodes[0],nodes[1]]
        return fields

class CDAMP4(DamperElement):
    type = 'CDAMP4'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)
            ## Value of the scalar damper (Real)
            self.b   = card.field(2)
            nids = [card.field(3),card.field(4)]
        else:
            self.eid = data[0]
            self.b   = data[1]
            nids     = [data[2],data[3]]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def B(self):
        return self.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        
    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['CDAMP4',self.eid,self.b,nodes[0],nodes[1]]
        return fields

class CDAMP5(DamperElement):
    type = 'CDAMP5'
    def __init__(self,card=None,data=None):
        DamperElement.__init__(self,card,data)
        
        if card:
            self.eid = card.field(1)
            ## Property ID
            self.pid = card.field(2)
            nids = [card.field(3),card.field(4)]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[3]]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['CDAMP5',self.eid,self.Pid(),nodes[0],nodes[1]]
        return fields

class CGAP(Element):
    type = 'CGAP'
    def __init__(self,card=None,data=None):
        Element.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            self.ga = card.field(3)
            self.gb = card.field(4)
            x1G0    = card.field(5)
            if isinstance(x1G0,int):
                self.g0 = x1G0
                self.x = None
            elif isinstance(x1G0,float):
                self.g0  = None
                x1  = x1G0
                x2  = card.field(6)
                x3  = card.field(7)
                self.x   = [x1,x2,x3]
                self.cid = card.field(8,0)
            else:
                raise Exception('invalid CGAP...x1/g0 = |%s|' %(x1G0))
            ###
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.ga  = data[2]
            self.gb  = data[3]
            self.g0  = data[4]
            x1  = data[5]
            x2  = data[6]
            x3  = data[7]
            self.x   = [x1,x2,x3]
            self.cid = data[8]
        ###

    def crossReference(self,mesh):
        self.ga = mesh.Node(self.ga)
        self.gb = mesh.Node(self.gb)
        if self.g0:
            self.g0 = mesh.Node(self.g0)
            self.x  = self.g0.Position()
        self.pid = mesh.Property(self.pid)
        self.cid = mesh.Coord(self.cid)

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        return self.cid.cid

    def rawFields(self):
        if self.g0 is not None:
            x = [self.g0,None,None]
        else:
            x = self.x
        fields = ['CGAP',self.eid,self.Pid(),self.ga,self.gb]+x+[self.Cid()]
        return fields

class CRAC2D(Element):
    type = 'CRAC2D'
    def __init__(self,card=None,data=None):
        Element.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,21) # caps at 18
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==18

    def rawFields(self):
        fields = ['CRAC2D',self.eid,self.Pid()]+self.nodes
        return fields

class CRAC3D(Element):
    type = 'CRAC3D'
    def __init__(self,card=None,data=None):
        Element.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,67) # cap at +3 = 67
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==64

    def rawFields(self):
        fields = ['CRAC3D',self.eid,self.Pid()]+self.nodes
        return fields

class PointElement(Element):
    def __init__(self,card=None,data=None):
        Element.__init__(self,card,data)
        
    def Mass(self):
        return self.mass

class CMASS1(PointElement):
    """
    Defines a scalar mass element.
    CMASS2 EID M   G1 C1 G2 C2
    CMASS1 EID PID G1 C1 G2 C2
    """
    type = 'CMASS2'
    def __init__(self,card=None,data=None):
        PointElement.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2,self.eid)
            self.g1 = card.field(3)
            self.c1 = card.field(4)
            self.g2 = card.field(5)
            self.c2 = card.field(6)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.g1  = data[2]
            self.c1  = data[3]
            self.g2  = data[4]
            self.c2  = data[5]
        ###

    def Mass(self):
        return self.pid.mass

    def crossReference(self,mesh):
        #self.g1 = mesh.Node(self.g1)
        #self.g2 = mesh.Node(self.g2)
        self.pid = mesh.Property(self.pid)

    def rawFields(self):
        fields = ['CMASS1',self.eid,self.Pid(),self.g1,self.c1,self.g2,self.c2]
        return fields

class CMASS2(PointElement):
    """
    Defines a scalar mass element without reference to a property entry.
    CMASS2 EID M G1 C1 G2 C2
    """
    type = 'CMASS2'
    def __init__(self,card=None,data=None):
        PointElement.__init__(self,card,data)
        if card:
            self.eid  = card.field(1)
            self.mass = card.field(2,0.)
            self.g1   = card.field(3)
            self.c1   = card.field(4)
            self.g2   = card.field(5)
            self.c2   = card.field(6)
        else:
            self.eid  = data[0]
            self.mass = data[1]
            self.g1   = data[2]
            self.c1   = data[3]
            self.g2   = data[4]
            self.c2   = data[5]
        ###

    def Mass(self):
        return self.mass

    def crossReference(self,mesh):
        """
        """
        #self.g1 = mesh.Node(self.g1)
        #self.g2 = mesh.Node(self.g2)
        pass

    def rawFields(self):
        fields = ['CMASS2',self.eid,self.mass,self.g1,self.c1,self.g2,self.c2]
        return fields

    def reprFields(self):
        mass = self.setBlankIfDefault(self.mass,0.0)
        fields = ['CMASS2',self.eid,mass,self.g1,self.c1,self.g2,self.c2]
        return fields

class CMASS3(PointElement):
    """
    Defines a scalar mass element that is connected only to scalar points.
    CMASS3 EID PID S1 S2
    """
    type = 'CMASS3'
    def __init__(self,card=None,data=None):
        PointElement.__init__(self,card,data)

        if card:
            self.eid  = card.field(1)
            self.pid  = card.field(2,self.eid)
            self.s1   = card.field(3)
            self.s2   = card.field(4)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.s1  = data[2]
            self.s2  = data[3]
        ###

    def Mass(self):
        return self.pid.mass

    def isSameCard(self,elem):
        if self.type!=elem.type:  return False
        fields1 = [self.eid,self.Pid(),self.s1,self.s2]
        fields2 = [elem.eid,elem.Pid(),elem.s1,elem.s2]
        return self.isSameFields(fields1,fields2)

    def crossReference(self,mesh):
        """
        links up the propertiy ID
        """
        #self.s1 = mesh.Node(self.s1)
        #self.s2 = mesh.Node(self.s2)
        self.pid = mesh.Property(self.pid)

    def rawFields(self):
        fields = ['CMASS3',self.eid,self.Pid(),self.s1,self.s2]
        return fields

class CMASS4(PointElement):
    """
    Defines a scalar mass element that is connected only to scalar points, without reference to a property entry
    CMASS4 EID M S1 S2
    """
    type = 'CMASS4'
    def __init__(self,card=None,data=None):
        PointElement.__init__(self,card,data)
        
        if card:
            self.eid  = card.field(1)
            self.mass = card.field(2,0.)
            self.s1 = card.field(3)
            self.s2 = card.field(4)
        else:
            self.eid  = data[0]
            self.mass = data[1]
            self.s1   = data[2]
            self.s2   = data[3]
        ###

    def Mass(self):
        return self.mass

    def isSameCard(self,elem):
        if self.type!=elem.type:  return False
        fields1 = [self.eid,self.Pid(),self.s1,self.s2]
        fields2 = [elem.eid,elem.Pid(),elem.s1,elem.s2]
        return self.isSameFields(fields1,fields2)

    def crossReference(self,mesh):
        """
        """
        #self.s1 = mesh.Node(self.s1)
        #self.s2 = mesh.Node(self.s2)
        pass

    def rawFields(self):
        fields = ['CMASS4',self.eid,self.mass,self.s1,self.s2]
        return fields

   
class CONM2(PointElement): # v0.1 not done
    """
    @todo support cid != 0
    """
    type = 'CONM2'
    # 'CONM2    501274  11064          132.274'
    def __init__(self,card=None,data=None):
        PointElement.__init__(self,card,data)
        if card:
            #self.nids  = [ card[1] ]
            #del self.nids
            #self.pid = None
            self.eid  = card.field(1)
            self.nid  = card.field(2)
            self.cid  = card.field(3,0)
            self.mass = card.field(4,0.0)
            self.X    = array(card.fields(5,8,[0.,0.,0.]))
            self.I    = card.fields(9,15,[0.]*6)
        else:
            self.eid  = data[0]
            self.nid  = data[1]
            self.cid  = data[2]
            self.mass = data[3]
            self.X    = data[4:7]
            self.I    = data[7:]
        ###

    def Mass(self):
        if self.cid==0:
            raise NotImplementedError('cid=0 is not supported for CONM2...')
        return self.mass

    def Centroid(self):
        if self.cid==0:
            raise NotImplementedError('cid=0 is not supported for CONM2...')
        if self.cid==-1:
            return self.X
        else:
            dx,matrix = self.cid.transformToGlobal(self.X)
            #print "dx = ",dx
            #print "X  = ",self.nid.Position()
            X2  = self.nid.Position()+dx
        return X2

    def crossReference(self,model):
        """
        @warning only supports cid=0
        """
        if self.cid!=-1:
            self.cid = model.Coord(self.cid)
        self.nid = model.Node(self.nid)

    def Nid(self):
        if isinstance(self.nid,int):
            return self.nid
        return self.nid.nid

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        return self.cid.cid

    def writeCodeAster(self):
        msg += "    DISCRET=_F(\n"
        msg += "             'CARA='M_T_D_N'\n"
        msg += "              NOEUD=N%s\n" %(self.Nid())
        msg += "              VALE=%g),\n" %(self.mass)
        return msg
    def rawFields(self):
        #print "X=%r" %(self.X)
        #print "I=%r" %(self.I)
        fields = ['CONM2',self.eid,self.Nid(),self.Cid(),self.mass]+list(self.X)+[None]+list(self.I)
        return fields

    def reprFields(self):
        I = []
        for i in self.I:
            if i==0.:
                I.append(None)
            else:
                I.append(i)
            ###
        ###
        X = []
        for x in self.X:
            if x==0.:
                X.append(None)
            else:
                X.append(x)
            ###
        ###

        cid  = self.setBlankIfDefault(self.Cid(),0)
        mass = self.setBlankIfDefault(self.mass,0.0)
        fields = ['CONM2',self.eid,self.Nid(),cid,mass]+X+[None]+I
        return fields

