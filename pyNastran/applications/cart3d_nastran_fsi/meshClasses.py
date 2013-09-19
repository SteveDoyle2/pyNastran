from numpy import array,cross,dot
from numpy.linalg import norm
import sys

# my code
from bdf import printCard
from mathFunctions import Area,Triangle_AreaCentroidNormal,Normal

def setBlankIfDefault(value,default):
    if value==default:
        return ''
    return value

def setDefaultIfBlank(value,default):
    if value=='':
        return default
    return value

class PARAM(object):
    def __init__(self,card):
        self.key   = card[1]
        self.value = card[2]

    def __repr__(self):
        fields = ['PARAM',self.key,self.value]
        return printCard(fields)

class GRID(object):
    def __init__(self,card):
        self.nid = int(card[1])
        self.id  = self.nid

        self.cid = setDefaultIfBlank(card[2],0)
        xyz = card[3:6]  # TODO:  is standard nastran???
        #displayCard(card)
        #print "xyz = ",xyz

        self.xyz = array([0.,0.,0.])
        for i,x in enumerate(xyz):
            self.xyz[i] = float(x)
        self.xyz = array(self.xyz)

    def Position(self):
        return self.xyzGlobal

    def crossReference(self,coord):
        #print str(self)
        self.xyzGlobal = coord.transformToGlobal(self.xyz)
        #return self.

    def __repr__(self):
        cid = setBlankIfDefault(self.cid,0)
        fields = ['GRID',self.nid,cid]+list(self.xyz)
        #print "fields = ",fields
        return printCard(fields)

class Constraint(object):
    def __init__(self,card):
        #self.type = card[0]
        self.lid  = card[1]

    def cleanNodes(self,nodes):
        nodes2 = []
        for node in nodes:
            if node=="":
                pass
            else:
                nodes2.append(node)
            ###
        ###
        #print "nodes2 = ",nodes2
        nodes = nodes2
        if len(nodes2)>1 and nodes2[1]=='THRU':
            nodes = [int(i) for i in range(nodes2[0],nodes2[2]+1)]

        self.nodes = nodes
        #print "*nodes = ",nodes
        #return nodes2

    def __repr__(self):
        fields = [self.type,self.cid]
        return printCard(fields)

class SPC1(Constraint):
#SPC1     3       246     209075  209096  209512  209513  209516
    type = 'SPC1'
    def __init__(self,card):
        Constraint.__init__(self,card)
        self.constrained = card[2]  # 246 = y; dx, dz dir
        nodes = card[3:]
        self.cleanNodes(nodes)
        #print "nodes = ",nodes

    def __repr__(self):
        test = [i for i in range(self.nodes[0],self.nodes[-1]+1)]
        #print "self.nodes = ",self.nodes
        #print "test       = ",test
        if self.nodes==test:
            nodes = [self.nodes[0],'THRU',self.nodes[-1]]
        else:
            nodes = [int(i) for i in self.nodes]
        fields = [self.type,self.lid,self.constrained]+nodes
        return printCard(fields)

class SPCADD(Constraint):
#SPCADD   2       1       3
    type = 'SPCADD'
    def __init__(self,card):
        Constraint.__init__(self,card)
        nodes = card[2:]
        self.cleanNodes(nodes)
        #print "self.nodes = ",self.nodes

    def __repr__(self):
        fields = [self.type,self.lid]+self.nodes
        return printCard(fields)


class Coord(object):
    def __init__(self,card):
        #self.type = card[0]
        self.cid  = card[1]
        self.dunno = card[2]
        #print "cid = ",self.cid

    def normalize(self,v):
        #print "v = ",v
        return v/norm(v)

    def __repr__(self):
        fields = [self.type,self.cid]
        return printCard(fields)

class CORD2R(Coord):  # working...
    type = 'CORD2R'
    def __init__(self,card=['CORD2R',0,0,  0.,0.,0.,  0.,0.,1., 1.,0.,0.]):
        Coord.__init__(self,card)
        
        #print card
        self.eo = array([card[3], card[4], card[5]] )
        self.ez = array([card[6], card[7], card[8]] )
        self.ex = array([card[9], card[10],card[11]])
        
        #print "eo = ",self.eo
        #print "ez = ",ez
        #print "ex = ",ex

        self.ez0 = self.normalize(self.ez-self.eo)
        self.ex0 = self.normalize(self.ex-self.eo)
        self.ey0 = cross(self.ez0,self.ex0)
        
        #print norm(self.ex)
        #print norm(self.ey)
        #print norm(self.ez)
        #print card
        #print str(self)
        #sys.exit()

        
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
        return printCard(fields)

class Material(object):
    def __init__(self,card):
        #self.type = card[0]
        self.mid  = card[1]
        
    def __repr__(self):
        fields = [self.type,self.mid]
        return printCard(fields)

class MAT1(Material):
    """MAT1     1      1.03+7  3.9615+6.3      .098"""
    type = 'MAT1'
    def __init__(self,card):
        Material.__init__(self,card)
        self.E     = card[2]
        self.G     = card[3]
        self.nu    = card[4]
        self.dunno = card[5]
        #print "card = ",card
        #print self
        #sys.exit()

    def __repr__(self):
        fields = [self.type,self.mid,self.mid,self.E,self.G,self.nu,self.dunno,]
        #print "hi...\n"
        return printCard(fields)

class Property(object):
    def __init__(self,card):
        #self.type = card[0]
        self.pid = card[1]
        self.mid = card[2]
        
    def __repr__(self):
        fields = [self.type,self.pid]
        return printCard(fields)

class PSHELL(Property):
    """PSHELL   41111   1      1.0000   1               1               0.02081"""
    type = 'PSHELL'
    def __init__(self,card):
        Property.__init__(self,card)
        self.dunno1 = card[3]
        self.dunno2 = card[4]
        self.dunno3 = card[5]
        self.dunno4 = card[6]
        self.dunno5 = card[7]
        self.thk    = card[8]
        #print "card = ",card
        #print self
        #sys.exit()

    def __repr__(self):
        fields = [self.type,self.pid,self.mid,self.dunno1,self.dunno2,self.dunno3,
                  self.dunno4,self.dunno5, self.thk]
        return printCard(fields)

class PSOLID(Property):
    """PSOLID   1       1       0"""
    type = 'PSOLID'
    def __init__(self,card):
        Property.__init__(self,card)
        self.dunno1 = card[3]
        #print "card = ",card
        #print self
        #sys.exit()

    def __repr__(self):
        fields = [self.type,self.pid,self.mid,self.dunno1]
        return printCard(fields)

class Element(object):
    def __init__(self,card):
        #displayCard(card)
        self.nodes = None
        self.eid = int(card[1])
        self.id = self.eid

        self.pid = None # CONM2s dont have a pid
        #self.nids = []
        pass

    def prepareNodeIDs(self,nids):
        self.nodes = []
        for nid in nids:
            self.nodes.append(int(nid))

    def Is(self,typeCheck):
        if self.type==typeCheck:
            return True
        return False

    def Centroid(self,nodes,debug=False):
        return None

    def getNodeIDs(self):
        nids = self.nodes
        #nids = []
        #for node in nodes:
        #    nids.append(node.nid)
        #print "nids[%s] = %s" %(self.eid,nids)
        return nids

    def getNodes(self,nodes):
        """
        returns nodes...???
        """
        #print "self.type = ",self.type
        nids = self.nodes
        #print "nids = ",self.nodes
        nNodes = len(self.nodes)
        nodesList = []
        for nid in range(nNodes):
            nodesList.append(nodes[nid])
        return nodesList

    #def Normal(self,a,b):
    #    """finds the unit normal vector of 2 vectors"""
    #    return Normal(a,b)

    def CentroidTriangle(self,n1,n2,n3,debug=False):
        if debug:
            print "n1=%s \nn2=%s \nn3=%s" %(n1,n2,n3)
        centroid = (n1+n2+n3)/3.
        return centroid

    #def Area(self,a,b):
    #    return 0.5*numpy.linalg.norm(numpy.cross(a,b))

    def __repr__(self):
        fields = [self.type,self.eid,self.pid]+self.nodes
        return printCard(fields)

class CTRIA3(Element):
    type = 'CTRIA3'
    def __init__(self,card):
        Element.__init__(self,card)
        self.pid = card[2]

        nids = card[3:6]
        self.prepareNodeIDs(nids)

        #print "self.xi = ",self.xi
        assert len(self.nodes)==3
        #raise

    def crossReference(self,nodes):
        pass

    def AreaCentroidNormal(self,nodes):
        return Triangle_AreaCentroidNormal(nodes)

    def Area(self,nodes):
        #(n0,n1,n2) = self.getNodes(nodes)
        (n0,n1,n2) = nodes
        a = n0-n1
        b = n0-n2
        area = Area(a,b)

    def Normal(self,nodes):
        #(n0,n1,n2) = self.getNodes(nodes)
        (n0,n1,n2) = nodes
        a = n0-n1
        b = n0-n2
        return Normal(a,b)

    def Centroid(self,nodes,debug=False):
        (n0,n1,n2) = self.getNodes(nodes)
        centroid = self.CentroidTriangle(n0,n1,n2,debug)
        return centroid

class CTETRA(Element):
#CTETRA   1       1       239     229     516     99      335     103
#         265     334     101     102

    type = 'CTETRA'
    def __init__(self,card):
        Element.__init__(self,card)
        self.pid = card[2]

        nids = card[3:13]
        self.prepareNodeIDs(nids)

class CQUAD4(Element):
    type = 'CQUAD4'
    def __init__(self,card):
        Element.__init__(self,card)
        self.pid = card[2]

        nids = card[3:7]
        self.prepareNodeIDs(nids)

        #print "self.xi = ",self.xi
        #print "nodes = ",self.nodes
        assert len(self.nodes)==4
        #for nid in nids:
        #    self.nodes.append(int(nid))

    def writeAsCTRIA3(self,newID):
        """
        triangle - 012
        triangle - 023
        """
        nodes1 = [self.nodes[0],self.nodes[1],self.nodes[2]]
        nodes2 = [self.nodes[0],self.nodes[2],self.nodes[3]]
        fields1 = ['CTRIA3',self.eid,self.mid]+nodes1
        fields2 = ['CTRIA3',newID,   self.mid]+nodes2
        return printCard(fields1)+printCard(fields2)
    
    def Normal(self,nodes):
        (n1,n2,n3,n4) = nodes
        #(n1,n2,n3,n4) = self.getNodes(nodes)
        a = n1-n3
        b = n2-n4
        return Normal(a,b)

    def AreaCentroidNormal(self,nodes):
        (area,centroid) = self.AreaCentroid(nodes)
        normal = self.Normal(nodes)
        return (area,centroid,normal)

    def AreaCentroid(self,nodes,debug=False):
        """
        1-----2
        |    /|
        | A1/ |
        |  /  |
        |/ A2 |
        4-----3
        
        centroid
           c = sum(ci*Ai)/sum(A)
           where:
             c=centroid
             A=area
        """
        if debug:
            print "nodes = ",nodes
        (n1,n2,n3,n4) = nodes
        #(n1,n2,n3,n4) = self.getNodes(nodes)
        a = n1-n2
        b = n2-n4
        area1 = Area(a,b)
        c1 = self.CentroidTriangle(n1,n2,n4)

        a = n2-n4
        b = n2-n3
        area2 = Area(a,b)
        c2 = self.CentroidTriangle(n2,n3,n4)
        
        area = area1+area2
        centroid = (c1*area1+c2*area2)/area
        if debug:
            print "c1=%s \n c2=%s \n a1=%s a2=%s" %(c1,c2,area1,area2)
            print "type(centroid=%s centroid=%s \n" %(type(centroid),centroid)
        return(area,centroid)
    ###

    def Centroid(self,nodes,debug=False):
        (area,centroid) = self.AreaCentroid(nodes,debug)
        return centroid

    def Area(self,nodes):
        """
        Area = 1/2*|a cross b|
        where a and b are the quad's cross node point vectors
        """
        (n1,n2,n3,n4) = self.getNodes(nodes)
        a = n1-n3
        b = n2-n4
        area = Area(a,b)
        return area
    ###

class CONM2(Element):
    type = 'CONM2'
    # 'CONM2    501274  11064          132.274'
    def __init__(self,card):
        Element.__init__(self,card)
        #self.nids  = [ card[1] ]
        #del self.nids
        self.dunno = card[2]
        self.blank = card[3]
        self.mass  = card[4]
        
        #print "nids       = ",self.nids
        #print 'self.dunno = ',self.dunno
        #print 'self.blank = ',self.blank
        #print "mass       = ",self.mass
        #print "card       = ",card
        #print str(self)
        #sys.exit()
    
    def __repr__(self):
        fields = [self.type,self.eid,self.dunno,self.blank,self.mass]
        #fields = [self.type,self.eid,self.blank,self.mass]
        return printCard(fields)

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

