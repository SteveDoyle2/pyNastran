from __future__ import print_function
import sys

from numpy import pi,matrix,zeros,ones,array,transpose,dot
from numpy.linalg import norm
from pyNastran.general.generalMath import printMatrix

from pyNastran.bdf.errors import *
from pyNastran.bdf.cards.baseCard import Element,Mid


class LineElement(Element):
    def __init__(self,card,data):
        Element.__init__(self,card,data)

    def C(self):
        """torsional constant"""
        return self.pid.C()

    def Area(self):
        """returns the area of the element face"""
        raise NotImplementedMethodError('implement self.Area() for %s' %(self.type))
    
    def E(self):
        """returns the Young's Modulus  \f$ E \f$"""
        return self.pid.mid.E()
     
    def G(self):
        """returns the Shear Modulus   \f$ G \f$"""
        return self.pid.mid.G()

    def J(self):
        """returns the Polar Moment of Inertia.   \f$ J \f$"""
        return self.pid.J()

    def I11(self):
        """returns the Moment of Inertia.   \f$ I_{11} \f$"""
        return self.pid.I11()

    def I22(self):
        """returns the Moment of Inertia.   \f$ I_{22} \f$"""
        return self.pid.I22()

    def I12(self):
        """returns the Moment of Inertia.   \f$ I_{12} \f$"""
        return self.pid.I12()

    def Nu(self):
        """returns Poisson's Ratio  \f$ \nu \f$"""
        return self.pid.mid.nu
    
    def Rho(self):
        """returns the material density  \f$ \rho \f$"""
        #print str(self.pid),type(self.pid)
        #raise NotImplementedMethodError('implement self.Rho() for %s' %(self.type))
        return self.pid.mid.rho

    def Nsm(self):
        """Placeholder method for the non-structural mass"""
        raise NotImplementedMethodError('implement self.Area() for %s' %(self.type))

    def MassPerLength(self):
        """Returns the mass per unit length"""
        return self.pid.MassPerLength()

    #def MassPerLength(self):
        """Returns the mass per unit length"""
        #massPerLength = self.Rho()*self.Area() + self.Nsm()

    def Mass(self):
        """
        returns the mass of a bar/beam/rod element
        
        \f[ \large  mass = \left( \rho A + nsm \right) L  \f]
        """
        L = self.Length()
        #print "type = ",self.type
        
        #try:
            #print "mat  = ",self.mid.type
        #except:
            #print "no material"
        #print "rho",self.Rho()
        #print "area",self.Area()
        #print "nsm",self.Nsm()
        #print "L",L
        mass = (self.Rho()*self.Area() + self.Nsm())*L
        #print "mass",mass
        #print "---"
        return mass

    def crossReference(self,mesh):
        self.nodes = mesh.Nodes(self.nodes)
        self.pid   = mesh.Property(self.pid)

    def Length_noXref(self,n1=None,n2=None):
        """
        Returns the length of a bar/rod/beam element
        \f[ \large \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  } \f]
        @param self the object pointer
        @param n1 a Node object (default=None)
        @param n2 a Node object (default=None)
        @note
            if n1 AND n2 are both none (the default), then the model must
            be cross-referenced already
        """
        #print self.type
        L = norm(n1.Position()-n2.Position())
        return L

    def Length(self):
        """
        Returns the length of a bar/rod/beam element
        \f[ \large \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  } \f]
        @param self the object pointer
        """
        #print self.type
        return self.Length_noXref(self.nodes[1],self.nodes[0])

    def k_Axial(self):
        """
        Returns the axial stiffness matrix.

        \f[ \large   k_{Axial} = \frac{AE}{2L} 
          \left[
          \begin{array}{cc}
              1 & -1 \\
             -1 &  1
          \end{array} \right]
        \f]
        """
        L = self.Length()
        E = self.E()
        A = self.Area()
        kMag = A*E/(2*L)
        M = Matrix(ones(1,1))
        M[0,1] = M[1,0] = -1
        return M

    def k_Torsion(self):  # not done
        """
        Returns the torsional stiffness matrix.

        \f[ \large   k_{Axial} = \frac{L}{GJ} 
          \left[
          \begin{array}{cc}
              1 & -1 \\
             -1 &  1
          \end{array} \right]
        \f]
        @warning formula not verified
        """
        L = self.Length()
        G = self.G()
        J = self.J()
        #A = self.Area()
        #kMag = A*E/(2*L)
        kMag = L/(G*J)
        M = Matrix(ones(1,1))
        M[0,1] = M[1,0] = -1
        return M

    def k_Bending(self):
        """
        Returns the bending stiffness matrix.

        \f[ \large  k_{Bending} = \frac{EI}{L^3} 
          \left[ 
          \begin{array}{cccc}
             12 &  6L   & -12 &  6L    \\
             6L &  4L^2 & -6L &  2L^2  \\
            -12 & -6L   &  12 & -6L    \\
             6L &  2L^2 & -6L &  4L^2
          \end{array} \right]
        \f] 
        """
        L = self.Length()
        E = self.E()
        I = self.I()
        kMag = E*I/LLL
        LL  = L*L
        LLL = L*LL
        sL = 6*L
        tLL = 2*LL
        fLL = 4*LL
        
        M = Matrix(zeros(4,4))
        matrix([[12.,  sL, -12,  sL],
                [sL,  fLL, -sL, tLL],
                [-12, -sL, 12., -sL],
                [sL,  tLL, -sL, fLL]])
        #M[1,0] = sL
        #M[2,0] = -12.
        #M[3,0] = sL

        #M[2,4] =  -sL
        #M[1,1] = M[3,3] = fLL
        
        return M

class CROD(LineElement):
    type = 'CROD'
    def __init__(self,card=None,data=None):
        LineElement.__init__(self,card,data)
        if card:
            self.eid = int(card.field(1))
            self.pid = int(card.field(2,self.eid))
            nids = card.fields(3,5)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        ###
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==2

    def Centroid(self):
        return (self.nodes[0].Position()+self.nodes[1].Position())/2.

    def Mid(self):
        return self.pid.Mid()

    def Area(self):
        return self.pid.A

    def Nsm(self):
        return self.pid.nsm

    def MassPerLength(self):
        #print "self.type = ",self.type
        massPerLength = self.pid.mid.rho*self.pid.A + self.pid.nsm

    def Rmatrix(self,model,is3D):
        """
        where   \f$ [R]_{ij} \f$ is the tranformation matrix
        \f[ \large  [R]_{ij} \left[ 
          \begin{array}{ccc}
              g_x \dot e_x & g_x \dot e_y &  g_x \dot e_z    \\
              g_y \dot e_x & g_y \dot e_y &  g_y \dot e_z    \\
              g_z \dot e_x & g_z \dot e_y &  g_z \dot e_z
          \end{array} \right]
        \f] 
        """
        (n1,n2) = self.nodeIDs()
        p1 = model.Node(n1).xyz
        p2 = model.Node(n2).xyz
        v1 = p2-p1
        v1 = v1/norm(v1)
        (l,m,n) = v1
        
        v1x = array([v1[0],0.,0.])
        v1y = array([0.,v1[1],0.])
        v1z = array([0.,0.,v1[2]])
        
        g1x = array([1.,0.,0.])
        g1y = array([0.,1.,0.])
        g1z = array([0.,0.,1.])

        if is3D:
            R = matrix([  #global rod
                        [dot(v1x,g1x),dot(v1y,g1x),dot(v1z,g1x)],
                        [dot(v1x,g1y),dot(v1y,g1y),dot(v1z,g1y)],
                        [dot(v1x,g1z),dot(v1y,g1z),dot(v1z,g1z)],
                      ]) # rod
            #R = matrix([
            #            [],
            #            [],
            #            [],
            #          ])
                      
        else:
            R = matrix([ # there can be no z component
                        [dot(v1x,g1x),dot(v1y,g1x)],
                        [dot(v1x,g1y),dot(v1y,g1y)],
                      ]) # rod
            #R = matrix([ # there can be no z component
            #            [dot(v1x,g1x),dot(v1y,g1x)],
            #            [dot(v1x,g1y),dot(v1y,g1y)],
            #          ]) # rod
        return R
        
    def Lambda(self,model):
        """
        2d  [l,m,0,0]
            [0,0,l,m]
        
        3d  [l,m,n,0,0,0]
            [0,0,0,l,m,n]
        """
        is3D = False
        #R = self.Rmatrix(model,is3D)
        
        (n1,n2) = self.nodeIDs()
        p1 = model.Node(n1).Position()
        p2 = model.Node(n2).Position()
        v1 = p2-p1
        print(v1)
        v1 = v1/norm(v1)
        (l,m,n) = v1
        if is3D:
            Lambda = matrix(zeros((2,6),'d')) # 3D
        else:
            Lambda = matrix(zeros((2,4),'d'))

        #print("R = \n",R)
        Lambda[0,0] = Lambda[1,2] = l
        Lambda[0,1] = Lambda[1,3] = m

        if is3D:
            Lambda[0,2] = Lambda[1,5] = n # 3D
        print("Lambda = \n",Lambda)
        return Lambda

    def Stiffness1(self,model):
        nIJV = [(nodes[0],1),(nodes[0],2),(nodes[1],1),]

    def Stiffness(self,model):
        Lambda = self.Lambda(model)
        #print "Lambda = \n",Lambda
        
        k = self.Stiffness1D(model) #/250000.
        #print R
        #print(k)
        K = dot(dot(transpose(Lambda),k),Lambda)
        #print(K)
        print("K[%s] = \n%s\n" %(self.eid,K))
        #sys.exit('stopping in ROD stiffness...elementsBars.py')

        nodes = self.nodeIDs()
        nIJV = [(nodes[0],1),(nodes[0],2),(nodes[1],1),(nodes[1],2),]

        return K,nIJV

    def displacementStress(self,model,q,dofs):
        (n1,n2) = self.nodeIDs()
        Lambda = self.Lambda(model)
        n11 = dofs[(n1,1)]
        n12 = dofs[(n1,2)]
        n21 = dofs[(n2,1)]
        n22 = dofs[(n2,2)]
        
        q2 = array( [q[n11],q[n12],q[n21],q[n22]] )
        print("q[%s] = %s" %(self.eid,q2))
        #print "Lambda = \n",Lambda
        
        #print "Lsize = ",Lambda.shape
        #print "qsize = ",q.shape
        u = dot(array(Lambda),q2)
        #L = self.Length()
        EL = self.E()/self.Length()
        
        stressX = -EL*u[0]+EL*u[1]
        print("stressX = %s [psi]" %(stressX))

    def Stiffness1D(self,model):
        """
        @todo remove this method after making sure it still works
        """
        #L = norm(r)
        (n1,n2) = self.nodeIDs()
        node1 = model.Node(n1)
        node2 = model.Node(n2)

        p1 = model.Node(n1).xyz
        p2 = model.Node(n2).xyz
        #print "node1 = ",node1
        #print "node2 = ",node2
        L = norm(p1-p2)
        
        if L==0.0:
            msg = 'invalid CROD length=0.0\n%s' %(self.__repr__())
            raise StiffnessMatrixError(msg)
        
        A = self.Area()
        #mat = self.mid
        E = self.E()
        print("A = ",A)
        print("E = ",E)
        print("L = ",L)
        f = A*E/L
        #f = 250000.
        K = 1*matrix([[1.,-1.],[-1.,1.]]) # rod
        
        print("A=%g E=%g L=%g  AE/L=%g" %(A,E,L,A*E/L))
        #print "K = \n",K
        return K

    def rawFields(self):
        fields = ['CROD',self.eid,self.Pid()]+self.nodeIDs()
        return fields

    def reprFields(self):
        return self.rawFields()

class CVISC(CROD):
    type = 'CVISC'
    def __init__(self,card=None,data=None):
        LineElement.__init__(self,card,data)
        if card:
            self.eid = int(card.field(1))
            self.pid = int(card.field(2,self.eid))
            nids = card.fields(3,5)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        ###
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==2
    ###
    def rawFields(self):  # not done...
        fields = ['CVISC',self.eid,self.Pid()]+self.nodeIDs()
        return fields

    def reprFields(self):
        return self.rawFields()
###

class CTUBE(CROD):
    type = 'CTUBE'
    def __init__(self,card=None,data=None):
        CROD.__init__(self,card,data)
    ###
    def Area(self):
        return self.pid.Area()

    def rawFields(self):
        fields = ['CTUBE',self.eid,self.Pid()]+self.nodeIDs()
        return fields
###

class CONROD(CROD):
    type = 'CONROD'
    def __init__(self,card=None,data=None):
        LineElement.__init__(self,card,data)
        if card:
            self.eid  = int(card.field(1))
            #print "self.eid = ",self.eid
            nids = card.fields(2,4)

            self.mid = int(card.field(4))
            self.A   = float(card.field(5))
            self.j   = float(card.field(6,0.0))
            self.c   = float(card.field(7,0.0))
            self.nsm = float(card.field(8,0.0))
        else:
            self.eid = data[0]
            nids = data[1:3]
            self.mid = data[3]
            self.A   = data[4]
            self.j   = data[5]
            self.c   = data[6]
            self.nsm = data[7]
        ###
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==2
        #print self.nodes

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.mid   = model.Material(self.mid)

    def Mid(self):
        return Mid(self)

    def Pid(self):
        return None

    def MassPerLength(self):
        massPerLength = self.mid.rho*self.A + self.nsm

    def C(self):
        """torsional constant"""
        return self.c

    def Area(self):
        return self.A
    
    def J(self):
        """returns the Polar Moment of Inertia.   \f$ J \f$"""
        return self.j

    def Nsm(self):
        """Placeholder method for the non-structural mass"""
        return self.nsm

    def E(self):
        """returns the Young's Modulus  \f$ E \f$"""
        return self.mid.E()
     
    def G(self):
        """returns the Shear Modulus   \f$ G \f$"""
        return self.mid.G()

    def I12(self):
        """the I12 for a rod is 0"""
        return 0.

    def I22(self):
        return self.I11()

    def I11(self):
        """returns the Moment of Inertia.   \f$ I \f$"""
        r = self.Radius()
        #A=piD2/4
        #I = piR^4/4 = pi*D^4/16*4 = pi*D^4/64 = A*D^2/16 = A^2/pi
        #  = A*R^2/4 = A^2/pi/4
        I = self.A**2/(4*pi)
        return I

    def Radius(self):
        """returns the Radius of the Rod.   \f$ R \f$"""
        R = (self.A/pi)**0.5
        return R

    def Nu(self):
        """returns Poisson's Ratio  \f$ \nu \f$"""
        return self.mid.nu
    
    def Rho(self):
        """returns the material density  \f$ \rho \f$"""
        return self.mid.rho

    def writeCodeAster(self):
        msg += "    POUTRE=_F(GROUP_MA='CONROD_%s',\n" %(self.eid)
        msg += "              SECTION='CERCLE',  # circular section\n"
        if thickness:
            msg += "              CARA=('R','EP'),   # radius, thickness\n"
            msg += "              VALE=(%g,%g),\n"  %(self.Radius(),self.Thickness())
        else:
            msg += "              CARA=('R')   # radius\n"
            msg += "              VALE=(%g),\n"  %(self.Radius())
        ###
        return msg

    def rawFields(self):
        fields = ['CONROD',self.eid]+self.nodeIDs()+[self.Mid(),self.A,self.j,self.c,self.nsm]
        return fields

    def reprFields(self):
        j   = self.setBlankIfDefault(self.j,  0.0)
        c   = self.setBlankIfDefault(self.c,  0.0)
        nsm = self.setBlankIfDefault(self.nsm,0.0)
        #print "nodes = ",self.nodeIDs()
        #print "mid   = ",self.Mid()
        #print "eid   = ",self.eid
        fields = ['CONROD',self.eid]+self.nodeIDs()+[self.Mid(),self.A,j,c,nsm]
        return fields

class CBAR(LineElement):
    """
    CBAR EID PID GA GB X1 X2 X3 OFFT
    PA PB W1A W2A W3A W1B W2B W3B
    or
    CBAR EID PID GA GB G0 OFFT
    PA PB W1A W2A W3A W1B W2B W3B

    CBAR       22062       4   21798   21799   0.0+0   0.0+0     -1.
                               0.0+0   0.0+0     -9.   0.0+0   0.0+0     -9.
    """
    type = 'CBAR'
    asterType = 'CBAR'
    def __init__(self,card=None,data=None):
        LineElement.__init__(self,card,data)
        if card:
            self.eid = int(card.field(1))
            self.pid = int(card.field(2,self.eid))
            self.ga  = int(card.field(3))
            self.gb  = int(card.field(4))
            self.initX_G0(card)

            self.offt = card.field(8,'GGG')
            #print 'self.offt = |%s|' %(self.offt)

            self.pa = card.field(9,0)
            self.pb = card.field(10,0)

            self.w1a = float(card.field(11,0.0))
            self.w2a = float(card.field(12,0.0))
            self.w3a = float(card.field(13,0.0))

            self.w1b = float(card.field(14,0.0))
            self.w2b = float(card.field(15,0.0))
            self.w3b = float(card.field(16,0.0))
        else: ## @todo verify
            #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,g0]]
            #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,x1,x2,x3]]

            main = data[0]

            flag = data[1][0]
            if flag in [0,1]:
                self.g0 = None
                self.x1 = data[1][1]
                self.x2 = data[1][2]
                self.x3 = data[1][3]
            else:
                self.g0 = data[1][1]
                self.x1 = None
                self.x2 = None
                self.x3 = None

            self.eid  = main[0]
            self.pid  = main[1]
            self.ga   = main[2]
            self.gb   = main[3]
            #self.offt = str(data[4]) # GGG
            self.offt = 'GGG'
            self.pa   = main[4]
            self.pb   = main[5]

            self.w1a = main[6]
            self.w2a = main[7]
            self.w3a = main[8]

            self.w1b = main[9]
            self.w2b = main[10]
            self.w3b = main[11]
        ###
        assert isinstance(self.offt,str),'invalid offt expected a string of length 3 offt=|%s|' %(self.offt)
        assert self.offt[0] in ['G','B','O'],'invalid offt parameter of %s...offt=%s' %(self.type,self.offt)
        assert self.offt[1] in ['G','B','O'],'invalid offt parameter of %s...offt=%s' %(self.type,self.offt)
        assert self.offt[2] in ['G','B','O'],'invalid offt parameter of %s...offt=%s' %(self.type,self.offt)

    def Mid(self):
        return self.pid.Mid()

    def Area(self):
        A = self.pid.Area()
        assert isinstance(A,float)
        return A

    def Length(self):
        L = self.Length_noXref(self.ga,self.gb)
        assert isinstance(L,float)
        return L

    def Nsm(self):
        nsm = self.pid.Nsm()
        assert isinstance(nsm,float)
        return nsm

    def Centroid(self):
        return (self.ga.Position()+self.gb.Position())/2.

    def initX_G0(self,card):
        field5 = card.field(5)
        if isinstance(field5,int):
            self.g0 = field5
            self.x1 = None
            self.x2 = None
            self.x3 = None
        elif isinstance(field5,float):
            self.g0 = None
            self.x1 = float(card.field(5,0.0))
            self.x2 = float(card.field(6,0.0))
            self.x3 = float(card.field(7,0.0))
        else:
            #msg = 'field5 on %s is the wrong type...id=%s field5=%s type=%s' %(self.type,self.eid,field5,type(field5))
            #raise InvalidFieldError(msg)
            self.g0 = None
            self.x1 = 0.
            self.x2 = 0.
            self.x3 = 0.
        #if self.eid==14100238:
            #print "g0=%s x1=%s x2=%s x3=%s" %(self.g0,self.x1,self.x2,self.x3)

    def crossReference(self,mesh):
        """
        set g0-ga to x1,x2,x3
        """
        #if self.g0:
        #    v = nodes[self.g0].Position()-nodes[self.ga].Position()
        #    self.x1 = v[0]
        #    self.x2 = v[1]
        #    self.x3 = v[2]
        ###
        self.ga  = mesh.Node(self.ga)
        self.gb  = mesh.Node(self.gb)
        self.pid = mesh.Property(self.pid)

    #def updateNodes(self,nodes):
    #    """@todo maybe improve"""
    #    self.crossReference(self,nodes)

    def Ga(self):
        if isinstance(self.ga,int):
            return self.ga
        else:
            return self.ga.nid
        ###

    def Gb(self):
        if isinstance(self.gb,int):
            return self.gb
        else:
            return self.gb.nid
        ###

    def getX_G0_defaults(self):
        if self.g0:
            return (self.g0,None,None)
        else:
            #x1 = self.setBlankIfDefault(self.x1,0.0)
            #x2 = self.setBlankIfDefault(self.x2,0.0)
            #x3 = self.setBlankIfDefault(self.x3,0.0)
            return (self.x1,self.x2,self.x3)
        ###

    def nodeIDs(self):
        return [self.Ga(),self.Gb()]

    def rawFields(self):
        """@todo not perfectly accurate"""
        (x1,x2,x3) = self.getX_G0_defaults()
        offt = self.setBlankIfDefault(self.offt,'GGG')
        fields = ['CBAR',self.eid,self.Pid(),self.Ga(),self.Gb(),x1,x2,x3,offt,
                         self.pa,self.pb,self.w1a,self.w2a,self.w3a,self.w1b,self.w2b,self.w3b]

        return fields

    def reprFields(self):
        pa = self.setBlankIfDefault(self.pa,0)
        pb = self.setBlankIfDefault(self.pb,0)

        w1a = self.setBlankIfDefault(self.w1a,0.0)
        w2a = self.setBlankIfDefault(self.w2a,0.0)
        w3a = self.setBlankIfDefault(self.w3a,0.0)
        w1b = self.setBlankIfDefault(self.w1b,0.0)
        w2b = self.setBlankIfDefault(self.w2b,0.0)
        w3b = self.setBlankIfDefault(self.w3b,0.0)
        (x1,x2,x3) = self.getX_G0_defaults()
        offt = self.setBlankIfDefault(self.offt,'GGG')
        fields = ['CBAR',self.eid,self.Pid(),self.Ga(),self.Gb(),x1,x2,x3,offt,
                         pa,pb,w1a,w2a,w3a,w1b,w2b,w3b]

        return fields

class CBEAM3(CBAR):
    """
    Defines a three-node beam element
    """
    type = 'CBEAM3'
    def __init__(self,card=None,data=None):
        LineElement.__init__(self,card,data)
        if card:
            self.eid = int(card.field(1))
            self.eid = int(card.field(1))
            self.pid = int(card.field(2,self.eid))
            self.ga  = int(card.field(3))
            self.gb  = int(card.field(4))
            self.gc  = int(card.field(5))

            self.initX_G0(card)

            self.w1a = float(card.field(9,0.0))
            self.w2a = float(card.field(10,0.0))
            self.w3a = float(card.field(11,0.0))

            self.w1b = float(card.field(12,0.0))
            self.w2b = float(card.field(13,0.0))
            self.w3b = float(card.field(14,0.0))

            self.w1c = float(card.field(15,0.0))
            self.w2c = float(card.field(16,0.0))
            self.w3c = float(card.field(17,0.0))

            self.twa = card.field(18,0.)
            self.twb = card.field(19,0.)
            self.twc = card.field(20,0.)

            self.sa = card.field(21)
            self.sb = card.field(22)
            self.sc = card.field(23)
        else:
            raise NotImplementedError(data)
        ###

    def crossReference(self,model):
        self.ga  = model.Node(self.ga)
        self.gb  = model.Node(self.gb)
        self.gc  = model.Node(self.gc)
        self.pid = model.Property(self.pid)

    def rawFields(self):
        (x1,x2,x3) = self.getX_G0_defaults()
        (ga,gb,gc) = self.nodeIDs()
        fields = ['CBEAM3',self.eid,self.Pid(),ga,gb,gc,x1,x2,x3,
                  self.w1a,self.w2a,self.w3a, self.w1b,self.w2b,self.w3b, self.w1c,self.w2c,self.w3c,
                  self.twa,self.twb,self.twc,self.sa,self.sb,self.sc]
        return fields

    def reprFields(self):
        w1a = self.setBlankIfDefault(self.w1a,0.0)
        w2a = self.setBlankIfDefault(self.w2a,0.0)
        w3a = self.setBlankIfDefault(self.w3a,0.0)
        w1b = self.setBlankIfDefault(self.w1b,0.0)
        w2b = self.setBlankIfDefault(self.w2b,0.0)
        w3b = self.setBlankIfDefault(self.w3b,0.0)
        w1c = self.setBlankIfDefault(self.w1c,0.0)
        w2c = self.setBlankIfDefault(self.w2c,0.0)
        w3c = self.setBlankIfDefault(self.w3c,0.0)

        twa = self.setBlankIfDefault(self.twa,0.0)
        twb = self.setBlankIfDefault(self.twb,0.0)
        twc = self.setBlankIfDefault(self.twc,0.0)

        (x1,x2,x3) = self.getX_G0_defaults()
        (ga,gb,gc) = self.nodeIDs()
        fields = ['CBEAM',self.eid,self.Pid(),ga,gb,x1,x2,x3,offt,
                  w1a,w2a,w3a, w1b,w2b,w3b, w1c,w2c,w3c,
                  twa,twb,twc,self.sa,self.sb,self.sc]
        return fields


class CBEAM(CBAR):
    """
    CBEAM EID PID GA GB X1 X2 X3 OFFT/BIT
    PA PB W1A W2A W3A W1B W2B W3B
    SA SB
    or
    CBEAM EID PID GA GB G0 - - OFFT/BIT
    PA PB W1A W2A W3A W1B W2B W3B
    SA SB
    """
    type = 'CBEAM'
    def __init__(self,card=None,data=None):
        LineElement.__init__(self,card,data)
        if card:
            self.eid = int(card.field(1))
            self.pid = int(card.field(2,self.eid))
            self.ga  = int(card.field(3))
            self.gb  = int(card.field(4))

            self.initX_G0(card)
            self.initOfftBit(card)
            self.pa = card.field(9)
            self.pb = card.field(10)

            self.w1a = float(card.field(11,0.0))
            self.w2a = float(card.field(12,0.0))
            self.w3a = float(card.field(13,0.0))

            self.w1b = float(card.field(14,0.0))
            self.w2b = float(card.field(15,0.0))
            self.w3b = float(card.field(16,0.0))

            self.sa = card.field(17)
            self.sb = card.field(18)

        else: ## @todo verify
            #data = [[eid,pid,ga,gb,sa,sb, pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,g0]]
            #data = [[eid,pid,ga,gb,sa,sb, pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,x1,x2,x3]]

            main = data[0]

            flag = data[1][0]
            if flag in [0,1]:
                self.g0 = None
                self.x1 = data[1][1]
                self.x2 = data[1][2]
                self.x3 = data[1][3]
            else:
                self.g0 = data[1][1]
                self.x1 = None
                self.x2 = None
                self.x3 = None

            self.eid = main[0]
            self.pid = main[1]
            self.ga  = main[2]
            self.gb  = main[3]
            self.sa  = main[4]
            self.sb  = main[5]

            self.isOfft = True ## @todo is this correct???
            #self.offt = str(data[6]) # GGG
            self.offt = 'GGG' ## @todo is this correct???

            self.pa   = main[6]
            self.pb   = main[7]

            self.w1a = main[8]
            self.w2a = main[9]
            self.w3a = main[10]

            self.w1b = main[11]
            self.w2b = main[12]
            self.w3b = main[13]
        ###

    def initOfftBit(self,card):
        field8 = card.field(8)
        if isinstance(field8,float):
            self.isOfft = False
            self.offt = None
            self.bit = field8
        elif field8 is None:
            self.isOfft = True
            self.offt = 'GGG' # default
            self.bit = None
        elif isinstance(field8,str):
            self.isOfft = True
            self.bit = None
            self.offt = field8
            #print "self.offt = ",self.offt
            assert self.offt[0] in ['G','B','O'],'invalid offt parameter of CBEAM...offt=%s' %(self.offt)
            assert self.offt[1] in ['G','B','O'],'invalid offt parameter of CBEAM...offt=%s' %(self.offt)
            assert self.offt[2] in ['G','B','O'],'invalid offt parameter of CBEAM...offt=%s' %(self.offt)
        else:
            msg = 'field8 on %s card is not a string(offt) or bit (float)...field8=%s\n' %(self.type,field8)
            raise CardInstantiationError(msg)
        ###

    def Mid(self):
        return self.pid.Mid()

    def Area(self):
        return self.pid.Area()

    def Nsm(self):
        #print "CBEAM pid = ",str(self.pid)
        return self.pid.Nsm()

    def getOfft_Bit_defaults(self):
        if self.isOfft:
            field8 = self.offt
        else:
            field8 = self.setBlankIfDefault(self.bit,0.0)
        return field8

    def crossReference(self,model):
        self.ga  = model.Node(self.ga)
        self.gb  = model.Node(self.gb)
        self.pid = model.Property(self.pid)

    def Stiffness(self,bdf,r,A,E,I):
        """
        from makeTruss???
        """
        Ke = matrix( zeros((6,6),'d'))
        L = r
        AE = A*E
        EI = E*I

        if 1:
            Ke[0,0] =  AE/L
            Ke[3,0] = -AE/L
            Ke[0,3] = -AE/L
            Ke[3,3] =  AE/L

            Ke[1,1] = 12*EI/L**3
            Ke[1,2] =  6*EI/L**2
            Ke[2,1] = Ke[1,2]  # 6*EI/L**2
            Ke[2,2] =  4*EI/L

            Ke[1,4] = -Ke[1,1] # -12*EI/L**3
            Ke[1,5] =  Ke[1,2] #  6*EI/L**2
            Ke[2,4] = -Ke[1,2] # -6*EI/L**2
            Ke[2,5] = 2*EI/L

            Ke[4,1] = -Ke[1,4] # -12*EI/L**3
            Ke[4,2] =  Ke[2,4] #  -6*EI/L**2
            Ke[5,1] =  Ke[1,2] #   6*EI/L**2
            Ke[5,2] =  Ke[2,5] #   2*EI/L

            Ke[4,4] = Ke[1,1] #  12*EI/L**3
            Ke[4,5] = Ke[2,4] #  -6*EI/L**2
            Ke[5,4] = Ke[2,4] #  -6*EI/L**2
            Ke[5,5] = Ke[2,2] #   4*EI/L
        else:
            Ke[0,0] =  AE
            Ke[3,0] = -AE
            Ke[0,3] = -AE
            Ke[3,3] =  AE

            Ke[1,1] = 12*EI/L**2
            Ke[1,2] =  6*EI/L
            Ke[2,1] = Ke[1,2]  # 6*EI/L**2
            Ke[2,2] =  4*EI

            Ke[1,4] = -Ke[1,1] # -12*EI/L**3
            Ke[1,5] =  Ke[1,2] #  6*EI/L**2
            Ke[2,4] = -Ke[1,2] # -6*EI/L**2
            Ke[2,5] = 2*EI

            Ke[4,1] = -Ke[1,4] # -12*EI/L**3
            Ke[4,2] =  Ke[2,4] #  -6*EI/L**2
            Ke[5,1] =  Ke[1,2] #   6*EI/L**2
            Ke[5,2] =  Ke[2,5] #   2*EI/L

            Ke[4,4] = Ke[1,1] #  12*EI/L**3
            Ke[4,5] = Ke[2,4] #  -6*EI/L**2
            Ke[5,4] = Ke[2,4] #  -6*EI/L**2
            Ke[5,5] = Ke[2,2] #   4*EI/L

            Ke = Ke/L
        ###
        return Ke

    def rawFields(self):
        (x1,x2,x3) = self.getX_G0_defaults()
        offt = self.getOfft_Bit_defaults()
        ga,gb = self.nodeIDs()
        fields = ['CBEAM',self.eid,self.Pid(),ga,gb,x1,x2,x3,offt,
                  self.pa,self.pb,self.w1a,self.w2a,self.w3a,self.w1b,self.w2b,self.w3b,
                  self.sa,self.sb]
        return fields

    def reprFields(self):
        w1a = self.setBlankIfDefault(self.w1a,0.0)
        w2a = self.setBlankIfDefault(self.w2a,0.0)
        w3a = self.setBlankIfDefault(self.w3a,0.0)
        w1b = self.setBlankIfDefault(self.w1b,0.0)
        w2b = self.setBlankIfDefault(self.w2b,0.0)
        w3b = self.setBlankIfDefault(self.w3b,0.0)
        (x1,x2,x3) = self.getX_G0_defaults()
        offt = self.getOfft_Bit_defaults()
        ga,gb = self.nodeIDs()
        fields = ['CBEAM',self.eid,self.Pid(),ga,gb,x1,x2,x3,offt,
                  self.pa,self.pb,w1a,w2a,w3a,w1b,w2b,w3b,
                  self.sa,self.sb]
        return fields

class CBEND(LineElement):
    type = 'CBEND'
    def __init__(self,card=None,data=None):
        LineElement.__init__(self,card,data)
        self.eid  = int(card.field(1))
        self.pid  = card.field(2,self.eid)
        self.ga = card.field(3)
        self.gb = card.field(4)
        x1Go = card.field(5)
        #self.initX_G0(card)
        if isinstance(x1Go,int):
            self.g0 = x1Go
            self.x1 = None
            self.x2 = None
            self.x3 = None
        elif isinstance(x1Go,float):
            self.g0 = None
            self.x1 = x1Go
            self.x2 = card.field(6)
            self.x3 = card.field(7)
        else:
            raise Exception('invalid x1Go=|%s| on CBEND' %(x1Go))
        self.geom = card.field(8)
        assert self.geom in [1,2,3,4],'geom is invalid geom=|%s|' %(self.geom)

    def Area(self):
        return self.pid.Area()

    def rawFields(self):
        (x1,x2,x3) = self.getX_G0_defaults()
        fields = ['CBEND',self.eid,self.Pid(),self.Ga(),self.Gb(),x1,x2,x3,self.geom]
        return fields
    
    def reprFields(self):
        return self.rawFields()
###

