from numpy import matrix,zeros,ones
from numpy.linalg import norm

class CardInstantiationError(RuntimeError):
    pass

from elements import Element

class LineElement(Element):
    def __init__(self,card):
        Element.__init__(self,card)

    def C(self):
        """torsional constant"""
        return self.pid.C()

    def area(self):
        raise Exception('implement self.area() for %s' %(self.type))
    
    def E(self):
        return self.pid.mid.E
     
    def G(self):
        return self.pid.mid.G

    def nu(self):
        return self.pid.mid.nu
    
    def rho(self):
        #print str(self.pid),type(self.pid)
        return self.pid.mid.rho

    def nsm(self):
        raise Exception('implement self.area() for %s' %(self.type))

    def mass(self):
        L = self.length()
        #print "rho",self.rho()
        #print "area",self.area()
        #print "nsm",self.nsm()
        #print "L",L
        mass = (self.rho()*self.area() + self.nsm())*L
        #print "mass",mass
        #print "---"
        return mass

    def crossReference(self,mesh):
        self.nids = mesh.Nodes(self.nodes)
        self.pid  = mesh.Property(self.pid)

    def length(self):
        #print self.type
        L = norm(self.nids[1].Position()-self.nids[0].Position())
        return L

    def length2(self,n1,n2):    # length refers to an alternate method
        #print self.type
        L = norm(n1.Position()-n2.Position())
        return L

    def k_Axial(self):
        L = self.length()
        E = self.E()
        A = self.Area()
        kMag = A*E/(2*L)
        M = Matrix(ones(1,1))
        M[0,1] = M[1,0] = -1
        return M

    def k_Bending(self):
        """
        kFac = EI/L^3
         k/kFac = 
        
         12   6L    -12   6L
         6L   4L^2  -6L   2L^2
        -12  -6L     12  -6L
         6L   2L^2  -6L   4L^2
        """
        L = self.length()
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
        M[1,0] = sL
        M[2,0] = -12.
        M[3,0] = sL

        M[2,4] =  -sL
        M[1,1] = M[3,3] = fLL
        
        return M

class CROD(LineElement):
    type = 'CROD'
    def __init__(self,card):
        LineElement.__init__(self,card)
        self.id  = card.field(1)
        self.pid = card.field(2,self.id)

        nids = card.fields(3,5)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==2

    def area(self):
        return self.pid.A

    def nsm(self):
        return self.pid.nsm

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodes
        return self.printCard(fields)


class CTUBE(CROD):
    type = 'CTUBE'
    def __init__(self,card):
        CROD.__init__(self,card)
    ###
    def area():
        return self.pid.area()
###

class CONROD(CROD):
    type = 'CONROD'
    def __init__(self,card):
        LineElement.__init__(self,card)
        self.id  = card.field(1)

        nids = card.fields(2,4)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==2
        #print self.nodes
        
        self.mid = card.field(4)
        self.A   = card.field(5)
        self.J   = card.field(6)
        self.c   = card.field(7,0.0)
        self.nsm = card.field(8,0.0)

    def Stiffness(self,bdf): # ,bdf,L,A,E
        #L = norm(r)
        (n1,n2) = self.getNodeIDs()
        node1 = bdf.Node(n1)
        node2 = bdf.Node(n2)

        p1 = bdf.Node(n1).xyz
        p2 = bdf.Node(n2).xyz
        #print "node1 = ",node1
        #print "node2 = ",node2
        L = norm(p1-p2)
        
        if L==0.0:
            msg = 'invalid CROD length=0.0\n%s' %(self.__repr__())
            raise RuntimeError(msg)
        
        A = self.A
        mat = bdf.Material(self.mid)
        E = mat.E
        #print "L = ",L
        K = A*E/L*matrix([[1.,-1.],[-1.,1.]]) # rod
        
        #print "K = \n",K
        return K

    def __repr__(self):
        c   = self.setBlankIfDefault(self.c,  0.0)
        nsm = self.setBlankIfDefault(self.nsm,0.0)
        #print "nodes",self.nids
        fields = [self.type,self.eid]+self.nids+[self.mid,self.A,self.J,self.c,self.nsm]
        return self.printCard(fields)


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
    def __init__(self,card):
        LineElement.__init__(self,card)
        self.pid = card.field(2)
        self.ga = card.field(3)
        self.gb = card.field(4)
        self.initX_G0(card)

        self.offt = card.field(8,'GGG')
        #print 'self.offt = |%s|' %(self.offt)
        assert self.offt[0] in ['G','B','O'],'invalid offt parameter of %s...offt=%s' %(self.type,self.offt)
        assert self.offt[1] in ['G','B','O'],'invalid offt parameter of %s...offt=%s' %(self.type,self.offt)
        assert self.offt[2] in ['G','B','O'],'invalid offt parameter of %s...offt=%s' %(self.type,self.offt)

        self.pa = card.field(9)
        self.pb = card.field(10)
        
        self.w1a = card.field(11,0.0)
        self.w2a = card.field(12,0.0)
        self.w3a = card.field(13,0.0)

        self.w1b = card.field(14,0.0)
        self.w2b = card.field(15,0.0)
        self.w3b = card.field(16,0.0)

    def area(self):
        return self.pid.A

    def length(self):
        L = self.length2(self.ga,self.gb)
        return L

    def nsm(self):
        return self.pid.nsm

    def initX_G0(self,card):
        field5 = card.field(5)
        if isinstance(field5,int):
            self.g0 = field5
            self.x1 = None
            self.x2 = None
            self.x3 = None
        elif isinstance(field5,float):
            self.g0 = None
            self.x1 = card.field(5,0.0)
            self.x2 = card.field(6,0.0)
            self.x3 = card.field(7,0.0)
        else:
            msg = 'field5 on %s is the wrong type...id=%s field5=%s' %(self.type,self.id,field5)
            raise RuntimeError(msg)
        #if self.id==14100238:
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
        self.ga = mesh.Node(self.ga)
        self.gb = mesh.Node(self.gb)
        self.pid  = mesh.Property(self.pid)

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

    def __repr__(self):
        w1a = self.setBlankIfDefault(self.w1a,0.0)
        w2a = self.setBlankIfDefault(self.w2a,0.0)
        w3a = self.setBlankIfDefault(self.w3a,0.0)
        w1b = self.setBlankIfDefault(self.w1b,0.0)
        w2b = self.setBlankIfDefault(self.w2b,0.0)
        w3b = self.setBlankIfDefault(self.w3b,0.0)
        (x1,x2,x3) = self.getX_G0_defaults()
        offt = self.setBlankIfDefault(self.offt,'GGG')
        fields = [self.type,self.eid,self.Pid(),self.Ga(),self.Gb(),x1,x2,x3,offt,
                  self.pa,self.pb,w1a,w2a,w3a,w1b,w2b,w3b]

        return self.printCard(fields)

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
    def __init__(self,card):
        LineElement.__init__(self,card)
        self.pid = card.field(2,self.id)
        self.ga = card.field(3)
        self.gb = card.field(4)

        self.initX_G0(card)
        self.initOfftBit(card)
        self.pa = card.field(9)
        self.pb = card.field(10)
        
        self.w1a = card.field(11,0.0)
        self.w2a = card.field(12,0.0)
        self.w3a = card.field(13,0.0)

        self.w1b = card.field(14,0.0)
        self.w2b = card.field(15,0.0)
        self.w3b = card.field(16,0.0)

        self.sa = card.field(17)
        self.sb = card.field(18)
    
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

    def getOfft_Bit_defaults(self):
        if self.isOfft:
            field8 = self.offt
        else:
            field8 = self.setBlankIfDefault(self.bit,0.0)
        return field8

    def crossReference(self,model):
        self.ga = mesh.Node(self.ga)
        self.gb = mesh.Node(self.gb)
        self.pid  = mesh.Property(self.pid)

    def Stiffness(self,bdf,r,A,E,I):
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

    def __repr__(self):
        w1a = self.setBlankIfDefault(self.w1a,0.0)
        w2a = self.setBlankIfDefault(self.w2a,0.0)
        w3a = self.setBlankIfDefault(self.w3a,0.0)
        w1b = self.setBlankIfDefault(self.w1b,0.0)
        w2b = self.setBlankIfDefault(self.w2b,0.0)
        w3b = self.setBlankIfDefault(self.w3b,0.0)
        (x1,x2,x3) = self.getX_G0_defaults()
        offt = self.getOfft_Bit_defaults()
        fields = [self.type,self.eid,self.Pid(),self.ga,self.gb,x1,x2,x3,offt,
                  self.pa,self.pb,w1a,w2a,w3a,w1b,w2b,w3b,
                  self.sa,self.sb]
        return self.printCard(fields)
