#import sys
#from numpy import array,cross,dot
#from numpy import array
#from numpy.linalg import norm
from numpy import zeros

# my code
from baseCard import Property

class PointProperty(Property):
    type = 'PointProperty'
    def __init__(self,card):
        Property.__init__(self,card)
        pass

class PMASS(PointProperty):
    def __init__(self,card):
        PointProperty.__init__(self,card,nOffset=0)
        
        nOffset *= 2
        self.pid = card.field(1+nOffset)
        self.Mass = card.field(2+nOffset,0.)

    def mass(self):
        return self.mass

    def __repr__(self):
        fields = ['PMASS',self.pid,self.Mass]

class SpringProperty(Property):
    type = 'SpringProperty'
    def __init__(self,card):
        Property.__init__(self,card)
        pass

class PELAS(SpringProperty):
    type = 'PELAS'
    def __init__(self,card,nPELAS=0):
        SpringProperty.__init__(self,card)
        self.pid = card.field(1+5*nPELAS) # 2 PELAS properties can be defined on 1 PELAS card
        self.k   = card.field(2+5*nPELAS) # these are split into 2 separate cards
        self.ge  = card.field(3+5*nPELAS)
        self.s   = card.field(4+5*nPELAS)

    def __repr__(self):
        fields = ['PELAS',self.pid,self.k,self.ge,self.s]
        return self.printCard(fields)

class LineProperty(Property):
    type = 'LineProperty'
    def __init__(self,card):
        Property.__init__(self,card)
        pass
    def D_bending(self):
        pass
    def D_axial(self):
        pass
    def D_thermal(self):
        pass
    def D_shear(self):
        pass
    
    def rho(self):
        return self.mid.rho

    def nsm(self):
        return self.nsm

    def E(self):
        return self.mid.E

    def G(self):
        return self.mid.G

    def nu(self):
        return self.mid.nu

class PROD(LineProperty):
    type = 'PROD'
    def __init__(self,card):
        LineProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)
        self.A   = card.field(3)
        self.J   = card.field(4)
        self.c   = card.field(5,0.0)
        self.nsm = card.field(6,0.0)

    def crossReference(self,mesh):
        self.mid = mesh.Material(self.mid)

    def __repr__(self):
        c   = self.setBlankIfDefault(self.c,0.0)
        nsm = self.setBlankIfDefault(self.nsm,0.0)
        fields = ['PROD',self.pid,self.Mid(),self.A,self.J,c,nsm]
        return self.printCard(fields)

class PTUBE(LineProperty):
    type = 'PTUBE'
    def __init__(self,card):
        LineProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)
        self.OD1 = card.field(3)
        self.t   = card.field(4,self.outerDiameter/2.)
        self.nsm = card.field(5,0.0)
        self.OD2 = card.field(6,self.outerDiameter)

    def __repr__(self):
        t   = self.setBlankIfDefault(self.t,self.OD1/2.)
        nsm = self.setBlankIfDefault(self.nsm,0.0)
        OD2 = self.setBlankIfDefault(self.OD2,self.OD1)
        fields = ['PTUBE',self.pid,self.Mid(),self.OD1,t,nsm,OD2]
        return self.printCard(fields)
    
    def area(self):
        """
        this shouldnt be so hard...
        """
        return (self.area1()+self.area2())/2.
        factor = pi()
        A = 1.
        #Di = Do-2t
        #Ri = Ro-t
        #A/pi = Ro^2-Ri^2 = Ro^2-(Ro-t)^2 = Ro^2-(Ro-t)(Ro-t)
        #A/pi = Ro^2-(Ro^2-2Ro*t + tt)
        #A/pi = 2Ro*t-tt = t*(2Ro-t)
        #R/2*(2R-R/2) = R/2*(3/2*R) = RR3/4
        
        #2A/pi = t*(2Ro1-t)+(2Ro2-t)*t = t*(2Ro1-2t+2Ro2)
        #factor = pi/2
        #A = pi*t*(Ro1+Ro2-t)
        #A = pi*t*(Do-t)
        
        
        #A/pi = DoDo/4 - (Do/2-t)*(Do/2-t)
        #A/pi = DoDo/4 - (DoDo/4 -2*Do/2*t+tt)
        #A/pi = DoDo/4 - (DoDo/4 -Do/t    +tt)
        #A/pi = Do*t - tt
        
        #A/pi    = t*(Do-t)
        #Aavg/pi = t*(Do1-t + Do2-t)/2.
        #Aavg/pi = t*(Do1+Do2 - 2*t)/2.
        factor = pi()/2.
        #A = self.t*(self.OD1+self.OD2 - 2*self.t)
        return A*factor

    def area1(self):
        Dout = self.OD
        Din  = Dout-self.t
        A = pi()/4.*(Dout*Dout-Din*Din)
        return A1

    def area2(self):
        Dout = self.OD2
        Din  = Dout-self.t
        A = pi()/4.*(Dout*Dout-Din*Din)
        return A2

    def massMatrix(self):
        """@todo not done"""
        m = zeros(6,6)
        m[0,0] = 1.
        return m

class PBAR(LineProperty):
    type = 'PBAR'
    def __init__(self,card):
        """
        @todo
            support solution 600 default
            do a check for mid -> MAT1      for structural
            do a check for mid -> MAT4/MAT5 for thermal
        """
        LineProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)
        self.A   = card.field(3,0.0)
        self.I1  = card.field(4,0.0)
        self.I2  = card.field(5,0.0)
        self.J   = card.field(6,0.0) # default=1/2(I1+I2) for SOL=600, otherwise 0.0
        self.nsm = card.field(7,0.0)

        self.C1  = card.field(9, 0.0)
        self.C2  = card.field(10,0.0)
        self.D1  = card.field(11,0.0)
        self.D2  = card.field(12,0.0)
        self.E1  = card.field(13,0.0)
        self.E2  = card.field(14,0.0)
        self.F1  = card.field(15,0.0)
        self.F2  = card.field(16,0.0)

        self.K1  = card.field(17)
        self.K2  = card.field(18)
        self.I12 = card.field(19,0.0)
        if self.A==0.0:
            assert self.K1==None
            assert self.K2==None

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        A   = self.setBlankIfDefault(self.A,0.0)
        I1  = self.setBlankIfDefault(self.I1,0.0)
        I2  = self.setBlankIfDefault(self.I2,0.0)
        I12 = self.setBlankIfDefault(self.I12,0.0)
        J   = self.setBlankIfDefault(self.J,0.0)
        nsm = self.setBlankIfDefault(self.nsm,0.0)
        
        C1  = self.setBlankIfDefault(self.C1,0.0)
        C2  = self.setBlankIfDefault(self.C2,0.0)

        D1  = self.setBlankIfDefault(self.D1,0.0)
        D2  = self.setBlankIfDefault(self.D2,0.0)

        E1  = self.setBlankIfDefault(self.E1,0.0)
        E2  = self.setBlankIfDefault(self.E2,0.0)

        F1  = self.setBlankIfDefault(self.F1,0.0)
        F2  = self.setBlankIfDefault(self.F2,0.0) # must have 1 on line, if line3 is not empty
        
        line3 = [self.K1,self.K2,I12]
        #print "line3 = ",line3
        
        line1 = ['PBAR',self.pid,self.Mid(),self.A,I1,I2,J,nsm,None]

        if line3==[None,None,None]:
            line2 = [C1,C2,D1,D2,E1,E2,F1,F2]
        else:
            line2 = [C1,C2,D1,D2,E1,E2,F1,self.F2]
        fields = line1+line2+line3

        return self.printCard(fields)

class PBARL(LineProperty): # not done, what if all of dim is blank and no nsm...
    type = 'PBARL'
    validTypes = ["ROD", "TUBE", "I", "CHAN", "T", "BOX", "BAR", "CROSS", "H",
    "T1", "I1", "CHAN1", "Z", "CHAN2", "T2", "BOX1", "HEXA", "HAT",
    "HAT1", "DBOX"] # for GROUP="MSCBML0"

    def __init__(self,card):
        LineProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)
        self.group = card.field(3,'MSCBMLO')
        self.type = card.field(4)
        assert type in self.validTypes

        self.dim = fields(9) # confusing entry...
        nDim = len(self.dim)-1
        if nDim>0:
            self.nsm = self.dim.pop()
        else:
            self.nsm = 0.0
        ###

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def nsm():
        return self.nsm

    def __repr__(self):
        fields = ['PBARL',self.pid,self.Mid(),group,type,None,None,None,None,
        ]+self.dim+[self.nsm]
        return self.printCard(fields)

class PBEAM(LineProperty):
    type = 'PBEAM'
    def __init__(self,card):
        """
        @todo fix 0th entry of self.so, self.xxb
        """
        LineProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)

        self.so  = [None] ## @todo what are these values???
        self.xxb = [None]
        self.A   = [card.field(3) ]
        self.i1  = [card.field(4) ]
        self.i2  = [card.field(5) ]
        self.i12 = [card.field(6) ]
        self.J   = [card.field(7) ]
        self.NSM = [card.field(8,0.0) ]
        self.c1  = [card.field(9) ]
        self.c2  = [card.field(10)]
        self.d1  = [card.field(11)]
        self.d2  = [card.field(12)]
        self.e1  = [card.field(13)]
        self.e2  = [card.field(14)]
        self.f1  = [card.field(15)]
        self.f2  = [card.field(16)]
        
        #fields = card.fields(0)
        #print "fieldsPBEAM = ",fields
        #fieldsMid = fields[16:]
        #print "fieldsMid = ",fieldsMid

        #fields = card.fields(9)
        nFields = card.nFields()-16 # 17+16 (leading + trailing fields)
        # counting continuation cards
        nMajor    = nFields/16
        nLeftover = nFields%16
        #print "nMajor=%s nLeftover=%s" %(nMajor,nLeftover)
        if nLeftover==0:
            nMajor-=1

        #print "nMajor = ",nMajor
        for nRepeated in range(1,nMajor+1):
            nStart = nRepeated*16+1  # field 17 is the first possible so
            propFields = card.fields(nStart,nStart+16)
            #print "propFields = ",propFields
            self.so.append( propFields[0])
            self.xxb.append(propFields[1])
            self.A.append(  propFields[2])
            self.i1.append( propFields[3])
            self.i2.append( propFields[4])
            self.i12.append(propFields[5])
            self.J.append(  propFields[6])
            self.NSM.append(self.setDefaultIfBlank(propFields[7],0.0))
            self.c1.append( propFields[8])
            self.c2.append( propFields[9])
            self.d1.append( propFields[10])
            self.d2.append( propFields[11])
            self.e1.append( propFields[12])
            self.e2.append( propFields[13])
            self.f1.append( propFields[14])
            self.f2.append( propFields[15])
        #print "nRepeated = ",nRepeated

        # footer fields
        x = (nMajor)*16+17
        self.k1   = card.field(x)
        self.k2   = card.field(x+1)
        self.s1   = card.field(x+2)
        self.s2   = card.field(x+3)
        self.nsia = card.field(x+4)
        self.nsib = card.field(x+5)
        self.cwa  = card.field(x+6)
        self.cwb  = card.field(x+7)

        self.m1a = card.field(x+8)
        self.m2a = card.field(x+9)
        self.m1b = card.field(x+10)
        self.m2b = card.field(x+11)
        self.n1a = card.field(x+12)
        self.n2a = card.field(x+13)
        self.n1b = card.field(x+14)
        self.n2b = card.field(x+15)

    def area(self):
        """@warning area field not supported fully on PBEAM card"""
        #raise Exception(self.A[0])
        return self.A[0]

    def nsm(self):
        """@warning nsm field not supported fully on PBEAM card"""
        #raise Exception(self.nsm[0])
        return self.NSM[0]

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        fields = ['PBEAM',self.pid,self.Mid()]
        #print "fieldsA = ",fields
        
        #print len(self.so)
        i = 0
        for (so,xxb,A,i1,i2,i12,J,NSM,c1,c2,d1,d2,e1,e2,f1,f2) in zip(
            self.so,self.xxb,self.A,self.i1,self.i2,self.i12,self.J,self.NSM,
            self.c1,self.c2,self.d1,self.d2,self.e1,self.e2,self.f1,self.f2):
            NSM = self.setBlankIfDefault(NSM,0.0)
            if i==0:
                fields +=        [A,i1,i2,i12,J,NSM,c1,c2,d1,d2,e1,e2,f1,f2] # the 1st 2 fields aren't written
            else:
                fields += [so,xxb,A,i1,i2,i12,J,NSM,c1,c2,d1,d2,e1,e2,f1,f2]
            i+=1
        fields += [self.k1,self.k2,self.s1,self.s2,self.nsia,self.nsib,self.cwa,self.cwb,
                   self.m1a,self.m2a,self.m1b,self.m2b,self.n1a,self.n2a,self.n1b,self.n2b]
        #print fields
        return self.printCard(fields)
        
#class PBEAML(LineProperty): #not done
class PBEAM3(LineProperty): # not done, cleanup
    type = 'PBEAM3'
    def __init__(self,card):
        LineProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)

        self.A   = card.field(3)
        self.Iz  = card.field(4)
        self.Iy  = card.field(5)
        self.Iyz = card.field(6,0.0)
        self.J   = card.field(7,self.Iy+self.Iz)
        self.nsm = card.field(8,0.0)

        self.cy = card.field(9)
        self.cz = card.field(10)
        self.dy = card.field(11)
        self.dz = card.field(12)
        self.ey = card.field(13)
        self.dz = card.field(14)
        self.fy = card.field(15)
        self.fz = card.field(16)
        # more...

    def nsm(self):
        """@warning nsm field not supported fully on PBEAM3 card"""
        #raise Exception(self.nsm[0])
        return self.nsm

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        raise Exception('not done...')
        fields = ['PBEAM3',self.pid,self.Mid(),] # other
        return self.printCard(fields)

class ShellProperty(Property):
    type = 'ShellProperty'
    def __init__(self,card):
        Property.__init__(self,card)
        pass
    

    def S(self):
        """"
        Calculates the compliance matrix for a lamina
        [Q] = [S]^-1
        @todo finish...if necessary...
        """
        pass
        #return Q.inv()

    def ABDH(self):
        """
        tranforms load to strain/bending curvature taken at \f$ z=0 \f$
        \f[ \large  [T] \left[ 
          \begin{array}{c}
              Nx  \\
              Ny  \\
              Nz  \\
              Mx  \\
              My  \\
              Mz  \\
          \end{array} \right] = 

          [ABBD] \left[ 
          \begin{array}{c}
              \epsilon_{xx}  \\
              \epsilon_{yy}  \\
                \gamma_{xy}  \\
                \kappa_{xx}  \\
                \kappa_{yy}  \\
                \kappa_{xy}  \\
          \end{array} \right]
        \f] 


        [Nx] = [            ] [ e_xx0    ]
        [Ny] = [  [A]   [B] ] [ e_yy0    ]
        [Nz] = [            ] [ gamma_xy0]
        [Mx] = [            ] [ k_xx0    ]
        [My] = [  [B]   [D] ] [ k_yy0    ]
        [Mz] = [            ] [ k_xy0    ]
        
        \f[ \large  A_{ij} = \Sigma_{k=1}^N (\overline{Q_{ij}})_k (z_k  -z_{k-1}  ) = \Sigma_{k=1}^N (Q_{ij})_k t_k                                     \f]
        \f[ \large  B_{ij} = \Sigma_{k=1}^N (\overline{Q_{ij}})_k (z_k^2-z_{k-1}^2) = \Sigma_{k=1}^N (Q_{ij})_k (\overline{z} t_k)                      \f]
        \f[ \large  D_{ij} = \Sigma_{k=1}^N (\overline{Q_{ij}})_k (z_k^3-z_{k-1}^3) = \Sigma_{k=1}^N (Q_{ij})_k         (\overline{z}^2 t_k + t_k^3/12) \f]
        \f[ \large  H_{ij} =                                                      \Sigma_{k=1}^N (Q_{ij})_k (t_k -4/t^2 (\overline{z}^2 t_k + t_k^3/12) \f]
        
        p. 138 of "Introduction to Composite Material Design"
        """
        A = zeros(9,'d')
        B = copy.deepcopy(A)
        D = copy.deepcopy(A)
        H = copy.deepcopy(A)

        for i,layer in enumerate(self.layers):
            t = layer.t
            z0 = layer.z
            #z1 = z0+t
            zbar = (2*z0+t)/2. # (z1+z0)/2. 
            Qraw   = layer.Q() # needs E11, E22, G12, G13, nu12, nu21
            qlayer = Qall(Qout,layer.thetad)
            A += qlayer*t
            B += qlayer*t*z
            
            Dfactor = t*zbar*zbar+t**3/12.
            D += qlayer*Dfactor
            H += qlayer*(t-4./t**2 * Dfactor)
        B = 0.5*B
        D = D/3.
        H = H*5./4.
            
            
    def Qall(self,thetad):
        """
        Caculates the laminate tranformation  stiffness \f$ [Q]_{all} \f$
        \f[ \large  [Q]_{all} = [T]^{-1} [Q] [R][T][R]^{-1}  \f]
        \f[ \large  [Q]_{all} = [T]^{-1} [Q] [T]^{-T}        \f]
        
        p. 123 of "Introduction to Composite Material Design"
        """
        theta = radians(thetad)
        ct  = cos(theta)
        c2t = ct*ct
        c3t = ct*c2t
        c4t = c2t*c2t

        st  = sin(theta)
        s2t = st*st
        s3t = st*s2t
        s4t = s2t*s2t
        
        s2c2t = s2t*c2t
        #s4tpc4t = s4t+c4t
        
        Q11a = Q11*c4t + 2*(Q12+2*Q66)*s2c2t + Q22*s4t
        Q12a = (Q11+Q22-4*Q66)*s2c2t + Q12(s4t+c4t)
        Q22a = Q11*s4t + 2*(Q12+2*Q66)*s2c2t+Q22*c4t
        Q16a = (Q11-Q12-2*Q66)*st*c3t + (Q12-Q22+2*Q66)*s3t*ct
        Q26a = (Q11-Q12-2*Q66)*s3t*ct + (Q12-Q22+2*Q66)*st*c3t
        Q66a = (Q11+Q22-2*Q12-2*Q66)*s2c2t + Q66(s4t+c4t)
        Q44a = Q44*c2t + Q55*s2t
        Q55a = Q44*s2t + Q55*c2t
        Q45a = (Q55-Q44)*st*ct
        return array([Q11a,Q12a,Q22a,Q16a,Q26a,Q66a,Q44a,Q55a,Q45a])

    def Q(self):
        """"
        Calculates the stiffness matrix \f$ [Q] \f$ for a lamina
        @todo is this done?
        p. 114 "Introduction to Composite Material Design"
        """
        delta = 1-nu12*nu21
        Q11 = E1/delta
        Q12 = nu12*E2/delta
        Q22 = E2/delta
        Q66 = G12
        Q44 = G23
        Q55 = G13
        Qout = (Q11,Q22,Q12,Q44,Q55,Q66)
        return Qout
        
    def T(self,theta):
        """
        calculates the Transformation matricies \$ [T] \$ and  \f$ [T]^{-1} \f$
        @param self           the object pointer
        @param theta          in radians...
        @retval Tinv          the inverse transformation matrix
        @retval TinvTranspose the transposed inverse transformation matrix
        @todo document better

        tranformation matrix  \f$ [T] \f$
        \f[ \large  [T] \left[ 
          \begin{array}{ccc}
              m^2 & n^2 &  2mn    \\
              n^2 & m^2 & -2mn    \\
              -mn & mn  & m^2-n^2 
          \end{array} \right]
        \f] 

                 [ m^2  n^2        2mn]
        [T]    = [ n^2  m^2       -2mn]   # transformation matrix
                 [ -mn   mn    m^2-n^2]


        inverse transformation matrix \f$ [T]^{-1} \f$
        \f[ \large  [T]^{-1} \left[ 
          \begin{array}{ccc}
              m^2 & n^2 & -2mn    \\
              n^2 & m^2 &  2mn    \\
              mn  & -mn & m^2-n^2 
          \end{array} \right]
        \f] 

                 [ m^2  n^2       -2mn]
        [T]^-1 = [ n^2  m^2        2mn]   # inverse transform
                 [ mn   -mn    m^2-n^2]
        
        \f[ \large  \left[ \sigma_{xx}, \sigma_{yy}, \sigma_{xy} \right]^T = [T]^{-1} [Q] [R][T] \left[ \epsilon_{xx}, \epsilon_{yy}, \frac{1}{2} \gamma_{xy} \right]^T \f]
        
        p.119 "Introduction to Composite Material Design"
        """
        n = cos(theta)
        m = sin(theta)
        Tinv  = zeros((6,6),'d')
        mm = m*m
        nn = n*n
        nm = n*m
        Tinv[0][0] = Tinv[1][1] = mm
        Tinv[1][0] = Tinv[0][1] = nn
        Tinv[0][2] = -2*mn
        Tinv[1][2] =  2*mn

        Tinv[2][0] =  mn
        Tinv[2][1] = -mn
        Tinv[2][2] =  mm-nn
        TinvT = numpy.transpose(Tinv)
        return (Tinv,TinvT)

#class PCOMPG(ShellProperty): # not done...
class PCOMP(ShellProperty):
    """
    PCOMP     701512   0.0+0 1.549-2                   0.0+0   0.0+0     SYM
              300704   3.7-2   0.0+0     YES  300704   3.7-2     45.     YES
              300704   3.7-2    -45.     YES  300704   3.7-2     90.     YES
              300705      .5   0.0+0     YES
    """
    type = 'PCOMP'
    def __init__(self,card): # not done, cleanup
        ShellProperty.__init__(self,card)
        
        #fields = card.fields(1)

        self.pid = card.field(1)
        self.z0   = card.field(2)
        self.nsm  = card.field(3,0.0)
        self.sb   = card.field(4,0.0)
        self.ft   = card.field(5)
        self.TRef = card.field(6,0.0)
        self.ge   = card.field(7,0.0)
        
        ## symmetric flag - default = No Symmetry (NO)
        self.lam  = card.field(8)
        #print "lam = ",self.lam
        
        nPlyFields = card.nFields()-8 # -8 for the first 8 fields (1st line)
        #plyCards = card.fields(9)
        
        # counting plies
        nMajor    = nPlyFields/4
        nLeftover = nPlyFields%4
        if nLeftover:
            nMajor+=1
        nplies = nMajor

        #iPly = 0
        plies = []
        for i in range(9,nplies*4,4):  # doesnt support single ply per line
            defaults = [None,None,0.0,'NO']
            (mid,t,theta,sout) = card.fields(i,i+4,defaults)
            ply = [mid,t,theta,sout]
            if ply!=defaults: # if they're not all defaults...
                plies.append(ply)
            #iPly +=1
        #print "nplies = ",nplies
        
        ## list of plies
        self.plies = plies
        
        #self.plies = []
        #if self.lam=='SYM':
        #    if nplies%2==1:  # 0th layer is the core layer
        #        plies[0][1] = plies[0][1]/2. # cut the thickness in half to make the ply have an even number of plies, is there a better way???
        #
        #    pliesLower = plies.reverse()
        #    self.plies = pliesLower+plies
        #    #print str(self)
        ###
        #sys.exit()

    def hasCoreLayer(self):
        """is there a center layer (matters most for a symmetrical ply)"""
        return self.nPlies%2==1 # True if has a core, False otherwise

    def nPlies(self):
        """
        how many plies are there?
        @note 
            returns half the number of actual plies if the ply is symmetrical (not considering the core)
        """
        return len(self.plies)

    def isSymmetrical(self):
        """is the ply a symmetrical ply"""
        if self.lam=='SYM':
            return True
        return False

    def crossReference(self,model):
        for iPly,ply in enumerate(self.plies):
            mid = self.plies[iPly][0]
            #print mid
            self.plies[iPly][0] = model.Material(mid)  # mid
        ###

    def material(self,iPly):
        Mid = self.plies[iPly][0]
        return Mid

    #def nsm(self,iPly):
    #    material = self.material(iPly)
    #    return material.nsm

    def nsm(self):
        return self.nsm

    def mid(self,iPly):
        Mid = self.material(iPly)
        if isinstance(Mid,int):
            return Mid
        return Mid.mid

    def theta(self,iPly):
        Theta = self.plies[iPly][2]
        return Theta

    def sout(self,iPly):
        Sout = self.plies[iPly][3]
        return Sout

    def massPerArea(self,iPly='all'):
        """
        mass = rho*A*t
        but area comes from the element
        mass/A = rho*t for the various layers
        the final mass calculation will be done later
        """
        if iPly=='all': # get all layers
            massPerArea = 0.
            for (iply,ply) in enumerate(self.plies):
                massPerArea += self.massPerArea(iply)
            ###
            if self.isSymmetrical():
                if self.isCoreLayer():
                    massPerArea -= self.massPerArea(0)/2.  # cut out the thickness of half the core layer
                    
                ###
                return massPerArea*2.
            ###
            return massPerArea
        ###
        else:
            rho = self.rho(iPly)
            t = self.plies[iPly][1]
            return rho*t+self.nsm
        ###

    def rho(self,iPly):
        mid = self.material(iPly)
        return mid.rho

    def thickness(self,iPly='all'):
        if iPly=='all': # get all layers
            t = 0.
            for (iply,ply) in enumerate(self.plies):
                t += self.thickness(iply)
            ###
            if self.isSymmetrical():
                if self.isCoreLayer():
                    t -= self.thickness(0)/2.  # cut out the thickness of half the core layer
                return t*2.
            return t
            ###
        ###
        else:
            t = self.plies[iPly][1]
            return t
        ###

    def __repr__(self):
        fields = ['PCOMP',self.pid,self.z0,self.nsm,self.sb,self.ft,self.TRef,self.ge,self.lam,]
        #print "plies = ",self.plies
        for (iPly,ply) in enumerate(self.plies):
            (_mid,t,theta,sout) = ply
            mid = self.mid(iPly)
            theta = self.setBlankIfDefault(theta,0.0)
            sout  = self.setBlankIfDefault(sout,'NO')
            fields += [mid,t,theta,sout]
        return self.printCard(fields)

class SolidProperty(Property):
    type = 'SolidProperty'
    def __init__(self,card):
        Property.__init__(self,card)
        pass

class PLSOLID(SolidProperty):
    """
    Defines a fully nonlinear (i.e., large strain and large rotation) hyperelastic solid
    element.
    PLSOLID PID MID STR
    PLSOLID 20 21
    """
    type = 'PLSOLID'
    def __init__(self,card):
        SolidProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)
        self.ge  = card.field(3)
        self.str = card.field(4,'GRID')
        assert self.str in ['GRID','GAUS'],'card=%s doesnt have a valid stress/strain output value set\n'

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        stressStrain = self.setDefaultIfNone(self.str,'GRID')
        fields = ['PLSOLID',self.pid,self.Mid(),stressStrain]
        return self.printCard(fields)

class PSHELL(ShellProperty):
    """
    PSHELL PID MID1 T MID2 12I/T**3 MID3 TS/T NSM
    Z1 Z2 MID4
    PSHELL   41111   1      1.0000   1               1               0.02081"""
    type = 'PSHELL'
    def __init__(self,card):
        ShellProperty.__init__(self,card)
        self.pid  = card.field(1)
        self.mid  = card.field(2)
        self.t    = card.field(3)
        
        ## Material identification number for bending
        self.mid2 = card.field(4)
        ## \f$ \frac{12I}{t^3} \f
        self.twelveIt3 = card.field(5,1.0)  # poor name
        self.mid3  = card.field(6)
        self.tst   = card.field(7,0.833333)
        self.nsm   = card.field(8,0.0)
        self.z1    = card.field(9)
        self.z2    = card.field(10)
        self.mid4  = card.field(11)

    def rho(self):
        return self.mid.rho

    def massPerArea(self):
        """
        calculates mass per area
        \f[ \large  \frac{mass}{A} = nsm + \rho t \f]
        """
        mPerA = self.nsm + self.mid.rho*self.t
        return mPerA

    def crossReference(self,mesh):
        self.mid  = mesh.Material(self.mid)
        #self.mid2 = mesh.Material(self.mid2)
        #self.mid3 = mesh.Material(self.mid3)
        #self.mid4 = mesh.Material(self.mid4)

    def __repr__(self):
        fields = ['PSHELL',self.pid,self.Mid(),self.t,self.mid2,self.twelveIt3,self.mid3,self.tst,self.nsm,
                  self.z1,self.z2,self.mid4]
        return self.printCard(fields)

class PSOLID(SolidProperty):
    """
    PSOLID PID MID CORDM IN STRESS ISOP FCTN
    PSOLID   1       1       0
    PSOLID 2 100 6 TWO GRID REDUCED
    """
    type = 'PSOLID'
    def __init__(self,card):
        SolidProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)
        self.cordm  = card.field(3,0)
        self.integ  = card.field(4)
        self.stress = card.field(5)
        self.isop   = card.field(6)
        self.fctn   = card.field(7,'SMECH')

    def __repr__(self):
        cordm = self.setBlankIfDefault(self.cordm,0)
        fctn  = self.setBlankIfDefault(self.fctn,'SMECH')
        fields = ['PSOLID',self.pid,self.Mid(),cordm,self.integ,self.stress,self.isop,fctn]
        return self.printCard(fields)

class PCONEAX(Property): #not done
    type = 'PCONEAX'
    def __init__(self,card):
        Property.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)
        self.group = card.field(3,'MSCBMLO')
        self.type = card.field(4)
        self.dim = [] # confusing entry...

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        fields = ['PCONEAX',self.pid,self.Mid()]
        return self.printCard(fields)
    
