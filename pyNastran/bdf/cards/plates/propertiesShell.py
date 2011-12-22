import sys
from numpy import zeros,pi

# my code
from .baseCard import Property

class ShellProperty(Property):
    type = 'ShellProperty'
    def __init__(self,card=None,data=None):
        Property.__init__(self,card,data)
        pass

    def S(self):
        """"
        Calculates the compliance matrix for a lamina
        \f[ \large [Q] = [S]^{-1}  \f]
        @todo finish...if necessary...
        """
        return Q.inv()

    def ABDH(self):
        """
        tranforms load to strain/bending curvature taken at \f$ z=0 \f$
        \f[ \large  \left[ 
          \begin{array}{c}
              Nx  \\
              Ny  \\
              Nz  \\
              Mx  \\
              My  \\
              Mz  \\
          \end{array} \right] = 

          \left[ 
          \begin{array}{cc}
              A & B  \\
              B & D  \\
          \end{array} \right]
          
          \left[ 
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
        
        \f[ \large  A_{ij} = \Sigma_{k=1}^N (\overline{Q_{ij}})_k \left( z_k  -z_{k-1}   \right) = \Sigma_{k=1}^N (Q_{ij})_k t_k                                                                                    \f]
        \f[ \large  B_{ij} = \Sigma_{k=1}^N (\overline{Q_{ij}})_k \left( z_k^2-z_{k-1}^2 \right) = \Sigma_{k=1}^N (Q_{ij})_k                           \left( \overline{z} t_k                      \right)         \f]
        \f[ \large  D_{ij} = \Sigma_{k=1}^N (\overline{Q_{ij}})_k \left( z_k^3-z_{k-1}^3 \right) = \Sigma_{k=1}^N (Q_{ij})_k                           \left( \overline{z}^2 t_k + \frac{t_k^3}{12} \right)         \f]
        \f[ \large  H_{ij} =                                                                       \Sigma_{k=1}^N (Q_{ij})_k \left( t_k -\frac{4}{t^2} \left( \overline{z}^2 t_k + \frac{t_k^3}{12} \right) \right) \f]
        
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
        \f[ \large  [T] = \left[ 
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
        \f[ \large  [T]^{-1} = \left[ 
          \begin{array}{ccc}
              m^2 & n^2 & -2mn    \\
              n^2 & m^2 &  2mn    \\
              mn  & -mn & m^2-n^2 
          \end{array} \right]
        \f] 

                 [ m^2  n^2       -2mn]
        [T]^-1 = [ n^2  m^2        2mn]   # inverse transform
                 [ mn   -mn    m^2-n^2]

        \f[ \large  \left[ 
            \begin{array}{c}
                 \sigma_{xx}   \\
                 \sigma_{yy}   \\
                 \sigma_{xy} 
            \end{array} \right]

             = [T]^{-1} [Q] [R][T]
          
            \left[ \begin{array}{c}
                 \epsilon_{xx} \\
                 \epsilon_{yy} \\
                 \frac{1}{2} \gamma_{xy}
            \end{array} \right]
        \f] 
        
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
    def __init__(self,card=None,data=None): # not done, cleanup
        ShellProperty.__init__(self,card,data)
        
        if card:
            #fields = card.fields(1)

            self.pid = card.field(1)
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

            iPly = 1
            plies = []
            for i in range(9,nplies*4,4):  # doesnt support single ply per line
                defaults = [None,None,0.0,'NO']
                (mid,t,theta,sout) = card.fields(i,i+4,defaults)
                ply = [mid,t,theta,sout]
                assert t>0.,'thickness of PCOMP layer is invalid iLayer=%s t=%s' %(iPly,t)
                if ply!=defaults: # if they're not all defaults...
                    plies.append(ply)
                iPly +=1
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
            self.z0 = card.field(2,-0.5*self.Thickness())
        else:
            #print "len(data) = ",len(data)
            self.pid  = data[0]
            self.z0   = data[1]
            self.nsm  = data[2]
            self.sb   = data[3]
            self.ft   = data[4]
            self.TRef = data[5]
            self.ge   = data[6]
            self.lam  = data[7]
            Mid       = data[8]
            T         = data[9]
            Theta     = data[10]
            Sout      = data[11]

            self.plies = []
            #ply = [mid,t,theta,sout]
            for (mid,t,theta,sout) in zip(Mid,T,Theta,Sout):
                if sout==1:  ## @todo not sure  0=NO,1=YES (most likely)
                    sout='YES'
                elif sout==0:
                    sout= 'NO'
                else:
                    raise Exception('unsupported sout...needs debugging...\nPCOMP = %s' %(data))
                self.plies.append([mid,t,theta,sout])
            ###
        ###

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

    #def Nsm(self,iPly):
    #    material = self.material(iPly)
    #    return material.nsm

    def Nsm(self):
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
            rho = self.Rho(iPly)
            t = self.plies[iPly][1]
            return rho*t+self.nsm
        ###

    def Rho(self,iPly):
        mid = self.material(iPly)
        return mid.rho

    def Thickness(self,iPly='all'):
        if iPly=='all': # get all layers
            t = 0.
            for (iply,ply) in enumerate(self.plies):
                t += self.Thickness(iply)
            ###
            if self.isSymmetrical():
                if self.isCoreLayer():
                    t -= self.Thickness(0)/2.  # cut out the thickness of half the core layer
                return t*2.
            return t
            ###
        ###
        else:
            t = self.plies[iPly][1]
            return t
        ###

    def __repr__(self):
        #print "t = ",self.Thickness()

        nsm  = self.setBlankIfDefault(self.nsm,0.0)
        sb   = self.setBlankIfDefault(self.sb,0.0)
        TRef = self.setBlankIfDefault(self.TRef,0.0)
        ge   = self.setBlankIfDefault(self.ge,0.0)
        z0 = self.setBlankIfDefault(self.z0,-0.5*self.Thickness())

        fields = ['PCOMP',self.pid,z0,nsm,sb,self.ft,TRef,ge,self.lam,]
        #print "plies = ",self.plies
        for (iPly,ply) in enumerate(self.plies):
            (_mid,t,theta,sout) = ply
            mid = self.mid(iPly)
            theta = self.setBlankIfDefault(theta,0.0)
            sout  = self.setBlankIfDefault(sout,'NO')
            fields += [mid,t,theta,sout]
        return self.printCard(fields)

class PSHELL(ShellProperty):
    """
    PSHELL PID MID1 T MID2 12I/T**3 MID3 TS/T NSM
    Z1 Z2 MID4
    PSHELL   41111   1      1.0000   1               1               0.02081"""
    type = 'PSHELL'
    def __init__(self,card=None,data=None):
        ShellProperty.__init__(self,card,data)
        if card:
            self.pid  = card.field(1)
            self.mid  = card.field(2)
            self.t    = card.field(3)

            ## Material identification number for bending
            self.mid2 = card.field(4)
            ## \f$ \frac{12I}{t^3} \f$
            self.twelveIt3 = card.field(5,1.0)  # poor name
            self.mid3  = card.field(6)
            self.tst   = card.field(7,0.833333)
            self.nsm   = card.field(8,0.0)

            tOver2 = self.t/2.
            self.z1    = card.field(9,-tOver2)
            self.z2    = card.field(10,tOver2)
            self.mid4  = card.field(11)
            #if self.mid2 is None:
            #    assert self.mid3 is None
            #else: # mid2 is defined
            #    #print "self.mid2 = ",self.mid2
            #    assert self.mid2 >= -1
            #    #assert self.mid3 >   0

            #if self.mid is not None and self.mid2 is not None:
            #    assert self.mid4==None
            ###
        else:
            self.pid       = data[0]
            self.mid       = data[1]
            self.t         = data[2]
            self.mid2      = data[3]
            self.twelveIt3 = data[4]
            self.mid3      = data[5]
            self.tst       = data[6]
            self.nsm       = data[7]
            self.z1        = data[8]
            self.z2        = data[9]
            self.mid4      = data[10]
            #maxMid = max(self.mid,self.mid2,self.mid3,self.mid4)
        ###

        assert self.t>0.0,'the thickness must be defined on the PSHELL card (Ti field not supported)'

    def Thickness(self):
        return self.t

    def Rho(self):
        return self.mid.rho

    def Nsm(self):
        return self.nsm

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
        twelveIt3 = self.setBlankIfDefault(self.twelveIt3,1.0)
        tst       = self.setBlankIfDefault(self.tst,0.833333)
        nsm       = self.setBlankIfDefault(self.nsm,0.0)

        tOver2 = self.t/2.
        z1 = self.setBlankIfDefault(self.z1,-tOver2)
        z2 = self.setBlankIfDefault(self.z2, tOver2)

        fields = ['PSHELL',self.pid,self.Mid(),self.t,self.mid2,twelveIt3,self.mid3,tst,nsm,
                           z1,z2,self.mid4]
        #print fields
        return self.printCard(fields)

