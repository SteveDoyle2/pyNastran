## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import sys
from numpy import zeros,pi,transpose

# pyNastran
from pyNastran.bdf.cards.baseCard import Property,Material

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
        TinvT = transpose(Tinv)
        return (Tinv,TinvT)

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
            ## Property ID
            self.pid  = card.field(1)
            ## Non-Structural Mass
            self.nsm  = card.field(3,0.0)
            self.sb   = card.field(4,0.0)
            self.ft   = card.field(5)
            self.TRef = card.field(6,0.0)
            self.ge   = card.field(7,0.0)

            ## symmetric flag - default = No Symmetry (NO)
            self.lam  = card.field(8)
            #print "lam = ",self.lam

            nPlyFields = card.nFields()-9 # -8 for the first 8 fields (1st line)
            #plyCards = card.fields(9)

            # counting plies
            nMajor    = nPlyFields//4
            nLeftover = nPlyFields%4
            if nLeftover:
                nMajor+=1
            nplies = nMajor
            #print "nplies = ",nplies

            iPly = 1
            plies = []
            midLast = None
            tLast = None

            ## supports single ply per line
            for i in range(9,9+nplies*4,4):
                defaults = [midLast,tLast,0.0,'NO']
                (mid,t,theta,sout) = card.fields(i,i+4,defaults)
                ply = [mid,t,theta,sout]
                
                assert t>0.,'thickness of PCOMP layer is invalid iLayer=%s t=%s ply=[mid,t,theta,sout]=%s' %(iPly,t,ply)
                if ply!=defaults: # if they're not all defaults...
                    plies.append(ply)
                midLast = mid
                tLast = t
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
        #print self
        #print "nPlies = %s" %(self.nPlies())
        #if self.pid==2058:
        #    sys.exit()

    def hasCoreLayer(self):
        """is there a center layer (matters most for a symmetrical ply)"""
        return self.nPlies()%2==1 # True if has a core, False otherwise

    def nPlies(self):
        """
        how many plies are there?
        @note 
            returns half the number of actual plies if the ply is symmetrical (not considering the core)
        """
        return len(self.plies)

    def isSymmetrical(self):
        """is the laminate symmetrical laminate"""
        if self.lam=='SYM':
            return True
        return False

    def isSameCard(self,prop):
        if self.type!=mat.type:  return False
        fields2 = [prop.nsm,prop.sb,prop.ft,prop.TRef,prop.ge,prop.lam]
        fields1 = [self.nsm,self.sb,self.ft,self.TRef,self.ge,self.lam]

        for ply in self.plies:
            fields1 += ply
        for ply in prop.plies:
            fields2 += ply

        for (field1,field2) in zip(fields1,fields2):
            if not self.isSame(field1,field2):
                return False
            ###
        ###
        return True

    def crossReference(self,model):
        """
        links the material ID to the materials
        @param self the object pointer
        @param model a BDF object
        """
        for iPly,ply in enumerate(self.plies):
            mid = self.plies[iPly][0]
            #print mid
            self.plies[iPly][0] = model.Material(mid)  # mid
        ###

    def Nsm(self):
        return self.nsm

    def Mid(self,iPly):
        """
        gets the material ID of the ith ply
        @param self the object pointer
        @param iPly the ply ID (starts from 0)
        """
        Mid = self.Material(iPly)
        if isinstance(Mid,int):
            return Mid
        return Mid.mid

    def Mids(self):
        """
        gets the material IDs of all the plies
        @param self the object pointer
        @retval mids the material IDs
        """
        mids = []
        for iPly in range(len(self.plies)):
            mids.append(self.Mid(iPly))
        return mids

    def Rho(self,iPly):
        """
        gets the density of the ith ply
        @param self the object pointer
        @param iPly the ply ID (starts from 0)
        """
        mid = self.Material(iPly)
        return mid.rho

    def Material(self,iPly):
        """
        gets the material of the ith ply (not the ID unless it's not cross-referenced)
        @param self the object pointer
        @param iPly the ply ID (starts from 0)
        """
        Mid = self.plies[iPly][0]
        return Mid

    def Thickness(self,iPly='all'):
        """
        gets the thickness of the ith ply
        @param self the object pointer
        @param iPly the string 'all' (default) or the mass per area of the ith ply
        """
        if iPly=='all': # get all layers
            t = 0.
            for (iply,ply) in enumerate(self.plies):
                t += self.Thickness(iply)
            ###
            if self.isSymmetrical():
                if self.hasCoreLayer():
                    t -= self.Thickness(0)/2.  # cut out the thickness of half the core layer
                return t*2.
            return t
            ###
        ###
        else:
            t = self.plies[iPly][1]
            return t
        ###

    def Theta(self,iPly):
        """
        gets the ply angle of the ith ply (not the ID)
        @param self the object pointer
        @param iPly the ply ID (starts from 0)
        """
        Theta = self.plies[iPly][2]
        return Theta

    def sout(self,iPly):
        Sout = self.plies[iPly][3]
        return Sout

    def MassPerArea(self,iPly='all'):
        """
        mass = rho*A*t
        but area comes from the element
        mass/A = rho*t for the various layers
        the final mass calculation will be done later
        @param self the object pointer
        @param iPly the string 'all' (default) or the mass per area of the ith ply
        """
        if iPly=='all': # get all layers
            massPerArea = 0.
            for (iply,ply) in enumerate(self.plies):
                massPerArea += self.MassPerArea(iply)
            ###
            if self.isSymmetrical():
                if self.hasCoreLayer():
                    massPerArea -= self.MassPerArea(0)/2.  # cut out the thickness of half the core layer
                    
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

    def D(self):
        D = zeros([3,3])
        isSym = self.isSymmetrical()
        for (iply,ply) in enumerate(self.plies):
            theta = self.Theta(iply)
            t     = self.Thickness(iply)
            mat   = self.Material(iply)
            Di = mat.Dplate()
            transform = self.T(theta)
            D += Di*transform
        return D

    def rawFields(self):
        fields = ['PCOMP',self.pid,self.z0,self.nsm,self.sb,self.ft,self.TRef,self.ge,self.lam,]
        #print "plies = ",self.plies
        for (iPly,ply) in enumerate(self.plies):
            (_mid,t,theta,sout) = ply
            mid = self.Mid(iPly)
            fields += [mid,t,theta,sout]
        return fields

    def reprFields(self):
        #print "t = ",self.Thickness()
        nsm  = self.setBlankIfDefault(self.nsm, 0.0)
        sb   = self.setBlankIfDefault(self.sb,  0.0)
        TRef = self.setBlankIfDefault(self.TRef,0.0)
        ge   = self.setBlankIfDefault(self.ge,  0.0)
        z0 = self.setBlankIfDefault(self.z0,-0.5*self.Thickness())

        fields = ['PCOMP',self.pid,z0,nsm,sb,self.ft,TRef,ge,self.lam,]
        #print "plies = ",self.plies
        for (iPly,ply) in enumerate(self.plies):
            (_mid,t,theta,sout) = ply
            mid = self.Mid(iPly)
            #theta = self.setBlankIfDefault(theta,0.0)
            sout  = self.setBlankIfDefault(sout,'NO')
            fields += [mid,t,theta,sout]
        return fields

class PCOMPG(PCOMP):  ## @todo check for bugs in ply parser
    type = 'PCOMPG'
    def __init__(self,card=None,data=None):
        ShellProperty.__init__(self,card,data) ## @todo doesnt support data
        self.pid  = card.field(1)
        #self.z0 = 
        self.nsm  = card.field(3,0.0)
        self.sb   = card.field(4,0.0)
        self.ft   = card.field(5)
        self.TRef = card.field(6,0.0)
        self.ge   = card.field(7,0.0)
        self.lam  = card.field(8)
        fields = card.fields(9)

        T = 0. # thickness
        midLast = None
        tLast = None
        self.plies = []
        for i in range(len(fields)):
            gPlyID    = card.field(i)
            mid       = card.field(i+1,midLast)
            thickness = float(card.field(i+2,tLast)) # can be blank 2nd time thru
            theta     = card.field(i+3,0.0)
            sout      = card.field(i+4,'NO')
            self.plies.append([mid,thickness,theta,sout,gPlyID])
            #[mid,t,theta,sout] # PCOMP
            assert mid       is not None
            assert thickness is not None
            midLast = mid
            tLast = thickness
            T+=thickness
        ###

    def GlobalPlyID(self,iPly):
        gPlyID = self.plies[iPly][4]
        return gPlyID
        
    def rawFields(self):
        fields = ['PCOMPG',self.pid,self.z0,self.nsm,self.sb,self.ft,TRef,self.ge,self.lam,]
        for (iPly,ply) in enumerate(self.plies):
            (_mid,t,theta,sout,gPlyID) = ply
            mid = self.Mid(iPly)
            fields += [mid,t,theta,sout,gPlyID,None,None,None]
        return fields

    def reprFields(self):
        nsm  = self.setBlankIfDefault(self.nsm, 0.0)
        sb   = self.setBlankIfDefault(self.sb,  0.0)
        TRef = self.setBlankIfDefault(self.TRef,0.0)
        ge   = self.setBlankIfDefault(self.ge,  0.0)
        z0 = self.setBlankIfDefault(self.z0,-0.5*self.Thickness())

        fields = ['PCOMPG',self.pid,z0,nsm,sb,self.ft,TRef,ge,self.lam,]

        for (iPly,ply) in enumerate(self.plies):
            (_mid,t,theta,sout,gPlyID) = ply
            mid = self.Mid(iPly)
            #theta = self.setBlankIfDefault(theta,0.0)
            sout  = self.setBlankIfDefault(sout,'NO')
            fields += [mid,t,theta,sout,gPlyID,None,None,None]
        return fields

class PSHEAR(ShellProperty):
    type = 'PSHEAR'
    def __init__(self,card=None,nPDAMP=0,data=None):
        """
        Defines the properties of a shear panel (CSHEAR entry).
        PSHEAR PID MID T NSM F1 F2
        """
        ShellProperty.__init__(self,card,data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            ## Material ID
            self.mid = card.field(2)
            self.t   = card.field(3)
            self.nsm = card.field(4,0.0)
            self.f1  = card.field(5,0.0)
            self.f2  = card.field(6,0.0)
        else:
            #(pid,mid,t,nsm,f1,f2) = out
            self.pid = data[0]
            self.mid = data[1]
            self.t   = data[2]
            self.nsm = data[3]
            self.f1  = data[4]
            self.f2  = data[5]
        ###

    def isSameCard(self,prop):
        if self.type!=prop.type:  return False
        fields1 = self.rawFields()
        fields2 = prop.rawFields()
        return self.isSameFields(fields1,fields2)

    def rawFields(self):
        fields = ['PSHEAR',self.pid,self.Mid(),self.t,self.nsm,self.f1,self.f2]
        return fields

class PSHELL(ShellProperty):
    """
    PSHELL PID MID1 T MID2 12I/T**3 MID3 TS/T NSM
    Z1 Z2 MID4
    PSHELL   41111   1      1.0000   1               1               0.02081"""
    type = 'PSHELL'
    def __init__(self,card=None,data=None):
        ShellProperty.__init__(self,card,data)
        if card:
            self.pid  = int(card.field(1))
            self.mid1 = card.field(2)
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
            self.mid1      = data[1]
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

    def mid(self):
        if isinstance(self.mid1,Material):
            return self.mid1
        return self.mid2

    def Mid(self):
        if isinstance(self.mid1,Material):
            return self.mid1.mid  # self.Mid1()
        return self.Mid2()
    
    def Mid1(self):
        if isinstance(self.mid1,Material):
            return self.mid1.mid
        return self.mid1

    def Mid2(self):
        if isinstance(self.mid2,Material):
            return self.mid2.mid
        return self.mid2

    def Mid3(self):
        if isinstance(self.mid3,Material):
            return self.mid3.mid
        return self.mid3

    def Mid4(self):
        if isinstance(self.mid4,Material):
            return self.mid4.mid
        return self.mid4

    def Thickness(self):
        return self.t

    def Rho(self):
        return self.mid().rho

    def Nsm(self):
        return self.nsm

    def MassPerArea(self):
        """
        calculates mass per area
        \f[ \large  \frac{mass}{A} = nsm + \rho t \f]
        """
        massPerArea = self.nsm + self.Rho()*self.t
        return massPerArea

    def D(self):
        return self.mid().Dplate()
        
    def crossReference(self,mesh):
        if self.mid1:
            self.mid1 = mesh.Material(self.mid1)
        if self.mid2:
            self.mid2 = mesh.Material(self.mid2)
        if self.mid3:
            self.mid3 = mesh.Material(self.mid3)
        if self.mid4:
            self.mid4 = mesh.Material(self.mid4)

    def writeCodeAster(self):
        """
        http://www.caelinux.org/wiki/index.php/Contrib:KeesWouters/shell/static
        http://www.caelinux.org/wiki/index.php/Contrib:KeesWouters/platedynamics
        # the angle_rep is a direction angle, use either angle(a,b) or vecteur(x,y,z)
        # coque_ncou is the number of gauss nodes along the thickness, for linear analysis one node is sufficient.
        """
        msg = ''
        msg += "    COQUE=_F(GROUP_MA='P%s', # COQUE=PSHELL\n" %(self.pid)
        msg += "              EPAIS=%g, # EPAIS=thickness\n" %(self.t)
        msg += "              ANGL_REP=(0.,90.),     # ???\n"  ## @todo what is this?
        #msg += "              VECTEUR=(1.0,0.0,0.0,) #  local coordinate system\n"
        msg += "              EXCENTREMENT=%g,       # offset-Z1\n" %(self.z1)
        msg += "              COQUE_NCOU=1,          # Number of Integration Layers\n"
        msg += "              CARA=('NSM'), # ???\n"  ## @todo check
        msg += "              VALE=(%g),),\n"  %(self.nsm)
        return msg

    def rawFields(self):
        fields = ['PSHELL',self.pid,self.Mid1(),self.t,self.Mid2(),self.twelveIt3,self.Mid3(),self.tst,self.nsm,
                           self.z1,self.z2,self.Mid4()]
        return fields

    def reprFields(self):
        twelveIt3 = self.setBlankIfDefault(self.twelveIt3,1.0)
        tst       = self.setBlankIfDefault(self.tst,0.833333)
        nsm       = self.setBlankIfDefault(self.nsm,0.0)

        tOver2 = self.t/2.
        z1 = self.setBlankIfDefault(self.z1,-tOver2)
        z2 = self.setBlankIfDefault(self.z2, tOver2)

        fields = ['PSHELL',self.pid,self.Mid1(),self.t,self.Mid2(),twelveIt3,self.Mid3(),tst,nsm,
                           z1,z2,self.Mid4()]
        #print fields
        return fields

