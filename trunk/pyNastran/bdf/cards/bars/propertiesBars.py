import sys
from numpy import zeros,pi

from pyNastran.general.generalMath import integrateLine,integratePositiveLine

# pyNastran code
from ..baseCard import Property


def IyyBeam(b,h):
    return 1/12.*b*h**3

def IBeam(b,h):
    f = 1/12.*b*h
    
    Iyy = f*h*h # 1/12.*b*h**3
    Izz = f*b*b # 1/12.*h*b**3
    Iyz = 0.
    return (Iyy,Izz,Iyz)

def IBeamOffset(b,h,y,z):
    A = b*h
    f = 1./12.*A
    
    Iyy = f*h*h # 1/12.*b*h**3
    Izz = f*b*b # 1/12.*h*b**3
    Iyz = 0.
    
    Iyy += A*y*y
    Izz += A*z*z
    Iyz += A*y*z
    return (Iyy,Izz,Iyz)


def getInertiaRectangular(sections):
    """
    calculates the moment of inertia for a section about the CG
    @param sections [[b,h,y,z]_1,...] y,z is the centroid (x in the direction of the beam, y right, z up)
    @retval Area,Iyy,Izz,Iyz
    @ see http://www.webs1.uidaho.edu/mindworks/Machine_Design/Posters/PDF/Moment%20of%20Inertia.pdf
    """
    As = []
    Ax=0.; Ay=0.
    for section in sections:
        (b,h,x,y) = section
        A = b*h
        As.append(A)
        Ax += A*x
        Ay += A*y

    xCG = Ax/A
    yCG = Ay/A
    Axx=0; Ayy=0; Axy=0.
    for i,section in enumerate(sections):
        (b,h,x,y) = section
        #A = b*h
        #As.append(A)
        Axx += As[i]*(x-xCG)**2
        Ayy += As[i]*(y-yCG)**2
        Axy += As[i]*(x-xCG)*(y-yCG)
    Ixx = Axx/A
    Iyy = Ayy/A
    Ixy = Axy/A
    return (A,Ixx,Iyy,Ixy)

       


class LineProperty(Property):
    type = 'LineProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass
    def D_bending(self):
        pass
    def D_axial(self):
        pass
    def D_thermal(self):
        pass
    def D_shear(self):
        pass
    
    def Rho(self):
        return self.mid.rho

    def Area(self):
        return self.A

    def Nsm(self):
        return self.nsm

    def J(self):
        return self.j

    def I11(self):
        return self.i1

    def I22(self):
        return self.i2

    def I12(self):
        return self.i12

    def E(self):
        return self.mid.E

    def G(self):
        return self.mid.G

    def Nu(self):
        return self.mid.nu

    def IAreaL(self,dim):
        if self.Type=='ROD':
            R = dim[0]
            A = pi*R**2
            Iyy = A*R**2/4.
            Izz = Iyy
            Iyz = 0.
        elif self.Type=='TUBE':
            R1 = dim[0]
            R2 = dim[1]
            A1   = pi*R1**2
            Iyy1 = A1*R1**2/4.
            A2   = pi*R2**2
            Iyy2 = A2*R2**2/4.
            A   = A1-A2
            Iyy = Iyy1-Iyy2
            Izz = Iyy
            Iyz = 0.

        elif self.Type=='I':
            h1 = dim[5]                      #        d2
            w1 = dim[2]                      # |  ------------
            y1 = dim[0]/2.-h1                # |  |    A     | d5
            sections.append([w1,h1,0.,y1])   # |  ------------
                                             # |      | |
                                             # |     >| |<--d3
            h3 = dim[4]                      # |      |B|
            w3 = dim[1]                      # | d1   | |
            y3 = -dim[0]/2.+h3               # |      | |
            sections.append([w3,h3,0.,y1])   # |   ----------
                                             # |   |   C    |  d5
            h2 = dim[0]-h1-h3                # |   ----------
            w2 = dim[3]                      #         d1
            sections.append([w2,h2,0.,0.])   #

            (A,Iyy,Izz,Iyz) = getInertiaRectangular(sections)
            assert Iyz == 0.
        else:
            raise Exception('Type=%s is not supported for %s class...' %(self.Type,self.type))
        return (A,Iyy,Izz,Iyz)
    
    def areaL(self,dim):
        """
        Area method for the PBARL and PBEAML classes (pronounced Area-L)
        @param self the object pointer
        @param dim a list of the dimensions associated with self.Type
        @retval Area of the given cross section
        """
        if   self.Type=='ROD':  A=pi*dim[0]**2
        elif self.Type=='TUBE': A=pi*(dim[0]**2-dim[1]**2)
        elif self.Type=='I':
            h1 = dim[5]
            w1 = dim[2]

            h3 = dim[4]
            w3 = dim[1]

            h2 = dim[0]-h1-h3
            w2 = dim[3]
            A = h1*w1+h2*w2+h3*w3
        elif self.Type=='CHAN':
            h1 = dim[3]
            w1 = dim[0]

            h3 = h1
            w3 = w1
            h2 = dim[1]-h1-h3
            w2 = dim[2]
            A = h1*w1+h2*w2+h3*w3
        elif self.Type=='T':
            h1 = dim[2]
            w1 = dim[0]

            h2 = dim[1]-h1
            w2 = dim[3]
            A = h1*w1+h2*w2
        elif self.Type=='BOX':
            h1 = dim[2]
            w1 = dim[0]

            h2 = dim[1]-2*h1
            w2 = dim[3]
            A = 2*(h1*w1+h2*w2)
        elif self.Type=='BAR':
            h1 = dim[1]
            w1 = dim[0]
            A = h1*w1
        elif self.Type=='CROSS':
            h1 = dim[2]
            w1 = dim[1]

            h2 = dim[3]
            w2 = dim[0]
            A = h1*w1+h2*w2
        elif self.Type=='H':
            h1 = dim[2]
            w1 = dim[1]

            h2 = dim[3]
            w2 = dim[0]
            A = h1*w1+h2*w2
        elif self.Type=='T1':
            h1 = dim[0]
            w1 = dim[2]

            h2 = dim[3]
            w2 = dim[1]
            A = h1*w1+h2*w2
        elif self.Type=='I1':
            h2 = dim[2]
            w2 = dim[1]

            h1 = dim[3]-h2
            w1 = dim[0]+w2
            A = h1*w1+h2*w2
        elif self.Type=='CHAN1':
            h2 = dim[2]
            w2 = dim[1]

            h1 = dim[3]-h2
            w1 = dim[0]+w2
            A = h1*w1+h2*w2
        elif self.Type=='Z':
            h2 = dim[2]
            w2 = dim[1]

            h1 = dim[4]-h2
            w1 = dim[0]
            A = h1*w1+h2*w2
        elif self.Type=='CHAN2':
            h2 = dim[1]
            w2 = dim[3]

            h1 = dim[2]-h2
            w1 = dim[0]*2
            A = h1*w1+h2*w2

        elif self.Type=='T2':
            h1 = dim[3]
            w1 = dim[1]

            h2 = h1-dim[2]
            w2 = dim[0]
            A = h1*w1+h2*w2
        elif self.Type=='BOX1':
            h1 = dim[2]  # top
            w1 = dim[0]
            
            h2 = sefl.dim[3] # btm
            A1 = (h1+h2)*w1
            
            h3 = dim[1]-h1-h2  # left
            w3 = dim[5]
            
            w4 = dim[4] # right
            A2 = h3*(w3+w4)
            A  = A1+A2
        elif self.Type=='HEXA':
            hBox = dim[2]
            wBox = dim[1]

            wTri = dim[0]
            A = hBox*wBox - wTri*hBox
        elif self.Type=='HAT':
            w  = dim[1]      # constant width (note h is sometimes w)
            h1 = w           # horizontal lower bar
            w1 = dim[3]
            
            h2 = dim[0]-2*w  # vertical bar
            w2 = w
            
            h3 = w           # half of top bar
            w3 = dim[2]/2.
            
            A = 2*(h1*w1+h2*w2+h3*w3)  # symmetrical box
        elif self.Type=='HAT1':
            w = dim[3]
            
            h0 = dim[4]         # btm bar
            w0 = dim[0]/2.

            h2 = dim[1]-h0-2*w  # vertical bar
            w2 = w

            h3 = w              # top bar
            w3 = dim[2]/2.

            h1 = w              # upper, horizontal lower bar (see HAT)
            w1 = w0-w3
            
            A = 2*(h0*w0+h1*w1+h2*w2+h3*w3)
            
        elif self.Type=='DBOX':
            #
            #  |--2------5----
            #  |     |       |
            #  1     3       6
            #  |     |       |
            #  |--4--|---7---|
            #
            
            0,1,2,6,11
            1,2,3,7,12
            
            hTotal = dim[11]
            wTotal = dim[0]

            h2 = dim[6]
            w2 = dim[3]

            h4 = dim[7]
            w4 = w2

            h1 = hTotal-h2-h4
            w1 = dim[3]
            
            h5 = dim[8]
            w5 = wTotal-w2

            h7 = dim[9]
            w7 = w5

            h6 = hTotal-h5-h7
            w6 = dim[5]

            h3 = (h1+h6)/2.
            w3 = dim[4]
            A = h1*w1 +h2*w2 +h3*w3 +h4*w4 +h5*w5 +h6*w6 +h7*w7 +h8*w8
        else:
            raise Exception('Type=%s is not supported for %s class...' %(self.Type,self.type))
            
        return A

class IntegratedLineProperty(LineProperty):
    type = 'IntegratedLineProperty'
    def __init__(self,card,data):
        LineProperty.__init__(self,card,data)

    def Area(self):
        A = integratePositiveLine(self.xxb,self.A)
        return A

    def J(self):
        J = integratePositiveLine(self.xxb,self.j)
        return J

    def I11(self):
        i1 = integratePositiveLine(self.xxb,self.i1)
        return i1

    def I22(self):
        i2 = integratePositiveLine(self.xxb,self.i2)
        return i2

    def I12(self):
        i12 = integrateLine(self.xxb,self.i12)
        return i12

    def Nsm(self):
        nsm = integratePositiveLine(self.xxb,self.nsm)
        return nsm

        
class PROD(LineProperty):
    type = 'PROD'
    def __init__(self,card=None,data=None):
        LineProperty.__init__(self,card,data)

        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.A   = card.field(3)
            self.j   = card.field(4)
            self.c   = card.field(5,0.0)
            self.nsm = card.field(6,0.0)
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.A   = data[2]
            self.j   = data[3]
            self.c   = data[4]
            self.nsm = data[5]
        ###

    def crossReference(self,mesh):
        self.mid = mesh.Material(self.mid)

    def I11(self):
        """@todo whats the proper formula to use for a ROD"""
        return self.j/2.

    def I22(self):
        """@todo whats the proper formula to use for a ROD"""
        return self.j/2.

    def I12(self):
        return 0.

    def __repr__(self):
        c   = self.setBlankIfDefault(self.c,0.0)
        nsm = self.setBlankIfDefault(self.nsm,0.0)
        fields = ['PROD',self.pid,self.Mid(),self.A,self.j,c,nsm]
        return self.printCard(fields)

class PTUBE(LineProperty):
    type = 'PTUBE'
    def __init__(self,card=None,data=None):
        LineProperty.__init__(self,card,data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.OD1 = card.field(3)
            self.t   = card.field(4,self.outerDiameter/2.)
            self.nsm = card.field(5,0.0)
            self.OD2 = card.field(6,self.outerDiameter)
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.OD1 = data[2]
            self.t   = data[3]
            self.nsm = data[4]
            self.OD2 = self.OD1
            #self.OD2 = data[5]  ## @note quirk to this one...

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        t   = self.setBlankIfDefault(self.t,self.OD1/2.)
        nsm = self.setBlankIfDefault(self.nsm,0.0)
        OD2 = self.setBlankIfDefault(self.OD2,self.OD1)
        fields = ['PTUBE',self.pid,self.Mid(),self.OD1,t,nsm,OD2]
        return self.printCard(fields)
    
    def Area(self):
        A = (self.area1()+self.area2())/2.
        
        #A1 = pi*D1^2/4 - pi*((D1-2t)^2)/4
        #A2 = pi*D2^2/4 - pi*((D2-2t)^2)/4
        #A = A1+A2
        
        #A = pi*D1*t/2 + pi*D2*t/2 - pi*t
        #A = pi*t*(D1/2 + D2/2 - t)
        #A = pi*t*( (D1+D2)/2.-t )
        
        A2 = pi*t*( (D1+D2)/2.-t )
        
        assert A == A2,'AREA method has problem in PTUBE Aold=%s Anew=%s' %(A,A2)
        return A2
        
    def area1(self):
        """@todo remove after verifying formula..."""
        Dout = self.OD
        Din  = Dout-2*self.t
        A = pi()/4.*(Dout*Dout-Din*Din)
        return A1

    def area2(self):
        """@todo remove after verifying formula..."""
        Dout = self.OD2
        Din  = Dout-2*self.t
        A = pi()/4.*(Dout*Dout-Din*Din)
        return A2

    def massMatrix(self):
        """@todo not done"""
        m = zeros(6,6)
        m[0,0] = 1.
        return m

class PBAR(LineProperty):
    type = 'PBAR'
    def __init__(self,card=None,data=None):
        """
        @todo
            support solution 600 default
            do a check for mid -> MAT1      for structural
            do a check for mid -> MAT4/MAT5 for thermal
        """
        LineProperty.__init__(self,card,data)
        if card:
            ## property ID -> use Pid()
            self.pid = card.field(1)
            ## material ID -> use Mid()
            self.mid = card.field(2)
            ## Area -> use Area()
            self.A   = card.field(3,0.0)
            ## Izz -> use Izz()
            self.i1  = card.field(4,0.0)
            ## Iyy -> use Iyy()
            self.i2  = card.field(5,0.0)
            ## Polar Moment of Inertia J -> use J()
            self.j   = card.field(6,0.0) # default=1/2(I1+I2) for SOL=600, otherwise 0.0
            ## nonstructral mass -> use Nsm()
            self.nsm = card.field(7,0.0)

            self.C1  = card.field(9, 0.0)
            self.C2  = card.field(10,0.0)
            self.D1  = card.field(11,0.0)
            self.D2  = card.field(12,0.0)
            self.E1  = card.field(13,0.0)
            self.E2  = card.field(14,0.0)
            self.F1  = card.field(15,0.0)
            self.F2  = card.field(16,0.0)

            self.K1  = card.field(17,1e8) ## default=infinite
            self.K2  = card.field(18,1e8) ## default=infinite
            self.i12 = card.field(19,0.0)
            if self.A==0.0:
                assert self.K1==None
                assert self.K2==None
            ###
        else:
            self.mid = data[0]
            self.pid = data[1]
            self.A   = data[2]
            self.i1  = data[3]
            self.i2  = data[4]
            self.j   = data[5]

            self.nsm = data[6]
            #self.fe  = data[7] ## @todo not documented....
            self.C1  = data[8]
            self.C2  = data[9]
            self.D1  = data[10]
            self.D2  = data[11]
            self.E1  = data[12]
            self.E2  = data[13]
            self.F1  = data[14]
            self.F2  = data[15]
            self.K1  = data[16]
            self.K2  = data[17]
            self.i12 = data[18]
        ###

    def MassPerLength(self):
        """
        \f[ \frac{m}{L} = \rho A+nsm \f]
        """
        A   = self.Area()
        rho = self.Rho()
        nsm = self.Nsm()
        return rho*A+nsm

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def Area(self):
        return self.A

    #def Nsm(self):
    #    return self.nsm

    #def J(self):
    #    return self.j

    #def Izz(self):
    #    return self.i1

    #def Iyy(self):
    #    return self.i2

    #def Iyz(self):
    #    return self.i12

    def __repr__(self):
        A   = self.setBlankIfDefault(self.A,0.0)
        i1  = self.setBlankIfDefault(self.i1,0.0)
        i2  = self.setBlankIfDefault(self.i2,0.0)
        i12 = self.setBlankIfDefault(self.i12,0.0)
        j   = self.setBlankIfDefault(self.j,0.0)
        nsm = self.setBlankIfDefault(self.nsm,0.0)
        
        C1  = self.setBlankIfDefault(self.C1,0.0)
        C2  = self.setBlankIfDefault(self.C2,0.0)

        D1  = self.setBlankIfDefault(self.D1,0.0)
        D2  = self.setBlankIfDefault(self.D2,0.0)

        E1  = self.setBlankIfDefault(self.E1,0.0)
        E2  = self.setBlankIfDefault(self.E2,0.0)

        F1  = self.setBlankIfDefault(self.F1,0.0)
        F2  = self.setBlankIfDefault(self.F2,0.0) # must have 1 on line, if line3 is not empty
        
        line3 = [self.K1,self.K2,i12]
        #print "line3 = ",line3
        
        line1 = ['PBAR',self.pid,self.Mid(),self.A,i1,i2,j,nsm,None]

        if line3==[None,None,None]:
            line2 = [C1,C2,D1,D2,E1,E2,F1,F2]
        else:
            line2 = [C1,C2,D1,D2,E1,E2,F1,self.F2]
        fields = line1+line2+line3

        return self.printCard(fields)

class PBARL(LineProperty):
    """
    @warning doesnt support user-defined types
    """
    type = 'PBARL'
    validTypes = {
        "ROD"   : 1, 
        "TUBE"  : 2,
        "I"     : 6,
        "CHAN"  : 4,
        "T"     : 4,
        "BOX"   : 4,
        "BAR"   : 2,
        "CROSS" : 4,
        "H"     : 4,
        "T1"    : 4,
        "I1"    : 4,
        "CHAN1" : 4,
        "Z"     : 4,
        "CHAN2" : 4,
        "T2"    : 4,
        "BOX1"  : 6,
        "HEXA"  : 3,
        "HAT"   : 4,
        "HAT1"  : 5,
        "DBOX"  : 10, # was 12
        } # for GROUP="MSCBML0"

    def __init__(self,card=None,data=None):
        LineProperty.__init__(self,card,data)
        if card:
            self.pid   = card.field(1)
            self.mid   = card.field(2)
            self.group = card.field(3,'MSCBMLO')
            self.Type  = card.field(4)

            #nDim = len(self.dim)-1
            nDim = self.validTypes[self.Type]
            self.dim = card.fields(9,9+nDim+1)
            self.nsm = card.field(9+nDim+1,0.0)

            #self.dim = fields(9)
            if nDim>0:
                self.nsm = self.setDefaultIfBlank(self.dim.pop(),0.0)
            else:
                self.nsm = 0.0
            ###
            assert isinstance(self.nsm,float)
        else:
            self.pid   = data[0]
            self.mid   = data[1]
            self.group = data[2].strip()
            self.Type  = data[3].strip()
            self.dim   = list(data[4:-1])
            self.nsm   = data[-1]
            #print "group = |%s|" %(self.group)
            #print "Type  = |%s|" %(self.Type)
            #print "dim = ",self.dim
            #print str(self)
            #print "*PBARL = ",data
            #raise Exception('not finished...')
        assert self.Type in self.validTypes,'Invalid PBARL Type, Type=%s validTypes=%s' %(self.Type,self.validTypes.keys())
        assert len(self.dim)==self.validTypes[self.Type],'dim=%s len(dim)=%s Type=%s len(dimType)=%s' %(self.dim,len(self.dim),self.Type,self.validTypes[self.Type])
        assert None not in self.dim
            
    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def Area(self):
        return self.areaL(self.dim)

    def Nsm(self):
        return self.nsm
    
    def MassPerLength(self):
        """
        mass = L*(Area*rho+nsm)
        mass/L = Area*rho+nsm
        """
        rho  = self.Rho()
        area = self.Area()
        nsm  = self.Nsm()
        return area*rho+nsm

    def __repr__(self):
        group = self.setBlankIfDefault(self.group,'MSCBMLO')
        fields = ['PBARL',self.pid,self.Mid(),group,self.Type,None,None,None,None,
        ]+self.dim+[self.nsm]
        return self.printCard(fields)

class PBEAM(IntegratedLineProperty):
    type = 'PBEAM'
    def __init__(self,card=None,data=None):
        """
        @todo fix 0th entry of self.so, self.xxb
        """
        IntegratedLineProperty.__init__(self,card,data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            #print "pid = ",self.pid

            nFields = card.nFields()-1 # -1 for PBEAM field
            fields = card.fields()
            #print "  fields = ",fields

            self.so  = ['YES'] ## @todo what are these values (so[0],xxb[0])???
            self.xxb = [0.]
            self.A   = [card.field(3) ]
            self.i1  = [card.field(4,0.0) ]
            self.i2  = [card.field(5,0.0) ]
            self.i12 = [card.field(6,0.0) ]
            self.j   = [card.field(7,0.0) ]
            self.nsm = [card.field(8,0.0) ]
            self.c1  = [card.field(9,0.0) ]
            self.c2  = [card.field(10,0.0)]
            self.d1  = [card.field(11,0.0)]
            self.d2  = [card.field(12,0.0)]
            self.e1  = [card.field(13,0.0)]
            self.e2  = [card.field(14,0.0)]
            self.f1  = [card.field(15,0.0)]
            self.f2  = [card.field(16,0.0)]

            #fields = card.fields(0)
            #print "fieldsPBEAM = ",fields
            #fieldsMid = fields[16:]
            #print "fieldsMid = ",fieldsMid

            #fields = card.fields(9)
            #print ""
            #print "  nFields = ",nFields
            #nFields = card.nFields()-16 # 17+16 (leading + trailing fields)
            # counting continuation cards
            nMajor    = nFields/16
            nLeftover = nFields%16
            #print "  nMajor=%s nLeftover=%s" %(nMajor,nLeftover)
            if nLeftover==0:
                nMajor-=1

            if nMajor==0:
                nMajor=1

            x = (nMajor)*16+1
            if card.field(x) in ['YES','YESA','NO']: # there is no footer
                nMajor +=1
                x+=16

            #print "  nMajor = ",nMajor
            for nRepeated in range(1,nMajor):
                #print "  adding a major"
                nStart = nRepeated*16+1  # field 17 is the first possible so
                propFields = card.fields(nStart,nStart+16)
                #print "propFields = ",propFields

                #print "  so = ",propFields[0]
                assert propFields[0] in [None,'YES','YESA','NO'],"SO=%s for PBEAM pid=%s" %(propFields[0],self.pid)
                self.so.append( propFields[0])
                self.xxb.append(propFields[1])
                self.A.append(  propFields[2])
                self.i1.append( self.setDefaultIfBlank(propFields[3],0.0))
                self.i2.append( self.setDefaultIfBlank(propFields[4],0.0))
                self.i12.append(self.setDefaultIfBlank(propFields[5],0.0))
                self.j.append(  self.setDefaultIfBlank(propFields[6],0.0))
                self.nsm.append(self.setDefaultIfBlank(propFields[7],0.0))
                self.c1.append( self.setDefaultIfBlank(propFields[8],0.0))
                self.c2.append( self.setDefaultIfBlank(propFields[9],0.0))
                self.d1.append( self.setDefaultIfBlank(propFields[10],0.0))
                self.d2.append( self.setDefaultIfBlank(propFields[11],0.0))
                self.e1.append( self.setDefaultIfBlank(propFields[12],0.0))
                self.e2.append( self.setDefaultIfBlank(propFields[13],0.0))
                self.f1.append( self.setDefaultIfBlank(propFields[14],0.0))
                self.f2.append( self.setDefaultIfBlank(propFields[15],0.0))
            #print "nRepeated = ",nRepeated

            # footer fields
            self.k1   = card.field(x,1.0)

            assert self.k1 not in ['YES','YESA','NO'],'error reading PBEAM card pid=%s' %(self.pid)
            #print "  k1 = ",self.k1
            self.k2   = card.field(x+1,1.0)
            self.s1   = card.field(x+2,0.0)
            self.s2   = card.field(x+3,0.0)
            self.nsia = card.field(x+4,0.0)
            self.nsib = card.field(x+5,self.nsia)
            self.cwa  = card.field(x+6,0.0)
            self.cwb  = card.field(x+7,self.cwa)

            self.m1a = card.field(x+8,0.0)
            self.m2a = card.field(x+9,self.m1a)
            self.m1b = card.field(x+10,0.0)
            self.m2b = card.field(x+11,self.m1b)
            self.n1a = card.field(x+12,0.0)
            self.n2a = card.field(x+13,self.n1a)
            self.n1b = card.field(x+14,0.0)
            self.n2b = card.field(x+15,self.n1b)
        else:
            raise Exception('not supported')
        ###

    #def Area(self):
    #    """@warning area field not supported fully on PBEAM card"""
    #    #raise Exception(self.A[0])
    #    return self.A[0]

    #def Nsm(self):
    #    """@warning nsm field not supported fully on PBEAM card"""
    #    #raise Exception(self.nsm[0])
    #    return self.nsm[0]

    def MassPerLength(self):
        """
        mass = L*(Area*rho+nsm)
        mass/L = Area*rho+nsm
        """
        rho  = self.Rho()
        massPerLs = []
        for area,nsm in zip(self.A,self.nsm):
            massPerLs.append(area*rho+nsm)
        massPerL = integratePositiveLine(self.xxb,massPerLs)
        return massPerL

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        fields = ['PBEAM',self.pid,self.Mid()]
        #print "fieldsA = ",fields
        
        #print len(self.so)
        i = 0
        #print "pid=%s" %(self.pid)
        for (so,xxb,A,i1,i2,i12,j,nsm,c1,c2,d1,d2,e1,e2,f1,f2) in zip(
            self.so,self.xxb,self.A,self.i1,self.i2,self.i12,self.j,self.nsm,
            self.c1,self.c2,self.d1,self.d2,self.e1,self.e2,self.f1,self.f2):

            i1  = self.setBlankIfDefault(i1,0.0)
            i2  = self.setBlankIfDefault(i2,0.0)
            i12 = self.setBlankIfDefault(i12,0.0)
            j   = self.setBlankIfDefault(j,0.0)

            nsm = self.setBlankIfDefault(nsm,0.0)
            #c1 = self.setBlankIfDefault(c1,0.0)
            d1 = self.setBlankIfDefault(d1,0.0)
            e1 = self.setBlankIfDefault(e1,0.0)
            f1 = self.setBlankIfDefault(f1,0.0)

            c2 = self.setBlankIfDefault(c2,0.0)
            d2 = self.setBlankIfDefault(d2,0.0)
            e2 = self.setBlankIfDefault(e2,0.0)
            f2 = self.setBlankIfDefault(f2,0.0)

            #print "  i = ",i
            if i==0:
                fields +=        [A,i1,i2,i12,j,nsm,c1,c2,d1,d2,e1,e2,f1,f2] # the 1st 2 fields aren't written
            else:
                fields += [so,xxb,A,i1,i2,i12,j,nsm,c1,c2,d1,d2,e1,e2,f1,f2]
            i+=1
        k1 = self.setBlankIfDefault(self.k1,1.0)
        k2 = self.setBlankIfDefault(self.k2,1.0)
        s1 = self.setBlankIfDefault(self.s1,0.0)
        s2 = self.setBlankIfDefault(self.s2,0.0)
        
        nsia = self.setBlankIfDefault(self.nsia,0.0)
        nsib = self.setBlankIfDefault(self.nsib,self.nsia)

        cwa = self.setBlankIfDefault(self.cwa,0.0)
        cwb = self.setBlankIfDefault(self.cwb,self.cwa)

        m1a = self.setBlankIfDefault(self.m1a,0.0)
        m2a = self.setBlankIfDefault(self.m2a,self.m1a)
        m1b = self.setBlankIfDefault(self.m1b,0.0)
        m2b = self.setBlankIfDefault(self.m2b,self.m1b)

        #print "m1a=%s m2a=%s" %(m1a,m2a)
        #print "m1b=%s m2b=%s" %(m1b,m2b)

        n1a = self.setBlankIfDefault(self.n1a,0.0)
        n2a = self.setBlankIfDefault(self.n2a,self.n1a)
        n1b = self.setBlankIfDefault(self.n1b,0.0)
        n2b = self.setBlankIfDefault(self.n2b,self.n1b)
        #print "n1a=%s n2a=%s" %(n1a,n2a)
        #print "n1b=%s n2b=%s" %(n1b,n2b)

        footer = [k1,k2,s1,s2,nsia,nsib,cwa,cwb,
                  m1a,m2a,m1b,m2b,n1a,n2a,n1b,n2b]
        
        #if footer == [self.k1,None,None,None,None,None,None,None,   None,None,None,None,None,None,None,None]:
        #    footer = []
        fields+=footer
        #print fields
        return self.printCard(fields)
        
class PBEAML(IntegratedLineProperty):
    type = 'PBEAML'
    validTypes = {
        "ROD"   : 1, 
        "TUBE"  : 2,
        "I"     : 6,
        "CHAN"  : 4,
        "T"     : 4,
        "BOX"   : 4,
        "BAR"   : 2,
        "CROSS" : 4,
        "H"     : 4,
        "T1"    : 4,
        "I1"    : 4,
        "CHAN1" : 4,
        "Z"     : 4,
        "CHAN2" : 4,
        "T2"    : 4,
        "BOX1"  : 6,
        "HEXA"  : 3,
        "HAT"   : 4,
        "HAT1"  : 5,
        "DBOX"  : 10, # was 12
        } # for GROUP="MSCBML0"
    def __init__(self,card=None,data=None):
        IntegratedLineProperty.__init__(self,card,data)
        if card:
            self.pid   = card.field(1)
            self.mid   = card.field(2)
            self.group = card.field(3,'MSCBMLO')
            self.Type  = card.field(4)

            nDim = self.validTypes[self.Type]
            nAll = nDim+1

            #nDim = len(self.dim)-1
            dimAll = card.fields(9)
            
            self.dim = []
            Dim = []
            j=0
            self.xxb = [0.]
            self.so = ['YES']
            self.nsm = []
            
            n = 0
            i=0
            for i,dim in dimAll:
                if j<nDim:
                    Dim.append(dim)
                    j+=1
                else:
                    self.nsm.append(dim)
                    if n>0:
                        self.so.append( card.field(i+1,'YES')) # dimAll[i+1]
                        self.xxb.append(card.field(i+2,1.0  ))  #dimAll[i+2]
                    j = 0
                    n+=1
                    i+=2
                    self.dim.append(Dim)
                    Dim = []
                ###
            ###
            if j<nDim: # if the last field is blank
                self.dim.append(Dim)
                self.nsm.append(card.field(i,0.0))
            ###
    
    def MassPerLength(self):
        """
        mass = L*(Area*rho+nsm)
        mass/L = Area*rho+nsm
        """
        rho  = self.Rho()
        massPerLs = []
        for dim,n in zip(self.dim,self.nsm):
            a = self.areaL(dim)
            massPerLs.append(a*rho+nsm)
        massPerL = integratePositiveLine(self.xxb,massPerLs)
        return massPerL

    def Area(self):
        Areas = []
        for dim in self.dim:
            Areas.append( self.areaL(dim) )
        A = integrateLine(self.xxb,Areas)
        return A

    #def Mid(self):
    #    return self.mid

    def crossReference(self,model):
        """
        @warning For structural problems, PBEAML entries must reference a MAT1 material entry
        @warning For heat-transfer problems, the MID must reference a MAT4 or MAT5 material entry.
        @note what happens when there are 2 subcases?
        """
        self.mid = model.Material(self.mid)
    
    def verify(self,model,iSubcase):
        if self.IsThermalSolution(iSubcase):
            assert self.mid.type in ['MAT4','MAT5']
        else:
            assert self.mid.type in ['MAT1']
        ###

    def __repr__(self):
        fields = ['PBEAML',self.Pid,self.Mid(),self.group,self.Type,None,None,None,None]
        
        for i,(xxb,so,dim) in enumerate(zip(self.xxb,self.so,self.dim)):
            if i==0:
                fields += dim
            else:
                fields += [xxb,so]+dim
            ###
        ###
        print self.printCard(fields)
        
        raise Exception('verify PBEAML...')
        return self.printCard(fields)

class PBEAM3(LineProperty): # not done, cleanup
    type = 'PBEAM3'
    def __init__(self,card=None,data=None):
        LineProperty.__init__(self,card,data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)

            self.A   = card.field(3)
            self.iz  = card.field(4)
            self.iy  = card.field(5)
            self.iyz = card.field(6,0.0)
            self.j   = card.field(7,self.iy+self.iz)
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
        ###
        else:
            raise Exception('not implemented...')
        ###

    def Nsm(self):
        """@warning nsm field not supported fully on PBEAM3 card"""
        #raise Exception(self.nsm[0])
        return self.nsm

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        raise Exception('not done...')
        fields = ['PBEAM3',self.pid,self.Mid(),] # other
        return self.printCard(fields)
