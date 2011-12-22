import sys
from numpy import zeros,pi

# my code
from ..baseCard import Property

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

    def Nsm(self):
        return self.nsm

    def E(self):
        return self.mid.E

    def G(self):
        return self.mid.G

    def nu(self):
        return self.mid.nu

class PROD(LineProperty):
    type = 'PROD'
    def __init__(self,card=None,data=None):
        LineProperty.__init__(self,card,data)

        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.A   = card.field(3)
            self.J   = card.field(4)
            self.c   = card.field(5,0.0)
            self.nsm = card.field(6,0.0)
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.A   = data[2]
            self.J   = data[3]
            self.c   = data[4]
            self.nsm = data[5]
        ###

    def crossReference(self,mesh):
        self.mid = mesh.Material(self.mid)

    def __repr__(self):
        c   = self.setBlankIfDefault(self.c,0.0)
        nsm = self.setBlankIfDefault(self.nsm,0.0)
        fields = ['PROD',self.pid,self.Mid(),self.A,self.J,c,nsm]
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
    def __init__(self,card=None,data=None):
        """
        @todo
            support solution 600 default
            do a check for mid -> MAT1      for structural
            do a check for mid -> MAT4/MAT5 for thermal
        """
        LineProperty.__init__(self,card,data)
        if card:
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

            self.K1  = card.field(17,1e8) ## default=infinite
            self.K2  = card.field(18,1e8) ## default=infinite
            self.I12 = card.field(19,0.0)
            if self.A==0.0:
                assert self.K1==None
                assert self.K2==None
            ###
        else:
            self.mid = data[0]
            self.pid = data[1]
            self.A   = data[2]
            self.I1  = data[3]
            self.I2  = data[4]
            self.J   = data[5]
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
            self.I12 = data[18]
        ###

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def Area(self):
        return self.A

    def Nsm(self):
        return self.nsm

    def J(self):
        return self.J

    def Izz(self):
        return self.I1

    def Iyy(self):
        return self.I2

    def Iyz(self):
        return self.I12

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
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.group = card.field(3,'MSCBMLO')
            self.Type = card.field(4)

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
        if   self.Type=='ROD':  A=pi*self.dim[0]**2
        elif self.Type=='TUBE': A=pi*(self.dim[0]**2-self.dim[1]**2)
        elif self.Type=='I':
            h1 = self.dim[5]
            w1 = self.dim[2]

            h3 = self.dim[4]
            w3 = self.dim[1]

            h2 = self.dim[0]-h1-h3
            w2 = self.dim[3]
            A = h1*w1+h2*w2+h3*w3
        elif self.Type=='CHAN':
            h1 = self.dim[3]
            w1 = self.dim[0]

            h3 = h1
            w3 = w1
            h2 = self.dim[1]-h1-h3
            w2 = self.dim[2]
            A = h1*w1+h2*w2+h3*w3
        elif self.Type=='T':
            h1 = self.dim[2]
            w1 = self.dim[0]

            h2 = self.dim[1]-h1
            w2 = self.dim[3]
            A = h1*w1+h2*w2
        elif self.Type=='BOX':
            h1 = self.dim[2]
            w1 = self.dim[0]

            h2 = self.dim[1]-2*h1
            w2 = self.dim[3]
            A = 2*(h1*w1+h2*w2)
        elif self.Type=='BAR':
            h1 = self.dim[1]
            w1 = self.dim[0]
            A = h1*w1
        elif self.Type=='CROSS':
            h1 = self.dim[2]
            w1 = self.dim[1]

            h2 = self.dim[3]
            w2 = self.dim[0]
            A = h1*w1+h2*w2
        elif self.Type=='H':
            h1 = self.dim[2]
            w1 = self.dim[1]

            h2 = self.dim[3]
            w2 = self.dim[0]
            A = h1*w1+h2*w2
        elif self.Type=='T1':
            h1 = self.dim[0]
            w1 = self.dim[2]

            h2 = self.dim[3]
            w2 = self.dim[1]
            A = h1*w1+h2*w2
        elif self.Type=='I1':
            h2 = self.dim[2]
            w2 = self.dim[1]

            h1 = self.dim[3]-h2
            w1 = self.dim[0]+w2
            A = h1*w1+h2*w2
        elif self.Type=='CHAN1':
            h2 = self.dim[2]
            w2 = self.dim[1]

            h1 = self.dim[3]-h2
            w1 = self.dim[0]+w2
            A = h1*w1+h2*w2
        elif self.Type=='Z':
            h2 = self.dim[2]
            w2 = self.dim[1]

            h1 = self.dim[4]-h2
            w1 = self.dim[0]
            A = h1*w1+h2*w2
        elif self.Type=='CHAN2':
            h2 = self.dim[1]
            w2 = self.dim[3]

            h1 = self.dim[2]-h2
            w1 = self.dim[0]*2
            A = h1*w1+h2*w2

        elif self.Type=='T2':
            h1 = self.dim[3]
            w1 = self.dim[1]

            h2 = h1-self.dim[2]
            w2 = self.dim[0]
            A = h1*w1+h2*w2
        elif self.Type=='BOX1':
            h1 = self.dim[2]  # top
            w1 = self.dim[0]
            
            h2 = sefl.dim[3] # btm
            A1 = (h1+h2)*w1
            
            h3 = self.dim[1]-h1-h2  # left
            w3 = self.dim[5]
            
            w4 = self.dim[4] # right
            A2 = h3*(w3+w4)
            A  = A1+A2
        elif self.Type=='HEXA':
            hBox = self.dim[2]
            wBox = self.dim[1]

            wTri = self.dim[0]
            A = hBox*wBox - wTri*hBox
        elif self.Type=='HAT':
            w = self.dim[1] # constant width (note h is sometimes w)
            h1 = w                # horizontal lower bar
            w1 = self.dim[3]
            
            h2 = self.dim[0]-2*w  # vertical bar
            w2 = w
            
            h3 = w                # half of top bar
            w3 = self.dim[2]/2.
            
            A = 2*(h1*w1+h2*w2+h3*w3)  # symmetrical box
        elif self.Type=='HAT1':
            w = self.dim[3]
            
            h0 = self.dim[4]         # btm bar
            w0 = self.dim[0]/2.

            h2 = self.dim[1]-h0-2*w  # vertical bar
            w2 = w

            h3 = w              # top bar
            w3 = self.dim[2]/2.

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
            
            hTotal = self.dim[11]
            wTotal = self.dim[0]

            h2 = self.dim[6]
            w2 = self.dim[3]

            h4 = self.dim[7]
            w4 = w2

            h1 = hTotal-h2-h4
            w1 = self.dim[3]
            
            h5 = self.dim[8]
            w5 = wTotal-w2

            h7 = self.dim[9]
            w7 = w5

            h6 = hTotal-h5-h7
            w6 = self.dim[5]

            h3 = (h1+h6)/2.
            w3 = self.dim[4]
            A = h1*w1 +h2*w2 +h3*w3 +h4*w4 +h5*w5 +h6*w6 +h7*w7 +h8*w8
        else:
            raise Exception('Type=%s is not supported...' %(self.Type))
            
        return A
        
    def Nsm(self):
        return self.nsm

    def __repr__(self):
        group = self.setBlankIfDefault(self.group,'MSCBMLO')
        fields = ['PBARL',self.pid,self.Mid(),group,self.Type,None,None,None,None,
        ]+self.dim+[self.nsm]
        return self.printCard(fields)

class PBEAM(LineProperty):
    type = 'PBEAM'
    def __init__(self,card=None,data=None):
        """
        @todo fix 0th entry of self.so, self.xxb
        """
        LineProperty.__init__(self,card,data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            #print "pid    = ",self.pid

            nFields = card.nFields()-1 # -1 for PBEAM field
            fields = card.fields()
            #print "  fields = ",fields

            self.so  = [None] ## @todo what are these values???
            self.xxb = [None]
            self.A   = [card.field(3) ]
            self.i1  = [card.field(4,0.0) ]
            self.i2  = [card.field(5,0.0) ]
            self.i12 = [card.field(6,0.0) ]
            self.J   = [card.field(7,0.0) ]
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
                self.J.append(  self.setDefaultIfBlank(propFields[6],0.0))
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

    def Area(self):
        """@warning area field not supported fully on PBEAM card"""
        #raise Exception(self.A[0])
        return self.A[0]

    def Nsm(self):
        """@warning nsm field not supported fully on PBEAM card"""
        #raise Exception(self.nsm[0])
        return self.nsm[0]

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        fields = ['PBEAM',self.pid,self.Mid()]
        #print "fieldsA = ",fields
        
        #print len(self.so)
        i = 0
        #print "pid=%s" %(self.pid)
        for (so,xxb,A,i1,i2,i12,J,nsm,c1,c2,d1,d2,e1,e2,f1,f2) in zip(
            self.so,self.xxb,self.A,self.i1,self.i2,self.i12,self.J,self.nsm,
            self.c1,self.c2,self.d1,self.d2,self.e1,self.e2,self.f1,self.f2):

            i1  = self.setBlankIfDefault(i1,0.0)
            i2  = self.setBlankIfDefault(i2,0.0)
            i12 = self.setBlankIfDefault(i12,0.0)
            J   = self.setBlankIfDefault(J,0.0)

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
                fields +=        [A,i1,i2,i12,J,nsm,c1,c2,d1,d2,e1,e2,f1,f2] # the 1st 2 fields aren't written
            else:
                fields += [so,xxb,A,i1,i2,i12,J,nsm,c1,c2,d1,d2,e1,e2,f1,f2]
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
        
#class PBEAML(LineProperty): #not done
class PBEAM3(LineProperty): # not done, cleanup
    type = 'PBEAM3'
    def __init__(self,card=None,data=None):
        LineProperty.__init__(self,card,data)
        if card:
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
