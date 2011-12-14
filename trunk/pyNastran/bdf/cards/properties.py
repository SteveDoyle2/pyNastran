#import sys
#from numpy import array,cross,dot
#from numpy import array
#from numpy.linalg import norm
from numpy import zeros

# my code
from baseCard import Property

class PointProperty(Property):
    type = 'PointProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass

class NSM(PointProperty):
    """
    Defines a set of non structural mass.
    """
    ## Set points to either Property entries or Element entries. Properties are:
    validProperties = [
        'PSHELL', 'PCOMP', 'PBAR',  'PBARL', 'PBEAM',  'PBEAML', 'PBCOMP', 'PROD',
        'CONROD', 'PBEND', 'PSHEAR','PTUBE', 'PCONEAX','PRAC2D']
    def __init__(self,card=None,nOffset=0,data=None):
        #Element.__init__(self,card,data)
        nOffset *= 2
        if card:
            self.sid   = card.field(1)
            self.Type  = card.field(2)
            self.id    = card.field(3+nOffset)
            self.value = card.field(4+nOffset)
        else:
            self.sid   = data[0]
            self.Type  = data[1]
            self.id    = data[2]
            self.value = data[3]
        ###
    def __repr__(self):
        nodes = self.nodeIDs()
        fields = ['NSM',self.sid,self.Type,self.id,self.value]
        return self.printCard(fields)

class PMASS(PointProperty):
    def __init__(self,card=None,nOffset=0,data=None):
        PointProperty.__init__(self,card,nOffset=0)
        
        nOffset *= 2
        self.pid  = card.field(1+nOffset)
        self.mass = card.field(2+nOffset,0.)

    def Mass(self):
        return self.mass

    def __repr__(self):
        fields = ['PMASS',self.pid,self.Mass]

class SpringProperty(Property):
    type = 'SpringProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass

class PELAS(SpringProperty):
    type = 'PELAS'
    def __init__(self,card=None,nPELAS=0,data=None):
        SpringProperty.__init__(self,card,data)
        nOffset = nPELAS*5
        if card:
            self.pid = card.field(1+nOffset) # 2 PELAS properties can be defined on 1 PELAS card
            self.k   = card.field(2+nOffset) # these are split into 2 separate cards
            self.ge  = card.field(3+nOffset)
            self.s   = card.field(4+nOffset)
        else:
            self.pid = data[0]
            self.k   = data[1]
            self.ge  = data[2]
            self.s   = data[3]
        ###

    def __repr__(self):
        fields = ['PELAS',self.pid,self.k,self.ge,self.s]
        return self.printCard(fields)

class DamperProperty(Property):
    type = 'DamperProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass

class PDAMP(DamperProperty):
    type = 'PDAMP'
    def __init__(self,card=None,nPDAMP=0,data=None):
        DamperProperty.__init__(self,card,data)
        nOffset = nPDAMP*2
        if card:
            ## Property ID
            self.pid = card.field(1+nOffset) # 3 PDAMP properties can be defined on 1 PDAMP card
            ## Force per unit velocity (Real)
            self.b   = card.field(2+nOffset) # these are split into 2 separate cards
        else:
            self.pid = data[0]
            self.b   = data[1]
        ###

    def __repr__(self):
        fields = ['PDAMP',self.pid,self.b]
        return self.printCard(fields)


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

            self.K1  = card.field(17)
            self.K2  = card.field(18)
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
    doesnt support user-defined types
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
            print "group = |%s|" %(self.group)
            print "Type  = |%s|" %(self.Type)
            print "dim = ",self.dim
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
        #elif self.Type=='Z':
        #elif self.Type=='CHAN2':
        elif self.Type=='T2':
            h1 = self.dim[3]
            w1 = self.dim[1]

            h2 = h1-self.dim[2]
            w2 = self.dim[0]
            A = h1*w1+h2*w2
        #elif self.Type=='BOX1':
        elif self.Type=='HEXA':
            hBox = self.dim[2]
            wBox = self.dim[1]

            wTri = self.dim[0]
            A = hBox*wBox - wTri*hBox
        #elif self.Type=='HAT':
        #elif self.Type=='HAT1':
        #elif self.Type=='DBOX':
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
        else:
            #print "len(data) = ",len(data)
            self.pid  = data[0]
            self.z0   = data[1]
            self.nsm  = data[2]
            self.sb   = data[3]
            self.ft   = data[4]
            self.TRef = data[5]
            self.ge   = data[6]

            self.lam  = None ## @todo No Symmetry - check if this is correct

            Mid   = data[7]
            T     = data[8]
            Theta = data[9]
            Sout  = data[10]

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

class SolidProperty(Property):
    type = 'SolidProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass

class PLSOLID(SolidProperty):
    """
    Defines a fully nonlinear (i.e., large strain and large rotation) hyperelastic solid
    element.
    PLSOLID PID MID STR
    PLSOLID 20 21
    """
    type = 'PLSOLID'
    def __init__(self,card=None,data=None):
        SolidProperty.__init__(self,card,data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.ge  = card.field(3)
            self.str = card.field(4,'GRID')
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.ge  = data[2]
            self.str = data[3]
            print "data = ",data
        ###
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

class PSOLID(SolidProperty):
    """
    PSOLID PID MID CORDM IN STRESS ISOP FCTN
    PSOLID   1       1       0
    PSOLID 2 100 6 TWO GRID REDUCED
    """
    type = 'PSOLID'
    def __init__(self,card=None,data=None):
        SolidProperty.__init__(self,card,data)
        if card:
            self.pid    = card.field(1)
            self.mid    = card.field(2)
            self.cordm  = card.field(3,0)
            self.integ  = card.field(4)
            self.stress = card.field(5)
            self.isop   = card.field(6)
            self.fctn   = card.field(7,'SMECH')
        else:
            self.pid    = data[0]
            self.mid    = data[1]
            self.cordm  = data[2]
            self.integ  = data[3]
            self.stress = data[4]
            self.isop   = data[5]
            self.fctn   = data[6]

            if self.fctn=='SMEC':
                self.fctn = 'SMECH'
        ###

    def __repr__(self):
        cordm = self.setBlankIfDefault(self.cordm,0)
        fctn  = self.setBlankIfDefault(self.fctn,'SMECH')
        fields = ['PSOLID',self.pid,self.Mid(),cordm,self.integ,self.stress,self.isop,fctn]
        return self.printCard(fields)

class PCONEAX(Property): #not done
    type = 'PCONEAX'
    def __init__(self,card=None,data=None):
        Property.__init__(self,card,data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.group = card.field(3,'MSCBMLO')
            self.Type = card.field(4)
            self.dim = [] # confusing entry...
        else:
            raise Exception('not supported')
        ###

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        fields = ['PCONEAX',self.pid,self.Mid(),sefl.group,self.Type,]
        return self.printCard(fields)
    
