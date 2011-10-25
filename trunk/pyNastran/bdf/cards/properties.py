#import sys
#from numpy import array,cross,dot
#from numpy import array
#from numpy.linalg import norm
from numpy import zeros

# my code
from baseCard import Property

class LineProperty(Property):
    type = 'LineProperty'
    def __init__(self,card):
        Property.__init__(self,card)
        pass

class ShellProperty(Property):
    type = 'ShellProperty'
    def __init__(self,card):
        Property.__init__(self,card)
        pass

class SolidProperty(Property):
    type = 'SolidProperty'
    def __init__(self,card):
        Property.__init__(self,card)
        pass

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
        print "line3 = ",line3
        
        line1 = ['PBAR',self.pid,self.mid,self.A,I1,I2,J,nsm,None]

        if line3==[None,None,None]:
            line2 = [C1,C2,D1,D2,E1,E2,F1,F2]
        else:
            line2 = [C1,C2,D1,D2,E1,E2,F1,self.F2]
        fields = line1+line2+line3

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

    def __repr__(self):
        fields = ['PCONEAX',self.pid,self.mid]
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

    def __repr__(self):
        fields = ['PBARL',self.pid,self.mid,group,type,None,None,None,None,
        ]+self.dim+[self.nsm]
        return self.printCard(fields)

class PBEAM(LineProperty):
    type = 'PBEAM'
    def __init__(self,card):
        """@todo cleanup entries"""
        LineProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)

        self.A   = card.field(3)
        self.I1  = card.field(4)
        self.I2  = card.field(5)
        self.I12 = card.field(6)
        self.J   = card.field(7)
        self.NSM = card.field(8)
        self.C1  = card.field(9)
        self.C2  = card.field(10)
        self.D1  = card.field(11)
        self.D2  = card.field(12)
        self.E1  = card.field(13)
        self.E2  = card.field(14)
        self.F1  = card.field(15)
        self.F2  = card.field(16)

        self.so  = []
        self.xxb = []
        self.a   = []
        self.i1  = []
        self.i2  = []
        self.i12 = []
        self.j   = []
        self.nsm = []
        self.c1 = []
        self.c2 = []
        self.d1 = []
        self.d2 = []
        self.e1 = []
        self.e2 = []
        self.f1 = []
        self.f2 = []
        
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
            self.a.append(  propFields[2])
            self.i1.append( propFields[3])
            self.i2.append( propFields[4])
            self.i12.append(propFields[5])
            self.j.append(  propFields[6])
            self.nsm.append(propFields[7])
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

    def __repr__(self):
        fields = ['PBEAM',self.pid,self.mid,self.A, self.I1,self.I2,self.I12,self.J, self.NSM,
                          self.C1, self.C2, self.D1,self.D2,self.E1,self.E2, self.F1,self.F2]
        #print "fieldsA = ",fields
        
        #print len(self.so)
        for (so,xxb,a,i1,i2,i12,j,nsm,c1,c2,d1,d2,e1,e2,f1,f2) in zip(
            self.so,self.xxb,self.a,self.i1,self.i2,self.i12,self.j,self.nsm,
            self.c1,self.c2,self.d1,self.d2,self.e1,self.e2,self.f1,self.f2):
            fields += [so,xxb,a,i1,i2,i12,j,nsm,c1,c2,d1,d2,e1,e2,f1,f2]
            #print "asdf = ",asdf
        fields += [self.k1,self.k2,self.s1,self.s2,self.nsia,self.nsib,self.cwa,self.cwb,
                   self.m1a,self.m2a,self.m1b,self.m2b,self.n1a,self.n2a,self.n1b,self.n2b]
        #print fields
        #print "asdf = ",asdf
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

    def __repr__(self):
        raise Exception('not done...')
        fields = ['PBEAM3',self.pid,self.mid,] # other
        return self.printCard(fields)

#class PCOMPG(ShellProperty): # not done...
#    pass
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
        self.lam  = card.field(8)  # symmetric flag???
        
        nPlyFields = card.nFields()-8 # -8 for the first 8 fields (1st line)
        #plyCards = card.fields(9)
        
        # counting plies
        nMajor    = nPlyFields/4
        nLeftover = nPlyFields%4
        if nLeftover:
            nMajor+=1
        self.nplies = nMajor

        iPly = 0
        self.plies = []
        for i in range(9,self.nplies*4,4):
            defaults = [None,None,0.0,'NO']
            (mid,t,theta,sout) = card.fields(i,i+4,defaults)
            self.plies.append([mid,t,theta,sout])
            iPly +=1
        #print "nplies = ",self.nplies
        #print str(self)

    def __repr__(self):
        fields = ['PCOMP',self.pid,self.z0,self.nsm,self.sb,self.ft,self.TRef,self.ge,self.lam,]
        #print "plies = ",self.plies
        for ply in self.plies:
            (mid,t,theta,sout) = ply
            theta = self.setBlankIfDefault(theta,0.0)
            sout  = self.setBlankIfDefault(sout,'NO')
            fields += [mid,t,theta,sout]
        return self.printCard(fields)

class PELAS(LineProperty):
    type = 'PELAS'
    def __init__(self,card,nPELAS=0):
        LineProperty.__init__(self,card)
        self.pid = card.field(1+5*nPELAS) # 2 PELAS properties can be defined on 1 PELAS card
        self.k   = card.field(2+5*nPELAS) # these are split into 2 separate cards
        self.ge  = card.field(3+5*nPELAS)
        self.s   = card.field(4+5*nPELAS)

    def __repr__(self):
        fields = ['PELAS',self.pid,self.k,self.ge,self.s]
        return self.printCard(fields)

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
        

    def __repr__(self):
        stressStrain = self.setDefaultIfNone(self.str,'GRID')
        fields = ['PLSOLID',self.pid,self.mid,stressStrain]
        return self.printCard(fields)

class PROD(LineProperty):
    type = 'PROD'
    def __init__(self,card):
        LineProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)
        self.A   = card.field(3)
        self.J   = card.field(4)
        self.c   = card.field(5,0.0)
        self.nsm = card.field(6)

    def __repr__(self):
        c  = self.setBlankIfDefault(self.c,0.0)
        fields = ['PROD',self.pid,self.mid,self.A,self.J,c,self.nsm]
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
        self.mid2 = card.field(4)
        self.twelveIt3 = card.field(5,1.0)  # poor name
        self.mid3  = card.field(6)
        self.tst   = card.field(7,0.833333)
        self.nsm   = card.field(8,0.0)
        self.z1    = card.field(9)
        self.z2    = card.field(10)
        self.mid4  = card.field(11)

    def __repr__(self):
        fields = ['PSHELL',self.pid,self.mid,self.t,self.mid2,self.twelveIt3,self.mid3,self.tst,self.nsm,
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
        fields = ['PSOLID',self.pid,self.mid,cordm,self.integ,self.stress,self.isop,fctn]
        return self.printCard(fields)

class PTUBE(LineProperty):
    type = 'PTUBE'
    def __init__(self,card):
        LineProperty.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)
        self.outerDiameter = card.field(3)
        self.t   = card.field(4,self.outerDiameter/2.)
        self.nsm = card.field(5,0.0)
        self.outerDiameter2 = card.field(6,self.outerDiameter)

    def __repr__(self):
        t   = self.setBlankIfDefault(self.t,self.outerDiameter/2.)
        nsm = self.setBlankIfDefault(self.nsm,0.0)
        outerDiameter2 = self.setBlankIfDefault(self.outerDiameter2,self.outerDiameter)
        fields = ['PTUBE',self.pid,self.mid,self.outerDiameter,t,nsm,outerDiameter2]
        return self.printCard(fields)
    
    def massMatrix():
        """@todo not done"""
        m = zeros(6,6)
        m[0,0] = 1.
        return m
