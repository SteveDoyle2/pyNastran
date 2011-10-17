#import sys
#from numpy import array,cross,dot
#from numpy import array
#from numpy.linalg import norm
from numpy import zeros

# my code
from baseCard import Property

class PBARL(Property): # not done
    type = 'PBARL'
    def __init__(self,card):
        Property.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)
        self.group = card.field(3,'MSCBMLO')
        self.type = card.field(4)
        self.dim = [] # confusing entry...

    def __repr__(self):
        fields = [self.type,self.pid,self.mid,group,type,None,None,None,None,
        ]+self.dim
        return self.printCard(fields)

class PBEAM(Property): # not done, cleanup
    def __init__(self,card):
    type = 'PBEAM'
        Property.__init__(self,card)
        self.pid = card.field(1)
        self.mid = card.field(2)

        self.A   = card.field(3)
        self.I1  = card.field(4)
        self.I2  = card.field(5)
        self.I12 = card.field(6)
        self.J   = card.field(7)
        self.nsm = card.field(8)

        fields = card.fields(9)
        nFields = len(fields)
        # counting continuation cards
        nMajor    = nFields/16
        nLeftover = nFields%16
        if nLeftover:
            nMajor+=1
        
        for nRepeated in range(nMajor-1): # the -1 is for the last group of lines
            nStart = nRepeated*16
            so  = fields[nStart+10] # field 10 is the first possible so
            xxb = fields[nStart+11]
            a   = fields[nStart+12]
            i1  = fields[nStart+13]
            i2  = fields[nStart+14]
            i12 = fields[nStart+15]
            j   = fields[nStart+16]
            nsm = fields[nStart+17]
            c1 = fields[nStart+18]
            c2 = fields[nStart+19]
            d1 = fields[nStart+20]
            d2 = fields[nStart+21]
            e1 = fields[nStart+22]
            e2 = fields[nStart+23]
            f1 = fields[nStart+24]
            f2 = fields[nStart+25]

        # missing repeated lines
        x = nStart+10+16
        self.k1 = card.field(x)
        self.k2 = card.field(x+1)
        self.s1 = card.field(x+2)
        self.s2 = card.field(x+3)
        self.nsia = card.field(x+4)
        self.nsib = card.field(x+5)
        self.cwa = card.field(x+6)
        self.cwb = card.field(x+7)

        self.m1a = card.field(x+8)
        self.m2a = card.field(x+9)
        self.m1b = card.field(x+10)
        self.m2b = card.field(x+11)
        self.n1a = card.field(x+12)
        self.n2a = card.field(x+13)
        self.n1b = card.field(x+14)
        self.n2b = card.field(x+15)

    def __repr__(self):
        raise Exception('not done...')
        fields = [self.type,self.pid,self.mid,] # other
        return self.printCard(fields)
        
#class PBEAML(Property): #not done

class PBEAM3(Property): # not done, cleanup
    def __init__(self,card):
    type = 'PBEAM3'
        Property.__init__(self,card)
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
        fields = [self.type,self.pid,self.mid,] # other
        return self.printCard(fields)

#class PCOMPG(Property): # not done...
#    pass
class PCOMP(Property):
    """
    PCOMP     701512   0.0+0 1.549-2                   0.0+0   0.0+0     SYM
              300704   3.7-2   0.0+0     YES  300704   3.7-2     45.     YES
              300704   3.7-2    -45.     YES  300704   3.7-2     90.     YES
              300705      .5   0.0+0     YES
    """
    type = 'PCOMP'
    def __init__(self,card): # not done, cleanup
        Property.__init__(self,card)
        
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
        fields = [self.type,self.pid,self.z0,self.nsm,self.sb,self.ft,self.TRef,self.ge,self.lam,]
        #print "plies = ",self.plies
        for ply in self.plies:
            (mid,t,theta,sout) = ply
            theta = self.setBlankIfDefault(theta,0.0)
            sout  = self.setBlankIfDefault(sout,'NO')
            fields += [mid,t,theta,sout]
        return self.printCard(fields)

class PELAS(Property):
    type = 'PELAS'
    def __init__(self,card,nPELAS=0):
        Property.__init__(self,card)
        self.pid = card.field(1+5*nPELAS) # 2 PELAS properties can be defined on 1 PELAS (not supported in this program)
        self.k   = card.field(2+5*nPELAS)
        self.ge  = card.field(3+5*nPELAS)
        self.s   = card.field(4+5*nPELAS)

    def __repr__(self):
        fields = [self.type,self.pid,self.mid,self.k,self.ge,self.s]
        return self.printCard(fields)

class PROD(Property):
    type = 'PROD'
    def __init__(self,card):
        Property.__init__(self,card)
        self.pid    = card.field(1)
        self.mid    = card.field(2)
        self.A   = card.field(3)
        self.J   = card.field(4)
        self.c   = card.field(5,0.0)
        self.nsm = card.field(6)

    def __repr__(self):
        c  = self.setBlankIfDefault(self.c,0.0)
        fields = [self.type,self.pid,self.mid,self.A,self.J,c,self.nsm]
        return self.printCard(fields)

class PSHELL(Property):
    """
    PSHELL PID MID1 T MID2 12I/T**3 MID3 TS/T NSM
    Z1 Z2 MID4
    PSHELL   41111   1      1.0000   1               1               0.02081"""
    type = 'PSHELL'
    def __init__(self,card):
        Property.__init__(self,card)
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
        fields = [self.type,self.pid,self.mid,self.t,self.mid2,self.twelveIt3,self.mid3,self.tst,self.nsm,
                  self.z1,self.z2,self.mid4]
        return self.printCard(fields)

class PSOLID(Property):
    """
    PSOLID PID MID CORDM IN STRESS ISOP FCTN
    PSOLID   1       1       0
    PSOLID 2 100 6 TWO GRID REDUCED
    """
    type = 'PSOLID'
    def __init__(self,card):
        Property.__init__(self,card)
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
        fields = [self.type,self.pid,self.mid,cordm,self.integ,self.stress,self.isop,fctn]
        return self.printCard(fields)

class PTUBE(Property):
    type = 'PTUBE'
    def __init__(self,card):
        Property.__init__(self,card)
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
        fields = [self.type,self.pid,self.mid,self.outerDiameter,t,nsm,outerDiameter2]
        return self.printCard(fields)
    
    def massMatrix():
        """@todo not done"""
        m = zeros(6,6))
        m[0,0] = 1.
        return m
