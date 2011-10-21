from baseCard import BaseCard

class FLFACT(BaseCard):
    """
    FLFACT SID F1 F2 F3 F4 F5 F6 F7
    F8 F9 -etc.-
    
    FLFACT 97 .3 .7 3.5
    
    FLFACT SID F1 THRU FNF NF FMID       # delta quantity approach
    FLFACT 201 .200 THRU .100 11 .133333
    """
    type = 'FLFACT'
    def __init__(self,card):
        #Material.__init__(self,card)
        self.sid     = card.field(1)
        self.factors = card.fields(2)
        if self.factors[1]=='THRU':
            raise Exception('embedded THRUs not supported yet on FLFACT card\n')
            #(a,thru,b,n,dn) = factors
            #for i in range(

    def __repr__(self):
        fields = ['FLFACT',self.sid]+self.factors
        return self.printCard(fields)

class GUST(BaseCard):
    """
    Defines a stationary vertical gust for use in aeroelastic response analysis.
    GUST SID DLOAD WG  X0   V
    GUST 133 61    1.0 0.   1.+4
    """
    type = 'GUST'
    def __init__(self,card):
        self.sid   = card.field(1)
        self.dload = card.field(2)
        self.wg    = card.field(3)
        self.x0    = card.field(4)
        self.V     = card.field(5)
        #angle = self.wg*self.t*(t-(x-self.x0)/self.V) # T is the tabular function

    def __repr__(self):
        fields = ['GUST',self.sid,self.dload,self.wg,self.x0,self.V]
        return self.printCard(fields)


class Aero(BaseCard):
    """Base class for AERO and AEROS cards."""
    def __init__(self,card):
        pass

    def isSymmetricalXY(self):
        if self.symXY==1:
            return True
        return False

    def isSymmetricalXZ(self):
        if self.symXZ==1:
            return True
        return False

    def setGroundEffect(self):
        self.symXY = -1

    def isAntiSymmetricalXY(self):
        if self.symXY==-1:
            return True
        return False

    def isAntiSymmetricalXZ(self):
        if self.symXY==-1:
            return True
        return False


class AERO(BaseCard):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.
    AERO ACSID VELOCITY REFC RHOREF SYMXZ SYMXY
    AERO 3     1.3+4    100.  1.-5  1     -1
    """
    type = 'AERO'
    def __init__(self,card):
        Aero.__init__(self,card)
        self.acsid    = card.field(1)
        self.velocity = card.field(2)
        self.cRef     = card.field(3)
        self.rhoRef   = card.field(4)
        self.symXZ    = card.field(5,0)
        self.symXY    = card.field(6,0)
        #angle = self.wg*self.t*(t-(x-self.x0)/self.V) # T is the tabular function

    def __repr__(self):
        symXZ = self.setBlankIfDefault(self.symXZ,0)
        symXY = self.setBlankIfDefault(self.symXY,0)
        fields = ['AERO',self.acsid,self.velocity,self.cRef,self.rhoRef,symXZ,symXY]
        return self.printCard(fields)

class AEROS(BaseCard):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.
    AEROS ACSID RCSID REFC REFB REFS SYMXZ SYMXY
    AEROS 10   20     10.  100. 1000. 1
    """
    type = 'AEROS'
    def __init__(self,card):
        Aero.__init__(self,card)
        self.acsid  = card.field(1)
        self.rcsid  = card.field(2)
        self.cRef   = card.field(3)
        self.bRef   = card.field(4)
        self.Sref   = card.field(4)
        self.symXZ  = card.field(5,0)
        self.symXY  = card.field(6,0)

    def __repr__(self):
        symXZ = self.setBlankIfDefault(self.symXZ,0)
        symXY = self.setBlankIfDefault(self.symXY,0)
        fields = ['AEROS',self.acsid,self.velocity,self.cRef,self.rhoRef,symXZ,symXY]
        return self.printCard(fields)

class GRAV(BaseCard):
    """
    Defines acceleration vectors for gravity or other acceleration loading
    GRAV SID CID A     N1  N2 N3    MB
    GRAV 1   3   32.2 0.0 0.0 -1.0
    """
    type = 'GRAV'
    def __init__(self,card):
        Aero.__init__(self,card)
        self.sid  = card.field(1)
        self.cid  = card.field(2,0)
        self.a    = card.field(3)
        self.N    = card.fields(4,7,[0.,0.,0.]) # measured in cid
        assert max(self.N) != min(self.N)
        self.mb     = card.field(7,0)

    def __repr__(self):
        N = []
        for n in self.N:
            N.append(self.setBlankIfDefault(n,0.0))
        fields = ['GRAV',self.sid,self.cid,self.a]+N+[self.mb]
        return self.printCard(fields)
