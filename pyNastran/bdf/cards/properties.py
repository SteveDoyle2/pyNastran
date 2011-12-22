#import sys

# my code
#from baseCard import Property
from    bars.propertiesBars    import *
from  plates.propertiesShell   import *
from springs.propertiesSprings import *

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
        if card:
            nOffset *= 2
            self.sid   = card.field(1)
            self.Type  = card.field(2)
            self.id    = card.field(3+nOffset)
            self.value = card.field(4+nOffset)
        else:
            self.sid   = data[0]
            #sid=9 propSet=PBEA ID=538976333 value=0.0
            #sid=10 propSet=PDUM ID=538976312 value=2.80259692865e-45
            #sid=10 propSet=ELEM ID=542395973 value=0.0
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
        PointProperty.__init__(self,card,data)
        
        nOffset *= 2
        self.pid  = card.field(1+nOffset)
        self.mass = card.field(2+nOffset,0.)

    def Mass(self):
        return self.mass

    def __repr__(self):
        fields = ['PMASS',self.pid,self.Mass]

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

class SolidProperty(Property):
    type = 'SolidProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass

    def Rho(self):
        self.mid.rho

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
        assert self.str in ['GRID','GAUS'],'STR=|%s| doesnt have a valid stress/strain output value set\n' %(self.str)

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def __repr__(self):
        stressStrain = self.setBlankIfDefault(self.str,'GRID')
        fields = ['PLSOLID',self.pid,self.Mid(),stressStrain]
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
    
