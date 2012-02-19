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
    def crossReference(self,MODEL):
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
        PointProperty.__init__(self,card,data)
        if card:
            nOffset *= 2
            self.sid   = card.field(1)
            self.Type  = card.field(2)
            self.id    = card.field(3+nOffset)
            self.value = card.field(4+nOffset)
        else:
            self.sid   = data[0]
            #sid=9  propSet=PBEA ID=538976333 value=0.0
            #sid=10 propSet=PDUM ID=538976312 value=2.80259692865e-45
            #sid=10 propSet=ELEM ID=542395973 value=0.0
            self.Type  = data[1]
            self.id    = data[2]
            self.value = data[3]
        ###

    def rawFields(self):
        nodes = self.nodeIDs()
        fields = ['NSM',self.sid,self.Type,self.id,self.value]
        return fields

    def reprFields(self):
        return self.rawFields()

class PMASS(PointProperty):
    def __init__(self,card=None,nOffset=0,data=None):
        PointProperty.__init__(self,card,data)
        nOffset *= 2
        ## Property ID
        self.pid  = card.field(1+nOffset)
        self.mass = card.field(2+nOffset,0.)

    def Mass(self):
        return self.mass

    def rawFields(self):
        fields = ['PMASS',self.pid,self.mass]
        return fields

    def reprFields(self):
        return self.rawFields()

class BushingProperty(Property):
    type = 'BushingProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass
    def crossReference(self,model):
        pass

class PBUSH(BushingProperty):
    type = 'PBUSH'
    def __init__(self,card=None,data=None):
        BushingProperty.__init__(self,card,data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            
            nFields = card.nFields()
            self.vars = []
            iStart = 0
            while iStart<nFields:
                pname = card.field(2)
                if   pname=='K':     self.getK(card,iStart)
                elif pname=='B':   self.getB(card,iStart)
                elif pname=='GE':  self.getGE(card,iStart)
                elif pname=='RCV': self.getRCV(card,iStart)
                else:
                    break
                iStart += 8
            ###
        else:
            self.pid = data[0]
            self.b   = data[1]
            raise NotImplementedError('PBUSH data...')
        ###

    def getK(self,card,iStart):
        ## Flag indicating that the next 1 to 6 fields are stiffness values in the element coordinate system.
        self.k = card.field(iStart)
        ## Nominal stiffness values in directions 1 through 6. See Remarks 2. and 3. (Real; Default = 0.0)
        self.Ki = card.fields(i=iStart+1,j=iStart+6)
        self.vars.append('K')

    def getB(self,card,iStart):
        ## Flag indicating that the next 1 to 6 fields are force-per-velocity damping.
        self.b = card.field(iStart)
        ## Force per unit velocity (Real)
        ## Nominal damping coefficients in direction 1 through 6 in units of force per
        ## unit velocity. See Remarks 2., 3., and 9. (Real; Default = 0.0)
        self.Bi = card.fields(i=iStart+1,j=iStart+6)
        self.vars.append('B')

    def getGE(self,card,iStart):
        ## Flag indicating that the next fields, 1 through 6 are structural damping constants. See Remark 7. (Character)
        self.ge = card.field(iStart)
        ## Nominal structural damping constant in directions 1 through 6. See
        ## Remarks 2. and 3. (Real; Default = 0.0)
        self.GEi = card.fields(i=iStart+1,j=iStart+6)
        self.vars.append('GE')

    def getRCV(self,card,iStart):
        ## Flag indicating that the next 1 to 4 fields are stress or strain coefficients. (Character)
        self.ge = card.field(iStart)
        self.sa = card.field(iStart+1,1.)
        self.sb = card.field(iStart+2,1.)
        self.ea = card.field(iStart+3,1.)
        self.et = card.field(iStart+4,1.)
        self.vars.append('RCV')

    def rawFields(self):
        fields = ['PBUSH',self.pid]
        for var in self.vars:
            if var=='K':
                fields += [None,'K']+self.Ki
            elif var=='B':
                fields += [None,'B']+self.Bi
            elif var=='GE':
                fields += [None,'GE']+self.GEi
            elif var=='RCV':
                fields += [None,'RCV',self.sa,self.st,self.ea,self.et]
            else:
                raise Exception('not supported PBUSH field...')
        return fields

    def reprFields(self):
        return self.rawFields()

class DamperProperty(Property):
    type = 'DamperProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass
    def crossReference(self,model):
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

    def rawFields(self):
        fields = ['PDAMP',self.pid,self.b]
        return fields

    def reprFields(self):
        return self.rawFields()

class PDAMP5(DamperProperty):
    type = 'PDAMP5'
    def __init__(self,card=None,data=None):
        """
        Defines the damping multiplier and references the material properties for damping. CDAMP5 is intended
        for heat transfer analysis only.
        """
        DamperProperty.__init__(self,card,data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            ## Material ID
            self.mid = card.field(2)
            ## Damping multiplier. (Real > 0.0)
            ## B is the mass that multiplies the heat capacity CP on the MAT4 or MAT5 entry.
            self.b   = card.field(3)
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.b   = data[2]
        ###

    def crossReference(self,model):
        self.mid = model.Material(self.mid)
    
    def Mid(self):
        if isinstance(self.mid,int):
            return self.mid
        return self.mid.mid
    
    def reprFields(self):
        return self.rawFields()

    def rawFields(self):
        fields = ['PDAMP5',self.pid,self.Mid(),self.b]
        return fields

class PDAMPT(DamperProperty):
    type = 'PDAMPT'
    def __init__(self,card=None,data=None):
        DamperProperty.__init__(self,card,data)
        if card:
            ## Property ID
            self.pid  = card.field(1)
            ## Identification number of a TABLEDi entry that defines the damping
            ## force per-unit velocity versus frequency relationship
            self.tbid = card.field(2,0)
        else:
            self.pid  = data[0]
            self.tbid = data[1]
        ###

    def crossReference(self,model):
        self.tbid = Table(self.tbid)
    
    def Tbid(self):
        if isinstance(self.tbid,int):
            return self.tbid
        return self.tbid.tid

    def reprFields(self):
        return self.rawFields()

    def rawFields(self):
        fields = ['PDAMPT',self.pid,self.Tbid()]
        return fields

class PGAP(Property):
    type = 'PGAP'
    def __init__(self,card=None,data=None):
        """
        Defines the properties of the gap element (CGAP entry).
        """
        Property.__init__(self,card,data)
        if card:
            ## Property ID
            self.pid   = card.field(1)
            ## initial gap opening
            self.u0    = card.field(2,0.)
            ## preload
            self.f0    = card.field(3,0.)
            ## axial stiffness of closed gap
            self.ka    = card.field(4,1.e8)
            ## axial stiffness of open gap
            self.kb    = card.field(5,1e-14*self.ka)
            ## static friction coeff
            self.mu1   = card.field(7,0.)
            ## transverse stiffness of closed gap
            self.kt    = card.field(6,self.mu1*self.ka)
            ## kinetic friction coeff
            self.mu2   = card.field(8,self.mu1)
            self.tmax  = card.field(9,0.)
            self.mar   = card.field(10,100.)
            self.trmin = card.field(11,0.001)
        else:
            self.pid   = data[0]
            self.u0    = data[1]
            self.f0    = data[2]
            self.ka    = data[3]
            self.kb    = data[4]
            self.kt    = data[5]
            self.mu1   = data[6]
            self.mu2   = data[7]
            self.tmax  = data[8]
            self.mar   = data[9]
            self.trmin = data[10]
        ###

    def rawFields(self):
        fields = ['PGAP',self.pid,self.u0,self.f0,self.ka,self.kb,self.kt,self.mu1,self.mu2,
                         self.tmax,self.mar,self.trmin]
        return fields

    def reprFields(self):
        u0    = self.setBlankIfDefault(self.u0,0.)
        f0    = self.setBlankIfDefault(self.f0,0.)
        ka    = self.setBlankIfDefault(self.ka,1.e8)
        kb    = self.setBlankIfDefault(self.kb,1e-14*self.ka)
        kt    = self.setBlankIfDefault(self.kt,self.mu1*self.ka)
        mu1   = self.setBlankIfDefault(self.mu1,0.)
        mu2   = self.setBlankIfDefault(self.mu2,self.mu1)
        tmax  = self.setBlankIfDefault(self.tmax,0.)
        mar   = self.setBlankIfDefault(self.mar,100.)
        trmin = self.setBlankIfDefault(self.trmin,0.001)
        
        fields = ['PGAP',self.pid,u0,f0,ka,kb,kt,mu1,mu2,
                         tmax,mar,trmin]
        return fields

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
            ## Property ID
            self.pid = card.field(1)
            ## Material ID
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

    def rawFields(self):
        stressStrain = self.setBlankIfDefault(self.str,'GRID')
        fields = ['PLSOLID',self.pid,self.Mid(),stressStrain]
        return fields

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
            ## Property ID
            self.pid    = card.field(1)
            ## Material ID
            self.mid    = card.field(2)
            self.cordm  = card.field(3,0)
            self.integ  = card.field(4)
            #validIntegration = ['THREE','TWO','FULL','BUBBLE',2,3,None,'REDUCED']
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

    def rawFields(self):
        fields = ['PSOLID',self.pid,self.Mid(),self.cordm,self.integ,self.stress,self.isop,self.fctn]
        return fields

    def reprFields(self):
        cordm = self.setBlankIfDefault(self.cordm,0)
        fctn  = self.setBlankIfDefault(self.fctn,'SMECH')
        fields = ['PSOLID',self.pid,self.Mid(),cordm,self.integ,self.stress,self.isop,fctn]
        return fields

class PCONEAX(Property): #not done
    type = 'PCONEAX'
    def __init__(self,card=None,data=None):
        Property.__init__(self,card,data)
        if card:
            ## Property ID
            self.pid   = card.field(1)
            ## Material ID
            self.mid   = card.field(2)
            self.group = card.field(3,'MSCBMLO')
            self.Type  = card.field(4)
            self.dim   = [] # confusing entry...
        else:
            raise NotImplementedError('not supported')
        ###

    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def rawFields(self):
        fields = ['PCONEAX',self.pid,self.Mid(),self.group,self.Type,]
        return fields
    
