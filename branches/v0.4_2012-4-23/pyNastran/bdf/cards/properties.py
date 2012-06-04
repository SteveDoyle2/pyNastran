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
#import sys

# pyNastran
#from baseCard import Property
from pyNastran.bdf.cards.bars.propertiesBars import *
from pyNastran.bdf.cards.plates.propertiesShell import *
from pyNastran.bdf.cards.springs.propertiesSprings import *
from pyNastran.bdf.cards.mass.propertiesMass import *


class BushingProperty(Property):
    type = 'BushingProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass
    def crossReference(self,model):
        pass

class PFAST(Property):
    type = 'PFAST'
    def __init__(self,card=None,data=None):
        Property.__init__(self,card,data)
        ## Property ID
        self.pid   = card.field(1)
        ## diameter of the fastener
        self.d     = card.field(2)
        ## Specifies the element stiffness coordinate system
        self.mcid  = card.field(3,-1)
        self.mflag = card.field(4,0) # 0-absolute 1-relative
        ## stiffness values in directions 1-3
        self.kt1  = card.field(5)
        self.kt2  = card.field(6)
        self.kt3  = card.field(7)
        ## Rotational stiffness values in directions 1-3
        self.kr1  = card.field(8, 0.0)
        self.kr2  = card.field(9, 0.0)
        self.kr3  = card.field(10,0.0)
        ## Lumped mass of fastener
        self.mass = card.field(11,0.0)
        ## Structural damping
        self.ge   = card.field(12,0.0)

    def crossReference(self,model):
        self.mcid = model.Coord(self.mcid)
    
    def Mcid(self):
        if isinstance(self.mcid,int):
            return self.mcid
        return self.mcid.cid

    def Mass(self):
        return self.mass

    def rawFields(self):
        fields = ['PFAST',self.pid,self.d,self.Mcid(),self.mflag,self.kt1,self.kt2,self.kt3,self.kr1,
                          self.kr2,self.kr2,self.mass,self.ge]
    
    def reprFields(self):
        return self.rawFields()

class PBUSH(BushingProperty):
    type = 'PBUSH'
    def __init__(self,card=None,data=None):
        BushingProperty.__init__(self,card,data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            
            nFields = card.nFields()
            self.vars = []
            iStart = 2
            while iStart<nFields:
                pname = card.field(iStart)
                if   pname=='K':   self.getK(card,iStart)
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
        #print self
        #sys.exit()

    def getK(self,card,iStart):
        ## Flag indicating that the next 1 to 6 fields are stiffness values in the element coordinate system.
        #self.k = card.field(iStart)
        ## Nominal stiffness values in directions 1 through 6. See Remarks 2. and 3. (Real; Default = 0.0)
        self.Ki = card.fields(i=iStart+1,j=iStart+6)
        #print "Ki = ",self.Ki
        self.vars.append('K')

    def getB(self,card,iStart):
        ## Flag indicating that the next 1 to 6 fields are force-per-velocity damping.
        #self.b = card.field(iStart)
        ## Force per unit velocity (Real)
        ## Nominal damping coefficients in direction 1 through 6 in units of force per
        ## unit velocity. See Remarks 2., 3., and 9. (Real; Default = 0.0)
        self.Bi = card.fields(i=iStart+1,j=iStart+6)
        self.vars.append('B')

    def getGE(self,card,iStart):
        ## Flag indicating that the next fields, 1 through 6 are structural damping constants. See Remark 7. (Character)
        #self.ge = card.field(iStart)
        ## Nominal structural damping constant in directions 1 through 6. See
        ## Remarks 2. and 3. (Real; Default = 0.0)
        self.GEi = card.fields(i=iStart+1,j=iStart+6)
        self.vars.append('GE')

    def getRCV(self,card,iStart):
        ## Flag indicating that the next 1 to 4 fields are stress or strain coefficients. (Character)
        #self.ge = card.field(iStart)
        self.sa = card.field(iStart+1,1.)
        self.st = card.field(iStart+2,1.)
        self.ea = card.field(iStart+3,1.)
        self.et = card.field(iStart+4,1.)
        self.vars.append('RCV')

    def rawFields(self):
        fields = ['PBUSH',self.pid]
        for var in self.vars:
            if var=='K':
                fields += ['K']+self.Ki
            elif var=='B':
                fields += ['B']+self.Bi
            elif var=='GE':
                fields += ['GE']+self.GEi
            elif var=='RCV':
                fields += ['RCV',self.sa,self.st,self.ea,self.et]
            else:
                raise Exception('not supported PBUSH field...')
            nSpaces = 8-(len(fields)-1)%8
            
            #print "nSpaces = ",nSpaces
            if nSpaces<8:
                fields += [None]*(nSpaces+1)
            ###
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
            #(pid,u0,f0,ka,kb,kt,mu1,mu2,tmax,mar,trmin) = out
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
    
