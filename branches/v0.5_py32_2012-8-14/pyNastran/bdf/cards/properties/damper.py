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
# pylint: disable=C0103,R0902,R0904,R0914



from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Property

class DamperProperty(Property):
    type = 'DamperProperty'
    def __init__(self, card, data):
        Property.__init__(self, card, data)
        pass
    def crossReference(self,model):
        pass

class PVISC(DamperProperty):
    type = 'PVISC'
    def __init__(self,card=None,nPVISC=0,data=None):
        if card:
            self.pid = card.field(1+4*nPVISC)
            self.ce = card.field(2+4*nPVISC)
            self.cr = card.field(3+4*nPVISC,0.)
        else:
            self.pid = data[0]
            self.ce = data[1]
            self.cr = data[2]
        ###
    
    def crossReference(self,model):
        pass

    def rawFields(self):
        fields = ['PVISC',self.pid,self.ce,self.cr]
        return fields
    
    def reprFields(self):
        cr = set_blank_if_default(self.cr,0.)
        fields = ['PVISC',self.pid,self.ce,cr]
        return fields

class PDAMP(DamperProperty):
    type = 'PDAMP'
    def __init__(self,card=None,nPDAMP=0,data=None):
        DamperProperty.__init__(self, card, data)
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
        DamperProperty.__init__(self, card, data)
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
        DamperProperty.__init__(self, card, data)
        if card:
            ## Property ID
            self.pid  = card.field(1)
            ## Identification number of a TABLEDi entry that defines the
            ## damping force per-unit velocity versus frequency relationship
            self.tbid = card.field(2,0)
        else:
            self.pid  = data[0]
            self.tbid = data[1]
        ###

    def crossReference(self,model):
        self.tbid = model.Table(self.tbid)
    
    def Tbid(self):
        if self.tbid==0:
            return None
        elif isinstance(self.tbid,int):
            return self.tbid
        return self.tbid.tid

    def reprFields(self):
        return self.rawFields()

    def rawFields(self):
        fields = ['PDAMPT',self.pid,self.Tbid()]
        return fields

