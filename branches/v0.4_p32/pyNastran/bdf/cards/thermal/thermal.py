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
from pyNastran.bdf.cards.baseCard import BaseCard

class ThermalCard(BaseCard):
    def __init__(self,card,data):
        pass
    def crossReference(self,model):
        raise Exception('%s has not defined the crossReference method' %(self.type))
    def isSameCard(self,obj):
        return False

class ThermalBC(ThermalCard):
    def __init__(self,card,data):
        pass

class ThermalElement(ThermalCard):
    pid = 0
    def __init__(self,card,data):
        pass

    def Pid(self):
        if isinstance(self.pid,int):
            return self.pid
        else:
            return self.pid.pid
        ###

class ThermalProperty(ThermalCard):
    def __init__(self,card,data):
        pass

class ThermalLoadDefault(ThermalCard):
    def __init__(self,card,data):
        pass

class ThermalLoad(ThermalCard):
    def __init__(self,card,data):
        pass

#-------------------------------------------------------
# Elements

class CHBDYE(ThermalElement):
    """
    Defines a boundary condition surface element with reference to a heat conduction element.
    """
    type = 'CHBDYE'

    hexMap = {1: [4,3,2,1],  # CHEXA8/CHEXA20
              2: [1,2,6,5],
              3: [2,3,7,6],
              4: [3,4,8,7],
              5: [4,1,5,8],
              6: [5,6,7,8],
             }

    pentMap = {1: [3,2,1],  # CPENTA
               2: [1,2,5,4],
               3: [2,3,6,5],
               4: [3,1,4,6],
               5: [4,5,6],
              }

    tetMap = {1: [1,3,2],  # CTETRA
              2: [1,2,4],
              3: [2,3,4],
              4: [3,1,4],
             }

    sideMaps = {'CHEXA': hexMap,'CPENTA':pentMap,'CTETRA':tetMap,'CTRIA3':[1,2,3],'CQUAD4':[1,2,3,4]}

    def __init__(self,card=None,data=None):
        ThermalElement.__init__(self,card,data)
        ## Surface element ID number for a side of an
        ## element. (0 < Integer < 100,000,000)
        self.eid  = card.field(1)
        
        ## A heat conduction element identification
        self.eid2 = card.field(2)

        ## A consistent element side identification number (1 < Integer < 6)
        self.side = card.field(3)
        assert 0<self.side<7

        ## A VIEW entry identification number for the front face
        self.iViewFront  = card.field(4,0)
        ## A VIEW entry identification number for the back face
        self.iViewBack   = card.field(5,0)

        ## RADM identification number for front face of surface element (Integer > 0)
        self.radMidFront = card.field(6,0)        
        ## RADM identification number for back face of surface element (Integer > 0)
        self.radMidBack  = card.field(7,0)
        
        self.grids = []

    #def crossReference(self,model):
    #    pass
    
    def sideToEIDs(self,eid):
        sideIDs = self.sideMaps[eid.type][self.side]
        ## [1,2,3]
        
        # id-1 is for the 0 based python index
        nodes = [enodes[id-1] for id in range(len(eid.nodes)) if id in sideIDs ]
        return side
        
    def rawFields(self):
        fields = ['CHBDYE',self.eid,self.eid2,self.side,self.iViewFront,self.iViewBack,self.radMidFront,self.radMidBack]
        return fields

    def reprFields(self):
        #eids = self.collapseThruBy(self.eids)  ## @todo is this done
        fields = ['CHBDYE',self.eid,self.eid2,self.side,self.iViewFront,self.iViewBack,self.radMidFront,self.radMidBack]
        return fields

class CHBDYG(ThermalElement):
    """
    Defines a boundary condition surface element without reference to a property entry.
    """
    type = 'CHBDYG'
    def __init__(self,card=None,data=None):
        ThermalElement.__init__(self,card,data)
        
        if card:
            ## Surface element ID
            self.eid  = card.field(1)
            # no field 2

            ## Surface type
            self.Type = card.field(3)
            assert self.Type in ['REV','AREA3','AREA4','AREA6','AREA8']

            ## A VIEW entry identification number for the front face
            self.iViewFront  = card.field(4,0)

            ## A VIEW entry identification number for the back face
            self.iViewBack  = card.field(8,0)

            ## RADM identification number for front face of surface element (Integer > 0)
            self.radMidFront = card.field(6,0)

            ## RADM identification number for back face of surface element (Integer > 0)
            self.radMidBack = card.field(7,0)
            # no field 8

            ## Grid point IDs of grids bounding the surface (Integer > 0)
            self.grids = card.fields(9)
        else:
            self.eid         = data[0]
            self.Type        = data[1]
            self.iViewFront  = data[2]
            self.iViewBack   = data[3]
            self.radMidFront = data[4]
            self.radMidBack  = data[5]
            self.grids       = data[6:14]
        ###

    def crossReference(self,mesh):
        pass
        #self.pid = mesh.Phbdy(self.pid)

    def rawFields(self):
        fields = ['CHBDYG',self.eid,None,self.Type,self.iViewFront,self.iViewBack,self.radMidFront,self.radMidBack,None,
        ]+self.grids
        return fields

    def reprFields(self):
        iViewFront  = self.setBlankIfDefault(self.iViewFront,  0)
        iViewBack   = self.setBlankIfDefault(self.iViewBack,   0)
        radMidFront = self.setBlankIfDefault(self.radMidFront, 0)
        radMidBack  = self.setBlankIfDefault(self.radMidBack,  0)

        fields = ['CHBDYG',self.eid,None,self.Type,iViewFront,iViewBack,radMidFront,radMidBack,None,
        ]+self.grids
        #print "chbdyg fields = %s" %(fields)
        return fields

class CHBDYP(ThermalElement):
    """
    Defines a boundary condition surface element with reference to a PHBDY entry
    """
    type = 'CHBDYP'
    def __init__(self,card=None,data=None):
        ThermalElement.__init__(self,card,data)
        if card:
            ## Surface element ID
            self.eid  = card.field(1)

            ## PHBDY property entry identification numbers. (Integer > 0)
            self.pid  = card.field(2)

            self.Type = card.field(3)
            #print "self.Type = ",self.Type
            #assert self.Type in ['POINT','LINE','ELCYL','FTUBE','TUBE'],'CHBDYP Type=|%s|' (self.Type)

            ## A VIEW entry identification number for the front face.
            self.iViewFront = card.field(4,0)

            ## A VIEW entry identification number for the back face.
            self.iViewBack  = card.field(5,0)

            ## Grid point identification numbers of grids bounding the surface. (Integer > 0)
            self.g1 = card.field(6)        
            ## Grid point identification numbers of grids bounding the surface. (Integer > 0)
            self.g2 = card.field(7)

            ## Orientation grid point. (Integer > 0; Default = 0)
            self.g0 = card.field(8,0)

            ## RADM identification number for front face of surface element. (Integer > 0)
            self.radMidFront = card.field(9,0)

            ## RADM identification number for back face of surface element. (Integer > 0)
            self.radMidBack = card.field(10,0)

            ## Grid point identification number of a midside node if it is used with the line type surface element.
            self.gmid = card.field(11)
            ## Coordinate system for defining orientation vector. (Integer > 0;Default = 0
            self.ce   = card.field(12,0)

            ## Components of the orientation vector in coordinate system CE. The origin of the orientation vector is grid point G1. (Real or blank)
            self.e1   = card.field(13)
            self.e2   = card.field(14)
            self.e3   = card.field(15)
        else:
            raise NotImplementedError()

    def crossReference(self,mesh):
        self.pid = mesh.Phbdy(self.pid)

    def rawFields(self):
        fields = ['CHBDYP',self.eid,self.Pid(),self.Type,self.iViewFront,self.iViewBack,self.g1,self.g2,self.g0,
                  self.radMidFront,self.radMidBack,self.gmid,self.ce,self.e1,self.e2,self.e3]
        return fields

    def reprFields(self):
        iViewFront  = self.setBlankIfDefault(self.iViewFront,  0)
        iViewBack   = self.setBlankIfDefault(self.iViewBack,   0)
        radMidFront = self.setBlankIfDefault(self.radMidFront, 0)
        radMidBack  = self.setBlankIfDefault(self.radMidBack,  0)

        g0 = self.setBlankIfDefault(self.g0, 0)
        ce = self.setBlankIfDefault(self.ce, 0)

        fields = ['CHBDYP',self.eid,self.Pid(),self.Type,iViewFront,iViewBack,self.g1,self.g2,g0,
                  radMidFront,radMidBack,self.gmid,ce,self.e1,self.e2,self.e3]
        return fields

# Elements
#-------------------------------------------------------
# Properties


class PCONV(ThermalProperty):
    """
    Specifies the free convection boundary condition properties of a boundary condition
    surface element used for heat transfer analysis.
    """
    type = 'PCONV'
    def __init__(self,card=None,data=None):
        ## Convection property identification number. (Integer > 0)
        self.pconid = card.field(1)
        ## Material property identification number. (Integer > 0)
        self.mid    = card.field(2)
        ## Type of formula used for free convection. (Integer 0, 1, 10, 11, 20, or 21)
        self.form   = card.field(3,0)
        assert self.form in [0,1,10,11,20,21]

        ## Free convection exponent as implemented within the context of the
        ## particular form that is chosen
        self.expf   = card.field(4,0.0)
        ## Formula type for various configurations of free convection
        self.ftype  = card.field(5,0)
        ## Identification number of a TABLEHT entry that specifies the twovariable
        ## tabular function of the free convection heat transfer coefficient
        self.tid    = card.field(6)
        #7
        #8
        ## Characteristic length
        self.chlen  = card.field(9)
        ## Grid ID of the referenced inlet point
        self.gidin  = card.field(10)
        ## Coordinate system for defining orientation vector. (Integer > 0;Default = 0
        self.ce     = card.field(11,0)
        ## Components of the orientation vector in coordinate system CE. The origin of the orientation vector is grid point G1. (Real or blank)
        self.e1     = card.field(12)
        self.e2     = card.field(13)
        self.e3     = card.field(14)

    #def crossReference(self,model):
    #    pass

    def rawFields(self):
        fields = ['PCONV',self.pconid,self.mid,self.form,self.expf,self.ftype,self.tid,None,None,
                  self.chlen,self.gidin,self.ce,self.e1,self.e2,self.e3]
        return fields

    def reprFields(self):
        form  = self.setBlankIfDefault(self.form,  0)
        expf  = self.setBlankIfDefault(self.expf,0.0)
        ftype = self.setBlankIfDefault(self.ftype, 0)
        ce    = self.setBlankIfDefault(self.ce,    0)
        fields = ['PCONV',self.pconid,self.mid,form,expf,ftype,self.tid,None,None,
                          self.chlen,self.gidin,ce,self.e1,self.e2,self.e3]
        return fields

class PCONVM(ThermalProperty):
    """
    Specifies the free convection boundary condition properties of a boundary condition
    surface element used for heat transfer analysis.
    """
    type = 'PCONVM'
    def __init__(self,card=None,data=None):
        ## Convection property identification number. (Integer > 0)
        self.pconid = card.field(1)
        ## Material property identification number. (Integer > 0)
        self.mid    = card.field(2)
        ## Type of formula used for free convection. (Integer 0, 1, 10, 11, 20, or 21)
        self.form   = card.field(3,0)
        assert self.form in [0,1,10,11,20,21]

        ## Flag for mass flow convection. (Integer = 0 or 1; Default = 0)
        self.flag  = card.field(4,0)
        ## Constant coefficient used for forced convection
        self.coef  = card.field(5)
        ## Reynolds number convection exponent. (Real > 0.0; Default = 0.0)
        self.expr  = card.field(6,0.0)
        ## Prandtl number convection exponent for heat transfer into the working
        ## fluid. (Real > 0.0; Default = 0.0)
        self.exppi = card.field(7,0.0)
        ## Prandtl number convection exponent for heat transfer into the working
        ## fluid. (Real > 0.0; Default = 0.0)
        self.exppo = card.field(8,0.0)

    #def crossReference(self,model):
    #    pass

    def rawFields(self):
        fields = ['PCONVM',self.pconid,self.mid,self.form,self.flag,self.coef,self.expr,self.exppi,self.exppo]
        return fields

    def reprFields(self):
        form  = self.setBlankIfDefault(self.form,   0)
        flag  = self.setBlankIfDefault(self.flag,   0)
        expr  = self.setBlankIfDefault(self.expr, 0.0)
        exppi = self.setBlankIfDefault(self.exppi,0.0)
        exppo = self.setBlankIfDefault(self.exppo,0.0)
        fields = ['PCONVM',self.pconid,self.mid,form,flag,self.coef,expr,exppi,exppo]
        return fields

class PHBDY(ThermalProperty):
    """
    A property entry referenced by CHBDYP entries to give auxiliary geometric
    information for boundary condition surface elements
    """
    type = 'PHBDY'
    def __init__(self,card=None,data=None):
        ## Property identification number. (Unique Integer among all PHBDY entries). (Integer > 0)
        self.pid = card.field(1)

        ## Area factor of the surface used only for CHBDYP element
        ## TYPE = 'POINT', TYPE = 'LINE', TYPE = 'TUBE', or
        ## TYPE = 'ELCYL'. For TYPE = 'TUBE', AF is the constant thickness of
        ## the hollow tube. (Real > 0.0 or blank)
        self.af = card.field(2)
        
        ## Diameters associated with the surface. Used with CHBDYP element TYPE='ELCYL','TUBE','FTUBE'
        self.d1 = card.field(3)
        self.d2 = card.field(4,self.d1)


    #def crossReference(self,model):
    #    pass

    def rawFields(self):
        fields = ['PHBDY',self.pid,self.af,self.d1,self.d2]
        return fields

    def reprFields(self):
        d2  = self.setBlankIfDefault(self.d2,self.d1)
        fields = ['PHBDY',self.pid,self.af,self.d1,d2]
        return fields


# Properties
#-------------------------------------------------------
# Boundary Conditions

class CONV(ThermalBC):
    """
    Specifies a free convection boundary condition for heat transfer analysis through
    connection to a surface element (CHBDYi entry).
    """
    type = 'CONV'
    def __init__(self,card=None,data=None):
        #ThermalBC.__init__(self,card,data)
        ## CHBDYG, CHBDYE, or CHBDYP surface element identification number. (Integer > 0)
        self.eid     = card.field(1)
        
        ## Convection property identification number of a PCONV entry
        self.pconID  = card.field(2)
        ## Point for film convection fluid property temperature
        self.flmnd   = card.field(3,0)
        ## Control point for free convection boundary condition.
        self.cntrlnd = card.field(4,0)

        TA1          = card.field(5)
        nFields = card.nFields()-1 # maybe off by 1...
        
        ## Ambient points used for convection 0's are allowed for TA2 and higher.
        self.ta = card.fields(5,card.nfields,[TA1]*nFields)

    def crossReference(self,model):
        self.eid = model.Element(self.eid)
        assert self.eid.type in ['CHBDYG','CHBDYE','CHBDYP']

    def TA(self,i=None):
        if i is None:
            return self.ta
        return self.ta[i]

    def rawFields(self):
        fields = ['CONV',self.eid,self.pconID,self.flmnd,self.cntrlnd]+self.ta
        return fields

    def reprFields(self):
        flmnd   = self.setBlankIfDefault(self.flmnd,  0)
        cntrlnd = self.setBlankIfDefault(self.cntrlnd,0)
        fields = ['CONV',self.eid,self.pconID,flmnd,cntrlnd]+self.ta
        return fields

class RADM(ThermalBC):
    """
    Defines the radiation properties of a boundary element for heat transfer analysis
    """
    type = 'RADM'
    def __init__(self,card=None,data=None):
        ThermalBC.__init__(self,card,data)
        ## Material identification number
        self.radmid = card.field(1)
        self.absorb = card.field(2)
        self.emissivity = card.fields(3)
        assert self.radmid > 0
        assert 0. <= self.absorb <= 1.0
        for e in self.emissivity:
            assert 0. <= e <= 1.0
        
    #def crossReference(self,model):
    #    pass

    def reprFields(self):
        fields = ['RADM',self.radmid,self.absorb] + self.emissivity
        return fields

class RADBC(ThermalBC):
    """
    Specifies an CHBDYi element face for application of radiation boundary conditions
    """
    type = 'RADBC'
    def __init__(self,card=None,data=None):
        ThermalBC.__init__(self,card,data)
        
        ## NODAMB Ambient point for radiation exchange. (Integer > 0)
        self.nodamb = card.field(1)
        ## Radiation view factor between the face and the ambient point. (Real > 0.0)
        self.famb   = card.field(2)
        ## Control point for thermal flux load. (Integer > 0; Default = 0)
        self.cntrlnd = card.field(3,0)
        
        eids = card.fields(4)
        ## CHBDYi element identification number
        self.eids = self.expandThruBy(eids)
        
    def crossReference(self,model):
        for i,eid in enumerate(self,eids):
            self.eids[i] = model.Element(eid)
        ###

    def Eids(self):
        eids = []
        for eid in self.eids:
            eids.append(self.Eid(eid))
        ###
        return eids

    def Eid(self,eid):
        if isinstance(eid,int):
            return eid
        return eid.eid

    def rawFields(self):
        fields = ['RADBC',self.nodamb,self.famb,self.cntrlnd]+self.Eids()
        return fields

    def reprFields(self):
        cntrlnd = self.setBlankIfDefault(self.cntrlnd,0)
        eids   = self.collapseThruBy(self.Eids())
        fields = ['RADBC',self.nodamb,self.famb,cntrlnd]+eids
        return fields

# Boundary Conditions
#-------------------------------------------------------
# Loads

class QBDY1(ThermalLoad):
    """
    Defines a uniform heat flux into CHBDYj elements.
    """
    type = 'QBDY1'
    def __init__(self,card=None,data=None):
        ThermalLoad.__init__(self,card,data)  # self.sid
        
        if card:
            ## Load set identification number. (Integer > 0)
            self.sid = card.field(1)

            ## Heat flux into element (FLOAT)
            self.qFlux = card.field(2)
            eids       = card.fields(3)
            ## CHBDYj element identification numbers (Integer)
            self.eids = self.expandThru(eids)  # @warning should this use expandThruBy ???
        else:
            self.sid   = data[0]
            self.qFlux = data[1]
            self.eids  = data[2:]
        ###

    def crossReference(self,model):
        self.eid = model.Element(eid)

    def Eid(self):
        if isinstance(self.eid,int):
            return self.eid
        return self.eid.eid

    def nQFluxTerms(self):
        return len(self.qFlux)

    def rawFields(self):
        fields = ['QBDY1',self.qFlux,self.sid]+list(self.eids)+[self.qFlux]
        return fields

    def reprFields(self):
        eids = self.collapseThruBy(self.eids)
        fields = ['QBDY1',self.qFlux,self.sid]+list(eids)+[self.qFlux]
        #print "FIELDS = ",fields
        return fields

class QBDY2(ThermalLoad): # not tested
    """
    Defines a uniform heat flux load for a boundary surface.
    """
    type = 'QBDY2'
    def __init__(self,card=None,data=None):
        ThermalLoad.__init__(self,card,data)

        if card:
            ## Load set identification number. (Integer > 0)
            self.sid   = card.field(1)
            ## Identification number of an CHBDYj element. (Integer > 0)
            self.eid   = card.field(2)
            ## Heat flux at the i-th grid point on the referenced CHBDYj element. (Real or blank)
            self.qFlux = self.removeTrailingNones(card.fields(3))
        else:
            self.sid = data[0]
            self.eid = data[1]
            self.qFlux = data[2]
        ###

    def crossReference(self,model):
        self.eid = model.Element(eid)

    def Eid(self):
        if isinstance(self.eid,int):
            return self.eid
        return self.eid.eid

    def nQFluxTerms(self):
        return len(self.qFlux)

    def rawFields(self):
        fields = ['QBDY2',self.sid,self.Eid(),self.qFlux]
        return fields

    def reprFields(self):
        return self.rawFields()

class QBDY3(ThermalLoad):
    """
    Defines a uniform heat flux load for a boundary surface.
    """
    type = 'QBDY3'
    def __init__(self,card=None,data=None):
        ThermalLoad.__init__(self,card,data)

        if card:
            ## Load set identification number. (Integer > 0)
            self.sid = card.field(1)
            ## Heat flux into element
            self.Q0      = card.field(2)
            ## Control point for thermal flux load. (Integer > 0; Default = 0)
            self.cntrlnd = card.field(3,0)
            eids         = card.fields(4)
            ## CHBDYj element identification numbers
            self.eids = self.expandThruBy(eids)
        else:
            self.sid     = data[0]
            self.Q0      = data[1]
            self.cntrlnd = data[2]
            self.eids    = list(data[3:])
        ###

    def crossReference(self,model):
        for i,eid in enumerate(self,eids):
            self.eids[i] = model.Element(eid)
        ###

    def Eids(self):
        eids = []
        for eid in self.eids:
            eids.append(self.Eid(eid))
        ###
        return eids

    def Eid(self,eid):
        if isinstance(eid,int):
            return eid
        return eid.eid

    def rawFields(self):
        eids = self.Eids()
        eids.sort()
        fields = ['QBDY3',self.sid,self.Q0,self.cntrlnd]+eids
        return fields

    def reprFields(self):
        cntrlnd = self.setBlankIfDefault(self.cntrlnd,0)
        eids = self.collapseThruBy(self.Eids())
        eids.sort()
        fields = ['QBDY3',self.sid,self.Q0,cntrlnd]+eids
        return fields

class QHBDY(ThermalLoad):
    """
    Defines a uniform heat flux into a set of grid points.
    """
    type = 'QHBDY'
    def __init__(self,card=None,data=None):
        ThermalLoad.__init__(self,card,data)
        
        if card:
            ## Load set identification number. (Integer > 0)
            self.sid = card.field(1)

            self.flag = card.field(2)
            assert self.flag in ['POINT','LINE','REV','AREA3','AREA4','AREA6','AREA8']

            ## Magnitude of thermal flux into face. Q0 is positive for heat into the surface. (Real)
            self.Q0   = card.field(3)

            ## Area factor depends on type. (Real > 0.0 or blank)
            self.af    = card.field(4)
            self.grids = card.fields(5)

            ## Grid point identification of connected grid points. (Integer > 0 or blank)
            self.grids = self.expandThruBy(self.grids)
        else:
            self.sid   = data[0]
            self.flag  = data[1]
            self.Q0    = data[2]
            self.af    = data[3]
            self.grids = data[4:]
        ###

    #def crossReference(self,model):
    #    pass

    def rawFields(self):
        fields = ['QHBDY',self.sid,self.flag,self.Q0,self.af]+self.grids
        return fields

    def reprFields(self):
        return self.rawFields()

class TEMP(ThermalLoad):
    """
    Defines temperature at grid points for determination of thermal loading,
    temperature-dependent material properties, or stress recovery.
    """
    type = 'TEMP'
    def __init__(self,card=None,data=None):
        ThermalLoad.__init__(self,card,data)
        
        if card:
            ## Load set identification number. (Integer > 0)
            self.sid = card.field(1)

            fields = card.fields(2)
            nFields = len(fields)
            assert nFields%2==0

            ## dictionary of temperatures where the key is the grid ID (Gi) and the value is the temperature (Ti)
            self.temperatures={}
            for i in range(0,nFields,2):
                self.temperatures[fields[i]] = fields[i+1]
            ###
        else:
            #print "TEMP data = ",data
            self.sid = data[0]
            self.temperatures = {data[1]: data[2]}
        ###

    def add(self,tempObj):
        assert self.sid==tempObj.sid
        for gid,temp in self.tempObj.temperatures.items():
            self.temperatures[gid] = temp
        ###

    def crossReference(self,model):
        pass

    def rawFields(self):
        """Writes the TEMP card"""
        fields = ['TEMP',self.sid]

        nTemps = len(self.temperatures)-1
        for i,(gid,temp) in enumerate(sorted(self.temperatures.items())):
            fields += [gid,temp]
            if i%3==2 and nTemps>i: # start a new TEMP card
                fields += [None,'TEMP',self.sid]
        return fields

    def reprFields(self):
        """Writes the TEMP card"""
        return self.rawFields()

# Loads
#-------------------------------------------------------
# Default Loads

class TEMPD(ThermalLoadDefault):
    """
    Defines a temperature value for all grid points of the structural model that have not
    been given a temperature on a TEMP entry
    """
    type = 'TEMPD'
    def __init__(self,card=None,data=None):
        ThermalLoadDefault.__init__(self,card,data)
        if card:

            fields = card.fields(1)
            nFields = len(fields)
            assert nFields%2==0

            ## dictionary of temperatures where the key is the set ID (SIDi) and the value is the temperature (Ti)
            self.temperatures={}
            for i in range(0,nFields,2):
                self.temperatures[fields[i]] = fields[i+1]
            ###
        else:
            self.temperatures = {data[0]: data[1] }
        ###

    def add(self,tempdObj):
        for sid,tempd in self.tempdObj.temperatures.items():
            self.temperatures[gid] = temp
        ###

    def crossReference(self,model):
        pass

    def reprFields(self):
        """Writes the TEMPD card"""
        fields = ['TEMPD']

        nTemps = len(self.temperatures)-1
        #print "self.temperatures = ",self.temperatures
        #print "nTemps = ",nTemps
        for i,(gid,temp) in enumerate(sorted(self.temperatures.items())):
            fields += [gid,temp]
            if i%4==3 and nTemps>i: # start a new TEMP card
                fields += ['TEMPD']
        return fields

