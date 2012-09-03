# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import (BaseCard, expand_thru_by,
                                          collapse_thru_by)


class ThermalCard(BaseCard):
    def __init__(self, card, data):
        pass

    def cross_reference(self, model):
        raise NotImplementedError('%s has not defined the cross_reference'
                                  'method' % (self.type))

    def isSameCard(self, obj, debug=False):
        return False


class ThermalBC(ThermalCard):
    def __init__(self, card, data):
        pass


class ThermalElement(ThermalCard):
    pid = 0

    def __init__(self, card, data):
        pass

    def nodeIDs(self):
        return []

    def Pid(self):
        if isinstance(self.pid, int):
            return self.pid
        else:
            return self.pid.pid


class ThermalProperty(ThermalCard):
    def __init__(self, card, data):
        pass

#-------------------------------------------------------
# Elements


class CHBDYE(ThermalElement):
    """
    Defines a boundary condition surface element with reference to a heat
    conduction element.
    """
    type = 'CHBDYE'

    hexMap = {1: [4, 3, 2, 1],  # CHEXA8/CHEXA20
              2: [1, 2, 6, 5],
              3: [2, 3, 7, 6],
              4: [3, 4, 8, 7],
              5: [4, 1, 5, 8],
              6: [5, 6, 7, 8],
              }

    pentMap = {1: [3, 2, 1],  # CPENTA
               2: [1, 2, 5, 4],
               3: [2, 3, 6, 5],
               4: [3, 1, 4, 6],
               5: [4, 5, 6],
               }

    tetMap = {1: [1, 3, 2],  # CTETRA
              2: [1, 2, 4],
              3: [2, 3, 4],
              4: [3, 1, 4],
              }

    sideMaps = {'CHEXA': hexMap, 'CPENTA': pentMap, 'CTETRA': tetMap,
                'CTRIA3': [1, 2, 3], 'CQUAD4': [1, 2, 3, 4]}

    def __init__(self, card=None, data=None):
        ThermalElement.__init__(self, card, data)
        ## Surface element ID number for a side of an
        ## element. (0 < Integer < 100,000,000)
        self.eid = card.field(1)

        ## A heat conduction element identification
        self.eid2 = card.field(2)

        ## A consistent element side identification number (1 < Integer < 6)
        self.side = card.field(3)
        assert 0 < self.side < 7

        ## A VIEW entry identification number for the front face
        self.iViewFront = card.field(4, 0)
        ## A VIEW entry identification number for the back face
        self.iViewBack = card.field(5, 0)

        ## RADM identification number for front face of surface element
        ## (Integer > 0)
        self.radMidFront = card.field(6, 0)
        ## RADM identification number for back face of surface element (Integer > 0)
        self.radMidBack = card.field(7, 0)

        self.grids = []

    #def cross_reference(self,model):
    #    pass

    def sideToEIDs(self, eid):
        sideIDs = self.sideMaps[eid.type][self.side]
        ## [1,2,3]

        # id-1 is for the 0 based python index
        nodes = [enodes[id - 1] for id in xrange(len(eid.nodes)) if id in sideIDs]
        return side

    def rawFields(self):
        fields = ['CHBDYE', self.eid, self.eid2, self.side, self.iViewFront,
                  self.iViewBack, self.radMidFront, self.radMidBack]
        return fields

    def reprFields(self):
        #eids = collapse_thru_by(self.eids)  ## @todo is this done
        fields = ['CHBDYE', self.eid, self.eid2, self.side, self.iViewFront,
                  self.iViewBack, self.radMidFront, self.radMidBack]
        return fields


class CHBDYG(ThermalElement):
    """
    Defines a boundary condition surface element without reference to a
    property entry.
    """
    type = 'CHBDYG'

    def __init__(self, card=None, data=None):
        ThermalElement.__init__(self, card, data)

        if card:
            ## Surface element ID
            self.eid = card.field(1)
            # no field 2

            ## Surface type
            self.Type = card.field(3)
            assert self.Type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8']

            ## A VIEW entry identification number for the front face
            self.iViewFront = card.field(4, 0)

            ## A VIEW entry identification number for the back face
            self.iViewBack = card.field(8, 0)

            ## RADM identification number for front face of surface element
            ## (Integer > 0)
            self.radMidFront = card.field(6, 0)

            ## RADM identification number for back face of surface element
            ## (Integer > 0)
            self.radMidBack = card.field(7, 0)
            # no field 8

            ## Grid point IDs of grids bounding the surface (Integer > 0)
            self.grids = card.fields(9)
        else:
            self.eid = data[0]
            self.Type = data[1]
            self.iViewFront = data[2]
            self.iViewBack = data[3]
            self.radMidFront = data[4]
            self.radMidBack = data[5]
            self.grids = data[6:14]

    def cross_reference(self, mesh):
        pass
        #self.pid = mesh.Phbdy(self.pid)

    def rawFields(self):
        fields = ['CHBDYG', self.eid, None, self.Type, self.iViewFront,
                  self.iViewBack, self.radMidFront, self.radMidBack, None,
                  ] + self.grids
        return fields

    def reprFields(self):
        iViewFront = set_blank_if_default(self.iViewFront, 0)
        iViewBack = set_blank_if_default(self.iViewBack, 0)
        radMidFront = set_blank_if_default(self.radMidFront, 0)
        radMidBack = set_blank_if_default(self.radMidBack, 0)

        fields = ['CHBDYG', self.eid, None, self.Type, iViewFront, iViewBack,
                  radMidFront, radMidBack, None, ] + self.grids
        return fields


class CHBDYP(ThermalElement):
    """
    Defines a boundary condition surface element with reference to a PHBDY entry
    """
    type = 'CHBDYP'

    def __init__(self, card=None, data=None):
        ThermalElement.__init__(self, card, data)
        if card:
            ## Surface element ID
            self.eid = card.field(1)

            ## PHBDY property entry identification numbers. (Integer > 0)
            self.pid = card.field(2)

            self.Type = card.field(3)
            #print "self.Type = ",self.Type
            #assert self.Type in ['POINT','LINE','ELCYL','FTUBE','TUBE'],'CHBDYP Type=|%s|' (self.Type)

            ## A VIEW entry identification number for the front face.
            self.iViewFront = card.field(4, 0)

            ## A VIEW entry identification number for the back face.
            self.iViewBack = card.field(5, 0)

            ## Grid point identification numbers of grids bounding the surface. (Integer > 0)
            self.g1 = card.field(6)
            ## Grid point identification numbers of grids bounding the surface. (Integer > 0)
            self.g2 = card.field(7)

            ## Orientation grid point. (Integer > 0; Default = 0)
            self.g0 = card.field(8, 0)

            ## RADM identification number for front face of surface element. (Integer > 0)
            self.radMidFront = card.field(9, 0)

            ## RADM identification number for back face of surface element. (Integer > 0)
            self.radMidBack = card.field(10, 0)

            ## Grid point identification number of a midside node if it is used with the line type surface element.
            self.gmid = card.field(11)
            ## Coordinate system for defining orientation vector. (Integer > 0;Default = 0
            self.ce = card.field(12, 0)

            ## Components of the orientation vector in coordinate system CE.
            ## The origin of the orientation vector is grid point G1.
            ## (Real or blank)
            self.e1 = card.field(13)
            self.e2 = card.field(14)
            self.e3 = card.field(15)
        else:
            raise NotImplementedError()

    def cross_reference(self, mesh):
        self.pid = mesh.Phbdy(self.pid)

    def rawFields(self):
        fields = ['CHBDYP', self.eid, self.Pid(), self.Type, self.iViewFront,
                  self.iViewBack, self.g1, self.g2, self.g0, self.radMidFront,
                  self.radMidBack, self.gmid, self.ce, self.e1, self.e2,
                  self.e3]
        return fields

    def reprFields(self):
        iViewFront = set_blank_if_default(self.iViewFront, 0)
        iViewBack = set_blank_if_default(self.iViewBack, 0)
        radMidFront = set_blank_if_default(self.radMidFront, 0)
        radMidBack = set_blank_if_default(self.radMidBack, 0)

        g0 = set_blank_if_default(self.g0, 0)
        ce = set_blank_if_default(self.ce, 0)

        fields = ['CHBDYP', self.eid, self.Pid(), self.Type, iViewFront,
                  iViewBack, self.g1, self.g2, g0, radMidFront, radMidBack,
                  self.gmid, ce, self.e1, self.e2, self.e3]
        return fields

# Elements
#-------------------------------------------------------
# Properties


class PCONV(ThermalProperty):
    """
    Specifies the free convection boundary condition properties of a boundary
    condition surface element used for heat transfer analysis.
    """
    type = 'PCONV'

    def __init__(self, card=None, data=None):
        ## Convection property identification number. (Integer > 0)
        self.pconid = card.field(1)
        ## Material property identification number. (Integer > 0)
        self.mid = card.field(2)
        ## Type of formula used for free convection. (Integer 0, 1, 10, 11, 20, or 21)
        self.form = card.field(3, 0)
        assert self.form in [0, 1, 10, 11, 20, 21]

        ## Free convection exponent as implemented within the context of the
        ## particular form that is chosen
        self.expf = card.field(4, 0.0)
        ## Formula type for various configurations of free convection
        self.ftype = card.field(5, 0)
        ## Identification number of a TABLEHT entry that specifies the twovariable
        ## tabular function of the free convection heat transfer coefficient
        self.tid = card.field(6)

        ## Characteristic length
        self.chlen = card.field(9)
        ## Grid ID of the referenced inlet point
        self.gidin = card.field(10)
        ## Coordinate system for defining orientation vector.
        ## (Integer > 0;Default = 0
        self.ce = card.field(11, 0)
        ## Components of the orientation vector in coordinate system CE. The
        ## origin of the orientation vector is grid point G1. (Real or blank)
        self.e1 = card.field(12)
        self.e2 = card.field(13)
        self.e3 = card.field(14)

    #def cross_reference(self,model):
    #    pass

    def rawFields(self):
        fields = ['PCONV', self.pconid, self.mid, self.form, self.expf,
                  self.ftype, self.tid, None, None, self.chlen, self.gidin,
                  self.ce, self.e1, self.e2, self.e3]
        return fields

    def reprFields(self):
        form = set_blank_if_default(self.form, 0)
        expf = set_blank_if_default(self.expf, 0.0)
        ftype = set_blank_if_default(self.ftype, 0)
        ce = set_blank_if_default(self.ce, 0)
        fields = ['PCONV', self.pconid, self.mid, form, expf, ftype, self.tid,
                  None, None, self.chlen, self.gidin, ce, self.e1, self.e2,
                  self.e3]
        return fields


class PCONVM(ThermalProperty):
    """
    Specifies the free convection boundary condition properties of a boundary
    condition surface element used for heat transfer analysis.
    """
    type = 'PCONVM'

    def __init__(self, card=None, data=None):
        ## Convection property identification number. (Integer > 0)
        self.pconid = card.field(1)
        ## Material property identification number. (Integer > 0)
        self.mid = card.field(2)
        ## Type of formula used for free convection.
        ## (Integer 0, 1, 10, 11, 20, or 21)
        self.form = card.field(3, 0)
        assert self.form in [0, 1, 10, 11, 20, 21]

        ## Flag for mass flow convection. (Integer = 0 or 1; Default = 0)
        self.flag = card.field(4, 0)
        ## Constant coefficient used for forced convection
        self.coef = card.field(5)
        ## Reynolds number convection exponent. (Real > 0.0; Default = 0.0)
        self.expr = card.field(6, 0.0)
        ## Prandtl number convection exponent for heat transfer into the working
        ## fluid. (Real > 0.0; Default = 0.0)
        self.exppi = card.field(7, 0.0)
        ## Prandtl number convection exponent for heat transfer into the working
        ## fluid. (Real > 0.0; Default = 0.0)
        self.exppo = card.field(8, 0.0)

    #def cross_reference(self,model):
    #    pass

    def rawFields(self):
        fields = ['PCONVM', self.pconid, self.mid, self.form,
                  self.flag, self.coef, self.expr, self.exppi, self.exppo]
        return fields

    def reprFields(self):
        form = set_blank_if_default(self.form, 0)
        flag = set_blank_if_default(self.flag, 0)
        expr = set_blank_if_default(self.expr, 0.0)
        exppi = set_blank_if_default(self.exppi, 0.0)
        exppo = set_blank_if_default(self.exppo, 0.0)
        fields = ['PCONVM', self.pconid, self.mid, form, flag,
                  self.coef, expr, exppi, exppo]
        return fields


class PHBDY(ThermalProperty):
    """
    A property entry referenced by CHBDYP entries to give auxiliary geometric
    information for boundary condition surface elements
    """
    type = 'PHBDY'

    def __init__(self, card=None, data=None):
        ## Property identification number. (Unique Integer among all PHBDY
        ## entries). (Integer > 0)
        self.pid = card.field(1)

        ## Area factor of the surface used only for CHBDYP element
        ## TYPE = 'POINT', TYPE = 'LINE', TYPE = 'TUBE', or
        ## TYPE = 'ELCYL'. For TYPE = 'TUBE', AF is the constant thickness of
        ## the hollow tube. (Real > 0.0 or blank)
        self.af = card.field(2)

        ## Diameters associated with the surface. Used with CHBDYP element
        ## TYPE='ELCYL','TUBE','FTUBE'
        self.d1 = card.field(3)
        self.d2 = card.field(4, self.d1)

    #def cross_reference(self,model):
    #    pass
    def rawFields(self):
        fields = ['PHBDY', self.pid, self.af, self.d1, self.d2]
        return fields

    def reprFields(self):
        d2 = set_blank_if_default(self.d2, self.d1)
        fields = ['PHBDY', self.pid, self.af, self.d1, d2]
        return fields


# Properties
#-------------------------------------------------------
# Boundary Conditions

class CONV(ThermalBC):
    """
    Specifies a free convection boundary condition for heat transfer analysis
    through connection to a surface element (CHBDYi entry).
    """
    type = 'CONV'

    def __init__(self, card=None, data=None):
        #ThermalBC.__init__(self, card, data)
        ## CHBDYG, CHBDYE, or CHBDYP surface element identification number.
        ## (Integer > 0)
        self.eid = card.field(1)

        ## Convection property identification number of a PCONV entry
        self.pconID = card.field(2)
        ## Point for film convection fluid property temperature
        self.flmnd = card.field(3, 0)
        ## Control point for free convection boundary condition.
        self.cntrlnd = card.field(4, 0)

        TA1 = card.field(5)
        nFields = card.nFields() - 1  # maybe off by 1...

        ## Ambient points used for convection 0's are allowed for TA2 and higher.
        self.ta = card.fields(5, card.nfields, [TA1] * nFields)

    def cross_reference(self, model):
        self.eid = model.Element(self.eid)
        assert self.eid.type in ['CHBDYG', 'CHBDYE', 'CHBDYP']

    def TA(self, i=None):
        if i is None:
            return self.ta
        return self.ta[i]

    def rawFields(self):
        fields = ['CONV', self.eid, self.pconID, self.flmnd,
                  self.cntrlnd] + self.ta
        return fields

    def reprFields(self):
        flmnd = set_blank_if_default(self.flmnd, 0)
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)
        fields = ['CONV', self.eid, self.pconID, flmnd, cntrlnd] + self.ta
        return fields

class CONVM(ThermalBC):
    """
    Specifies a forced convection boundary condition for heat transfer analysis through connection to a
    surface element (CHBDYi entry).
    """
    type = 'CONV'

    def __init__(self, card=None, data=None):
        self.eid = card.field(1)
        self.pconvmID = card.field(2)
        self.filmNode = card.field(3, 0)
        self.cntmdot = card.field(4)
        self.ta1 = card.field(5)
        self.ta2 = card.field(6,self.ta1)
        self.mdot = card.field(7, 1.0)

    def cross_reference(self, model):
        self.eid = model.CYBDY(self.eid)
        self.pconvmID = model.PCONV(self.pconvmID)
        self.filmNode = model.Grid(self.filmNode)

    def film_node(self):
        if isinstance(self.filmNode, int):
            return self.filmNode
        return self.filmNode.nid

    def rawFields(self):
        fields = ['CONVM', self.eid, self.pconvmID, self.filmNode, self.cntmdot,
                  self.ta1, self.ta2, self.mdot]
        return fields

    def reprFields(self):
        filmNode = set_blank_if_default(self.filmNode, 0)
        ta2 = set_blank_if_default(self.ta2, self.ta1)
        mdot = set_blank_if_default(self.mdot, 1.0)
        fields = ['CONVM', self.eid, self.pconvmID, filmNode, self.cntmdot,
                  self.ta1, ta2, mdot]
        return fields


class RADM(ThermalBC):
    """
    Defines the radiation properties of a boundary element for heat transfer
    analysis
    """
    type = 'RADM'

    def __init__(self, card=None, data=None):
        ThermalBC.__init__(self, card, data)
        ## Material identification number
        self.radmid = card.field(1)
        self.absorb = card.field(2)
        self.emissivity = card.fields(3)
        assert self.radmid > 0
        assert 0. <= self.absorb <= 1.0
        for e in self.emissivity:
            assert 0. <= e <= 1.0

    #def cross_reference(self,model):
    #    pass

    def reprFields(self):
        fields = ['RADM', self.radmid, self.absorb] + self.emissivity
        return fields


class RADBC(ThermalBC):
    """
    Specifies an CHBDYi element face for application of radiation boundary
    conditions
    """
    type = 'RADBC'

    def __init__(self, card=None, data=None):
        ThermalBC.__init__(self, card, data)

        ## NODAMB Ambient point for radiation exchange. (Integer > 0)
        self.nodamb = card.field(1)
        ## Radiation view factor between the face and the ambient point.
        ## (Real > 0.0)
        self.famb = card.field(2)
        ## Control point for thermal flux load. (Integer > 0; Default = 0)
        self.cntrlnd = card.field(3, 0)

        eids = card.fields(4)
        ## CHBDYi element identification number
        self.eids = expand_thru_by(eids)

    def cross_reference(self, model):
        for i, eid in enumerate(self.eids):
            self.eids[i] = model.Element(eid)

    def Eids(self):
        eids = []
        for eid in self.eids:
            eids.append(self.Eid(eid))
        return eids

    def Eid(self, eid):
        if isinstance(eid, int):
            return eid
        return eid.eid

    def rawFields(self):
        fields = ['RADBC', self.nodamb, self.famb, self.cntrlnd] + self.Eids()
        return fields

    def reprFields(self):
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)
        eids = collapse_thru_by(self.Eids())
        fields = ['RADBC', self.nodamb, self.famb, cntrlnd] + eids
        return fields

# Boundary Conditions
#-------------------------------------------------------
