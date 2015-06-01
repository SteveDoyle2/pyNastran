# pylint: disable=C0103,R0902,R0904,R0914,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import integer_types
from six.moves import range

from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.baseCard import (BaseCard, expand_thru_by,
                                          collapse_thru_by)
from pyNastran.bdf.bdfInterface.assign_type import (fields, integer, double,
    integer_or_blank, double_or_blank, integer_or_string, string, blank)

class ThermalCard(BaseCard):
    def __init__(self, card, data):
        pass

    def cross_reference(self, model):
        raise NotImplementedError('%s has not defined the cross_reference '
                                  'method' % self.type)

    def _is_same_card(self, obj, debug=False):
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
        if isinstance(self.pid, integer_types):
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

    hexMap = {
        1: [4, 3, 2, 1],
        2: [1, 2, 6, 5],
        3: [2, 3, 7, 6],
        4: [3, 4, 8, 7],
        5: [4, 1, 5, 8],
        6: [5, 6, 7, 8],
    }

    pentMap = {
        1: [3, 2, 1],
        2: [1, 2, 5, 4],
        3: [2, 3, 6, 5],
        4: [3, 1, 4, 6],
        5: [4, 5, 6],
    }

    tetMap = {
        1: [1, 3, 2],
        2: [1, 2, 4],
        3: [2, 3, 4],
        4: [3, 1, 4],
    }

    sideMaps = {
        'CHEXA': hexMap,
        'CPENTA': pentMap,
        'CTETRA': tetMap,
        'CTRIA3': [1, 2, 3],
        'CQUAD4': [1, 2, 3, 4],
    }

    def __init__(self, card=None, data=None, comment=''):
        ThermalElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Surface element ID number for a side of an
            #: element. (0 < Integer < 100,000,000)
            self.eid = integer(card, 1, 'eid')

            #: A heat conduction element identification
            self.eid2 = integer(card, 2, 'eid2')

            #: A consistent element side identification number
            #: (1 < Integer < 6)
            self.side = integer(card, 3, 'side')
            assert 0 < self.side < 7

            #: A VIEW entry identification number for the front face
            self.iViewFront = integer_or_blank(card, 4, 'iViewFront', 0)
            #: A VIEW entry identification number for the back face
            self.iViewBack = integer_or_blank(card, 5, 'iViewBack', 0)

            #: RADM identification number for front face of surface element
            #: (Integer > 0)
            self.radMidFront = integer_or_blank(card, 6, 'radMidFront', 0)
            #: RADM identification number for back face of surface element
            #: (Integer > 0)
            self.radMidBack = integer_or_blank(card, 7, 'radMidBack', 0)
            assert len(card) <= 8, 'len(CHBDYE card) = %i' % len(card)
        else:
            raise NotImplementedError(data)
        self.grids = []

    def cross_reference(self, model):
        pass

    def _verify(self, xref=False):
        eid = self.Eid()
        eid2 = self.Eid2()
        pid = self.Pid()
        assert isinstance(eid, int)
        assert isinstance(eid2, int)
        assert isinstance(pid, int)

    # def sideToEIDs(self, eid):
    #     sideIDs = self.sideMaps[eid.type][self.side]
    #     # [1,2,3]
    #
    #     # id-1 is for the 0 based python index
    #     nodes = [enodes[id - 1] for id in range(len(eid.nodes))
    #              if id in sideIDs]
    #     return side

    def Eid(self):
        return self.eid

    def Eid2(self):
        return self.eid2

    def raw_fields(self):
        list_fields = ['CHBDYE', self.eid, self.eid2, self.side,
                       self.iViewFront, self.iViewBack, self.radMidFront,
                       self.radMidBack]
        return list_fields

    def repr_fields(self):
        """
        .. todo:: is this done
        """
        #eids = collapse_thru_by(self.eids)
        list_fields = ['CHBDYE', self.eid, self.eid2, self.side,
                       self.iViewFront, self.iViewBack, self.radMidFront,
                       self.radMidBack]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class CHBDYG(ThermalElement):
    """
    Defines a boundary condition surface element without reference to a
    property entry.
    """
    type = 'CHBDYG'

    def __init__(self, card=None, data=None, comment=''):
        ThermalElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Surface element ID
            self.eid = integer(card, 1, 'eid')
            # no field 2

            #: Surface type
            self.Type = string(card, 3, 'Type')
            assert self.Type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8']

            #: A VIEW entry identification number for the front face
            self.iViewFront = integer_or_blank(card, 4, 'iViewFront', 0)

            #: A VIEW entry identification number for the back face
            self.iViewBack = integer_or_blank(card, 8, 'iViewBack', 0)

            #: RADM identification number for front face of surface element
            #: (Integer > 0)
            self.radMidFront = integer_or_blank(card, 6, 'radMidFront', 0)

            #: RADM identification number for back face of surface element
            #: (Integer > 0)
            self.radMidBack = integer_or_blank(card, 7, 'radMidBack', 0)
            # no field 8

            #: Grid point IDs of grids bounding the surface (Integer > 0)
            self.grids = []
            n = 1
            for i in range(9, len(card)):
                grid = integer_or_blank(card, i, 'grid%i' % n)
                if grid is not None:
                    self.grids.append(grid)
            assert len(self.grids) > 0, 'card=%s' % card
        else:
            self.eid = data[0]
            self.Type = data[1]
            self.iViewFront = data[2]
            self.iViewBack = data[3]
            self.radMidFront = data[4]
            self.radMidBack = data[5]
            self.grids = data[6:14]

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        assert isinstance(eid, int)
        assert isinstance(pid, int)

    def cross_reference(self, model):
        pass
        #msg = ' which is required by CHBDYG eid=%s' % self.eid
        #self.pid = model.Phbdy(self.pid, msg=msg)
        #self.grids

    def Eid(self):
        return self.eid

    def raw_fields(self):
        list_fields = (['CHBDYG', self.eid, None, self.Type, self.iViewFront,
                        self.iViewBack, self.radMidFront, self.radMidBack, None,] +
                       self.grids)
        return list_fields

    def repr_fields(self):
        iViewFront = set_blank_if_default(self.iViewFront, 0)
        iViewBack = set_blank_if_default(self.iViewBack, 0)
        radMidFront = set_blank_if_default(self.radMidFront, 0)
        radMidBack = set_blank_if_default(self.radMidBack, 0)

        list_fields = (['CHBDYG', self.eid, None, self.Type, iViewFront,
                        iViewBack, radMidFront, radMidBack, None, ] + self.grids)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class CHBDYP(ThermalElement):
    """
    Defines a boundary condition surface element with reference to a PHBDY
    entry
    """
    type = 'CHBDYP'

    def __init__(self, card=None, data=None, comment=''):
        ThermalElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Surface element ID
            self.eid = integer(card, 1, 'eid')

            #: PHBDY property entry identification numbers. (Integer > 0)
            self.pid = integer(card, 2, 'pid')
            assert self.pid > 0

            self.Type = string(card, 3, 'Type')
            #print("self.Type = %s" % self.Type)
            # msg = 'CHBDYP Type=%r' % self.Type
            #assert self.Type in ['POINT','LINE','ELCYL','FTUBE','TUBE'], msg

            #: A VIEW entry identification number for the front face.
            self.iViewFront = integer_or_blank(card, 4, 'iViewFront', 0)

            #: A VIEW entry identification number for the back face.
            self.iViewBack = integer_or_blank(card, 5, 'iViewBack', 0)

            #: Grid point identification numbers of grids bounding the surface.
            #: (Integer > 0)
            self.g1 = integer(card, 6, 'g1')
            #: Grid point identification numbers of grids bounding the surface.
            #: (Integer > 0)
            if self.Type != 'POINT':
                self.g2 = integer(card, 7, 'g2')
            else:
                self.g2 = blank(card, 7, 'g2')

            #: Orientation grid point. (Integer > 0; Default = 0)
            self.g0 = integer_or_blank(card, 8, 'g0', 0)

            #: RADM identification number for front face of surface element.
            #: (Integer > 0)
            self.radMidFront = integer_or_blank(card, 9, 'radMidFront', 0)

            #: RADM identification number for back face of surface element.
            #: (Integer > 0)
            self.radMidBack = integer_or_blank(card, 10, 'radMidBack', 0)

            #: Grid point identification number of a midside node if it is used
            #: with the line type surface element.
            self.gmid = integer_or_blank(card, 11, 'gmid')
            #: Coordinate system for defining orientation vector.
            #: (Integer > 0; Default = 0
            self.ce = integer_or_blank(card, 12, 'ce', 0)

            #: Components of the orientation vector in coordinate system CE.
            #: The origin of the orientation vector is grid point G1.
            #: (Real or blank)
            self.e1 = double_or_blank(card, 13, 'e3')
            self.e2 = double_or_blank(card, 14, 'e3')
            self.e3 = double_or_blank(card, 15, 'e3')
            assert len(card) <= 16, 'len(CHBDYP card) = %i' % len(card)
        else:
            raise NotImplementedError(data)

    def cross_reference(self, model):
        msg = ' which is required by CHBDYP pid=%s' % self.pid
        self.pid = model.Phbdy(self.pid, msg=msg)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        assert isinstance(eid, int)
        assert isinstance(pid, int)

    def Eid(self):
        return self.eid

    def raw_fields(self):
        list_fields = ['CHBDYP', self.eid, self.Pid(), self.Type,
                       self.iViewFront, self.iViewBack, self.g1, self.g2, self.g0,
                       self.radMidFront, self.radMidBack, self.gmid, self.ce,
                       self.e1, self.e2, self.e3]
        return list_fields

    def repr_fields(self):
        iViewFront = set_blank_if_default(self.iViewFront, 0)
        iViewBack = set_blank_if_default(self.iViewBack, 0)
        radMidFront = set_blank_if_default(self.radMidFront, 0)
        radMidBack = set_blank_if_default(self.radMidBack, 0)

        g0 = set_blank_if_default(self.g0, 0)
        ce = set_blank_if_default(self.ce, 0)

        list_fields = ['CHBDYP', self.eid, self.Pid(), self.Type, iViewFront,
                       iViewBack, self.g1, self.g2, g0, radMidFront, radMidBack,
                       self.gmid, ce, self.e1, self.e2, self.e3]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)

# Elements
#-------------------------------------------------------
# Properties


class PCONV(ThermalProperty):
    """
    Specifies the free convection boundary condition properties of a boundary
    condition surface element used for heat transfer analysis.
    """
    type = 'PCONV'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Convection property identification number. (Integer > 0)
            self.pconid = integer(card, 1, 'pconid')
            assert self.pconid > 0

            #: Material property identification number. (Integer > 0)
            self.mid = integer(card, 2, 'mid')
            assert self.mid > 0

            #: Type of formula used for free convection.
            #: (Integer 0, 1, 10, 11, 20, or 21)
            self.form = integer_or_blank(card, 3, 'form', 0)
            assert self.form in [0, 1, 10, 11, 20, 21]

            #: Free convection exponent as implemented within the context of the
            #: particular form that is chosen
            self.expf = double_or_blank(card, 4, 'expf', 0.0)

            #: Formula type for various configurations of free convection
            self.ftype = integer_or_blank(card, 5, 'ftype', 0)

            #: Identification number of a TABLEHT entry that specifies the two
            #: variable tabular function of the free convection heat transfer
            #: coefficient
            self.tid = integer_or_blank(card, 6, 'tid')

            #: Characteristic length
            self.chlen = double_or_blank(card, 9, 'chlen')

            #: Grid ID of the referenced inlet point
            self.gidin = double_or_blank(card, 10, 'gidin')

            #: Coordinate system for defining orientation vector.
            #: (Integer > 0;Default = 0
            self.ce = integer_or_blank(card, 11, 'ce', 0)

            #: Components of the orientation vector in coordinate system CE. The
            #: origin of the orientation vector is grid point G1. (Real or blank)
            self.e1 = double_or_blank(card, 12, 'e1')
            self.e2 = double_or_blank(card, 13, 'e2')
            self.e3 = double_or_blank(card, 14, 'e3')
            assert len(card) <= 15, 'len(PCONV card) = %i' % len(card)
        else:
            raise NotImplementedError(data)

    #def cross_reference(self,model):
    #    pass

    def raw_fields(self):
        list_fields = ['PCONV', self.pconid, self.mid, self.form, self.expf,
                       self.ftype, self.tid, None, None, self.chlen, self.gidin,
                       self.ce, self.e1, self.e2, self.e3]
        return list_fields

    def repr_fields(self):
        form = set_blank_if_default(self.form, 0)
        expf = set_blank_if_default(self.expf, 0.0)
        ftype = set_blank_if_default(self.ftype, 0)
        ce = set_blank_if_default(self.ce, 0)
        list_fields = ['PCONV', self.pconid, self.mid, form, expf, ftype, self.tid,
                       None, None, self.chlen, self.gidin, ce, self.e1, self.e2,
                       self.e3]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class PCONVM(ThermalProperty):
    """
    Specifies the free convection boundary condition properties of a boundary
    condition surface element used for heat transfer analysis.
    """
    type = 'PCONVM'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Convection property identification number. (Integer > 0)
            self.pconid = integer(card, 1, 'pconid')
            assert self.pconid > 0

            #: Material property identification number. (Integer > 0)
            self.mid = integer(card, 2, 'mid')
            assert self.mid > 0

            #: Type of formula used for free convection.
            #: (Integer 0, 1, 10, 11, 20, or 21)
            self.form = integer_or_blank(card, 3, 'form', 0)
            assert self.form in [0, 1, 10, 11, 20, 21]

            #: Flag for mass flow convection. (Integer = 0 or 1; Default = 0)
            self.flag = integer_or_blank(card, 4, 'flag', 0)

            #: Constant coefficient used for forced convection
            self.coef = double(card, 5, 'coef')

            #: Reynolds number convection exponent. (Real > 0.0; Default = 0.0)
            self.expr = double_or_blank(card, 6, 'expr', 0.0)

            #: Prandtl number convection exponent for heat transfer into the
            #: workingfluid. (Real > 0.0; Default = 0.0)
            self.exppi = double_or_blank(card, 7, 'exppi', 0.0)

            #: Prandtl number convection exponent for heat transfer into the
            #: working fluid. (Real > 0.0; Default = 0.0)
            self.exppo = double_or_blank(card, 8, 'exppo', 0.0)
            assert len(card) <= 9, 'len(PCONVM card) = %i' % len(card)
        else:
            raise NotImplementedError(data)
    #def cross_reference(self,model):
    #    pass

    def raw_fields(self):
        list_fields = ['PCONVM', self.pconid, self.mid, self.form,
                       self.flag, self.coef, self.expr, self.exppi, self.exppo]
        return list_fields

    def repr_fields(self):
        form = set_blank_if_default(self.form, 0)
        flag = set_blank_if_default(self.flag, 0)
        expr = set_blank_if_default(self.expr, 0.0)
        exppi = set_blank_if_default(self.exppi, 0.0)
        exppo = set_blank_if_default(self.exppo, 0.0)
        list_fields = ['PCONVM', self.pconid, self.mid, form, flag,
                       self.coef, expr, exppi, exppo]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class PHBDY(ThermalProperty):
    """
    A property entry referenced by CHBDYP entries to give auxiliary geometric
    information for boundary condition surface elements
    """
    type = 'PHBDY'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Property identification number. (Unique Integer among all PHBDY
            #: entries). (Integer > 0)
            self.pid = integer(card, 1, 'pid')
            assert self.pid > 0

            #: Area factor of the surface used only for CHBDYP element
            #: TYPE = 'POINT', TYPE = 'LINE', TYPE = 'TUBE', or
            #: TYPE = 'ELCYL'. For TYPE = 'TUBE', AF is the constant thickness
            #: of the hollow tube. (Real > 0.0 or blank)
            self.af = double_or_blank(card, 2, 'af')

            #: Diameters associated with the surface. Used with CHBDYP element
            #: TYPE='ELCYL','TUBE','FTUBE'
            self.d1 = double_or_blank(card, 3, 'd1')
            self.d2 = double_or_blank(card, 4, 'd2', self.d1)
            assert len(card) <= 5, 'len(PHBDY card) = %i' % len(card)
        else:
            raise NotImplementedError(data)

    #def cross_reference(self,model):
    #    pass

    def raw_fields(self):
        list_fields = ['PHBDY', self.pid, self.af, self.d1, self.d2]
        return list_fields

    def repr_fields(self):
        d2 = set_blank_if_default(self.d2, self.d1)
        list_fields = ['PHBDY', self.pid, self.af, self.d1, d2]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


# Properties
#-------------------------------------------------------
# Boundary Conditions

class CONV(ThermalBC):
    """
    Specifies a free convection boundary condition for heat transfer analysis
    through connection to a surface element (CHBDYi entry).
    """
    type = 'CONV'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        #ThermalBC.__init__(self, card, data)
        if card:
            #: CHBDYG, CHBDYE, or CHBDYP surface element identification number.
            #: (Integer > 0)
            self.eid = integer(card, 1, 'eid')
            assert self.eid > 0

            #: Convection property identification number of a PCONV entry
            self.pconID = integer(card, 2, 'pconID')

            #: Point for film convection fluid property temperature
            self.flmnd = integer_or_blank(card, 3, 'flmnd', 0)

            #: Control point for free convection boundary condition.
            self.cntrlnd = integer_or_blank(card, 4, 'cntrlnd', 0)

            TA1 = integer(card, 5, 'TA1')
            assert TA1 > 0

            #: Ambient points used for convection 0's are allowed for TA2 and
            #: higher.  (Integer > 0 for TA1 and Integer > 0 for TA2 through TA8;
            #: Default for TA2 through TA8 is TA1.)
            TA2 = integer_or_blank(card, 6, 'ta2', TA1)
            TA3 = integer_or_blank(card, 7, 'ta3', TA1)
            TA4 = integer_or_blank(card, 8, 'ta4', TA1)
            TA5 = integer_or_blank(card, 9, 'ta5', TA1)
            TA6 = integer_or_blank(card, 10, 'ta6', TA1)
            TA7 = integer_or_blank(card, 11, 'ta7', TA1)
            TA8 = integer_or_blank(card, 12, 'ta8', TA1)
            self.ta = [TA1, TA2, TA3, TA4, TA5, TA6, TA7, TA8]
            assert len(card) <= 13, 'len(CONV card) = %i' % len(card)
        else:
            raise NotImplementedError(data)

    def cross_reference(self, model):
        msg = ' which is required by CONV eid=%s' % self.eid
        self.eid = model.Element(self.eid, msg=msg)
        if model._xref == 1:  # True
            assert self.eid.type in ['CHBDYG', 'CHBDYE', 'CHBDYP']

    def TA(self, i=None):
        if i is None:
            return self.ta
        return self.ta[i]

    def raw_fields(self):
        list_fields = ['CONV', self.eid, self.pconID, self.flmnd,
                       self.cntrlnd] + self.ta
        return list_fields

    def repr_fields(self):
        flmnd = set_blank_if_default(self.flmnd, 0)
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)

        ta0 = self.ta[0]
        ta = [ta0]
        for tai in self.ta[1:]:
            ta.append(set_blank_if_default(tai, ta0))

        list_fields = ['CONV', self.eid, self.pconID, flmnd, cntrlnd] + ta
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class CONVM(ThermalBC):
    """
    Specifies a forced convection boundary condition for heat transfer analysis
    through connection to a surface element (CHBDYi entry).
    """
    type = 'CONV'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pconvmID = integer(card, 2, 'pconvmID')

            self.filmNode = integer_or_blank(card, 3, 'filmNode', 0)
            assert self.filmNode >= 0

            self.cntmdot = integer(card, 4, 'cntmdot')
            assert self.cntmdot > 0

            self.ta1 = integer(card, 5, 'ta1')
            self.ta2 = integer_or_blank(card, 6, 'ta2', self.ta1)

            self.mdot = double_or_blank(card, 7, 'mdot', 1.0)
            assert len(card) <= 8, 'len(CONVM card) = %i' % len(card)
        else:
            raise NotImplementedError(data)

    def cross_reference(self, model):
        msg = ' which is required by CONVM eid=%s' % self.eid
        self.eid = model.CYBDY(self.eid, msg=msg)
        self.pconvmID = model.PCONV(self.pconvmID, msg=msg)
        self.filmNode = model.Grid(self.filmNode, msg=msg)

    def film_node(self):
        if isinstance(self.filmNode, integer_types):
            return self.filmNode
        return self.filmNode.nid

    def raw_fields(self):
        list_fields = ['CONVM', self.eid, self.pconvmID, self.filmNode,
                       self.cntmdot, self.ta1, self.ta2, self.mdot]
        return list_fields

    def repr_fields(self):
        filmNode = set_blank_if_default(self.filmNode, 0)
        ta2 = set_blank_if_default(self.ta2, self.ta1)
        mdot = set_blank_if_default(self.mdot, 1.0)
        list_fields = ['CONVM', self.eid, self.pconvmID, filmNode,
                       self.cntmdot, self.ta1, ta2, mdot]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class RADM(ThermalBC):
    """
    Defines the radiation properties of a boundary element for heat transfer
    analysis
    """
    type = 'RADM'

    def __init__(self, card=None, data=None, comment=''):
        ThermalBC.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Material identification number
            self.radmid = integer(card, 1, 'radmid')
            assert self.radmid > 0

            self.absorb = double(card, 2, 'absorb')
            assert 0. <= self.absorb <= 1.0

            nfields = card.nfields

            self.emissivity = fields(double, card, 'emissivity', i=3, j=nfields)
        else:
            raise NotImplementedError(data)
        for e in self.emissivity:
            assert 0. <= e <= 1.0

    #def cross_reference(self,model):
    #    pass

    def repr_fields(self):
        list_fields = ['RADM', self.radmid, self.absorb] + self.emissivity
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class RADBC(ThermalBC):
    """
    Specifies an CHBDYi element face for application of radiation boundary
    conditions
    """
    type = 'RADBC'

    def __init__(self, card=None, data=None, comment=''):
        ThermalBC.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            #: NODAMB Ambient point for radiation exchange. (Integer > 0)
            self.nodamb = integer(card, 1, 'nodamb')
            assert self.nodamb > 0

            #: Radiation view factor between the face and the ambient point.
            #: (Real > 0.0)
            self.famb = double(card, 2, 'famb')
            assert self.famb > 0.0

            #: Control point for thermal flux load. (Integer > 0; Default = 0)
            self.cntrlnd = integer_or_blank(card, 3, 'cntrlnd', 0)
            assert self.cntrlnd >= 0

            nfields = card.nfields
            eids = fields(integer_or_string, card, 'eid', i=4, j=nfields)
            #: CHBDYi element identification number
            self.eids = expand_thru_by(eids)
        else:
            raise NotImplementedError(data)

        min_eid = min(self.eids)
        if min_eid < 1:
            msg = 'min(eids)=%i' % min_eid
            raise ValueError(msg)

    def cross_reference(self, model):
        msg = ' which is required by RADBC pid=%s' % self.nodamb
        for i, eid in enumerate(self.eids):
            self.eids[i] = model.Element(eid, msg=msg)

    def Eids(self):
        eids = []
        for eid in self.eids:
            eids.append(self.Eid(eid))
        return eids

    def Eid(self, eid):
        if isinstance(eid, integer_types):
            return eid
        return eid.eid

    def raw_fields(self):
        list_fields = (['RADBC', self.nodamb, self.famb, self.cntrlnd] +
                       self.Eids())
        return list_fields

    def repr_fields(self):
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)
        eids = collapse_thru_by(self.Eids())
        list_fields = ['RADBC', self.nodamb, self.famb, cntrlnd] + eids
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)

# Boundary Conditions
#-------------------------------------------------------
