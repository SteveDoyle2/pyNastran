# pylint: disable=R0902,R0904,R0914,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import range

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import (BaseCard, expand_thru_by,
                                           _node_ids)
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.assign_type import (
    fields, integer, double, integer_or_blank, double_or_blank,
    integer_or_string, string, blank)


class ThermalCard(BaseCard):
    def __init__(self):
        pass

    def cross_reference(self, model):
        raise NotImplementedError('%s has not defined the cross_reference '
                                  'method' % self.type)

    #def _is_same_card(self, obj, debug=False):
        #return False


class ThermalBC(ThermalCard):
    def __init__(self):
        ThermalCard.__init__(self)


class ThermalElement(ThermalCard):
    pid = 0

    def __init__(self):
        ThermalCard.__init__(self)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return []

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        else:
            return self.pid_ref.pid


class ThermalProperty(ThermalCard):
    def __init__(self):
        ThermalCard.__init__(self)

#-------------------------------------------------------
# Elements


class CHBDYE(ThermalElement):
    """
    Defines a boundary condition surface element with reference to a heat
    conduction element.

    +--------+-----+------+------+--------+--------+---------+---------+
    | CHBDYE | EID | EID2 | SIDE | IVIEWF | IVIEWB | RADMIDF | RADMIDB |
    +--------+-----+------+------+--------+--------+---------+---------+
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

    side_maps = {
        'CHEXA': hexMap,
        'CPENTA': pentMap,
        'CTETRA': tetMap,
        'CTRIA3': [1, 2, 3],
        'CQUAD4': [1, 2, 3, 4],
    }

    def __init__(self, eid, eid2, side, iViewFront=0, iViewBack=0,
                 radMidFront=0, radMidBack=0, comment=''):
        ThermalElement.__init__(self)
        if comment:
            self._comment = comment
        #: Surface element ID number for a side of an
        #: element. (0 < Integer < 100,000,000)
        self.eid = eid

        #: A heat conduction element identification
        self.eid2 = eid2

        #: A consistent element side identification number
        #: (1 < Integer < 6)
        self.side = side
        assert 0 < side < 7

        #: A VIEW entry identification number for the front face
        self.iViewFront = iViewFront

        #: A VIEW entry identification number for the back face
        self.iViewBack = iViewBack

        #: RADM identification number for front face of surface element
        #: (Integer > 0)
        self.radMidFront = radMidFront
        #: RADM identification number for back face of surface element
        #: (Integer > 0)
        self.radMidBack = radMidBack
        self.grids = []

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        eid2 = integer(card, 2, 'eid2')
        side = integer(card, 3, 'side')

        iViewFront = integer_or_blank(card, 4, 'iViewFront', 0)
        iViewBack = integer_or_blank(card, 5, 'iViewBack', 0)
        radMidFront = integer_or_blank(card, 6, 'radMidFront', 0)
        radMidBack = integer_or_blank(card, 7, 'radMidBack', 0)
        assert len(card) <= 8, 'len(CHBDYE card) = %i\ncard=%s' % (len(card), card)
        return CHBDYE(eid, eid2, side, iViewFront, iViewBack,
                      radMidFront, radMidBack, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid, eid2, side, iviewf, iviewb, radmidf, radmidb = data
        return CHBDYE(eid, eid2, side, iviewf, iviewb,
                      radmidf, radmidb, comment=comment)

    def cross_reference(self, model):
        pass

    def safe_cross_reference(self):
        pass

    def uncross_reference(self):
        pass

    @property
    def nodes(self):
        return []

    @property
    def node_ids(self):
        return _node_ids(self, nodes=self.grids, allow_empty_nodes=False, msg='')

    def get_edge_ids(self):
        # TODO: not implemented
        return []

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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class CHBDYG(ThermalElement):
    """
    Defines a boundary condition surface element without reference to a
    property entry.

    +--------+-----+----+------+--------+--------+---------+---------+-----+
    | CHBDYG | EID |    | TYPE | IVIEWF | IVIEWB | RADMIDF | RADMIDB |     |
    +--------+-----+----+------+--------+--------+---------+---------+-----+
    |        | G1  | G2 |  G3  |   G4   |   G5   |   G6    |   G7    |  G8 |
    +--------+-----+----+------+--------+--------+---------+---------+-----+
    """
    type = 'CHBDYG'

    def __init__(self, eid, Type, nodes, iViewFront=0, iViewBack=0,
                 radMidFront=0, radMidBack=0, comment=''):
        ThermalElement.__init__(self)
        if comment:
            self._comment = comment
        #: Surface element ID
        self.eid = eid

        #: Surface type
        self.Type = Type

        #: A VIEW entry identification number for the front face
        self.iViewFront = iViewFront

        #: A VIEW entry identification number for the back face
        self.iViewBack = iViewBack

        #: RADM identification number for front face of surface element
        #: (Integer > 0)
        self.radMidFront = radMidFront

        #: RADM identification number for back face of surface element
        #: (Integer > 0)
        self.radMidBack = radMidBack

        #: Grid point IDs of grids bounding the surface (Integer > 0)
        self.nodes = nodes

        assert self.Type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], 'Type=%r' % Type
        assert len(nodes) > 0, nodes

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        # no field 2

        Type = string(card, 3, 'Type')
        iViewFront = integer_or_blank(card, 4, 'iViewFront', 0)
        iViewBack = integer_or_blank(card, 8, 'iViewBack', 0)
        radMidFront = integer_or_blank(card, 6, 'radMidFront', 0)
        radMidBack = integer_or_blank(card, 7, 'radMidBack', 0)
        # no field 8

        n = 1
        nodes = []
        for i in range(9, len(card)):
            grid = integer_or_blank(card, i, 'grid%i' % n)
            if grid is not None:
                nodes.append(grid)
        assert len(nodes) > 0, 'card=%s' % card
        return CHBDYG(eid, Type, nodes,
                      iViewFront=iViewFront, iViewBack=iViewBack,
                      radMidFront=radMidFront, radMidBack=radMidBack,
                      comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        Type = data[1]
        iViewFront = data[2]
        iViewBack = data[3]
        radMidFront = data[4]
        radMidBack = data[5]
        nodes = [datai for datai in data[6:14] if datai > 0]
        if Type == 5:
            Type = 'AREA4'
        elif Type == 4:
            Type = 'AREA3'
        #elif Type == 7:
            #Type = 'AREA6'
        elif Type == 9:
            Type = 'AREA8'
        assert Type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], 'Type=%r data=%s' % (Type, data)
        return CHBDYG(eid, Type, nodes,
                      iViewFront=iViewFront, iViewBack=iViewBack,
                      radMidFront=radMidFront, radMidBack=radMidBack,
                      comment=comment)

    def _verify(self, xref=False):
        eid = self.Eid()
        assert isinstance(eid, int)

    @property
    def node_ids(self):
        # TODO: is this correct?
        return _node_ids(self, nodes=self.nodes, allow_empty_nodes=False, msg='')

    def get_edge_ids(self):
        # TODO: not implemented
        return []

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CHBDYG eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=False, msg=msg)
        self.nodes_ref = self.nodes

    def safe_cross_reference(self, model):
        msg = ' which is required by CHBDYG eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=False, msg=msg)
        self.nodes_ref = self.nodes

    def uncross_reference(self):
        self.nodes = self.node_ids
        del self.nodes_ref

    def Eid(self):
        return self.eid

    def raw_fields(self):
        list_fields = (
            ['CHBDYG', self.eid, None, self.Type, self.iViewFront,
             self.iViewBack, self.radMidFront, self.radMidBack, None,] +
            self.node_ids)
        return list_fields

    def repr_fields(self):
        iViewFront = set_blank_if_default(self.iViewFront, 0)
        iViewBack = set_blank_if_default(self.iViewBack, 0)
        radMidFront = set_blank_if_default(self.radMidFront, 0)
        radMidBack = set_blank_if_default(self.radMidBack, 0)

        list_fields = (['CHBDYG', self.eid, None, self.Type, iViewFront,
                        iViewBack, radMidFront, radMidBack, None, ] + self.node_ids)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class CHBDYP(ThermalElement):
    """
    Defines a boundary condition surface element with reference to a PHBDY
    entry

    +--------+---------+---------+------+--------+--------+----+----+----+
    | CHBDYP |   EID   |   PID   | TYPE | IVIEWF | IVIEWB | G1 | G2 | G0 |
    +--------+---------+---------+------+--------+--------+----+----+----+
    |        | RADMIDF | RADMIDB | GMID |   CE   |   E1   | E2 | E3 |    |
    +--------+---------+---------+------+--------+--------+----+----+----+
    """
    type = 'CHBDYP'

    def __init__(self, eid, pid, Type, g1, g2, g0=0, gmid=None, ce=0,
                 iViewFront=0, iViewBack=0,
                 radMidFront=0, radMidBack=0,
                 e1=None, e2=None, e3=None, comment=''):
        ThermalElement.__init__(self)
        if comment:
            self._comment = comment
        #: Surface element ID
        self.eid = eid

        #: PHBDY property entry identification numbers. (Integer > 0)
        self.pid = pid
        self.Type = Type

        #: A VIEW entry identification number for the front face.
        self.iViewFront = iViewFront

        #: A VIEW entry identification number for the back face.
        self.iViewBack = iViewBack

        #: Grid point identification numbers of grids bounding the surface.
        #: (Integer > 0)
        self.g1 = g1
        #: Grid point identification numbers of grids bounding the surface.
        #: (Integer > 0)
        self.g2 = g2

        #: Orientation grid point. (Integer > 0; Default = 0)
        self.g0 = g0

        #: RADM identification number for front face of surface element.
        #: (Integer > 0)
        self.radMidFront = radMidFront

        #: RADM identification number for back face of surface element.
        #: (Integer > 0)
        self.radMidBack = radMidBack

        #: Grid point identification number of a midside node if it is used
        #: with the line type surface element.
        self.gmid = gmid

        #: Coordinate system for defining orientation vector.
        #: (Integer > 0; Default = 0
        self.ce = ce

        #: Components of the orientation vector in coordinate system CE.
        #: The origin of the orientation vector is grid point G1.
        #: (Real or blank)
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3

        assert self.pid > 0

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')

        Type = string(card, 3, 'Type')
        #print("self.Type = %s" % self.Type)
        # msg = 'CHBDYP Type=%r' % self.Type
        #assert self.Type in ['POINT','LINE','ELCYL','FTUBE','TUBE'], msg

        iViewFront = integer_or_blank(card, 4, 'iViewFront', 0)
        iViewBack = integer_or_blank(card, 5, 'iViewBack', 0)
        g1 = integer(card, 6, 'g1')

        if Type != 'POINT':
            g2 = integer(card, 7, 'g2')
        else:
            g2 = blank(card, 7, 'g2')

        g0 = integer_or_blank(card, 8, 'g0', 0)
        radMidFront = integer_or_blank(card, 9, 'radMidFront', 0)
        radMidBack = integer_or_blank(card, 10, 'radMidBack', 0)
        gmid = integer_or_blank(card, 11, 'gmid')
        ce = integer_or_blank(card, 12, 'ce', 0)
        e1 = double_or_blank(card, 13, 'e3')
        e2 = double_or_blank(card, 14, 'e3')
        e3 = double_or_blank(card, 15, 'e3')
        assert len(card) <= 16, 'len(CHBDYP card) = %i\ncard=%s' % (len(card), card)
        return CHBDYP(eid, pid, Type, g1, g2, g0=g0, gmid=gmid, ce=ce,
                      iViewFront=iViewFront, iViewBack=iViewBack,
                      radMidFront=radMidFront, radMidBack=radMidBack,
                      e1=e1, e2=e2, e3=e3, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        [eid, pid, Type, iviewf, iviewb, g1, g2, g0, radmidf, radmidb,
         dislin, ce, e1, e2, e3] = data
        #eid = data[0]
        #Type = data[1]
        #iViewFront = data[2]
        #iViewBack = data[3]
        #radMidFront = data[4]
        #radMidBack = data[5]
        #nodes = [datai for datai in data[6:14] if datai > 0]
        if Type == 1:
            Type = 'POINT'
        elif Type == 2:
            Type = 'LINE'
        elif Type == 7:
            Type = 'FTUBE'
        else:
            raise NotImplementedError('CHBDYP Type=%r data=%s' % (Type, data))
        #assert Type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], 'Type=%r data=%s' % (Type, data)

        assert dislin == 0, 'CHBDYP dislin=%r data=%s' % (dislin, data)
        gmid = dislin
        return CHBDYP(eid, pid, Type, g1, g2, g0=g0, gmid=gmid, ce=ce,
                      iViewFront=iviewf, iViewBack=iviewb,
                      radMidFront=radmidf, radMidBack=radmidb,
                      e1=e1, e2=e2, e3=e3, comment=comment)

    @property
    def nodes(self):
        return [self.g1, self.g2, self.g0, self.gmid]

    @nodes.setter
    def nodes(self, nodes):
        self.g1 = nodes[0]
        self.g2 = nodes[1]
        self.g0 = nodes[2]
        self.gmid = nodes[3]
        assert len(nodes) == 4, len(nodes)

    @property
    def node_ids(self):
        return _node_ids(self, nodes=self.nodes, allow_empty_nodes=True, msg='')

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CHBDYP pid=%s' % self.pid
        self.pid = model.Phbdy(self.pid, msg=msg)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.pid_ref = self.pid
        self.nodes_ref = self.nodes

    def safe_cross_reference(self, model):
        msg = ' which is required by CHBDYP pid=%s' % self.pid
        self.pid = model.Phbdy(self.pid, msg=msg)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.pid_ref = self.pid
        self.nodes_ref = self.nodes

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        assert isinstance(eid, int)
        assert isinstance(pid, int)

    def Eid(self):
        return self.eid

    def raw_fields(self):
        (g1, g2, g0, gmid) = self.node_ids
        list_fields = ['CHBDYP', self.eid, self.Pid(), self.Type,
                       self.iViewFront, self.iViewBack, g1, g2, g0,
                       self.radMidFront, self.radMidBack, gmid, self.ce,
                       self.e1, self.e2, self.e3]
        return list_fields

    def repr_fields(self):
        iViewFront = set_blank_if_default(self.iViewFront, 0)
        iViewBack = set_blank_if_default(self.iViewBack, 0)
        radMidFront = set_blank_if_default(self.radMidFront, 0)
        radMidBack = set_blank_if_default(self.radMidBack, 0)

        (g1, g2, g0, gmid) = self.node_ids
        g0 = set_blank_if_default(g0, 0)
        ce = set_blank_if_default(self.ce, 0)

        list_fields = ['CHBDYP', self.eid, self.Pid(), self.Type, iViewFront,
                       iViewBack, g1, g2, g0, radMidFront, radMidBack,
                       gmid, ce, self.e1, self.e2, self.e3]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

# Elements
#-------------------------------------------------------
# Properties


class PCONV(ThermalProperty):
    """
    Specifies the free convection boundary condition properties of a boundary
    condition surface element used for heat transfer analysis.
    """
    type = 'PCONV'

    def __init__(self, pconid, mid, form, expf, ftype, tid, chlen, gidin, ce,
                 e1, e2, e3, comment=''):
        ThermalProperty.__init__(self)
        if comment:
            self._comment = comment
        #: Convection property identification number. (Integer > 0)
        self.pconid = pconid

        #: Material property identification number. (Integer > 0)
        self.mid = mid

        #: Type of formula used for free convection.
        #: (Integer 0, 1, 10, 11, 20, or 21)
        self.form = form

        #: Free convection exponent as implemented within the context of the
        #: particular form that is chosen
        self.expf = expf

        #: Formula type for various configurations of free convection
        self.ftype = ftype

        #: Identification number of a TABLEHT entry that specifies the two
        #: variable tabular function of the free convection heat transfer
        #: coefficient
        self.tid = tid

        #: Characteristic length
        self.chlen = chlen

        #: Grid ID of the referenced inlet point
        self.gidin = gidin

        #: Coordinate system for defining orientation vector.
        #: (Integer > 0;Default = 0
        self.ce = ce

        #: Components of the orientation vector in coordinate system CE. The
        #: origin of the orientation vector is grid point G1. (Real or blank)
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3
        assert self.pconid > 0
        assert self.mid > 0
        assert self.form in [0, 1, 10, 11, 20, 21]

    @classmethod
    def add_card(cls, card, comment=''):
        pconid = integer(card, 1, 'pconid')
        mid = integer(card, 2, 'mid')
        form = integer_or_blank(card, 3, 'form', 0)
        expf = double_or_blank(card, 4, 'expf', 0.0)
        ftype = integer_or_blank(card, 5, 'ftype', 0)
        tid = integer_or_blank(card, 6, 'tid')
        chlen = double_or_blank(card, 9, 'chlen')
        gidin = double_or_blank(card, 10, 'gidin')
        ce = integer_or_blank(card, 11, 'ce', 0)
        e1 = double_or_blank(card, 12, 'e1')
        e2 = double_or_blank(card, 13, 'e2')
        e3 = double_or_blank(card, 14, 'e3')
        assert len(card) <= 15, 'len(PCONV card) = %i\ncard=%s' % (len(card), card)
        return PCONV(pconid, mid, form, expf, ftype, tid, chlen, gidin, ce,
                     e1, e2, e3, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        (pconid, mid, form, expf, ftype, tid, chlen, gidin, ce, e1, e2, e3) = data
        return PCONV(pconid, mid, form, expf, ftype, tid, chlen, gidin, ce,
                     e1, e2, e3, comment=comment)

    #def cross_reference(self, model):
        #pass

    def uncross_reference(self):
        pass

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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PCONVM(ThermalProperty):
    """
    Specifies the free convection boundary condition properties of a boundary
    condition surface element used for heat transfer analysis.
    """
    type = 'PCONVM'

    def __init__(self, pconid, mid, form, flag, coef, expr, exppi, exppo,
                 comment=''):
        ThermalProperty.__init__(self)
        if comment:
            self._comment = comment
        #: Convection property identification number. (Integer > 0)
        self.pconid = pconid
        assert self.pconid > 0

        #: Material property identification number. (Integer > 0)
        self.mid = mid
        assert self.mid > 0

        #: Type of formula used for free convection.
        #: (Integer 0, 1, 10, 11, 20, or 21)
        self.form = form
        assert self.form in [0, 1, 10, 11, 20, 21]

        #: Flag for mass flow convection. (Integer = 0 or 1; Default = 0)
        self.flag = flag

        #: Constant coefficient used for forced convection
        self.coef = coef

        #: Reynolds number convection exponent. (Real > 0.0; Default = 0.0)
        self.expr = expr

        #: Prandtl number convection exponent for heat transfer into the
        #: workingfluid. (Real > 0.0; Default = 0.0)
        self.exppi = exppi

        #: Prandtl number convection exponent for heat transfer into the
        #: working fluid. (Real > 0.0; Default = 0.0)
        self.exppo = exppo

    @classmethod
    def add_card(cls, card, comment=''):
        pconid = integer(card, 1, 'pconid')
        mid = integer(card, 2, 'mid')
        form = integer_or_blank(card, 3, 'form', 0)
        flag = integer_or_blank(card, 4, 'flag', 0)
        coef = double(card, 5, 'coef')
        expr = double_or_blank(card, 6, 'expr', 0.0)
        exppi = double_or_blank(card, 7, 'exppi', 0.0)
        exppo = double_or_blank(card, 8, 'exppo', 0.0)
        assert len(card) <= 9, 'len(PCONVM card) = %i\ncard=%s' % (len(card), card)
        return PCONVM(pconid, mid, form, flag, coef, expr, exppi, exppo,
                      comment=comment)

    #def cross_reference(self, model):
        #pass

    def uncross_reference(self):
        pass

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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PHBDY(ThermalProperty):
    """
    A property entry referenced by CHBDYP entries to give auxiliary geometric
    information for boundary condition surface elements
    """
    type = 'PHBDY'

    def __init__(self, pid, af, d1, d2, comment=''):
        ThermalProperty.__init__(self)
        if comment:
            self._comment = comment

        #: Property identification number. (Unique Integer among all PHBDY
        #: entries). (Integer > 0)
        self.pid = pid
        assert self.pid > 0

        #: Area factor of the surface used only for CHBDYP element
        #: TYPE = 'POINT', TYPE = 'LINE', TYPE = 'TUBE', or
        #: TYPE = 'ELCYL'. For TYPE = 'TUBE', AF is the constant thickness
        #: of the hollow tube. (Real > 0.0 or blank)
        self.af = af

        #: Diameters associated with the surface. Used with CHBDYP element
        #: TYPE='ELCYL','TUBE','FTUBE'
        self.d1 = d1
        self.d2 = d2

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        af = double_or_blank(card, 2, 'af')
        d1 = double_or_blank(card, 3, 'd1')
        d2 = double_or_blank(card, 4, 'd2', d1)
        assert len(card) <= 5, 'len(PHBDY card) = %i\ncard=%s' % (len(card), card)
        return PHBDY(pid, af, d1, d2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        pid = data[0]
        af = data[1]
        d1 = data[2]
        d2 = data[3]
        assert len(data) == 4, 'data = %s' % data
        return PHBDY(pid, af, d1, d2, comment=comment)

    #def cross_reference(self, model):
        #pass

    def uncross_reference(self):
        pass

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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


# Properties
#-------------------------------------------------------
# Boundary Conditions

class CONV(ThermalBC):
    """
    Specifies a free convection boundary condition for heat transfer analysis
    through connection to a surface element (CHBDYi entry).
    """
    type = 'CONV'

    def __init__(self, eid, pconid, film_node, cntrlnd, ta, comment=''):
        ThermalBC.__init__(self)
        if comment:
            self._comment = comment

        #: CHBDYG, CHBDYE, or CHBDYP surface element identification number.
        #: (Integer > 0)
        self.eid = eid

        #: Convection property identification number of a PCONV entry
        self.pconid = pconid

        #: Point for film convection fluid property temperature
        self.film_node = film_node

        #: Control point for free convection boundary condition.
        self.cntrlnd = cntrlnd

        #: Ambient points used for convection 0's are allowed for TA2 and
        #: higher.  (Integer > 0 for TA1 and Integer > 0 for TA2 through TA8;
        #: Default for TA2 through TA8 is TA1.)
        self.ta = ta
        assert self.eid > 0, 'eid=%s\n%s' % (eid, str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pconid = integer(card, 2, 'pconid')
        film_node = integer_or_blank(card, 3, 'film_node', 0)
        cntrlnd = integer_or_blank(card, 4, 'cntrlnd', 0)

        ta1 = integer(card, 5, 'TA1')
        assert ta1 > 0, ta1

        ta2 = integer_or_blank(card, 6, 'ta2', ta1)
        ta3 = integer_or_blank(card, 7, 'ta3', ta1)
        ta4 = integer_or_blank(card, 8, 'ta4', ta1)
        ta5 = integer_or_blank(card, 9, 'ta5', ta1)
        ta6 = integer_or_blank(card, 10, 'ta6', ta1)
        ta7 = integer_or_blank(card, 11, 'ta7', ta1)
        ta8 = integer_or_blank(card, 12, 'ta8', ta1)
        ta = [ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8]
        assert len(card) <= 13, 'len(CONV card) = %i\ncard=%s' % (len(card), card)
        return CONV(eid, pconid, film_node, cntrlnd, ta, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        #data_in = [eid, pconid, flmnd, cntrlnd,
                   #[ta1, ta2, ta3, ta5, ta6, ta7, ta8],
                   #[wt1, wt2, wt3, wt5, wt6, wt7, wt8]]
        ## weights are unique to MSC nastran
        eid, pconid, film_node, cntrlnd, ta, weights = data
        #ta1, ta2, ta3, ta5, ta6, ta7, ta8 = ta
        #wt1, wt2, wt3, wt5, wt6, wt7, wt8 = aft

        return CONV(eid, pconid, film_node, cntrlnd, ta, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CONV eid=%s' % self.eid
        ## TODO: eid???
        self.eid_ref = model.Element(self.eid, msg=msg)
        if model._xref == 1:  # True
            assert self.eid_ref.type in ['CHBDYG', 'CHBDYE', 'CHBDYP']

    def uncross_reference(self):
        del self.eid_ref

    def TA(self, i=None):
        if i is None:
            return self.ta
        return self.ta[i]

    def raw_fields(self):
        list_fields = ['CONV', self.eid, self.pconid, self.film_node,
                       self.cntrlnd] + self.ta
        return list_fields

    def repr_fields(self):
        film_node = set_blank_if_default(self.film_node, 0)
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)

        ta0 = self.ta[0]
        ta = [ta0]
        for tai in self.ta[1:]:
            ta.append(set_blank_if_default(tai, ta0))
        list_fields = ['CONV', self.eid, self.pconid, film_node, cntrlnd] + ta
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class CONVM(ThermalBC):
    """
    Specifies a forced convection boundary condition for heat transfer analysis
    through connection to a surface element (CHBDYi entry).
    """
    type = 'CONVM'

    def __init__(self, eid, pconvm, film_node, cntmdot, ta1, ta2, mdot,
                 comment=''):
        ThermalBC.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pconvm = pconvm
        self.film_node = film_node
        self.cntmdot = cntmdot
        self.ta1 = ta1
        self.ta2 = ta2
        self.mdot = mdot
        assert film_node >= 0
        if self.cntmdot == 0:
            if self.mdot is None:
                self.mdot = 1.0
        else:
            assert self.cntmdot > 0

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pconvm = integer(card, 2, 'pconvm')
        film_node = integer_or_blank(card, 3, 'film_node', 0)
        cntmdot = integer_or_blank(card, 4, 'cntmdot', 0)
        ta1 = integer(card, 5, 'ta1')
        ta2 = integer_or_blank(card, 6, 'ta2', ta1)
        mdot = double_or_blank(card, 7, 'mdot', 1.0)
        assert len(card) <= 8, 'len(CONVM card) = %i\ncard=%s' % (len(card), card)
        return CONVM(eid, pconvm, film_node, cntmdot, ta1, ta2, mdot,
                     comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        (eid, pconvm_id, film_node, cntrlnd, ta1, ta2, mdot) = data
        return CONVM(eid, pconvm_id, film_node, cntrlnd, ta1, ta2, mdot,
                     comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CONVM eid=%s' % self.eid
        self.eid = model.CYBDY(self.eid, msg=msg)
        self.pconvm = model.PCONV(self.pconvm, msg=msg)
        self.film_node = model.Grid(self.film_node, msg=msg)
        self.eid_ref = self.eid
        self.pconvm_ref = self.pconvm
        self.film_node_ref = self.film_node

    @property
    def film_node_id(self):
        if isinstance(self.film_node, integer_types):
            return self.film_node
        return self.film_node.nid

    def raw_fields(self):
        list_fields = ['CONVM', self.eid, self.pconvm, self.film_node_id,
                       self.cntmdot, self.ta1, self.ta2, self.mdot]
        return list_fields

    def repr_fields(self):
        film_node_id = set_blank_if_default(self.film_node_id, 0)
        ta2 = set_blank_if_default(self.ta2, self.ta1)
        mdot = set_blank_if_default(self.mdot, 1.0)
        list_fields = ['CONVM', self.eid, self.pconvm, film_node_id,
                       self.cntmdot, self.ta1, ta2, mdot]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RADM(ThermalBC):
    """
    Defines the radiation properties of a boundary element for heat transfer
    analysis
    """
    type = 'RADM'

    def __init__(self, radmid, absorb, emissivity, comment=''):
        ThermalBC.__init__(self)
        if comment:
            self._comment = comment

        #: Material identification number
        self.radmid = radmid

        self.absorb = absorb
        self.emissivity = emissivity

        assert self.radmid > 0, str(self)
        assert 0. <= self.absorb <= 1.0, str(self)
        for e in self.emissivity:
            assert 0. <= e <= 1.0, str(self)

    @classmethod
    def add_card(cls, card, comment=''):
        nfields = card.nfields
        radmid = integer(card, 1, 'radmid')
        absorb = double(card, 2, 'absorb')
        emissivity = fields(double, card, 'emissivity', i=3, j=nfields)
        return RADM(radmid, absorb, emissivity, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        radmid, absorb = data[:2]
        emissivity  = data[2:]
        return RADM(radmid, absorb, emissivity, comment=comment)

    #def cross_reference(self, model):
        #pass

    def repr_fields(self):
        list_fields = ['RADM', self.radmid, self.absorb] + self.emissivity
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RADBC(ThermalBC):
    """
    Specifies an CHBDYi element face for application of radiation boundary
    conditions
    """
    type = 'RADBC'

    def __init__(self, card=None, data=None, comment=''):
        ThermalBC.__init__(self)
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
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by RADBC pid=%s' % self.nodamb
        for i, eid in enumerate(self.eids):
            self.eids[i] = model.Element(eid, msg=msg)
        self.eids_ref = self.eids

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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

# Boundary Conditions
#-------------------------------------------------------
