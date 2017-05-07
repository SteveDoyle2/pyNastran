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


class ThermalBC(ThermalCard):
    def __init__(self):
        ThermalCard.__init__(self)


class ThermalElement(ThermalCard):
    pid = 0

    def __init__(self):
        ThermalCard.__init__(self)

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
    |   1    |  2  |   3  |  4   |   5    |    6   |    7    |    8    |
    +========+=====+======+======+========+========+=========+=========+
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

    def __init__(self, eid, eid2, side, iview_front=0, ivew_back=0,
                 rad_mid_front=0, rad_mid_back=0, comment=''):
        """
        Creates a CHBDYE card

        Parameters
        ----------
        eid : int
            surface element ID number for a side of an element
        eid2: int
            a heat conduction element identification
        side: int
            a consistent element side identification number (1-6)
        iview_front: int; default=0
            a VIEW entry identification number for the front face
        ivew_back: int; default=0
            a VIEW entry identification number for the back face
        rad_mid_front: int; default=0
            RADM identification number for front face of surface element
        rad_mid_back: int; default=0
            RADM identification number for back face of surface element
        comment : str; default=''
            a comment for the card
        """
        ThermalElement.__init__(self)
        if comment:
            self.comment = comment
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
        self.iview_front = iview_front

        #: A VIEW entry identification number for the back face
        self.ivew_back = ivew_back

        #: RADM identification number for front face of surface element
        #: (Integer > 0)
        self.rad_mid_front = rad_mid_front
        #: RADM identification number for back face of surface element
        #: (Integer > 0)
        self.rad_mid_back = rad_mid_back
        self.grids = []

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CHBDYE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        eid2 = integer(card, 2, 'eid2')
        side = integer(card, 3, 'side')

        iview_front = integer_or_blank(card, 4, 'iview_front', 0)
        ivew_back = integer_or_blank(card, 5, 'ivew_back', 0)
        rad_mid_front = integer_or_blank(card, 6, 'rad_mid_front', 0)
        rad_mid_back = integer_or_blank(card, 7, 'rad_mid_back', 0)
        assert len(card) <= 8, 'len(CHBDYE card) = %i\ncard=%s' % (len(card), card)
        return CHBDYE(eid, eid2, side, iview_front, ivew_back,
                      rad_mid_front, rad_mid_back, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CHBDYE card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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
        assert isinstance(eid, integer_types)
        assert isinstance(eid2, integer_types)
        assert isinstance(pid, integer_types)

    #def side_to_eids(self, eid):
        #side_ids = self.side_maps[eid.type][self.side]
        ## [1,2,3]

        ## id-1 is for the 0 based python index
        #nodes = [enodes[id - 1] for id in range(len(eid.nodes))
                 #if id in side_ids]
        #return side

    def Eid(self):
        return self.eid

    def Eid2(self):
        return self.eid2

    def raw_fields(self):
        list_fields = ['CHBDYE', self.eid, self.eid2, self.side,
                       self.iview_front, self.ivew_back, self.rad_mid_front,
                       self.rad_mid_back]
        return list_fields

    def repr_fields(self):
        """
        .. todo:: is this done
        """
        #eids = collapse_thru_by(self.eids)
        list_fields = ['CHBDYE', self.eid, self.eid2, self.side,
                       self.iview_front, self.ivew_back, self.rad_mid_front,
                       self.rad_mid_back]
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
    |    1   |  2  |  3 |   4  |    5   |    6   |    7    |    8    |  9  |
    +========+=====+====+======+========+========+=========+=========+=====+
    | CHBDYG | EID |    | TYPE | IVIEWF | IVIEWB | RADMIDF | RADMIDB |     |
    +--------+-----+----+------+--------+--------+---------+---------+-----+
    |        | G1  | G2 |  G3  |   G4   |   G5   |   G6    |   G7    |  G8 |
    +--------+-----+----+------+--------+--------+---------+---------+-----+
    """
    type = 'CHBDYG'

    def __init__(self, eid, Type, nodes, iview_front=0, ivew_back=0,
                 rad_mid_front=0, rad_mid_back=0, comment=''):
        ThermalElement.__init__(self)
        if comment:
            self.comment = comment
        #: Surface element ID
        self.eid = eid

        #: Surface type
        self.Type = Type

        #: A VIEW entry identification number for the front face
        self.iview_front = iview_front

        #: A VIEW entry identification number for the back face
        self.ivew_back = ivew_back

        #: RADM identification number for front face of surface element
        #: (Integer > 0)
        self.rad_mid_front = rad_mid_front

        #: RADM identification number for back face of surface element
        #: (Integer > 0)
        self.rad_mid_back = rad_mid_back

        #: Grid point IDs of grids bounding the surface (Integer > 0)
        self.nodes = nodes

        assert self.Type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], 'Type=%r' % Type
        assert len(nodes) > 0, nodes

    def validate(self):
        #assert self.Type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], 'Type=%r' % self.Type
        if self.Type != 'REV':
            nnodes = int(self.Type[-1])
            if len(self.nodes) != nnodes:
                msg = 'nnodes=%s Type=%r' % (len(self.nodes), self.Type)
                raise ValueError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CHBDYG card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        # no field 2

        Type = string(card, 3, 'Type')
        i_view_front = integer_or_blank(card, 4, 'iview_front', 0)
        i_view_back = integer_or_blank(card, 8, 'ivew_back', 0)
        rad_mid_front = integer_or_blank(card, 6, 'rad_mid_front', 0)
        rad_mid_back = integer_or_blank(card, 7, 'rad_mid_back', 0)
        # no field 8

        n = 1
        nodes = []
        for i in range(9, len(card)):
            grid = integer_or_blank(card, i, 'grid%i' % n)
            nodes.append(grid)  # used to have a None option
        assert len(nodes) > 0, 'card=%s' % card
        return CHBDYG(eid, Type, nodes,
                      iview_front=i_view_front, ivew_back=i_view_back,
                      rad_mid_front=rad_mid_front, rad_mid_back=rad_mid_back,
                      comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CHBDYG card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        Type = data[1]
        i_view_front = data[2]
        i_view_back = data[3]
        rad_mid_front = data[4]
        rad_mid_back = data[5]
        nodes = [datai for datai in data[6:14] if datai > 0]
        if Type == 3:
            Type = 'REV'
        elif Type == 4:
            Type = 'AREA3'
        elif Type == 5:
            Type = 'AREA4'
        #elif Type == 7: # ???
            #Type = 'AREA6'
        elif Type == 8:
            Type = 'AREA6'
        elif Type == 9:
            Type = 'AREA8'
        else:
            raise NotImplementedError('eid=%s Type=%r' % (eid, Type))

        assert Type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], 'Type=%r data=%s' % (Type, data)
        return CHBDYG(eid, Type, nodes,
                      iview_front=i_view_front, ivew_back=i_view_back,
                      rad_mid_front=rad_mid_front, rad_mid_back=rad_mid_back,
                      comment=comment)

    def _verify(self, xref=False):
        eid = self.Eid()
        assert isinstance(eid, integer_types)

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
            ['CHBDYG', self.eid, None, self.Type, self.iview_front,
             self.ivew_back, self.rad_mid_front, self.rad_mid_back, None,] +
            self.node_ids)
        return list_fields

    def repr_fields(self):
        i_view_front = set_blank_if_default(self.iview_front, 0)
        i_view_back = set_blank_if_default(self.ivew_back, 0)
        rad_mid_front = set_blank_if_default(self.rad_mid_front, 0)
        rad_mid_back = set_blank_if_default(self.rad_mid_back, 0)

        list_fields = (['CHBDYG', self.eid, None, self.Type, i_view_front,
                        i_view_back, rad_mid_front, rad_mid_back, None, ] + self.node_ids)
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
    |    1   |    2    |    3    |   4  |    5   |    6   |  7 |  8 |  9 |
    +========+=========+=========+======+========+========+====+====+====+
    | CHBDYP |   EID   |   PID   | TYPE | IVIEWF | IVIEWB | G1 | G2 | G0 |
    +--------+---------+---------+------+--------+--------+----+----+----+
    |        | RADMIDF | RADMIDB | GMID |   CE   |   E1   | E2 | E3 |    |
    +--------+---------+---------+------+--------+--------+----+----+----+
    """
    type = 'CHBDYP'

    def __init__(self, eid, pid, Type, g1, g2, g0=0, gmid=None, ce=0,
                 iview_front=0, ivew_back=0,
                 rad_mid_front=0, rad_mid_back=0,
                 e1=None, e2=None, e3=None, comment=''):
        """
        Creates a CHBDYP card

        Parameters
        ----------
        eid : int
            Surface element ID
        pid : int
            PHBDY property entry identification numbers. (Integer > 0)
        Type : str
            Surface type
            Must be {POINT, LINE, ELCYL, FTUBE, TUBE}
        iview_front : int; default=0
            A VIEW entry identification number for the front face.
        ivew_back : int; default=0
            A VIEW entry identification number for the back face.
        g1 / g2 : int
            Grid point identification numbers of grids bounding the surface
        g0 : int; default=0
            Orientation grid point
        rad_mid_front : int
            RADM identification number for front face of surface element
        rad_mid_back : int
            RADM identification number for back face of surface element.
        gmid : int
            Grid point identification number of a midside node if it is used
            with the line type surface element.
        ce : int; default=0
            Coordinate system for defining orientation vector
        e1 / e2 / e3 : float; default=None
            Components of the orientation vector in coordinate system CE.
            The origin of the orientation vector is grid point G1.
        comment : str; default=''
            a comment for the card
        """
        ThermalElement.__init__(self)
        if comment:
            self.comment = comment
        #: Surface element ID
        self.eid = eid

        #: PHBDY property entry identification numbers. (Integer > 0)
        self.pid = pid
        self.Type = Type

        #: A VIEW entry identification number for the front face.
        self.iview_front = iview_front

        #: A VIEW entry identification number for the back face.
        self.ivew_back = ivew_back

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
        self.rad_mid_front = rad_mid_front

        #: RADM identification number for back face of surface element.
        #: (Integer > 0)
        self.rad_mid_back = rad_mid_back

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
        """
        Adds a CHBDYP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')

        Type = string(card, 3, 'Type')
        #print("self.Type = %s" % self.Type)
        # msg = 'CHBDYP Type=%r' % self.Type
        #assert self.Type in ['POINT','LINE','ELCYL','FTUBE','TUBE'], msg

        iview_front = integer_or_blank(card, 4, 'iview_front', 0)
        ivew_back = integer_or_blank(card, 5, 'ivew_back', 0)
        g1 = integer(card, 6, 'g1')

        if Type != 'POINT':
            g2 = integer(card, 7, 'g2')
        else:
            g2 = blank(card, 7, 'g2')

        g0 = integer_or_blank(card, 8, 'g0', 0)
        rad_mid_front = integer_or_blank(card, 9, 'rad_mid_front', 0)
        rad_mid_back = integer_or_blank(card, 10, 'rad_mid_back', 0)
        gmid = integer_or_blank(card, 11, 'gmid')
        ce = integer_or_blank(card, 12, 'ce', 0)
        e1 = double_or_blank(card, 13, 'e3')
        e2 = double_or_blank(card, 14, 'e3')
        e3 = double_or_blank(card, 15, 'e3')
        assert len(card) <= 16, 'len(CHBDYP card) = %i\ncard=%s' % (len(card), card)
        return CHBDYP(eid, pid, Type, g1, g2, g0=g0, gmid=gmid, ce=ce,
                      iview_front=iview_front, ivew_back=ivew_back,
                      rad_mid_front=rad_mid_front, rad_mid_back=rad_mid_back,
                      e1=e1, e2=e2, e3=e3, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CHBDYP card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        [eid, pid, Type, iviewf, iviewb, g1, g2, g0, radmidf, radmidb,
         dislin, ce, e1, e2, e3] = data
        #eid = data[0]
        #Type = data[1]
        #iview_front = data[2]
        #ivew_back = data[3]
        #rad_mid_front = data[4]
        #rad_mid_back = data[5]
        #nodes = [datai for datai in data[6:14] if datai > 0]
        if Type == 1:
            Type = 'POINT'
        elif Type == 2:
            Type = 'LINE'
        elif Type == 6:
            Type = 'ELCYL'
        elif Type == 7:
            Type = 'FTUBE'
        elif Type == 10:
            Type = 'TUBE'
        else:
            raise NotImplementedError('CHBDYP Type=%r data=%s' % (Type, data))
        #assert Type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], 'Type=%r data=%s' % (Type, data)

        #assert dislin == 0, 'CHBDYP dislin=%r data=%s' % (dislin, data)
        if dislin  == 0:
            gmid = None
        else:
            gmid = dislin
        return CHBDYP(eid, pid, Type, g1, g2, g0=g0, gmid=gmid, ce=ce,
                      iview_front=iviewf, ivew_back=iviewb,
                      rad_mid_front=radmidf, rad_mid_back=radmidb,
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
        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)

    def Eid(self):
        return self.eid

    def raw_fields(self):
        (g1, g2, g0, gmid) = self.node_ids
        list_fields = ['CHBDYP', self.eid, self.Pid(), self.Type,
                       self.iview_front, self.ivew_back, g1, g2, g0,
                       self.rad_mid_front, self.rad_mid_back, gmid, self.ce,
                       self.e1, self.e2, self.e3]
        return list_fields

    def repr_fields(self):
        iview_front = set_blank_if_default(self.iview_front, 0)
        ivew_back = set_blank_if_default(self.ivew_back, 0)
        rad_mid_front = set_blank_if_default(self.rad_mid_front, 0)
        rad_mid_back = set_blank_if_default(self.rad_mid_back, 0)

        (g1, g2, g0, gmid) = self.node_ids
        g0 = set_blank_if_default(g0, 0)
        ce = set_blank_if_default(self.ce, 0)

        list_fields = ['CHBDYP', self.eid, self.Pid(), self.Type, iview_front,
                       ivew_back, g1, g2, g0, rad_mid_front, rad_mid_back,
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

    Format (MSC 2005.2)
    +-------+--------+-------+-------+-------+-------+-----+----+----+
    |   1   |    2   |   3   |   4   |   5   |   6   |  7  | 8  |  9 |
    +=======+========+=======+=======+=======+=======+=====+====+====+
    | PCONV | PCONID |  MID  | FORM  | EXPF  | FTYPE | TID |    |    |
    +-------+--------+-------+-------+-------+-------+-----+----+----+
    |       | CHLEN  | GIDIN |  CE   |  E1   |   E2  |  E3 |    |    |
    +-------+--------+-------+-------+-------+-------+-----+----+----+
    | PCONV |   38   |  21   |   2   |  54   |       |     |    |    |
    +-------+--------+-------+-------+-------+-------+-----+----+----+
    |       |   2.0  |  235  |   0   |  1.0  |  0.0  | 0.0 |    |    |
    +-------+--------+-------+-------+-------+-------+-----+----+----+

    Alternate format (MSC 2005.2):
    +-------+--------+-------+-------+-------+-------+-----+----+----+
    |   1   |    2   |   3   |   4   |   5   |   6   |  7  | 8  |  9 |
    +=======+========+=======+=======+=======+=======+=====+====+====+
    | PCONV | PCONID |  MID  | FORM  | EXPF  |   3   | H1  | H2 | H3 |
    +-------+--------+-------+-------+-------+-------+-----+----+----+
    |       |   H4   |  H5   |  H6   |  H7   |  H8   |     |    |    |
    +-------+--------+-------+-------+-------+-------+-----+----+----+
    | PCONV |   7    |   3   | 10.32 | 10.05 | 10.09 |     |    |    |
    +-------+--------+-------+-------+-------+-------+-----+----+----+
    |       | 10.37  |       |       |       |       |     |    |    |
    +-------+--------+-------+-------+-------+-------+-----+----+----+

    .. todo:: alternate format is not supported; NX not checked
    """
    type = 'PCONV'

    def __init__(self, pconid, mid, form=0, expf=0.0, ftype=0, tid=None,
                 chlen=None, gidin=None, ce=0,
                 e1=None, e2=None, e3=None, comment=''):
        """
        Creates a PCONV card

        Parameters
        ----------
        pconid : int
            Convection property ID
        mid : int
            Material ID
        form : int; default=0
            Type of formula used for free convection
            Must be {0, 1, 10, 11, 20, or 21}
        expf : float; default=0.0
            Free convection exponent as implemented within the context
            of the particular form that is chosen
        ftype : int; default=0
            Formula type for various configurations of free convection
        tid : int; default=None
            Identification number of a TABLEHT entry that specifies the
            two variable tabular function of the free convection heat
            transfer coefficient
        chlen : float; default=None
            Characteristic length
        gidin : int; default=None
            Grid ID of the referenced inlet point
        ce : int; default=0
            Coordinate system for defining orientation vector.
        e1 / e2 / e3 : List[float]; default=None
            Components of the orientation vector in coordinate system CE.
            The origin of the orientation vector is grid point G1
        comment : str; default=''
            a comment for the card
        """
        ThermalProperty.__init__(self)
        if comment:
            self.comment = comment
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
        """
        Adds a PCONV card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        return PCONV(pconid, mid,
                     form=form, expf=expf, ftype=ftype,
                     tid=tid, chlen=chlen, gidin=gidin,
                     ce=ce, e1=e1, e2=e2, e3=e3, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PCONV card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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

    +--------+--------+-----+------+------+-------+------+-------+-------+
    |    1   |    2   |  3  |   4  |   5  |   6   |   7  |   8   |   9   |
    +========+========+=====+======+======+=======+======+=======+=======+
    | PCONVM | PCONID | MID | FORM | FLAG | COEF  | EXPR | EXPPI | EXPPO |
    +--------+--------+-----+------+------+-------+------+-------+-------+
    | PCONVM |    3   |  2  |   1  |   1  | 0.023 | 0.80 | 0.40  | 0.30  |
    +--------+--------+-----+------+------+-------+------+-------+-------+
    """
    type = 'PCONVM'

    def __init__(self, pconid, mid, coef, form=0, flag=0,
                 expr=0.0, exppi=0.0, exppo=0.0, comment=''):
        """
        Creates a PCONVM card

        Parameters
        ----------
        pconid : int
            Convection property ID
        mid: int
            Material ID
        coef: float
            Constant coefficient used for forced convection
        form: int; default=0
            Type of formula used for free convection
            Must be {0, 1, 10, 11, 20, or 21}
        flag: int; default=0
            Flag for mass flow convection
        expr: float; default=0.0
            Reynolds number convection exponent
        exppi: float; default=0.0
            Prandtl number convection exponent for heat transfer into
            the working fluid
        exppo: float; default=0.0
            Prandtl number convection exponent for heat transfer out of
            the working fluid
        comment : str; default=''
            a comment for the card
        """
        ThermalProperty.__init__(self)
        if comment:
            self.comment = comment
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
        #: working fluid. (Real > 0.0; Default = 0.0)
        self.exppi = exppi

        #: Prandtl number convection exponent for heat transfer out of the
        #: working fluid. (Real > 0.0; Default = 0.0)
        self.exppo = exppo

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PCONVM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pconid = integer(card, 1, 'pconid')
        mid = integer(card, 2, 'mid')
        form = integer_or_blank(card, 3, 'form', 0)
        flag = integer_or_blank(card, 4, 'flag', 0)
        coef = double(card, 5, 'coef')
        expr = double_or_blank(card, 6, 'expr', 0.0)
        exppi = double_or_blank(card, 7, 'exppi', 0.0)
        exppo = double_or_blank(card, 8, 'exppo', 0.0)
        assert len(card) <= 9, 'len(PCONVM card) = %i\ncard=%s' % (len(card), card)
        return PCONVM(pconid, mid, coef, form=form, flag=flag,
                      expr=expr, exppi=exppi, exppo=exppo, comment=comment)

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

    +-------+-----+------+-----+-----+
    |   1   |  2  |   3  |  4  | 5   |
    +=======+=====+======+=====+=====+
    | PHBDY | PID |  AF  | D1  | D2  |
    +-------+-----+------+-----+-----+
    | PHBDY |  2  | 0.02 | 1.0 | 1.0 |
    +-------+-----+------+-----+-----+
    """
    type = 'PHBDY'

    def __init__(self, pid, af=None, d1=None, d2=None, comment=''):
        """
        Creates a PHBDY card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id
        af : int
            Area factor of the surface used only for CHBDYP element
            Must be {POINT, LINE, TUBE, ELCYL}
            TUBE : constant thickness of hollow tube
        d1, d2 : float; default=None
            Diameters associated with the surface
            Used with CHBDYP [ELCYL, TUBE, FTUBE] surface elements
        comment : str; default=''
            a comment for the card
        """
        ThermalProperty.__init__(self)
        if comment:
            self.comment = comment
        if d2 is None:
            d2 = d1

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
        """
        Adds a PHBDY card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        af = double_or_blank(card, 2, 'af')
        d1 = double_or_blank(card, 3, 'd1')
        d2 = double_or_blank(card, 4, 'd2', d1)
        assert len(card) <= 5, 'len(PHBDY card) = %i\ncard=%s' % (len(card), card)
        return PHBDY(pid, af, d1, d2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PHBDY card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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

    def __init__(self, eid, pconid, ta, film_node=0, cntrlnd=0, comment=''):
        """
        Creates a CONV card

        Parameters
        ----------
        eid : int
            element id
        pconid : int
            Convection property ID
        mid : int
            Material ID
        ta : List[int]
            Ambient points used for convection 0's are allowed for TA2
            and higher
        film_node : int; default=0
            Point for film convection fluid property temperature
        cntrlnd : int; default=0
            Control point for free convection boundary condition
        comment : str; default=''
            a comment for the card
        """
        ThermalBC.__init__(self)
        if comment:
            self.comment = comment

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
        if isinstance(ta, integer_types):
            self.ta = [ta]
        else:
            self.ta = ta
        assert self.eid > 0, 'eid=%s\n%s' % (eid, str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CONV card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        return CONV(eid, pconid, ta, film_node, cntrlnd, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CONV card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        #data_in = [eid, pconid, flmnd, cntrlnd,
                   #[ta1, ta2, ta3, ta5, ta6, ta7, ta8],
                   #[wt1, wt2, wt3, wt5, wt6, wt7, wt8]]
        ## weights are unique to MSC nastran
        eid, pconid, film_node, cntrlnd, ta, weights = data
        #ta1, ta2, ta3, ta5, ta6, ta7, ta8 = ta
        #wt1, wt2, wt3, wt5, wt6, wt7, wt8 = aft

        return CONV(eid, pconid, ta, film_node, cntrlnd, comment=comment)

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

    +-------+-----+--------+-------+---------+-----+-----+------+
    |   1   |  2  |    3   |   4   |    5    |  6  |  7  |   8  |
    +=======+=====+========+=======+=========+=====+=====+======+
    | CONVM | EID | PCONID | FLMND | CNTMDOT | TA1 | TA2 | Mdot |
    +-------+-----+--------+-------+---------+-----+-----+------+
      CONVM | 101 |    1   |  201  |   301   |  20 |  21 |      |
    +-------+-----+--------+-------+---------+-----+-----+------+
    """
    type = 'CONVM'

    def __init__(self, eid, pconvm, ta1, film_node=0, cntmdot=0,
                 ta2=None, mdot=1.0, comment=''):
        """
        Creates a CONVM card

        Parameters
        ----------
        eid : int
            element id (CHBDYP)
        pconid : int
            property ID (PCONVM)
        mid : int
            Material ID
        ta1 : int
            ambient point for convection
        ta2 : int; default=None
            None : ta1
            ambient point for convection
        film_node : int; default=0
        cntmdot : int; default=0
            control point used for controlling mass flow
            0/blank is only allowed when mdot > 0
        mdot : float; default=1.0
            a multiplier for the mass flow rate in case there is no
            point associated with the CNTRLND field
            required if cntmdot = 0
        comment : str; default=''
            a comment for the card
        """
        ThermalBC.__init__(self)
        if comment:
            self.comment = comment
        if ta2 is None:
            ta2 = ta1
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
        """
        Adds a CONVM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pconvm = integer(card, 2, 'pconvm')
        film_node = integer_or_blank(card, 3, 'film_node', 0)
        cntmdot = integer_or_blank(card, 4, 'cntmdot', 0)
        ta1 = integer(card, 5, 'ta1')
        ta2 = integer_or_blank(card, 6, 'ta2', ta1)
        mdot = double_or_blank(card, 7, 'mdot', 1.0)
        assert len(card) <= 8, 'len(CONVM card) = %i\ncard=%s' % (len(card), card)
        return CONVM(eid, pconvm, ta1, film_node=film_node, cntmdot=cntmdot,
                     ta2=ta2, mdot=mdot, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CONVM card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        (eid, pconvm_id, film_node, cntrlnd, ta1, ta2, mdot) = data
        return CONVM(eid, pconvm_id, ta1, film_node, cntrlnd, ta2, mdot,
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
            self.comment = comment

        #: Material identification number
        self.radmid = radmid

        self.absorb = absorb
        if isinstance(emissivity, float):
            self.emissivity = [emissivity]
        else:
            self.emissivity = emissivity

        assert self.radmid > 0, str(self)
        assert 0. <= self.absorb <= 1.0, str(self)
        for e in self.emissivity:
            assert 0. <= e <= 1.0, str(self)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RADM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        nfields = card.nfields
        radmid = integer(card, 1, 'radmid')
        absorb = double(card, 2, 'absorb')
        emissivity = fields(double, card, 'emissivity', i=3, j=nfields)
        return RADM(radmid, absorb, emissivity, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a RADM card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        radmid, absorb = data[:2]
        emissivity  = data[2:]
        return RADM(radmid, absorb, emissivity, comment=comment)

    #def cross_reference(self, model):
        #pass

    def raw_fields(self):
        list_fields = ['RADM', self.radmid, self.absorb] + self.emissivity
        return list_fields

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

    def __init__(self, nodamb, famb, cntrlnd, eids, comment=''):
        ThermalBC.__init__(self)
        if comment:
            self.comment = comment

        #: NODAMB Ambient point for radiation exchange. (Integer > 0)
        self.nodamb = nodamb

        #: Radiation view factor between the face and the ambient point.
        #: (Real > 0.0)
        self.famb = famb

        #: Control point for thermal flux load. (Integer > 0; Default = 0)
        self.cntrlnd = cntrlnd

        #: CHBDYi element identification number
        self.eids = expand_thru_by(eids)

        assert self.nodamb > 0
        assert self.famb > 0.0
        assert self.cntrlnd >= 0
        min_eid = min(self.eids)
        if min_eid < 1:
            msg = 'min(eids)=%i' % min_eid
            raise ValueError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RADBC card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        nodamb = integer(card, 1, 'nodamb')
        famb = double(card, 2, 'famb')
        cntrlnd = integer_or_blank(card, 3, 'cntrlnd', 0)

        nfields = card.nfields
        eids = fields(integer_or_string, card, 'eid', i=4, j=nfields)
        return RADBC(nodamb, famb, cntrlnd, eids, comment=comment)

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
