# pylint: disable=R0902,R0904,R0914,C0111
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import BaseCard, _node_ids
#from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank,
    string, blank, string_or_blank)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class ThermalCard(BaseCard):
    def __init__(self):
        BaseCard.__init__(self)

    def cross_reference(self, model: BDF) -> None:
        raise NotImplementedError('%s has not defined the cross_reference '
                                  'method' % self.__class__.__name__)


class ThermalBC(ThermalCard):
    def __init__(self):
        ThermalCard.__init__(self)


class ThermalElement(ThermalCard):

    def __init__(self):
        ThermalCard.__init__(self)

    #def Centroid(self):
        #return np.zeros(3)

    #def center_of_mass(self):
        #return self.Centroid()

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
    _properties = ['hex_map', 'pent_map', 'tet_map', 'side_maps']
    hex_map = {
        1: [4, 3, 2, 1],
        2: [1, 2, 6, 5],
        3: [2, 3, 7, 6],
        4: [3, 4, 8, 7],
        5: [4, 1, 5, 8],
        6: [5, 6, 7, 8],
    }

    pent_map = {
        1: [3, 2, 1],
        2: [1, 2, 5, 4],
        3: [2, 3, 6, 5],
        4: [3, 1, 4, 6],
        5: [4, 5, 6],
    }

    tet_map = {
        1: [1, 3, 2],
        2: [1, 2, 4],
        3: [2, 3, 4],
        4: [3, 1, 4],
    }

    side_maps = {
        'CHEXA': hex_map,
        'CPENTA': pent_map,
        'CTETRA': tet_map,
        'CTRIA3': [1, 2, 3],
        'CQUAD4': [1, 2, 3, 4],
    }

    #pid = 0

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        eid2 = 2
        side = 1
        return CHBDYE(eid, eid2, side, iview_front=0, iview_back=0,
                      rad_mid_front=0, rad_mid_back=0, comment='')

    def __init__(self, eid, eid2, side, iview_front=0, iview_back=0,
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
        iview_back: int; default=0
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
        self.iview_back = iview_back

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
        iview_back = integer_or_blank(card, 5, 'iview_back', 0)
        rad_mid_front = integer_or_blank(card, 6, 'rad_mid_front', 0)
        rad_mid_back = integer_or_blank(card, 7, 'rad_mid_back', 0)
        assert len(card) <= 8, 'len(CHBDYE card) = %i\ncard=%s' % (len(card), card)
        return CHBDYE(eid, eid2, side, iview_front, iview_back,
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

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, xref_errors):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    @property
    def nodes(self):
        return []

    @property
    def node_ids(self):
        return _node_ids(self, nodes=self.grids, allow_empty_nodes=False, msg='')

    def Pid(self):
        return 0
        #if self.pid_ref is not None:
            #return self.pid_ref.pid
        #return self.pid

    def get_edge_ids(self):
        # TODO: not implemented
        return []

    def _verify(self, xref):
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
                       self.iview_front, self.iview_back, self.rad_mid_front,
                       self.rad_mid_back]
        return list_fields

    def repr_fields(self):
        """
        .. todo:: is this done
        """
        #eids = collapse_thru_by(self.eids)
        list_fields = ['CHBDYE', self.eid, self.eid2, self.side,
                       self.iview_front, self.iview_back, self.rad_mid_front,
                       self.rad_mid_back]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
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
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        surface_type = 'AREA3'
        nodes = [1, 2]
        return CHBDYG(eid, surface_type, nodes, iview_front=0, iview_back=0,
                      rad_mid_front=0, rad_mid_back=0, comment='')

    def __init__(self, eid, surface_type, nodes, iview_front=0, iview_back=0,
                 rad_mid_front=0, rad_mid_back=0, comment=''):
        ThermalElement.__init__(self)
        if comment:
            self.comment = comment
        #: Surface element ID
        self.eid = eid

        #: Surface type
        self.surface_type = surface_type

        #: A VIEW entry identification number for the front face
        self.iview_front = iview_front

        #: A VIEW entry identification number for the back face
        self.iview_back = iview_back

        #: RADM identification number for front face of surface element
        #: (Integer > 0)
        self.rad_mid_front = rad_mid_front

        #: RADM identification number for back face of surface element
        #: (Integer > 0)
        self.rad_mid_back = rad_mid_back

        #: Grid point IDs of grids bounding the surface (Integer > 0)
        self.nodes = nodes

        assert self.surface_type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], 'surface_type=%r' % surface_type
        assert len(nodes) > 0, nodes
        self.nodes_ref = None

    def validate(self):
        #assert self.surface_type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], 'surface_type=%r' % self.surface_type
        if self.surface_type == 'REV':
            assert len(self.nodes) in [2, 3], 'CHBDYG: REV; nodes=%s' % str(self.nodes)
        elif self.surface_type == 'REV1':
            assert len(self.nodes) == 3, 'CHBDYG: REV; nodes=%s' % str(self.nodes)
        else:
            if self.surface_type == 'AREA3':
                nnodes_required = 3
                nnodes_allowed = 3
            elif self.surface_type == 'AREA4':
                nnodes_required = 4
                nnodes_allowed = 4
            elif self.surface_type == 'AREA6':
                nnodes_required = 3
                nnodes_allowed = 6
            elif self.surface_type == 'AREA8':
                nnodes_required = 4
                nnodes_allowed = 8
            else:
                raise RuntimeError('CHBDYG surface_type=%r' % self.surface_type)

            if len(self.nodes) < nnodes_required:
                msg = 'nnodes=%s nnodes_required=%s; surface_type=%r' % (
                    len(self.nodes), nnodes_required, self.surface_type)
                raise ValueError(msg)
            if len(self.nodes) > nnodes_allowed:
                msg = 'nnodes=%s nnodes_allowed=%s; surface_type=%r' % (
                    len(self.nodes), nnodes_allowed, self.surface_type)
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

        surface_type = string(card, 3, 'Type')
        iview_front = integer_or_blank(card, 4, 'iview_front', 0)
        iview_back = integer_or_blank(card, 8, 'iview_back', 0)
        rad_mid_front = integer_or_blank(card, 6, 'rad_mid_front', 0)
        rad_mid_back = integer_or_blank(card, 7, 'rad_mid_back', 0)
        # no field 8

        n = 1
        nodes = []
        for i in range(9, len(card)):
            grid = integer_or_blank(card, i, 'grid%i' % n)
            nodes.append(grid)  # used to have a None option
        assert len(nodes) > 0, 'card=%s' % card
        return CHBDYG(eid, surface_type, nodes,
                      iview_front=iview_front, iview_back=iview_back,
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
        surface_type = data[1]
        i_view_front = data[2]
        i_view_back = data[3]
        rad_mid_front = data[4]
        rad_mid_back = data[5]
        nodes = [datai for datai in data[6:14] if datai > 0]
        surface_type_int_to_str = {
            3 : 'REV',
            4 : 'AREA3',
            5 : 'AREA4',
            #7: 'AREA6',# ???
            8 : 'AREA6',
            9 : 'AREA8',
        }
        try:
            surface_type = surface_type_int_to_str[surface_type]
        except KeyError:
            raise NotImplementedError('eid=%s surface_type=%r' % (eid, surface_type))

        assert surface_type in ['REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], 'surface_type=%r data=%s' % (surface_type, data)
        return CHBDYG(eid, surface_type, nodes,
                      iview_front=i_view_front, iview_back=i_view_back,
                      rad_mid_front=rad_mid_front, rad_mid_back=rad_mid_back,
                      comment=comment)

    def _verify(self, xref):
        eid = self.Eid()
        assert isinstance(eid, integer_types)

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        # TODO: is this correct?
        return _node_ids(self, nodes=self.nodes_ref, allow_empty_nodes=True, msg='')

    def get_edge_ids(self):
        # TODO: not implemented
        return []

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CHBDYG eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by CHBDYG eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def Eid(self):
        return self.eid

    def raw_fields(self):
        list_fields = (
            ['CHBDYG', self.eid, None, self.surface_type, self.iview_front,
             self.iview_back, self.rad_mid_front, self.rad_mid_back, None,] +
            self.node_ids)
        return list_fields

    def repr_fields(self):
        i_view_front = set_blank_if_default(self.iview_front, 0)
        i_view_back = set_blank_if_default(self.iview_back, 0)
        rad_mid_front = set_blank_if_default(self.rad_mid_front, 0)
        rad_mid_back = set_blank_if_default(self.rad_mid_back, 0)

        list_fields = (['CHBDYG', self.eid, None, self.surface_type, i_view_front,
                        i_view_back, rad_mid_front, rad_mid_back, None, ] + self.node_ids)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
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
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 2
        surface_type = 'POINT'
        g1 = 1
        g2 = None
        return CHBDYP(eid, pid, surface_type, g1, g2, g0=0, gmid=None, ce=0,
                      iview_front=0, iview_back=0,
                      rad_mid_front=0, rad_mid_back=0,
                      e1=None, e2=None, e3=None, comment='')

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        if isinstance(self.nodes, np.ndarray):
            self.nodes = self.nodes.tolist()
        self.nodes = [None if (nid is None or np.isnan(nid)) else nid
                      for nid in self.nodes]

    def __init__(self, eid, pid, surface_type, g1, g2, g0=0, gmid=None, ce=0,
                 iview_front=0, iview_back=0,
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
        surface_type : str
            Surface type
            Must be {POINT, LINE, ELCYL, FTUBE, TUBE}
        iview_front : int; default=0
            A VIEW entry identification number for the front face.
        iview_back : int; default=0
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
        self.surface_type = surface_type

        #: A VIEW entry identification number for the front face.
        self.iview_front = iview_front

        #: A VIEW entry identification number for the back face.
        self.iview_back = iview_back

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
        self.nodes_ref = None
        self.pid_ref = None
        self.ce_ref = None
        assert surface_type in ['POINT', 'LINE', 'ELCYL', 'FTUBE', 'TUBE'], surface_type

    @property
    def Type(self):
        return self.surface_type

    @Type.setter
    def Type(self, surface_type):
        self.surface_type = surface_type

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
        surface_type = string(card, 3, 'Type')

        iview_front = integer_or_blank(card, 4, 'iview_front', 0)
        iview_back = integer_or_blank(card, 5, 'iview_back', 0)
        g1 = integer(card, 6, 'g1')

        if surface_type != 'POINT':
            g2 = integer(card, 7, 'g2')
        else:
            g2 = blank(card, 7, 'g2')

        g0 = integer_or_blank(card, 8, 'g0', 0)
        rad_mid_front = integer_or_blank(card, 9, 'rad_mid_front', 0)
        rad_mid_back = integer_or_blank(card, 10, 'rad_mid_back', 0)
        gmid = integer_or_blank(card, 11, 'gmid')
        ce = integer_or_blank(card, 12, 'ce', 0)
        e1 = double_or_blank(card, 13, 'e1')
        e2 = double_or_blank(card, 14, 'e2')
        e3 = double_or_blank(card, 15, 'e3')
        assert len(card) <= 16, 'len(CHBDYP card) = %i\ncard=%s' % (len(card), card)
        return CHBDYP(eid, pid, surface_type, g1, g2, g0=g0, gmid=gmid, ce=ce,
                      iview_front=iview_front, iview_back=iview_back,
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
        [eid, pid, surface_type, iview_front, iview_back, g1, g2, g0, radmidf, radmidb,
         dislin, ce, e1, e2, e3] = data
        #eid = data[0]
        #surface_type = data[1]
        #iview_front = data[2]
        #iview_back = data[3]
        #rad_mid_front = data[4]
        #rad_mid_back = data[5]
        #nodes = [datai for datai in data[6:14] if datai > 0]
        surface_type_int_to_str = {
            1 : 'POINT',
            2 : 'LINE',
            6 : 'ELCYL',
            7 : 'FTUBE',
            10 : 'TUBE',
        }
        try:
            surface_type = surface_type_int_to_str[surface_type]
        except KeyError:  # pragma: no cover
            raise NotImplementedError('CHBDYP surface_type=%r data=%s' % (surface_type, data))

        #assert dislin == 0, 'CHBDYP dislin=%r data=%s' % (dislin, data)
        if dislin == 0:
            gmid = None
        else:
            gmid = dislin
        return CHBDYP(eid, pid, surface_type, g1, g2, g0=g0, gmid=gmid, ce=ce,
                      iview_front=iview_front, iview_back=iview_back,
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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        try:
            msg = ', which is required by CHBDYP pid=%s' % self.pid
            self.pid_ref = model.Phbdy(self.pid, msg=msg)
            self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
            self.ce_ref = model.Coord(self.ce, msg)
        except KeyError:
            print(self.get_stats())
            raise

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by CHBDYP pid=%s' % self.pid
        self.pid_ref = model.Phbdy(self.pid, msg=msg)
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.ce_ref = model.safe_coord(self.ce, self.pid, xref_errors, msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.ce = self.Ce()
        self.nodes_ref = None
        self.pid_ref = None
        self.ce_ref = None

    def _verify(self, xref):
        eid = self.Eid()
        pid = self.Pid()
        ce = self.Ce()
        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        assert isinstance(ce, integer_types)

    def Eid(self):
        return self.eid

    def Pid(self):
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    def Ce(self):
        """gets the coordinate system, CE"""
        if self.ce_ref is not None:
            return self.ce_ref.cid
        return self.ce

    def raw_fields(self):
        (g1, g2, g0, gmid) = self.node_ids
        list_fields = ['CHBDYP', self.eid, self.Pid(), self.surface_type,
                       self.iview_front, self.iview_back, g1, g2, g0,
                       self.rad_mid_front, self.rad_mid_back, gmid, self.Ce(),
                       self.e1, self.e2, self.e3]
        return list_fields

    def repr_fields(self):
        iview_front = set_blank_if_default(self.iview_front, 0)
        iview_back = set_blank_if_default(self.iview_back, 0)
        rad_mid_front = set_blank_if_default(self.rad_mid_front, 0)
        rad_mid_back = set_blank_if_default(self.rad_mid_back, 0)

        (g1, g2, g0, gmid) = self.node_ids
        g0 = set_blank_if_default(g0, 0)
        ce = set_blank_if_default(self.Ce(), 0)

        list_fields = ['CHBDYP', self.eid, self.Pid(), self.surface_type, iview_front,
                       iview_back, g1, g2, g0, rad_mid_front, rad_mid_back,
                       gmid, ce, self.e1, self.e2, self.e3]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
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

    @classmethod
    def _init_from_empty(cls):
        pconid = 1
        #mid = 1
        return PCONV(pconid, mid=None, form=0, expf=0.0, ftype=0, tid=None,
                     chlen=None, gidin=None, ce=0, e1=None, e2=None, e3=None, comment='')

    def __init__(self, pconid, mid=None, form=0, expf=0.0, ftype=0, tid=None,
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
        assert mid is None or self.mid > 0
        assert self.form in [0, 1, 10, 11, 20, 21]
        self.ce_ref = None
        self.gidin_ref = None

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
        mid = integer_or_blank(card, 2, 'mid')
        form = integer_or_blank(card, 3, 'form', 0)
        expf = double_or_blank(card, 4, 'expf', 0.0)
        ftype = integer_or_blank(card, 5, 'ftype', 0)
        tid = integer_or_blank(card, 6, 'tid')
        chlen = double_or_blank(card, 9, 'chlen')
        gidin = integer_or_blank(card, 10, 'gidin')
        ce = integer_or_blank(card, 11, 'ce', 0)
        e1 = double_or_blank(card, 12, 'e1')
        e2 = double_or_blank(card, 13, 'e2')
        e3 = double_or_blank(card, 14, 'e3')
        assert len(card) <= 15, 'len(PCONV card) = %i\ncard=%s' % (len(card), card)
        return PCONV(pconid, mid=mid,
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

    def Ce(self):
        """gets the coordinate system, CE"""
        if self.ce_ref is not None:
            return self.ce_ref.cid
        return self.ce

    def Gidin(self):
        """gets the grid input node, gidin"""
        if self.gidin_ref is not None:
            return self.gidin_ref.nid
        return self.gidin

    def cross_reference(self, model: BDF) -> None:
        msg = 'which is required by PCONV pconid=%s' % self.pconid
        self.ce_ref = model.Coord(self.ce, msg)
        if self.gidin is not None:
            self.gidin_ref = model.Node(self.gidin, msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.ce = self.Ce()
        self.ce_ref = None
        self.gidin = self.Gidin()
        self.gidin_ref = None

    def raw_fields(self):
        list_fields = ['PCONV', self.pconid, self.mid, self.form, self.expf,
                       self.ftype, self.tid, None, None, self.chlen, self.Gidin(),
                       self.Ce(), self.e1, self.e2, self.e3]
        return list_fields

    def repr_fields(self):
        form = set_blank_if_default(self.form, 0)
        expf = set_blank_if_default(self.expf, 0.0)
        ftype = set_blank_if_default(self.ftype, 0)
        ce = set_blank_if_default(self.Ce(), 0)
        list_fields = ['PCONV', self.pconid, self.mid, form, expf, ftype, self.tid,
                       None, None, self.chlen, self.Gidin(), ce, self.e1, self.e2,
                       self.e3]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
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

    @classmethod
    def _init_from_empty(cls):
        pconid = 1
        mid = 1
        coeff = 0.1
        return PCONVM(pconid, mid, coeff, form=0, flag=0,
                      expr=0.0, exppi=0.0, exppo=0.0, comment='')

    def __init__(self, pconid, mid, coeff, form=0, flag=0,
                 expr=0.0, exppi=0.0, exppo=0.0, comment=''):
        """
        Creates a PCONVM card

        Parameters
        ----------
        pconid : int
            Convection property ID
        mid: int
            Material ID
        coeff: float
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
        self.coef = coeff

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

    #def cross_reference(self, model: BDF) -> None:
        #pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
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

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        return PHBDY(pid, af=None, d1=None, d2=None, comment='')

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

    #def cross_reference(self, model: BDF) -> None:
        #pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        list_fields = ['PHBDY', self.pid, self.af, self.d1, self.d2]
        return list_fields

    def repr_fields(self):
        d2 = set_blank_if_default(self.d2, self.d1)
        list_fields = ['PHBDY', self.pid, self.af, self.d1, d2]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
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

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        #nodamb = 1
        pconid = 2
        ta = 1.0
        return CONV(eid, pconid, ta, film_node=0, cntrlnd=0, comment='')

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
        self.eid_ref = None

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
        del weights
        #ta1, ta2, ta3, ta5, ta6, ta7, ta8 = ta
        #wt1, wt2, wt3, wt5, wt6, wt7, wt8 = aft

        return CONV(eid, pconid, ta, film_node, cntrlnd, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CONV eid=%s' % self.eid
        ## TODO: eid???
        self.eid_ref = model.Element(self.eid, msg=msg)
        if model._xref == 1:  # True
            assert self.eid_ref.type in ['CHBDYG', 'CHBDYE', 'CHBDYP']

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.eid_ref = None

    def Eid(self):
        if self.eid_ref is not None:
            return self.eid_ref.eid
        return self.eid

    def TA(self, i=None):
        if i is None:
            return self.ta
        return self.ta[i]

    def raw_fields(self):
        list_fields = ['CONV', self.Eid(), self.pconid, self.film_node,
                       self.cntrlnd] + self.ta
        return list_fields

    def repr_fields(self):
        film_node = set_blank_if_default(self.film_node, 0)
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)

        ta0 = self.ta[0]
        ta = [ta0]
        for tai in self.ta[1:]:
            ta.append(set_blank_if_default(tai, ta0))
        list_fields = ['CONV', self.Eid(), self.pconid, film_node, cntrlnd] + ta
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class TEMPBC(ThermalBC):
    type = 'TEMPBC'
    _properties = ['eid']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        Type = 1
        nodes = [1, 2]
        temps = [10., 20.]
        return TEMPBC(sid, Type, nodes, temps, comment='')

    def __init__(self, sid, Type, nodes, temps, comment=''):
        ThermalBC.__init__(self)
        if comment:
            self.comment = comment

        self.sid = sid
        self.Type = Type
        self.nodes = nodes
        self.temps = temps

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TEMPBC card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        Type = string_or_blank(card, 2, 'Type', 'STAT')
        nfields_left = len(card) - 3
        assert nfields_left > 0, card
        assert nfields_left % 2 == 0, card

        temps = []
        nodes = []
        for i in range(nfields_left // 2):
            ifield = 3 + i*2
            temp = double(card, ifield, 'temp_%i'%  ((i+1)))
            nid = integer(card, ifield+1, 'temp_%i'%  ((i+1)))
            temps.append(temp)
            nodes.append(nid)
        return TEMPBC(sid, Type, nodes, temps, comment=comment)

    @property
    def eid(self):
        return self.sid

    def raw_fields(self):
        list_fields = ['TEMPBC', self.sid, self.Type]
        for temp, node in zip(self.temps, self.nodes):
            list_fields.extend([temp, node])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
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
    | CONVM | 101 |    1   |  201  |   301   |  20 |  21 |      |
    +-------+-----+--------+-------+---------+-----+-----+------+

    """
    type = 'CONVM'
    _properties = ['film_node_id', 'pconvm_id']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pconvm = 2
        ta1 = 1.0
        return CONVM(eid, pconvm, ta1, film_node=0, cntmdot=0,
                     ta2=None, mdot=1.0, comment='')

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
        self.eid_ref = None
        self.pconvm_ref = None
        self.film_node_ref = None

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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CONVM eid=%s' % self.eid
        self.eid_ref = model.CYBDY(self.eid, msg=msg)
        self.pconvm_ref = model.PCONV(self.pconvm, msg=msg)
        self.film_node_ref = model.Grid(self.film_node, msg=msg)

    def Eid(self):
        if self.eid_ref is not None:
            return self.eid_ref.eid
        return self.eid

    @property
    def film_node_id(self):
        if self.film_node_ref is not None:
            return self.film_node_ref.nid
        return self.film_node

    @property
    def pconvm_id(self):
        if self.pconvm_ref is not None:
            return self.pconvm_ref.pconid
        return self.pconvm

    def raw_fields(self):
        list_fields = ['CONVM', self.Eid(), self.pconvm_id, self.film_node_id,
                       self.cntmdot, self.ta1, self.ta2, self.mdot]
        return list_fields

    def repr_fields(self):
        film_node_id = set_blank_if_default(self.film_node_id, 0)
        ta2 = set_blank_if_default(self.ta2, self.ta1)
        mdot = set_blank_if_default(self.mdot, 1.0)
        list_fields = ['CONVM', self.Eid(), self.pconvm_id, film_node_id,
                       self.cntmdot, self.ta1, ta2, mdot]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

# Boundary Conditions
#-------------------------------------------------------
