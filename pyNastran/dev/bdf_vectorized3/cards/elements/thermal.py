from __future__ import annotations
from itertools import zip_longest
from typing import TYPE_CHECKING
import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8  # , print_float_8, print_field_8
# from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
# from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type_force import force_double, force_double_or_blank
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double,  # string,
    blank,
    integer_or_blank, double_or_blank, string_or_blank,
)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check, find_missing
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.cards.elements.bars import set_blank_if_default
#from pyNastran.bdf.cards.properties.bars import _bar_areaL # PBARL as pbarl, A_I1_I2_I12

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, parse_check,
    make_idim, hslice_by_idim, save_ifile_comment,
    remove_unused_primary, remove_unused_duplicate)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_int, array_default_float,
    array_float, get_print_card_size)

from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    # from pyNastran.dev.bdf_vectorized3.bdf import BDF


class BDYOR:
    """

    +--------+------+--------+--------+---------+---------+----+-----+-----+
    |    1   |   2  |    3   |    4   |    5    |    6    | 7  |  8  |  9  |
    +--------+------+--------+--------+---------+---------+----+-----+-----+
    | BDYOR  | TYPE | IVIEWF | IVIEWB | RADMIDF | RADMIDB |    | PID |  GO |
    +--------+------+--------+--------+---------+---------+----+-----+-----+
    |        |  CE  |   E1   |   E2   |    E3   |         |    |     |     |
    +--------+------+--------+--------+---------+---------+----+-----+-----+

    """
    type = 'BDYOR'

    #: allows the get_field method and update_field methods to be used
    #_field_map = {2:'cp', 6:'cd', 7:'ps', 8:'seid'}

    @classmethod
    def _init_from_empty(cls):
        surface_type = ''
        iview_front = 0
        iview_back = 0
        rad_mid_front = 0
        rad_mid_back = 0
        pid = -1
        g0 = 0
        ce = 0
        e123 = [0., 0., 0.]
        return BDYOR(surface_type,
                     iview_front, iview_back,
                     rad_mid_front, rad_mid_back,
                     pid, g0,
                     ce, e123, comment='')

    def __init__(self, surface_type: str,
                 iview_front: int, iview_back: int,
                 rad_mid_front: int, rad_mid_back: int,
                 pid: int, g0: int,
                 ce: int, e123: list[float],
                 ifile: int=0,
                 comment: str=''):
        """
        Creates the BDYOR card

        Parameters
        ----------
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        self.ifile = ifile
        self.surface_type = surface_type
        self.iview_front = iview_front
        self.iview_back = iview_back
        self.rad_mid_front = rad_mid_front
        self.rad_mid_back = rad_mid_back
        self.pid = pid
        self.g0 = g0
        self.ce = ce
        self.e123 = e123

    @classmethod
    def add_card(cls, card: BDFCard, ifile: int, comment: str=''):
        """
        Adds a BYDOR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        surface_type = string_or_blank(card, 1, 'chbdy_type', default='')
        iview_front = integer_or_blank(card, 2, 'iview_front', default=0)
        iview_back = integer_or_blank(card, 3, 'iview_back', default=0)
        rad_mid_front = integer_or_blank(card, 4, 'rad_mid_front', default=0)
        rad_mid_back = integer_or_blank(card, 5, 'rad_mid_back', default=0)
        #
        pid = integer_or_blank(card, 7, 'pid', default=-1)
        g0 = integer_or_blank(card, 8, 'g0', default=-1)
        ce = integer_or_blank(card, 9, 'ce', default=0)
        e123 = [
            double_or_blank(card, 10, 'e1', default=np.nan),
            double_or_blank(card, 11, 'e2', default=np.nan),
            double_or_blank(card, 12, 'e3', default=np.nan),
        ]
        assert len(card) <= 13, f'len(BDYOR card) = {len(card):d}\ncard={card}'
        return BDYOR(surface_type, iview_front, iview_back,
                     rad_mid_front, rad_mid_back,
                     pid, g0,
                     ce, e123, ifile=ifile, comment=comment)

    def __deepcopy__(self, memo_dict):
        raw_fields = self.raw_fields()
        card = BDFCard(raw_fields)
        return self.add_card(card, self.ifile, comment=self.comment)

    def raw_fields(self):
        g0 = '' if self.g0 == -1 else self.g0
        iview_front = '' if self.iview_front == 0 else self.iview_front
        iview_back = '' if self.iview_back == 0 else self.iview_back
        rad_mid_front = '' if self.rad_mid_front == 0 else self.rad_mid_front
        rad_mid_back = '' if self.rad_mid_back == 0 else self.rad_mid_back
        ce = '' if self.ce == 0 else self.ce

        list_fields = ['BDYOR', self.surface_type,
                       iview_front, iview_back,
                       rad_mid_front, rad_mid_back, None,
                       self.pid, g0,
                       ce] + self.e123
        return list_fields

    def write(self, size: int=8):
        list_fields = self.raw_fields()
        assert len(list_fields) > 0
        return print_card_8(list_fields)


class ThermalElement(VectorizedBaseCard):
    pass
    # def __init__(self, model: BDF):
    #     super().__init__(model)
    #     self.element_id = np.array([], dtype='int32')


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
    _id_name = 'element_id'

    def clear(self) -> None:
        self.element_id = np.zeros((0, 2), dtype='int32')
        self.n = 0

    def add(self, eid: int, eid2: int, side: int,
            iview_front: int=0, iview_back: int=0,
            rad_mid_front: int=0, rad_mid_back: int=0,
            comment: str='') -> int:
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
        self.cards.append((eid, eid2, side, iview_front, iview_back,
                           rad_mid_front, rad_mid_back, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        eid2 = integer(card, 2, 'eid2')
        side = integer(card, 3, 'side')

        iview_front = integer_or_blank(card, 4, 'iview_front', default=0)
        iview_back = integer_or_blank(card, 5, 'iview_back', default=0)
        rad_mid_front = integer_or_blank(card, 6, 'rad_mid_front', default=0)
        rad_mid_back = integer_or_blank(card, 7, 'rad_mid_back', default=0)
        assert len(card) <= 8, f'len(CHBDYE card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, eid2, side, iview_front, iview_back,
                           rad_mid_front, rad_mid_back, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self):
        ncards = len(self.cards)

        #: eid1: Surface element ID number for a side of an element. (0 < Integer < 100,000,000)
        #: eid2: A heat conduction element identification
        element_id = np.zeros((ncards, 2), dtype='int32')

        #: A consistent element side identification number
        #: (1 < Integer < 6)
        side = np.zeros(ncards, dtype='int32')

        #: A VIEW entry identification number for the front face
        #: A VIEW entry identification number for the back face
        iview = np.zeros((ncards, 2), dtype='int32')

        #: RADM identification number for front face of surface element (Integer > 0)
        #: RADM identification number for back face of surface element (Integer > 0)
        rad_mid = np.zeros((ncards, 2), dtype='int32')

        #grids = []
        for icard, card in enumerate(self.cards):
            (eid, eid2, sidei, iview_front, iview_back,
             rad_mid_front, rad_mid_back, comment) = card
            element_id[icard, :] = [eid, eid2]

            side[icard] = sidei
            iview[icard] = [iview_front, iview_back]
            rad_mid[icard, :] = [rad_mid_front, rad_mid_back]
            assert 0 < sidei < 7, sidei
        self._save(element_id, side, iview, rad_mid)
        apply_bydor_default(self)
        self.cards = []

    def _save(self, element_id, side, iview, rad_mid):
        assert len(self.element_id) == 0
        self.element_id = element_id
        self.side = side
        self.iview = iview
        self.rad_mid = rad_mid
        self.n = len(element_id)

    def __apply_slice__(self, load: CHBDYE, i: np.ndarray) -> None:
        load.element_id = self.element_id[i]
        load.side = self.side[i]
        load.iview = self.iview[i, :]
        load.rad_mid = self.rad_mid[i, :]
        load.n = len(i)

    @property
    def iview_front(self) -> np.ndarray:
        return self.iview[:, 0]

    @property
    def iview_back(self) -> np.ndarray:
        return self.iview[:, 1]

    @property
    def rad_mid_front(self) -> np.ndarray:
        return self.rad_mid[:, 0]

    @property
    def rad_mid_back(self) -> np.ndarray:
        return self.rad_mid[:, 1]

    @iview_front.setter
    def iview_front(self, iview_front: np.ndarray) -> None:
        self.iview[:, 0] = iview_front

    @iview_back.setter
    def iview_back(self, iview_back: np.ndarray) -> None:
        self.iview[:, 1] = iview_back

    @rad_mid_front.setter
    def rad_mid_front(self, rad_mid_front: np.ndarray) -> None:
        self.rad_mid[:, 0] = rad_mid_front

    @rad_mid_back.setter
    def rad_mid_back(self, rad_mid_back: np.ndarray) -> None:
        self.rad_mid[:, 1] = rad_mid_back

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.iview.max(), self.rad_mid.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        # TODO: default for iview, rad_mid is specified on BDYOR
        element_ids = self.element_id[:, 0]
        element_id2s = self.element_id[:, 1]
        sides = array_str(self.side, size=size)
        iview_fronts   = array_default_int(self.iview[:, 0],   default=0, size=size)
        iview_backs    = array_default_int(self.iview[:, 1],   default=0, size=size)
        rad_mid_fronts = array_default_int(self.rad_mid[:, 0], default=0, size=size)
        rad_mid_backs  = array_default_int(self.rad_mid[:, 1], default=0, size=size)
        for eid, eid2, side, iview_front, iview_back, rad_mid_front, rad_mid_back in zip(
                element_ids, element_id2s, sides,
                iview_fronts, iview_backs,
                rad_mid_fronts, rad_mid_backs):
            assert eid != 0, eid
            assert eid2 != 0, eid2
            #assert iview_front != 0, iview_front
            #assert iview_back != 0, iview_back
            #assert rad_mid_front != 0, rad_mid_front
            #assert rad_mid_back != 0, rad_mid_back
            list_fields = ['CHBDYE', eid, eid2, side,
                           iview_front, iview_back, rad_mid_front,
                           rad_mid_back]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.prod]
                if prop.n > 0]


class CONV(VectorizedBaseCard):
    """
    Specifies a free convection boundary condition for heat transfer analysis
    through connection to a surface element (CHBDYi entry).

    """
    _id_name = 'element_id'

    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.n = 0

    def add(self, eid: int, pconid: int,
            ta: list[int], film_node: int=0, control_node: int=0,
            comment: str='') -> int:
        """
        Creates a CONV card

        Parameters
        ----------
        eid : int
            element id
        pconid : int
            Convection property ID
        ta : list[int]
            Ambient points used for convection 0's are allowed for TA2
            and higher
        film_node : int; default=0
            Point for film convection fluid property temperature
        control_node : int; default=0
            Control point for free convection boundary condition
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pconid, film_node, control_node, ta, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pconid = integer(card, 2, 'pconid')
        film_node = integer_or_blank(card, 3, 'film_node', default=0)
        control_node = integer_or_blank(card, 4, 'control_node', default=0)

        ta1 = integer(card, 5, 'TA1')
        assert ta1 > 0, ta1

        ta2 = integer_or_blank(card, 6, 'ta2', default=ta1)
        ta3 = integer_or_blank(card, 7, 'ta3', default=ta1)
        ta4 = integer_or_blank(card, 8, 'ta4', default=ta1)
        ta5 = integer_or_blank(card, 9, 'ta5', default=ta1)
        ta6 = integer_or_blank(card, 10, 'ta6', default=ta1)
        ta7 = integer_or_blank(card, 11, 'ta7', default=ta1)
        ta8 = integer_or_blank(card, 12, 'ta8', default=ta1)
        ta = [ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8]
        assert len(card) <= 13, f'len(CONV card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pconid, film_node, control_node, ta, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype

        #: CHBDYG, CHBDYE, or CHBDYP surface element identification number.
        #: (Integer > 0)
        element_id = np.zeros(ncards, dtype=idtype)

        #: Convection property identification number of a PCONV entry
        pconv_id = np.zeros(ncards, dtype='int32')

        #: Point for film convection fluid property temperature
        film_node = np.zeros(ncards, dtype='int32')

        #: Control point for free convection boundary condition.
        control_node = np.zeros(ncards, dtype='int32')

        #: Ambient points used for convection 0's are allowed for TA2 and
        #: higher.  (Integer > 0 for TA1 and Integer > 0 for TA2 through TA8;
        #: Default for TA2 through TA8 is TA1.)
        temp_ambient = np.zeros((ncards, 8), dtype=idtype)

        #grids = []
        for icard, card in enumerate(self.cards):
            (eid, pconid, film_nodei, control_nodei, ta, comment) = card
            element_id[icard] = eid
            pconv_id[icard] = pconid
            film_node[icard] = film_nodei
            control_node[icard] = control_nodei
            temp_ambient[icard, :] = ta
        self._save(element_id, pconv_id, film_node, control_node, temp_ambient)
        self.sort()
        self.cards = []

    def _save(self, element_id, pconv_id, film_node, control_node, temp_ambient):
        nelements = len(element_id)
        assert element_id.min() > 0, element_id
        assert pconv_id.min() > 0, pconv_id
        assert film_node.min() >= 0, film_node
        assert control_node.min() >= 0, control_node
        self.element_id = element_id
        self.pconv_id = pconv_id
        self.film_node = film_node
        self.control_node = control_node
        self.temp_ambient = temp_ambient
        self.n = nelements

    def __apply_slice__(self, load: CONV, i: np.ndarray) -> None:
        load.element_id = self.element_id[i]
        load.pconv_id = self.pconv_id[i]
        load.film_node = self.film_node[i]
        load.control_node = self.control_node[i]
        load.temp_ambient = self.temp_ambient[i]
        load.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        node_id = np.unique(np.hstack([self.film_node, self.control_node]))
        geom_check(self,
                   missing,
                   node=(nid, node_id), filter_node0=True,
                   pconv_id=(self.model.pconv.pconv_id, self.pconv_id),
                   )

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.pconv_id.max(), self.control_node.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        #pconv_ids = array_str(self.pconv_id, size=size)
        film_nodes = array_default_int(self.film_node, default=0, size=size)
        control_nodes = array_default_int(self.control_node, default=0, size=size)
        temp_ambients = array_default_int(self.temp_ambient, default=0, size=size)
        #self.temp_ambient[icard] = ta
        for eid, pconv_id, film_node, control_node, tas in zip(element_ids, self.pconv_id,
                                                              film_nodes, control_nodes, temp_ambients):
            assert pconv_id > 0, pconv_id

            ta0 = tas[0]
            ta = [ta0]
            for tai in tas[1:]:
                ta.append(set_blank_if_default(tai, ta0))
            list_fields = ['CONV', eid, pconv_id, film_node, control_node] + ta
            bdf_file.write(print_card(list_fields))
        return


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
    _id_name = 'element_id'

    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.n = 0

    def add(self, eid: int, surface_type: str, nodes: list[int],
            iview_front: int=0, iview_back: int=0,
            rad_mid_front: int=0, rad_mid_back: int=0,
            comment: str='') -> int:
        """Creates a CHBDYG card"""
        self.cards.append((eid, surface_type, [iview_front, iview_back],
                           [rad_mid_front, rad_mid_back], nodes, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        # no field 2

        surface_type = string_or_blank(card, 3, 'Type', default='')
        assert len(surface_type) <= 5, surface_type
        iview_front = integer_or_blank(card, 4, 'iview_front', default=0)
        iview_back = integer_or_blank(card, 8, 'iview_back', default=0)
        rad_mid_front = integer_or_blank(card, 6, 'rad_mid_front', default=0)
        rad_mid_back = integer_or_blank(card, 7, 'rad_mid_back', default=0)
        # no field 8

        n = 1
        nodes = []
        for i in range(9, len(card)):
            grid = integer_or_blank(card, i, 'grid%d' % n, default=0)
            nodes.append(grid)  # used to have a None option
        assert len(nodes) > 0, 'card=%s' % card
        self.cards.append((eid, surface_type, [iview_front, iview_back],
                           [rad_mid_front, rad_mid_back], nodes, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype

        #: eid1: Surface element ID number for a side of an element. (0 < Integer < 100,000,000)
        #: eid2: A heat conduction element identification
        element_id = np.zeros(ncards, dtype=idtype)

        # : A consistent element side identification number
        # : (1 < Integer < 6)
        # self.side = np.zeros(ncards, dtype='int32')
        # if self.surface_type == 'AREA3':
        #     nnodes_required = 3
        #     nnodes_allowed = 3
        # elif self.surface_type == 'AREA4':
        #     nnodes_required = 4
        #     nnodes_allowed = 4
        # elif self.surface_type == 'AREA6':
        #     nnodes_required = 3
        #     nnodes_allowed = 6
        # elif self.surface_type == 'AREA8':
        #     nnodes_required = 4
        #     nnodes_allowed = 8
        surface_type = np.zeros(ncards, dtype='|U5')

        #: A VIEW entry identification number for the front face
        #: A VIEW entry identification number for the back face
        iview = np.zeros((ncards, 2), dtype='int32')

        #: RADM identification number for front face of surface element (Integer > 0)
        #: RADM identification number for back face of surface element (Integer > 0)
        rad_mid = np.zeros((ncards, 2), dtype='int32')

        nnode = np.zeros(ncards, dtype='int32')
        grid_list = []
        for icard, card in enumerate(self.cards):
            (eid, surface_typei, iviewi, rad_midi, nodesi, comment) = card
            element_id[icard] = eid
            surface_type[icard] = surface_typei
            iview[icard] = iviewi
            rad_mid[icard] = rad_midi
            nnode[icard] = len(nodesi)
            grid_list.extend(nodesi)
        grid = np.array(grid_list, dtype=idtype)
        self._save(element_id, surface_type, iview, rad_mid, nnode, grid)
        apply_bydor_default(self)
        self.cards = []

    def _save(self, element_id, surface_type, iview, rad_mid, nnode, grid):
        self.element_id = element_id
        self.surface_type = surface_type
        self.iview = iview
        self.rad_mid = rad_mid
        self.nnode = nnode
        self.grid = grid
        self.n = len(element_id)

    def __apply_slice__(self, load: CHBDYG, i: np.ndarray) -> None:
        load.element_id = self.element_id[i]
        load.surface_type = self.surface_type[i]
        load.iview = self.iview[i, :]
        load.rad_mid = self.rad_mid[i, :]
        load.grid = hslice_by_idim(i, self.igrid, self.grid)
        load.nnode = self.nnode[i]
        load.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['node_id'].append(self.grid)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        node_id = used_dict['node_id']
        ncards_removed = remove_unused_primary(
            self, node_id, self.grid, 'node_id')
        return ncards_removed

    @property
    def igrid(self) -> np.ndarray:
        return make_idim(self.n, self.nnode)

    @property
    def iview_front(self) -> np.ndarray:
        return self.iview[:, 0]

    @property
    def iview_back(self) -> np.ndarray:
        return self.iview[:, 1]

    @property
    def rad_mid_front(self) -> np.ndarray:
        return self.rad_mid[:, 0]

    @property
    def rad_mid_back(self) -> np.ndarray:
        return self.rad_mid[:, 1]

    @iview_front.setter
    def iview_front(self, iview_front: np.ndarray) -> None:
        self.iview[:, 0] = iview_front

    @iview_back.setter
    def iview_back(self, iview_back: np.ndarray) -> None:
        self.iview[:, 1] = iview_back

    @rad_mid_front.setter
    def rad_mid_front(self, rad_mid_front: np.ndarray) -> None:
        self.rad_mid[:, 0] = rad_mid_front

    @rad_mid_back.setter
    def rad_mid_back(self, rad_mid_back: np.ndarray) -> None:
        self.rad_mid[:, 1] = rad_mid_back

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.igrid.max(),
                   self.iview.max(), self.rad_mid.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        # TODO: default for iview, rad_mid is specified on BDYOR
        element_ids = array_str(self.element_id, size=size)
        iview_fronts   = array_default_int(self.iview[:, 0],   default=0, size=size)
        iview_backs    = array_default_int(self.iview[:, 1],   default=0, size=size)
        rad_mid_fronts = array_default_int(self.rad_mid[:, 0], default=0, size=size)
        rad_mid_backs  = array_default_int(self.rad_mid[:, 1], default=0, size=size)

        #i_view_front = set_blank_if_default(self.iview_front, 0)
        #i_view_back = set_blank_if_default(self.iview_back, 0)
        #rad_mid_front = set_blank_if_default(self.rad_mid_front, 0)
        #rad_mid_back = set_blank_if_default(self.rad_mid_back, 0)

        for eid, (inode0, inode1), surface_type, iview_front, iview_back, \
            rad_mid_front, rad_mid_back in zip(element_ids, self.igrid, self.surface_type,
                                               iview_fronts, iview_backs,
                                               rad_mid_fronts, rad_mid_backs):
            nodes = self.grid[inode0:inode1].tolist()
            list_fields = (['CHBDYG', eid, None, surface_type, iview_front,
                            iview_back, rad_mid_front, rad_mid_back, None, ] + nodes)
            bdf_file.write(print_card(list_fields))
        return


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
    _id_name = 'element_id'

    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.n = 0

    def add(self, eid: int, pid: int, surface_type: str,
            g1: int, g2: int,
            g0: int=0, gmid: int=0,
            ce: int=0,
            iview_front: int=0, iview_back: int=0,
            rad_mid_front: int=0, rad_mid_back: int=0,
            e1=None, e2=None, e3=None,
            comment: str='') -> int:
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
        g2 = 0 if g2 is None else g2
        gmid = 0 if gmid is None else gmid
        self.cards.append((eid, pid, surface_type,
                           [g1, g2, gmid], g0,
                           [iview_front, iview_back],
                           [rad_mid_front, rad_mid_back],
                           ce, [e1, e2, e3], comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        fdouble_or_blank = double_or_blank if self.model.is_strict_card_parser else force_double_or_blank
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=-1)
        surface_type = string_or_blank(card, 3, 'Type', default='')

        iview_front = integer_or_blank(card, 4, 'iview_front', default=0)
        iview_back = integer_or_blank(card, 5, 'iview_back', default=0)
        g1 = integer(card, 6, 'g1')

        if surface_type == 'POINT':
            g2 = blank(card, 7, 'g2')
            g2 = 0
        else:
            g2 = integer(card, 7, 'g2')

        g0 = integer_or_blank(card, 8, 'g0', default=0)
        rad_mid_front = integer_or_blank(card, 9, 'rad_mid_front', default=0)
        rad_mid_back = integer_or_blank(card, 10, 'rad_mid_back', default=0)
        gmid = integer_or_blank(card, 11, 'gmid', default=0)
        ce = integer_or_blank(card, 12, 'ce', default=-1)
        e1 = fdouble_or_blank(card, 13, 'e1', default=np.nan)
        e2 = fdouble_or_blank(card, 14, 'e2', default=np.nan)
        e3 = fdouble_or_blank(card, 15, 'e3', default=np.nan)
        assert len(card) <= 16, f'len(CHBDYP card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, surface_type,
                           [g1, g2, gmid], g0,
                           [iview_front, iview_back],
                           [rad_mid_front, rad_mid_back],
                           ce, [e1, e2, e3], comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype

        #: Surface element ID
        element_id = np.zeros(ncards, dtype=idtype)

        #: PHBDY property entry identification numbers. (Integer > 0)
        property_id = np.zeros(ncards, dtype='int32')
        surface_type = np.zeros(ncards, dtype='|U8')

        #: A VIEW entry identification number for the front/back face.
        iview = np.zeros((ncards, 2), dtype='int32')

        #: Grid point identification numbers of grids bounding the surface.
        #: (Integer > 0)
        # 0/1: g1/g2: Grid point identification numbers of grids bounding the surface
        # 2:   gmid:  Grid point identification number of a midside node if it is used
        #:            with the line type surface element.
        nodes = np.zeros((ncards, 3), dtype=idtype)

        #: Orientation grid point. (Integer > 0; Default = 0)
        g0 = np.zeros(ncards, dtype=idtype)

        #: RADM identification number for front/back face of surface element.
        #: (Integer > 0)
        rad_mid = np.zeros((ncards, 2), dtype='int32')

        #: Coordinate system for defining orientation vector.
        #: (Integer > 0; Default = 0
        ce = np.zeros(ncards, dtype='int32')

        #: Components of the orientation vector in coordinate system CE.
        #: The origin of the orientation vector is grid point G1.
        #: (Real or blank)
        ce_orientation = np.zeros((ncards, 3), dtype='float64')

        #self.nnode = np.zeros(ncards, dtype='int32')
        #grids = []
        for icard, card in enumerate(self.cards):
            (eid, pid, surface_typei, nodesi, g0i,
             iviewi, rad_midi,
             cei, ce_orientationi, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid

            nodes[icard, :] = nodesi
            if nodesi[1] == 0:
                assert surface_typei == 'POINT', f'eid={eid} surface_type={surface_typei!r} has g2=0'
            g0[icard] = g0i
            surface_type[icard] = surface_typei
            iview[icard] = iviewi
            rad_mid[icard, :] = rad_midi
            ce[icard] = cei
            ce_orientation[icard, :] = ce_orientationi
        self._save(element_id, property_id, nodes, g0,
                   surface_type, iview, rad_mid, ce, ce_orientation)
        apply_bydor_default(self)
        self.cards = []

    def _save(self, element_id, property_id, nodes, g0,
              surface_type, iview, rad_mid, ce, ce_orientation) -> None:
        assert isinstance(surface_type[0], str), surface_type.dtype.name
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.g0 = g0
        self.surface_type = surface_type
        self.iview = iview
        self.rad_mid = rad_mid
        self.ce = ce
        self.ce_orientation = ce_orientation
        self.n = len(element_id)

    def __apply_slice__(self, load: CHBDYP, i: np.ndarray) -> None:
        load.element_id = self.element_id[i]
        load.property_id = self.property_id[i]
        load.nodes = self.nodes[i, :]
        load.g0 = self.g0[i]
        load.surface_type = self.surface_type[i]
        load.iview = self.iview[i, :]
        load.rad_mid = self.rad_mid[i, :]
        load.ce = self.ce[i]
        load.ce_orientation = self.ce_orientation[i]
        load.n = len(i)

    @property
    def iview_front(self) -> np.ndarray:
        return self.iview[:, 0]

    @property
    def iview_back(self) -> np.ndarray:
        return self.iview[:, 1]

    @property
    def rad_mid_front(self) -> np.ndarray:
        return self.rad_mid[:, 0]

    @property
    def rad_mid_back(self) -> np.ndarray:
        return self.rad_mid[:, 1]

    @iview_front.setter
    def iview_front(self, iview_front: np.ndarray) -> None:
        self.iview[:, 0] = iview_front

    @iview_back.setter
    def iview_back(self, iview_back: np.ndarray) -> None:
        self.iview[:, 1] = iview_back

    @rad_mid_front.setter
    def rad_mid_front(self, rad_mid_front: np.ndarray) -> None:
        self.rad_mid[:, 0] = rad_mid_front

    @rad_mid_back.setter
    def rad_mid_back(self, rad_mid_back: np.ndarray) -> None:
        self.rad_mid[:, 1] = rad_mid_back

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max(),
                   self.iview.max(), self.rad_mid.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        g0s = array_str(self.g0, size=size)
        nodes = array_default_int(self.nodes, default=0, size=size)
        iview_fronts   = array_default_int(self.iview[:, 0],   default=0, size=size)
        iview_backs    = array_default_int(self.iview[:, 1],   default=0, size=size)
        rad_mid_fronts = array_default_int(self.rad_mid[:, 0], default=0, size=size)
        rad_mid_backs  = array_default_int(self.rad_mid[:, 1], default=0, size=size)
        ces  = array_default_int(self.ce, default=0, size=size)

        for eid, pid, surface_type, (g1, g2, gmid), g0, iview_front, iview_back, \
            rad_mid_front, rad_mid_back, ce, (e1, e2, e3) in zip_longest(
                element_ids, property_ids, self.surface_type, nodes, g0s,
                iview_fronts, iview_backs,
                rad_mid_fronts, rad_mid_backs,
                ces, self.ce_orientation):
            #(g1, g2, g0, gmid) = self.node_ids

            list_fields = ['CHBDYP', eid, pid, surface_type, iview_front,
                           iview_back, g1, g2, g0, rad_mid_front, rad_mid_back,
                           gmid, ce, e1, e2, e3]
            bdf_file.write(print_card(list_fields))
        return


def apply_bydor_default(chbdyx: CHBDYE | CHBDYG | CHBDYP) -> None:
    card_type = chbdyx.type
    model = chbdyx.model
    data_temp_default = []
    SURFACE_TYPE_DEFAULT = ''
    PROPERTY_ID_DEFAULT = -1
    CE_DEFAULT = -1
    #E123_DEFAULT = 0
    RAD_MID_DEFAULT = 0
    IVIEW_DEFAULT = 0

    bdyor = model.bdyor
    if bdyor is not None:
        #print(bdyor)
        pid_ = bdyor.pid
        surface_type_ = bdyor.surface_type
        #print('surface_type_ =', surface_type_)
        iview_front_ = bdyor.iview_front
        iview_back_ = bdyor.iview_back
        rad_mid_front_ = bdyor.rad_mid_front
        rad_mid_back_ = bdyor.rad_mid_back
        ce_ = bdyor.ce
        e1_, e2_, e3_ = bdyor.e123
        #model.log.info(f'setting default g0 as {g0_}')
        #model.log.info(f'setting default x_ as {x_}')
        #model.log.info(f'setting default offt_ as {offt_}')
    else:
        pid_ = 0
        surface_type_ = None
        iview_front_ = 0
        iview_back_ = 0
        rad_mid_front_ = 0
        rad_mid_back_ = 0
        ce_ = 0
        e1_, e2_, e3_ = 0., 0., 0.

    #PROPERTY_ID_DEFAULT = None
    #SURFACE_TYPE_DEFAULT = None
    if card_type == 'CHBDYP':
        data_temp_default.extend([
            # name,         data,               temp_value           default_value
            ('property_id', chbdyx.property_id, PROPERTY_ID_DEFAULT, pid_),
        ])
    if card_type in {'CHBDYG', 'CHBDYP'}:
        data_temp_default.append(
            ('surface_type', chbdyx.surface_type, SURFACE_TYPE_DEFAULT, surface_type_))

    data_temp_default += [
        # array to update,    placeholder_value,    default_value
        ('iview_front', chbdyx.iview_front, IVIEW_DEFAULT, iview_front_),
        ('iview_back', chbdyx.iview_back, IVIEW_DEFAULT, iview_back_),
        ('rad_mid_front', chbdyx.rad_mid_front, RAD_MID_DEFAULT, rad_mid_front_),
        ('rad_mid_back', chbdyx.rad_mid_back, RAD_MID_DEFAULT, rad_mid_back_),
    ]
    if card_type == 'CHBDYP':
        data_temp_default += [
            ('ce', chbdyx.ce, CE_DEFAULT, ce_),
            #(chbdyx.ce[:, 0], np.nan, ce_),
            ('ce_orientation1', chbdyx.ce_orientation[:, 0], np.nan, e1_),
            ('ce_orientation2', chbdyx.ce_orientation[:, 1], np.nan, e2_),
            ('ce_orientation3', chbdyx.ce_orientation[:, 2], np.nan, e3_),
        ]

    for name, data, temp_value, default_value in data_temp_default:
        #print(f'{name!r} temp_value={temp_value} default_value={default_value}')
        if isinstance(temp_value, (int, str)):
            ibad = np.where(data == temp_value)[0]
        elif np.isnan(temp_value):
            ibad = np.where(np.isnan(data))[0]
        else:
            ibad = np.where(data == temp_value)[0]
        if len(ibad):
            if isinstance(default_value, np.ndarray):
                data[ibad] = default_value[ibad]
            else:
                data[ibad] = default_value

    if card_type in {'CHBDYG', 'CHBDYP'}:
        empty_surface_type = (chbdyx.surface_type == '')
        assert not np.any(empty_surface_type), f'{card_type} surface_type={chbdyx.surface_type}'
    if card_type == 'CHBDYP':
        assert chbdyx.property_id.min() >= 0, f'{card_type} property_id={chbdyx.property_id}'
        assert chbdyx.ce.min() >= 0, f'{card_type} ce={chbdyx.ce}'
        assert np.abs(chbdyx.ce_orientation).max() >= 0.0, f'{card_type} ce_orientation={chbdyx.ce_orientation}'


class PHBDY(VectorizedBaseCard):
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
    _id_name = 'property_id'

    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.n = 0

    def add(self, pid: int,
            area_factor: float=None,
            d1=None, d2=None,
            comment: str='') -> int:
        """
        Creates a PHBDY card

        Parameters
        ----------
        pid : int
            property id
        area_factor : int
            Area factor of the surface used only for CHBDYP element
            Must be {POINT, LINE, TUBE, ELCYL}
            TUBE : constant thickness of hollow tube
        d1, d2 : float; default=None
            Diameters associated with the surface
            Used with CHBDYP [ELCYL, TUBE, FTUBE] surface elements
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, area_factor, d1, d2, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a PHBDY card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        fdouble_or_blank = double_or_blank if self.model.is_strict_card_parser else force_double_or_blank
        pid = integer(card, 1, 'pid')
        af = fdouble_or_blank(card, 2, 'af', default=np.nan)
        d1 = fdouble_or_blank(card, 3, 'd1', default=np.nan)
        d2 = fdouble_or_blank(card, 4, 'd2', default=d1)
        assert len(card) <= 5, f'len(PHBDY card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, af, d1, d2, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)

        #: Property identification number. (Unique Integer among all PHBDY
        #: entries). (Integer > 0)
        property_id = np.zeros(ncards, dtype='int32')

        #: Area factor of the surface used only for CHBDYP element
        #: TYPE = 'POINT', TYPE = 'LINE', TYPE = 'TUBE', or
        #: TYPE = 'ELCYL'. For TYPE = 'TUBE', AF is the constant thickness
        #: of the hollow tube. (Real > 0.0 or blank)
        area_factor = np.zeros(ncards, dtype='float64')

        #: Diameters associated with the surface. Used with CHBDYP element
        #: TYPE='ELCYL','TUBE','FTUBE'
        diameter =  np.zeros((ncards, 2), dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, af, d1, d2, comment) = card
            property_id[icard] = pid
            area_factor[icard] = af
            d2 = d1 if d2 is None else d2
            diameter[icard, :] = [d1, d2]
        self._save(property_id, area_factor, diameter)
        self.cards = []

    def _save(self, property_id, area_factor, diameter):
        assert len(self.property_id) == 0
        nproperties = len(property_id)
        self.property_id = property_id
        self.area_factor = area_factor
        self.diameter = diameter
        self.n = nproperties

    def __apply_slice__(self, prop: PHBDY, i: np.ndarray) -> None:
        prop.property_id = self.property_id[i]
        prop.area_factor = self.area_factor[i]
        prop.diameter = self.diameter[i, :]
        prop.n = len(i)

    @property
    def max_id(self) -> int:
        return self.property_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_ids = array_str(self.property_id, size=size)
        for pid, af, d1, d2 in zip(property_ids, self.area_factor,
                                   self.diameter[:, 0], self.diameter[:, 1]):
            d2 = set_blank_if_default(d2, d1)
            list_fields = ['PHBDY', pid, af, d1, d2]
            bdf_file.write(print_card(list_fields))
        return

    # @property
    # def allowed_materials(self) -> list[Any]:
    #     return [mat for mat in [self.model.mat1] if mat.n > 0]
    #
    # def mass_per_length(self) -> np.ndarray:
    #     return line_mid_mass_per_length(self.material_id, self.nsm, self.A,
    #                                     self.allowed_materials)
    #
    # def area(self) -> np.ndarray:
    #     return self.A


# class CONV(VectorizedBaseCard):
#     """
#     Specifies a free convection boundary condition for heat transfer analysis
#     through connection to a surface element (CHBDYi entry).
#
#     """
#     def clear(self) -> None:
#         self.element_id = np.array([], dtype='int32')
#         self.n = 0
#
#     @VectorizedBaseCard.parse_cards_check
#     def parse_cards(self) -> None:
#         ncards = len(self.cards)
#         if ncards == 0:
#             return
#
#         #: CHBDYG, CHBDYE, or CHBDYP surface element identification number.
#         #: (Integer > 0)
#         self.element_id = np.zeros(ncards, dtype='int32')
#
#         #: Convection property identification number of a PCONV entry
#         self.pconv_id = np.zeros(ncards, dtype='int32')
#
#         #: Point for film convection fluid property temperature
#         self.film_node = np.zeros(ncards, dtype='int32')
#
#         #: Control point for free convection boundary condition.
#         self.control_node = np.zeros(ncards, dtype='int32')
#
#         #: Ambient points used for convection 0's are allowed for TA2 and
#         #: higher.  (Integer > 0 for TA1 and Integer > 0 for TA2 through TA8;
#         #: Default for TA2 through TA8 is TA1.)
#         self.temp_ambient = np.zeros((ncards, 8), dtype='int32')
#
#         #grids = []
#         for icard, card_comment in enumerate(self.cards):
#             card, comment = card_comment
#
#             eid = integer(card, 1, 'eid')
#             pconid = integer(card, 2, 'pconid')
#             film_node = integer_or_blank(card, 3, 'film_node', 0)
#             control_node = integer_or_blank(card, 4, 'control_node', 0)
#
#             ta1 = integer(card, 5, 'TA1')
#             assert ta1 > 0, ta1
#
#             ta2 = integer_or_blank(card, 6, 'ta2', ta1)
#             ta3 = integer_or_blank(card, 7, 'ta3', ta1)
#             ta4 = integer_or_blank(card, 8, 'ta4', ta1)
#             ta5 = integer_or_blank(card, 9, 'ta5', ta1)
#             ta6 = integer_or_blank(card, 10, 'ta6', ta1)
#             ta7 = integer_or_blank(card, 11, 'ta7', ta1)
#             ta8 = integer_or_blank(card, 12, 'ta8', ta1)
#             ta = [ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8]
#
#             self.element_id[icard] = eid
#             self.pconv_id[icard] = pconid
#             self.film_node[icard] = film_node
#             self.control_node[icard] = control_node
#             self.temp_ambient[icard] = ta
#             assert len(card) <= 13, f'len(CONV card) = {len(card):d}\ncard={card}'
#
#     def write(self, size: int=8) -> str:
#         if len(self.element_id) == 0:
#             return ''
#         lines = []
#         if size == 8:
#             print_card = print_card_8
#
#         element_ids = array_str(self.element_id, size=size)
#         #pconv_ids = array_str(self.pconv_id, size=size)
#         film_nodes = array_default_int(self.film_node, default=0, size=size)
#         control_nodes = array_default_int(self.control_node, default=0, size=size)
#         temp_ambients = array_default_int(self.temp_ambient, default=0, size=size)
#         #self.temp_ambient[icard] = ta
#         for eid, pconv_id, film_node, control_node, tas in zip(element_ids, self.pconv_id,
#                                                               film_nodes, control_nodes, temp_ambients):
#             assert pconv_id > 0, pconv_id
#
#             ta0 = tas[0]
#             ta = [ta0]
#             for tai in tas[1:]:
#                 ta.append(set_blank_if_default(tai, ta0))
#             list_fields = ['CONV', eid, pconv_id, film_node, control_node] + ta
#
#             lines.append(print_card(list_fields))
#         return ''.join(lines)


class PCONV(VectorizedBaseCard):
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
    _id_name = 'pconv_id'

    def clear(self) -> None:
        self.pconv_id = np.array([], dtype='int32')
        self.n = 0

    def add(self, pconv_id: int, mid: int=None,
            form: int=0, exponent_free_convection=0.0,
            free_convection_type: int=0,
            table_id: int=0,
            chlen: float=None, gidin: int=None, ce: int=0,
            e: list[float]=None,
            comment: str='') -> int:
        """
        Creates a PCONV card

        Parameters
        ----------
        pconv_id : int
            Convection property ID
        mid : int; default=None
            Material ID
        form : int; default=0
            Type of formula used for free convection
            Must be {0, 1, 10, 11, 20, or 21}
        exponent_free_convection : float; default=0.0
            Free convection exponent as implemented within the context
            of the particular form that is chosen
        free_convection_type : int; default=0
            Formula type for various configurations of free convection
        table_id : int; default=None
            Identification number of a TABLEHT entry that specifies the
            two variable tabular function of the free convection heat
            transfer coefficient
        chlen : float; default=None
            Characteristic length
        gidin : int; default=None
            Grid ID of the referenced inlet point
        ce : int; default=0
            Coordinate system for defining orientation vector.
        e : list[float]; default=None
            Components of the orientation vector in coordinate system CE.
            The origin of the orientation vector is grid point G1
        comment : str; default=''
            a comment for the card

        """
        mid = 0 if mid is None else mid
        table_id = 0 if table_id is None else table_id
        gidin = 0 if gidin is None else gidin
        coord_e = ce
        self.cards.append((pconv_id, mid, form, exponent_free_convection,
                           free_convection_type, table_id,
                           chlen, gidin, coord_e, e, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        fdouble_or_blank = double_or_blank if self.model.is_strict_card_parser else force_double_or_blank
        pconv_id = integer(card, 1, 'pconv_id')
        mid = integer_or_blank(card, 2, 'mid')
        form = integer_or_blank(card, 3, 'form', default=0)
        exponent_free_convection = fdouble_or_blank(card, 4, 'expf', default=0.0)
        free_convection_type = integer_or_blank(card, 5, 'ftype', default=0)
        table_id = integer_or_blank(card, 6, 'tid', default=0)
        chlen = fdouble_or_blank(card, 9, 'chlen', default=np.nan)
        gidin = integer_or_blank(card, 10, 'gidin', default=0)
        coord_e = integer_or_blank(card, 11, 'ce', default=0)
        e = [
            fdouble_or_blank(card, 12, 'e1', default=np.nan),
            fdouble_or_blank(card, 13, 'e2', default=np.nan),
            fdouble_or_blank(card, 14, 'e3', default=np.nan),
        ]
        assert len(card) <= 15, f'len(PCONV card) = {len(card):d}\ncard={card}'
        self.cards.append((pconv_id, mid, form, exponent_free_convection,
                           free_convection_type, table_id,
                           chlen, gidin, coord_e, e, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)

        #: Convection property identification number. (Integer > 0)
        pconv_id = np.zeros(ncards, dtype='int32')

        #: Material property identification number. (Integer > 0)
        material_id = np.zeros(ncards, dtype='int32')

        #: Type of formula used for free convection.
        #: (Integer 0, 1, 10, 11, 20, or 21)
        form = np.zeros(ncards, dtype='int32')

        #: Free convection exponent as implemented within the context of the
        #: particular form that is chosen
        exponent_free_convection = np.zeros(ncards, dtype='float64')

        #: Formula type for various configurations of free convection
        free_convection_type = np.zeros(ncards, dtype='int32')

        #: Identification number of a TABLEHT entry that specifies the two
        #: variable tabular function of the free convection heat transfer
        #: coefficient
        table_id = np.zeros(ncards, dtype='int32')

        #: Characteristic length
        characteristic_length = np.zeros(ncards, dtype='float64')

        #: Grid ID of the referenced inlet point
        grid_inlet = np.zeros(ncards, dtype='int32')

        #: Coordinate system for defining orientation vector.
        #: (Integer > 0;Default = 0
        coord_e = np.zeros(ncards, dtype='int32')

        #: Components of the orientation vector in coordinate system CE. The
        #: origin of the orientation vector is grid point G1. (Real or blank)
        e = np.zeros((ncards, 3), dtype='float64')

        #assert self.pconid > 0
        #assert mid is None or self.mid > 0
        #assert self.form in [0, 1, 10, 11, 20, 21]

        for icard, card in enumerate(self.cards):
            (pconid, mid, formi,
             exponent_free_convectioni, free_convection_typei,
             table_idi, chleni, gidini, coord_ei, e123i, comment) = card

            pconv_id[icard] = pconid
            material_id[icard] = mid
            form[icard] = formi
            exponent_free_convection[icard] = exponent_free_convectioni
            free_convection_type[icard] = free_convection_typei
            table_id[icard] = table_idi
            characteristic_length[icard] = chleni
            grid_inlet[icard] = gidini
            coord_e[icard] = coord_ei
            if e123i is not None:
                e[icard] = e123i
        self._save(pconv_id, material_id, form,
                   exponent_free_convection,
                   free_convection_type,
                   table_id, characteristic_length,
                   grid_inlet, coord_e, e)
        self.cards = []

    def _save(self, pconv_id, material_id, form,
              exponent_free_convection,
              free_convection_type,
              table_id, characteristic_length,
              grid_inlet, coord_e, e) -> None:
        assert len(self.pconv_id) == 0
        self.pconv_id = pconv_id
        self.material_id = material_id
        self.form = form
        self.exponent_free_convection = exponent_free_convection

        nproperties = len(pconv_id)
        self.free_convection_type = free_convection_type
        self.table_id = table_id
        self.characteristic_length = characteristic_length
        self.grid_inlet = grid_inlet
        self.coord_e = coord_e
        self.e = e
        self.n = nproperties

    def _save_nx(self, pconv_id, material_id, form, exponent_free_convection):
        assert len(self.pconv_id) == 0
        self.pconv_id = pconv_id
        self.material_id = material_id
        self.form = form
        self.exponent_free_convection = exponent_free_convection

        nproperties = len(pconv_id)
        self.free_convection_type = np.zeros(nproperties, dtype='int32')
        self.table_id = np.zeros(nproperties, dtype='int32')
        self.characteristic_length = np.zeros(nproperties, dtype='float64')
        self.grid_inlet = np.zeros(nproperties, dtype='int32')
        self.coord_e = np.zeros(nproperties, dtype='int32')
        self.e = np.zeros((nproperties, 3), dtype='float64')
        self.n = nproperties

    def __apply_slice__(self, prop: PCONV, i: np.ndarray) -> None:
        prop.pconv_id = self.pconv_id[i]
        prop.material_id = self.material_id[i]
        prop.form = self.form[i]
        prop.exponent_free_convection = self.exponent_free_convection[i]

        prop.free_convection_type = self.free_convection_type[i]
        prop.table_id = self.table_id[i]
        prop.characteristic_length = self.characteristic_length[i]
        prop.grid_inlet = self.grid_inlet[i]
        prop.coord_e = self.coord_e[i]
        prop.e = self.e[i]
        prop.n = len(i)

    @property
    def max_id(self) -> int:
        return max(self.pconv_id.max(), self.material_id.max(), self.table_id.max(),
                   self.grid_inlet.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        pconv_ids = array_str(self.pconv_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        forms = array_default_int(self.form, default=0.0, size=size)
        exponent_free_convections = array_default_float(self.exponent_free_convection, default=0, size=size)
        free_convection_types = array_default_int(self.free_convection_type, default=0, size=size)
        coord_ids = array_default_int(self.coord_e, default=0, size=size)
        grid_inlets = array_str(self.grid_inlet, size=size)
        table_ids = array_str(self.table_id, size=size)
        e1s = self.e[:, 0]
        e2s = self.e[:, 1]
        e3s = self.e[:, 2]
        for pconv_id, mid, form, expf, \
            ftype, table_id, characteristic_length, grid_inlet, \
            coord_e, e1, e2, e3 in zip(pconv_ids, material_ids, forms, exponent_free_convections,
                                       free_convection_types, table_ids, self.characteristic_length,
                                       grid_inlets, coord_ids, e1s, e2s, e3s):
            list_fields = ['PCONV', pconv_id, mid, form, expf, ftype, table_id,
                           None, None, characteristic_length, grid_inlet,
                           coord_e, e1, e2, e3]
            bdf_file.write(print_card(list_fields))
        return

    # @property
    # def allowed_materials(self) -> list[Any]:
    #     return [mat for mat in [self.model.mat1] if mat.n > 0]
    #
    # def mass_per_length(self) -> np.ndarray:
    #     return line_mid_mass_per_length(self.material_id, self.nsm, self.A,
    #                                     self.allowed_materials)
    #
    # def area(self) -> np.ndarray:
    #     return self.A


class CONVM(VectorizedBaseCard):
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
    _id_name = 'element_id'

    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.n = 0

    def add(self, eid: int, pconvm_id: int, ta1: int,
            film_node: int=0, cntmdot: int=0,
            ta2: int=0, mdot: float=1.0,
            ifile: int=0, comment: str='') -> int:
        """
        Creates a CONVM card

        Parameters
        ----------
        eid : int
            element id (CHBDYP)
        pconvm_id : int
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
            point associated with the CONTROL_NODE field
            required if cntmdot = 0
        comment : str; default=''
            a comment for the card

        """
        ta2 = 0 if ta2 is None else ta2
        self.cards.append((eid, pconvm_id, film_node, cntmdot, [ta1, ta2], mdot, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a CONVM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        fdouble_or_blank = double_or_blank if self.model.is_strict_card_parser else force_double_or_blank
        eid = integer(card, 1, 'eid')
        pconvm = integer(card, 2, 'pconvm')
        film_node = integer_or_blank(card, 3, 'film_node', default=0)
        cntmdot = integer_or_blank(card, 4, 'cntmdot', default=0)
        ta1 = integer(card, 5, 'ta1')
        ta2 = integer_or_blank(card, 6, 'ta2', default=ta1)
        mdot = fdouble_or_blank(card, 7, 'mdot', default=1.0)
        assert len(card) <= 8, f'len(CONVM card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pconvm, film_node, cntmdot, [ta1, ta2], mdot, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)

        ifile = np.zeros(ncards, dtype='int32')
        # Parameters
        # ----------
        # eid : int
        #     element id (CHBDYP)
        element_id = np.zeros(ncards, dtype='int32')

        # pconid : int
        #     property ID (PCONVM)
        pconvm_id = np.zeros(ncards, dtype='int32')

        # ta1 : int
        #     ambient point for convection
        # ta2 : int; default=None
        #     None : ta1
        #     ambient point for convection
        temp_ambient = np.zeros((ncards, 2), dtype='int32')

        # film_node : int; default=0
        film_node = np.zeros(ncards, dtype='int32')

        # cntmdot : int; default=0
        #     control point used for controlling mass flow
        #     0/blank is only allowed when mdot > 0
        control_node_mdot = np.zeros(ncards, dtype='int32')

        # mdot : float; default=1.0
        #     a multiplier for the mass flow rate in case there is no
        #     point associated with the CONTROL_NODE field
        #     required if cntmdot = 0
        mdot = np.zeros(ncards, dtype='float64')

        # comment : str; default=''
        #     a comment for the card

        #grids = []
        comment = {}
        for icard, card in enumerate(self.cards):
            (eid, pconvm, film_nodei, cntmdot, ta, mdoti, ifilei, commenti) = card
            ifile[icard] = ifilei
            element_id[icard] = eid
            pconvm_id[icard] = pconvm
            film_node[icard] = film_nodei
            #control_node[icard] = control_node
            temp_ambient[icard] = ta
            control_node_mdot[icard] = cntmdot
            mdot[icard] = mdoti
        self._save(element_id, pconvm_id, film_node, temp_ambient, control_node_mdot, mdot,
                   ifile, comment)
        self.cards = []

    def _save(self, element_id, pconvm_id, film_node, temp_ambient, control_node_mdot, mdot,
              ifile=None, comment=None):
        ncards = len(element_id)
        if mdot is None:
            mdot = np.ones(ncards, dtype='float64')
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        assert len(self.element_id) == 0
        save_ifile_comment(self, ifile, comment)
        self.element_id = element_id
        self.pconvm_id = pconvm_id
        self.film_node = film_node
        self.temp_ambient = temp_ambient
        self.control_node_mdot = control_node_mdot
        self.mdot = mdot
        self.n = len(ifile)

    def __apply_slice__(self, elem: CONVM, i: np.ndarray) -> None:
        elem.element_id = self.element_id[i]
        elem.pconvm_id = self.pconvm_id[i]
        elem.film_node = self.film_node[i]
        elem.temp_ambient = self.temp_ambient[i]
        elem.control_node_mdot = self.control_node_mdot[i]
        elem.mdot = self.mdot[i]
        elem.n = len(i)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.pconvm_id.max(),
                   self.film_node.max(), self.control_node_mdot.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        pconvm_ids = array_str(self.pconvm_id, size=size)
        film_nodes = array_default_int(self.film_node, default=0, size=size)
        control_node_mdots = array_default_int(self.control_node_mdot, default=0, size=size)
        temp_ambients = array_default_int(self.temp_ambient, default=0, size=size)
        #self.temp_ambient[icard] = ta
        for eid, pconvm_id, film_node, control_node_mdot, \
            (ta1, ta2), mdot in zip(element_ids, pconvm_ids, film_nodes,
                                    control_node_mdots, temp_ambients, self.mdot):
            #assert pconvm_id > 0, pconvm_id

            ta2 = set_blank_if_default(ta2, ta1)
            #mdot = set_blank_if_default(self.mdot, 1.0)
            list_fields = ['CONVM', eid, pconvm_id, film_node,
                           control_node_mdot, ta1, ta2, mdot]
            bdf_file.write(print_card(list_fields))
        return


class PCONVM(VectorizedBaseCard):
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
    _id_name = 'pconvm_id'

    def clear(self) -> None:
        self.pconvm_id = np.array([], dtype='int32')
        self.n = 0

    def add(self, pconvm_id: int, mid: int, coeff: float,
            form: int=0, flag: int=0,
            expr: float=0.0,
            exppi: float=0.0, exppo: float=0.0,
            comment: str='') -> int:
        """
        Creates a PCONVM card

        Parameters
        ----------
        pconvm_id : int
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
        self.cards.append((pconvm_id, mid, form, flag, coeff, expr, exppi, exppo, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        fdouble = double if self.model.is_strict_card_parser else force_double
        fdouble_or_blank = double_or_blank if self.model.is_strict_card_parser else force_double_or_blank
        pconvm_id = integer(card, 1, 'pconid')
        mid = integer(card, 2, 'mid')
        form = integer_or_blank(card, 3, 'form', default=0)
        flag = integer_or_blank(card, 4, 'flag', default=0)
        coeff = fdouble(card, 5, 'coef')
        expr = fdouble_or_blank(card, 6, 'expr', default=0.0)
        exppi = fdouble_or_blank(card, 7, 'exppi', default=0.0)
        exppo = fdouble_or_blank(card, 8, 'exppo', default=0.0)
        assert len(card) <= 9, f'len(PCONVM card) = {len(card):d}\ncard={card}'
        self.cards.append((pconvm_id, mid, form, flag, coeff, expr, exppi, exppo, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)

        # Parameters
        # ----------
        # pconid : int
        #     Convection property ID
        pconvm_id = np.zeros(ncards, dtype='int32')

        # mid: int
        #     Material ID
        material_id = np.zeros(ncards, dtype='int32')

        # coeff: float
        #     Constant coefficient used for forced convection
        coefficient = np.zeros(ncards, dtype='float64')

        # form: int; default=0
        #     Type of formula used for free convection
        #     Must be {0, 1, 10, 11, 20, or 21}
        form = np.zeros(ncards, dtype='int32')

        # flag: int; default=0
        #     Flag for mass flow convection
        flag = np.zeros(ncards, dtype='int32')

        # expr: float; default=0.0
        #     Reynolds number convection exponent
        expr = np.zeros(ncards, dtype='float64')

        # exppi: float; default=0.0
        #     Prandtl number convection exponent for heat transfer into
        #     the working fluid
        exppi = np.zeros(ncards, dtype='float64')

        # exppo: float; default=0.0
        #     Prandtl number convection exponent for heat transfer out of
        #     the working fluid
        exppo = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pconid, mid, formi, flagi, coeffi, expri, exppii, exppoi, comment) = card
            pconvm_id[icard] = pconid
            material_id[icard] = mid
            form[icard] = formi
            flag[icard] = flagi
            coefficient[icard] = coeffi
            expr[icard] = expri
            exppi[icard] = exppii
            exppo[icard] = exppoi
        self._save(pconvm_id, material_id, form, flag, coefficient, expr, exppi, exppo)
        self.cards = []

    def _save(self, pconvm_id, material_id, form, flag, coefficient, expr, exppi, exppo):
        self.pconvm_id = pconvm_id
        self.material_id = material_id
        self.form = form
        self.flag = flag
        self.coefficient = coefficient
        self.expr = expr
        self.exppi = exppi
        self.exppo = exppo
        self.n = len(pconvm_id)

    def __apply_slice__(self, prop: PCONVM, i: np.ndarray) -> None:
        prop.pconvm_id = self.pconvm_id[i]
        prop.material_id = self.material_id[i]
        prop.form = self.form[i]
        prop.flag = self.flag[i]
        prop.coefficient = self.coefficient[i]
        prop.expr = self.expr[i]
        prop.exppi = self.exppi[i]
        prop.exppo = self.exppo[i]
        prop.n = len(i)

    @property
    def max_id(self) -> int:
        return max(self.pconvm_id.max(), self.material_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        pconvm_ids = array_str(self.pconvm_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        forms = array_default_int(self.form, default=0, size=size)
        flags = array_default_int(self.flag, default=0, size=size)

        coefficients = array_float(self.coefficient, size=size, is_double=False)
        exprs = array_default_float(self.expr, default=0.0, size=size, is_double=False)
        exppis = array_default_float(self.exppi, default=0.0, size=size, is_double=False)
        exppos = array_default_float(self.exppo, default=0.0, size=size, is_double=False)
        #form = set_blank_if_default(self.form, 0)
        #flag = set_blank_if_default(self.flag, 0)
        for pconvm_id, mid, form, flag, \
            coeff, expr, exppi, exppo in zip_longest(pconvm_ids, material_ids, forms, flags,
                                                     coefficients, exprs, exppis, exppos):
            assert pconvm_id != '0', pconvm_id
            list_fields = ['PCONVM', pconvm_id, mid, form, flag,
                           coeff, expr, exppi, exppo]
            bdf_file.write(print_card(list_fields))
        return
