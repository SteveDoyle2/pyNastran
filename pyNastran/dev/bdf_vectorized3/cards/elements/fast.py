from __future__ import annotations
from itertools import zip_longest
from typing import Any, TYPE_CHECKING
import numpy as np
#from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string, # blank,
    integer_or_blank, double_or_blank, # string_or_blank,
    #integer_double_or_blank,
    #fields,
)
#from pyNastran.bdf.cards.elements.bars import set_blank_if_default

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, get_print_card_8_16,
    parse_element_check, parse_property_check)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int
#from .rod import line_length, line_centroid, line_centroid_with_spoints
#from .utils import get_mass_from_property
if TYPE_CHECKING:
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class CFAST(Element):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')

    def add(self, eid: int, pid: int, fast_type: str,
            ids: list[int],
            gs=None, ga=None, gb=None,
            xs=None, ys=None, zs=None, comment: str='') -> int:
        """Creates a CFAST card"""
        ga = 0 if ga is None else ga
        gb = 0 if gb is None else gb
        gs = 0 if gs is None else gs
        self.cards.append((eid, pid, fast_type, ids, [ga, gb], gs, [xs, ys, zs], comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CFAST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        fast_type = string(card, 3, 'Type')
        ida = integer(card, 4, 'ida')
        idb = integer(card, 5, 'idb')
        gs = integer_or_blank(card, 6, 'gs', default=0)
        ga = integer_or_blank(card, 7, 'ga', default=0)
        gb = integer_or_blank(card, 8, 'gb', default=0)
        xs = double_or_blank(card, 9, 'xs')
        ys = double_or_blank(card, 10, 'ys')
        zs = double_or_blank(card, 11, 'zs')
        assert len(card) <= 12, f'len(CFAST card) = {len(card):d}\ncard={card}'

        #if self.fast_type not in ['PROP', 'ELEM']:
            #msg = f'CFAST; eid={self.eid} Type={self.Type!r} must be in [PROP, ELEM]'
            #raise TypeError(msg)

        #gab_is_none = self.ga is None or self.gb is None
        #xyz_is_none = self.xs is None or self.ys is None or self.zs is None
        #if self.gs is None and gab_is_none and xyz_is_none:
            #msg = ('CFAST; eid=%s; gs=%s is not an integer or\n'
                   #'              [ga=%s, gb=%s] are not integers or \n'
                   #'              [xs=%s, ys=%s, zs=%s] are not floats' % (
                       #self.eid, self.gs, self.ga, self.gb, self.xs, self.ys, self.zs))
            #raise ValueError(msg)
        #if self.Type=='PROP': # PSHELL/PCOMP  ida & idb
        #return CFAST(eid, pid, Type, ida, idb, gs=gs, ga=ga, gb=gb,
                     #xs=xs, ys=ys, zs=zs, comment=comment)
        self.cards.append((eid, pid, fast_type, [ida, idb],
                           [ga, gb], gs, [xs, ys, zs], comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)

        # If both GA and GB are specified, they must lie on or at least
        # have projections onto surface patches A and B respectively.
        # If GA and GB are both specified, GS is ignored. By default,
        # the locations of user specified GA and GB will not be changed.
        # If the user specifies "SWLDPRM, MOVGAB, 1,", then the locations
        # will be corrected so that they lie on the surface patches A and
        # B within machine precision. The length of the fastener is the
        # final distance between GA and GB. If the length is zero, the
        # normal to patch A is used to define the axis of the fastener.
        nodes = np.zeros((ncards, 2), dtype=idtype)

        fast_type = np.zeros(ncards, dtype='|U4')
        ids = np.zeros((ncards, 2), dtype='int32')
        fastener_node = np.zeros(ncards, dtype='int32')
        fastener_xyz = np.zeros((ncards, 3), dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, pid, fast_typei, idsi, nodesi, gs, xyz, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nodesi
            fast_type[icard] = fast_typei
            ids[icard, :] = idsi
            fastener_node[icard] = gs
            fastener_xyz[icard, :] = xyz
        self._save(element_id, property_id, nodes, fast_type, ids, fastener_node, fastener_xyz)
        self.cards = []

    def _save(self, element_id, property_id, nodes, fast_type, ids, fastener_node, fastener_xyz):
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.fast_type = fast_type
        self.ids = ids
        self.fastener_node = fastener_node
        self.fastener_xyz = fastener_xyz

    def __apply_slice__(self, elem: CFAST, i: np.ndarray) -> None:
        self._slice_comment(elem, i)

        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.fast_type = self.fast_type[i]
        elem.ids = self.ids[i, :]
        elem.fastener_node = self.fastener_node[i]
        elem.fastener_xyz = self.fastener_xyz[i, :]
        elem.n = len(i)

    def set_from_op2(self, element_id, property_id, gs, elem_grid_flag, nodes):
        assert element_id.min() > 0, element_id
        assert property_id.min() > 0, property_id
        #assert nodes.min() > 0, nodes
        nelement = len(element_id)

        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.fast_type = elem_grid_flag

        self.ids = np.zeros((nelement, 2), dtype='int32')
        self.fastener_node = gs
        self.fastener_xyz = np.zeros((nelement, 3), dtype='float64')
        self.n = nelement

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        gab_nodes = array_default_int(self.nodes, default=0, size=size)
        fastener_node = array_default_int(self.fastener_node, default=0, size=size)
        for eid, pid, fast_type, (ida, idb), (ga, gb), gs, (xs, ys, zs) in zip_longest(
                element_ids, property_ids, self.fast_type, self.ids,
                gab_nodes, fastener_node, self.fastener_xyz):

            xs = None if np.isnan(xs) else xs
            ys = None if np.isnan(ys) else ys
            zs = None if np.isnan(zs) else zs
            list_fields = ['CFAST', eid, pid, fast_type, ida, idb,
                           gs, ga, gb, xs, ys, zs]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pbush]
                if prop.n > 0]

    def length(self) -> np.ndarray:
        self.model.log.warning('CFAST length is not supported')
        length = np.zeros(len(self.element_id), dtype='float64')
        return length


class PFAST(Property):
    """
    +------+-----+-----+------+-------+---------+--------+--------+-----+
    |   1  |  2  |  3  |   4  |   5   |    6    |    7   |    8   |  9  |
    +======+=====+=====+======+=======+=========+========+========+=====+
    |PFAST | PID |  D  | MCID | MFLAG |   KT1   |   KT2  |   KT3  | KR1 |
    +------+-----+-----+------+-------+---------+--------+--------+-----+
    |      | KR2 | KR3 | MASS |   GE  |         |        |        |     |
    +------+-----+-----+------+-------+---------+--------+--------+-----+
    |PFAST |  7  | 1.1 | 70   |       | 100000. | 46000. | 12300. |     |
    +------+-----+-----+------+-------+---------+--------+--------+-----+
    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.zeros(0, dtype='int32')
        self.diameter = np.zeros(0, dtype='float64')
        self.kt = np.zeros((0, 3), dtype='float64')
        self.kr = np.zeros((0, 3), dtype='float64')
        self._mass = np.zeros(0, dtype='float64')
        self.ge = np.zeros(0, dtype='float64')
        self.mflag = np.zeros(0, dtype='int32')
        self.coord_id = np.zeros(0, dtype='int32')

    def add(self, pid: int, d: int,
            kt1: float, kt2: float, kt3: float,
            mcid: int=-1, mflag: int=0,
            kr1: float=0., kr2: float=0., kr3: float=0.,
            mass: float=0., ge: float=0.,
            comment: str='') -> int:
        """
        Creates a PAST card

        Parameters
        ----------
        pid : int
            property id
        d : int
            diameter of the fastener
        kt1, kt2, kt3 : float
            stiffness values in directions 1-3
        mcid : int; default=01
            specifies the element stiffness coordinate system
        mflag : int; default=0
            0-absolute; 1-relative
        kr1, kr2, kr3 : float; default=0.0
            rotational stiffness values in directions 1-3
        mass : float; default=0.0
            lumped mass of the fastener
        ge : float; default=0.0
            structural damping
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, d,
                           kt1, kt2, kr3,
                           kr1, kr2, kr3,
                           mass, ge,
                           mcid, mflag, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PFAST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        d = double(card, 2, 'd')
        mcid = integer_or_blank(card, 3, 'mcid', default=-1)
        mflag = integer_or_blank(card, 4, 'mflag', default=0)

        kt1 = double(card, 5, 'kt1')
        kt2 = double(card, 6, 'kt2')
        kt3 = double(card, 7, 'kt3')

        kr1 = double_or_blank(card, 8, 'kr1', default=0.0)
        kr2 = double_or_blank(card, 9, 'kr2', default=0.0)
        kr3 = double_or_blank(card, 10, 'kr3', default=0.0)
        mass = double_or_blank(card, 11, 'mass', default=0.0)
        ge = double_or_blank(card, 12, 'ge', default=0.0)
        assert len(card) <= 13, f'len(PFAST card) = {len(card):d}\ncard={card}'
        #return PFAST(pid, d, kt1, kt2, kt3, mcid=mcid, mflag=mflag,
                     #kr1=kr1, kr2=kr2, kr3=kr3, mass=mass, ge=ge,
                     #comment=comment)
        self.cards.append((pid, d,
                           kt1, kt2, kt3,
                           kr1, kr2, kr3,
                           mass, ge,
                           mcid, mflag, comment))
        self.n += 1
        return self.n

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        diameter = np.zeros(ncards, dtype='float64')
        kt = np.zeros((ncards, 3), dtype='float64')
        kr = np.zeros((ncards, 3), dtype='float64')
        mass = np.zeros(ncards, dtype='float64')
        ge = np.zeros(ncards, dtype='float64')
        mflag = np.zeros(ncards, dtype='int32')
        coord_id = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (pid, d, kt1, kt2, kt3, kr1, kr2, kr3,
             massi, gei, mcid, mflagi, comment) = card

            property_id[icard] = pid
            diameter[icard] = d
            coord_id[icard] = mcid
            mflag[icard] = mflagi
            kt[icard, :] = [kt1, kt2, kt3]
            kr[icard, :] = [kr1, kr2, kr3]
            mass[icard] = massi
            ge[icard] = gei
        self._save(property_id, diameter, kt, kr, mass, coord_id, mflag, ge)
        self.sort()
        self.cards = []

    def _save(self, property_id, diameter, kt, kr, mass, coord_id, mflag, ge):
        self.property_id = property_id
        self.diameter = diameter
        self.kt = kt
        self.kr = kr
        self._mass = mass
        self.coord_id = coord_id
        self.mflag = mflag
        self.ge = ge

    def __apply_slice__(self, prop: PFAST, i: np.ndarray) -> None:
        self._slice_comment(prop, i)

        prop.property_id = self.property_id[i]
        prop.diameter = self.diameter[i]
        prop.kt = self.kt[i]
        prop.kr = self.kr[i]
        prop.coord_id = self.coord_id[i]
        prop.coord_id = self.coord_id[i]
        prop.mflag = self.mflag[i]
        prop.ge = self.ge[i]
        prop.n = len(i)

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        for pid, diameter, (kt1, kt2, kt3), (kr1, kr2, kr3), mcid, mflag, mass, ge in zip_longest(
                self.property_id, self.diameter, self.kt, self.kr,
                self.coord_id, self.mflag, self._mass, self.ge):

            #mass = None if mass == 0. else mass
            list_fields = ['PFAST', pid, diameter, mcid, mflag, kt1,
                           kt2, kt3, kr1, kr2, kr3, mass,
                           ge]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def all_materials(self) -> list[Any]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[Any]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def mass(self) -> np.ndarray:
        mass = self._mass
        return mass
