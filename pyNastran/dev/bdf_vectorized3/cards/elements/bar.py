from __future__ import annotations
from copy import deepcopy
from itertools import count
from typing import Union, Optional, Any, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16 # , print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string, integer_or_double,
    integer_or_blank, double_or_blank, string_or_blank,
    integer_double_or_blank, integer_string_or_blank,
    blank)
from pyNastran.bdf.bdf_interface.assign_type_force import force_double_or_blank, lax_double_or_blank
from pyNastran.bdf.cards.elements.bars import set_blank_if_default # init_x_g0,
from pyNastran.bdf.cards.properties.bars import (
    _bar_areaL, to_fields, get_beam_sections, parse_pbrsect_options)
# PBARL as pbarl, A_I1_I2_I12

from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, make_idim, hslice_by_idim,
    searchsorted_filter, parse_element_check, parse_property_check)
from pyNastran.dev.bdf_vectorized3.cards.elements.rod import (
    line_mid_mass_per_length, line_length, line_vector_length, line_centroid,
    e_g_nu_from_property_id,
)
from pyNastran.dev.bdf_vectorized3.cards.elements.utils import (
    get_density_from_material, basic_mass_material_id)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_int, array_default_str, array_default_float,
    array_default_float_nan, get_print_card_size)
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.cards.elements.beam import CBEAM, BEAMOR
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from ..coord import COORD
    #from pyNastran.dev.bdf_vectorized3.cards.elements.beam import BEAMOR, CBEAM


class BAROR(BaseCard):
    """
    +-------+---+-----+---+---+-------+-----+-------+------+
    |   1   | 2 |  3  | 4 | 5 |   6   |  7  |   8   |  9   |
    +=======+===+=====+===+===+=======+=====+=======+======+
    | BAROR |   | PID |   |   | G0/X1 |  X2 |  X3   | OFFT |
    +-------+---+-----+---+---+-------+-----+-------+------+
    | BAROR |   | 39  |   |   |  0.6  | 2.9 | -5.87 | GOG  |
    +-------+---+-----+---+---+-------+-----+-------+------+

    """
    type = 'BAROR'

    #@classmethod
    #def _init_from_empty(cls):
        #pid = 1
        #is_g0 = True
        #g0 = 1
        #x = None
        #return BAROR(pid, is_g0, g0, x, offt='GGG', comment='')

    def __init__(self, pid: int=0, g0: int=0, x=None, offt: str='GGG', comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        if x is None:
            x = np.array([0., 0., 0.])
        self.pid = pid
        self.g0 = g0
        self.x = x
        self.offt = offt
        #if isinstance(offt, integer_types):
            #raise NotImplementedError('the integer form of offt is not supported; offt=%s' % offt)

    def __deepcopy__(self, memo: dict[str, Any]):
        copy = type(self)()
        memo[id(self)] = copy
        if self.comment:
            copy.comment = self.comment
        copy.pid = deepcopy(self.pid, memo)
        copy.g0 = deepcopy(self.g0, memo)
        copy.x = deepcopy(self.x, memo)
        copy.offt = deepcopy(self.offt, memo)
        return copy

    @classmethod
    def add_card(cls, model: BDF, card: BDFCard, comment: str=''):
        PROPERTY_ID_DEFAULT = 0
        GO_X_DEFAULT = 0
        OFFT_DEFAULT = ''

        pid = integer_or_blank(card, 2, 'pid', default=PROPERTY_ID_DEFAULT)
        fdouble_or_blank = force_double_or_blank if model.is_lax_parser else double_or_blank
        # x / g0
        field5 = integer_double_or_blank(card, 5, 'g0_x1', default=0.)
        if isinstance(field5, integer_types):
            g0 = field5
            x = [np.nan, np.nan, np.nan]
            blank(card, 6, 'x2')
            blank(card, 7, 'x3')
        elif isinstance(field5, float):
            g0 = GO_X_DEFAULT
            x = np.array([field5,
                          fdouble_or_blank(card, 6, 'x2', default=0.),
                          fdouble_or_blank(card, 7, 'x3', default=0.)],
                          dtype='float64')
        else:
            raise NotImplementedError('BAROR field5 = %r' % field5)
        offt = integer_string_or_blank(card, 8, 'offt', default=OFFT_DEFAULT)
        assert len(card) <= 9, f'len(BAROR card) = {len(card):d}\ncard={card}'
        return BAROR(pid, g0, x, offt=offt, comment=comment)

    def raw_fields(self):
        """
        Gets the fields of the card in their full form
        """
        list_fields = ['BAROR', None, None] + self.x.tolist() + [self.offt]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


def init_x_g0(card: BDFCard, eid: int):
    """common method to read the x/g0 field for the CBAR, CBEAM, CBEAM3"""
    field5 = integer_double_or_blank(card, 5, 'g0_x1', default=np.nan)

    #GO_DEFAULT = -100
    if isinstance(field5, integer_types):
        g0 = field5
        x = [np.nan, np.nan, np.nan]
    elif isinstance(field5, float):
        g0 = 0
        x = np.array([field5,
                      double_or_blank(card, 6, 'x2', default=np.nan),
                      double_or_blank(card, 7, 'x3', default=np.nan)], dtype='float64')
        #if norm(x) == 0.0:
            #msg = 'G0 vector defining plane 1 is not defined.\n'
            #msg += 'G0 = %s\n' % g0
            #msg += 'X  = %s\n' % x
            #raise RuntimeError(msg)
    else:
        msg = ('field5 on %s (G0/X1) is the wrong type...id=%s field5=%s '
               'type=%s' % (card.field(0), eid, field5, type(field5)))
        raise RuntimeError(msg)
    return x, g0

def split_offt_vector(offt: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    neids = len(offt)
    offt_vector = np.full(neids, '', dtype='|U1')
    offt_end_a = np.full(neids, '', dtype='|U1')
    offt_end_b = np.full(neids, '', dtype='|U1')
    for i, (offt_vectori, offt_end_ai, offt_end_bi) in enumerate(offt):
        offt_vector[i] = offt_vectori
        offt_end_a[i] = offt_end_ai
        offt_end_b[i] = offt_end_bi
    return offt_vector, offt_end_a, offt_end_b

def get_bar_vector(elem, xyz1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    #if self.g0:
        #v = xyz[self.g0] - xyz[self.Ga()]
    #else:
        #v = self.x

    # get the vector v, which defines the projection on to the elemental
    # coordinate frame
    is_g0 = elem.is_g0
    is_x = elem.is_x
    v = np.full(elem.x.shape, np.nan, dtype=elem.x.dtype)
    if np.any(is_g0):
        grid_g0 = elem.model.grid.slice_card_by_node_id(elem.g0[is_g0], sort_ids=False)
        n0 = grid_g0.xyz_cid0()
        v[is_g0, :] = n0 - xyz1[is_g0, :]

    grid = elem.model.grid
    # get the cd frames for the nodes

    #n1 = self.nodes[:, 0]
    #in1 = np.searchsorted(grid.node_id, n1)
    #assert np.array_equal(grid.node_id[in1], n1)
    #cd = grid.cd[in1]

    # get the cd frames for nodes A/B (1/2)
    inode = np.searchsorted(grid.node_id, elem.nodes)
    assert np.array_equal(grid.node_id[inode], elem.nodes)
    cd = grid.cd[inode]
    cd1 = cd[:, 0]

    if np.any(is_x):
        coords = elem.model.coord
        is_cd1 = (cd1 > 0) & (is_x)

        v[is_x, :] = elem.x[is_x, :]
        if np.any(is_cd1):
            #cd1_ref = coords.slice_card_by_id(cd1, sort_ids=False)
            #print(cd1_ref.get_stats())
            #print(cd1_ref.get_object_methods())
            #v[is_cd1, :] = cd1_ref.transform_local_xyz_to_global(self.x)
            x_vector = elem.x[is_cd1, :]
            cd1i = cd1[is_cd1]

            #vi = coords.transform_local_xyz_to_global_coords(x_vector, cd1i)  # old
            # transform_node_to_global
            vi = coords.transform_offset_xyz_to_global_xyz(x_vector, cd1i)
            v[is_cd1, :] = vi

            #v[is_cd1, :] = coords.transform_local_xyz_to_global(x_vector, cd1[is_cd1])
            #cd1_ref = model.Coord(cd1)
            #cd2_ref = model.Coord(cd2)
    vmax = v.max(axis=1)
    inan = np.isnan(vmax)
    if np.any(inan):
        element_id = elem.element_id[inan]
        elem.model.log.error(f'{elem.type} element_id={element_id} has nan vectors')
    return v, cd


class CBAR(Element):
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='int32')
        self.offt = np.array([], dtype='|U3')
        self.g0 = np.array([], dtype='int32')
        self.x = np.full((0, 3), np.nan, dtype='float64')

        self.pa = np.array([], dtype='int32')
        self.pb = np.array([], dtype='int32')
        self.wa = np.zeros((0, 3), dtype='float64')
        self.wb = np.zeros((0, 3), dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int],
            x: Optional[list[float]], g0: Optional[int],
            offt: str='GGG', pa: int=0, pb: int=0,
            wa: Optional[list[float]]=None, wb: Optional[list[float]]=None,
            comment: str='', validate: bool=False) -> int:
        assert x is None or g0 is None, f'pid={pid} x={x} g0={g0}'
        self.cards.append((eid, pid, nids, x, g0, offt, pa, pb, wa, wb, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        PROPERTY_ID_DEFAULT = 0
        OFFT_DEFAULT = ''

        eid = integer(card, 1, 'eid')

        #pid_default = eid
        #if baror is not None:
            #if baror.pid is not None:
                #pid_default = baror.pid
            #if baror.x is None:
                #x1_default = baror.g0
                #x2_default = None
                #x3_default = None
            #else:
                #x1_default, x2_default, x3_default = baror.x
            #offt_default = baror.offt

        pid = integer_or_blank(card, 2, 'pid', default=PROPERTY_ID_DEFAULT)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        x, g0 = init_x_g0(card, eid)

        # doesn't exist in NX nastran
        offt = integer_string_or_blank(card, 8, 'offt', default=OFFT_DEFAULT)
        #print('cls.offt = %r' % (cls.offt))

        pa = integer_or_blank(card, 9, 'pa', default=0)
        pb = integer_or_blank(card, 10, 'pb', default=0)

        wa = np.array([double_or_blank(card, 11, 'w1a', 0.0),
                       double_or_blank(card, 12, 'w2a', 0.0),
                       double_or_blank(card, 13, 'w3a', 0.0)], dtype='float64')

        wb = np.array([double_or_blank(card, 14, 'w1b', 0.0),
                       double_or_blank(card, 15, 'w2b', 0.0),
                       double_or_blank(card, 16, 'w3b', 0.0)], dtype='float64')
        assert len(card) <= 17, f'len(CBAR card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, [ga, gb], x, g0, offt, pa, pb, wa, wb, comment))
        self.n += 1
        return self.n - 1

    def __apply_slice__(self, elem: CBAR, i: np.ndarray) -> None:  # ignore[override]
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.offt = self.offt[i]
        elem.g0 = self.g0[i]
        elem.x = self.x[i, :]
        elem.pa = self.pa[i]
        elem.pb = self.pb[i]
        elem.wa = self.wa[i, :]
        elem.wb = self.wb[i, :]
        elem.n = len(i)

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 2), dtype=idtype)
        offt = np.full(ncards, '', dtype='|U3')
        g0 = np.zeros(ncards, dtype=idtype)
        x = np.full((ncards, 3), np.nan, dtype='float64')

        pa = np.zeros(ncards, dtype='int32')
        pb = np.zeros(ncards, dtype='int32')
        wa = np.zeros((ncards, 3), dtype='float64')
        wb = np.zeros((ncards, 3), dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, xi, g0i, offti, pai, pbi, wai, wbi, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids

            if g0i in {None, 0, -1}:
                g0i = 0
            else:
                assert g0i > 0, g0i
                xi = [np.nan, np.nan, np.nan]
            g0[icard] = g0i
            x[icard, :] = xi
            offt[icard] = offti

            pa[icard] = pai
            pb[icard] = pbi
            if wai is not None:
                wa[icard, :] = wai
            if wbi is not None:
                wb[icard, :] = wbi
        self._save(element_id, property_id, nodes, g0, x,
                   offt, pa, pb, wa, wb)
        baror = self.model.baror
        apply_bar_default(self, baror)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, nodes, g0, x,
              offt, pa, pb, wa, wb) -> None:
        #assert len(element_id) == len(property_id), 'A1'
        #assert len(self.element_id) == len(self.property_id), 'A2'
        if len(self.element_id):
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            nodes = np.vstack([self.nodes, nodes])
            offt = np.hstack([self.offt, offt])
            g0 = np.hstack([self.g0, g0])
            x = np.vstack([self.x, x])

            pa = np.hstack([self.pa, pa])
            pb = np.hstack([self.pb, pb])
            wa = np.vstack([self.wa, wa])
            wb = np.vstack([self.wb, wb])

        assert len(element_id) == len(property_id), 'B'
        self.element_id = element_id
        self.property_id = property_id
        #print(element_id)
        #print(property_id)
        self.nodes = nodes
        self.offt = offt
        self.g0 = g0
        self.x = x

        self.pa = pa
        self.pb = pb
        self.wa = wa
        self.wb = wb

    def convert(self, xyz_scale: float=1.0,
                mass_scale: float=1.0, **kwargs):
        ## TODO: probably wrong for CD=1
        self.x *= xyz_scale
        self.wa *= xyz_scale
        self.wb *= xyz_scale

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['element_id'].append(self.element_id)
        used_dict['property_id'].append(self.property_id)
        used_dict['node_id'].append(self.nodes.ravel())
        g0 = self.g0[self.is_g0]
        if len(g0):
            used_dict['node_id'].append(g0)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(),
                   self.nodes.max(), self.g0.max())

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)
        pas = array_default_int(self.pa, default=0, size=size)
        pbs = array_default_int(self.pb, default=0, size=size)
        was = array_default_float(self.wa, default=0.0, size=size, is_double=False)
        wbs = array_default_float(self.wb, default=0.0, size=size, is_double=False)
        for eid, pid, nodes, g0, x, offt, pa, pb, wa, wb in zip(element_ids, property_ids, nodes,
                                                                self.g0, self.x, self.offt, pas, pbs, was, wbs):
            n1, n2 = nodes
            w1a, w2a, w3a = wa
            #w1a = set_blank_if_default(wa[0], 0.0)
            #w2a = set_blank_if_default(wa[1], 0.0)
            #w3a = set_blank_if_default(wa[2], 0.0)

            w1b, w2b, w3b = wb
            #w1b = set_blank_if_default(wb[0], 0.0)
            #w2b = set_blank_if_default(wb[1], 0.0)
            #w3b = set_blank_if_default(wb[2], 0.0)
            if g0 in {-1, 0}:
                x1, x2, x3 = x # self.get_x_g0_defaults()
            else:
                x1 = g0
                x2 = ''
                x3 = ''

            # offt doesn't exist in NX nastran
            offt = set_blank_if_default(offt, 'GGG')

            list_fields = ['CBAR', eid, pid, n1, n2,
                           x1, x2, x3, offt, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
            #print(list_fields)
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pbar, self.model.pbarl]
                if prop.n > 0]

    def mass_material_id(self) -> np.ndarray:
        material_id = basic_mass_material_id(self.property_id, self.allowed_properties, 'CBAR')
        return material_id

    def mass(self) -> np.ndarray:
        pid = self.property_id
        mass_per_length = np.full(len(pid), np.nan, dtype='float64')
        for prop in self.allowed_properties:
            i_lookup, i_all = searchsorted_filter(prop.property_id, pid, msg='')
            if len(i_lookup) == 0:
                continue

            # we're at least using some properties
            mass_per_lengthi = prop.mass_per_length()
            mass_per_length[i_lookup] = mass_per_lengthi[i_all]
        length = self.length()
        mass = mass_per_length * length
        return mass

    def mass_breakdown(self) -> np.ndarray:
        """
        [L, rho, A, nsm, mpl, mass]
        """
        pid = self.property_id
        rho = np.full(len(pid), np.nan, dtype='float64')
        area = np.full(len(pid), np.nan, dtype='float64')
        nsm = np.full(len(pid), np.nan, dtype='float64')
        mass_per_length = np.full(len(pid), np.nan, dtype='float64')
        for prop in self.allowed_properties:
            i_lookup, i_all = searchsorted_filter(prop.property_id, pid, msg='')
            if len(i_lookup) == 0:
                continue

            # we're at least using some properties
            breakdowni = prop.mass_per_length_breakdown() # [rho, A, nsm, mpl]

            rho[i_lookup] = breakdowni[i_all, 0]
            area[i_lookup] = breakdowni[i_all, 1]
            nsm[i_lookup] = breakdowni[i_all, 2]
            mass_per_length[i_lookup] = breakdowni[i_all, 3]
        length = self.length()
        mass = mass_per_length * length
        breakdown = np.column_stack([length, rho, area, nsm, mass_per_length, mass])
        return breakdown

    def area(self) -> np.ndarray:
        pid = self.property_id
        area = np.full(len(pid), np.nan, dtype='float64')
        for prop in self.allowed_properties:
            i_lookup, i_all = searchsorted_filter(prop.property_id, pid, msg='')
            if len(i_lookup) == 0:
                continue

            # we're at least using some properties
            areai = prop.area()
            area[i_lookup] = areai[i_all]
        return area

    def inertia(self) -> np.ndarray:
        inertia = inertia_from_property_id(self.property_id,
                                           self.allowed_properties)
        return inertia

    def k(self) -> np.ndarray:
        k1_k2 = k_from_property_id(self.property_id,
                                   self.allowed_properties)
        return k1_k2

    def e_g_nu(self) -> np.ndarray:
        e_g_nu = e_g_nu_from_property_id(self.property_id, self.allowed_properties)
        return e_g_nu

    def line_vector_length(self) -> tuple[np.ndarray, np.ndarray]:
        line_vector, length = line_vector_length(self.model, self.nodes)
        return line_vector, length

    def length(self) -> np.ndarray:
        length = line_length(self.model, self.nodes)
        return length

    def centroid(self) -> np.ndarray:
        centroid = line_centroid(self.model, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def volume(self) -> np.ndarray:
        A = self.area()
        L = self.length()
        return A * L

    @property
    def is_x(self) -> np.ndarray:
        return (self.g0 == 0)

    @property
    def is_g0(self) -> np.ndarray:
        return ~self.is_x

    def get_xyz(self) -> tuple[np.ndarray, np.ndarray]:
        #neids = len(self.element_id)
        grid = self.model.grid
        xyz = grid.xyz_cid0()
        nid = grid.node_id
        inode = np.searchsorted(nid, self.nodes)
        assert np.array_equal(nid[inode], self.nodes)
        in1 = inode[:, 0]
        in2 = inode[:, 1]
        xyz1 = xyz[in1, :]
        xyz2 = xyz[in2, :]
        return xyz1, xyz2

    def get_bar_vector(self, xyz1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        v, cd = get_bar_vector(self, xyz1)
        return v, cd

    def get_axes(self, xyz1: np.ndarray, xyz2: np.ndarray,
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                            np.ndarray, np.ndarray, np.ndarray]:
        log = self.model.log
        coords: COORD = self.model.coord
        #xyz1, xyz2 = self.get_xyz()

        neids = xyz1.shape[0]
        i = xyz2 - xyz1
        ihat_norm = np.linalg.norm(i, axis=1)
        assert len(ihat_norm) == neids
        if min(ihat_norm) == 0.:
            msg = 'xyz1=%s xyz2=%s\n%s' % (xyz1, xyz2, self)
            log.error(msg)
            raise ValueError(msg)
        i_offset = i / ihat_norm[:, np.newaxis]

        #log.info(f'x =\n{self.x}')
        #log.info(f'g0   = {self.g0}')
        v, cd = self.get_bar_vector(xyz1)
        cd1 = cd[:, 0]
        cd2 = cd[:, 1]
        vnorm1 = np.linalg.norm(v, axis=1)
        if np.any(vnorm1 == 0.):
            raise RuntimeError(f'vnorm1 = {vnorm1}')

        offt_vector, offt_end_a, offt_end_b = split_offt_vector(self.offt)
        is_rotate_v_g = (offt_vector == 'G')
        is_rotate_wa_g = (offt_end_a == 'G')
        is_rotate_wb_g = (offt_end_b == 'G')

        #is_rotate_v_b = (offt_vector == 'B')
        #is_rotate_wa_b = (offt_end_a == 'B')
        #is_rotate_wb_b = (offt_end_b == 'B')

        is_rotate_wa_o = (offt_end_a == 'O')
        is_rotate_wb_o = (offt_end_b == 'O')

        uofft_vector = np.unique(offt_vector)
        uofft_end_a = np.unique(offt_end_a)
        uofft_end_b = np.unique(offt_end_b)

        msg = ''
        for i, offt_vectori in enumerate(uofft_vector):
            if offt_vectori not in 'GB':
                msg += f'OFFT field[0]={offt_vectori} and must be G/B; offt={self.offt[i]}\n'
        for i, offt_end_ai in enumerate(uofft_end_a):
            if offt_end_ai not in 'GBO':
                msg += f'OFFT field[1]={offt_end_ai} and must be G/B/O; offt={self.offt[i]}\n'
        for i, offt_end_bi in enumerate(uofft_end_b):
            if offt_end_bi not in 'GBO':
                msg += f'OFFT field[2]={offt_end_bi} and must be G/B/O; offt={self.offt[i]}\n'
        if msg:
            log.error(msg)
            raise ValueError(msg)

        #--------------------------------------------------------------------------
        # rotate v
        #log.info(f'offt = {self.offt}')
        #log.info(f'v0 =\n{v}')
        #log.info(f'cd =\n{cd}')

        if np.any(is_rotate_v_g):
            # end A
            # global - cid != 0
            icd1_v_vector = (is_rotate_v_g) & (cd1 != 0)
            cd1_v_vector = cd1[icd1_v_vector]
            eid_v_vector = self.element_id[icd1_v_vector]
            if np.any(cd1_v_vector):
                cd1_ref: COORD = coords.slice_card_by_id(cd1_v_vector, sort_ids=False)
                v1v = v[icd1_v_vector, :]
                v2 = cd1_ref.transform_xyz_to_global_assuming_rectangular(v1v)
                v[icd1_v_vector, :] = v2
                del v1v, v2
            del icd1_v_vector, cd1_v_vector, eid_v_vector

        #elif offt_vector == 'B':
            # basic - cid = 0
            #pass

        if np.any(np.isnan(v.max(axis=1))):
            raise RuntimeError(f'v = {v}')

        vnorm = np.linalg.norm(v, axis=1)
        if np.any(vnorm == 0.):
            raise RuntimeError(f'vnorm = {vnorm}')

        #--------------------------------------------------------------------------
        # determine the bar vectors
        #log.info(f'v =\n{v}')
        #log.info(f'ihat =\n{i_offset}')
        ihat = i_offset
        #vhat = safe_normalize(v)

        vnorm = np.linalg.norm(v, axis=1)
        vhat = np.full(v.shape, np.nan, dtype=v.dtype)
        izero = (vnorm > 0)
        if np.any(izero):
            vhat[izero] = v[izero, :] / vnorm[izero, np.newaxis]
            #vhat = v / vnorm[:, np.newaxis] # j
        z = np.cross(ihat, vhat) # k
        norm_z = np.linalg.norm(z, axis=1)
        assert len(norm_z) == neids

        #if np.any(np.isnan(zhat.max(axis=1))):
        #print(f'norm_z = {norm_z}')

        zhat = safe_normalize(z)

        #zhat = z / norm_z[:, np.newaxis]
        yhat = np.cross(zhat, ihat) # j
        norm_i = np.linalg.norm(ihat, axis=1)
        norm_yhat = np.linalg.norm(yhat, axis=1)
        xform_offset = np.dstack([ihat, yhat, zhat]) # 3x3 unit matrix

        #if np.any(np.isnan(yhat.max(axis=1))):
            #print(f'norm_yhat = {norm_yhat}')

        del norm_i, norm_z, norm_yhat
        #--------------------------------------------------------------------------
        # rotate wa
        # wa defines the offset at end A
        wa = self.wa.copy()  # we're going to be inplace hacking it, so copy :)
        assert not np.isnan(np.max(wa)), wa

        if np.any(is_rotate_wa_g):
            icd1_vector = (is_rotate_wa_g) & (cd1 != 0)
            cd1_vector = cd1[icd1_vector]
            if np.any(icd1_vector):
                cd1_ref = coords.slice_card_by_id(cd1_vector, sort_ids=False)
                wai1 = wa[icd1_vector, :]
                wai2 = cd1_ref.transform_xyz_to_global_assuming_rectangular(wai1)
                #print('eids.shape =', self.element_id.shape)
                #print('len(cd1_vector) =', len(cd1_vector))
                #print('icd1_vector.shape =', icd1_vector.shape)
                #print('is_rotate_wa.shape =', is_rotate_wa.shape)
                #print('wai1.shape =', wai1.shape)
                #print('wai2.shape =', wai2.shape)
                #print('wa.shape =', wa.shape)
                wa[icd1_vector, :] = wai2
            del cd1_vector, icd1_vector
        #elif offt_end_a == 'B':
            #pass
        if np.any(is_rotate_wa_o):
            # rotate point wa from the local frame to the global frame
            #wa = wa @ xform_offset
            wao1 = wa[is_rotate_wa_o, :]
            To = xform_offset[is_rotate_wa_o, :, :]
            wao = np.einsum('ni,nij->nj', wao1, To)
            wa[is_rotate_wa_o, :] = wao
            del wao1, To, wao

        assert not np.isnan(np.max(wa)), wa

        #--------------------------------------------------------------------------
        # rotate wb
        # wb defines the offset at end B
        wb = self.wb.copy()  # we're going to be inplace hacking it, so copy :)
        if np.any(is_rotate_wb_g):
            icd2_vector = (is_rotate_wb_g) & (cd2 != 0)
            cd2_vector = cd2[icd2_vector]
            #cd2_vector = cd2[is_rotate_wb]
            #icd2_vector = (cd2_vector != 0)
            if np.any(icd2_vector):
                # MasterModelTaxi
                #wb = cd2_ref.transform_node_to_global_assuming_rectangular(wb)
                cd2_ref = coords.slice_card_by_id(cd2_vector, sort_ids=False)
                wbi1 = wb[icd2_vector, :]
                wbi2 = cd2_ref.transform_xyz_to_global_assuming_rectangular(wbi1)
                wb[icd2_vector, :] = wbi2
            del cd2_vector, icd2_vector
        #elif offt_end_b == 'B':
            #pass

        if np.any(is_rotate_wb_o):
            # rotate point wb from the local frame to the global frame

            wbo1 = wb[is_rotate_wb_o, :]
            To = xform_offset[is_rotate_wb_o, :, :]
            wbo = np.einsum('ni,nij->nj', wbo1, To)
            wb[is_rotate_wb_o, :] = wbo
            del wbo1, To, wbo
            #wb = wb @ xform_offset
            #ib = n2 + wb

        assert not np.isnan(np.max(wb)), wb

        #ihat = xform[0, :]
        #yhat = xform[1, :]
        #zhat = xform[2, :]
        #wa, wb, _ihat, jhat, khat = out

        # we finally have the nodal coordaintes!!!! :)
        return v, ihat, yhat, zhat, wa, wb

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no bar properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

def safe_normalize(vector: np.ndarray, axis: int=1) -> np.ndarray:
    norm_vector = np.linalg.norm(vector, axis=1)
    vector_hat = np.full(vector.shape, np.nan, dtype=vector.dtype)
    izero = (norm_vector > 0)
    if np.any(izero):
        vector_hat[izero] = vector[izero, :] / norm_vector[izero, np.newaxis]
    return vector_hat

PBARL_MSG = '\n' + """
+-------+------+------+-------+------+------+------+------+------+
|   1   |   2  |   3  |   4   |   5  |   6  |   7  |   8  |   9  |
+=======+======+======+=======+======+======+======+======+======+
| PBARL | PID  | MID  | GROUP | TYPE |      |      |      |      |
+-------+------+------+-------+------+------+------+------+------+
|       | DIM1 | DIM2 | DIM3  | DIM4 | DIM5 | DIM6 | DIM7 | DIM8 |
+-------+------+------+-------+------+------+------+------+------+
|       | DIM9 | etc. |  NSM  |      |      |      |      |      |
+-------+------+------+-------+------+------+------+------+------+""".strip()

class PBAR(Property):
    """
    Defines the properties of a simple beam element (CBAR entry).

    +------+-----+-----+-----+----+----+----+-----+-----+
    |   1  |  2  |  3  |  4  |  5 |  6 |  7 |  8  |  9  |
    +======+=====+=====+=====+====+====+====+=====+=====+
    | PBAR | PID | MID |  A  | I1 | I2 | J  | NSM |     |
    +------+-----+-----+-----+----+----+----+-----+-----+
    |      | C1  | C2  | D1  | D2 | E1 | E2 | F1  | F2  |
    +------+-----+-----+-----+----+----+----+-----+-----+
    |      | K1  | K2  | I12 |    |    |    |     |     |
    +------+-----+-----+-----+----+----+----+-----+-----+

    .. todo::
        support solution 600 default
        do a check for mid -> MAT1      for structural
        do a check for mid -> MAT4/MAT5 for thermal
    """
    _show_attributes = ['property_id', 'material_id', 'A', 'J',
                        'c', 'd', 'e', 'f', 'I', 'k', 'nsm']
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')

        self.A = np.array([], dtype='float64')
        self.J = np.array([], dtype='float64')

        self.c = np.zeros((0, 2), dtype='float64')
        self.d = np.zeros((0, 2), dtype='float64')
        self.e = np.zeros((0, 2), dtype='float64')
        self.f = np.zeros((0, 2), dtype='float64')

        self.I = np.zeros((0, 3), dtype='float64')
        self.k = np.zeros((0, 2), dtype='float64')
        self.nsm = np.array([], dtype='float64')

    def add(self, pid: int, mid: int, A: float=0.,
                i1: float=0., i2: float=0., i12: float=0., j: float=0.,
                nsm: float=0.,
                c1: float=0., c2: float=0.,
                d1: float=0., d2: float=0.,
                e1: float=0., e2: float=0.,
                f1: float=0., f2: float=0.,
                k1: float=1.e8, k2: float=1.e8, comment: str='') -> int:
        self.cards.append((pid, mid, A, i1, i2, i12, j, nsm, c1, c2, d1, d2,
                           e1, e2, f1, f2, k1, k2, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        A = double_or_blank(card, 3, 'A', default=0.0)
        i1 = double_or_blank(card, 4, 'I1', default=0.0)
        i2 = double_or_blank(card, 5, 'I2', default=0.0)

        j = double_or_blank(card, 6, 'J', default=0.0)
        nsm = double_or_blank(card, 7, 'nsm', default=0.0)

        c1 = double_or_blank(card, 9, 'C1', default=0.0)
        c2 = double_or_blank(card, 10, 'C2', default=0.0)
        d1 = double_or_blank(card, 11, 'D1', default=0.0)
        d2 = double_or_blank(card, 12, 'D2', default=0.0)
        e1 = double_or_blank(card, 13, 'E1', default=0.0)
        e2 = double_or_blank(card, 14, 'E2', default=0.0)
        f1 = double_or_blank(card, 15, 'F1', default=0.0)
        f2 = double_or_blank(card, 16, 'F2', default=0.0)

        i12 = double_or_blank(card, 19, 'I12', default=0.0)

        if A == 0.0:
            blank(card, 17, 'K1')
            blank(card, 18, 'K2')
            k1 = np.nan
            k2 = np.nan
        elif i12 != 0.0:
            # K1 / K2 are ignored
            k1 = np.nan
            k2 = np.nan
        else:
            #: default=infinite; assume 1e8
            k1 = double_or_blank(card, 17, 'K1', 1e8)
            #: default=infinite; assume 1e8
            k2 = double_or_blank(card, 18, 'K2', 1e8)

        assert len(card) <= 20, f'len(PBAR card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, mid, A, i1, i2, i12, j, nsm, c1, c2, d1, d2,
                           e1, e2, f1, f2, k1, k2, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        #Type = np.full(ncards, '', dtype='|U8')
        #group = np.full(ncards, '', dtype='|U8')

        #material_id[icard] = mid
        #group[icard] = group
        A = np.zeros(ncards, dtype='float64')
        J = np.zeros(ncards, dtype='float64')

        c = np.zeros((ncards, 2), dtype='float64')
        d = np.zeros((ncards, 2), dtype='float64')
        e = np.zeros((ncards, 2), dtype='float64')
        f = np.zeros((ncards, 2), dtype='float64')

        I = np.zeros((ncards, 3), dtype='float64')
        k = np.zeros((ncards, 2), dtype='float64')
        nsm = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, mid, Ai, i1, i2, i12, j, nsmi, c1, c2, d1, d2,
             e1, e2, f1, f2, k1, k2, comment) = card
            property_id[icard] = pid
            material_id[icard] = mid
            #group[icard] = group
            A[icard] = Ai
            I[icard, :] = [i1, i2, i12]
            J[icard] = j

            c[icard, :] = [c1, c2]
            d[icard, :] = [d1, d2]
            e[icard, :] = [e1, e2]
            f[icard, :] = [f1, f2]

            k[icard, :] = [k1, k2]
            nsm[icard] = nsmi
        self._save(property_id, material_id, A, J, c, d, e, f, I, k, nsm)
        self.sort()
        self.cards = []

    def _save(self, property_id, material_id, A, J, c, d, e, f, I, k, nsm):
        if len(self.property_id):
            property_id = np.hstack([self.property_id, property_id])
            material_id = np.hstack([self.property_id, material_id])
            A = np.hstack([self.A, A])
            J = np.hstack([self.J, J])

            c = np.vstack([self.c, c])
            d = np.vstack([self.d, d])
            e = np.vstack([self.e, e])
            f = np.vstack([self.f, f])
            I = np.vstack([self.I, I])
            k = np.vstack([self.k, k])
            nsm = np.hstack([self.nsm, nsm])

        self.property_id = property_id
        self.material_id = material_id
        #Type = np.full(ncards, '', dtype='|U8')
        #group = np.full(ncards, '', dtype='|U8')

        #material_id[icard] = mid
        #group[icard] = group
        self.A = A
        self.I = I
        self.J = J

        self.c = c
        self.d = d
        self.e = e
        self.f = f

        self.k = k
        self.nsm = nsm

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    def convert(self, xyz_scale: float=1.0,
                area_scale: float=1.0,
                area_inertia_scale: float=1.0,
                nsm_per_length_scale: float=1.0, **kwargs):
        self.A *= area_scale
        self.I *= area_inertia_scale
        self.J *= area_inertia_scale

        self.c *= xyz_scale
        self.d *= xyz_scale
        self.e *= xyz_scale
        self.f *= xyz_scale
        self.nsm *= nsm_per_length_scale

        ## TODO: probably wrong
        ##self.k

    def __apply_slice__(self, prop: PBAR, i: np.ndarray) -> None:  # ignore[override]
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        #self.Type = np.full(ncards, '', dtype='|U8')

        #prop.Type = Type[i]
        #prop.group = self.group[i]
        prop.A = self.A[i]
        prop.I = self.I[i, :]
        prop.J = self.J[i]

        prop.c = self.c[i, :]
        prop.d = self.d[i, :]
        prop.e = self.e[i, :]
        prop.f = self.f[i, :]

        prop.k = self.k[i, :]
        prop.nsm = self.nsm[i]
        prop.n = len(i)

    def validate(self) -> None:
        if np.any(self.i1 < 0.):
            raise ValueError('I1=%r must be greater than or equal to 0.0' % self.i1)
        if np.any(self.i2 < 0.):
            raise ValueError('I2=%r must be greater than or equal to 0.0' % self.i2)
        if np.any(self.j < 0.):
            raise ValueError('J=%r must be greater than or equal to 0.0' % self.j)

    @property
    def i1(self) -> np.ndarray:
        return self.I[:, 0]
    @property
    def i2(self) -> np.ndarray:
        return self.I[:, 1]
    @property
    def i12(self) -> np.ndarray:
        return self.I[:, 2]
    @property
    def j(self) -> np.ndarray:
        return self.J
    @property
    def k1(self) -> np.ndarray:
        return self.k[:, 0]
    @property
    def k2(self) -> np.ndarray:
        return self.k[:, 1]

    @property
    def c1(self) -> np.ndarray:
        return self.c[:, 0]
    @property
    def c2(self) -> np.ndarray:
        return self.c[:, 1]
    @property
    def d1(self) -> np.ndarray:
        return self.d[:, 0]
    @property
    def d2(self) -> np.ndarray:
        return self.d[:, 1]
    @property
    def e1(self) -> np.ndarray:
        return self.e[:, 0]
    @property
    def e2(self) -> np.ndarray:
        return self.e[:, 1]
    @property
    def f1(self) -> np.ndarray:
        return self.f[:, 0]
    @property
    def f2(self) -> np.ndarray:
        return self.f[:, 1]

    @i1.setter
    def i1(self, i1: np.ndarray) -> None:
        self.I[:, 0] = i1
    @i2.setter
    def i2(self, i2: np.ndarray) -> None:
        self.I[:, 1] = i2
    @i12.setter
    def i12(self, i12: np.ndarray) -> None:
        self.I[:, 2] = i12
    @j.setter
    def j(self, j: np.ndarray) -> None:
        self.J = j

    @k1.setter
    def k1(self, k1: np.ndarray) -> None:
        self.k[:, 0] = k1
    @k2.setter
    def k2(self, k2: np.ndarray) -> None:
        self.k[:, 1] = k2

    @c1.setter
    def c1(self, c1: np.ndarray) -> None:
        self.c[:, 0] = c1
    @c2.setter
    def c2(self, c2: np.ndarray) -> None:
        self.c[:, 1] = c2
    @d1.setter
    def d1(self, d1: np.ndarray) -> None:
        self.d[:, 0] = d1
    @d2.setter
    def d2(self, d2: np.ndarray) -> None:
        self.d[:, 1] = d2
    @e1.setter
    def e1(self, e1: np.ndarray) -> None:
        self.e[:, 0] = e1
    @e2.setter
    def e2(self, e2: np.ndarray) -> None:
        self.e[:, 1] = e2
    @f1.setter
    def f1(self, f1: np.ndarray) -> None:
        self.f[:, 0] = f1
    @f2.setter
    def f2(self, f2: np.ndarray) -> None:
        self.f[:, 1] = f2

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          msg=f'no bar materials for {self.type}')
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max())

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_ids = array_str(self.property_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        iis = array_default_float(self.I, default=0.0, size=size)
        areas = array_default_float(self.A, default=0.0, size=size)
        js = array_default_float(self.J, default=0.0, size=size)
        nsms = array_default_float(self.nsm, default=0.0, size=size)
        cs = array_default_float(self.c, default=0.0, size=size)
        ds = array_default_float(self.d, default=0.0, size=size)
        es = array_default_float(self.e, default=0.0, size=size)
        fs = array_default_float(self.f, default=0.0, size=size)
        ks = array_default_float_nan(self.k, default=1e8, size=size)
        for pid, mid, A, I, j, nsm, c, d, e, f, k in zip(property_ids, material_ids,
                                                         areas, iis, js,
                                                         nsms,
                                                         cs, ds, es, fs,
                                                         ks):
            c1, c2 = c
            d1, d2 = d
            e1, e2 = e
            f1, f2 = f
            i1, i2, i12 = I
            k1, k2 = k

            #i1 = set_blank_if_default(i1, 0.0)
            #i2 = set_blank_if_default(i2, 0.0)
            #i12 = set_blank_if_default(i12, 0.0)
            #j = set_blank_if_default(j, 0.0)
            #nsm = set_blank_if_default(nsm, 0.0)

            #c1 = set_blank_if_default(c1, 0.0)
            #c2 = set_blank_if_default(c2, 0.0)

            #d1 = set_blank_if_default(d1, 0.0)
            #d2 = set_blank_if_default(d2, 0.0)

            #e1 = set_blank_if_default(e1, 0.0)
            #e2 = set_blank_if_default(e2, 0.0)

            #f1 = set_blank_if_default(f1, 0.0)
            #f2 = set_blank_if_default(f2, 0.0)

            #k1 = set_blank_if_default(k1, 1e8)
            #k2 = set_blank_if_default(k2, 1e8)

            list_fields = ['PBAR', pid, mid, A, i1, i2, j, nsm,
                           None, c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, i12]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def all_materials(self) -> list[Any]:
        if self.model.is_thermal:
            materials = [self.model.mat4, self.model.mat5]
        else:
            materials = [self.model.mat1]
        return materials

    @property
    def allowed_materials(self) -> list[Any]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def mass_per_length(self) -> np.ndarray:
        mass_per_length = line_mid_mass_per_length(self.material_id, self.nsm, self.A,
                                                   self.allowed_materials)
        return mass_per_length

    def area(self) -> np.ndarray:
        return self.A

    def inertia(self) -> np.ndarray:
        inertia = np.column_stack([self.i1, self.i2, self.i12, self.j])
        return inertia

    def e_g_nu(self) -> np.ndarray:
        """calculates E, G, nu"""
        e_g_nu = e_g_nu_from_isotropic_material(self.material_id, self.allowed_materials)
        return e_g_nu

def e_g_nu_from_isotropic_material(material_id: np.ndarray,
                                   allowed_materials: list[Material]) -> np.ndarray:
    """calculates E, G, nu"""

    assert len(allowed_materials) > 0, allowed_materials
    nmaterials = len(material_id)
    if nmaterials == 0:
        raise RuntimeError(f'material_id={material_id}')

    material_id_check = np.zeros(nmaterials, dtype='int32')
    e_g_nu = np.full((nmaterials, 3), np.nan, dtype='float64')
    for mat in allowed_materials:
        mat_material_ids = mat.material_id

        i_lookup, i_all = searchsorted_filter(mat_material_ids, material_id)
        if len(i_all) == 0:
            continue

        material_id_check[i_lookup] = mat_material_ids[i_all]
        assert mat.type == 'MAT1', mat.type
        e_g_nu[i_lookup, 0] = mat.E[i_all]
        e_g_nu[i_lookup, 1] = mat.G[i_all]
        e_g_nu[i_lookup, 2] = mat.nu[i_all]

    assert len(e_g_nu) == nmaterials
    return e_g_nu


class PBARL(Property):
    valid_types = {
        "ROD": 1,
        "TUBE": 2,
        "TUBE2": 2,
        "I": 6,
        "CHAN": 4,
        "T": 4,
        "BOX": 4,
        "BAR": 2,
        "CROSS": 4,
        "H": 4,
        "T1": 4,
        "I1": 4,
        "CHAN1": 4,
        "Z": 4,
        "CHAN2": 4,
        "T2": 4,
        "BOX1": 6,
        "HEXA": 3,
        "HAT": 4,
        "HAT1": 5,
        "DBOX": 10,  # was 12

        # approximate
        #'I', 'CHAN', 'T', 'CHAN1', 'T1', 'CHAN2', 'T2', 'L' and 'BOX1'.
        'L' : 4,
    }  # for GROUP="MSCBML0"

    _skip_equality_check = True  # assume unequal
    _show_attributes = ['property_id', 'material_id', 'ndim', 'Type',
                        'group', 'nsm', 'dims']
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.ndim = np.array([], dtype='int32')
        self.Type = np.array([], dtype='|U8')
        self.group = np.array([], dtype='|U8')
        self.nsm = np.array([], dtype='float64')
        self.dims = np.array([], dtype='float64')

    #def slice_card_by_property_id(self, property_id: np.ndarray) -> PBARL:
        #"""uses a node_ids to extract PBARLs"""
        #iprop = self.index(property_id)
        #prop = self.slice_card_by_index(iprop)
        #return prop

    def validate(self) -> None:
        utypes = np.unique(self.Type)
        for utype in utypes:
            if utype not in self.valid_types:
                raise ValueError(f'PBARL Type={utype!r} is not valid')

        #try:
            #ndim = self.valid_types[self.Type]
        #except KeyError:
            #allowed = list(self.valid_types.keys())
            #msg = f'PBARL pid={self.pid}; Type={self.Type}; allowed={allowed}'
            #raise KeyError(msg)

        #print(self.valid_types)
        #print(self.Type)
        #print(self.property_id)
        ndim = [self.valid_types[Type] for Type in self.Type]
        #print(ndim)
        #if not isinstance(self.dim, list):
            #msg = 'PBARL pid=%s; dim must be a list; type=%r' % (self.pid, type(self.dim))
            #raise TypeError(msg)
        #if len(self.dim) != ndim:
            #msg = 'dim=%s len(dim)=%s Type=%s len(dimType)=%s' % (
                #self.dim, len(self.dim), self.Type,
                #self.valid_types[self.Type])
            #raise RuntimeError(msg)

        #assert len(self.dim) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(self.dim))
        #if not isinstance(self.group, str):
            #raise TypeError('Invalid group; pid=%s group=%r' % (self.pid, self.group))

    def __apply_slice__(self, prop: PBARL, i: np.ndarray) -> None:  # ignore[override]
        assert self.ndim.sum() == len(self.dims)
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        prop.Type = self.Type[i]
        prop.group = self.group[i]
        prop.nsm = self.nsm[i]

        idim = self.idim # [i, :]
        prop.dims = hslice_by_idim(i, idim, self.dims)

        prop.ndim = self.ndim[i]
        prop.n = len(i)
        assert isinstance(prop.ndim, np.ndarray), prop.ndim
        assert prop.ndim.sum() == len(prop.dims), f'prop.ndim={prop.ndim} len(prop.dims)={len(prop.dims)}'

    def add(self, pid: int, mid: int, bar_type: str, dim: list[float],
            group: str='MSCBML0', nsm: float=0., comment: str='') -> int:
        if isinstance(dim, integer_types):
            dim = [float(dim)]
        elif isinstance(dim, float_types):
            dim = [dim]

        if not isinstance(group, str):
            msg = f'PBARL: property_id={pid:d} group={group!r}'
            raise TypeError(msg)

        ndim = len(dim)
        assert ndim > 0, f'PBARL: property_id={pid} dims={dim}'

        ndim = self.valid_types[bar_type]
        if len(dim) != ndim:
            raise ValueError(f'PBARL pid={pid:d} bar_type={bar_type} ndim={ndim:d} len(dims)={dim}')

        self.cards.append((pid, mid, group, bar_type, dim, nsm, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        group = string_or_blank(card, 3, 'group', default='MSCBML0')
        bar_type = string(card, 4, 'Type')

        try:
            ndim = self.valid_types[bar_type]
        except KeyError:
            keys = list(self.valid_types.keys())
            raise KeyError('%r is not a valid PBARL type\nallowed_types={%s}' % (
                bar_type, ', '.join(sorted(keys))))

        # PBARL
        # 9. For DBOX section, the default value for DIM5 to DIM10 are
        #    based on the following rules:
        #     a. DIM5, DIM6, DIM7 and DIM8 have a default value of
        #        DIM4 if not provided.
        #     b. DIM9 and DIM10 have a default value of DIM6 if not
        #        provided.
        dim = []
        if bar_type == 'DBOX':
            for ioffset in range(ndim):
                if ioffset in {4, 5, 6, 7}:
                    dim4 = dim[3]
                    dimi = double_or_blank(card, 9 + ioffset, f'ndim={ndim}; dim{ioffset+1}',
                                           default=dim4, end=PBARL_MSG)
                elif ioffset in {8, 9}:
                    dim6 = dim[5]
                    dimi = double_or_blank(card, 9 + ioffset, f'ndim={ndim}; dim{ioffset+1}',
                                           default=dim6, end=PBARL_MSG)
                else:
                    dimi = double(card, 9 + ioffset, f'ndim={ndim}; dim{ioffset+1}', end=PBARL_MSG)
                dim.append(dimi)
        else:
            for ioffset in range(ndim):
                dimi = double(card, 9 + ioffset, f'ndim={ndim}; dim{ioffset+1}', end=PBARL_MSG)
                dim.append(dimi)

        #: dimension list
        assert len(dim) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(dim))
        #assert len(dims) == len(self.dim), 'PBARL ndim=%s len(dims)=%s' % (ndim, len(self.dim))

        ndim = len(dim)
        assert ndim > 0, f'PBARL: property_id={pid} dims={dim}'

        nsm = double_or_blank(card, 9 + ndim, 'nsm', default=0.0)
        assert len(dim) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(dim))
        self.cards.append((pid, mid, group, bar_type, dim, nsm, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        ndim = np.zeros(ncards, dtype='int32')
        Type = np.full(ncards, '', dtype='|U8')
        group = np.full(ncards, '', dtype='|U8')
        nsm = np.zeros(ncards, dtype='float64')
        dims = []

        for icard, card in enumerate(self.cards):
            (pid, mid, groupi, bar_type, dim, nsmi, comment) = card
            ndimi = len(dim)
            property_id[icard] = pid
            material_id[icard] = mid
            group[icard] = groupi
            Type[icard] = bar_type
            ndim[icard] = ndimi
            nsm[icard] = nsmi
            assert ndimi > 0, f'PBARL: property_id={pid} dims={dims}'
            dims.extend(dim)
        dims_array = np.array(dims, dtype='float64')
        self._save(property_id, material_id, ndim, Type, group, nsm, dims_array)
        self.sort()
        self.cards = []
        assert len(self.dims) == self.ndim.sum()
        assert isinstance(self.ndim, np.ndarray), self.ndim

    def _save(self, property_id, material_id, ndim, Type, group, nsm, dims) -> None:
        if len(self.property_id):
            property_id = np.hstack([self.property_id, property_id])
            material_id = np.hstack([self.material_id, material_id])
            ndim = np.hstack([self.ndim, ndim])
            Type = np.hstack([self.Type, Type])
            group = np.hstack([self.group, group])
            nsm = np.hstack([self.nsm, nsm])
            dims = np.hstack([self.dims, dims])
        self.property_id = property_id
        self.material_id = material_id
        self.ndim = ndim
        self.Type = Type
        self.group = group
        self.nsm = nsm
        self.dims = dims
        self.n = len(property_id)
        assert len(dims) == self.ndim.sum()
        assert isinstance(self.ndim, np.ndarray), self.ndim

        for pid, beam_type, ndimi, idim in zip(
            self.property_id, self.Type, self.ndim, self.idim):
            idim0, idim1 = idim
            dim = self.dims[idim0 : idim1].tolist()
            ndim = self.valid_types[beam_type]
            assert len(dim) == ndim, f'PBARL pid={pid:d} bar_type={beam_type} ndim={ndim:d} len(dims)={dim}'

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    def convert(self, xyz_scale: float=1.0,
                nsm_per_length_scale: float=1.0, **kwargs):
        self.dims *= xyz_scale
        self.nsm *= nsm_per_length_scale

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          msg=f'no bar materials for {self.type}')
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max())

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        assert isinstance(self.ndim, np.ndarray), self.ndim
        property_ids = array_str(self.property_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        groups = array_default_str(self.group, default='MSCBML0', size=size)
        nsms = array_default_float(self.nsm, default=0.0, size=size)
        for pid, mid, beam_type, ndimi, idim, group, nsm in zip(property_ids, material_ids,
                                                                self.Type, self.ndim,
                                                                self.idim,
                                                                groups, nsms):
            #nsm = set_blank_if_default(nsm, 0.)
            idim0, idim1 = idim
            dim = self.dims[idim0 : idim1].tolist()
            ndim = self.valid_types[beam_type]
            assert len(dim) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(dim))
            list_fields = ['PBARL', pid, mid, group, beam_type, None,
                           None, None, None] + dim + [nsm]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def idim(self) -> np.ndarray:
        idim = make_idim(self.n, self.ndim)
        return idim

    @property
    def all_materials(self) -> list[Any]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[Any]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials), f'{self.type}: all_materials={all_materials}'
        return materials

    def mass_per_length(self) -> np.ndarray:
        assert isinstance(self.ndim, np.ndarray), self.ndim
        nsm = self.nsm
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        if rho.max() == 0. and rho.min() == 0. and nsm.max() == 0. and nsm.min() == 0.:
            return np.zeros(len(rho), dtype=rho.dtype)

        nproperties = len(self.property_id)
        area = self.area()
        mass_per_length = rho * area + nsm

        assert len(mass_per_length) == nproperties
        return mass_per_length

    def mass_per_length_breakdown(self) -> np.ndarray:
        """[rho, A, nsm, mpl]"""
        assert isinstance(self.ndim, np.ndarray), self.ndim
        nsm = self.nsm
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        if rho.max() == 0. and rho.min() == 0. and nsm.max() == 0. and nsm.min() == 0.:
            return np.zeros(len(rho), dtype=rho.dtype)

        nproperties = len(self.property_id)
        area = self.area()
        mass_per_length = rho * area + nsm

        assert len(mass_per_length) == nproperties
        breakdown = np.column_stack([rho, area, nsm, mass_per_length])
        return breakdown

    def area(self) -> np.ndarray:
        nproperties = len(self.property_id)
        area = np.zeros(nproperties, dtype='float64')
        for i, beam_type, idim in zip(count(), self.Type, self.idim):
            idim0, idim1 = idim
            dim = self.dims[idim0:idim1]
            #prop = pbarl(self.property_id[i], self.material_id[i], beam_type, dim)
            #A, I1, I2, I12 = A_I1_I2_I12(prop, beam_type, dim)
            A = _bar_areaL('PBARL', beam_type, dim, self)[0]
            area[i] = A
        return area

    def I(self):
        """A, I1, I2, I12"""
        nprop = len(self.property_id)
        I = np.full((nprop, 4), np.nan, dtype='float64')
        for i, pid, beam_type, idim in zip(count(), self.property_id,
                                           self.Type, self.idim):
            idim0, idim1 = idim
            dim = self.dims[idim0 : idim1].tolist()
            #from pyNastran.bdf.cards.properties.bars import A_I1_I2_I12
            I[i] = _bar_areaL('PBARL', beam_type, dim, self)
        return I

    def i1(self):
        i1 = self.I()[:, 1]
        return i1
    def i2(self):
        i1 = self.I()[:, 2]
        return i1
    def i12(self):
        i12 = self.I()[:, 3]
        return i12


class PBRSECT(Property):
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.form = np.array([], dtype='|U8')

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PBRSECT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : list[str]
            this card is special and is not a ``BDFCard`` like other cards
        comment : str; default=''
            a comment for the card

        """
        line0 = card[0]
        if '\t' in line0:
            line0 = line0.expandtabs()

        bdf_card = BDFCard(to_fields([line0], 'PBMSECT'))
        unused_line0_eq = line0[16:]
        lines_joined = ','.join(card[1:]).replace(' ', '').replace(',,', ',')

        if lines_joined:
            fields = get_beam_sections(lines_joined)
            options = [field.split('=', 1) for field in fields]
            #C:\MSC.Software\MSC.Nastran\msc20051\nast\tpl\zbr3.dat
            #options = [
                #[u'OUTP', u'201'],
                #[u'T', u'1.0'],
                #[u'BRP', u'202'],
                #[u'T(11)', u'[1.2'],
                #[u'PT', u'(202'], [u'224)]'],
                #[u'T(12)', u'[1.2'],
                #[u'PT', u'(224'],
                #[u'205)]'],
            #]
        else:
            options = []

        pid = integer(bdf_card, 1, 'pid')
        mid = integer(bdf_card, 2, 'mid')
        form = string_or_blank(bdf_card, 3, 'form', default='')

        nsm, brps, inps, outp, ts = parse_pbrsect_options(pid, options)
        if options != []:
            print(card)
            raise RuntimeError(f'PBRSECT pid={pid:d}; nsm={nsm} brps={brps} inps={inps} outp={outp} ts={ts}')
        self.cards.append((pid, mid, form, nsm, brps, inps, outp, ts, comment))
        self.n += 1
        #return PBRSECT(pid, mid, form, options, comment=comment)
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        #self.Type = np.full(ncards, '', dtype='|U8')
        #self.group = np.full(ncards, '', dtype='|U8')

        #self.material_id[icard] = mid
        #self.group[icard] = group
        form = np.zeros(ncards, dtype='|U8')
        #self.A = np.zeros(ncards, dtype='float64')
        #self.J = np.zeros(ncards, dtype='float64')

        #self.c = np.zeros((ncards, 2), dtype='float64')
        #self.d = np.zeros((ncards, 2), dtype='float64')
        #self.e = np.zeros((ncards, 2), dtype='float64')
        #self.f = np.zeros((ncards, 2), dtype='float64')

        #self.I = np.zeros((ncards, 3), dtype='float64')
        #self.k = np.zeros((ncards, 2), dtype='float64')
        #self.nsm = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, mid, form, nsm, brps, inps, outp, ts, comment) = card
            assert form in {'', 'GS'}, f'PBRSECT pid={pid} form={form}'
            property_id[icard] = pid
            material_id[icard] = mid
            form[icard] = form
            #self.I[icard, :] = [i1, i2, i12]
            #self.J[icard] = j

            #self.c[icard, :] = [c1, c2]
            #self.d[icard, :] = [d1, d2]
            #self.e[icard, :] = [e1, e2]
            #self.f[icard, :] = [f1, f2]

            #self.k[icard, :] = [k1, k2]
            #self.nsm[icard] = nsm

        self.property_id = property_id
        self.material_id = material_id
        self.form = form
        self.sort()
        self.cards = []

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        #raise RuntimeError('PBRSECT')
        #return ''
        return


class CBARAO(Element):
    """
    Per MSC 2016.1
    +--------+------+-------+------+-----+--------+-----+----+----+
    |   1    |  2   |   3   |  4   |  5  |    6   |  7  | 8  |  9 |
    +========+======+=======+======+=====+========+=====+====+====+
    | CBARAO | EID  | SCALE |  X1  | X2  |  X3    | X4  | X5 | X6 |
    +--------+------+-------+------+-----+--------+-----+----+----+
    | CBARAO | 1065 |  FR   | 0.2  | 0.4 |  0.6   | 0.8 |    |    |
    +--------+------+-------+------+-----+--------+-----+----+----+

    Alternate form (not supported):
    +--------+------+-------+------+-----+--------+-----+----+----+
    |   1    |  2   |   3   |  4   |  5  |    6   |  7  | 8  |  9 |
    +========+======+=======+======+=====+========+=====+====+====+
    | CBARAO | EID  | SCALE | NPTS | X1  | DELTAX |     |    |    |
    +--------+------+-------+------+-----+--------+-----+----+----+
    | CBARAO | 1065 |  FR   |  4   | 0.2 |  0.2   |     |    |    |
    +--------+------+-------+------+-----+--------+-----+----+----+

    """
    @Element.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.scale = np.array([], dtype='|U3')
        self.g0 = np.array([], dtype='int32')
        #self.x = np.array([], dtype='float64')
        self.station = np.array([], dtype='float64')

    def add(self, eid: int, scale: str, x: list[float], comment: str='') -> int:
        """
        Creates a CBARAO card, which defines additional output locations
        for the CBAR card.

        It also changes the OP2 element type from a CBAR-34 to a CBAR-100.
        However, it is ignored if there are no PLOAD1s in the model.
        Furthermore, the type is changed for the whole deck, regardless of
        whether there are PLOAD1s in the other load cases.

        Parameters
        ----------
        eid : int
            element id
        scale : str
            defines what x means
            LE : x is in absolute coordinates along the bar
            FR : x is in fractional
        x : list[float]
            the additional output locations
            len(x) <= 6
        comment : str; default=''
            a comment for the card

        Notes
        -----
        MSC only

        """
        assert isinstance(scale, str), scale
        self.cards.append((eid, scale, x, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CBARAO card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        scale = string(card, 2, 'scale')
        x1_npoints = integer_or_double(card, 3, 'x1/npoints')
        if isinstance(x1_npoints, integer_types):
            npoints = x1_npoints
            assert 0 < npoints < 7, 'CBARAO npoints=%r must be 1-6' % npoints
            x1 = double(card, 4, 'x1')
            delta_x = double(card, 5, 'delta_x')
            x = np.linspace(x1, x1 + delta_x * (npoints-1), num=npoints)
            assert len(x) == npoints, x
        else:
            x = [
                x1_npoints,
                double_or_blank(card, 4, 'x2'),
                double_or_blank(card, 5, 'x3'),
                double_or_blank(card, 6, 'x4'),
                double_or_blank(card, 7, 'x5'),
                double_or_blank(card, 8, 'x6'),
            ]
            x = [xi for xi in x if xi is not None]
        assert len(card) <= 9, f'len(CBARAO card) = {len(card):d}\ncard={card}'
        nstation = len(x)
        assert nstation > 0, x
        self.cards.append((eid, scale, x, comment))
        self.n += 1
        return self.n - 1

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['element_id'].append(self.element_id)

    def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        #LE : x is in absolute coordinates along the bar
        #FR : x is in fractional
        is_le = (self.scale == 'LE')
        nle = is_le.sum()
        if nle:
            for ile, (istation0, istation1) in zip(is_le, self.istation):
                if not ile:
                    continue
                self.station[istation0:istation1] *= xyz_scale

    def __apply_slice__(self, elem: CBARAO, i: np.ndarray) -> None:
        elem.element_id = self.element_id[i]
        elem.scale = self.scale[i]
        istation = self.istation # [i, :]
        elem.station = hslice_by_idim(i, istation, self.station)
        elem.nstation = self.nstation[i]
        elem.n = len(i)

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        station: list[float] = []

        scale = np.zeros(ncards, dtype='|U2')
        nstation = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (eid, scalei, stationi, comment) = card
            element_id[icard] = eid
            assert scalei in {'FR', 'LE'}, (eid, scalei)
            scale[icard] = scalei
            nstationi = len(stationi)
            nstation[icard] = nstationi
            assert nstationi > 0, stationi
            station.extend(stationi)

        self._save(element_id, scale, nstation, station)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray,
              scale: np.ndarray,
              nstation: np.ndarray,
              station: list[float]) -> None:
        if len(self.element_id) != 0:
            raise NotImplementedError()
        self.element_id = element_id
        self.scale = scale
        self.nstation = nstation
        self.station = np.array(station, dtype=self.model.fdtype)


    @property
    def istation(self) -> np.ndarray:
        istation = make_idim(self.n, self.nstation)
        return istation

    @property
    def max_id(self) -> int:
        return self.element_id.max()

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)


        element_ids = array_str(self.element_id, size=size)

        for eid, scale, (istation0, istation1) in zip(element_ids, self.scale, self.istation):
            station = self.station[istation0:istation1].tolist()
            assert len(station) > 0, station
            list_fields = ['CBARAO', eid, scale] + station
            #print(list_fields)
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pbar, self.model.pbarl]
                if prop.n > 0]

    def mass_material_id(self) -> np.ndarray:
        material_id = basic_mass_material_id(self.property_id, self.allowed_properties, 'CBAR')
        return material_id

    def geom_check(self, missing: dict[str, np.ndarray]):
        eid = self.model.cbar.element_id
        #pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          #msg=f'no bar properties for {self.type}')
        eid.sort()
        geom_check(self,
                   missing,
                   element_id=(eid, self.element_id),
                   )

#class CBARAO(BaseCard):

    #def __init__(self, eid, scale, x, comment=''):
        #"""
        #Creates a CBARAO card, which defines additional output locations
        #for the CBAR card.

        #It also changes the OP2 element type from a CBAR-34 to a CBAR-100.
        #However, it is ignored if there are no PLOAD1s in the model.
        #Furthermore, the type is changed for the whole deck, regardless of
        #whether there are PLOAD1s in the other load cases.

        #Parameters
        #----------
        #eid : int
            #element id
        #scale : str
            #defines what x means
            #LE : x is in absolute coordinates along the bar
            #FR : x is in fractional
        #x : list[float]
            #the additional output locations (doesn't include the end points)
            #len(x) <= 6
        #comment : str; default=''
            #a comment for the card

        #MSC only

        #"""
        #if comment:
            #self.comment = comment
        #self.eid = eid
        #self.scale = scale
        #self.x = np.unique(x).tolist()

def apply_bar_default(bar: Union[CBAR, CBEAM],
                      baror: Union[BAROR, BEAMOR]) -> None:
    model = bar.model
    data_temp_default = []
    PROPERTY_ID_DEFAULT = 0
    G0_DEFAULT = 0
    OFFT_DEFAULT = ''

    if baror is not None:
        data_temp_default.append((bar.property_id, PROPERTY_ID_DEFAULT, baror.pid))
        g0_ = baror.g0
        x_ = baror.x
        offt_ = baror.offt
        model.log.info(f'setting default g0 as {g0_}')
        model.log.info(f'setting default x_ as {x_}')
        model.log.info(f'setting default offt_ as {offt_}')
    else:
        offt_ = 'GGG'
        g0_ = 0
        x_ = np.zeros(3, dtype='float64')

    data_temp_default += [
        #(bar.property_id, PROPERTY_ID_DEFAULT, baror.pid),
        (bar.property_id, PROPERTY_ID_DEFAULT, bar.element_id),
    ]
    if baror is None or baror.type == 'BAROR':
        data_temp_default.append((bar.offt, OFFT_DEFAULT, offt_))
    else:
        if isinstance(baror.offt, str):
            data_temp_default.append((bar.offt, OFFT_DEFAULT, offt_))
        else:
            data_temp_default.append((bar.bit, -1, baror.offt))

    #print(bar.x)
    #print(bar.g0)
    is_x_nan = np.isnan(bar.x)
    x_any_nan = ~np.any(is_x_nan, axis=1) # x is blank
    x_all_nan = ~np.all(is_x_nan, axis=1) # x is blank
    g0_blank = (bar.g0 == G0_DEFAULT)     # g0 is blank

    # we need to fix things that have:
    #   - x = nan (for all values)
    #   - g0 = 0
    #
    # (x is nan) or (g0 == 0)
    is_blank = x_all_nan | g0_blank
    #print('x_all_nan =', x_all_nan)
    #print('g0_blank =', g0_blank)
    #print('is_blank =', is_blank)

    #g0     : 0
    #offt   : 1
    #pid    : 1
    #type   : 'BAROR'
    #x      : array([1., 0., 0.])
    #print(baror.get_stats())
    if len(is_blank):
        if g0_ == 0:
            # use x
            ix = is_x_nan[:, 0]
            iy = is_x_nan[:, 1]
            iz = is_x_nan[:, 2]
            bar.x[ix, 0] = x_[0]
            bar.x[iy, 1] = x_[1]
            bar.x[iz, 2] = x_[2]
            #model.log.info(f'x as {bar.x}')
        else:
            # use g0
            bar.g0[g0_blank] = g0_
        # no value was set for g0 or x
        # use 0
    else:
        # no 100% blank data
        pass
    is_x_nan = ~np.isnan(bar.x)
    x_any_nan = ~np.any(is_x_nan, axis=1) # x is blank
    g0_blank = (bar.g0 == G0_DEFAULT)     # g0 is blank
    is_blank = x_any_nan & g0_blank
    assert is_blank.sum() == 0, is_blank

    #if g0_ == 0:
        # use x default
    for data, temp_value, default_value in data_temp_default:
        ibad = np.where(data == temp_value)[0]
        if len(ibad):
            if isinstance(default_value, np.ndarray):
                data[ibad] = default_value[ibad]
            else:
                data[ibad] = default_value
    assert bar.property_id.min() > 0, bar.property_id


def inertia_from_property_id(property_id: np.ndarray,
                             allowed_properties: list[Property]) -> np.ndarray:
    npid = len(property_id)
    inertia = np.full((npid, 4), np.nan, dtype='float64')
    for prop in allowed_properties:
        i_lookup, i_all = searchsorted_filter(prop.property_id, property_id, msg='')
        if len(i_lookup) == 0:
            continue

        # we're at least using some properties
        inertiai = prop.inertia()
        inertia[i_lookup] = inertiai[i_all]
    return inertia

def k_from_property_id(property_id: np.ndarray,
                       allowed_properties: list[Property]) -> np.ndarray:
    npid = len(property_id)
    k1_k2 = np.full((npid, 2), np.nan, dtype='float64')
    for prop in allowed_properties:
        i_lookup, i_all = searchsorted_filter(prop.property_id, property_id, msg='')
        if len(i_lookup) == 0:
            continue

        # we're at least using some properties
        if prop.type in {'PBAR', 'PBEAM'}:
            k1_k2i = prop.k
        else:
            k1_k2i = prop.k()
        k1_k2[i_lookup] = k1_k2i[i_all]
    return k1_k2

