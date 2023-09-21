from __future__ import annotations
from itertools import count
from typing import Union, Optional, Any, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16 # , print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string, integer_or_double,
    integer_or_blank, double_or_blank, string_or_blank,
    integer_double_or_blank, integer_string_or_blank,
    blank)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default # init_x_g0,
from pyNastran.bdf.cards.properties.bars import _bar_areaL, to_fields, get_beam_sections, parse_pbrsect_options
# PBARL as pbarl, A_I1_I2_I12

from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import Element, Property, make_idim, hslice_by_idim, searchsorted_filter, get_print_card_8_16
from pyNastran.dev.bdf_vectorized3.cards.elements.rod import line_mid_mass_per_length, line_length, line_vector_length, line_centroid
from pyNastran.dev.bdf_vectorized3.cards.elements.utils import get_density_from_material, basic_mass_material_id
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int, array_default_str
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg, cast_int_array
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
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

    @classmethod
    def add_card(cls, card, comment=''):
        PROPERTY_ID_DEFAULT = 0
        GO_X_DEFAULT = 0
        OFFT_DEFAULT = ''

        pid = integer_or_blank(card, 2, 'pid', default=PROPERTY_ID_DEFAULT)

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
                          double_or_blank(card, 6, 'x2', default=0.),
                          double_or_blank(card, 7, 'x3', default=0.)],
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
        g0 = -1
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


class CBAR(Element):
    def __init__(self, model: BDF):
        super().__init__(model)
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
        return self.n

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

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        element_id = []
        property_id = np.zeros(ncards, dtype='int32')
        nodes = []
        offt = np.full(ncards, '', dtype='|U3')
        g0 = np.zeros(ncards, dtype='int32')
        x = np.full((ncards, 3), np.nan, dtype='float64')

        pa = np.zeros(ncards, dtype='int32')
        pb = np.zeros(ncards, dtype='int32')
        wa = np.zeros((ncards, 3), dtype='float64')
        wb = np.zeros((ncards, 3), dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, xi, g0i, offti, pai, pbi, wai, wbi, comment) = card
            element_id.append(eid)
            property_id[icard] = pid
            nodes.append(nids)

            if g0i is None:
                g0i = 0
            else:
                xi = [np.nan, np.nan, np.nan]
            g0[icard] = g0i
            x[icard, :] = xi
            offt[icard] = offti

            pa[icard] = pai
            pb[icard] = pbi
            wa[icard, :] = wai
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
        self.element_id = cast_int_array(element_id)
        self.property_id = property_id
        #print(element_id)
        #print(property_id)
        self.nodes = cast_int_array(nodes)
        self.offt = offt
        self.g0 = g0
        self.x = x

        self.pa = pa
        self.pb = pb
        self.wa = wa
        self.wb = wb

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.element_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)
        pas = array_default_int(self.pa, default=0, size=size)
        pbs = array_default_int(self.pb, default=0, size=size)
        for eid, pid, nodes, g0, x, offt, pa, pb, wa, wb in zip(element_ids, property_ids, nodes,
                                                                self.g0, self.x, self.offt, pas, pbs, self.wa, self.wb):
            n1, n2 = nodes
            w1a = set_blank_if_default(wa[0], 0.0)
            w2a = set_blank_if_default(wa[1], 0.0)
            w3a = set_blank_if_default(wa[2], 0.0)

            w1b = set_blank_if_default(wb[0], 0.0)
            w2b = set_blank_if_default(wb[1], 0.0)
            w3b = set_blank_if_default(wb[2], 0.0)
            if g0 == -1:
                x1, x2, x3 = x # self.get_x_g0_defaults()
            else:
                x1 = g0
                x2 = ''
                x3 = ''

            # offt doesn't exist in NX nastran
            offt = set_blank_if_default(offt, 'GGG')

            list_fields = ['CBAR', eid, pid, n1, n2,
                           x1, x2, x3, offt, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
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

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no bar properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

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
    #def slice_card_by_property_id(self, property_id: np.ndarray) -> PBAR:
        #"""uses a node_ids to extract PBARs"""
        #iprop = self.index(property_id)
        #prop = self.slice_card_by_index(iprop)
        #return prop

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
        return self.n

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
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
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

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.property_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        property_ids = array_str(self.property_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        for pid, mid, A, I, j, nsm, c, d, e, f, k in zip(property_ids, material_ids,
                                                         self.A, self.I, self.J,
                                                         self.nsm,
                                                         self.c, self.d, self.e, self.f,
                                                         self.k):
            c1, c2 = c
            d1, d2 = d
            e1, e2 = e
            f1, f2 = f
            i1, i2, i12 = I
            k1, k2 = k

            i1 = set_blank_if_default(i1, 0.0)
            i2 = set_blank_if_default(i2, 0.0)
            i12 = set_blank_if_default(i12, 0.0)
            j = set_blank_if_default(j, 0.0)
            nsm = set_blank_if_default(nsm, 0.0)

            c1 = set_blank_if_default(c1, 0.0)
            c2 = set_blank_if_default(c2, 0.0)

            d1 = set_blank_if_default(d1, 0.0)
            d2 = set_blank_if_default(d2, 0.0)

            e1 = set_blank_if_default(e1, 0.0)
            e2 = set_blank_if_default(e2, 0.0)

            f1 = set_blank_if_default(f1, 0.0)
            f2 = set_blank_if_default(f2, 0.0)

            k1 = set_blank_if_default(k1, 1e8)
            k2 = set_blank_if_default(k2, 1e8)

            list_fields = ['PBAR', pid, mid, A, i1, i2, j, nsm,
                           None, c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, i12]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def all_materials(self) -> list[Any]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[Any]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.materials}'
        return materials

    def mass_per_length(self) -> np.ndarray:
        mass_per_length = line_mid_mass_per_length(self.material_id, self.nsm, self.A,
                                                   self.allowed_materials)
        return mass_per_length

    def area(self) -> np.ndarray:
        return self.A


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
    def __init__(self, model: BDF):
        super().__init__(model)
        #self.model = model
        self.material_id = np.array([], dtype='int32')
        self.ndim = np.array([], dtype='int32')
        self.Type = np.array([], dtype='|U8')
        self.group = np.array([], dtype='|U8')
        self.nsm = np.array([], dtype='float64')
        self.dims = np.array([], dtype='float64')

    def slice_card_by_property_id(self, property_id: np.ndarray) -> PBARL:
        """uses a node_ids to extract PBARLs"""
        iprop = self.index(property_id)
        prop = self.slice_card_by_index(iprop)
        return prop

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
        self.write()
        assert self.ndim.sum() == len(self.dims)
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        prop.Type = self.Type[i]
        prop.group = self.group[i]
        prop.nsm = self.nsm[i]

        idim = self.idim # [i, :]
        prop.dims = hslice_by_idim(i, idim, self.dims)

        prop.ndim = self.ndim[i]
        assert isinstance(prop.ndim, np.ndarray), prop.ndim
        assert prop.ndim.sum() == len(prop.dims), f'prop.ndim={prop.ndim} len(prop.dims)={len(prop.dims)}'
        self.write()

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
        return self.n

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
        assert len(card) <= 17, f'len(CBAR card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, mid, group, bar_type, dim, nsm, comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
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

        self.write()

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          msg=f'no bar materials for {self.type}')
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.property_id) == 0:
            return
        print_card = get_print_card_8_16(size)
        assert isinstance(self.ndim, np.ndarray), self.ndim
        property_ids = array_str(self.property_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        groups = array_default_str(self.group, default='MSCBML0', size=size)
        for pid, mid, beam_type, ndimi, idim, group, nsm in zip(property_ids, material_ids,
                                                                self.Type, self.ndim,
                                                                self.idim,
                                                                groups, self.nsm):
            nsm = set_blank_if_default(nsm, 0.)
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

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        self.property_id = np.zeros(ncards, dtype='int32')
        self.material_id = np.zeros(ncards, dtype='int32')
        #self.Type = np.full(ncards, '', dtype='|U8')
        #self.group = np.full(ncards, '', dtype='|U8')

        #self.material_id[icard] = mid
        #self.group[icard] = group
        self.form = np.zeros(ncards, dtype='|U8')
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
            self.property_id[icard] = pid
            self.material_id[icard] = mid
            self.form[icard] = form
            #self.I[icard, :] = [i1, i2, i12]
            #self.J[icard] = j

            #self.c[icard, :] = [c1, c2]
            #self.d[icard, :] = [d1, d2]
            #self.e[icard, :] = [e1, e2]
            #self.f[icard, :] = [f1, f2]

            #self.k[icard, :] = [k1, k2]
            #self.nsm[icard] = nsm
        self.sort()
        self.cards = []

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if self.n == 0:
            return
        #raise RuntimeError('PBRSECT')
        #return ''


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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')
        self.scale = np.array([], dtype='|U3')
        self.g0 = np.array([], dtype='int32')
        self.x = np.array([], dtype='float64')

    #def add(self, eid: int, pid: int, nids: list[int],
            #x: Optional[list[float]], g0: Optional[int],
            #offt: str='GGG', pa: int=0, pb: int=0,
            #wa: Optional[list[float]]=None, wb: Optional[list[float]]=None,
            #comment: str='', validate: bool=False):
        #assert x is None or g0 is None, f'pid={pid} x={x} g0={g0}'
        #self.cards.append((eid, pid, nids, x, g0, offt, pa, pb, wa, wb, comment))
        #self.n += 1

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
        self.cards.append((eid, scale, x, comment))
        self.n += 1
        return self.n

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
        self.cards.append((eid, scale, x, comment))
        self.n += 1
        return self.n

    def __apply_slice__(self, elem: CBARAO, i: np.ndarray) -> None:
        elem.element_id = self.element_id[i]
        elem.scale = self.scale[i]
        istation = self.istation # [i, :]
        elem.station = hslice_by_idim(i, istation, self.station)
        elem.nstation = self.nstation[i]
        elem.n = len(i)

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        element_id = []
        station = []

        scale = np.zeros(ncards, dtype='|U2')
        nstation = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (eid, scalei, stationi, comment) = card
            element_id.append(eid)
            assert scalei in ['FR'], (eid, scalei)
            scale[icard] = scalei
            nstation[icard] = len(station)
            station.extend(stationi)
        self._save(element_id, scale, nstation, station)
        self.sort()
        self.cards = []

    def _save(self, element_id, scale, nstation, station) -> None:
        if len(self.element_id) != 0:
            raise NotImplementedError()
        self.element_id = cast_int_array(element_id)
        self.scale = scale
        self.nstation = nstation
        self.station = np.array(station, dtype=self.model.fdtype)

    @property
    def istation(self) -> np.ndarray:
        istation = make_idim(self.n, self.nstation)
        return istation

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.element_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        element_ids = array_str(self.element_id, size=size)

        for eid, scale, (istation0, istation1) in zip(element_ids, self.scale, self.istation):
            station = self.station[istation0:istation1]
            list_fields = ['CBARAO', eid, self.scale] + station.tolist()
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
        (bar.offt, OFFT_DEFAULT, offt_),
    ]
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
