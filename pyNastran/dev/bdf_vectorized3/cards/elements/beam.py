from __future__ import annotations
from itertools import count, zip_longest
from typing import Union, Optional, Any, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import zip_strict, integer_types, float_types
from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16 # , print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string, blank,
    integer_or_double,
    integer_or_blank, double_or_blank, string_or_blank,
    integer_double_or_blank, integer_string_or_blank, double_string_or_blank)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default # init_x_g0,
from pyNastran.bdf.cards.properties.bars import _bar_areaL # PBARL as pbarl, A_I1_I2_I12
from pyNastran.utils.mathematics import integrate_positive_unit_line # integrate_unit_line,

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, make_idim, hslice_by_idim, searchsorted_filter, # vslice_by_idim,
    parse_element_check, parse_property_check,
    get_print_card_8_16
)
from .rod import line_pid_mass_per_length, line_length, line_vector_length, line_centroid, e_g_nu_from_property_id
from .bar import (apply_bar_default, init_x_g0, get_bar_vector, split_offt_vector,
                  inertia_from_property_id, k_from_property_id,
                  e_g_nu_from_isotropic_material)
from .utils import get_density_from_material
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_int, array_default_float, array_default_str)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from ..coord import COORD
    from ..materials import MAT1


class BEAMOR(BaseCard):
    """
    +--------+-----+---+---+---+-------+-----+-------+------+
    |    1   |  2  | 3 | 4 | 5 |   6   |  7  |   8   |  9   |
    +========+=====+===+===+===+=======+=====+=======+======+
    | BEAMOR | PID |   |   |   | G0/X1 |  X2 |  X3   | OFFT |
    +--------+-----+---+---+---+-------+-----+-------+------+
    | BEAMOR | 39  |   |   |   |  0.6  | 2.9 | -5.87 | GOG  |
    +--------+-----+---+---+---+-------+-----+-------+------+

    """
    type = 'BEAMOR'
    def __init__(self, pid, g0, x, offt='GGG', comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.g0 = g0
        self.x = x
        self.offt = offt

    #@classmethod
    #def _init_from_empty(cls):
        #pid = 1
        #is_g0 = True
        #g0 = 1
        #x = None
        #return BEAMOR(pid, is_g0, g0, x, offt='GGG', comment='')

    @classmethod
    def add_card(cls, card, comment=''):
        PROPERTY_ID_DEFAULT = 0
        GO_X_DEFAULT = 0
        OFFT_DEFAULT = ''
        pid = integer_or_blank(card, 2, 'pid', default=PROPERTY_ID_DEFAULT)

        # x / g0
        field5 = integer_double_or_blank(card, 5, 'g0_x1', default=0.0)
        if isinstance(field5, integer_types):
            #is_g0 = True
            g0 = field5
            x = [np.nan, np.nan, np.nan]
            blank(card, 6, 'x2')
            blank(card, 7, 'x3')
        elif isinstance(field5, float):
            #is_g0 = False
            g0 = GO_X_DEFAULT
            x = np.array([field5,
                          double_or_blank(card, 6, 'x2', default=0.0),
                          double_or_blank(card, 7, 'x3', default=0.0)],
                         dtype='float64')
        else:
            raise NotImplementedError('BEAMOR field5 = %r' % field5)
        offt = integer_string_or_blank(card, 8, 'offt', OFFT_DEFAULT)
        assert len(card) <= 9, f'len(BEAMOR card) = {len(card):d}\ncard={card}'
        return BEAMOR(pid, g0, x, offt=offt, comment=comment)

    def raw_fields(self):
        return ['BEAMOR', None, self.pid, None, None] + list(self.x) + [self.offt]

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class CBEAM(Element):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |    9     |
    +=======+=====+=====+=====+=====+=====+=====+=====+==========+
    | CBEAM | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |     |     |     |     |     |          |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |    9     |
    +=======+=====+=====+=====+=====+=====+=====+=====+==========+
    | CBEAM | EID | PID | GA  | GB  | G0  |     |     | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |     |     |     |     |     |          |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+

    bit is an MSC specific field
    NX 2020 added offt
    """

    @Element.clear_check
    def clear(self) -> None:
        self.element_id: np.array = np.array([], dtype='int32')
        self.property_id: np.array = np.array([], dtype='int32')
        self.nodes: np.array = np.zeros((0, 2), dtype='int32')
        self.offt: np.array = np.array([], dtype='|U3')
        self.bit: np.array = np.array([], dtype='int32')
        self.g0: np.array = np.array([], dtype='int32')
        self.x: np.array = np.zeros((0, 3), dtype='float64')

        # pin flags
        self.pa: np.array = np.array([], dtype='int32')
        self.pb: np.array = np.array([], dtype='int32')

        # offset vectors at A/B
        self.wa: np.array = np.zeros((0, 3), dtype='float64')
        self.wb: np.array = np.zeros((0, 3), dtype='float64')

        # scalar points at end A/B for warping
        self.sa: np.array = np.zeros([], dtype='int32')
        self.sb: np.array = np.zeros([], dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int],
            x: Optional[list[float]], g0: Optional[int],
            offt: str='GGG', bit=None,
            pa: int=0, pb: int=0,
            wa=None, wb=None,
            sa: int=0, sb: int=0, comment: str='') -> int:
        if wa is None:
            wa = [0., 0., 0.]
        if wb is None:
            wb = [0., 0., 0.]
        if bit:
            assert isinstance(bit, integer_types), f'offt/bit={bit!r} and should be an integer'
            self.cards.append((eid, pid, nids, g0, x, bit, [pa, pb], wa, wb, sa, sb, comment))
        else:
            self.cards.append((eid, pid, nids, g0, x, offt, [pa, pb], wa, wb, sa, sb, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        PROPERTY_ID_DEFAULT = 0
        OFFT_DEFAULT = ''
        eid = integer(card, 1, 'eid')

        pid = integer_or_blank(card, 2, 'pid', default=PROPERTY_ID_DEFAULT)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        #x, g0 = init_x_g0(card, eid, x1_default, x2_default, x3_default)
        x, g0 = init_x_g0(card, eid)
        # doesn't exist in NX nastran
        offt = integer_string_or_blank(card, 8, 'offt', default=OFFT_DEFAULT)

        pa = integer_or_blank(card, 9, 'pa', default=0)
        pb = integer_or_blank(card, 10, 'pb', default=0)

        wa = np.array([double_or_blank(card, 11, 'w1a', default=0.0),
                       double_or_blank(card, 12, 'w2a', default=0.0),
                       double_or_blank(card, 13, 'w3a', default=0.0)], dtype='float64')

        wb = np.array([double_or_blank(card, 14, 'w1b', 0.0),
                       double_or_blank(card, 15, 'w2b', 0.0),
                       double_or_blank(card, 16, 'w3b', 0.0)], dtype='float64')
        sa = integer_or_blank(card, 17, 'sa', default=0)
        sb = integer_or_blank(card, 18, 'sb', default=0)

        assert len(card) <= 19, f'len(CBEAM card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, [ga, gb], g0, x, offt, [pa, pb], wa, wb, sa, sb, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        offt = np.full(ncards, '', dtype='|U3')
        bit = np.full(ncards, -1, dtype='int32')
        g0 = np.zeros(ncards, dtype=idtype)
        x = np.full((ncards, 3), np.nan, dtype='float64')

        pa = np.zeros(ncards, dtype='int32')
        pb = np.zeros(ncards, dtype='int32')
        wa = np.zeros((ncards, 3), dtype='float64')
        wb = np.zeros((ncards, 3), dtype='float64')
        sa = np.zeros(ncards, dtype='int32')
        sb = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, g0i, xi, offti, [pai, pbi],
             wai, wbi, sai, sbi, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            if g0i in {None, 0, -1}:
                x[icard, :] = xi
            else:
                assert g0i > 0, card
                g0[icard] = g0i
            if isinstance(offti, str):
                offt[icard] = offti
            else:
                assert isinstance(offti, integer_types), f'offt/bit={offti!r} and should be an integer'
                bit[icard] = offti
            pa[icard] = pai
            pb[icard] = pbi
            wa[icard, :] = wai
            wb[icard, :] = wbi
            sa[icard] = sai
            sb[icard] = sbi

        self._save(element_id, property_id, nodes,
                   offt, bit,
                   g0, x,
                   pa, pb, wa, wb, sa, sb)
        beamor = self.model.beamor
        apply_bar_default(self, beamor)
        self.cards = []

    def _save(self, element_id, property_id, nodes,
              offt, bit,
              g0, x,
              pa, pb, wa, wb, sa, sb) -> None:
        if len(self.element_id) != 0:
            asdf
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.offt = offt
        self.bit = bit
        self.g0 = g0
        self.x = x

        self.pa = pa
        self.pb = pb
        self.wa = wa
        self.wb = wb
        self.sa = sa
        self.sb = sb
        self.n = len(property_id)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['element_id'].append(self.element_id)
        used_dict['property_id'].append(self.property_id)
        used_dict['node_id'].append(self.nodes.ravel())
        g0 = self.g0[self.is_g0]
        if len(g0):
            used_dict['node_id'].append(g0)

        sa = self.sa[self.sa != 0]
        sb = self.sb[self.sb != 0]
        if len(sa):
            used_dict['spoint_id'].append(sa)
        if len(sb):
            used_dict['spoint_id'].append(sb)

    def convert(self, xyz_scale: float=1.0,
                mass_scale: float=1.0, **kwargs):
        # easy
        self.wa *= xyz_scale
        self.wb *= xyz_scale

        ## TODO: probably wrong for CD=1
        self.x *= xyz_scale

    def __apply_slice__(self, elem: CBEAM, i: np.ndarray) -> None:
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.offt = self.offt[i]
        elem.bit = self.bit[i]
        elem.g0 = self.g0[i]
        elem.x = self.x[i, :]
        elem.pa = self.pa[i]
        elem.pb = self.pb[i]
        elem.wa = self.wa[i, :]
        elem.wb = self.wb[i, :]
        elem.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no beam properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes_ = array_str(self.nodes, size=size)
        offts = array_default_str(self.offt, default='GGG', size=size)
        bits = array_str(self.bit, size=size)
        ibit = (self.bit != -1)
        offts[ibit] = bits[ibit]
        pas = array_default_int(self.pa, default=0, size=size)
        pbs = array_default_int(self.pb, default=0, size=size)
        was = array_default_float(self.wa, default=0, size=size, is_double=False)
        wbs = array_default_float(self.wb, default=0, size=size, is_double=False)
        for eid, pid, nodes, g0, x, is_g0, offt, pa, pb, wa, wb in zip_longest(
            element_ids, property_ids, nodes_,
            self.g0, self.x, self.is_g0, offts,
            pas, pbs, was, wbs):

            n1, n2 = nodes
            w1a, w2a, w3a = wa
            w1b, w2b, w3b = wb
            if is_g0:
                x1 = g0
                x2 = ''
                x3 = ''
            else:
                x1, x2, x3 = x # self.get_x_g0_defaults()

            # offt doesn't exist in NX nastran
            #offt = set_blank_if_default(offt, 'GGG')

            list_fields = ['CBEAM', eid, pid, n1, n2,
                           x1, x2, x3, offt, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def is_x(self) -> np.ndarray:
        return (self.g0 == 0)

    @property
    def is_g0(self) -> np.ndarray:
        return ~self.is_x

    @property
    def all_properties(self) -> list[Union[PBEAM, PBEAML, PBCOMP]]:
        model = self.model
        return [model.pbeam, model.pbeaml, model.pbcomp]

    @property
    def allowed_properties(self):
        all_properties = self.all_properties
        props = [prop for prop in all_properties if prop.n > 0]
        assert len(props) > 0, f'{self.type}: all_props={all_properties}'
        return props

    def mass(self) -> np.ndarray:
        #pid = self.property_id
        mass_per_length = line_pid_mass_per_length(self.property_id, self.allowed_properties)
        length = self.length()
        mass = mass_per_length * length
        return mass

    def line_vector_length(self) -> tuple[np.ndarray, np.ndarray]:
        line_vector, length = line_vector_length(self.model, self.nodes)
        return line_vector, length

    def length(self) -> np.ndarray:
        missing = np.setdiff1d(self.nodes.flatten(), self.model.grid.node_id)
        if len(missing):
            raise RuntimeError(missing)
        length = line_length(self.model, self.nodes)
        inan = np.isnan(length)
        if np.any(inan):
            msg = 'CBEAM has nan length\n'
            msg += f'eids={self.element_id[inan]}\n'
            msg += f'nid1={self.nodes[inan,0]}\n'
            msg += f'nid2={self.nodes[inan,1]}\n'
            raise RuntimeError(msg)
        return length

    def centroid(self) -> np.ndarray:
        centroid = line_centroid(self.model, self.nodes)
        return centroid

    def e_g_nu(self) -> np.ndarray:
        e_g_nu = e_g_nu_from_property_id(self.property_id, self.allowed_properties)
        return e_g_nu

    def inertia(self) -> np.ndarray:
        inertia = inertia_from_property_id(self.property_id,
                                           self.allowed_properties)
        return inertia

    def k(self) -> np.ndarray:
        k1_k2 = k_from_property_id(self.property_id,
                                   self.allowed_properties)
        return k1_k2

    def get_bar_vector(self, xyz1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        v, cd = get_bar_vector(self, xyz1)
        return v, cd

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

    def get_axes(self, xyz1: np.ndarray, xyz2: np.ndarray,
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                            np.ndarray, np.ndarray, np.ndarray]:
        log = self.model.log
        coords = self.model.coord
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
            if np.any(cd1_v_vector):
                #v[icd1_vector, :] = np.nan
                cd1_ref: COORD = coords.slice_card_by_id(cd1_v_vector)
                v1v = v[icd1_v_vector, :]
                v[icd1_v_vector, :] = cd1_ref.transform_xyz_to_global_assuming_rectangular(v1v)
                del v1v
            del icd1_v_vector, cd1_v_vector

        #elif offt_vector == 'B':
            # basic - cid = 0
            #pass

        if np.any(np.isnan(v.max(axis=1))):
            raise RuntimeError(f'v = {v}')

        #--------------------------------------------------------------------------
        # determine the bar vectors
        #log.info(f'v =\n{v}')
        #log.info(f'ihat =\n{i_offset}')
        ihat = i_offset

        vnorm = np.linalg.norm(v, axis=1)

        #if np.any(np.isnan(v.max(axis=1))):
        #print(f'vnorm = {vnorm}')

        vhat = v / vnorm[:, np.newaxis] # j
        z = np.cross(ihat, vhat) # k
        norm_z = np.linalg.norm(z, axis=1)
        assert len(norm_z) == neids

        #if np.any(np.isnan(zhat.max(axis=1))):
        #print(f'norm_z = {norm_z}')

        zhat = z / norm_z[:, np.newaxis]
        yhat = np.cross(zhat, ihat) # j
        norm_i = np.linalg.norm(ihat, axis=1)
        norm_yhat = np.linalg.norm(yhat, axis=1)
        xform_offset = np.dstack([ihat, yhat, zhat]) # 3x3 unit matrix
        #del ihat, yhat, zhat, norm_z, norm_yhat

        if np.any(np.isnan(yhat.max(axis=1))):
            self.model.log.error(f'norm_yhat = {norm_yhat}')

        del norm_i, norm_z, norm_yhat
        #aaa
        #--------------------------------------------------------------------------
        # rotate wa
        # wa defines the offset at end A
        wa = self.wa.copy()  # we're going to be inplace hacking it, so copy :)

        if np.any(is_rotate_wa_g):
            icd1_vector = (is_rotate_wa_g) & (cd1 != 0)
            cd1_vector = cd1[icd1_vector]
            if np.any(icd1_vector):
                cd1_ref = coords.slice_card_by_id(cd1_vector)
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
                cd2_ref = coords.slice_card_by_id(cd2_vector)
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

    #def check_missing_ids(self, property_id: np.ndarray):
        #missing_coords = np.setdiff1d(coord_id, self.coord_id)
        #if len(missing_coords):
            #raise RuntimeError(f'coords={missing_coords} not found in {self.coord_id}')

    def center_of_mass(self) -> np.ndarray:
        #self.check_missing(self.property_id)
        log = self.model.log

        xyz1, xyz2 = self.get_xyz()
        neids = xyz1.shape[0]
        centroid = (xyz1 + xyz2) / 2.
        assert centroid.shape[0] == self.nodes.shape[0]
        assert not np.isnan(np.max(xyz1)), xyz1
        assert not np.isnan(np.max(xyz2)), xyz2

        v, ihat, jhat, khat, wa, wb = self.get_axes(xyz1, xyz2)

        # we finally have the nodal coordaintes!!!! :)
        p1 = xyz1 + wa
        p2 = xyz2 + wb
        # ----------------------------------
        # now some mass properties :(
        mass_per_length = np.full(neids, np.nan, dtype='float64')
        nsm_per_length = np.full(neids, np.nan, dtype='float64')
        nsm_centroid = np.full((neids, 3), np.nan, dtype='float64')

        #log.debug(f'property_id = {self.property_id}')
        for prop in self.allowed_properties:
            pids_common = np.intersect1d(prop.property_id, self.property_id)
            #ind = prop.property_id[ipid]
            if len(pids_common) == 0:
                log.debug(f'  skipping {prop.type}; pids={prop.property_id}')
                continue

            if 0:
                ipid = np.searchsorted(prop.property_id, self.property_id)
                ipid = ipid[ipid < len(prop.property_id)]
                if len(ipid) == 0:
                    log.warning(f'skipping {prop.type}; pids={prop.property_id}')
                    continue
                is_valid = (prop.property_id[ipid] == self.property_id)
                ipid = ipid[is_valid]
            else:
                ipid = np.array([i for i, pid  in enumerate(self.property_id)
                                 if pid in prop.property_id])
                if len(ipid) == 0:
                    log.warning(f'skipping {prop.type}; pids={prop.property_id}')
                    continue

            prop2 = prop.slice_card_by_property_id(pids_common)
            log.info(f'running...{prop.type}: pids={prop.property_id}')
            if prop.type == 'PBEAM':
                #ipid = prop.index(self.property_id)
                #ipid = np.array([pid for pid in self.property_id
                                 #if pid in pids_common])
                #ipid = prop.index(pids_common, assume_sorted=True,
                                  #inverse=False)
                #ipidrev = prop.index(pids_common, assume_sorted=True,
                                     #inverse=True)
                #prop2 = prop.slice_card_by_id(ipid)
                m1a = prop2.m1a
                m1b = prop2.m1b
                m2a = prop2.m2a
                m2b = prop2.m2b
                #rho = prop2.rho()

                # we don't call the MassPerLength method so we can put the NSM centroid
                # on a different axis (the PBEAM is weird)
                mass_per_lengths_, nsm_per_lengths_ = prop2.rhoarea_nsm()

                #ipidrev2 = np.zeros(neids, )
                for jpid, pid, mpl, nsmpl in zip(count(), pids_common,
                                                 mass_per_lengths_, nsm_per_lengths_):
                    iipid = np.where(self.property_id == pid)[0]
                    if len(iipid) == 0:
                        log.warning(f'  skipping {prop.type}; pid={pid}')
                        continue
                    #ipidrev2.append(_pid)
                    #log.debug(f'iipid = {iipid}')

                    nsm_n1 = (p1[iipid, :] + jhat[iipid, :] * m1a[jpid] + khat[iipid, :] * m2a[jpid])
                    nsm_n2 = (p2[iipid, :] + jhat[iipid, :] * m1b[jpid] + khat[iipid, :] * m2b[jpid])

                    #log.debug(f'  jhat={jhat}')
                    #log.debug(f'  khat={khat}')
                    #log.debug(f'  m1a={m1a}')
                    #log.debug(f'  m2a={m2a}')
                    #log.debug(f'  m1b={m1b}')
                    #log.debug(f'  m2b={m2b}')

                    #print("nsm_per_length=%s" % nsm_per_length)
                    #log.debug('  nsm_n1=%s' % nsm_n1)
                    #log.debug('  nsm_n2=%s' % nsm_n2)
                    nsm_centroid[iipid] = (nsm_n1 + nsm_n2) / 2.
                    mass_per_length[iipid] = mpl
                    nsm_per_length[iipid] = nsmpl
                    #nsm_centroid[iipid] = 0.

                #if nsm != 0.:
                    #p1_nsm = p1 + prop.ma
                    #p2_nsm = p2 + prop.mb
            elif prop.type == 'PBEAML':
                prop2 = prop.slice_card_by_property_id(pids_common)
                mass_per_lengths_ = prop2.mass_per_length()
                for jpid, pid, mpl in zip(count(), pids_common, mass_per_lengths_):
                    iipid = np.where(self.property_id == pid)[0]
                    if len(iipid) == 0:
                        log.warning(f'  skipping {prop.type}; pid={pid}')
                        continue
                    mass_per_length[iipid] = mpl

                    #mass_per_length = prop.MassPerLength() # includes simplified nsm
                    nsm_centroid[iipid, :] = (p1[iipid, :] + p2[iipid, :]) / 2.

                    # mass_per_length already includes nsm
                    nsm_per_length[iipid] = 0.

                #print('mass_per_lengths=%s nsm_per_lengths=%s' % (
                    #mass_per_lengths, nsm_per_lengths))
                #print('mass_per_length=%s nsm_per_length=%s' % (
                    #mass_per_length, nsm_per_length))
            elif prop.type == 'PBCOMP':
                prop2 = prop.slice_card_by_property_id(pids_common)
                mass_per_lengths_ = prop2.mass_per_length()
                m1s = prop2.m1
                m2s = prop2.m2
                for jpid, pid, mpl, m1, m2 in zip(count(), pids_common,
                                                  mass_per_lengths_, m1s, m2s):
                    iipid = np.where(self.property_id == pid)[0]
                    if len(iipid) == 0:
                        log.warning(f'  skipping {prop.type}; pid={pid}')
                        continue
                    mass_per_length[iipid] = mpl

                    # already accounted for in mass_per_length
                    nsm_per_length[iipid] = 0.0

                    nsm_n1 = (p1[iipid, :] + jhat[iipid, :] * m1 + khat[jpid, :] * m2)
                    nsm_n2 = (p2[iipid, :] + jhat[iipid, :] * m1 + khat[jpid, :] * m2)
                    nsm_centroid[iipid, :] = (nsm_n1 + nsm_n2) / 2.
            #elif prop.type == 'PBMSECT':
                #continue
                #mass_per_length = prop.MassPerLength()
                #m = mass_per_length * length
                #nsm = prop.nsm
            elif prop.type == 'PBMSECT':
                raise RuntimeError(prop)
                #mass_per_length = 0.  ## TODO: fix me
                #nsm_per_length = prop.nsm
                #nsm_centroid = (p1 + p2) / 2.
            else:
                raise NotImplementedError(prop.type)

            assert isinstance(mass_per_length, np.ndarray), mass_per_length
            assert isinstance(nsm_per_length, np.ndarray), nsm_per_length
            assert nsm_centroid.shape == (self.n, 3), nsm_centroid.shape

        if np.isnan(nsm_centroid.max()):
            inan = np.isnan(nsm_centroid.max(axis=1))
            assert len(inan) == len(self.element_id)
            #eid = self.element_id[inan]
            pid = self.property_id[inan]
            upid = np.unique(pid)
            #eid_pid = np.column_stack([eid, pid])
            #'[eid,pid]={eid_pid}\n'
            raise RuntimeError(f'nan nsm_centroids for upids={upid}')

        assert not np.isnan(mass_per_length.max()), mass_per_length
        assert not np.isnan(nsm_per_length.max()), nsm_per_length
        total_mass = mass_per_length + nsm_per_length

        assert centroid.shape == (self.n, 3), centroid.shape
        if np.abs(total_mass).sum() == 0.0:
            return centroid

        center_of_mass = (centroid * mass_per_length[:, np.newaxis] + nsm_centroid * nsm_per_length[:, np.newaxis]) / total_mass[:, np.newaxis]
        assert not np.isnan(center_of_mass.max()), center_of_mass
        assert mass_per_length.shape == (self.n, ), mass_per_length.shape
        assert nsm_per_length.shape == (self.n, ), nsm_per_length.shape
        assert total_mass.shape == (self.n, ), total_mass.shape
        assert nsm_centroid.shape == (self.n, 3), nsm_centroid.shape
        assert center_of_mass.shape == (self.n, 3), center_of_mass.shape
        return center_of_mass
        #return self.centroid()

    def area(self) -> np.ndarray:
        pid = self.property_id
        area = np.full(len(pid), np.nan, dtype='float64')
        log = self.model.log
        for prop in self.allowed_properties:
            i_lookup, i_all = searchsorted_filter(prop.property_id, pid, msg='')
            if len(i_lookup) == 0:
                continue
            # we're at least using some properties
            areai = prop.area()
            area_all = areai[i_all]
            inan = np.isnan(area_all)
            if np.any(inan):
                msg = f'{prop.type} has nan area for property_ids={prop.property_id[inan]}\n'
                log.warning(msg)

            area[i_lookup] = areai[i_all]

        inan = np.isnan(area)
        if np.any(inan):
            msg = 'CBEAM has nan area\n'
            msg += f'eids={self.element_id[inan]}\n'
            msg += f'pid={self.property_id[inan]}\n'
            msg += f'all_properties={self.all_properties}'
            #msg += f'As={self.nodes[inan]}\n'
            raise RuntimeError(msg)
        return area

    def volume(self) -> np.ndarray:
        A = self.area()
        L = self.length()
        return A * L

    @property
    def is_offt(self) -> np.ndarray:
        is_offt = (self.bit == -1)
        return is_offt

    @property
    def is_bit(self) -> np.ndarray:
        return not self.is_offt


class PBEAM(Property):
    """
    Defines the properties of a beam element (CBEAM entry). This element may be
    used to model tapered beams.


    +-------+-------+-------+-------+-------+-------+--------+-------+--------+
    | PBEAM |  PID  |  MID  | A(A)  | I1(A) | I2(A) | I12(A) | J(A)  | NSM(A) |
    +-------+-------+-------+-------+-------+-------+--------+-------+--------+
    |       | C1(A) | C2(A) | D1(A) | D2(A) | E1(A) | E2(A)  | F1(A) | F2(A)  |
    +-------+-------+-------+-------+-------+-------+--------+-------+--------+

    The next two continuations are repeated for each intermediate station as
    described in Remark 5. and SO and X/XB must be specified.

    +----+------+----+----+----+-----+----+-----+
    | SO | X/XB | A  | I1 | I2 | I12 | J  | NSM |
    +----+------+----+----+----+-----+----+-----+
    | C1 |  C2  | D1 | D2 | E1 | E2  | F1 | F2  |
    +----+------+----+----+----+-----+----+-----+

    The last two continuations are:
    +-------+-------+-------+-------+--------+--------+-------+-------+
    |   K1  |   K2  |   S1  |   S2  | NSI(A) | NSI(B) | CW(A) | CW(B) |
    +-------+-------+-------+-------+--------+--------+-------+-------+
    | M1(A) | M2(A) | M1(B) | M2(B) | N1(A)  | N2(A)  | N1(B) | N2(B) |
    +-------+-------+-------+-------+--------+--------+-------+-------+
    """
    @Property.clear_check
    def clear(self) -> None:
        self.property_id: np.array = np.array([], dtype='int32')
        self.material_id: np.array = np.array([], dtype='int32')
        self.nstation = np.array([], dtype='int32')

        self.s1 = np.array([], dtype='float64')
        self.s2 = np.array([], dtype='float64')
        self.k1 = np.array([], dtype='float64')
        self.k2 = np.array([], dtype='float64')

        self.nsia = np.array([], dtype='float64')
        self.nsib = np.array([], dtype='float64')
        self.cwa = np.array([], dtype='float64')
        self.cwb = np.array([], dtype='float64')

        self.m1a = np.array([], dtype='float64')
        self.m2a = np.array([], dtype='float64')
        self.m1b = np.array([], dtype='float64')
        self.m2b = np.array([], dtype='float64')
        self.n1a = np.array([], dtype='float64')
        self.n2a = np.array([], dtype='float64')
        self.n1b = np.array([], dtype='float64')
        self.n2b = np.array([], dtype='float64')

        self.xxb = np.array([], dtype='float64')
        self.so = np.array([], dtype='|U4')
        self.A = np.array([], dtype='float64')
        self.J = np.array([], dtype='float64')
        self.I1 = np.array([], dtype='float64')
        self.I2 = np.array([], dtype='float64')
        self.I12 = np.array([], dtype='float64')
        self.nsm = np.array([], dtype='float64')

        self.c1 = np.array([], dtype='float64')
        self.c2 = np.array([], dtype='float64')
        self.d1 = np.array([], dtype='float64')
        self.d2 = np.array([], dtype='float64')
        self.e1 = np.array([], dtype='float64')
        self.e2 = np.array([], dtype='float64')
        self.f1 = np.array([], dtype='float64')
        self.f2 = np.array([], dtype='float64')

    def add(self, pid, mid, xxb, so, area, i1, i2, i12, j, nsm=None,
            c1=None, c2=None, d1=None, d2=None,
            e1=None, e2=None, f1=None, f2=None,
            k1=1., k2=1., s1=0., s2=0.,
            nsia=0., nsib=None, cwa=0., cwb=None,
            m1a=0., m2a=0., m1b=None, m2b=None,
            n1a=0., n2a=0., n1b=None, n2b=None,
            comment='') -> int:
        self.cards.append((pid, mid,
                           xxb, so, area, j, i1, i2, i12, nsm,
                           c1, c2, d1, d2, e1, e2, f1, f2,
                           s1, s2, k1, k2,
                           nsia, nsib, cwa, cwb,
                           m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b,
                           comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'property_id')
        mid = integer(card, 2, 'material_id')

        area0 = double(card, 3, 'Area')
        i1a  = double_or_blank(card, 4, 'I1',  default=0.0)
        i2a  = double_or_blank(card, 5, 'I2',  default=0.0)
        i12a = double_or_blank(card, 6, 'I12', default=0.0)
        ja   = double_or_blank(card, 7, 'J',   default=0.0)
        nsma = double_or_blank(card, 8, 'nsm', default=0.0)
        area = [area0]
        i1 = [i1a]
        i2 = [i2a]
        i12 = [i12a]
        j = [ja]
        nsm = [nsma]

        assert area[0] >= 0., 'PBEAM pid=%s area=%s' % (pid, area)
        assert i1[0] >= 0., 'PBEAM pid=%s i1=%s' % (pid, i1)
        assert i2[0] >= 0., 'PBEAM pid=%s i2=%s' % (pid, i2)

        # we'll do a check for warping later; cwa/cwb -> j > 0.0
        assert j[0] >= 0., 'PBEAM pid=%s j=%s' % (pid, j)

        if i1a * i2a - i12a ** 2 <= 0.:
            msg = 'I1 * I2 - I12^2=0 and must be greater than 0.0 at End A\n'
            msg += 'i1=%s i2=%s i12=%s i1*i2-i12^2=%s'  % (i1a, i2a, i12a, i1a*i2a-i12a**2)
            raise ValueError(msg)

        # TODO: can you have a single lined PBEAM...I think so...
        # the second line is blank, so all values would be None (End A)
        # the NO xxb would be implicitly defined as 1.0
        # the End B values would try to use End A values, but because they're not set,
        # the defaults would get applied
        # the final 2 lines will default
        # finally, there would be no output at End A, but there would be output at End A.
        ifield = 9
        field9 = double_string_or_blank(card, 9, 'field9', default=0.0)
        if isinstance(field9, float):
            # C/D/E/F
            c1a = double_or_blank(card, 9,  'c1', default=0.0)
            c2a = double_or_blank(card, 10, 'c2', default=0.0)
            d1a = double_or_blank(card, 11, 'd1', default=0.0)
            d2a = double_or_blank(card, 12, 'd2', default=0.0)
            e1a = double_or_blank(card, 13, 'e1', default=0.0)
            e2a = double_or_blank(card, 14, 'e2', default=0.0)
            f1a = double_or_blank(card, 15, 'f1', default=0.0)
            f2a = double_or_blank(card, 16, 'f2', default=0.0)
            c1 = [c1a]
            c2 = [c2a]
            d1 = [d1a]
            d2 = [d2a]
            e1 = [e1a]
            e2 = [e2a]
            f1 = [f1a]
            f2 = [f2a]
            so = ['YES']
            ifield += 8 # 9 + 8 = 17
        else:
            c1a = c2a = d1a = d2a = e1a = e2a = f1a = f2a = 0.0
            c1 = [None]
            c2 = [None]
            d1 = [None]
            d2 = [None]
            e1 = [None]
            e2 = [None]
            f1 = [None]
            f2 = [None]
            so = ['NO']
            if field9 not in ['YES', 'YESA', 'NO']:
                msg = ('field9=%r on the PBEAM pid=%s must be [YES, YESA, NO] '
                       'because C/D/E/F at A is not specified' % field9)
                raise ValueError(field9)
        xxb = [0.]

        # --------------------------------------------------------
        irow = 0
        nrows_max = 10
        for irow in range(nrows_max):
            nrepeated = irow + 1
            SOi_k1 = double_string_or_blank(card, ifield, 'SO_%d/K1' % nrepeated)
            if isinstance(SOi_k1, float) or SOi_k1 is None:
                # we found K1
                break
            else:
                soi = string(card, ifield, 'SO%i' % nrepeated)
                xxbi = double(card, ifield + 1, 'x/xb%i' % nrepeated)
                if xxbi == 1.0:
                    # these have already been checked such that they're greater than 0
                    # so when we interpolate, our values will be correct
                    areai = double_or_blank(card, ifield + 2, 'Area%d' % nrepeated, default=area0)
                    i1i   = double_or_blank(card, ifield + 3, 'I1 %d'  % nrepeated, default=i1a)
                    i2i   = double_or_blank(card, ifield + 4, 'I2 %d'  % nrepeated, default=i2a)
                    i12i  = double_or_blank(card, ifield + 5, 'I12 %d' % nrepeated, default=i12a)
                    ji    = double_or_blank(card, ifield + 6, 'J%i' % nrepeated, ja)
                    nsmi  = double_or_blank(card, ifield + 7, 'nsm%i' % nrepeated, nsma)

                    assert areai >= 0., areai
                    assert i1i >= 0., i1i
                    assert i2i >= 0., i2i
                    assert ji >= 0., ji

                    # we'll do a check for warping later; cwa/cwb -> j > 0.0
                    #assert j[-1] >= 0., j
                    if i1i * i2i - i12i ** 2 <= 0.:
                        msg = 'I1 * I2 - I12^2=0 and must be greater than 0.0 at End B\n'
                        msg += 'xxb=1.0 i1=%s i2=%s i12=%s'  % (i1i, i2i, i12i)
                        raise ValueError(msg)
                else:
                    ## 0.0 < x/xb < 1.0
                    # we'll go through and do linear interpolation afterwards
                    areai = double_or_blank(card, ifield + 2, 'Area%d' % nrepeated, default=np.nan)
                    i1i   = double_or_blank(card, ifield + 3, 'I1 %d'  % nrepeated, default=np.nan)
                    i2i   = double_or_blank(card, ifield + 4, 'I2 %d'  % nrepeated, default=np.nan)
                    i12i  = double_or_blank(card, ifield + 5, 'I12 %d' % nrepeated, default=np.nan)
                    ji    = double_or_blank(card, ifield + 6, 'J%d'    % nrepeated, default=np.nan)
                    nsmi  = double_or_blank(card, ifield + 7, 'nsm%d'  % nrepeated, default=np.nan)
                    #assert areai >= 0., areai
                    #assert i1i >= 0., i1i
                    #assert i2i >= 0., i2i

                so.append(soi)
                xxb.append(xxbi)
                area.append(areai)
                i1.append(i1i)
                i2.append(i2i)
                i12.append(i12i)
                j.append(ji)
                nsm.append(nsmi)

                if soi == 'YES':
                    c1i = double_or_blank(card, ifield + 8,  'c1 %d' % nrepeated, default=np.nan)
                    c2i = double_or_blank(card, ifield + 9,  'c2 %d' % nrepeated, default=np.nan)
                    d1i = double_or_blank(card, ifield + 10, 'd1 %d' % nrepeated, default=np.nan)
                    d2i = double_or_blank(card, ifield + 11, 'd2 %d' % nrepeated, default=np.nan)
                    e1i = double_or_blank(card, ifield + 12, 'e1 %d' % nrepeated, default=np.nan)
                    e2i = double_or_blank(card, ifield + 13, 'e2 %d' % nrepeated, default=np.nan)
                    f1i = double_or_blank(card, ifield + 14, 'f1 %d' % nrepeated, default=np.nan)
                    f2i = double_or_blank(card, ifield + 15, 'f2 %d' % nrepeated, default=np.nan)
                    ifield += 16
                elif soi == 'YESA':
                    c1i = c1a
                    c2i = c2a
                    d1i = d1a
                    d2i = d2a
                    e1i = e1a
                    e2i = e2a
                    f1i = f1a
                    f2i = f2a
                    ifield += 8
                elif soi == 'NO':
                    c1i = None
                    c2i = None
                    d1i = None
                    d2i = None
                    e1i = None
                    e2i = None
                    f1i = None
                    f2i = None
                    ifield += 8
                else:
                    raise RuntimeError(f'PBEAM: pid={pid} so={soi!r} and not [YES, YESA, NO]')
                c1.append(c1i)
                c2.append(c2i)
                d1.append(d1i)
                d2.append(d2i)
                e1.append(e1i)
                e2.append(e2i)
                f1.append(f1i)
                f2.append(f2i)
        if irow != 0:
            assert min(xxb) == 0.0, 'pid=%s x/xb=%s' % (pid, xxb)
            assert max(xxb) == 1.0, 'pid=%s x/xb=%s' % (pid, xxb)
            assert len(xxb) == len(np.unique(xxb)), xxb

        # calculate:
        #    k1, k2, s1, s2
        #    m1a, m2a, n1a, n2a, etc.

        # footer fields
        #: Shear stiffness factor K in K*A*G for plane 1.
        k1 = double_or_blank(card, ifield, 'k1', default=1.0)
        #: Shear stiffness factor K in K*A*G for plane 2.
        k2 = double_or_blank(card, ifield + 1, 'k2', default=1.0)

        #: Shear relief coefficient due to taper for plane 1.
        s1 = double_or_blank(card, ifield + 2, 's1', default=0.0)
        #: Shear relief coefficient due to taper for plane 2.
        s2 = double_or_blank(card, ifield + 3, 's2', default=0.0)

        #: non structural mass moment of inertia per unit length
        #: about nsm center of gravity at Point A.
        nsia = double_or_blank(card, ifield + 4, 'nsia', default=0.0)
        #: non structural mass moment of inertia per unit length
        #: about nsm center of gravity at Point B.
        nsib = double_or_blank(card, ifield + 5, 'nsib', default=nsia)

        #: warping coefficient for end A.
        cwa = double_or_blank(card, ifield + 6, 'cwa', default=0.0)
        #: warping coefficient for end B.
        cwb = double_or_blank(card, ifield + 7, 'cwb', default=cwa)

        #: y coordinate of center of gravity of
        #: nonstructural mass for end A.
        m1a = double_or_blank(card, ifield + 8, 'm1a', default=0.0)
        #: z coordinate of center of gravity of
        #: nonstructural mass for end A.
        m2a = double_or_blank(card, ifield + 9, 'm2a', default=0.0)

        #: y coordinate of center of gravity of
        #: nonstructural mass for end B.
        m1b = double_or_blank(card, ifield + 10, 'm1b', default=m1a)
        #: z coordinate of center of gravity of
        #: nonstructural mass for end B.
        m2b = double_or_blank(card, ifield + 11, 'm2b', default=m2a)

        #: y coordinate of neutral axis for end A.
        n1a = double_or_blank(card, ifield + 12, 'n1a', default=0.0)
        #: z coordinate of neutral axis for end A.
        n2a = double_or_blank(card, ifield + 13, 'n2a', default=0.0)


        #: y coordinate of neutral axis for end B.
        n1b = double_or_blank(card, ifield + 14, 'n1a', default=n1a)
        #: z coordinate of neutral axis for end B.
        n2b = double_or_blank(card, ifield + 15, 'n2b', default=n2a)


        ifield += 16
        if len(card) > ifield:
            msg = 'len(card)=%s is too long; max=%s\n' % (len(card), ifield)
            msg += 'You probably have a empty line after the YESA/NO line.\n'
            msg += 'The next line must have K1.\n'
            msg += 'pid = %s\n' % pid
            msg += 'mid = %s\n' % mid
            msg += 's0 = %s\n' % so
            msg += 'xxb = %s\n' % xxb

            msg += 'A = %s\n' % area
            msg += 'i1 = %s\n' % i1
            msg += 'i2 = %s\n' % i2
            msg += 'i12 = %s\n' % i12
            msg += 'j = %s\n' % j
            msg += 'nsm = %s\n\n' % nsm

            msg += 'c1 = %s\n' % c1
            msg += 'c2 = %s\n' % c2
            msg += 'd1 = %s\n' % d1
            msg += 'd2 = %s\n' % d2
            msg += 'e1 = %s\n' % e1
            msg += 'e2 = %s\n' % e2
            msg += 'f1 = %s\n' % f1
            msg += 'f2 = %s\n\n' % f2

            msg += 'k1 = %s\n' % k1
            msg += 'k2 = %s\n' % k2
            msg += 's1 = %s\n' % s1
            msg += 's2 = %s\n' % s2
            msg += 'nsia = %s\n' % nsia
            msg += 'nsib = %s\n\n' % nsib

            msg += 'cwa = %s\n' % cwa
            msg += 'cwb = %s\n' % cwb
            msg += 'm1a = %s\n' % m1a
            msg += 'm2a = %s\n' % m2a
            msg += 'mb1 = %s\n' % m1b
            msg += 'm2b = %s\n' % m2b
            msg += 'n1a = %s\n' % n1a
            msg += 'n2a = %s\n' % n2a
            msg += 'n1b = %s\n' % n1b
            msg += 'n2b = %s\n' % n2b
            raise RuntimeError(msg)


        # ------------------------------
        # we now need to sort the xxb values
        out = _sort_pbeam(
            pid, xxb, so, area, i1, i2, i12, j, nsm,
            c1, c2, d1, d2, e1, e2, f1, f2)
        xxb, so, area, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2, f1, f2 = out

        # fill in any values of 0.0 on area, i1, ..., with the linearly interpolated value
        out = _interpolate_pbeam_sections(
            pid, xxb, so, area, i1, i2, i12, j, nsm,
            c1, c2, d1, d2, e1, e2, f1, f2,)
        xxb, so, area, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2, f1, f2 = out

        if cwa or cwb:  # if either is non-zero
            for i, xxbi, ji in zip(count(), xxb, j):
                if ji < 0.:
                    msg = 'Warping Check Error; j[%i] must be greater than 0.0' % i
                    msg += '  cwa=%s cwb=%s\n' % (cwa, cwb)
                    msg += '  i=%s xxb=%s j=%s; j[%i]=%s\n' % (i, xxbi, j, i, ji)
                    raise ValueError(msg)

        assert nsm is not None, nsm
        self.cards.append((pid, mid,
                           xxb, so, area, j, i1, i2, i12, nsm,
                           c1, c2, d1, d2, e1, e2, f1, f2,
                           s1, s2, k1, k2,
                           nsia, nsib, cwa, cwb,
                           m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b,
                           comment))

        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        property_id = np.zeros(ncards, dtype=idtype)
        material_id = np.zeros(ncards, dtype=idtype)
        xxb_list = []
        so_list = []
        A_list = []
        J_list = []
        I1_list = []
        I2_list = []
        I12_list = []
        nsm_list = []
        c1_list = []
        c2_list = []
        d1_list = []
        d2_list = []
        e1_list = []
        e2_list = []
        f1_list = []
        f2_list = []
        nstation = np.zeros(ncards, dtype='int32')

        s1 = np.zeros(ncards, dtype='float64')
        s2 = np.zeros(ncards, dtype='float64')
        k1 = np.zeros(ncards, dtype='float64')
        k2 = np.zeros(ncards, dtype='float64')

        nsia = np.zeros(ncards, dtype='float64')
        nsib = np.zeros(ncards, dtype='float64')
        cwa = np.zeros(ncards, dtype='float64')
        cwb = np.zeros(ncards, dtype='float64')

        m1a = np.zeros(ncards, dtype='float64')
        m2a = np.zeros(ncards, dtype='float64')
        m1b = np.zeros(ncards, dtype='float64')
        m2b = np.zeros(ncards, dtype='float64')
        n1a = np.zeros(ncards, dtype='float64')
        n2a = np.zeros(ncards, dtype='float64')
        n1b = np.zeros(ncards, dtype='float64')
        n2b = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, mid,
             xxbi, soi, areai, ji, i1i, i2i, i12i, nsmi,
             c1i, c2i, d1i, d2i, e1i, e2i, f1i, f2i,
             s1i, s2i, k1i, k2i,
             nsiai, nsibi, cwai, cwbi,
             m1ai, m2ai, m1bi, m2bi,
             n1ai, n2ai, n1bi, n2bi,
             comment) = card
            nstations = len(areai)
            if nsmi is None:
                nsmi = np.zeros(nstations)
            if c1i is None:
                c1i = np.zeros(nstations)
            if c2i is None:
                c2i = np.zeros(nstations)
            if d1i is None:
                d1i = np.zeros(nstations)
            if d2i is None:
                d2i = np.zeros(nstations)
            if e1i is None:
                e1i = np.zeros(nstations)
            if e2i is None:
                e2i = np.zeros(nstations)
            if f1i is None:
                f1i = np.zeros(nstations)
            if f2i is None:
                f2i = np.zeros(nstations)
            #if nsmi is None:
                #nsmi = np.zeros(nstations)

            if m1ai is None:
                m1ai = 0.0
            if m2ai is None:
                m2ai = 0.0

            if m1bi is None:
                m1bi = m1ai
            if m2bi is None:
                m2bi = m2ai

            if n1ai is None:
                n1ai = 0.0
            if n2ai is None:
                n2ai = 0.0

            if n1bi is None:
                n1bi = n1ai
            if n2bi is None:
                n2bi = n2ai

            nstationi = len(xxbi)
            property_id[icard] = pid
            material_id[icard] = mid
            nstation[icard] = nstationi

            s1[icard] = s1i
            s2[icard] = s2i
            k1[icard] = k1i
            k2[icard] = k2i

            xxb_list.extend(xxbi)
            so_list.extend(soi)
            A_list.extend(areai)
            J_list.extend(ji)
            I1_list.extend(i1i)
            I2_list.extend(i2i)
            I12_list.extend(i12i)
            nsm_list.extend(nsmi)

            c1_list.extend(c1i)
            c2_list.extend(c2i)
            d1_list.extend(d1i)
            d2_list.extend(d2i)
            e1_list.extend(e1i)
            e2_list.extend(e2i)
            f1_list.extend(f1i)
            f2_list.extend(f2i)

            nsia[icard] = nsiai
            nsib[icard] = nsibi
            cwa[icard] = cwai
            cwb[icard] = cwbi

            m1a[icard] = m1ai
            m2a[icard] = m2ai
            m1b[icard] = m1bi
            m2b[icard] = m2bi
            n1a[icard] = n1ai
            n2a[icard] = n2ai
            n1b[icard] = n1bi
            n2b[icard] = n2bi

            #self.area[istation0:istation1] = area
            #return PBEAM(
                #pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                #c1, c2, d1, d2, e1, e2, f1, f2,
                #k1, k2, s1, s2,
                #nsia, nsib, cwa, cwb, m1a,
                #m2a, m1b, m2b, n1a, n2a, n1b, n2b,
                #comment=comment)
        xxb = np.array(xxb_list, dtype='float64')
        so = np.array(so_list, dtype='|U4')
        A = np.array(A_list, dtype='float64')
        J = np.array(J_list, dtype='float64')
        I1 = np.array(I1_list, dtype='float64')
        I2 = np.array(I2_list, dtype='float64')
        I12 = np.array(I12_list, dtype='float64')
        nsm = np.array(nsm_list, dtype='float64')

        c1 = np.array(c1_list, dtype='float64')
        c2 = np.array(c2_list, dtype='float64')
        d1 = np.array(d1_list, dtype='float64')
        d2 = np.array(d2_list, dtype='float64')
        e1 = np.array(e1_list, dtype='float64')
        e2 = np.array(e2_list, dtype='float64')
        f1 = np.array(f1_list, dtype='float64')
        f2 = np.array(f2_list, dtype='float64')
        self._save(property_id, material_id,
                   nstation, xxb, so,
                   A, J, I1, I2, I12, nsm,
                   c1, c2, d1, d2, e1, e2, f1, f2,
                   s1, s2, k1, k2,
                   nsia, nsib, cwa, cwb,
                   m1a, m2a, m1b, m2b,
                   n1a, n2a, n1b, n2b)
        self.sort()

    def _save(self, property_id, material_id,
              nstation, xxb, so,
              A, J, I1, I2, I12, nsm,
              c1, c2, d1, d2, e1, e2, f1, f2,
              s1, s2, k1, k2,
              nsia, nsib, cwa, cwb,
              m1a, m2a, m1b, m2b,
              n1a, n2a, n1b, n2b) -> None:
        self.property_id = property_id
        self.material_id = material_id

        self.nstation = nstation
        self.xxb = xxb
        self.so = so
        self.A = A
        self.J = J
        self.I1 = I1
        self.I2 = I2
        self.I12 = I12
        self._nsm = nsm

        self.c1 = c1
        self.c2 = c2
        self.d1 = d1
        self.d2 = d2
        self.e1 = e1
        self.e2 = e2
        self.f1 = f1
        self.f2 = f2

        self.s1 = s1
        self.s2 = s2
        self.k1 = k1
        self.k2 = k2

        self.nsia = nsia
        self.nsib = nsib
        self.cwa = cwa
        self.cwb = cwb
        self.m1a = m1a
        self.m2a = m2a
        self.m1b = m1b
        self.m2b = m2b
        self.n1a = n1a
        self.n2a = n2a
        self.n1b = n1b
        self.n2b = n2b
        self.n = len(property_id)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    def convert(self, xyz_scale: float=1.0,
                area_scale: float=1.0,
                area_inertia_scale:float=1.0,
                nsm_per_length_scale: float=1.0, **kwargs):
        # easy
        #self.xxb # percent of length
        self.A *= area_scale
        self.J *= area_inertia_scale
        self.I1 *= area_inertia_scale
        self.I2 *= area_inertia_scale
        self.I12 *= area_inertia_scale

        self.c1 *= xyz_scale
        self.c2 *= xyz_scale
        self.d1 *= xyz_scale
        self.d2 *= xyz_scale
        self.e1 *= xyz_scale
        self.e2 *= xyz_scale
        self.f1 *= xyz_scale
        self.f2 *= xyz_scale
        self._nsm *= nsm_per_length_scale

        # ???
        self.m1a *= xyz_scale
        self.m1b *= xyz_scale
        self.n2a *= xyz_scale
        self.n2b *= xyz_scale

        #self.k1 = k1
        #self.k2 = k2

        #self.nsia = nsia
        #self.nsib = nsib
        #self.cwa = cwa
        #self.cwb = cwb

    def __apply_slice__(self, prop: PBEAM, i: np.ndarray) -> None:
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]

        prop.s1 = self.s1[i]
        prop.s2 = self.s2[i]
        prop.k1 = self.k1[i]
        prop.k2 = self.k2[i]

        prop.nsia = self.nsia[i]
        prop.nsib = self.nsib[i]
        prop.cwa = self.cwa[i]
        prop.cwb = self.cwb[i]
        prop.m1a = self.m1a[i]
        prop.m2a = self.m2a[i]
        prop.m1b = self.m1b[i]
        prop.m2b = self.m2b[i]
        prop.n1a = self.n1a[i]
        prop.n2a = self.n2a[i]
        prop.n1b = self.n1b[i]
        prop.n2b = self.n2b[i]

        istation = self.istation
        prop.xxb = hslice_by_idim(i, istation, self.xxb)
        prop.so = hslice_by_idim(i, istation, self.so)
        prop.A = hslice_by_idim(i, istation, self.A)
        prop.J = hslice_by_idim(i, istation, self.J)
        prop.I1 = hslice_by_idim(i, istation, self.I1)
        prop.I2 = hslice_by_idim(i, istation, self.I2)
        prop.I12 = hslice_by_idim(i, istation, self.I12)
        prop._nsm = hslice_by_idim(i, istation, self._nsm)

        prop.c1 = hslice_by_idim(i, istation, self.c1)
        prop.c2 = hslice_by_idim(i, istation, self.c2)
        prop.d1 = hslice_by_idim(i, istation, self.d1)
        prop.d2 = hslice_by_idim(i, istation, self.d2)
        prop.e1 = hslice_by_idim(i, istation, self.e1)
        prop.e2 = hslice_by_idim(i, istation, self.e2)
        prop.f1 = hslice_by_idim(i, istation, self.f1)
        prop.f2 = hslice_by_idim(i, istation, self.f2)
        prop.nstation = self.nstation[i]
        prop.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        materials = self.allowed_materials
        mids = hstack_msg([prop.material_id for prop in materials],
                          msg=f'no materials for {self.type}; {self.all_materials}')
        mids = np.unique(mids)
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    @property
    def all_materials(self) -> list[MAT1]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[MAT1]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def rho(self) -> np.ndarray:
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        return rho

    def mass_per_length(self) -> np.ndarray:
        nproperties = len(self.property_id)
        rho = self.rho()
        mass_per_length = np.zeros(nproperties, dtype='float64')
        if rho.max() == 0. and rho.min() == 0. and self._nsm.max() == 0. and self._nsm.min() == 0.:
            return mass_per_length

        for i, rhoi, istation in zip(count(), rho, self.istation):
            istation0, istation1 = istation
            assert istation1 > istation0
            xxb = self.xxb[istation0:istation1]
            area = self.A[istation0:istation1]
            nsm = self._nsm[istation0:istation1]
            mass_per_lengths = rhoi * area + nsm
            mass_per_lengthi = integrate_positive_unit_line(xxb, mass_per_lengths)
            mass_per_length[i] = mass_per_lengthi

        assert len(mass_per_length) == nproperties
        return mass_per_length

    def rhoarea_nsm(self) -> np.ndarray:
        nproperties = len(self.property_id)
        rho = self.rho()

        nsm_per_length = np.zeros(nproperties, dtype='float64')
        rho_area = np.zeros(nproperties, dtype='float64')
        if rho.max() == 0. and rho.min() == 0. and self._nsm.max() == 0. and self._nsm.min() == 0.:
            return rho_area, nsm_per_length

        #assert len(rho) == len(self.A)
        assert len(rho) == len(nsm_per_length)
        for i, istation in zip(count(), self.istation):
            istation0, istation1 = istation
            assert istation1 > istation0
            areai = self.A[istation0:istation1]
            xxbi = self.xxb[istation0:istation1]
            nsmi = self._nsm[istation0:istation1]
            nsm_per_length[i] = integrate_positive_unit_line(xxbi, nsmi)

            rho_areas = rho[i] * areai # + nsmi
            #mass_per_lengths = rho[i] * areai + nsmi
            rho_area[i] = integrate_positive_unit_line(xxbi, rho_areas)
            #mass_per_length[i] = integrate_positive_unit_line(xxbi, mass_per_lengths)
            #rho_area[i] = mass_per_lengthi

        assert len(nsm_per_length) == nproperties
        return rho_area, nsm_per_length

    @property
    def is_small_field(self) -> bool:
        return max(self.property_id.max(),
                   self.material_id.max(),
                   self.s1.max(), self.s2.max()) < 99_999_999

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:

        if size == 8 and self.is_small_field:
            print_card = print_card_8
        else:
            print_card = print_card_16

        for pid, mid, (istation0, istation1), nstation, \
            s1, s2, k1, k2, \
            nsia, nsib, cwa, cwb, \
            m1a, m2a, m1b, m2b, \
            n1a, n2a, n1b, n2b in zip(self.property_id, self.material_id,
                                      self.istation, self.nstation,
                                      self.s1, self.s2, self.k1, self.k2,
                                      self.nsia, self.nsib, self.cwa, self.cwb,
                                      self.m1a, self.m2a, self.m1b, self.m2b,
                                      self.n1a, self.n2a, self.n1b, self.n2b,):

            xxb_ = self.xxb[istation0:istation1]
            so_ = self.so[istation0:istation1].tolist()  # need to fix the type
            A_ = self.A[istation0:istation1]
            j_ = self.J[istation0:istation1]
            i1_ = self.I1[istation0:istation1]
            i2_ = self.I2[istation0:istation1]
            i12_ = self.I12[istation0:istation1]
            nsm_ = self._nsm[istation0:istation1]
            c1_ = self.c1[istation0:istation1]
            c2_ = self.c2[istation0:istation1]
            d1_ = self.d1[istation0:istation1]
            d2_ = self.d2[istation0:istation1]
            e1_ = self.e1[istation0:istation1]
            e2_ = self.e2[istation0:istation1]
            f1_ = self.f1[istation0:istation1]
            f2_ = self.f2[istation0:istation1]
            assert nstation > 0

            write_so = True
            if len(xxb_) == 2 and xxb_[0] == 0. and xxb_[1] == 1.0:
                is_same = all([value[0] == value[1] for value in (
                    A_, j_, i1_, i2_, i12_, nsm_,
                    c1_, c2_, d1_, d2_, e1_, e2_, f1_, f2_)])
                write_so = not is_same

            # still need to save these
            #nsia = nsib = 0.
            #cwa = cwb = 0.
            #k1 = k2 = 1.0
            #s1 = s2 = 0.0
            #m1a = m1b = 0.
            #m2a = m2b = 0.
            #n1a = n1b = 0.
            #n2a = n2b = 0.
            list_fields = ['PBEAM', pid, mid]

            i = 0
            for (so, xxb, A, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2, f1,
                 f2) in zip_longest(so_, xxb_, A_, i1_, i2_, i12_,
                            j_, nsm_, c1_, c2_, d1_, d2_,
                            e1_, e2_, f1_, f2_):
                if not self.write_default_fields:
                    i1 = set_blank_if_default(i1, 0.0)
                    i2 = set_blank_if_default(i2, 0.0)
                    i12 = set_blank_if_default(i12, 0.0)
                    j = set_blank_if_default(j, 0.0)

                    nsm = set_blank_if_default(nsm, 0.0)
                    c1 = set_blank_if_default(c1, 0.0)
                    d1 = set_blank_if_default(d1, 0.0)
                    e1 = set_blank_if_default(e1, 0.0)
                    f1 = set_blank_if_default(f1, 0.0)

                    c2 = set_blank_if_default(c2, 0.0)
                    d2 = set_blank_if_default(d2, 0.0)
                    e2 = set_blank_if_default(e2, 0.0)
                    f2 = set_blank_if_default(f2, 0.0)

                if i == 0:  # the first 2 fields aren't written
                    list_fields += [A, i1, i2, i12, j, nsm,
                                    c1, c2, d1, d2, e1, e2, f1, f2]
                else:
                    if so in ['YES']:
                        list_fields += ['YES', xxb, A, i1, i2, i12, j, nsm,
                                        c1, c2, d1, d2, e1, e2, f1, f2]
                    elif so in ['NO']:
                        list_fields += ['NO', xxb, A, i1, i2, i12, j, nsm]
                    elif so in ['YESA']:
                        list_fields += ['YESA', xxb, A, i1, i2, i12, j, nsm]
                    else:
                        raise RuntimeError('so=%r type(so)=%s' % (so, type(so)))
                i += 1
                if not write_so:
                    break

            if not self.write_default_fields:
                k1 = set_blank_if_default(k1, 1.0)
                k2 = set_blank_if_default(k2, 1.0)

                s1 = set_blank_if_default(s1, 0.0)
                s2 = set_blank_if_default(s2, 0.0)
                #k1 = self.k1
                #k2 = self.k2
                #s1 = self.s1
                #s2 = self.s2
                nsib = set_blank_if_default(nsib, nsia)
                nsia = set_blank_if_default(nsia, 0.0)

                cwb = set_blank_if_default(cwb, cwa)
                cwa = set_blank_if_default(cwa, 0.0)

                #m1a = self.m1a
                #m2a = self.m2a
                #m1b = self.m1b
                #m2b = self.m2b
                # Point A/B
                # Directions 1/2
                m1b = set_blank_if_default(m1b, m1a)
                m2b = set_blank_if_default(m2b, m2a)
                m1a = set_blank_if_default(m1a, 0.0)
                m2a = set_blank_if_default(m2a, 0.0)

                n1b = set_blank_if_default(n1b, n1a)
                n2b = set_blank_if_default(n2b, n2b)
                n1a = set_blank_if_default(n1a, 0.0)
                n2a = set_blank_if_default(n2a, 0.0)

            footer = [k1, k2, s1, s2, nsia, nsib, cwa, cwb,
                      m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b]
            if footer != [None] * len(footer):
                list_fields += footer
            bdf_file.write(print_card(list_fields))
        return

    @property
    def istation(self) -> np.ndarray:
        return make_idim(self.n, self.nstation)

    @property
    def k(self) -> np.ndarray:
        return np.column_stack([self.k1, self.k2])

    def e_g_nu(self) -> np.ndarray:
        """calculates E, G, nu"""
        e_g_nu = e_g_nu_from_isotropic_material(self.material_id, self.allowed_materials)
        return e_g_nu

    def area(self) -> np.ndarray:
        nproperties = len(self.property_id)
        areas = np.zeros(nproperties, dtype='float64')
        for i, istation in zip(count(), self.istation):
            istation0, istation1 = istation
            assert istation1 > istation0
            xxb = self.xxb[istation0:istation1]
            area = self.A[istation0:istation1]
            areasi = integrate_positive_unit_line(xxb, area)
            areas[i] = areasi
        assert len(areas) == nproperties
        return areas

    def inertia(self) -> np.ndarray:
        """i1, i2, i12, j"""
        nproperties = len(self.property_id)
        inertias = np.zeros((nproperties, 4), dtype='float64')
        for i, istation in zip(count(), self.istation):
            istation0, istation1 = istation
            assert istation1 > istation0
            xxb = self.xxb[istation0:istation1]
            i1si = self.I1[istation0:istation1]
            i2si = self.I2[istation0:istation1]
            i12si = self.I12[istation0:istation1]
            jsi = self.J[istation0:istation1]

            i1i = integrate_positive_unit_line(xxb, i1si)
            i2i = integrate_positive_unit_line(xxb, i2si)
            i12i = integrate_positive_unit_line(xxb, i12si)
            ji = integrate_positive_unit_line(xxb, jsi)

            inertias[i, :] = [i1i, i2i, i12i, ji]
        assert len(inertias) == nproperties
        return inertias

    def to_old_card(self) -> list[Any]:
        from pyNastran.bdf.bdf import BDF
        model = BDF()
        card_name = self.type
        cards = []
        dicti = model.properties
        for icard in range(self.n):
            card_obj = self.slice_card_by_index(icard)
            card_lines = card_obj.write(size=16).split('\n')
            unused_bdf_card = model.add_card(
                card_lines, card_name, comment='',
                ifile=None, is_list=False, has_none=False)
            eid = self.property_id[icard]
            card = dicti[eid]
            cards.append(card)
        return cards


class PBEAML(Property):
    valid_types = {
        "ROD": 1,
        "TUBE": 2,
        "TUBE2": 2,
        "L": 4,
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
        "DBOX": 10,  # TODO: was 12???
    }  # for GROUP="MSCBML0"

    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.Type = np.array([], dtype='|U8')
        self.group = np.array([], dtype='|U8')
        self._nsm = np.array([], dtype='float64')
        self.dims = np.zeros([], dtype='float64')
        self.ndim = np.array([], dtype='int32')
        self.nstation = np.array([], dtype='int32')

    def slice_card_by_property_id(self, property_id: np.ndarray) -> PBEAML:
        """uses a node_ids to extract GRIDs"""
        iprop = self.index(property_id)
        #assert len(self.node_id) > 0, self.node_id
        #i = np.searchsorted(self.node_id, node_id)
        prop = self.slice_card_by_index(iprop)
        return prop

    def add(self, pid: int, mid: int, beam_type: str,
            xxb: list[float], dims: list[list[float]],
            so=None, nsm=None,
            group: str='MSCBML0', comment: str='') -> int:
        nxxb = len(xxb)
        self.model.log.info(f'pid={pid} so0={so} beam_type={beam_type!r} dims={dims} nsm0={nsm} xxb={xxb}')
        if so is None:
            so = ['YES'] * nxxb
        elif isinstance(so, str):
            so = [so] * nxxb

        if nsm is None:
            nsm = [0.] * nxxb
        elif isinstance(nsm, float_types):
            nsm = [nsm] * nxxb
        ndim = self.valid_types[beam_type]

        nstation = len(xxb)
        self.model.log.info(f'  nstation={nstation} ndim={ndim}')
        assert nstation == len(dims), f'pid={pid} nstation={nstation} dims={dims}'
        assert nstation == len(nsm), f'pid={pid} nstation={nstation} nsm={nsm}'
        self.cards.append((pid, mid, beam_type, group, xxb, so, nsm, ndim, dims, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        group = string_or_blank(card, 3, 'group', default='MSCBML0')
        beam_type = string(card, 4, 'Type')

        # determine the number of required dimensions on the PBEAM
        ndim = self.valid_types[beam_type]

        #: dimension list
        dims: list[list[float]] = []
        dim: list[float] = []

        #: Section position
        xxb: list[float] = [0.]

        #: Output flag
        so: list[str] = ['YES']  # station 0

        #: non-structural mass :math:`nsm`
        nsm: list[float] = []

        ioffset = 9
        n = 0

        #n_so = (len(card) - 9) // (ndim + 2) #- 1
        #n_extra = (len(card) - 9) % (ndim + 2)
        xxbi = 0.0
        while ioffset < len(card):
            if n > 0:
                soi = string_or_blank(card, ioffset, 'so_n=%d' % n, default='YES')
                xxbi = double_or_blank(card, ioffset + 1, 'xxb_n=%d' % n, default=1.0)
                so.append(soi)
                xxb.append(xxbi)
                ioffset += 2

            # PBARL
            # 9. For DBOX section, the default value for DIM5 to DIM10 are
            #    based on the following rules:
            #     a. DIM5, DIM6, DIM7 and DIM8 have a default value of
            #        DIM4if not provided.
            #     b. DIM9 and DIM10 have a default value of DIM6 if not
            #        provided.

            #If any of the fields NSM(B), DIMi(B) are blank on the
            #continuation entry for End B, the values are set to the
            #values given for end A. For the continuation entries that
            #have values of X(j)/XB between 0.0 and 1.0 and use the
            #default option (blank field), a linear interpolation between
            #the values at ends A and B is performed to obtain the
            #missing field.
            dim = []
            if beam_type == 'DBOX':
                for ii in range(ndim):
                    field_name = 'istation=%s; ndim=%s; dim%i' % (n, ndim, ii+1)
                    if ii in [4, 5, 6, 7]:
                        dim4 = dim[3]
                        dimi = double_or_blank(card, ioffset, field_name, default=dim4)
                    elif ii in [8, 9]:
                        dim6 = dim[5]
                        dimi = double_or_blank(card, ioffset, field_name, default=dim6)
                    else:
                        dimi = double(card, ioffset, field_name)
                    dim.append(dimi)
                    ioffset += 1
            else:
                for ii in range(ndim):
                    field_name = 'istation=%s; ndim=%s; dim%d' % (n, ndim, ii+1)
                    if xxbi == 0.0:
                        dimi = double(card, ioffset, field_name)
                    elif xxbi == 1.0:
                        dims0 = dims[0]
                        dimi = double_or_blank(card, ioffset, field_name, dims0[ii])
                    else:
                        ## TODO: use linear interpolation
                        dimi = double(card, ioffset, field_name)

                    dim.append(dimi)
                    ioffset += 1
            dims.append(dim)

            nsmi = double_or_blank(card, ioffset, 'nsm_n=%d' % n, 0.0)
            nsm.append(nsmi)
            n += 1
            ioffset += 1
        assert len(card) > 5, card

        nstationi = len(xxb)
        assert nstationi == len(nsm), f'pid={pid} nstation={nstationi} nsm={nsm}'
        self.cards.append((pid, mid, beam_type, group, xxb, so, nsm, ndim, dims, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        property_id = np.zeros(ncards, dtype=idtype)
        material_id = np.zeros(ncards, dtype=idtype)

        #idim = np.zeros((ncards, 2), dtype='int32')  # for all properties
        ndim = np.zeros(ncards, dtype='int32')

        #istation = np.zeros((ncards, 2), dtype='int32')  # for all properties
        nstation = np.zeros(ncards, dtype='int32')

        Type = np.full(ncards, '', dtype='|U8')
        group = np.full(ncards, '', dtype='|U8')
        #nsm = np.zeros(ncards, dtype='float64')

        all_xxb = []
        all_dims = []
        all_so = []
        all_nsm = []

        #idim0 = 0
        #istation0 = 0
        for icard, card in enumerate(self.cards):
            (pid, mid, beam_type, groupi, xxbi, soi, nsmi, ndimi, dims, comment) = card

            nstationi = len(xxbi)
            all_xxb.extend(xxbi)
            for dim in dims:
                _bar_areaL('PBEAML', beam_type, dim, self)
                all_dims.extend(dim)
            all_so.extend(soi)
            all_nsm.extend(nsmi)
            #station.extend(xxbi)

            #idim1 = idim0 + ndimi * nstationi
            #istation1 = istation0 + nstationi
            nstation[icard] = nstationi
            property_id[icard] = pid
            material_id[icard] = mid
            group[icard] = groupi
            Type[icard] = beam_type
            ndim[icard] = ndimi
            #assert nstationi >= 2, f'pid={pid} mid={mid} beam_type={beam_type!r} xxb={xxbi} dims={dims}'
            #print(f'pid={pid} mid={mid} beam_type={beam_type!r} xxb={xxbi} dims={dims} nsm={nsmi}')
            assert nstationi == len(nsmi), f'pid={pid} nstation={nstationi} nsmi={nsmi}'
            #idim[icard, :] = [idim0, idim1]
            #istation[icard, :] = [istation0, istation1]

            #self.nsm[icard] = nsm
            #idim0 = idim1
            #istation0 = istation1

        xxb = np.array(all_xxb, dtype='float64')
        dims = np.array(all_dims, dtype='float64')
        so = np.array(all_so, dtype='|U4')
        nsm = np.array(all_nsm, dtype='float64')
        #ndim_total = self.ndim.sum()
        self._save(
            property_id, material_id,
            ndim,
            nstation,
            Type, group,
            xxb, dims, so, nsm)
        nstation_total = self.nstation.sum()
        self.sort()

        assert self.xxb.shape[0] == nstation_total
        #assert len(self.dims) == ndim_total
        assert len(self.so) == nstation_total
        assert len(self._nsm) == nstation_total, f'nsm={self._nsm}; nstations_total={nstation_total}'
        self.cards = []

    @property
    def istation(self) -> np.ndarray:
        nprops = len(self.property_id)
        istation = np.zeros((nprops, 2), dtype='int32')
        csum = np.cumsum(self.nstation)
        istation[:, 0] = np.hstack([0, csum[:-1]])
        istation[:, 1] = csum
        return istation

    @property
    def idim(self) -> np.ndarray:
        nprops = len(self.property_id)
        idim = np.zeros((nprops, 2), dtype='int32')
        csum = np.cumsum(self.ndim * self.nstation)
        idim[:, 0] = np.hstack([0, csum[:-1]])
        idim[:, 1] = csum
        return idim

    #@property
    #def idim(self) -> np.ndarray:
        #return make_idim(self.n, self.ndim)

    def _save(self, property_id: np.ndarray, material_id: np.ndarray,
              ndim: np.ndarray,
              nstation: np.ndarray,
              Type: np.ndarray,
              group: np.ndarray,
              xxb: np.ndarray, dims: np.ndarray,
              so: np.ndarray, nsm: np.ndarray):
        self.property_id = property_id
        self.material_id = material_id

        assert ndim.min() >= 1, ndim
        self.ndim = ndim

        #assert istation.ndim == 2, istation.shape
        #assert istation.min() == 0, istation
        #assert nstation.min() >= 2, nstation
        assert nstation.min() >= 1, nstation
        #self.istation = istation
        self.nstation = nstation

        self.Type = Type
        self.group = group

        self._nsm = nsm
        self.xxb = xxb
        self.dims = dims
        self.so = so
        self._nsm = nsm
        #ndims = ndim.sum()
        nstations = self.nstation.sum()
        assert nstations == len(self.xxb), (self.xxb, self.xxb)
        assert nstations == len(self.so), (self.xxb, self.so)
        assert nstations == len(self._nsm), (self.xxb, self._nsm)
        #assert nx == len(self.so), (self.xxb, self.so)


    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    def convert(self, xyz_scale: float=1.0,
                nsm_per_length_scale: float=1.0, **kwargs):
        self.xxb *= xyz_scale
        self.dims *= xyz_scale
        self._nsm *= nsm_per_length_scale

    def __apply_slice__(self, prop: PBEAML, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]

        #self.istation = hslice_by_idim(i, idim, elements)
        prop.Type = self.Type[i]
        prop.group = self.group[i]
        #prop._nsm = self._nsm[i]

        idim = self.idim
        istation = self.istation
        prop.dims = hslice_by_idim(i, idim, self.dims)
        prop.xxb = hslice_by_idim(i, istation, self.xxb)
        prop._nsm = hslice_by_idim(i, istation, self._nsm)

        #self.idim = np.zeros((ncards, 2), dtype='int32')  # for all properties
        #prop.idim = self.idim[i, :]
        prop.ndim = self.ndim[i]

        #self.istation = istation
        #prop.istation = self.istation[i, :]
        prop.nstation = self.nstation[i]

        prop.so = self.so[i]
        prop.n = len(i)

        nproperties = len(prop.property_id)
        nstations = prop.nstation.sum()
        assert len(prop.xxb) == nstations
        assert len(prop.group) == nproperties
        assert len(prop._nsm) == nstations

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          msg=f'no materials for {self.type}; {self.all_materials}')
        mids = np.unique(mids)
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        assert len(self.property_id) == len(self.material_id)
        assert len(self.property_id) == len(self.ndim)
        assert len(self.property_id) == len(self.idim)
        assert len(self.property_id) == len(self.nstation)
        assert len(self.property_id) == len(self.istation)
        assert len(self.property_id) == len(self.group)
        for pid, mid, beam_type, ndim, idim, nstation, istation, group in zip_strict(
                self.property_id, self.material_id, self.Type,
                self.ndim, self.idim,
                self.nstation, self.istation,
                self.group):

            idim0, idim1 = idim
            istation0, istation1 = istation
            xxb = self.xxb[istation0:istation1]
            nsm = self._nsm[istation0:istation1]
            so = self.so[istation0:istation1]
            dims = self.dims[idim0 : idim1].reshape(nstation, ndim)

            group = set_blank_if_default(group, 'MSCBML0')
            ndim = self.valid_types[beam_type]
            assert len(dims[0, :]) == ndim, 'PBEAML ndim=%s len(dims)=%s' % (ndim, len(dims))

            list_fields = ['PBEAML', pid, mid, group, beam_type,
                           None, None, None, None]
            #print("xxb=%s so=%s dim=%s nsm=%s" % (
                #xxb, so, dims, nsm))

            dims_equal = True
            dim0 = dims[0, :]
            for dim in dims[1:, :]:
                if not np.array_equal(dim0, dim):
                    dims_equal = False
                    break

            if dims_equal and len(xxb) == 2 and so[0] == so[1] and len(nsm) == 2 and nsm[0] == nsm[1]:
                list_fields += dims[0].tolist() + [nsm[0]]
            else:
                for (i, xxbi, soi, dimi, nsmi) in zip(count(), xxb, so, dims, nsm):
                    if i == 0:
                        list_fields += dimi.tolist() + [nsmi]
                    else:
                        list_fields += [soi, xxbi] + dimi.tolist() + [nsmi]
            bdf_file.write(print_card(list_fields))
        return

    def area(self) -> np.ndarray:
        nproperties = len(self.property_id)
        area = np.zeros(nproperties, dtype='float64')
        for i, pid, beam_type, ndim, idim, nstation, istation in zip(count(), self.property_id, self.Type,
                                                                     self.ndim, self.idim,
                                                                     self.nstation, self.istation,):
            idim0, idim1 = idim
            istation0, istation1 = istation
            xxb = self.xxb[istation0:istation1]
            #nsm = self.nsm[istation0:istation1]
            #so = self.so[istation0:istation1]
            dims = self.dims[idim0 : idim1].reshape(nstation, ndim)
            dims_str = str(dims).replace('\n', '')
            #self.model.log.info(f'pid={pid} beam_type={beam_type!r} dims={dims_str} idim0={idim0} idim1={idim1}')

            #prop = pbarl(self.property_id[i], self.material_id[i], beam_type, dim)
            #A, I1, I2, I12 = A_I1_I2_I12(prop, beam_type, dim)
            areasi = []
            for dim in dims:
                area_i1_i2_i12 = _bar_areaL('PBEAML', beam_type, dim, self)
                areasi.append(area_i1_i2_i12[0])
            if len(xxb) == 1:
                areaii = areasi[0]
                xxb = [0., 1.]
                areasi = [areaii, areaii]
            A = integrate_positive_unit_line(xxb, areasi)
            area[i] = A
        return area

    def rho(self) -> np.ndarray:
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        return rho

    def nsm(self) -> np.ndarray:
        nproperties = len(self.property_id)
        nsm = np.zeros(nproperties, dtype='float64')
        for i, beam_type, ndim, idim, nstation, istation in zip(count(), self.Type,
                                                                self.ndim, self.idim,
                                                                self.nstation, self.istation,):
            #idim0, idim1 = idim
            istation0, istation1 = istation
            xxb = self.xxb[istation0:istation1]
            nsms = self._nsm[istation0:istation1]
            #so = self.so[istation0:istation1]
            #dims = self.dims[idim0 : idim1].reshape(nstation, ndim)

            #prop = pbarl(self.property_id[i], self.material_id[i], beam_type, dim)
            #A, I1, I2, I12 = A_I1_I2_I12(prop, beam_type, dim)
            #areasi = []
            #for dim in dims:
                #area_i1_i2_i12 = _bar_areaL('PBEAML', beam_type, dim, self)
                #areasi.append(area_i1_i2_i12[0])

            assert len(xxb) == len(nsms)
            if len(xxb) == 1:
                nsmi = nsms[0]
                xxb = [0., 1.]
                nsms = [nsmi, nsmi]
            nsmi = integrate_positive_unit_line(xxb, nsms)
            nsm[i] = nsmi
        return nsm

    def rho(self) -> np.ndarray:
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        return rho

    @property
    def all_materials(self) -> list[MAT1]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[MAT1]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def mass_per_length(self) -> np.ndarray:
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        A = self.area()
        nsm = self.nsm()
        return rho * A + nsm


class PBCOMP(Property):
    """
    +--------+------+-----+-----+------+----+-----+--------+-----+
    |   1    |   2  |  3  |  4  |   5  |  6 |  7  |   8    |  9  |
    +========+======+=====+=====+======+====+=====+========+=====+
    | PBCOMP | PID  | MID | A   |  I1  | I2 | I12 |   J    | NSM |
    +--------+------+-----+-----+------+----+-----+--------+-----+
    |        |  K1  | K2  | M1  |  M2  | N1 | N2  | SYMOPT |     |
    +--------+------+-----+-----+------+----+-----+--------+-----+
    |        |  Y1  | Z1  | C1  | MID1 |    |     |        |     |
    +--------+------+-----+-----+------+----+-----+--------+-----+
    |        |  Y2  | Z2  | C2  | MID2 |    |     |        |     |
    +--------+------+-----+-----+------+----+-----+--------+-----+
    |        | ...  | ... | ... |      |    |     |        |     |
    +--------+------+-----+-----+------+----+-----+--------+-----+
    """
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self._area = np.array([], dtype='float64')
        self.j = np.array([], dtype='float64')
        self.i1 = np.array([], dtype='float64')
        self.i2 = np.array([], dtype='float64')
        self.i12 = np.array([], dtype='float64')
        self.nsm = np.array([], dtype='float64')
        self.k = np.array([], dtype='float64')
        self.m1 = np.array([], dtype='float64')
        self.n1 = np.array([], dtype='float64')
        self.m2 = np.array([], dtype='float64')
        self.n2 = np.array([], dtype='float64')
        self.y = np.array([], dtype='float64')
        self.z = np.array([], dtype='float64')
        self.c = np.array([], dtype='float64')
        self.material_ids = np.array([], dtype='int32')

    def add(self, pid: int, mid: int, y: list[float], z: list[float],
            c: list[float], mids: list[int],
            area: float=0.0, i1: float=0.0, i2: float=0.0, i12: float=0.0,
            j: float=0.0, nsm: float=0.0,
            k1: float=1.0, k2: float=1.0,
            m1: float=0.0, m2: float=0.0,
            n1: float=0.0, n2: float=0.0,
            symopt: int=0, comment: str='') -> int:
        """
        Creates a PBCOMP card

        Parameters
        ----------
        pid : int
            Property ID
        mid : int
            Material ID
        mids : list[int]
            Material ID for the i-th integration point
        y / z : list[float]
            The (y,z) coordinates of the lumped areas in the element
            coordinate system
        c : list[float]; default=0.0
            Fraction of the total area for the i-th lumped area
            default not supported...
        area : float
            Area of beam cross section
        i1 / i2 : float; default=0.0
            Area moment of inertia about plane 1/2 about the neutral axis
        i12 : float; default=0.0
           area product of inertia
        j : float; default=0.0
            Torsional moment of interia
        nsm : float; default=0.0
            Nonstructural mass per unit length
        k1 / k2 : float; default=1.0
            Shear stiffness factor K in K*A*G for plane 1/2
        m1 / m2 : float; default=0.0
            The (y,z) coordinates of center of gravity of nonstructural mass
        n1 / n2 : float; default=0.0
            The (y,z) coordinates of neutral axis
        symopt : int; default=0
            Symmetry option to input lumped areas for the beam cross section
            0 < Integer < 5
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, mid, y, z, c, mids,
                           area, i1, i2, i12, j, nsm,
                           k1, k2, m1, m2, n1, n2,
                           symopt, y, z, c, mids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        area = double_or_blank(card, 3, 'Area', default=0.0)
        i1 = double_or_blank(card, 4, 'I1', default=0.0)
        i2 = double_or_blank(card, 5, 'I2', default=0.0)
        i12 = double_or_blank(card, 6, 'I12', default=0.0)
        j = double_or_blank(card, 7, 'J', default=0.0)
        nsm = double_or_blank(card, 8, 'nsm', default=0.0)
        k1 = double_or_blank(card, 9, 'k1', default=1.0)
        k2 = double_or_blank(card, 10, 'k2', default=1.0)
        m1 = double_or_blank(card, 11, 'm1', default=0.0)
        m2 = double_or_blank(card, 12, 'm2', default=0.0)
        n1 = double_or_blank(card, 13, 'n1', default=0.0)
        n2 = double_or_blank(card, 14, 'n2', default=0.0)
        symopt = integer_or_blank(card, 15, 'symopt', default=0)
        y = []
        z = []
        c = []
        mids = []

        nfields = len(card) - 17
        nrows = nfields // 8
        if nfields % 8 > 0:
            nrows += 1

        for row in range(nrows):
            i = 8 * row + 17
            yi = double_or_blank(card, i, 'y' + str(row), default=0.0)
            zi = double_or_blank(card, i + 1, 'z' + str(row), default=0.0)
            ci = double_or_blank(card, i + 2, 'c' + str(row), default=0.0)
            mid = integer_or_blank(card, i + 3, 'mid' + str(row), mid)
            y.append(yi)
            z.append(zi)
            c.append(ci)
            mids.append(mid)
        if symopt != 0:
            assert len(y) > 0, f'y={y} symopt={symopt} card={card}'
        else:
            assert len(y) == 0, f'y={y} symopt={symopt} card={card}'

        self.cards.append((pid, mid, y, z, c, mids,
                           area, i1, i2, i12, j, nsm,
                           k1, k2, m1, m2, n1, n2,
                           symopt, y, z, c, mids, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        self.property_id = np.zeros(ncards, dtype='int32')
        self.material_id = np.zeros(ncards, dtype='int32')
        self.symopt = np.zeros(ncards, dtype='int32')
        self._area = np.zeros(ncards, dtype='float64')
        self.j = np.zeros(ncards, dtype='float64')
        self.inertia = np.zeros((ncards, 3), dtype='float64')
        #self.i1 = np.zeros(ncards, dtype='float64')
        #self.i2 = np.zeros(ncards, dtype='float64')
        #self.i12 = np.zeros(ncards, dtype='float64')
        self.nsm = np.zeros(ncards, dtype='float64')
        self.k = np.zeros((ncards, 2), dtype='float64')

        self.m1 = np.zeros(ncards, dtype='float64')
        self.m2 = np.zeros(ncards, dtype='float64')
        self.n1 = np.zeros(ncards, dtype='float64')
        self.n2 = np.zeros(ncards, dtype='float64')
        #self.y = np.array([], dtype='float64')
        #self.z = np.array([], dtype='float64')
        #self.c = np.array([], dtype='float64')
        #self.material_ids = np.array([], dtype='int32')

        self.nstation = np.zeros(ncards, dtype='int32')

        self.k1 = np.zeros(ncards, dtype='float64')
        self.k2 = np.zeros(ncards, dtype='float64')

        y_list = []
        z_list = []
        c_list = []
        mids_list = []
        for icard, card in enumerate(self.cards):
            (pid, mid, y, z, c, mids,
             area, i1, i2, i12, j, nsm,
             k1, k2, m1, m2, n1, n2,
             symopt, y, z, c, mids, comment) = card

            nstation = len(mids)
            self.property_id[icard] = pid
            self.material_id[icard] = mid
            self.nstation[icard] = nstation
            self.symopt[icard] = symopt

            #self.s1[icard] = s1
            #self.s2[icard] = s2
            self._area[icard] = area
            self.inertia[icard, :] = [i1, i2, i12]
            self.j[icard] = j
            self.nsm[icard] = nsm
            self.k[icard, :] = [k1, k2]
            self.m1[icard] = m1
            self.n1[icard] = n1
            self.m2[icard] = m2
            self.n2[icard] = n2
            #assert len(y) > 0

            y_list.extend(y)
            z_list.extend(z)
            c_list.extend(c)
            mids_list.extend(mids)

        self.y = np.array(y_list, dtype='float64')
        self.z = np.array(z_list, dtype='float64')
        self.c = np.array(c_list, dtype='float64')
        #assert len(y_list) > 0, y_list
        self.material_ids = np.array(mids_list, dtype='int32')
        self.sort()

    def __apply_slice__(self, prop: PBCOMP, i: np.ndarray) -> None:
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]

        prop.k = self.k[i, :]
        prop.m1 = self.m1[i]
        prop.m2 = self.m2[i]
        prop.n1 = self.n1[i]
        prop.n2 = self.n2[i]
        prop.nsm = self.nsm[i]

        istation = self.istation
        prop.y = hslice_by_idim(i, istation, self.y)
        prop.z = hslice_by_idim(i, istation, self.z)
        prop.c = hslice_by_idim(i, istation, self.c)
        prop.material_ids = hslice_by_idim(i, istation, self.material_ids)
        prop.nstation = self.nstation[i]
        prop.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        materials = self.allowed_materials
        mids = hstack_msg([mat.material_id for mat in materials],
                          msg=f'no materials for {self.type}; {self.all_materials}')
        #mids2 = hstack_msg([prop.material_ids for mat in materials],
                           #msg=f'no materials for {self.type}; {self.all_materials}')
        #material_ids = np.hstack([mids, mids2])
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_ids.ravel()))

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)
        if len(self.material_ids):
            used_dict['material_id'].append(self.material_ids.ravel())

    @property
    def all_materials(self) -> list[MAT1]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[MAT1]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def mass_per_length(self) -> np.ndarray:
        #return self.nsm
        nproperties = len(self.property_id)
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        if rho.max() == 0. and rho.min() == 0. and self.nsm.max() == 0. and self.nsm.min() == 0.:
            return np.zeros(nproperties, dtype=rho.dtype)

        mass_per_length = rho * self._area + self.nsm
        self.model.log.debug(f'pids={self.property_id} mass/length={mass_per_length}')
        #for i, rhoi, istation in zip(count(), rho, self.istation):
            #istation0, istation1 = istation
            #assert istation1 > istation0
            #xxb = self.xxb[istation0:istation1]
            #area = self.A[istation0:istation1]
            #nsm = self.nsm[istation0:istation1]
            #mass_per_lengths = rhoi * area + nsm
            #mass_per_lengthi = integrate_positive_unit_line(xxb, mass_per_lengths)
            #mass_per_length[i] = mass_per_lengthi

        assert len(mass_per_length) == nproperties
        return mass_per_length


    @property
    def istation(self) -> np.ndarray:
        idim = make_idim(self.n, self.nstation)
        return idim

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_ids = array_str(self.property_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        for pid, mid, A, i1, i2, i12, j, nsm, \
            (k1, k2), m1, m2, n1, n2, symopt, (istation0, istation1) in zip_longest(
            property_ids, material_ids, self._area, self.i1, self.i2, self.i12, self.j, self.nsm,
            self.k, self.m1, self.m2, self.n1, self.n2, self.symopt, self.istation):

            area = set_blank_if_default(A, default=0.0)
            j = set_blank_if_default(j, default=0.0)
            i1 = set_blank_if_default(i1, default=0.0)
            i2 = set_blank_if_default(i2, default=0.0)
            i12 = set_blank_if_default(i12, default=0.0)
            nsm = set_blank_if_default(nsm, default=0.0)

            k1 = set_blank_if_default(k1, default=1.0)
            k2 = set_blank_if_default(k2, default=1.0)

            m1 = set_blank_if_default(m1, default=0.0)
            m2 = set_blank_if_default(m2, default=0.0)

            n1 = set_blank_if_default(n1, default=0.0)
            n2 = set_blank_if_default(n2, default=0.0)

            symopt = set_blank_if_default(symopt, default=0)

            list_fields = ['PBCOMP', pid, mid, area, i1, i2, i12, j,
                           nsm, k1, k2, m1, m2, n1, n2, symopt, None]

            y = self.y[istation0:istation1]
            z = self.z[istation0:istation1]
            c = self.c[istation0:istation1]
            mids = self.material_ids[istation0:istation1]
            for (yi, zi, ci, mid) in zip_longest(y, z, c, mids):
                ci = set_blank_if_default(ci, default=0.0)
                list_fields += [yi, zi, ci, mid, None, None, None, None]
            #assert len(y) > 0, list_fields
            bdf_file.write(print_card(list_fields))
        return

    @property
    def istation(self) -> np.ndarray:
        return make_idim(self.n, self.nstation)

    def area(self) -> np.ndarray:
        return self._area

def _sort_pbeam(pid, xxb, so, area, i1, i2, i12, j, nsm,
                c1, c2, d1, d2, e1, e2, f1, f2, ensure_xxb_1_section=True):
    nxxb = len(xxb)
    # sort xxb
    ixxb = np.argsort(xxb)
    duplicate_xxb = (len(xxb) == 1)
    if duplicate_xxb and ensure_xxb_1_section:
        ixxb = np.array([0, 0])
    self_so = np.array(so, dtype='|U8')[ixxb]
    self_xxb = np.array(xxb, dtype='float64')[ixxb]
    #print('ixxb = %s' % ixxb)
    #print('i12 = %s' % i12)

    assert len(area) == nxxb, 'pid=%s len(xxb)=%s len(A =)=%s' % (pid, nxxb, len(area))
    assert len(i1) == nxxb, 'pid=%s len(xxb)=%s len(i1 )=%s' % (pid, nxxb, len(i1))
    assert len(i2) == nxxb, 'pid=%s len(xxb)=%s len(i2 )=%s' % (pid, nxxb, len(i2))
    assert len(i12) == nxxb, 'pid=%s len(xxb)=%s len(i12)=%s' % (pid, nxxb, len(i12))
    assert len(j) == nxxb, 'pid=%s len(xxb)=%s len(j =)=%s' % (pid, nxxb, len(j))
    assert len(nsm) == nxxb, 'pid=%s len(xxb)=%s len(nsm)=%s' % (pid, nxxb, len(nsm))

    self_A = np.array(area, dtype='float64')[ixxb]
    self_i1 = np.array(i1, dtype='float64')[ixxb]
    self_i2 = np.array(i2, dtype='float64')[ixxb]
    self_i12 = np.array(i12, dtype='float64')[ixxb]
    self_j = np.array(j, dtype='float64')[ixxb]
    self_nsm = np.array(nsm, dtype='float64')[ixxb]

    #assert len(area) == nxxb, 'pid=%s len(xxb)=%s len(area)=%s' % (nxxb, len(area))
    assert len(c1) == nxxb, 'pid=%s len(xxb)=%s len(c1)=%s' % (pid, nxxb, len(c1))
    assert len(c2) == nxxb, 'pid=%s len(xxb)=%s len(c2)=%s' % (pid, nxxb, len(c2))
    assert len(d1) == nxxb, 'pid=%s len(xxb)=%s len(d1)=%s' % (pid, nxxb, len(d1))
    assert len(d2) == nxxb, 'pid=%s len(xxb)=%s len(d2)=%s' % (pid, nxxb, len(d2))
    assert len(e1) == nxxb, 'pid=%s len(xxb)=%s len(e1)=%s' % (pid, nxxb, len(e1))
    assert len(e2) == nxxb, 'pid=%s len(xxb)=%s len(e2)=%s' % (pid, nxxb, len(e2))
    assert len(f1) == nxxb, 'pid=%s len(xxb)=%s len(f1)=%s' % (pid, nxxb, len(f1))
    assert len(f2) == nxxb, 'pid=%s len(xxb)=%s len(f2)=%s' % (pid, nxxb, len(f2))

    self_c1 = np.array(c1, dtype='float64')[ixxb]
    self_c2 = np.array(c2, dtype='float64')[ixxb]
    self_d1 = np.array(d1, dtype='float64')[ixxb]
    self_d2 = np.array(d2, dtype='float64')[ixxb]
    self_e1 = np.array(e1, dtype='float64')[ixxb]
    self_e2 = np.array(e2, dtype='float64')[ixxb]
    self_f1 = np.array(f1, dtype='float64')[ixxb]
    self_f2 = np.array(f2, dtype='float64')[ixxb]

    if duplicate_xxb:
        self_xxb = np.array([0., 1.], dtype='float64')

    out = (
        self_xxb, self_so, self_A, self_i1, self_i2, self_i12, self_j, self_nsm,
        self_c1, self_c2, self_d1, self_d2, self_e1, self_e2, self_f1, self_f2,
    )
    return out

def _interpolate_pbeam_sections(pid, xxb, so, area, i1, i2, i12, j, nsm,
                                c1, c2, d1, d2, e1, e2, f1, f2):
    """
    now we interpolate to fix up missing data
    from arrays that were potentially out of order
    (they're sorted now)

    we've also already checked xxb=0.0 and xxb=1.0 for I1, I2, I12, J

    If any fields 4 through 9 are blank on the continuation with the value of
    X/XB =1.0, then the values for A, I1, I2, I12, J and NSM are set to the values
    given for end A.

    For the continuations that have intermediate values of X/XB between 0.0
    and 1.0 and use the default option (any of the fields 4 through 9 are blank), a
    linear interpolation between the values at ends A and B is performed to obtain
    the missing section properties.
    """
    if len(xxb) == 2:
        out = (xxb, so, area, i1, i2, i12, j, nsm,
               c1, c2, d1, d2, e1, e2, f1, f2)
        return out
    assert len(xxb) > 2, xxb

    def _linear_interpolate2(xxb, area):
        isnan = np.isnan(area[1:-1])
        inan = np.arange(1, len(area) - 1)[isnan]
        m = (area[-1] - area[0]) / (xxb[-1] - xxb[0])
        area[inan] = area[0] + m * xxb[inan]
        return area
    area = _linear_interpolate2(xxb, area)
    i1 = _linear_interpolate2(xxb, i1)
    i2 = _linear_interpolate2(xxb, i2)
    i12 = _linear_interpolate2(xxb, i12)
    j = _linear_interpolate2(xxb, j)
    nsm = _linear_interpolate2(xxb, nsm)
    c1 = _linear_interpolate2(xxb, c1)
    c2 = _linear_interpolate2(xxb, c2)
    d1 = _linear_interpolate2(xxb, d1)
    d2 = _linear_interpolate2(xxb, d2)
    e1 = _linear_interpolate2(xxb, e1)
    e2 = _linear_interpolate2(xxb, e2)
    f1 = _linear_interpolate2(xxb, f1)
    f2 = _linear_interpolate2(xxb, f2)

    for i, xxbi in enumerate(xxb):
        #i, xxb, a, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2, f1, f2 = interp_data
        #if xxbi not in [0., 1.]:
        assert area[i] >= 0., area
        assert i1[i] >= 0., i1
        assert i2[i] >= 0., i2
        assert j[i] >= 0., j  # we check warping later
        if i1[i] * i2[i] - i12[i] ** 2 <= 0.:
            msg = 'I1 * I2 - I12^2=0 and must be greater than 0.0 at End B\n'
            msg += f'pid={pid} xxb={xxbi} i1={i1[i]} i2={i2[i]} i12={i12[i]}'
            raise ValueError(msg)

    out = (xxb, so, area, i1, i2, i12, j, nsm,
           c1, c2, d1, d2, e1, e2, f1, f2)
    return out

def _linearly_interpolate(i: int, x: np.ndarray, y: np.ndarray):
    """
    For the continuations that have intermediate values of X/XB between 0.0
    and 1.0 and use the default option (any of the fields 4 through 9 are blank), a
    linear interpolation between the values at ends A and B is performed to obtain
    the missing section properties.
    """
    inonzero = np.where(y != 0.)[0]
    if len(inonzero) == 0:
        return 0.
    idiff = (inonzero - i)
    ineg = (idiff < 0)
    ipos = (idiff > 0)
    ilow = i + idiff[ineg]
    ihigh = i + idiff[ipos]

    xlow = x[ilow]
    ylow = y[ilow]
    xhigh = x[ihigh]
    yhigh = y[ihigh]
    yi = y[ilow] + (yhigh - ylow) / (xhigh - xlow) * x[i]
    assert isinstance(yi[0], float), yi
    return yi


class CBEND(Element):
    """
    NX 2020.1

    Defines a curved beam, curved pipe, or elbow element.

    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |    9     |
    +=======+=====+=====+=====+=====+=====+=====+=====+==========+
    | CBEND | EID | PID | GA  | GB  | X1  | X2  | X3  |   GEOM   |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+

    """

    @Element.clear_check
    def clear(self) -> None:
        self.element_id: np.array = np.array([], dtype='int32')
        self.property_id: np.array = np.array([], dtype='int32')
        self.nodes: np.array = np.zeros((0, 2), dtype='int32')
        self.g0: np.array = np.array([], dtype='int32')
        self.x: np.array = np.zeros((0, 3), dtype='float64')
        self.geom_flag: np.array = np.array([], dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int],
            g0: Optional[int], x: Optional[list[float]],
            geom: str='GGG', comment: str='') -> int:
        assert g0 is None or x is None, (g0, x)
        if g0 is None and not isinstance(x, (list, tuple, np.ndarray)):
            raise TypeError(f'PBEND eid={eid} x={x} and should be a list[float]')
        if x is None and not isinstance(g0, integer_types):
            raise TypeError(f'PBEND eid={eid} g0={g0} and should be an integer')

        self.cards.append((eid, pid, nids, g0, x, geom, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CBEND card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        x1_g0 = integer_double_or_blank(card, 5, 'x1_g0', 0.0)
        if isinstance(x1_g0, integer_types):
            g0 = x1_g0
            x = None
        elif isinstance(x1_g0, float):
            g0 = None
            x = np.array([double_or_blank(card, 5, 'x1', 0.0),
                          double_or_blank(card, 6, 'x2', 0.0),
                          double_or_blank(card, 7, 'x3', 0.0)], dtype='float64')
            if np.linalg.norm(x) == 0.0:
                msg = 'G0 vector defining plane 1 is not defined.\n'
                msg += 'G0 = %s\n' % g0
                msg += 'X  = %s\n' % x
                raise RuntimeError(msg)
        else:
            raise ValueError('invalid x1/g0=%r on CBEND' % x1_g0)
        geom = integer(card, 8, 'geom')

        assert len(card) == 9, f'len(CBEND card) = {len(card):d}\ncard={card}'
        #return CBEND(eid, pid, [ga, gb], g0, x, geom, comment=comment)
        self.cards.append((eid, pid, [ga, gb], g0, x, geom, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        g0 = np.zeros(ncards, dtype=idtype)
        x = np.full((ncards, 3), np.nan, dtype='float64')
        geom_flag = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, g0i, xi, geom_flagi, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            if g0i in {None, 0, -1}:
                x[icard, :] = xi
            else:
                assert g0i > 0, card
                g0[icard] = g0i

            geom_flag[icard] = geom_flagi

        self._save(element_id, property_id, nodes,
                   g0, x, geom_flag)
        self.cards = []

    def _save(self, element_id, property_id, nodes,
              g0, x, geom_flag) -> None:
        if len(self.element_id) != 0:
            asdf
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.g0 = g0
        self.x = x

        self.geom_flag = geom_flag
        self.n = len(property_id)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['element_id'].append(self.element_id)
        used_dict['property_id'].append(self.property_id)
        used_dict['node_id'].append(self.nodes.ravel())
        g0 = self.g0[self.is_g0]
        if len(g0):
            used_dict['node_id'].append(g0)

    def convert(self, xyz_scale: float=1.0,
                mass_scale: float=1.0, **kwargs):
        ## TODO: probably wrong for CD=1
        self.x *= xyz_scale

    def __apply_slice__(self, elem: CBEND, i: np.ndarray) -> None:
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.g0 = self.g0[i]
        elem.x = self.x[i, :]
        elem.geom_flag = self.geom_flag[i]
        elem.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no bend properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes_ = array_str(self.nodes, size=size)
        geom_flag = array_str(self.geom_flag, size=size)
        for eid, pid, nodes, g0, x, is_g0, geom_flagi in zip_longest(
            element_ids, property_ids, nodes_,
            self.g0, self.x, self.is_g0, geom_flag):

            n1, n2 = nodes
            if is_g0:
                x1 = g0
                x2 = ''
                x3 = ''
            else:
                x1, x2, x3 = x # self.get_x_g0_defaults()

            list_fields = ['CBEND', eid, pid, n1, n2,
                           x1, x2, x3, geom_flagi]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def is_x(self) -> np.ndarray:
        return (self.g0 == 0)

    @property
    def is_g0(self) -> np.ndarray:
        return ~self.is_x

    @property
    def all_properties(self) -> list[PBEND]:
        model = self.model
        return [model.pbend]

    @property
    def allowed_properties(self):
        all_properties = self.all_properties
        props = [prop for prop in all_properties if prop.n > 0]
        assert len(props) > 0, f'{self.type}: all_props={all_properties}'
        return props

    def mass(self) -> np.ndarray:
        #pid = self.property_id
        mass_per_length = line_pid_mass_per_length(self.property_id, self.allowed_properties)
        length = self.length()
        mass = mass_per_length * length
        return mass

    def line_vector_length(self) -> tuple[np.ndarray, np.ndarray]:
        line_vector, length = line_vector_length(self.model, self.nodes)
        return line_vector, length

    def length(self) -> np.ndarray:
        nelement = self.n
        length = np.full(nelement, np.nan, dtype='float64')

        igeom_flag1 = (self.geom_flag == 1)
        igeom_flag2 = (self.geom_flag == 2)
        igeom_flag3 = (self.geom_flag == 3)
        igeom_flag4 = (self.geom_flag == 4)
        igeom_flag_no = ~(igeom_flag1 | igeom_flag2 | igeom_flag3 | igeom_flag4)
        if igeom_flag_no.sum():
            i = np.where(igeom_flag_no)[0]
            cbends = self.slice_card_by_index(i).write()
            raise RuntimeError(f'Invalid geom flags (should be [1, 2, 3, 4])\n{cbends}')

        grid = self.model.grid
        xyz_cid0 = grid.xyz_cid0()

        i = np.full((nelement, 3), np.nan, dtype='float64')
        j = np.full((nelement, 3), np.nan, dtype='float64')
        k = np.full((nelement, 3), np.nan, dtype='float64')
        inode = grid.index(self.nodes)

        inode1 = inode[:, 0]
        inode2 = inode[:, 1]

        xyz1 = xyz_cid0[inode1, :]
        xyz2 = xyz_cid0[inode2, :]

        v, cd = get_bar_vector(self, xyz1)
        cd1 = cd[:, 0]
        cd1_ref = self.model.coord.slice_card_by_id(cd1, sort_ids=False)
        xyz0 = cd1_ref.transform_xyz_to_global_assuming_rectangular(v)

        print('xyz0 =', xyz0)
        length12 = np.linalg.norm(xyz2 - xyz1, axis=1)
        if igeom_flag1.sum():
            xyz1i = xyz1[igeom_flag1, :]
            xyz2i = xyz2[igeom_flag1, :]
            xyz0i = xyz0[igeom_flag1, :]

            dxyz10 = xyz1i - xyz0i
            radius = np.linalg.norm(dxyz10, axis=1)
            theta = np.arctan2(xyz2i - xyz0i, xyz1i - xyz0i)
            theta_deg = np.degrees(theta)
            #iprime1 = xyz1_1 - xyz0_1
            #j1 = xyz2_1 - xyz0_1
            #j1_norm = np.linalg.norm(j1, axis=1)
            #j1 /= j1_norm[:, np.newaxis]
            #k1 = np.cross(iprime1, j1, axis=1)
            #k1_norm = np.linalg.norm(k1, axis=1)
            #k1 /= k1_norm[:, np.newaxis]
            #i1 = np.cross(k1, j1, axis=1)

            #i[igeom_flag1, :] = i1
            #j[igeom_flag1, :] = j1
            #k[igeom_flag1, :] = k1
            #coord0 = self.model.coord.slice_card_by_id(0)
            ##for xyz0_1i, xyz1_1i, xyz2_1i in zip(xyz0i, xyz1i, xyz2i):
                ##coord0.add_cord2c(1, xyz2_1i, xyz0_1i, xyz1_1i)
            #coord0.parse_cards()
            #print(i1, j1, k1)
            chord_sq = radius**2 - length12[igeom_flag1]**2
            chord = np.full(self.element_id.shape, np.nan, dtype=xyz0.dtype)
            ipos = chord_sq > 0
            chord = 2 * np.sqrt(chord_sq[ipos])
            length[igeom_flag1] = chord[ipos]


        if igeom_flag2.sum():
            i = np.where(igeom_flag2)[0]
            cbends = self.slice_card_by_index(i).write()
            raise RuntimeError(f'Invalid geom flags (should be [1, 2, 3, 4])\n{cbends}')
        #if igeom_flag3.sum():
            #i = np.where(igeom_flag2)[0]
            #cbends = self.slice_card_by_index(i).write()
            #raise RuntimeError(f'Invalid geom flags (should be [1, 2, 3, 4])\n{cbends}')
        if igeom_flag4.sum():
            i = np.where(igeom_flag3)[0]
            cbends = self.slice_card_by_index(i).write()
            raise RuntimeError(f'Invalid geom flags (should be [1, 2, 3, 4])\n{cbends}')

        #length = line_length(self.model, self.nodes)
        #inan = np.isnan(length)
        #if np.any(inan):
            #msg = 'CBEAM has nan length\n'
            #msg += f'eids={self.element_id[inan]}\n'
            #msg += f'nid1={self.nodes[inan,0]}\n'
            #msg += f'nid2={self.nodes[inan,1]}\n'
            #raise RuntimeError(msg)
        print(f'length = {length}')
        return length

    def centroid(self) -> np.ndarray:
        centroid = line_centroid(self.model, self.nodes)
        return centroid

    #def get_bar_vector(self, xyz1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        #v, cd = get_bar_vector(self, xyz1)
        #return v, cd

    def get_xyz(self) -> tuple[np.ndarray, np.ndarray]:
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

    def get_axes(self, xyz1: np.ndarray, xyz2: np.ndarray,
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                            np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError('get_axes')

    def center_of_mass(self) -> np.ndarray:
        """
        Geom flag
         - 1: Curve from 1->2 (A->B) with center at O
         -    The center of curvature lies on the line AO (or its extension) or vector v.

              O
             / \
            /   \
           /     \
          1       2

        i' = n1 - n0
        j = n2 - n0
        k = i' x j
        i = k x j

         - 2: The tangent of centroid arc at end A is parallel to line AO or vector v.
              Point O (or vector v) and the arc AB must be on the same side of the chord AB.
            - Let point O to be point T
            - Let point O be the center as before

          https://www.varsitytutors.com/hotmath/hotmath_help/topics/tangent-to-a-circle
          AT = given
          OA = given
          (OA)^2 + (AT)^2 = (OT)^2
          OT = sqrt(AT^2 - OA^2)
          i' = n1 - n0
          j' = n2 - n0
          k = i' x j'
          j = k x i' ??? should be close, might have flipped sign...make a simple test that includes plotting

        """
        #self.check_missing(self.property_id)
        #log = self.model.log

        xyz1, xyz2 = self.get_xyz()
        #neids = xyz1.shape[0]
        centroid = (xyz1 + xyz2) / 2.

        ## TODO: doesn't consider curvature...
        return centroid

    def area(self) -> np.ndarray:
        pid = self.property_id
        area = np.full(len(pid), np.nan, dtype='float64')
        log = self.model.log
        for prop in self.allowed_properties:
            i_lookup, i_all = searchsorted_filter(prop.property_id, pid, msg='')
            if len(i_lookup) == 0:
                continue
            # we're at least using some properties
            areai = prop.area()
            area_all = areai[i_all]
            inan = np.isnan(area_all)
            if np.any(inan):
                msg = f'{prop.type} has nan area for property_ids={prop.property_id[inan]}\n'
                log.warning(msg)

            area[i_lookup] = areai[i_all]

        inan = np.isnan(area)
        if np.any(inan):
            msg = 'CBEND has nan area\n'
            msg += f'eids={self.element_id[inan]}\n'
            msg += f'pid={self.property_id[inan]}\n'
            msg += f'all_properties={self.all_properties}'
            #msg += f'As={self.nodes[inan]}\n'
            raise RuntimeError(msg)
        return area

    def volume(self) -> np.ndarray:
        A = self.area()
        L = self.length()
        return A * L


class PBEND(Property):
    """
    MSC/NX Option A

    +-------+------+-------+-----+----+----+--------+----+--------+
    |   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
    +=======+======+=======+=====+====+====+========+====+========+
    | PBEND | PID  |  MID  | A   | I1 | I2 |   J    | RB | THETAB |
    +-------+------+-------+-----+----+----+--------+----+--------+
    |       |  C1  |  C2   | D1  | D2 | E1 |   E2   | F1 |   F2   |
    +-------+------+-------+-----+----+----+--------+----+--------+
    |       |  K1  |  K2   | NSM | RC | ZC | DELTAN |    |        |
    +-------+------+-------+-----+----+----+--------+----+--------+

    MSC Option B

    +-------+------+-------+-----+----+----+--------+----+--------+
    |   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
    +=======+======+=======+=====+====+====+========+====+========+
    | PBEND | PID  |  MID  | FSI | RM | T  |   P    | RB | THETAB |
    +-------+------+-------+-----+----+----+--------+----+--------+
    |       |      |       | NSM | RC | ZC |        |    |        |
    +-------+------+-------+-----+----+----+--------+----+--------+

    NX Option B

    +-------+------+-------+-----+----+----+--------+----+--------+
    |   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
    +=======+======+=======+=====+====+====+========+====+========+
    | PBEND | PID  |  MID  | FSI | RM | T  |   P    | RB | THETAB |
    +-------+------+-------+-----+----+----+--------+----+--------+
    |       | SACL | ALPHA | NSM | RC | ZC | FLANGE |    |        |
    +-------+------+-------+-----+----+----+--------+----+--------+
    |       |  KX  |  KY   | KZ  |    | SY |   SZ   |    |        |
    +-------+------+-------+-----+----+----+--------+----+--------+
    """
    @Property.clear_check
    def clear(self) -> None:
        self.property_id: np.array = np.array([], dtype='int32')
        self.material_id: np.array = np.array([], dtype='int32')
        self.beam_type: np.array = np.array([], dtype='int32')

        self.nsm = np.array([], dtype='float64')

        #: Bend radius of the line of centroids.
        self.rb = np.array([], dtype='float64')

        #------------------------------------------------------
        ## beam_type = 1

        #: Arc angle of element
        self.theta_b = np.array([], dtype='float64')

        #: Shear stiffness factor K in K*A*G for plane 1 and plane 2.
        self.k1 = np.array([], dtype='float64')
        self.k2 = np.array([], dtype='float64')

        self.A = np.array([], dtype='float64')
        self.J = np.array([], dtype='float64')
        self.I1 = np.array([], dtype='float64')
        self.I2 = np.array([], dtype='float64')

        #: The r,z locations from the geometric centroid for stress data recovery.
        self.c1 = np.array([], dtype='float64')
        self.c2 = np.array([], dtype='float64')
        self.d1 = np.array([], dtype='float64')
        self.d2 = np.array([], dtype='float64')
        self.e1 = np.array([], dtype='float64')
        self.e2 = np.array([], dtype='float64')
        self.f1 = np.array([], dtype='float64')
        self.f2 = np.array([], dtype='float64')

        #------------------------------------------------------
        ## beam_type = 2

        # Flag selecting the flexibility and stress intensification factors.
        # See "Flexibility and stress intensification factors" in the Simcenter
        # Nastran Element Library. (Integer = 1-6)
        self.fsi = np.array([], dtype='int32')

        #: Internal pressure.
        self.p = np.array([], dtype='float64')

        #: Wall thickness of the curved pipe.
        self.t = np.array([], dtype='float64')

        #: Mean cross-sectional radius of the curved pipe.
        self.rm = np.array([], dtype='float64')

        self.coincedent_spacing = np.array([], dtype='float64')

        #: For FSI=6, the user defined flexibility factor for the:
        #  - X: torsional moment
        #  - Y: out-of-plane bending moment
        #  - Z: in-plane bending moment
        self.kx = np.array([], dtype='float64')
        self.ky = np.array([], dtype='float64')
        self.kz = np.array([], dtype='float64')

        #: For FSI=6, the user defined stress intensificatation factor for the:
        #  - Y: out-of-plane bending moment
        #  - Z: in-plane bending moment
        self.sy = np.array([], dtype='float64')
        self.sz = np.array([], dtype='float64')

    #def add(self, pid, mid, xxb, so, area, i1, i2, i12, j, nsm=None,
            #c1=None, c2=None, d1=None, d2=None,
            #e1=None, e2=None, f1=None, f2=None,
            #k1=1., k2=1., s1=0., s2=0.,
            #nsia=0., nsib=None, cwa=0., cwb=None,
            #m1a=0., m2a=0., m1b=None, m2b=None,
            #n1a=0., n2a=0., n1b=None, n2b=None,
            #comment='') -> int:
        #self.cards.append((pid, mid,
                           #xxb, so, area, j, i1, i2, i12, nsm,
                           #c1, c2, d1, d2, e1, e2, f1, f2,
                           #s1, s2, k1, k2,
                           #nsia, nsib, cwa, cwb,
                           #m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b,
                           #comment))
        #self.n += 1
        #return self.n - 1

    def add_beam_type_1(self, pid, mid,
                        A, i1, i2, j,
                        rb=None, theta_b=None,
                        c1=0., c2=0., d1=0., d2=0., e1=0., e2=0., f1=0., f2=0.,
                        k1=None, k2=None,
                        nsm=0., rc=0., zc=0., delta_n=0., comment=''):
        """
        +-------+------+-------+-----+----+----+--------+----+--------+
        |   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
        +=======+======+=======+=====+====+====+========+====+========+
        | PBEND | PID  |  MID  | A   | I1 | I2 |   J    | RB | THETAB |
        +-------+------+-------+-----+----+----+--------+----+--------+
        |       |  C1  |  C2   | D1  | D2 | E1 |   E2   | F1 |   F2   |
        +-------+------+-------+-----+----+----+--------+----+--------+
        |       |  K1  |  K2   | NSM | RC | ZC | DELTAN |    |        |
        +-------+------+-------+-----+----+----+--------+----+--------+

        Parameters
        ----------
        A : float
            cross-sectional area
        i1, i2 : float
            area moments of inertia for plane 1/2
        j : float
            torsional stiffness
        rb : float; default=None
            bend radius of the line of centroids
        theta_b : float; default=None
            arc angle of element (degrees)
        c1, c2, d1, d2, e1, e2, f1, f2 : float; default=0.0
            the r/z locations from the geometric centroid for stress recovery
        k1, k2 : float; default=None
            Shear stiffness factor K in K*A*G for plane 1 and plane 2
        nsm : float; default=0.
            nonstructural mass per unit length???
        zc : float; default=None
            Offset of the geometric centroid in a direction perpendicular to
            the plane of points GA and GB and vector v.
        delta_n : float; default=None
            Radial offset of the neutral axis from the geometric centroid,
            positive is toward the center of curvature
        """
        beam_type = 1
        fsi = None
        rm = None
        t = None
        p = None
        centerline_spacing = None
        alpha = None
        flange = None
        kx = None
        ky = None
        kz = None
        sy = None
        sz = None
        #return PBEND(pid, mid, beam_type, A, i1, i2, j,
                     #c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                     #nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b, comment=comment)

        self.cards.append((pid, mid, beam_type,
                           # type=1
                           A, i1, i2, j,
                           c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                           nsm,
                           # common
                           rc, zc,
                           # beam_type=1
                           delta_n,
                           # beam_type=2
                           fsi, rm, t, p,
                           rb, theta_b,
                           centerline_spacing, alpha, flange, kx, ky, kz, sy, sz,
                           comment))
        self.n += 1
        return self.n - 1

    def add_beam_type_2(self, pid, mid,
                        fsi, rm, t, p=None, rb=None, theta_b=None,
                        nsm=0., rc=0., zc=0., comment=''):
        """
        +-------+------+-------+-----+----+----+--------+----+--------+
        |   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
        +=======+======+=======+=====+====+====+========+====+========+
        | PBEND | PID  |  MID  | FSI | RM | T  |   P    | RB | THETAB |
        +-------+------+-------+-----+----+----+--------+----+--------+
        |       |      |       | NSM | RC | ZC |        |    |        |
        +-------+------+-------+-----+----+----+--------+----+--------+

        Parameters
        ----------
        fsi : int
            Flag selecting the flexibility and stress intensification
            factors. See Remark 3. (Integer = 1, 2, or 3)
        rm : float
            Mean cross-sectional radius of the curved pipe
        t : float
            Wall thickness of the curved pipe
        p : float; default=None
            Internal pressure
        rb : float; default=None
            bend radius of the line of centroids
        theta_b : float; default=None
            arc angle of element (degrees)
        nsm : float; default=0.
            nonstructural mass per unit length???
        rc : float; default=None
            Radial offset of the geometric centroid from points GA and GB.
        zc : float; default=None
            Offset of the geometric centroid in a direction perpendicular
            to the plane of points GA and GB and vector v
        """
        beam_type = 2
        A = None
        i1 = None
        i2 = None
        j = None
        c1 = None
        c2 = None
        d1 = None
        d2 = None
        e1 = None
        e2 = None
        f1 = None
        f2 = None
        k1 = None
        k2 = None
        delta_n = None
        centerline_spacing = 0.1
        alpha = 1.0
        flange = 5

        kx = 1.0
        ky = 1.0
        kz = 1.0
        sy = 1.0
        sz = 1.0

        #return PBEND(pid, mid, beam_type, A, i1, i2, j,
                     #c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                     #nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b, comment=comment)
        self.cards.append((pid, mid, beam_type,
                           # type=1
                           A, i1, i2, j,
                           c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                           nsm,
                           # common
                           rc, zc,
                           # beam_type=1
                           delta_n,
                           # beam_type=2
                           fsi, rm, t, p,
                           rb, theta_b,
                           centerline_spacing, alpha, flange, kx, ky, kz, sy, sz,
                           comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PBEND card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')

        value3 = integer_or_double(card, 3, 'Area/FSI')
        #print("PBEND: area/fsi=%s" % value3)

        # MSC/NX option A
        A = None
        i1 = None
        i2 = None
        j = None
        c1 = None
        c2 = None
        d1 = None
        d2 = None
        e1 = None
        e2 = None
        f1 = None
        f2 = None
        k1 = None
        k2 = None
        delta_n = None

        # MSC option B
        rm = None
        t = None
        p = None

        # NX option B
        centerline_spacing = None
        alpha = None
        flange = None
        kx = None
        ky = None
        kz = None
        sy = None
        sz = None
        if isinstance(value3, float):
            fsi = 0
            beam_type = 1
            #: Area of the beam cross section
            A = double(card, 3, 'A')

            #: Area moments of inertia in planes 1 and 2.
            i1 = double(card, 4, 'I1')
            i2 = double(card, 5, 'I2')

            #: Torsional stiffness :math:`J`
            j = double(card, 6, 'J')

            # line2
            #: The r,z locations from the geometric centroid for stress
            #: data recovery.
            c1 = double_or_blank(card, 9, 'c1', default=0.)
            c2 = double_or_blank(card, 10, 'c2', default=0.)
            d1 = double_or_blank(card, 11, 'd1', default=0.)
            d2 = double_or_blank(card, 12, 'd2', default=0.)
            e1 = double_or_blank(card, 13, 'e1', default=0.)
            e2 = double_or_blank(card, 14, 'e2', default=0.)
            f1 = double_or_blank(card, 15, 'f1', default=0.)
            f2 = double_or_blank(card, 16, 'f2', default=0.)

            # line 3
            #: Shear stiffness factor K in K*A*G for plane 1.
            k1 = double_or_blank(card, 17, 'k1')
            #: Shear stiffness factor K in K*A*G for plane 2.
            k2 = double_or_blank(card, 18, 'k2')

            #: Nonstructural mass per unit length.
            nsm = double_or_blank(card, 19, 'nsm', default=0.)

            #: Radial offset of the geometric centroid from points GA and GB.
            rc = double_or_blank(card, 20, 'rc', default=0.)

            #: Offset of the geometric centroid in a direction perpendicular
            #: to the plane of points GA and GB and vector v
            zc = double_or_blank(card, 21, 'zc', default=0.)

            #: Radial offset of the neutral axis from the geometric
            #: centroid, positive is toward the center of curvature
            delta_n = double_or_blank(card, 22, 'delta_n', default=0.)

        elif isinstance(value3, int):  # alternate form
            beam_type = 2
            #: Flag selecting the flexibility and stress intensification
            #: factors. See Remark 3. (Integer = 1, 2, or 3)
            fsi = integer(card, 3, 'fsi')
            if fsi in [1, 2, 3]:
                # assuming MSC
                #: Mean cross-sectional radius of the curved pipe
                rm = double(card, 4, 'rm')

                #: Wall thickness of the curved pipe
                t = double(card, 5, 't')

                #: Internal pressure
                p = double_or_blank(card, 6, 'p')

                # line3
                # Non-structural mass :math:`nsm`
                nsm = double_or_blank(card, 11, 'nsm', default=0.)
                rc = double_or_blank(card, 12, 'rc', default=0.)
                zc = double_or_blank(card, 13, 'zc', default=0.)
            elif fsi in [4, 5, 6]:
                # Non-structural mass :math:`nsm`
                nsm = double_or_blank(card, 11, 'nsm', default=0.)
                rc = double_or_blank(card, 12, 'rc', default=0.)
                zc = double_or_blank(card, 13, 'zc', default=0.)

                #sacl = double_or_blank(card, 9, 'sacl')
                #alpha = double_or_blank(card, 10, 'alpha', 0.)
                #flange = integer_or_blank(card, 15, 'flange', 0)
                #kx = double_or_blank(card, 18, 'kx', 1.0)
                #ky = double_or_blank(card, 19, 'ky', 1.0)
                #kz = double_or_blank(card, 20, 'kz', 1.0)
                #sy = double_or_blank(card, 22, 'sy', 1.0)
                #sz = double_or_blank(card, 23, 'sz', 1.0)
            else:
                assert fsi in [1, 2, 3, 4, 5, 6], 'pid=%s fsi=%s\ncard:%s' % (pid, fsi, card)
        else:
            raise RuntimeError('Area/FSI on CBEND must be defined...')
        assert fsi in [0, 1, 2, 3, 4, 5, 6], 'pid=%s fsi=%s\ncard:%s' % (pid, fsi, card)

        #: Bend radius of the line of centroids
        rb = double_or_blank(card, 7, 'rb')

        #: Arc angle :math:`\theta_B` of element  (optional)
        theta_b = double_or_blank(card, 8, 'thetab')
        assert len(card) <= 23, f'len(PBEND card) = {len(card):d}\ncard={card}'
        #return PBEND(pid, mid, beam_type, A, i1, i2, j, c1, c2, d1, d2,
                     #e1, e2, f1, f2, k1, k2, nsm,
                     #rc, zc, delta_n, fsi, rm, t,
                     #p, rb, theta_b, comment=comment)
        self.cards.append((pid, mid, beam_type,
                           # type=1
                           A, i1, i2, j,
                           c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                           nsm,
                           # common
                           rc, zc,
                           # beam_type=1
                           delta_n,
                           # beam_type=2
                           fsi, rm, t, p,
                           rb, theta_b,
                           centerline_spacing, alpha, flange, kx, ky, kz, sy, sz,
                           comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype

        # common
        property_id = np.zeros(ncards, dtype=idtype)
        material_id = np.zeros(ncards, dtype=idtype)
        nsm = np.zeros(ncards, dtype='float64')
        beam_type = np.zeros(ncards, dtype='int32')
        fsi = np.zeros(ncards, dtype='int32')
        rb = np.zeros(ncards, dtype='float64')
        theta_b = np.zeros(ncards, dtype='float64')

        # beam_type = 1
        A = np.zeros(ncards, dtype='float64')
        J = np.zeros(ncards, dtype='float64')
        I1 = np.zeros(ncards, dtype='float64')
        I2 = np.zeros(ncards, dtype='float64')

        k1 = np.zeros(ncards, dtype='float64')
        k2 = np.zeros(ncards, dtype='float64')

        c1 = np.zeros(ncards, dtype='float64')
        c2 = np.zeros(ncards, dtype='float64')
        d1 = np.zeros(ncards, dtype='float64')
        d2 = np.zeros(ncards, dtype='float64')
        e1 = np.zeros(ncards, dtype='float64')
        e2 = np.zeros(ncards, dtype='float64')
        f1 = np.zeros(ncards, dtype='float64')
        f2 = np.zeros(ncards, dtype='float64')

        rc = np.zeros(ncards, dtype='float64')
        zc = np.zeros(ncards, dtype='float64')
        delta_n = np.zeros(ncards, dtype='float64')

        # beam_type = 2
        p = np.full(ncards, np.nan, dtype='float64')
        t = np.full(ncards, np.nan, dtype='float64')

        rm = np.full(ncards, np.nan, dtype='float64')
        kx = np.full(ncards, np.nan, dtype='float64')
        ky = np.full(ncards, np.nan, dtype='float64')
        kz = np.full(ncards, np.nan, dtype='float64')
        alpha = np.full(ncards, np.nan, dtype='float64')
        centerline_spacing = np.full(ncards, np.nan, dtype='float64')
        flange = np.full(ncards, np.nan, dtype='float64')
        sy = np.full(ncards, np.nan, dtype='float64')
        sz = np.full(ncards, np.nan, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, mid, beam_typei, areai, i1i, i2i, ji,
             c1i, c2i, d1i, d2i, e1i, e2i, f1i, f2i,
             k1i, k2i, nsmi,
             rci, zci, delta_ni, fsii, rmi, ti,
             pi, rbi, theta_bi,
             centerline_spacingi, alphai, flangei, kxi, kyi, kzi, syi, szi,
             comment) = card

            #assert beam_typei == 1, card
            # common
            property_id[icard] = pid
            material_id[icard] = mid
            nsm[icard] = nsmi
            beam_type[icard] = beam_typei
            rc[icard] = rci
            zc[icard] = zci
            rb[icard] = rbi

            # beam_type = 1
            if beam_typei == 1:
                A[icard] = areai
                J[icard] = ji
                I1[icard] = i1i
                I2[icard] = i2i
                theta_b[icard] = theta_bi

                k1[icard] = k1i
                k2[icard] = k2i

                c1[icard] = c1i
                c2[icard] = c2i
                d1[icard] = d1i
                d2[icard] = d2i
                e1[icard] = e1i
                e2[icard] = e2i
                f1[icard] = f1i
                f2[icard] = f2i
            elif beam_typei == 2:
                fsi[icard] = fsii
                rm[icard] = rmi
                p[icard] = pi
                t[icard] = ti

                delta_n[icard] = delta_ni

                # beam_type = 2
                rm[icard] = rmi

                flange[icard] = flangei
                alpha[icard] = alphai
                centerline_spacing[icard] = centerline_spacingi
                kx[icard] = kxi
                ky[icard] = kyi
                kz[icard] = kzi
                sy[icard] = syi
                sz[icard] = szi

        self._save(property_id, material_id, beam_type, nsm,
                   rb, theta_b, rc, zc,
                   # beam_type = 1
                   A, J, I1, I2,
                   c1, c2, d1, d2, e1, e2, f1, f2,
                   k1, k2, delta_n,
                   # beam_type = 2
                   fsi, rm, t, p,
                   centerline_spacing, alpha, flange,
                   kx, ky, kz, sy, sz,
                   )
        self.sort()

    def _save(self, property_id, material_id, beam_type, nsm,
              rb, theta_b, rc, zc,
              # beam_type = 1
              A, J, I1, I2,
              c1, c2, d1, d2, e1, e2, f1, f2,
              k1, k2, delta_n,
              # beam_type = 2
              fsi, rm, t, p,
              centerline_spacing, alpha, flange,
              kx, ky, kz, sy, sz) -> None:
        self.property_id = property_id
        self.material_id = material_id
        self.beam_type = beam_type
        self.rb = rb
        self.theta_b = theta_b
        self.rc = rc
        self.zc = zc
        self._nsm = nsm

        # beam_type = 1
        self.A = A
        self.J = J
        self.I1 = I1
        self.I2 = I2

        self.c1 = c1
        self.c2 = c2
        self.d1 = d1
        self.d2 = d2
        self.e1 = e1
        self.e2 = e2
        self.f1 = f1
        self.f2 = f2

        self.k1 = k1
        self.k2 = k2
        self.delta_n = delta_n

        # beam_type = 2
        self.fsi = fsi
        self.p = p
        self.t = t
        self.rm = rm
        self.centerline_spacing = centerline_spacing
        self.alpha = alpha
        self.flange = flange
        self.kx = kx
        self.ky = ky
        self.kz = kz
        self.sy = sy
        self.sz = sz
        self.n = len(property_id)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    def convert(self, xyz_scale: float=1.0,
                area_scale: float=1.0,
                area_inertia_scale:float=1.0,
                nsm_per_length_scale: float=1.0,
                pressure_scale: float=1.0,
                **kwargs):
        # common (easy)
        self.rc *= xyz_scale
        self.zc *= xyz_scale
        self.rb *= xyz_scale
        self._nsm *= nsm_per_length_scale

        # beam_type = 1 (easy)
        self.A *= area_scale
        self.J *= area_inertia_scale
        self.I1 *= area_inertia_scale
        self.I2 *= area_inertia_scale

        self.c1 *= xyz_scale
        self.c2 *= xyz_scale
        self.d1 *= xyz_scale
        self.d2 *= xyz_scale
        self.e1 *= xyz_scale
        self.e2 *= xyz_scale
        self.f1 *= xyz_scale
        self.f2 *= xyz_scale
        self.delta_n *= xyz_scale

        # beam_type = 2 (easy)
        self.p *= pressure_scale
        self.t *= xyz_scale
        self.centerline_spacing *= xyz_scale

        # ???
        #self.k1 = k1
        #self.k2 = k2

    def __apply_slice__(self, prop: PBEND, i: np.ndarray) -> None:
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        prop._nsm = self._nsm[i]
        prop.rc = self.rc[i]
        prop.zc = self.zc[i]

        prop.rb = self.rb[i]
        prop.theta_b = self.theta_b[i]

        # beam_type = 1
        prop.c1 = self.c1[i]
        prop.c2 = self.c2[i]
        prop.d1 = self.d1[i]
        prop.d2 = self.d2[i]
        prop.e1 = self.e1[i]
        prop.e2 = self.e2[i]
        prop.f1 = self.f1[i]
        prop.f2 = self.f2[i]

        prop.k1 = self.k1[i]
        prop.k2 = self.k2[i]
        prop.delta_n = self.delta_n[i]

        # beam_type = 2
        prop.fsi = self.fsi[i]
        prop.p = self.p[i]
        prop.t = self.t[i]

        prop.centerline_spacing = self.centerline_spacing[i]
        prop.alpha = self.alpha[i]
        prop.kx = self.kx[i]
        prop.ky = self.ky[i]
        prop.kz = self.kz[i]
        prop.sy = self.sy[i]
        prop.sz = self.sz[i]
        prop.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        materials = self.allowed_materials
        mids = hstack_msg([prop.material_id for prop in materials],
                          msg=f'no materials for {self.type}; {self.all_materials}')
        mids = np.unique(mids)
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    @property
    def all_materials(self) -> list[MAT1]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[MAT1]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def rho(self) -> np.ndarray:
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        return rho

    def mass_per_length(self) -> np.ndarray:
        mpl = np.full(self.n, np.nan, dtype='float64')
        raise NotImplementedError('PBEND')
        return mpl

    @property
    def is_small_field(self) -> bool:
        return max(self.property_id.max(),
                   self.material_id.max(), ) < 99_999_999

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:

        if size == 8 and self.is_small_field:
            print_card = print_card_8
        else:
            print_card = print_card_16

        ps = array_default_float(self.p, default=0., size=size, is_double=False)
        rbs = array_default_float(self.rb, default=0., size=size, is_double=False)
        for (pid, mid, beam_type, nsm, rc, zc,
             # beam_type=1
             area, j, i1, i2, rb, theta_b,
             c1, c2, d1, d2, e1, e2, f1, f2,
             k1, k2, delta_n,
             # beam_type=2
             fsi, rm, t, p, rb, theta_b,
             centerline_spacing, alpha, flange,
             kx, ky, kz, sy, sz) in zip_longest(self.property_id, self.material_id, self.beam_type,
                                                self._nsm, self.rc, self.zc,
                                                # beam_type=1
                                                self.A, self.J, self.I1, self.I2,
                                                self.rb, self.theta_b,
                                                self.c1, self.c2, self.d1, self.d2,
                                                self.e1, self.e2, self.f1, self.f2,
                                                self.k1, self.k2,  self.delta_n,
                                                # beam_type=2
                                                self.fsi, self.rm, self.t, ps, rbs, self.theta_b,
                                                self.centerline_spacing, self.alpha, self.flange,
                                                self.kx, self.ky, self.kz, self.sy, self.sz):

            list_fields = ['PBEND', pid, mid, ]  # other
            if beam_type == 1:
                list_fields += [
                    area, i1, i2, j, rb, theta_b,
                    c1, c2, d1, d2, e1, e2, f1, f2,
                    k1, k2, nsm, rc, zc, delta_n]
                #print("beam_type=0 I1=%s I2=%s; J=%s RM=%s T=%s P=%s" % (
                    #i1, i2, j, rm, t, p), list_fields)
            elif beam_type == 2:
                list_fields += [fsi, rm, t, p, rb, theta_b,
                                centerline_spacing, alpha, nsm, rc, zc, flange, None, None,
                                kx, ky, kz, None, sy, sz]
                #print(f'pid={pid} fsi={fsi} rm={rm} t={t} p={p.strip()} rb={rb.strip()} theta_b={theta_b} nsm={nsm} rc={rc} zc={zc}')
            elif beam_type == 0:
                raise RuntimeError(beam_type)
                # dunno
                #list_fields += [A, i1, i2, j, rb,
                                #theta_b, c1, c2, d1, d2,
                                #e1, e2, f1, f2, k1, k2,
                                #nsm, rc, zc, delta_n]

            bdf_file.write(print_card(list_fields))
        #af
        return

    def area(self) -> np.ndarray:
        area = np.full(self.n, np.nan, dtype='float64')
        ibeamtype1 = (self.beam_type == 1)
        area[ibeamtype1] = self.A
        return area
