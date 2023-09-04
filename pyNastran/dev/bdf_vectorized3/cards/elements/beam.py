from __future__ import annotations
from itertools import count, zip_longest
from typing import Optional, Any, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import zip_strict, integer_types, float_types
from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16 # , print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string, blank,
    integer_or_blank, double_or_blank, string_or_blank,
    integer_double_or_blank, integer_string_or_blank, double_string_or_blank)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default # init_x_g0,
from pyNastran.bdf.cards.properties.bars import _bar_areaL # PBARL as pbarl, A_I1_I2_I12
from pyNastran.utils.mathematics import integrate_positive_unit_line # integrate_unit_line,

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, make_idim, hslice_by_idim, searchsorted_filter # vslice_by_idim,
)
from .rod import line_pid_mass_per_length, line_length, line_vector_length, line_centroid
from .bar import apply_bar_default, init_x_g0
from .utils import get_density_from_material
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran2.utils import hstack_msg

if TYPE_CHECKING:
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


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
    def add(self, eid: int, pid: int, nids: list[int],
            x: Optional[list[float]], g0: Optional[int],
            offt: str='GGG', bit=None,
            pa: int=0, pb: int=0,
            wa=None, wb=None,
            sa: int=0, sb: int=0, comment: str='') -> None:
        if wa is None:
            wa = [0., 0., 0.]
        if wb is None:
            wb = [0., 0., 0.]
        self.cards.append((eid, pid, nids, g0, x, offt, [pa, pb], wa, wb, sa, sb, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
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

    def parse_cards(self):
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 2), dtype='int32')
        offt = np.full(ncards, '', dtype='|U3')
        g0 = np.full(ncards, -1, dtype='int32')
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
            if g0i is None:
                x[icard, :] = xi
            else:
                g0[icard] = g0i
            offt[icard] = offti
            pa[icard] = pai
            pb[icard] = pbi
            wa[icard, :] = wai
            wb[icard, :] = wbi
            sa[icard] = sai
            sb[icard] = sbi

        self._save(element_id, property_id, nodes, offt, g0, x,
                   pa, pb, wa, wb, sa, sb)
        beamor = self.model.beamor
        apply_bar_default(self, beamor)
        self.cards = []

    def _save(self, element_id, property_id, nodes, offt, g0, x,
              pa, pb, wa, wb, sa, sb):
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.offt = offt
        self.g0 = g0
        self.x = x

        self.pa = pa
        self.pb = pb
        self.wa = wa
        self.wb = wb
        self.sa = sa
        self.sb = sb
        self.n = len(property_id)

    def __apply_slice__(self, elem: CBEAM, i: np.ndarray) -> None:
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

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no beam properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    def write(self, size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> str:
        if len(self.element_id) == 0:
            return ''

        lines = []
        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes_ = array_str(self.nodes, size=size)
        pas = array_default_int(self.pa, default=0, size=size)
        pbs = array_default_int(self.pb, default=0, size=size)
        for eid, pid, nodes, g0, x, offt, pa, pb, wa, wb in zip_longest(element_ids, property_ids, nodes_,
                                                                        self.g0, self.x, self.offt,
                                                                        pas, pbs, self.wa, self.wb):
            n1, n2 = nodes
            w1a = set_blank_if_default(wa[0], default=0.0)
            w2a = set_blank_if_default(wa[1], default=0.0)
            w3a = set_blank_if_default(wa[2], default=0.0)

            w1b = set_blank_if_default(wb[0], default=0.0)
            w2b = set_blank_if_default(wb[1], default=0.0)
            w3b = set_blank_if_default(wb[2], default=0.0)
            if g0 == -1:
                x1, x2, x3 = x # self.get_x_g0_defaults()
            else:
                x1 = g0
                x2 = ''
                x3 = ''

            # offt doesn't exist in NX nastran
            offt = set_blank_if_default(offt, 'GGG')

            list_fields = ['CBEAM', eid, pid, n1, n2,
                           x1, x2, x3, offt, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
            lines.append(print_card_8(list_fields))
        return ''.join(lines)

    @property
    def all_properties(self):
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

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

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

    def is_offt(self) -> np.ndarray:
        is_offt = (self.g0 == -1)
        return is_offt

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
        return self.n

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

        self.cards.append((pid, mid,
                           xxb, so, area, j, i1, i2, i12, nsm,
                           c1, c2, d1, d2, e1, e2, f1, f2,
                           s1, s2, k1, k2,
                           nsia, nsib, cwa, cwb,
                           m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b,
                           comment))

        self.n += 1
        return self.n

    def parse_cards(self):
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
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
        self.nsm = nsm

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
        prop.nsm = hslice_by_idim(i, istation, self.nsm)

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
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

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
        nproperties = len(self.property_id)
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        if rho.max() == 0. and rho.min() == 0. and self.nsm.max() == 0. and self.nsm.min() == 0.:
            return np.zeros(nproperties, dtype=rho.dtype)

        mass_per_length = np.zeros(nproperties, dtype='float64')
        for i, rhoi, istation in zip(count(), rho, self.istation):
            istation0, istation1 = istation
            assert istation1 > istation0
            xxb = self.xxb[istation0:istation1]
            area = self.A[istation0:istation1]
            nsm = self.nsm[istation0:istation1]
            mass_per_lengths = rhoi * area + nsm
            mass_per_lengthi = integrate_positive_unit_line(xxb, mass_per_lengths)
            mass_per_length[i] = mass_per_lengthi

        assert len(mass_per_length) == nproperties
        return mass_per_length

    @property
    def is_small_field(self):
        return max(self.property_id.max(),
                   self.material_id.max(),
                   self.s1.max(), self.s2.max()) < 99_999_999

    def write(self, size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> str:
        if len(self.property_id) == 0:
            return ''
        lines = []
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
            nsm_ = self.nsm[istation0:istation1]
            c1_ = self.c1[istation0:istation1]
            c2_ = self.c2[istation0:istation1]
            d1_ = self.d1[istation0:istation1]
            d2_ = self.d2[istation0:istation1]
            e1_ = self.e1[istation0:istation1]
            e2_ = self.e2[istation0:istation1]
            f1_ = self.f1[istation0:istation1]
            f2_ = self.f2[istation0:istation1]
            assert nstation > 0

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
            lines.append(print_card(list_fields))
        return ''.join(lines)

    @property
    def istation(self) -> np.ndarray:
        return make_idim(self.n, self.nstation)

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

    def __init__(self, model: BDF):
        super().__init__(model)
        #prop.property_id = self.property_id[i]
        self.material_id = np.array([], dtype='int32')

        #self.istation = hslice_by_idim(i, idim, elements)
        #idim = self.idim
        self.Type = np.array([], dtype='|U8')
        self.group = np.array([], dtype='|U8')
        self.nsm = np.array([], dtype='float64')

        self.dims = np.zeros([], dtype='float64')

        self.idim = np.zeros((0, 2), dtype='int32')  # for all properties
        self.ndim = np.array([], dtype='int32')

        self.istation = np.zeros((0, 2), dtype='int32')
        self.nstation = np.array([], dtype='int32')
        #self.ndim =

    def slice_card_by_property_id(self, property_id: np.ndarray) -> PBEAML:
        """uses a node_ids to extract GRIDs"""
        iprop = self.index(property_id)
        #assert len(self.node_id) > 0, self.node_id
        #i = np.searchsorted(self.node_id, node_id)
        prop = self.slice_card_by_index(iprop)
        return prop

    def add(self, pid: int, mid: int, beam_type: str,
            xxb, dims, so=None, nsm=None,
            group: str='MSCBML0', comment: str='') -> None:
        nxxb = len(xxb)
        if so is None:
            so = ['YES'] * nxxb
        elif isinstance(so, str):
            so = [so] * nxxb

        if nsm is None:
            nsm = [0.] * nxxb
        elif isinstance(nsm, float_types):
            nsm = [nsm] * nxxb
        ndim = self.valid_types[beam_type]
        self.cards.append((pid, mid, beam_type, group, xxb, so, nsm, ndim, dims, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        group = string_or_blank(card, 3, 'group', default='MSCBML0')
        beam_type = string(card, 4, 'Type')

        # determine the number of required dimensions on the PBEAM
        ndim = self.valid_types[beam_type]

        #: dimension list
        dims = []
        dim = []

        #: Section position
        xxb = [0.]

        #: Output flag
        so = ['YES']  # station 0

        #: non-structural mass :math:`nsm`
        nsm = []

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
        self.cards.append((pid, mid, beam_type, group, xxb, so, nsm, ndim, dims, comment))
        self.n += 1

    def parse_cards(self):
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return

        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')

        idim = np.zeros((ncards, 2), dtype='int32')  # for all properties
        ndim = np.zeros(ncards, dtype='int32')

        istation = np.zeros((ncards, 2), dtype='int32')  # for all properties
        nstation = np.zeros(ncards, dtype='int32')

        Type = np.full(ncards, '', dtype='|U8')
        group = np.full(ncards, '', dtype='|U8')
        #nsm = np.zeros(ncards, dtype='float64')

        all_xxb = []
        all_dims = []
        all_so = []
        all_nsm = []

        idim0 = 0
        istation0 = 0
        for icard, card in enumerate(self.cards):
            (pid, mid, beam_type, groupi, xxbi, soi, nsmi, ndimi, dims, comment) = card

            nstationi = len(xxbi)
            all_xxb.extend(xxbi)
            for dim in dims:
                all_dims.extend(dim)
            all_so.extend(soi)
            all_nsm.extend(nsmi)
            #station.extend(xxbi)

            idim1 = idim0 + ndimi * nstationi
            istation1 = istation0 + nstationi
            nstation[icard] = nstationi
            property_id[icard] = pid
            material_id[icard] = mid
            group[icard] = groupi
            Type[icard] = beam_type
            ndim[icard] = ndimi

            idim[icard, :] = [idim0, idim1]
            istation[icard, :] = [istation0, istation1]

            #self.nsm[icard] = nsm
            idim0 = idim1
            istation0 = istation1

        xxb = np.array(all_xxb, dtype='float64')
        dims = np.array(all_dims, dtype='float64')
        so = np.array(all_so, dtype='|U4')
        nsm = np.array(all_nsm, dtype='float64')
        #ndim_total = self.ndim.sum()
        self._save(property_id, material_id, idim, ndim, istation, nstation, Type, group,
                   xxb, dims, so, nsm)
        nstation_total = self.nstation.sum()
        self.sort()

        assert self.xxb.shape[0] == nstation_total
        #assert len(self.dims) == ndim_total
        assert len(self.so) == nstation_total
        assert len(self.nsm) == nstation_total
        self.cards = []

    def _save(self, property_id, material_id, idim, ndim, istation, nstation, Type, group,
              xxb, dims, so, nsm):
        self.property_id = property_id
        self.material_id = material_id

        assert idim.ndim == 2, idim.shape
        assert idim.min() == 0, idim
        assert ndim.min() >= 1, idim
        self.idim = idim
        self.ndim = ndim

        assert istation.ndim == 2, istation.shape
        assert istation.min() == 0, istation
        assert nstation.min() >= 2, nstation
        self.istation = istation
        self.nstation = nstation

        self.Type = Type
        self.group = group

        self.nsm = nsm
        self.xxb = xxb
        self.dims = dims
        self.so = so
        self.nsm = nsm

    def __apply_slice__(self, prop: PBEAML, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]

        #self.istation = hslice_by_idim(i, idim, elements)
        prop.Type = self.Type[i]
        prop.group = self.group[i]
        prop.nsm = self.nsm[i]

        idim = self.idim
        istation = self.istation
        prop.dims = hslice_by_idim(i, idim, self.dims)
        prop.xxb = hslice_by_idim(i, istation, self.xxb)

        #self.idim = np.zeros((ncards, 2), dtype='int32')  # for all properties
        prop.idim = self.idim[i, :]
        prop.ndim = self.ndim[i]

        #self.istation = istation
        prop.istation = self.istation[i, :]
        prop.nstation = self.nstation[i]

        prop.so = self.so[i]
        prop.nsm = self.nsm[i]

    #@property
    #def idim(self) -> np.ndarray:
        #return make_idim(self.n, self.ndim)

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          msg=f'no materials for {self.type}; {self.all_materials}')
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    def write(self, size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> str:
        if len(self.property_id) == 0:
            return ''
        lines = []
        if size == 8:
            print_card = print_card_8

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
            nsm = self.nsm[istation0:istation1]
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
            lines.append(print_card(list_fields))
        assert len(lines) > 0, lines
        return ''.join(lines)

    def area(self) -> np.ndarray:
        nproperties = len(self.property_id)
        area = np.zeros(nproperties, dtype='float64')
        for i, beam_type, ndim, idim, nstation, istation in zip(count(), self.Type,
                                                                self.ndim, self.idim,
                                                                self.nstation, self.istation,):
            idim0, idim1 = idim
            istation0, istation1 = istation
            xxb = self.xxb[istation0:istation1]
            #nsm = self.nsm[istation0:istation1]
            #so = self.so[istation0:istation1]
            dims = self.dims[idim0 : idim1].reshape(nstation, ndim)

            #prop = pbarl(self.property_id[i], self.material_id[i], beam_type, dim)
            #A, I1, I2, I12 = A_I1_I2_I12(prop, beam_type, dim)
            areasi = []
            for dim in dims:
                areai = _bar_areaL('PBEAML', beam_type, dim, self)
                areasi.append(areai)
            A = integrate_positive_unit_line(xxb, areasi)
            area[i] = A
        return area

    def rho(self) -> np.ndarray:
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        return rho

    @property
    def all_materials(self) -> list[Any]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[Any]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.materials}'
        return materials

    def mass_per_length(self):
        #nproperties = len(self.property_id)
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        A = self.area()
        nsm = 0.
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
    def __init__(self, model: BDF):
        super().__init__(model)
        #self.property_id = np.array([], dtype='int32')
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
            symopt: int=0, comment: str='') -> PBCOMP:
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

    def add_card(self, card: BDFCard, comment: str='') -> None:
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
            yi = double(card, i, 'y' + str(row))
            zi = double(card, i + 1, 'z' + str(row))
            ci = double_or_blank(card, i + 2, 'c' + str(row), 0.0)
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

    def parse_cards(self):
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return

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

        #self.Type = np.full(ncards, '', dtype='|U8')
        #self.group = np.full(ncards, '', dtype='|U8')
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
        self.property_id = self.property_id[i]
        self.material_id = self.material_id[i]

        self.Type = self.Type[i]
        self.group = self.group[i]
        self.k = self.k[i, :]
        self.m1 = self.m1[i]
        self.m2 = self.m2[i]
        self.n1 = self.n1[i]
        self.n2 = self.n2[i]

        istation = self.istation
        self.y = hslice_by_idim(i, istation, self.y)
        self.z = hslice_by_idim(i, istation, self.z)
        self.c = hslice_by_idim(i, istation, self.c)
        self.material_ids = hslice_by_idim(i, istation, self.material_ids)

        self.c1 = hslice_by_idim(i, istation, self.c1)
        self.c2 = hslice_by_idim(i, istation, self.c2)
        self.d1 = hslice_by_idim(i, istation, self.d1)
        self.d2 = hslice_by_idim(i, istation, self.d2)
        self.e1 = hslice_by_idim(i, istation, self.e1)
        self.e2 = hslice_by_idim(i, istation, self.e2)
        self.f1 = hslice_by_idim(i, istation, self.f1)
        self.f2 = hslice_by_idim(i, istation, self.f2)
        self.nstation = self.nstation[i, :]

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

    def write(self, size: int=8):
        if len(self.property_id) == 0:
            return ''
        lines = []
        if size == 8:
            print_card = print_card_8
        else:
            print_card = print_card_16


        property_ids = array_str(self.property_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        for pid, mid, A, i1, i2, i12, j, nsm, \
            (k1, k2), m1, m2, n1, n2, symopt, (istation0, istation1) in zip_longest(
            property_ids, material_ids, self._area, self.i1, self.i2, self.i12, self.j, self.nsm,
            self.k, self.m1, self.m2, self.n1, self.n2, self.symopt, self.istation):

            area = set_blank_if_default(A, 0.0)
            j = set_blank_if_default(j, 0.0)
            i1 = set_blank_if_default(i1, 0.0)
            i2 = set_blank_if_default(i2, 0.0)
            i12 = set_blank_if_default(i12, 0.0)
            nsm = set_blank_if_default(nsm, 0.0)

            k1 = set_blank_if_default(k1, 1.0)
            k2 = set_blank_if_default(k2, 1.0)

            m1 = set_blank_if_default(m1, 0.0)
            m2 = set_blank_if_default(m2, 0.0)

            n1 = set_blank_if_default(n1, 0.0)
            n2 = set_blank_if_default(n2, 0.0)

            symopt = set_blank_if_default(symopt, 0)

            list_fields = ['PBCOMP', pid, mid, area, i1, i2, i12, j,
                           nsm, k1, k2, m1, m2, n1, n2, symopt, None]

            y = self.y[istation0:istation1]
            z = self.z[istation0:istation1]
            c = self.c[istation0:istation1]
            mids = self.material_ids[istation0:istation1]
            for (yi, zi, ci, mid) in zip_longest(y, z, c, mids):
                ci = set_blank_if_default(ci, 0.0)
                list_fields += [yi, zi, ci, mid, None, None, None, None]
            #assert len(y) > 0, list_fields
            lines.append(print_card(list_fields))
        return ''.join(lines)

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
