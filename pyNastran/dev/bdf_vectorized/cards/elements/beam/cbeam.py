"""
http://www.ce.memphis.edu/7117/notes/presentations/chapter_04b.pdf
"""
#import numpy as np
from numpy import array, arange, zeros, unique, searchsorted, nan, full
from numpy.linalg import norm  # type: ignore

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_double_or_blank,
    double_string_or_blank)

from pyNastran.dev.bdf_vectorized.cards.elements.element import Element


class CBEAM(Element):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    | CBEAM | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |
    +-------+-----+-----+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    | CBEAM | EID | PID | GA  | GB  | G0  |     |     | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |
    +-------+-----+-----+

    """
    type = 'CBEAM'
    def __init__(self, model):
        """
        Defines the CBEAM object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        Element.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 2), 'int32')
            self.is_g0 = zeros(ncards, 'bool')
            self.g0 = full(ncards, nan, 'int32')
            self.x = full((ncards, 3), nan, float_fmt)
            self.is_offt = zeros(ncards, 'bool')
            self.bit = full(ncards, nan, 'int32')
            self.offt = full(ncards, nan, '|U3')
            self.pin_flags = zeros((ncards, 2), 'int32')
            self.wa = zeros((ncards, 3), float_fmt)
            self.wb = zeros((ncards, 3), float_fmt)
            self.sa = zeros(ncards, 'int32')
            self.sb = zeros(ncards, 'int32')

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i
        eid = integer(card, 1, 'element_id')
        self.element_id[i] = eid
        self.property_id[i] = integer_or_blank(card, 2, 'property_id', eid)
        self.node_ids[i] = [integer(card, 3, 'GA'),
                            integer(card, 4, 'GB')]

        #---------------------------------------------------------
        # x / g0
        field5 = integer_double_or_blank(card, 5, 'g0_x1', 0.0)
        if isinstance(field5, integer_types):
            self.is_g0[i] = True
            self.g0[i] = field5
        elif isinstance(field5, float):
            self.is_g0[i] = False
            x = array([field5,
                       double_or_blank(card, 6, 'x2', 0.0),
                       double_or_blank(card, 7, 'x3', 0.0)])
            self.x[i, :] = x
            if norm(x) == 0.0:
                msg = 'G0 vector defining plane 1 is not defined on CBEAM eid=%s.\n' % eid
                msg += 'G0 = %s\n' % field5
                msg += 'X  = %s\n' % x
                msg += '%s' % card
                raise RuntimeError(msg)
        else:
            msg = ('field5 on CBEAM (G0/X1) is the wrong type...id=%s field5=%s '
                   'type=%s' % (self.eid, field5, type(field5)))
            raise RuntimeError(msg)

        #---------------------------------------------------------
        # offt/bit
        field8 = double_string_or_blank(card, 8, 'offt/bit', 'GGG')
        if isinstance(field8, float):
            self.is_offt[i] = False
            self.bit[i] = field8
        elif isinstance(field8, str):
            self.is_offt[i] = True
            offt = field8
            msg = 'invalid offt parameter of CBEAM...offt=%s' % offt
            assert offt[0] in ['G', 'B', 'O', 'E'], msg
            assert offt[1] in ['G', 'B', 'O', 'E'], msg
            assert offt[2] in ['G', 'B', 'O', 'E'], msg
            self.offt[i] = offt
        else:
            msg = ('field8 on %s (offt/bit) is the wrong type...id=%s field5=%s '
                   'type=%s' % (self.type, self.eid, field8, type(field8)))
            raise RuntimeError(msg)

        self.pin_flags[i, :] = [integer_or_blank(card, 9, 'pa', 0),
                                integer_or_blank(card, 10, 'pb', 0)]

        self.wa[i, :] = [double_or_blank(card, 11, 'w1a', 0.0),
                         double_or_blank(card, 12, 'w2a', 0.0),
                         double_or_blank(card, 13, 'w3a', 0.0),]

        self.wb[i, :] = [double_or_blank(card, 14, 'w1b', 0.0),
                         double_or_blank(card, 15, 'w2b', 0.0),
                         double_or_blank(card, 16, 'w3b', 0.0),]
        self.sa[i] = integer_or_blank(card, 17, 'sa', 0)
        self.sb[i] = integer_or_blank(card, 18, 'sb', 0)
        assert len(card) <= 19, 'len(CBEAM card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]

            self.is_g0 = self.is_g0[i]
            self.g0 = self.g0[i]
            self.x = self.x[i, :]

            self.is_offt = self.is_offt[i]
            self.offt = self.offt[i]
            self.bit = self.bit[i]

            self.pin_flags = self.pin_flags[i, :]
            self.wa = self.wa[i, :]
            self.wb = self.wb[i, :]

            self.sa = self.sa[i]
            self.sb = self.sb[i]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CBEAM IDs...')
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'node_id' : nid_map,
            'property' : pid_map,
        }
        """
        if self.n:
            eid_map = maps['element']
            nid_map = maps['node']
            pid_map = maps['property']
            for i, (eid, pid, nids) in enumerate(zip(self.element_id, self.property_id, self.node_ids)):
                self.element_id[i] = eid_map[eid]
                self.property_id[i] = pid_map[pid]
                self.node_ids[i, 0] = nid_map[nids[0]]
                self.node_ids[i, 1] = nid_map[nids[1]]

    #=========================================================================
    def get_mass_by_element_id(self, grid_cid0=None, total=False):
        """
        mass = rho * A * L + nsm
        """
        if self.n == 0:
            return 0.0
        return [0.0]
        #if grid_cid0 is None:
            #grid_cid0 = self.model.grid.get_position_by_node_index()
        #p1 = grid_cid0[self.node_ids[:, 0]]
        #p2 = grid_cid0[self.node_ids[:, 1]]
        #L = p2 - p1
        #i = self.model.properties_bar.get_index(self.property_id)
        #A = self.model.properties_bar.get_Area[i]
        #material_id = self.model.properties_bar.material_id[i]

        #rho, E, J = self.model.Materials.get_rho_E_J(material_id)
        #rho = self.model.Materials.get_rho(self.mid)
        #E = self.model.Materials.get_E(self.mid)
        #J = self.model.Materials.get_J(self.mid)

        #mass = norm(L, axis=1) * A * rho + self.nsm
        #if total:
            #return mass.sum()
        #else:
            #return mass

    #=========================================================================
    def write_card(self, bdf_file, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, self.element_id)

            for (eid, pid, n, is_g0, g0, x, is_offt, offt, bit, pin, wa, wb, sa, sb) in zip(
                    self.element_id[i], self.property_id[i], self.node_ids[i],
                    self.is_g0[i], self.g0[i], self.x[i],
                    self.is_offt[i], self.offt[i], self.bit[i],
                    self.pin_flags[i], self.wa[i], self.wb[i], self.sa[i], self.sb[i]):

                x1 = g0 if is_g0 else x[0]
                x2 = 0 if is_g0 else x[1]
                x3 = 0 if is_g0 else x[2]
                offt_bit = offt if is_offt else bit
                #print('is_offt=%s offt=%s bit=%s offt_bit=%s' % (is_offt, offt, bit, offt_bit))

                pa = set_blank_if_default(pin[0], 0)
                pb = set_blank_if_default(pin[1], 0)

                w1a = set_blank_if_default(wa[0], 0.0)
                w2a = set_blank_if_default(wa[1], 0.0)
                w3a = set_blank_if_default(wa[2], 0.0)
                w1b = set_blank_if_default(wb[0], 0.0)
                w2b = set_blank_if_default(wb[1], 0.0)
                w3b = set_blank_if_default(wb[2], 0.0)

                sa = set_blank_if_default(sa, 0)
                sb = set_blank_if_default(sb, 0)
                card = ['CBEAM', eid, pid, n[0], n[1], x1, x2, x3, offt_bit,
                        pa, pb, w1a, w2a, w3a, w1b, w2b, w3b,
                        sa, sb]
                bdf_file.write(print_card_8(card))

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        obj = CBEAM(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.property_id = self.property_id[i]
        obj.node_ids = self.node_ids[i, :]
        obj.is_g0 = self.is_g0[i]
        obj.g0 = self.g0[i]
        obj.x = self.x[i, :]
        obj.is_offt = self.is_offt[i]
        obj.offt = self.offt[i]
        obj.bit = self.bit[i]
        obj.pin_flags = self.pin_flags[i]
        obj.wa = self.wa[i]
        obj.wb = self.wb[i]
        obj.sa = self.sa[i]
        obj.sb = self.sb[i]
        return obj

    #def get_stiffness_matrix(self, model, node_ids, index0s, fnorm=1.0):
        #K = np.zeros((12, 12), dtype='float64')
        #kaxial = E * A / L
        #ktorsion = G * J / L
        #K[0, 0] = K[6, 6] = kaxial
        #K[0, 6] = K[6, 0] = -kaxial

        #K[3, 3] = K[9, 9] = ktorsion
        #K[3, 9] = K[9, 3] = ktorsion

        ## Iy, L^3
        #eiy_3 = E * Iy / L**3
        #eiy_2 = E * Iy / L**2
        #eiy_1 = E * Iy / L

        #eiz_3 = E * Iz / L**3
        #eiz_2 = E * Iz / L**2
        #eiz_1 = E * Iz / L

        #K[1, 1] = K[7, 7] = 12 * eiz_3
        #K[2, 2] = K[8, 8] = 12 * eiy_3

        #K[4, 4] = K[10, 10] = 4 * eiy_1
        #K[5, 5] = K[11, 11] = 4 * eiz_1

        #K[4, 10] = K[10, 4] = 2 * eiy_1
        #K[5, 11] = K[11, 5] = 2 * eiz_1

        #K[4, 2] = K[8, 10] = K[10, 8] = K[4, 8] = K[8, 4] = 6 * eiy_2
        #K[2, 4] = K[2, 10] = K[10, 2] = -6 * eiy_2

        #K[1, 5] = K[1, 11] = K[5, 1] = K[11, 1] = 6 * eiz_2
        #K[5, 7] = K[11, 7] = K[7, 5] = K[7, 11] = -6 * eiz_2

        #K[1, 1] = K[7, 7] = 12 * eiz_3
        #K[2, 2] = K[8, 8] = 12 * eiy_3
        #K[2, 8] = K[8, 2] = -12 * eiy_3
        #K[1, 7] = K[7, 1] = -12 * eiz_3

        #return K, dofs, n_ijv
