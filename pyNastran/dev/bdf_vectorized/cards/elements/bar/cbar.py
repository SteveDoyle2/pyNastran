from numpy import array, arange, zeros, unique, searchsorted, full, nan
from numpy.linalg import norm  # type: ignore

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, string_or_blank)
from pyNastran.bdf.cards.elements.bars import BAROR

from pyNastran.bdf.field_writer_8 import set_string8_blank_if_default
from pyNastran.dev.bdf_vectorized.cards.elements.element import Element


class CBAR(Element):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    | CBAR  | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B  |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    | CBAR  | EID | PID | GA  | GB  | G0  |     |     | OFFT |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B  |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+

    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    |  CBAR | 2     |  39 | 7     | 6     |  105   |       |       |  GGG  |
    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    |       |       | 513 | 0.0+0 | 0.0+0 |    -9. | 0.0+0 | 0.0+0 |   -9. |
    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    """
    type = 'CBAR'
    def __init__(self, model):
        """
        Defines the CBAR object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        Element.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        self.n = ncards
        if self.n:
            assert isinstance(ncards, int), ncards
            float_fmt = self.model.float_fmt
            #: Element ID
            self.element_id = zeros(ncards, dtype='int32')
            #: Property ID
            self.property_id = zeros(ncards, dtype='int32')
            self.node_ids = zeros((ncards, 2), dtype='int32')
            self.is_g0 = zeros(ncards, dtype='bool')
            self.g0 = full(ncards, nan, dtype='int32')
            self.x = full((ncards, 3), nan, dtype=float_fmt)
            self.offt = full(ncards, nan,dtype= '|U3')
            self.pin_flags = zeros((ncards, 2), dtype='int32')
            self.wa = zeros((ncards, 3), dtype=float_fmt)
            self.wb = zeros((ncards, 3), dtype=float_fmt)

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i

        if 0 and self.model.cbaror.n > 0:
            cbaror = self.model.cbaror
            pid_default = cbaror.property_id
            is_g0_default = cbaror.is_g0
            x1_default = cbaror.x[0]
            x2_default = cbaror.x[1]
            x3_default = cbaror.x[2]
            g0_default = cbaror.g0
            offt_default = cbaror.offt
        else:
            pid_default = None
            is_g0_default = None
            x1_default = 0.0
            x2_default = 0.0
            x3_default = 0.0
            g0_default = None
            offt_default = 'GGG'

        eid = integer(card, 1, 'element_id')
        self.element_id[i] = eid
        if pid_default is not None:
            self.property_id[i] = integer_or_blank(card, 2, 'property_id', pid_default)
        else:
            self.property_id[i] = integer_or_blank(card, 2, 'property_id', eid)
        self.node_ids[i] = [integer(card, 3, 'GA'),
                            integer(card, 4, 'GB')]

        #---------------------------------------------------------
        # x / g0
        if g0_default is not None:
            field5 = integer_double_or_blank(card, 5, 'g0_x1', g0_default)
        else:
            field5 = integer_double_or_blank(card, 5, 'g0_x1', x1_default)

        if isinstance(field5, integer_types):
            self.is_g0[i] = True
            self.g0[i] = field5
        elif isinstance(field5, float):
            self.is_g0[i] = False
            x = array([field5,
                       double_or_blank(card, 6, 'x2', x2_default),
                       double_or_blank(card, 7, 'x3', x3_default)], dtype='float64')
            self.x[i, :] = x
            if norm(x) == 0.0:
                msg = 'G0 vector defining plane 1 is not defined on CBAR %s.\n' % eid
                msg += 'G0 = %s\n' % field5
                msg += 'X  = %s\n' % x
                msg += '%s' % card
                raise RuntimeError(msg)
        else:
            msg = ('field5 on CBAR (G0/X1) is the wrong type...id=%s field5=%s '
                   'type=%s' % (self.eid, field5, type(field5)))
            raise RuntimeError(msg)

        #---------------------------------------------------------
        # offt
        # bit doesn't exist on the CBAR
        offt = string_or_blank(card, 8, 'offt', offt_default)

        msg = 'invalid offt parameter of CBEAM...offt=%s' % offt
        assert offt[0] in ['G', 'B', 'O', 'E'], msg
        assert offt[1] in ['G', 'B', 'O', 'E'], msg
        assert offt[2] in ['G', 'B', 'O', 'E'], msg
        self.offt[i] = offt

        self.pin_flags[i, :] = [integer_or_blank(card, 9, 'pa', 0),
                                integer_or_blank(card, 10, 'pb', 0)]

        self.wa[i, :] = [double_or_blank(card, 11, 'w1a', 0.0),
                         double_or_blank(card, 12, 'w2a', 0.0),
                         double_or_blank(card, 13, 'w3a', 0.0),]

        self.wb[i, :] = [double_or_blank(card, 14, 'w1b', 0.0),
                         double_or_blank(card, 15, 'w2b', 0.0),
                         double_or_blank(card, 16, 'w3b', 0.0),]
        assert len(card) <= 17, 'len(CBAR card) = %i\ncard=%s' % (len(card), card)

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

            self.offt = self.offt[i]

            self.pin_flags = self.pin_flags[i, :]
            self.wa = self.wa[i, :]
            self.wb = self.wb[i, :]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CBAR IDs...')
            self._cards = []
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
            for i, (eid, pid, nids) in enumerate(zip(self.element_id, self.property_id,
                                                     self.node_ids)):
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
        if grid_cid0 is None:
            grid_cid0 = self.model.grid.get_position_by_node_index()
        p1 = grid_cid0[self.node_ids[:, 0]]
        p2 = grid_cid0[self.node_ids[:, 1]]
        L = p2 - p1
        i = self.model.properties_bar.get_index(self.property_id)
        A = self.model.properties_bar.get_Area[i]
        material_id = self.model.properties_bar.material_id[i]

        rho, E, J = self.model.Materials.get_rho_E_J(material_id)
        rho = self.model.Materials.get_rho(self.mid)
        E = self.model.Materials.get_E(self.mid)
        J = self.model.Materials.get_J(self.mid)

        mass = norm(L, axis=1) * A * rho + self.nsm
        if total:
            return mass.sum()
        else:
            return mass

    #=========================================================================
    def write_card(self, bdf_file, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, self.element_id)

            for (eid, pid, n, is_g0, g0, x, offt, pin, wa, wb) in zip(
                self.element_id[i], self.property_id[i], self.node_ids[i],
                self.is_g0[i], self.g0[i], self.x[i],
                self.offt[i],
                self.pin_flags[i], self.wa[i], self.wb[i]):

                pa = set_blank_if_default(pin[0], 0)
                pb = set_blank_if_default(pin[1], 0)

                w1a = set_blank_if_default(wa[0], 0.0)
                w2a = set_blank_if_default(wa[1], 0.0)
                w3a = set_blank_if_default(wa[2], 0.0)
                w1b = set_blank_if_default(wb[0], 0.0)
                w2b = set_blank_if_default(wb[1], 0.0)
                w3b = set_blank_if_default(wb[2], 0.0)
                x1 = g0 if is_g0 else x[0]
                x2 = 0 if is_g0 else x[1]
                x3 = 0 if is_g0 else x[2]

                offt = set_string8_blank_if_default(offt, 'GGG')
                card = ['CBAR', eid, pid, n[0], n[1], x1, x2, x3, offt,
                        pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        obj = CBAR(self.model)
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
        obj.offt = self.offt[i]
        obj.pin_flags = self.pin_flags[i]
        obj.wa = self.wa[i]
        obj.wb = self.wb[i]
        return obj

    #def get_stiffness_matrix(self, model, node_ids, index0s, fnorm=1.0):
        #return K, dofs, n_ijv
