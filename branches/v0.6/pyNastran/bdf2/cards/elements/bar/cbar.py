from numpy import array, dot, arange, zeros, unique, searchsorted

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, string_or_blank, blank)


class CBAR(object):
    type = 'CBAR'
    def __init__(self, model):
        """
        Defines the CBAR object.

        :param self: the CBAR object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment=None):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        """
        :param self: the CBAR object
        """
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float

            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 2), 'int32')
            #self.offt = array(ncards, '|S3')
            self.pin_flags = zeros((ncards, 2), 'int32')
            self.wa = zeros((ncards, 3), float_fmt)
            self.wb = zeros((ncards, 3), float_fmt)

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'element_id')
                self.property_id[i] = integer_or_blank(card, 2, 'property_id', self.element_id[i])
                self.node_ids[i] = [integer(card, 3, 'GA'),
                                    integer(card, 4, 'GB')]

                #self.initX_G0(card)

                #self.offt[i, :] = string_or_blank(card, 8, 'offt', 'GGG')
                #print 'self.offt = |%s|' % self.offt

                self.pin_flags[i, :] = [integer_or_blank(card, 9, 'pa', 0),
                                        integer_or_blank(card, 10, 'pb', 0)]

                self.wa[i, :] = [double_or_blank(card, 11, 'w1a', 0.0),
                                 double_or_blank(card, 12, 'w2a', 0.0),
                                 double_or_blank(card, 13, 'w3a', 0.0),]

                self.wb[i, :] = [double_or_blank(card, 14, 'w1b', 0.0),
                                 double_or_blank(card, 15, 'w2b', 0.0),
                                 double_or_blank(card, 16, 'w3b', 0.0),]
                assert len(card) <= 17, 'len(CBAR card) = %i' % len(card)

            i = self.element_id.argsort()
            print "i", i, type(i)
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            #self.offt = self.offt[i, :]
            self.pin_flags = self.pin_flags[i, :]
            self.wa = self.wa[i, :]
            self.wb = self.wb[i, :]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CBAR IDs...')
            self._cards = []
            self._comments = []

    #=========================================================================
    def get_mass(self, grid_cid0=None, total=False):
        """
        mass = rho * A * L + nsm
        """
        if self.n == 0:
            return 0.0

        if grid_cid0 is None:
            grid_cid0 = self.model.grid.position()
        p1 = grid_cid0[self.node_ids[:, 0]]
        p2 = grid_cid0[self.node_ids[:, 1]]
        L = p2 - p1
        i = self.model.properties_bar.get_index(self.property_id)
        A = self.model.properties_bar.get_Area[i]
        material_id = self.model.properties_bar.material_id[i]

        rho, E, J = self.model.Materials.get_rho_E_J(material_id)
        rho = self.model.Materials.get_rho(self.mid)
        E   = self.model.Materials.get_E(self.mid)
        J   = self.model.Materials.get_J(self.mid)

        mass = norm(L, axis=1) * A * rho + self.nsm
        if total:
            return mass.sum()
        else:
            return mass

    #=========================================================================
    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('CBAR', self.n))
        return msg

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, self.element_id)

            offt = 0
            ga = None
            gb = None
            x1 = None
            x2 = None
            x3 = None
            for (eid, pid, n, pin, wa, wb) in zip(self.element_id[i], self.property_id[i], self.node_ids[i], 
                    self.pin_flags[i], self.wa[i], self.wb[i]):

                pa = set_blank_if_default(pin[0], 0)
                pb = set_blank_if_default(pin[1], 0)

                w1a = set_blank_if_default(wa[0], 0.0)
                w2a = set_blank_if_default(wa[1], 0.0)
                w3a = set_blank_if_default(wa[2], 0.0)
                w1b = set_blank_if_default(wb[0], 0.0)
                w2b = set_blank_if_default(wb[1], 0.0)
                w3b = set_blank_if_default(wb[2], 0.0)
                #(x1, x2, x3) = self.getX_G0_defaults()
                #offt = set_blank_if_default(offt, 'GGG')
                #list_fields = ['CBAR', self.eid, self.Pid(), self.Ga(), self.Gb(), x1, x2,
                #          x3, offt, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
                card = ['CBAR', eid, pid, n[0], n[1], x1, x2, x3, offt,
                        pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
                f.write(print_card(card))


    def get_stiffness(self, model, node_ids, index0s, fnorm=1.0):
        return(K, dofs, nIJV)