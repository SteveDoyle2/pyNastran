from numpy import array, dot, arange, zeros, unique, searchsorted

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)


class CONM2(object):
    type = 'CONM2'
    def __init__(self, model):
        """
        Defines the CONM2 object.

        :param self: the CONM2 object
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
        :param self: the CONM2 object
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
            self.node_id = zeros(ncards, 'int32')
            self.coord_id = zeros(ncards, 'int32')
            self.mass = zeros(ncards, float_fmt)
            self.X = zeros((ncards, 3), float_fmt)
            self.I = zeros((ncards, 6), float_fmt)

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'element_id')
                self.node_id[i] = integer(card, 2, 'node_id')
                self.coord_id[i] = integer_or_blank(card, 3, 'coord_id', 0)
                self.mass[i] = double_or_blank(card, 4, 'mass', 0.)
                self.X[i, :] = [double_or_blank(card, 5, 'x1', 0.0),
                                double_or_blank(card, 6, 'x2', 0.0),
                                double_or_blank(card, 7, 'x3', 0.0)]

                self.I[i, :] = [double_or_blank(card, 9, 'I11', 0.0),
                                double_or_blank(card, 10, 'I21', 0.0),
                                double_or_blank(card, 11, 'I22', 0.0),
                                double_or_blank(card, 12, 'I31', 0.0),
                                double_or_blank(card, 13, 'I32', 0.0),
                                double_or_blank(card, 14, 'I33', 0.0)]
                assert len(card) <= 15, 'len(CONM2 card) = %i' % len(card)


                self.node_ids[i] = [integer(card, 3, 'n1'),
                                    integer(card, 4, 'n2')]
                assert len(card) == 5, 'len(CONM2 card) = %i' % len(card)

            i = self.element_id.argsort()
            print "i", i, type(i)
            self.element_id = self.element_id[i]
            self.node_id = self.node_id[i]
            self.coord_id = self.coord_id[i]
            self.mass = self.mass[i]
            self.X = self.X[i, :]
            self.I = self.I[i, :]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CONM2 IDs...')
            self._cards = []
            self._comments = []

    def mass(self, total=False):
        """
        mass = rho * A * L + nsm
        """
        asdf
        grid_cid0 = self.grid.position()
        p1 = grid_cid0[self.node_ids[:, 0]]
        p2 = grid_cid0[self.node_ids[:, 1]]

        mass = norm(L, axis=1) * A * rho + self.nsm
        if total:
            return mass.sum()
        else:
            return mass

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('CONM2', self.n))
        return msg

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, self.element_id)

            cid = [cid if cid != 0 else '' for cid in self.coord_id]

            X0 = [x if x != 0.0 else '' for x in self.x[i, 0]]
            X1 = [x if x != 0.0 else '' for x in self.x[i, 1]]
            X2 = [x if x != 0.0 else '' for x in self.x[i, 2]]
            Mass = [x if x != 0.0 else '' for x in self.mass[i]]

            I0 = [x if x != 0.0 else '' for x in self.I[i, 0]]
            I1 = [x if x != 0.0 else '' for x in self.I[i, 1]]
            I2 = [x if x != 0.0 else '' for x in self.I[i, 2]]
            I3 = [x if x != 0.0 else '' for x in self.I[i, 3]]
            I= [x if x != 0.0 else '' for x in self.I[i, 4]]
            I5 = [x if x != 0.0 else '' for x in self.I[i, 5]]
            for (eid, n, cid, mass, x0, x1, x2, i0, i1, i2, i3, i4, i5) in zip(self.element_id[i], self.node_id[i],
                    cid, Mass, X0, X1, X2, I0, I1, I2, I3, I4, I5):
                card = ['CONM2', eid, pid, n, cid, mass, x0, x1, x2,
                        None, i0, i1, i2, i3, i4, i5]
                f.write(print_card(card))

    def get_stiffness(self, model, node_ids, index0s, fnorm=1.0):
        return(K, dofs, nIJV)