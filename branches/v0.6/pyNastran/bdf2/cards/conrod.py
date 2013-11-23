from numpy import zeros, unique

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string)

class CONROD(object):
    type = 'CONROD'
    def __init__(self, model):
        """
        Defines the CONROD object.

        :param self: the CONROD object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._conrod = []
        self._conrod_comment = []

    def add(self, card, comment=None):
        self._conrod.append(card)
        self._conrod_comment.append(comment)

    def build(self):
        """
        :param self: the CONROD object
        :param cards: the list of CONROD cards
        """
        cards = self._conrod
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            #: Property ID
            self.eid = zeros(ncards, 'int32')
            self.mid = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 2), 'int32')
            self.A = zeros(ncards, float_fmt)
            self.J = zeros(ncards, float_fmt)
            self.c = zeros(ncards, float_fmt)
            self.nsm = zeros(ncards, float_fmt)

            for i, card in enumerate(cards):
                self.eid[i] = integer(card, 1, 'eid')
                self.node_ids[i] = [integer(card, 2, 'n1'),
                                    integer(card, 3, 'n2')]

                self.mid[i] = integer(card, 4, 'mid')
                self.A[i] = double(card, 5, 'A')
                self.J[i] = double_or_blank(card, 6, 'j', 0.0)
                self.c[i] = double_or_blank(card, 7, 'c', 0.0)
                self.nsm[i] = double_or_blank(card, 8, 'nsm', 0.0)
                assert len(card) <= 9, 'len(CONROD card) = %i' % len(card)

            i = self.eid.argsort()
            self.eid = self.eid[i]
            self.node_ids = self.node_ids[i, :]
            self.mid = self.mid[i]
            self.A = self.A[i]
            self.J = self.J[i]
            self.c = self.c[i]
            self.nsm = self.nsm[i]
        
            unique_eids = unique(self.eid)
            if len(unique_eids) != len(self.eid):
                raise RuntimeError('There are duplicate CROD IDs...')

    def displacement_stress(self, q):
        grid_cid0 = self.grid.position()
        p1 = grid_cid0[self.node_ids[:, 0]]
        p2 = grid_cid0[self.node_ids[:, 1]]
        L = p2 - p1
        n = norm(L, axis=1)
        
        v = divide_2d_array_by_column_vector(L, n)
        mass = n * self.A * rho + self.nsm


        A = self.A
        J = self.J
        c = self.c
        ki_torsion = G * J / L
        ki_axial   = A * E / L
        E = self.Materials.get_E(self.mid)
        k = array([[1, -1],
                   [-1, 1]])
        Kall = []
        for (vx, vy) in v[:, 1:2]:
            L = array([[vx, vy,  0,  0],
                         0,  0, vx, vy])
            K = dot(transpose(L), dot(k, L))

        #nijv = zeros((len(self.eid), 4)
        #nijv = [1, 2, 3, 4, 5, 6]
        #nijv[dof1, dof2, dof3, dof4, dof5, dof6]
        
        return K

    def mass(self, total=False):
        """
        mass = rho * A * L + nsm
        """
        grid_cid0 = self.grid.position()
        p1 = grid_cid0[self.node_ids[:, 0]]
        p2 = grid_cid0[self.node_ids[:, 1]]
        L = p2 - p1
        rho = self.model.Materials.get_rho(self.mid)
        mass = norm(L, axis=1) * self.A * rho + self.nsm
        if total:
            return mass.sum()
        else:
            return mass
        
    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('CROD', self.n))
        return msg

    def write_bdf(self, f, size=8, eids=None):
        if self.n:
            for (eid, n12, mid, A, J, c, nsm) in zip(
                 self.eid, self.node_ids, self.mid, self.A, self.J,
                 self.c, self.nsm):

                #self.mid = integer(card, 4, 'mid')
                #self.A = double(card, 5, 'A')
                #self.j = double_or_blank(card, 6, 'j', 0.0)
                #self.c = double_or_blank(card, 7, 'c', 0.0)
                #self.nsm = double_or_blank(card, 8, 'nsm', 0.0)

                card = ['CONROD', eid, n12[0], n12[1], mid, A, J, c, nsm ]
                f.write(print_card(card))