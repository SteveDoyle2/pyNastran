from numpy import dot, arange, zeros, unique, searchsorted

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)


class CROD(object):
    type = 'CROD'
    def __init__(self, model):
        """
        Defines the CROD object.

        :param self: the CROD object
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
        :param self: the CROD object
        """
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            #: Property ID
            self.element_id = zeros(ncards, 'int32')
            self.property_id = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 2), 'int32')

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'element_id')
                self.property_id[i] = integer_or_blank(card, 2, 'property_id', self.element_id[i])
                self.node_ids[i] = [integer(card, 3, 'n1'),
                                    integer(card, 4, 'n2')]
                assert len(card) == 5, 'len(CROD card) = %i' % len(card)

            i = self.element_id.argsort()
            print "i", i, type(i)
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CROD IDs...')
            self._cards = []
            self._comments = []

    def mass(self, total=False):
        """
        mass = rho * A * L + nsm
        """
        grid_cid0 = self.grid.position()
        p1 = grid_cid0[self.node_ids[:, 0]]
        p2 = grid_cid0[self.node_ids[:, 1]]
        L = p2 - p1
        i = self.model.prod.get_index(pid)
        A = self.model.prod.A[i]
        mid = self.model.prod.mid[i]

        rho, E, J = self.model.Materials.get_rho_E_J(self.mid)
        rho = self.model.Materials.get_rho(self.mid)
        E   = self.model.Materials.get_E(self.mid)
        J   = self.model.Materials.get_J(self.mid)

        mass = norm(L, axis=1) * A * rho + self.nsm
        if total:
            return mass.sum()
        else:
            return mass

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('CROD', self.n))
        return msg

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, self.element_id)

            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i]):
                card = ['CROD', eid, pid, n[0], n[1] ]
                f.write(print_card(card))

    def get_stiffness(self, model, node_ids, index0s, fnorm=1.0):
        #print("----------------")
        A = self.Area()
        E = self.E()
        G = self.G()
        J = self.J()

        #========================
        #(n1, n2) = self.nodeIDs()
        n0, n1 = self.nodeIDs()

        i0, i1 = index0s
        node0 = self.nodes[0]
        node1 = self.nodes[1]


        p0 = model.Node(n0).xyz
        p1 = model.Node(n1).xyz
        L = norm(p0 - p1)
        if L == 0.0:
            msg = 'invalid CROD length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)
        #========================
        print("A=%g E=%g G=%g J=%g L=%g" % (A, E, G, J, L))
        k_axial = A * E / L
        k_torsion = G * J / L
        #k_axial = 1.0
        #k_torsion = 2.0

        k = matrix([[1., -1.], [-1., 1.]])  # 1D rod

        Lambda = self.Lambda(model)
        K = dot(dot(transpose(Lambda), k), Lambda)
        Ki, Kj = K.shape

        # for testing
        #K = ones((Ki, Ki), 'float64')

        K2 = zeros((Ki*2, Kj*2), 'float64')
        if k_axial == 0.0 and k_torsion == 0.0:
            dofs = []
            nIJV = []
            K2 = []
        elif k_torsion == 0.0: # axial; 2D or 3D
            K2 = K * k_axial
            dofs = array([
                i0, i0+1, i0+2,
                i1, i1+1, i1+2,
            ], 'int32')
            nIJV = [
                # axial
                (n0, 1), (n0, 2), (n0, 3),
                (n1, 1), (n1, 2), (n1, 3),
            ]
        elif k_axial == 0.0: # torsion; assume 3D
            K2 = K * k_torsion
            dofs = array([
                i0+3, i0+4, i0+5,
                i1+3, i1+4, i1+5,
            ], 'int32')
            nIJV = [
                # torsion
                (n0, 4), (n0, 5), (n0, 6),
                (n1, 4), (n1, 5), (n1, 6),
            ]

        else:  # axial + torsion; assume 3D
            # u1fx, u1fy, u1fz, u2fx, u2fy, u2fz
            K2[:Ki, :Ki] = K * k_axial

            # u1mx, u1my, u1mz, u2mx, u2my, u2mz
            K2[Ki:, Ki:] = K * k_torsion

            dofs = array([
                i0, i0+1, i0+2,
                i1, i1+1, i1+2,

                i0+3, i0+4, i0+5,
                i1+3, i1+4, i1+5,
            ], 'int32')
            nIJV = [
                # axial
                (n0, 1), (n0, 2), (n0, 3),
                (n1, 1), (n1, 2), (n1, 3),

                # torsion
                (n0, 4), (n0, 5), (n0, 6),
                (n1, 4), (n1, 5), (n1, 6),
            ]

        #Fg = dot(dot(transpose(Lambda), grav), Lambda)
        #print("K=\n", K / fnorm)
        #print("K2=\n", K2 / fnorm)

        #========================

        #print(K / fnorm)
        #print("K[%s] = \n%s\n" % (self.eid, list_print(K/fnorm)))

        print('dofs =', dofs)
        print('K =\n', list_print(K / fnorm))

        return(K2, dofs, nIJV)