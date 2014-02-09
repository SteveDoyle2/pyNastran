from numpy import array, dot, arange, zeros, unique, searchsorted, transpose
from numpy.linalg import norm

from pyNastran.bdf2.cards.elements.rod.conrod import _Lambda
from pyNastran.utils import list_print

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

    #=========================================================================
    def get_Area(self, property_ids):
        A = self.model.prod.get_Area(property_ids)
        return A

    def get_E(self, property_ids):
        E = self.model.prod.get_E(property_ids)
        return E

    def get_G(self, property_ids=None):
        G = self.model.prod.get_G(property_ids)
        return G

    def get_J(self, property_ids=None):
        J = self.model.prod.get_J(property_ids)
        return J

    def get_c(self, property_ids=None):
        c = self.model.prod.get_c(property_ids)
        return c

    def get_mass(self, element_ids=None, total=False):
        """
        mass = rho * A * L + nsm
        """
        if self.n == 0:
            return 0.0

        assert element_ids is None
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

    #=========================================================================
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

    #=========================================================================
    def get_stiffness(self, i, model, positions, index0s, knorm=1.0):
        #print("----------------")
        pid = self.property_id[i]
        assert isinstance(pid, int), pid
        A = self.get_Area(pid)
        E = self.get_E(pid)
        G = self.get_G(pid)
        J = self.get_J(pid)
        print('A=%s E=%s G=%s J=%s' % (A, E, G, J))

        #========================
        #(n1, n2) = self.node_ids()
        n0 = self.node_ids[i, 0]
        n1 = self.node_ids[i, 1]

        i0 = index0s[n0]
        i1 = index0s[n1]

        print("n0", n0)
        print("n1", n1)
        p0 = positions[n0]
        p1 = positions[n1]
        #p1 = model.Node(n1).xyz

        v1 = p0 - p1
        L = norm(v1)
        if L == 0.0:
            msg = 'invalid CROD length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)
        #========================
        print("A=%g E=%g G=%g J=%g L=%g" % (A, E, G, J, L))
        k_axial = A * E / L
        k_torsion = G * J / L
        #k_axial = 1.0
        #k_torsion = 2.0

        k = array([[1., -1.], [-1., 1.]])  # 1D rod

        Lambda = _Lambda(v1, debug=True)
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
        #print("K=\n", K / knorm)
        #print("K2=\n", K2 / knorm)

        #========================

        #print(K / knorm)
        #print("K[%s] = \n%s\n" % (self.eid, list_print(K/knorm)))

        print('dofs =', dofs)
        print('K =\n', list_print(K / knorm))

        return(K2, dofs, nIJV)

    def displacement_stress(self, model, positions, q, dofs):
        n = self.n
        o1 = zeros(n, 'float64')
        e1 = zeros(n, 'float64')
        f1 = zeros(n, 'float64')

        o4 = zeros(n, 'float64')
        e4 = zeros(n, 'float64')
        f4 = zeros(n, 'float64')


        As = self.get_Area(self.property_id)
        Es = self.get_E(self.property_id)
        Gs = self.get_G(self.property_id)
        Js = self.get_J(self.property_id)
        Cs = self.get_c(self.property_id)

        for i in xrange(n):
            A = As[i]
            E = Es[i]
            G = Gs[i]
            E = Es[i]
            J = Js[i]
            C = Cs[i]
            n1, n2 = self.node_ids[i, :]


            p1 = positions[n1]
            p2 = positions[n2]

            v1 = p1 - p2
            L = norm(p1 - p2)
            if L == 0.0:
                msg = 'invalid CROD length=0.0\n%s' % (self.__repr__())
                raise ZeroDivisionError(msg)

            #========================
            #mat = self.get_material_from_index(i)
            #jmat = searchsorted(mat.material_id, self.material_id[i])

            #E = mat.E[jmat]
            #G = mat.G[jmat]
            #G = self.G()

            #print("A=%g E=%g G=%g J=%g L=%g" % (A, E, G, J, L))
            k_axial = A * E / L
            k_torsion = G * J / L
            #k_axial = 1.0
            #k_torsion = 2.0

            #k = array([[1., -1.], [-1., 1.]])  # 1D rod

            Lambda = _Lambda(v1, debug=False)

            #print("**dofs =", dofs)
            n11 = dofs[(n1, 1)]
            n21 = dofs[(n2, 1)]

            n12 = dofs[(n1, 2)]
            n22 = dofs[(n2, 2)]

            n13 = dofs[(n1, 3)]
            n23 = dofs[(n2, 3)]

            # moments
            n14 = dofs[(n1, 4)]
            n24 = dofs[(n2, 4)]

            n15 = dofs[(n1, 5)]
            n25 = dofs[(n2, 5)]

            n16 = dofs[(n1, 6)]
            n26 = dofs[(n2, 6)]

            q_axial = array([
                q[n11], q[n12], q[n13],
                q[n21], q[n22], q[n23]
            ])
            q_torsion = array([
                q[n14], q[n15], q[n16],
                q[n24], q[n25], q[n26]
            ])
            #print("type=%s n1=%s n2=%s" % (self.type, n1, n2))
            #print("n11=%s n12=%s n21=%s n22=%s" %(n11,n12,n21,n22))

            #print("q2[%s] = %s" % (self.eid, q2))
            #print("Lambda = \n"+str(Lambda))

            #print "Lsize = ",Lambda.shape
            #print "qsize = ",q.shape
            u_axial = dot(array(Lambda), q_axial)
            du_axial = -u_axial[0] + u_axial[1]
            u_torsion = dot(array(Lambda), q_torsion)
            du_torsion = -u_torsion[0] + u_torsion[1]

            #L = self.Length()
            #E = self.E()
            #A = self.Area()

            #C = self.C()
            #J = self.J()
            #G = self.G()

            axial_strain = du_axial / L
            torsional_strain = du_torsion * C / L

            axial_stress = E * axial_strain
            torsional_stress = G * torsional_strain

            axial_force = axial_stress * A
            torsional_moment = du_torsion * G * J / L
            #print("axial_strain = %s [psi]" % axial_strain)
            #print("axial_stress = %s [psi]" % axial_stress)
            #print("axial_force  = %s [lb]\n" % axial_force)
            o1[i] = axial_stress
            o4[i] = torsional_stress

            e1[i] = axial_strain
            e4[i] = torsional_strain

            f1[i] = axial_force
            f4[i] = torsional_moment

        return (e1, e4,
                o1, o4,
                f1, f4)