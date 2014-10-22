from itertools import izip

from numpy import arange, zeros, unique
from numpy.linalg import norm
from numpy import dot, searchsorted, array, transpose

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string)

from pyNastran.bdf.dev_vectorized.cards.elements.rod.rod_element import RodElement

def _Lambda(v1, debug=True):
    """
    ::
      3d  [l,m,n,0,0,0]  2x6
          [0,0,0,l,m,n]
    """
    #R = self.Rmatrix(model,is3D)

    #p1 = model.Node(n1).Position()
    #p2 = model.Node(n2).Position()
    #v1 = p2 - p1
    if debug:
        print("v1=%s" % v1)
    n = norm(v1)
    if n == 0:
        raise ZeroDivisionError(v1)
    v1 = v1 / n
    (l, m, n) = v1
    #l = 1
    #m = 2
    #n = 3
    Lambda = zeros((2, 6), 'd')
    Lambda[0, 0] = Lambda[1, 3] = l
    Lambda[0, 1] = Lambda[1, 4] = m
    Lambda[0, 2] = Lambda[1, 5] = n

    #print("R = \n",R)
    #debug = True
    if debug:
        print("Lambda = \n" + str(Lambda))
        #sys.exit('asdf')
    return Lambda


class CONROD(RodElement):
    type = 'CONROD'
    def __init__(self, model):
        """
        Defines the CONROD object.

        :param self: the CONROD object
        :param model: the BDF object
        """
        RodElement.__init__(self, model)

    def allocate(self, ncards):
        float_fmt = self.model.float
        self.element_id = zeros(ncards, 'int32')
        self.material_id = zeros(ncards, 'int32')
        self.node_ids = zeros((ncards, 2), 'int32')
        self.A = zeros(ncards, float_fmt)
        self.J = zeros(ncards, float_fmt)
        self.c = zeros(ncards, float_fmt)
        self.nsm = zeros(ncards, float_fmt)

    def build(self):
        """
        :param self: the CONROD object
        :param cards: the list of CONROD cards
        """
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            #: Property ID
            self.element_id = zeros(ncards, 'int32')
            self.material_id = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 2), 'int32')
            self.A = zeros(ncards, float_fmt)
            self.J = zeros(ncards, float_fmt)
            self.c = zeros(ncards, float_fmt)
            self.nsm = zeros(ncards, float_fmt)

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'element_id')
                self.node_ids[i] = [integer(card, 2, 'node_1'),
                                    integer(card, 3, 'node_2')]

                self.material_id[i] = integer(card, 4, 'material_id')
                self.A[i] = double(card, 5, 'Area')
                self.J[i] = double_or_blank(card, 6, 'J', 0.0)
                self.c[i] = double_or_blank(card, 7, 'c', 0.0)
                self.nsm[i] = double_or_blank(card, 8, 'non_structural_mass', 0.0)
                assert len(card) <= 9, 'len(CONROD card) = %i' % len(card)

            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.node_ids = self.node_ids[i, :]
            self.material_id = self.material_id[i]
            self.A = self.A[i]
            self.J = self.J[i]
            self.c = self.c[i]
            self.nsm = self.nsm[i]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CONROD IDs...')
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    #=========================================================================
    def get_area_from_index(self, i):
        return self.A[i]

    #def get_E_from_index(self, i):
        #return self.a[i]

    def get_material_from_index(self, i):
        return self.model.materials.mat1

    def get_area(self, property_ids=None):
        A = self.A
        return A

    def get_E(self, property_ids=None):
        mat = self.model.materials.mat1[self.material_id[i]]
        E = mat.E()
        G = mat.G()
        return E

    def get_G(self, property_ids=None):
        mat = self.model.materials.mat1[self.material_id[i]]
        E = mat.E()
        G = mat.G()
        return G

    def get_J(self, property_ids=None):
        J = self.model.prod.get_J(property_ids)
        return J

    def get_c(self, property_ids=None):
        c = self.model.prod.get_c(property_ids)
        return c

    def get_mass(self, total=False):
        """
        mass = rho * A * L + nsm
        """
        if self.n == 0:
            return 0.0

        grid_cid0 = self.model.grid.get_positions()
        p1 = grid_cid0[self.node_ids[:, 0]]
        p2 = grid_cid0[self.node_ids[:, 1]]
        L = p2 - p1
        rho = self.model.materials.get_density(self.material_id)
        mass = norm(L, axis=1) * self.A * rho + self.nsm
        if total:
            return mass.sum()
        else:
            return mass

    #=========================================================================

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            for (eid, n12, mid, A, J, c, nsm) in izip(
                 self.element_id, self.node_ids, self.material_id, self.A, self.J,
                 self.c, self.nsm):

                card = ['CONROD', eid, n12[0], n12[1], mid, A, J, c, nsm ]
                f.write(print_card(card))

    def get_mass_matrix(self, i, model, positions, index0s, knorm=1.0):  # CROD/CONROD
        """
        Lumped:
        =======
          mi = 1/2 * rho * A * L
                 [ 1  0 ]
          M = mi [ 0  1 ]

        Consistent:
        ===========
          mi = 1/6 * rho * A * L
                 [ 2 1 ]
          M = mi [ 1 2 ]
        """
        A = self.get_area_from_index(i)
        mat = self.model.materials.mat1[self.material_id[i]]
        mid = self.material_id[i]
        rho = mat.get_density()
        #========================
        n0, n1 = self.node_ids[i, :]

        i0 = index0s[n0]
        i1 = index0s[n1]

        p0 = positions[n0]
        p1 = positions[n1]
        v1 = p0 - p1
        L = norm(v1)
        if L == 0.0:
            msg = 'invalid CROD length=0.0\n%s' % self.__repr__()
            raise ZeroDivisionError(msg)
        #========================
        mi = (rho * A * L + self.nsm[i]) / 6.
        m = array([[2., 1.],
                   [1., 2.]])  # 1D rod

        Lambda = _Lambda(v1, debug=False)
        M = dot(dot(transpose(Lambda), m), Lambda)
        Mi, Mj = M.shape
        dofs = array([
            i0, i0+1, i0+2,
            i1, i1+1, i1+2,
        ], 'int32')
        nIJV = [
            # axial
            (n0, 1), (n0, 2), (n0, 3),
            (n1, 1), (n1, 2), (n1, 3),

            # torsion -> NA
        ]

        print('dofs =', dofs)
        return(M, dofs, nIJV)

    def get_stiffness_matrix(self, i, model, positions, index0s, knorm=1.0):  # CROD/CONROD
        #print("----------------")
        A = self.get_area_from_index(i)
        #mat = self.get_material_from_index(i)
        #print mat
        #print mat.material_id[se]
        mat = self.model.materials.mat1[self.material_id[i]]
        E = mat.E()
        G = mat.G()
        #G = self.G()
        J = self.J[i]
        #J = self.J()

        #========================
        #(n1, n2) = self.nodeIDs()
        n0, n1 = self.node_ids[i, :]

        i0 = index0s[n0]
        i1 = index0s[n1]
        #node0 = self.nodes[0]
        #node1 = self.nodes[1]


        #p0 = model.Node(n0).xyz
        #p1 = model.Node(n1).xyz
        p0 = positions[n0]
        p1 = positions[n1]

        v1 = p0 - p1
        L = norm(v1)
        if L == 0.0:
            msg = 'invalid CROD length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)
        #========================
        #print("A=%r E=%r G=%r J=%r L=%r" % (A, E, G, J, L))
        k_axial = A * E / L
        k_torsion = G * J / L
        #k_axial = 1.0
        #k_torsion = 2.0

        k = array([[1., -1.],
                   [-1., 1.]])  # 1D rod

        Lambda = _Lambda(v1, debug=False)
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
        #print('K =\n', list_print(K / knorm))

        return(K2, dofs, nIJV)

    def displacement_stress(self, model, positions, q, dofs):
        n = self.n
        o1 = zeros(n, 'float64')
        e1 = zeros(n, 'float64')
        f1 = zeros(n, 'float64')

        o4 = zeros(n, 'float64')
        e4 = zeros(n, 'float64')
        f4 = zeros(n, 'float64')

        for i in xrange(n):
            n1, n2 = self.node_ids[i, :]


            p1 = positions[n1]
            p2 = positions[n2]

            v1 = p1 - p2
            L = norm(p1 - p2)
            if L == 0.0:
                msg = 'invalid CROD length=0.0\n%s' % (self.__repr__())
                raise ZeroDivisionError(msg)
            #========================
            A = self.get_area_from_index(i)

            mat = self.model.materials.mat1[self.material_id[i]]
            E = mat.E()
            G = mat.G()
            #mat = self.get_material_from_index(i)
            #jmat = searchsorted(mat.material_id, self.material_id[i])

            #E = mat.E[jmat]
            #G = mat.G[jmat]
            #G = self.G()
            J = self.J[i]
            C = self.c[i]

            print("A=%g E=%g G=%g J=%g L=%g" % (A, E, G, J, L))
            k_axial = A * E / L
            k_torsion = G * J / L

            #k = array([[1., -1.],
                       #[-1., 1.]])  # 1D rod

            Lambda = _Lambda(v1, debug=True)

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