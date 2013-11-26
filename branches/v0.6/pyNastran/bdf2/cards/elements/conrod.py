from numpy import zeros, unique
from numpy.linalg import norm
from numpy import dot, searchsorted, array, transpose

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string)

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
            self.element_id = zeros(ncards, 'int32')
            self.material_id = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 2), 'int32')
            self.A = zeros(ncards, float_fmt)
            self.J = zeros(ncards, float_fmt)
            self.c = zeros(ncards, float_fmt)
            self.nsm = zeros(ncards, float_fmt)

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'eid')
                self.node_ids[i] = [integer(card, 2, 'n1'),
                                    integer(card, 3, 'n2')]

                self.material_id[i] = integer(card, 4, 'mid')
                self.A[i] = double(card, 5, 'A')
                self.J[i] = double_or_blank(card, 6, 'j', 0.0)
                self.c[i] = double_or_blank(card, 7, 'c', 0.0)
                self.nsm[i] = double_or_blank(card, 8, 'nsm', 0.0)
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
                raise RuntimeError('There are duplicate CROD IDs...')

    def _displacement_stress(self, q):
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
        E = self.Materials.get_E(self.material_id)
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

    def get_area_from_index(self, i):
        return self.A[i]

    #def get_E_from_index(self, i):
        #return self.a[i]

    def get_material_from_index(self, i):
        return self.model.materials.mat1

    def get_stiffness(self, i, model, positions, index0s, fnorm=1.0):  # CROD/CONROD
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
        L = norm(p0 - p1)
        if L == 0.0:
            msg = 'invalid CROD length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)
        #========================
        print("A=%r E=%r G=%r J=%r L=%r" % (A, E, G, J, L))
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
        #print("K=\n", K / fnorm)
        #print("K2=\n", K2 / fnorm)

        #========================

        #print(K / fnorm)
        #print("K[%s] = \n%s\n" % (self.eid, list_print(K/fnorm)))

        print('dofs =', dofs)
        #print('K =\n', list_print(K / fnorm))

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
            #j = searchsorted(mat.material_id, self.material_id[i])

            #E = mat.E[j]
            G = mat.G[j]
            #G = self.G()
            J = self.J[i]
            C = self.c[i]

            print("A=%g E=%g G=%g J=%g L=%g" % (A, E, G, J, L))
            k_axial = A * E / L
            k_torsion = G * J / L
            #k_axial = 1.0
            #k_torsion = 2.0

            #k = array([[1., -1.], [-1., 1.]])  # 1D rod

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


