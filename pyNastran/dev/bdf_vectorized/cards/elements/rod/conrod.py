#from numpy import arange, zeros, unique, dot, array, transpose
import numpy as np
from numpy.linalg import norm  # type: ignore

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import integer, double, double_or_blank

from pyNastran.dev.bdf_vectorized.cards.elements.rod.rod_element import RodElement

def _Lambda(v1, debug=True):
    """
    ::
      3d  [l,m,n,0,0,0]  2x6
          [0,0,0,l,m,n]
    """
    #R = self.Rmatrix(model,is3D)

    #xyz1 = model.Node(n1).get_position()
    #xyz2 = model.Node(n2).get_position()
    #v1 = xyz2 - xyz1
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
    Lambda = np.zeros((2, 6), 'd')
    Lambda[0, 0] = Lambda[1, 3] = l
    Lambda[0, 1] = Lambda[1, 4] = m
    Lambda[0, 2] = Lambda[1, 5] = n

    #print("R = \n",R)
    #debug = True
    if debug:
        print("Lambda = \n" + str(Lambda))
    return Lambda


class CONROD(RodElement):
    type = 'CONROD'
    def __init__(self, model):
        """
        Defines the CONROD object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        RodElement.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            self.element_id = np.zeros(ncards, 'int32')
            self.material_id = np.zeros(ncards, 'int32')
            self.node_ids = np.zeros((ncards, 2), 'int32')
            self.A = np.zeros(ncards, float_fmt)
            self.J = np.zeros(ncards, float_fmt)
            self.c = np.zeros(ncards, float_fmt)
            self.nsm = np.zeros(ncards, float_fmt)

    def add_card(self, card, comment=''):
        self.model.log.debug('  adding CONROD')
        i = self.i
        eid = integer(card, 1, 'element_id')
        self.element_id[i] = eid
        self.node_ids[i] = [integer(card, 2, 'node_1'),
                            integer(card, 3, 'node_2')]

        self.material_id[i] = integer(card, 4, 'material_id')
        self.A[i] = double(card, 5, 'Area')
        self.J[i] = double_or_blank(card, 6, 'J', 0.0)
        self.c[i] = double_or_blank(card, 7, 'c', 0.0)
        self.nsm[i] = double_or_blank(card, 8, 'non_structural_mass', 0.0)
        assert len(card) <= 9, 'len(CONROD card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        """
        Parameters
        ----------
        cards : list[BDFCard(), ...]
            the list of CONROD cards
        """
        assert self.n != 0, self.n
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.node_ids = self.node_ids[i, :]
            self.material_id = self.material_id[i]
            self.A = self.A[i]
            self.J = self.J[i]
            self.c = self.c[i]
            self.nsm = self.nsm[i]

            unique_eids = np.unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CONROD IDs...')
        else:
            self.element_id = np.array([], dtype='int32')
            self.property_id = np.array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'property' : pid_map,
            'material' : mid_map,
        }
        """
        if self.n:
            nid_map = maps['node']
            mid_map = maps['material']
            for i, (nids, mids) in enumerate(zip(self.node_id, self.material_ids)):
                self.node_ids[i, 0] = nid_map[nids[0]]
                self.node_ids[i, 1] = nid_map[nids[1]]
                self.material_id[i] = mid_map[pid]

    #=========================================================================
    def get_area_from_index(self, i):
        return self.A[i]

    #def get_E_from_index(self, i):
        #return self.a[i]

    def get_length_by_element_id(self, element_id=None, grid_cid0=None):
        i = self.get_element_index_by_element_id(element_id)
        return self.get_length_by_element_index(i, grid_cid0)

    def get_length_by_element_index(self, i=None, grid_cid0=None):
        if i is None:
            i = np.arange(self.n)

        if grid_cid0 is None:
            grid_cid0 = self.model.grid.get_position_by_node_id()
        assert grid_cid0 is not None
        #if positions is None:

        self.model.log.debug("i = %s" % i)
        n = len(i)
        if n == 0:
            i = i[0]
        n1 = self.node_ids[i, 0]
        n2 = self.node_ids[i, 1]
        n1i = self.model.grid.get_node_index_by_node_id(n1)
        n2i = self.model.grid.get_node_index_by_node_id(n2)

        #n1, n2 = self.node_ids[i, :]
        self.model.log.debug('grids\n%s' % grid_cid0)
        p1 = grid_cid0[n1i, :]
        p2 = grid_cid0[n2i, :]
        v1 = p1 - p2
        L = norm(p1 - p2, axis=1)
        assert L.shape == (n,), L.shape
        return L

    def get_material_id_by_element_id(self, element_id=None):
        i = self.get_element_index_by_element_id(element_id)
        return self.get_material_id_by_element_index(i)

    def get_material_id_by_element_index(self, i=None):
        return self.material_id[i]

    def get_material_from_index(self, i):
        return self.model.materials.mat1

    def get_area_by_element_id(self, element_id=None):
        i = self.get_element_index_by_element_id(element_id, msg='')
        return self.get_area_by_element_index(i)

    def get_area_by_element_index(self, i=None):
        A = self.A[i]
        return A

    def get_E_by_element_id(self, element_id=None):
        i = self.get_element_index_by_element_id(element_id, msg='')
        mat = self.model.materials.mat1.slice_by_material_id(self.material_id[i])
        E = mat.E()
        G = mat.G()
        return E

    def get_G_by_element_id(self, element_id=None):
        i = self.get_element_index_by_element_id(element_id, msg='')
        mat = self.model.materials.mat1.slice_by_material_id(self.material_id[i])
        E = mat.E()
        G = mat.G()
        return G

    def get_J_by_element_id(self, element_id=None):
        J = self.model.prod.get_J_by_property_id(property_id)
        return J

    def get_c_by_element_id(self, element_id=None):
        c = self.model.prod.get_c_by_property_id(property_id)
        return c

    def get_node_indicies_by_element_index(self, i=None):
        if i is None:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 1])
        else:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 1])
        return i1, i2

    def get_density_by_element_id(self, element_id=None):
        i = self.get_element_index_by_element_id(element_id)
        return self.get_density_by_element_index(i)

    def get_density_by_element_index(self, i=None):
        j = self.material_id[i]
        rho = self.model.materials.get_density_by_material_id(j)
        return rho

    def get_non_structural_mass_by_element_id(self, element_id=None):
        i = self.get_element_index_by_element_id(element_id)
        return self.get_non_structural_mass_by_element_index(i)

    def get_non_structural_mass_by_element_index(self, i=None):
        return self.nsm[i]

    #def get_mass_by_element_index(self, total=False, i=None):
    def get_mass_by_element_index(self, total=False, i=None):
        """
        mass = rho * A * L + nsm
        """
        if self.n == 0:
            return 0.0

        #i = self.model.grid.get_node_index_by_node_id()
        #grid_cid0 = self.model.grid.get_positions_by_index(i)


        # doesn't support indicies
        #self.model.log.debug('grid_cid0 = %s' % grid_cid0)
        #self.model.log.debug('node0 = %s' % self.node_ids[:, 0])
        #p1 = grid_cid0[self.node_ids[:, 0]]
        #p2 = grid_cid0[self.node_ids[:, 1]]


        #grid_cid0 = self.model.grid.get_positions_by_index(self.model.grid.get_node_index_by_node_id())
        try:
            msg = ', which is required by CONROD get_mass'
            i1, i2 = self.get_node_indicies_by_element_index(i)
            #n1 = xyz_cid0[i1, :]
            #n2 = xyz_cid0[i2, :]
            p1 = self.model.grid.get_position_by_node_index(i1)
            p2 = self.model.grid.get_position_by_node_index(i2)
        except RuntimeError:
            for eid, (n1, n2) in zip(self.element_id, self.node_ids):
                msg = ', which is required by CONROD element_id=%s node1=%s\n' % (eid, n1)
                i1 = self.model.grid.get_node_index_by_node_id([n1], msg=msg)
                msg = ', which is required by CONROD element_id=%s node2=%s\n' % (eid, n2)
                i2 = self.model.grid.get_node_index_by_node_id([n2], msg=msg)
                p1 = self.model.grid.get_position_by_node_index(i1)
                p2 = self.model.grid.get_position_by_node_index(i2)

        L = norm(p2 - p1, axis=1)
        rho = self.model.materials.get_density_by_material_id(self.material_id[i])
        #rho = self.get_element_id_by_element_index(self.material_id)

        self.model.log.debug('*L = %s' % L)
        self.model.log.debug('*rho = %s' % rho)
        self.model.log.debug('*A = %s' % self.A[i])
        self.model.log.debug('*nsm = %s' % self.nsm[i])
        mass = L * (rho * self.A[i] + self.nsm[i])
        self.model.log.debug('*mass = %s' % mass)
        if total:
            return mass.sum()
        else:
            return mass

    def write_card(self, bdf_file, size=8, element_id=None):
        if self.n:
            if element_id is None:
                i = np.arange(self.n)
            else:
                i = np.searchsorted(self.element_id, element_id)
            if isinstance(i, (int, np.int64)):
                i = np.array([i])

            for (eid, n12, mid, A, J, c, nsm) in zip(
                 self.element_id, self.node_ids, self.material_id, self.A, self.J,
                 self.c, self.nsm):

                card = ['CONROD', eid, n12[0], n12[1], mid, A, J, c, nsm]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    #=========================================================================

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
        rho = mat.get_density_by_material_id()

        n1, n2 = self.node_ids[i, :]

        i1 = index0s[n1]
        i2 = index0s[n2]

        xyz1 = positions[n1]
        xyz2 = positions[n2]
        v1 = xyz1 - xyz2
        L = norm(v1)
        if L == 0.0:
            msg = 'invalid CONROD length=0.0\n%s' % self.__repr__()
            raise ZeroDivisionError(msg)
        #========================
        mi = (rho * A * L + self.nsm[i]) / 6.
        m = np.array([[2., 1.],
                      [1., 2.]])  # 1D rod

        Lambda = _Lambda(v1, debug=False)
        M = (Lambda.T @ m) @ Lambda

        Mi, Mj = M.shape
        dofs = np.array([
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,
        ], 'int32')
        n_ijv = [
            # axial
            (n1, 1), (n1, 2), (n1, 3),
            (n2, 1), (n2, 2), (n2, 3),

            # torsion -> NA
        ]
        self.model.log.info('dofs = %s' % dofs)
        return(M, dofs, n_ijv)

    def get_stiffness_matrix(self, i, model, positions, index0s, knorm=1.0):  # CROD/CONROD
        #print("----------------")
        A = self.get_area_from_index(i)
        #mat = self.get_material_from_index(i)
        #print(mat)
        #print(mat.material_id[se])
        mid = self.material_id[i]
        #self.model.log.info('mid = %s' % mid)
        #print(self.model.materials.mat1.print_card())
        imid = self.model.materials.mat1.get_material_index_by_material_id(mid)
        mat = self.model.materials.mat1[imid]
        E = mat.E[0]
        G = mat.G[0]
        #G = self.G()
        J = self.J[i]
        #J = self.J()

        #========================
        #(n1, n2) = self.node_ids
        n1, n2 = self.node_ids[i, :]

        i1 = index0s[n1]
        i2 = index0s[n2]

        xyz1 = positions[n1]
        xyz2 = positions[n2]
        #p1 = model.Node(n1).xyz

        dxyz12 = xyz1 - xyz2
        L = norm(dxyz12)
        if L == 0.0:
            msg = 'invalid CONROD length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)
        #========================
        #print("A=%r E=%r G=%r J=%r L=%r" % (A, E, G, J, L))
        k_axial = A * E / L
        k_torsion = G * J / L
        #k_axial = 1.0
        #k_torsion = 2.0

        k = np.array([[1., -1.],
                      [-1., 1.]])  # 1D rod

        Lambda = _Lambda(dxyz12, debug=False)
        K = (Lambda.T @ k) @ Lambda
        Ki, Kj = K.shape

        K2 = np.zeros((Ki*2, Kj*2), 'float64')
        if k_axial == 0.0 and k_torsion == 0.0:
            dofs = []
            n_ijv = []
            K2 = []
        elif k_torsion == 0.0: # axial; 2D or 3D
            K2 = K * k_axial
            dofs = np.array([
                i1, i1+1, i1+2,
                i2, i2+1, i2+2,
            ], 'int32')
            n_ijv = [
                # axial
                (n1, 1), (n1, 2), (n1, 3),
                (n2, 1), (n2, 2), (n2, 3),
            ]
        elif k_axial == 0.0: # torsion; assume 3D
            K2 = K * k_torsion
            dofs = np.array([
                i1+3, i1+4, i1+5,
                i2+3, i2+4, i2+5,
            ], 'int32')
            n_ijv = [
                # torsion
                (n1, 4), (n1, 5), (n1, 6),
                (n2, 4), (n2, 5), (n2, 6),
            ]

        else:  # axial + torsion; assume 3D
            # u1fx, u1fy, u1fz, u2fx, u2fy, u2fz
            K2[:Ki, :Ki] = K * k_axial

            # u1mx, u1my, u1mz, u2mx, u2my, u2mz
            K2[Ki:, Ki:] = K * k_torsion

            dofs = np.array([
                i1, i1+1, i1+2,
                i2, i2+1, i2+2,

                i1+3, i1+4, i1+5,
                i2+3, i2+4, i2+5,
            ], 'int32')
            n_ijv = [
                # axial
                (n1, 1), (n1, 2), (n1, 3),
                (n2, 1), (n2, 2), (n2, 3),

                # torsion
                (n1, 4), (n1, 5), (n1, 6),
                (n2, 4), (n2, 5), (n2, 6),
            ]

        #Fg = (Lambda.T @ grav) @ Lambda
        #print("K=\n", K / knorm)
        #print("K2=\n", K2 / knorm)

        #========================

        #print(K / knorm)
        #print("K[%s] = \n%s\n" % (self.eid, list_print(K/knorm)))

        self.model.log.info('dofs = %s' % dofs)
        self.model.log.debug('K =\n\n%s' % (K / knorm))

        return(K2, dofs, n_ijv)

    def displacement_stress(self, model, positions, q, dofs):
        n = self.n
        o1 = np.zeros(n, 'float64')
        e1 = np.zeros(n, 'float64')
        f1 = np.zeros(n, 'float64')

        o4 = np.zeros(n, 'float64')
        e4 = np.zeros(n, 'float64')
        f4 = np.zeros(n, 'float64')

        for i in range(n):
            J = self.J[i]
            C = self.c[i]
            n1, n2 = self.node_ids[i, :]


            xyz1 = positions[n1]
            xyz2 = positions[n2]

            v1 = xyz1 - xyz2
            L = norm(xyz1 - xyz2)
            if L == 0.0:
                msg = 'invalid CONROD length=0.0\n%s' % (self.__repr__())
                raise ZeroDivisionError(msg)

            A = self.get_area_from_index(i)
            mid = self.material_id[i]
            #mat = self.model.materials.mat1[mid]
            imid = self.model.materials.mat1.get_material_index_by_material_id(mid)
            mat = self.model.materials.mat1[imid]

            E = mat.E
            G = mat.G
            #========================
            #mat = self.get_material_from_index(i)
            #jmat = searchsorted(mat.material_id, self.material_id[i])

            #E = mat.E[jmat]
            #G = mat.G[jmat]
            #G = self.G()

            #print("A=%g E=%g G=%g J=%g L=%g" % (A, E, G, J, L))
            k_axial = A * E / L
            k_torsion = G * J / L

            #k = array([[1., -1.],
                       #[-1., 1.]])  # 1D rod

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

            q_axial = np.array([
                q[n11], q[n12], q[n13],
                q[n21], q[n22], q[n23]
            ])
            q_torsion = np.array([
                q[n14], q[n15], q[n16],
                q[n24], q[n25], q[n26]
            ])
            #print("type=%s n1=%s n2=%s" % (self.type, n1, n2))
            #print("n11=%s n12=%s n21=%s n22=%s" %(n11,n12,n21,n22))

            #print("q2[%s] = %s" % (self.eid, q2))
            #print("Lambda = \n"+str(Lambda))

            #print("Lsize = ",Lambda.shape)
            #print("qsize = ",q.shape)
            u_axial = np.array(Lambda) @ q_axial
            du_axial = u_axial[0] - u_axial[1]
            u_torsion = np.array(Lambda) @ q_torsion
            du_torsion = u_torsion[0] - u_torsion[1]

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

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        obj = CONROD(self.model)
        n = len(i)
        obj.n = n
        obj.i = n
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.node_ids = self.node_ids[i, :]
        obj.material_id = self.material_id[i]
        obj.A = self.A[i]
        obj.J = self.J[i]
        obj.c = self.c[i]
        obj.nsm = self.nsm[i]
        return obj
