from numpy import array, arange, zeros, unique, searchsorted, int64
from numpy.linalg import norm  # type: ignore

from pyNastran.dev.bdf_vectorized.cards.elements.rod.conrod import _Lambda
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.utils.dev import list_print

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import integer, integer_or_blank
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

from pyNastran.dev.bdf_vectorized.cards.elements.rod.rod_element import RodElement

class CTUBE(RodElement):
    type = 'CTUBE'
    def __init__(self, model):
        """
        Defines the CTUBE object.

        Parameters
        ----------
        model : BDF
            the BDF object

        """
        RodElement.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.model.log.debug('  nCTUBE=%s' % ncards)
            self.n = ncards
            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            #: Node IDs
            self.node_ids = zeros((ncards, 2), 'int32')

    def add_card(self, card: BDFCard, comment: str=''):
        self.model.log.debug('  adding CTUBE')
        i = self.i
        eid = integer(card, 1, 'element_id')
        #if comment:
            #self.set_comment(eid, comment)


        self.element_id[i] = integer(card, 1, 'element_id')
        self.property_id[i] = integer_or_blank(card, 2, 'property_id', self.element_id[i])
        self.node_ids[i] = [integer(card, 3, 'n1'),
                            integer(card, 4, 'n2')]
        assert len(card) == 5, 'len(CTUBE card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CTUBE IDs...')
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'property' : pid_map,
            'material' : mid_map,
        }
        """
        if self.n:
            nid_map = maps['node']
            pid_map = maps['property']
            for i, (nids, pid) in enumerate(zip(self.node_id, self.property_id)):
                self.node_ids[i, 0] = nid_map[nids[0]]
                self.node_ids[i, 1] = nid_map[nids[1]]
                self.property_id[i] = pid_map[pid]

    #=========================================================================
    def get_property_id_by_element_index(self, i=None):
        return self.property_id[i]

    def get_property_id_by_element_id(self, element_id):
        i = self.get_element_index_by_element_id(element_id)
        pid = self.get_property_id_by_element_index(i)
        return pid

    def get_E_by_element_id(self, element_id):
        pid = self.get_property_id_by_element_id(element_id)
        E = self.get_E_by_property_id(pid)
        return E

    def get_G_by_element_id(self, element_id=None):
        pid = self.get_property_id_by_element_id(element_id)
        G = self.get_G_by_property_id(pid)
        return G

    def get_J_by_element_id(self, element_id=None):
        pid = self.get_property_id_by_element_id(element_id)
        J = self.get_J_by_property_id(pid)
        return J

    def get_c_by_element_id(self, element_id=None):
        pid = self.get_property_id_by_element_id(element_id)
        c = self.get_c_by_property_id(pid)
        return c

    def get_area_by_element_id(self, element_id=None):
        i = self.get_element_index_by_element_id(element_id)
        return self.get_area_by_element_index(i)

    def get_area_by_element_index(self, i=None):
        pid = self.get_property_id_by_element_index(i)
        area = self.get_area_by_property_id(pid)
        return area

    def get_non_structural_mass_by_element_id(self, element_id=None):
        pid = self.get_property_id_by_element_id(element_id)
        nsm = self.get_non_structural_mass_by_property_id(pid)
        return nsm

    def get_material_id_by_element_id(self, element_id=None):
        pid = self.get_property_id_by_element_id(element_id)
        mid = self.get_material_id_by_property_id(pid)
        return mid

    def get_density_by_element_id(self, element_id=None):
        pid = self.get_property_id_by_element_id(element_id)
        density = self.get_density_by_property_id(pid)
        return density

    #=========================================================================

    def get_area_by_property_id(self, property_id=None):
        A = self.model.ptube.get_area_by_property_id(property_id)
        return A

    def get_E_by_property_id(self, property_id=None):
        E = self.model.ptube.get_E_by_property_id(property_id)
        return E

    def get_G_by_property_id(self, property_id=None):
        G = self.model.ptube.get_G_by_property_id(property_id)
        return G

    def get_J_by_property_id(self, property_id=None):
        J = self.model.ptube.get_J_by_property_id(property_id)
        return J

    def get_c_by_property_id(self, property_ids=None):
        c = self.model.ptube.get_c(property_ids)
        return c

    def get_non_structural_mass_by_property_id(self, property_id=None):
        c = self.model.ptube.get_non_structural_mass_by_property_id(property_id)
        return c

    def get_density_by_property_id(self, property_id=None):
        density = self.model.ptube.get_density_by_property_id(property_id)
        return density

    #=========================================================================
    def get_material_id_by_element_index(self, i=None):
        pid = self.get_property_id_by_element_index(i)
        mid = self.model.ptube.get_material_id_by_property_id(pid)
        return mid

    def get_material_id_by_property_id(self, property_id=None):
        mid = self.model.ptube.get_material_id_by_property_id(property_id)
        return mid

    def get_length_by_element_id(self, element_id=None, grid_cid0=None):
        i = self.get_element_index_by_element_id(element_id=None, msg='')
        L = self.get_length_by_element_index(i, grid_cid0)
        return L

    def get_length_by_element_index(self, i=None, grid_cid0=None):
        #if i is None:
            #i = arange(self.n)

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
    #=========================================================================

    def get_node_indicies(self, i=None):
        if i is None:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 1])
        else:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 1])
        return i1, i2

    def _node_locations(self, xyz_cid0, i=None):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_node_index()
        i1, i2 = self.get_node_indicies(i)
        n1 = xyz_cid0[i1, :]
        n2 = xyz_cid0[i2, :]
        return n1, n2

    def get_mass_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        if isinstance(element_id, integer_types):
            assert element_id > 0, element_id
        elif element_id is None:
            pass
        elif isinstance(element_id, list):
            assert min(element_id) > 0, element_id
        else:
            assert element_id.min() > 0, element_id

        i = self.get_element_index_by_element_id(element_id)
        return self.get_mass_by_element_index(i, xyz_cid0=None, total=False)

    def get_mass_by_element_index(self, i=None, xyz_cid0=None, total=False):
        """
        mass = rho * A * L + nsm
        """
        if self.n == 0:
            return 0.0
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_node_index()

        #assert i is None, i
        n1, n2 = self._node_locations(xyz_cid0, i)
        L = norm(n2 - n1, axis=1)
        i = self.model.ptube.get_property_index_by_property_id(self.property_id)
        A = self.model.ptube.get_area_by_property_index(i)
        mid = self.model.ptube.material_id[i]
        J = self.model.ptube.get_J_by_property_id(self.property_id)

        #rho, E = self.model.materials.get_density_E(mid)
        rho = self.model.materials.get_density_by_material_id(mid)
        #E = self.model.materials.get_E(mid)
        #return 0. if total else [0.]

        nsm = self.get_non_structural_mass_by_property_id()
        mass = L * (A * rho + nsm)
        if total:
            return mass.sum()
        else:
            return mass

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
        i = self.model.ptube.get_index(self.property_id)
        A = self.model.ptube.A[i]
        mid = self.model.ptube.material_id[i]
        rho = self.model.materials.get_density_by_material_id(mid)
        #========================
        xyz_cid0 = None
        #xyz1, xyz2 = self._node_locations(xyz_cid0)
        if self.n == 1:
            n1, n2 = self.node_ids[0, :]
        else:
            n1, n2 = self.node_ids[i, :]

        i1 = index0s[n1]
        i2 = index0s[n2]

        p1 = positions[n1]
        p2 = positions[n2]
        v1 = p1 - p2
        L = norm(v1)
        if L == 0.0:
            msg = 'invalid CTUBE length=0.0\n%s' % self.__repr__()
            raise ZeroDivisionError(msg)
        #========================
        nsm = self.get_non_structural_mass(self.property_id[i])
        mi = (rho * A * L + nsm) / 6.
        m = array([[2., 1.],
                   [1., 2.]])  # 1D rod

        Lambda = _Lambda(v1, debug=False)
        M = (Lambda.T @ m) @ Lambda

        Mi, Mj = M.shape
        dofs = array([
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
        return M, dofs, n_ijv

    #=========================================================================
    def write_card(self, bdf_file, size=8, is_double=False, element_id=None):
        if self.n:
            #print('eid %s' % element_id)
            if element_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_id)
            if isinstance(i, (int, int64)):
                i = array([i])

            #print('i =', i, type(i))
            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i]):
                card = ['CTUBE', eid, pid, n[0], n[1]]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    #=========================================================================
    def get_stiffness_matrix(self, i, model, positions, index0s, knorm=1.0):
        #print("----------------")
        pid = self.property_id[i]
        assert isinstance(pid, int), pid
        A = self.get_area_by_element_id(pid)
        E = self.get_E_by_element_id(pid)
        G = self.get_G_by_element_id(pid)
        J = self.get_J_by_element_id(pid)
        #print('A=%s E=%s G=%s J=%s' % (A, E, G, J))

        #========================
        #(n1, n2) = self.node_ids()
        n1 = self.node_ids[i, 0]
        n2 = self.node_ids[i, 1]

        i1 = index0s[n1]
        i2 = index0s[n2]

        #print("n0", n0)
        #print("n1", n1)
        n1 = positions[n1]
        n2 = positions[n2]
        #p1 = model.Node(n1).xyz

        v1 = n1 - n2
        L = norm(v1)
        if L == 0.0:
            msg = 'invalid CROD length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)
        #========================
        #print("A=%g E=%g G=%g J=%g L=%g" % (A, E, G, J, L))
        k_axial = A * E / L
        k_torsion = G * J / L
        #k_axial = 1.0
        #k_torsion = 2.0

        k = array([[1., -1.],
                   [-1., 1.]])  # 1D rod

        Lambda = _Lambda(v1, debug=True)
        K = Lambda.T @ k @ Lambda
        Ki, Kj = K.shape

        # for testing
        #K = ones((Ki, Ki), 'float64')

        K2 = zeros((Ki*2, Kj*2), 'float64')
        if k_axial == 0.0 and k_torsion == 0.0:
            dofs = []
            n_ijv = []
            K2 = []
        elif k_torsion == 0.0: # axial; 2D or 3D
            K2 = K * k_axial
            dofs = array([
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
            dofs = array([
                i1+3, i1+4, i1+5,
                i2+3, i2+4, i2+5,
            ], 'int32')
            n_ijv = [
                # torsion
                (n1, 4), (n1, 5), (n2, 6),
                (n2, 4), (n2, 5), (n1, 6),
            ]

        else:  # axial + torsion; assume 3D
            # u1fx, u1fy, u1fz, u2fx, u2fy, u2fz
            K2[:Ki, :Ki] = K * k_axial

            # u1mx, u1my, u1mz, u2mx, u2my, u2mz
            K2[Ki:, Ki:] = K * k_torsion

            dofs = array([
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
        self.model.log.info('K =\n%s' % list_print(K / knorm))
        return K2, dofs, n_ijv

    def displacement_stress(self, model, positions, q, dofs):
        n = self.n
        o1 = zeros(n, 'float64')
        e1 = zeros(n, 'float64')
        f1 = zeros(n, 'float64')

        o4 = zeros(n, 'float64')
        e4 = zeros(n, 'float64')
        f4 = zeros(n, 'float64')


        As = self.get_area_by_element_id(self.property_id)
        Es = self.get_E_by_element_id(self.property_id)
        Gs = self.get_G_by_element_id(self.property_id)
        Js = self.get_J_by_element_id(self.property_id)
        Cs = self.get_c_by_element_id(self.property_id)

        for i in range(n):
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
                msg = 'invalid CTUBE length=0.0\n%s' % (self.__repr__())
                raise ZeroDivisionError(msg)

            #========================
            #mat = self.get_material_from_index(i)
            #jmat = searchsorted(mat.material_id, self.material_id[i])

            #E = mat.E[jmat]
            #G = mat.G[jmat]
            #G = self.G()

            #print("A=%g E=%g G=%g J=%g L=%g" % (A, E, G, J, L))
            #k_axial = A * E / L
            #k_torsion = G * J / L
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

            #print("Lsize = ", Lambda.shape)
            #print("qsize = ", q.shape)
            u_axial = array(Lambda) @ q_axial
            du_axial = -u_axial[0] + u_axial[1]
            u_torsion = array(Lambda) @ q_torsion
            du_torsion = -u_torsion[0] + u_torsion[1]

            #L = self.Length()
            #E = self.E()
            #A = self.area()

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
        obj = CTUBE(self.model)
        n = len(i)
        obj.n = n
        obj.i = n
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.property_id = self.property_id[i]
        obj.node_ids = self.node_ids[i, :]
        return obj
