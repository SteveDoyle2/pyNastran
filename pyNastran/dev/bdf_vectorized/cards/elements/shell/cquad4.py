"""
http://ssl.mit.edu/files/website/theses/SM-2011-McLinkoRyan.pdf
"""
import numpy as np
from numpy import array, zeros, arange, searchsorted, unique, cross
from numpy.linalg import norm  # type: ignore

from pyNastran.dev.bdf_vectorized.cards.elements.shell.shell_element import ShellElement

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_double_or_blank)
from pyNastran.utils.dev import list_print


class CQUAD4(ShellElement):
    type = 'CQUAD4'
    def __init__(self, model):
        ShellElement.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            #: Node IDs
            self.node_ids = zeros((ncards, 4), 'int32')
            self.zoffset = zeros(ncards, 'int32')

            self.theta = np.full(ncards, np.nan, 'float32')
            self.mcid = np.full(ncards, np.nan, 'int32')
            self.is_theta = zeros(ncards, 'bool')

            self.t_flag = zeros(ncards, 'int32')
            self.thickness = zeros((ncards, 4), float_fmt)

    def add_card(self, card, comment=''):
        i = self.i
        self.element_id[i] = integer(card, 1, 'eid')
        self.property_id[i] = integer(card, 2, 'pid')

        self.node_ids[i, :] = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer(card, 6, 'n4')
        ]

        theta_mcid = integer_double_or_blank(card, 7, 'theta_mcid', 0.0)

        if isinstance(theta_mcid, float):
            self.is_theta[i] = 1
            self.theta[i] = theta_mcid
        else:
            self.is_theta[i] = 0
            self.mcid[i] = theta_mcid


        #self.theta_mcid =
        #self.zoffset = double_or_blank(card, 7, 'zoffset', 0.0)
        #blank(card, 8, 'blank')
        #blank(card, 9, 'blank')

        self.t_flag[i] = integer_or_blank(card, 10, 'tflag', 0)
        self.thickness[i, :] = [
            double_or_blank(card, 11, 'T1', 1.0),
            double_or_blank(card, 12, 'T2', 1.0),
            double_or_blank(card, 13, 'T3', 1.0),
            double_or_blank(card, 14, 'T4', 1.0),
        ]
        self.i += 1

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self.is_theta = self.is_theta[i]
            self.theta = self.theta[i]
            self.mcid = self.mcid[i]

            self.thickness = self.thickness[i, :]
            self.t_flag = self.t_flag[i]
            assert self.node_ids.min() > 0
        else:
            self.element_id = array([], 'int32')
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
            for i, (eid, pid, nids) in enumerate(zip(self.element_id, self.property_id, self.node_ids)):
                self.element_id[i] = eid_map[eid]
                self.property_id[i] = pid_map[pid]
                self.node_ids[i, 0] = nid_map[nids[0]]
                self.node_ids[i, 1] = nid_map[nids[1]]
                self.node_ids[i, 2] = nid_map[nids[2]]
                self.node_ids[i, 3] = nid_map[nids[3]]

    def get_node_indicies(self, i=None):
        if i is None:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 1])
            i3 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 2])
            i4 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 3])
        else:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 1])
            i3 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 2])
            i4 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 3])
        return i1, i2, i3, i4

    #=========================================================================
    def _node_locations(self, xyz_cid0, i=None):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_node_index()
        i1, i2, i3, i4 = self.get_node_indicies(i)
        n1 = xyz_cid0[i1, :]
        n2 = xyz_cid0[i2, :]
        n3 = xyz_cid0[i3, :]
        n4 = xyz_cid0[i4, :]
        return n1, n2, n3, n4

    def _mass_area_normal(self, element_id=None, node_ids=None, xyz_cid0=None,
                          calculate_mass=True, calculate_area=True,
                          calculate_normal=True):
        """
        Gets the mass, area, and normals of the CQUAD4s on a per
        element basis.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : (nnodes, 3) float ndarray; default=None -> calculate
            the GRIDs in CORD2R=0
        calculate_mass : bool; default=True
            should the mass be calculated
        calculate_area : bool; default=True
            should the area be calculated
        calculate_normal : bool; default=True
            should the normals be calculated

        .. note:: If node_ids is None, the positions of all the GRID cards
                  must be calculated
        """
        if element_id is None:
            element_id = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
            property_id = self.property_id[i]

        n1, n2, n3, n4 = self._node_locations(xyz_cid0, i)
        if calculate_mass:
            calculate_area = True
        normal, A = _cquad4_normal_A(n1, n2, n3, n4,
                                     calculate_area=calculate_area, normalize=True)

        massi = None
        if calculate_mass:
            #t = self.model.properties_shell.get_thickness_by_property_id(property_id)  # PSHELL
            #nsm = self.model.properties_shell.get_nonstructural_mass_by_property_id(property_id)  # PSHELL
            #rho = self.model.properties_shell.get_density_by_property_id(property_id)  # MAT1
            mpa = self.model.properties_shell.get_mass_per_area_by_property_id(property_id)
            assert mpa is not None
            #assert t is not None
            #assert nsm is not None
            #assert rho is not None
            #print(type(nsm)
            #print("A  =%s" % A)
            #print("rho=%s" % rho)
            #print("t  =%s" % t)
            #print("nsm=%s" % nsm)
            #massi = rho * A * t + nsm
            massi = mpa * A
        #print("massi =", massi)
        return massi, A, normal

    def get_centroid_by_element_id(self, element_id=None, node_ids=None, xyz_cid0=None):
        if element_id is None:
            element_id = self.element_id
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
        return self.get_centroid_by_element_index(i, xyz_cid0)

    def get_centroid_by_element_index(self, i=None, xyz_cid0=None):
        n1, n2, n3, n4 = self._node_locations(xyz_cid0, i)
        return (n1 + n2 + n3 + n4) / 4.

    def get_area_by_element_index(self, i=None, xyz_cid0=None):
        n1, n2, n3, n4 = self._node_locations(xyz_cid0, i)
        normal, A = _cquad4_normal_A(n1, n2, n3, n4, calculate_area=True, normalize=False)
        return A

    def get_normal_by_element_index(self, i=None, xyz_cid0=None):
        n1, n2, n3, n4 = self._node_locations(xyz_cid0, i)
        normal, A = _cquad4_normal_A(n1, n2, n3, n4, calculate_area=False, normalize=True)
        return normal

    #=========================================================================
    def write_card(self, bdf_file, element_id=None, size=8):
        if self.n:
            #print('    self.n = %s' % self.n)
            if element_id is None:
                i = arange(self.n)
            else:
                if isinstance(element_id, int):
                    element_id = [element_id]
                #assert len(unique(element_id)) == len(element_id), unique(element_id)
                i = searchsorted(self.element_id, element_id)
                i = np.asarray(i)
            #if len(i) == 1:
                #i = i[0]
            #print('    iwrite = %s' % i)
            #print('    all_nodes = %s' % str(self.node_ids))
            self.write_card_by_index(bdf_file, i, size=size)

    def write_card_by_index(self, bdf_file, i, size=8):
        for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i]):
            card = ['CQUAD4', eid, pid, n[0], n[1], n[2], n[3]]
            bdf_file.write(print_card_8(card))

    def _verify(self, xref=True):
        self.get_mass_by_element_id()
        self.get_area_by_element_id()
        self.get_normal_by_element_id()

    def rebuild(self):
        raise NotImplementedError()

    def _positions(self, nids_to_get):
        """
        Gets the positions of a list of nodes

        Parameters
        ----------
        :param nids_to_get:  the node IDs to get as an NDARRAY
        :param node_ids:     the node IDs that contains all the nids_to_get
                             as an NDARRAY
        :param grids_cid_0:  the GRIDs as an (N, )  NDARRAY

        Returns
        -------
        grids2_cid_0 : (nnodes, 3) float ndarray
            the corresponding positions of the requested GRIDs
        """
        positions = self.model.grid.get_position_by_node_id(nids_to_get)
        #grids2_cid_0 = grids_cid0[searchsorted(node_ids, nids_to_get), :]
        #return grids2_cid_0
        return positions

    #def slice_by_index(self, i):
        #i = self._validate_slice(i)
        #obj = CQUAD4(self.model)
        #obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        #obj.element_id = self.element_id[i]
        #obj.property_id = self.property_id[i]
        #obj.node_ids = self.node_ids[i, :]
        #obj.zoffset = self.zoffset[i]
        #obj.t_flag = self.t_flag[i]
        #obj.thickness = self.thickness[i, :]
        #return obj

    def get_stiffness_matrix(self, i, model, positions, index0s, knorm=1.0):  # CROD/CONROD
        """gets the stiffness matrix for CQUAD4/PSHELL"""
        area = self.get_area_by_element_index(i)
        #print('self.thickness =', self.thickness)
        #print('self.i =', i)

        n1, n2, n3, n4 = self.node_ids[i, :]
        pid = self.property_id[i]
        prop = self.model.properties_shell.get_property_by_property_id(pid)
        #prop.get_
        #print(prop)

        # TODO: nodal thickness not supported
        #thickness = self.thickness.flatten()[i]

        thickness = prop.thickness#.flatten()[i]
        #print('thickness =', thickness)

        mid1 = prop.material_id[0]
        mat = self.model.materials.get_shell_material(mid1)
        assert prop.material_id2[0] == -1, prop.material_id2
        assert prop.material_id3[0] == -1, prop.material_id3
        assert prop.material_id4[0] == -1, prop.material_id4
        #print(mat)
        E = mat.E[0]
        nu = mat.nu[0]
        #print('area=%s thickness=%s E=%e nu=%s' % (area, thickness, E, nu))
        #sdd
        i1 = index0s[n1]
        i2 = index0s[n2]
        i3 = index0s[n3]
        i4 = index0s[n4]

        xyz1 = positions[n1]
        xyz2 = positions[n2]
        xyz3 = positions[n3]
        xyz4 = positions[n4]
        xy = np.vstack([
            xyz1,
            xyz2,
            xyz3,
            xyz4,
        ])[:, :2]

        #print(xy)

        centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.

        #normal = self.Normal()
        is_theta = self.is_theta[i]
        if is_theta:
            theta = self.theta[i]

            if theta != 0.0:
                self.model.log.debug(self.model.coords)
                mcid_ref = self.model.coords.get_coord_index_by_coord_id(0)
                #print('mcid_ref\n', mcid_ref)
                #assert self.theta == 0.0, self.theta
                raise NotImplementedError('theta=%r' % theta)
        else:
            mcid = self.mcid[i]
            #print(mcid)
            #print(self.model.coords)
            mcid_ref = self.model.coords.get_coord_index_by_coord_id(mcid)
            i = mcid_ref.i
            jmat = np.cross(normal, i) # k x i
            jmat /= np.linalg.norm(jmat)
            imat = np.cross(jmat, normal)
            T = np.vstack([imat, jmat, normal])
            self.model.log.debug(T)
        if 0:
            xyz = np.vstack([
                xyz1,
                xyz2,
                xyz3,
                xyz4,
            ]).dot(T)

        dofs = array([
            i1, i1+1,
            i2, i2+1,
            i3, i3+1,
            i4, i4+1,
        ], 'int32')

        n_ijv = [
            # axial
            (n1, 1), (n1, 2),
            (n2, 1), (n2, 2),
            (n3, 1), (n3, 2),
            (n4, 1), (n4, 2),
        ]

        #n1 = 0.25 * (1 - u) * (1 - v)
        #n2 = 0.25 * (1 + u) * (1 - v)
        #n3 = 0.25 * (1 + u) * (1 + v)
        #n4 = 0.25 * (1 - u) * (1 + v)
        #wti = -0.57735
        #wtj = -0.57735

        is_heat_transfer = False
        if is_heat_transfer:
            u = wti
            v = wtj
            Ji = array([
                [v - 1.0, -v + 1.0, v + 1.0, -v - 1.0],
                [u - 1.0, -u - 1.0, u + 1.0, -u + 1.0],
            ]) / 4.
            J = Ji.dot(xy)
            Jinv = np.linalg.inv(J)
            detJ = np.linalg.det(J)
            darea = detJ
            #print('Ji*4 =\n', Ji*4.)
            #print('J =\n', J)
            #print('Jinv =\n', Jinv)

            #B2 = array([
                #[-1.0 - wti, 1.0 - wti, 1.0 - wti, -1.0 + wti],
                #[-1.0 - wtj, -1.0 + wtj, 1.0 - wtj, 1.0 + wtj],
            #])
            #print('B2 =\n', B2)
            #B = B2.dot(xy)
            B = Jinv.dot(Ji)
            #print('B =\n', B)
            k = 1. * darea
            #print('dA =', darea)

            K = k * B.T.dot(B)
            #print('K =\n', K)
        else:
            K = np.zeros((8, 8), dtype='float64')
            #u = wti
            #v = wtj
            wts = [-0.57735, 0.57735]
            for u in wts:
                for v in wts:
                    Ji = array([
                        [v - 1.0, -v + 1.0, v + 1.0, -v - 1.0],
                        [u - 1.0, -u - 1.0, u + 1.0, -u + 1.0],
                    ]) / 4.
                    J = Ji.dot(xy)
                    Jinv = np.linalg.inv(J)
                    det_j = np.linalg.det(J)
                    darea = det_j


                    B1 = Jinv.dot(Ji)
                    #print('B1 =\n', B1)
                    N1x, N2x, N3x, N4x = B1[0, :]
                    N1y, N2y, N3y, N4y = B1[1, :]
                    #print('Nix =', B1[0, :])


                    #N1x, N2x, N3x, N4x = v - 1.0, -v + 1.0, v + 1.0, -v - 1.0
                    #N1y, N2y, N3y, N4y = u - 1.0, -u - 1.0, u + 1.0, -u + 1.0
                    B = array([
                        [N1x, 0., N2x, 0., N3x, 0., N4x, 0.],
                        [0., N1y, 0., N2y, 0., N3y, 0., N4y],
                        [N1y, N1x, N2y, N2x, N3y, N3x, N4y, N4x]
                    ])
                    #print('B =\n', B)

                    #E = 1.0
                    #nu = 0.25
                    denom = 1 - nu**2
                    #C = E/(1 - (poisson^2))*[1 poisson 0; poisson 1 0;0 0 ((1-poisson)/2)];
                    C = np.array([
                        [E/denom, nu*E/denom, 0.],
                        [nu*E/denom, E/denom, 0.],
                        [0., 0., E/(2.0*(1.0+nu))],
                    ], dtype='float64')
                    #print('C =\n', C)
                    #print('thickness =', thickness)
                    Ki = (B.T @ C @ B) * (thickness * darea)
                    #print('Ki(%s,%s) =%s\n' % (u, v, Ki))
                    #print('Ki(%s,%s) =\n%s\n' % (u, v, list_print(Ki, '%.4e')))
                    K += Ki

        #K *= (thickness * darea)

        knorm = 1.0
        self.model.log.info("K_norm / %s = \n" % knorm + list_print(K / knorm, float_fmt='%-4.4f'))
        return(K, dofs, n_ijv)

    def displacement_stress(self, model, positions, q, dofs):
        n = self.n
        stress = zeros(n, 'float64')
        strain = zeros(n, 'float64')
        force = zeros(n, 'float64')

        i = self.get_element_index_by_element_id(self.element_id)
        As = self.get_area_by_element_index(i)
        #Gs = self.model.prod.get_G_by_property_id(self.property_id)
        #Es = self.model.prod.get_E_by_property_id(self.property_id)
        #Js = self.model.prod.get_J_by_property_id(self.property_id)
        #Cs = self.model.prod.get_c_by_property_id(self.property_id)

        for i in range(n):
            A = As[i]
            n1, n2, n3, n4 = self.node_ids[i, :]


            xyz1 = positions[n1]
            xyz2 = positions[n2]
            xyz3 = positions[n3]
            xyz4 = positions[n4]

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
                q[n21], q[n22], q[n23],
            ])
            q_torsion = array([
                q[n14], q[n15], q[n16],
                q[n24], q[n25], q[n26]
            ])
            #v1 = xyz1 - xyz2
            #L = norm(xyz1 - xyz2)
            #if L == 0.0:
                #msg = 'invalid CROD length=0.0\n%s' % (self.__repr__())
                #raise ZeroDivisionError(msg)

            #========================

            #print("A=%g E=%g G=%g J=%g L=%g" % (A, E, G, J, L))
            k_axial = A * E / L
            k_torsion = G * J / L
            #k_axial = 1.0
            #k_torsion = 2.0

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
            du_axial = u_axial[0] - u_axial[1]
            u_torsion = array(Lambda) @ q_torsion
            du_torsion = u_torsion[0] - u_torsion[1]

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



def _cquad4_normal_A(n1, n2, n3, n4, calculate_area=True, normalize=True):
    v13 = n1 - n3
    v24 = n2 - n4
    normal = cross(v13, v24)

    A = None
    if calculate_area:
        n = norm(normal, axis=1)
        A = 0.5 * n
    elif normalize:
        n = norm(normal, axis=1)

    if normalize:
        #n = len(_norm)
        #for i in range(n):
            #normal[i] /= _norm[i]
        #normal /= n
        nx = len(n)
        normal /= n.reshape(nx, 1)

    return normal, A
