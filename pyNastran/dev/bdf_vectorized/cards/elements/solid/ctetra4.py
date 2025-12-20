from itertools import count

import numpy as np
from numpy import zeros, arange, dot, cross, searchsorted, array, eye, ones

#from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import integer
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

from pyNastran.dev.bdf_vectorized.cards.elements.solid.solid_element import SolidElement

def volume4(xyz1, xyz2, xyz3, xyz4):
    r"""
    Gets the volume, :math:`V`, of the tetrahedron.

    .. math:: V = \frac{(a-d) \cdot \left( (b-d) \times (c-d) \right) }{6}
    """
    V = -dot((xyz1 - xyz4), cross(xyz2 - xyz4, xyz3 - xyz4)) / 6.
    #V = 1/6. * np.det(
        #np.hstack(
            #[1., 1., 1., 1.],
            #np.vstack(n1, n2, n3, n4).T,
        #),
    #)
    return V


class CTETRA4(SolidElement):
    type = 'CTETRA4'
    nnodes = 4
    def __init__(self, model):
        """
        Defines the CTETRA object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        SolidElement.__init__(self, model)

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i

        eid = integer(card, 1, 'element_id')
        if comment:
            self.set_comment(eid, comment)

        #: Element ID
        self.element_id[i] = eid
        #: Property ID
        self.property_id[i] = integer(card, 2, 'property_id')
        #: Node IDs
        nids = array([
            integer(card, 3, 'node_id_1'),
            integer(card, 4, 'node_id_2'),
            integer(card, 5, 'node_id_3'),
            integer(card, 6, 'node_id_4'),
        ], dtype='int32')
        assert 0 not in nids, '%s\n%s' % (nids, card)
        self.node_ids[i, :] = nids
        assert len(card) == 7, 'len(CTETRA4 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

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
                print(self.print_card(i))
                self.element_id[i] = eid_map[eid]
                self.property_id[i] = pid_map[pid]
                self.node_ids[i, 0] = nid_map[nids[0]]
                self.node_ids[i, 1] = nid_map[nids[1]]
                self.node_ids[i, 2] = nid_map[nids[2]]
                self.node_ids[i, 3] = nid_map[nids[3]]

    def get_mass_matrix(self, i, model, positions, index0s):
        r"""
        A mass matrix is a discrete representation of a continuous mass distribution.
        To compute our mass matrix for a tetrahedral element with linear shape
        functions we need the formula (pp. 266 in Cook)

                                                         a!b!c!d!
         \int_V N_1^a  N_2^b N_3^c N_4^d dV = 6V --------------------------     (**)
                                                    (3 + a + b +c + d)!

        A consistent element mass matrix (pp. 376 Cook) is defined as

          m = \int_V \rho N^T N dV           (***)

        This equation can be derived from work balance, the details of which is unimportant
        here (look Cook pp. 375-376 for details).
        Assumping \rho is constant over each tetrahedral element and using the linear shape
        functions the above definition (***) results in

                           |N_1|
        m =  \rho \int_V   |N_2|  |N_1 N_2 N_3 N_4| dV
                           |N_3|
                           |N_4|

                              |(N_1 N_1)   (N_1 N_2)   (N_1 N_3)   (N_1 N_4)|
        m =   \rho    \int_V  |(N_2 N_1)   (N_2 N_2)   (N_2 N_3)   (N_2 N_4)| dV
                              |(N_3 N_1)   (N_3 N_2)   (N_3 N_3)   (N_3 N_4)|
                              |(N_4 N_1)   (N_4 N_2)   (N_4 N_3)   (N_4 N_4)|

        by (**)

                        | 2 1 1 1|
        m = \rho  V/20  | 1 2 1 1|          (****)
                        | 1 1 2 1|
                        | 1 1 1 2|

                      V
        m_ij =  \rho --- (1+delta_ij)
                      20

        in 3D this means that for the tetrahedral element

                       | 2 2 2  1 1 1  1 1 1  1 1 1 |
                       | 2 2 2  1 1 1  1 1 1  1 1 1 |
                       | 2 2 2  1 1 1  1 1 1  1 1 1 |
                       |                            |
                       | 1 1 1  2 2 2  1 1 1  1 1 1 |
                       | 1 1 1  2 2 2  1 1 1  1 1 1 |
                   V   | 1 1 1  2 2 2  1 1 1  1 1 1 |
        Me = \rho ---  |                            |
                   20  | 1 1 1  1 1 1  2 2 2  1 1 1 |
                       | 1 1 1  1 1 1  2 2 2  1 1 1 |
                       | 1 1 1  1 1 1  2 2 2  1 1 1 |
                       |                            |
                       | 1 1 1  1 1 1  1 1 1  2 2 2 |
                       | 1 1 1  1 1 1  1 1 1  2 2 2 |
                       | 1 1 1  1 1 1  1 1 1  2 2 2 |

        Notice that in order to obtain the global/system mass matrix an assembly similar to the
        stiffness matrix assembly must be carried out. Further, the global M matrix will
        have the same sub-block pattern as the global K matrix.

        A consistent mass matrix is often not used in computer graphics. Instead and
        ad-hoc approach named ``lumped'' mass matrix is applied.
        The lumped mass matrix is obtained by placing particle masses at the nodes.
        This corresponds to shifting all the masses in the rows of (****) onto the
        diagonal. In 3D this yields the element mass matrix

                       | 1 0 0  0 0 0  0 0 0  0 0 0 |
                       | 0 1 0  0 0 0  0 0 0  0 0 0 |
                       | 0 0 1  0 0 0  0 0 0  0 0 0 |
                       |                            |
                       | 0 0 0  1 0 0  0 0 0  0 0 0 |
                       | 0 0 0  0 1 0  0 0 0  0 0 0 |
                   V   | 0 0 0  0 0 1  0 0 0  0 0 0 |
        Me = \rho ---  |                            |
                   4   | 0 0 0  0 0 0  1 0 0  0 0 0 |
                       | 0 0 0  0 0 0  0 1 0  0 0 0 |
                       | 0 0 0  0 0 0  0 0 1  0 0 0 |
                       |                            |
                       | 0 0 0  0 0 0  0 0 0  1 0 0 |
                       | 0 0 0  0 0 0  0 0 0  0 1 0 |
                       | 0 0 0  0 0 0  0 0 0  0 0 1 |

        Thus a lumped mass matrix is diagonal whereas a consistent mass matrix
        is not. Observe that the global mass matrix would also diagonal and the
        assembly is simplified to an iteration over all tetrahedra, while
        incementing the nodal mass by one fourth of the tetrahedral mass.

         for each node n
           mass(n) = 0
         next n
         for each tetrahedron e
           mass(n_i) += \rho_e Ve / 4
           mass(n_j) += \rho_e Ve / 4
           mass(n_k) += \rho_e Ve / 4
           mass(n_m) += \rho_e Ve / 4
         next e

        where n_i,n_j,n_k and n_m are the four nodes of the e'th tetrahedron.

        The advantage of lumping is less storage and higher performance. On the downside
        lumping introduces a discontinouty in the displacement field.

        Obrien.shen state that the errors in lumping is negligeble for small-size course
        meshes used in computer graphics. However, for finer meshes the errors becomes
        noticeable.

        There do exist other approaches for computing mass matrices, even methods which
        combine other methods. We refer the interested reader to Cook for more details. Here
        we have limited our selfes to the two most common methods.

        It is worthwhile to notice that under the reasonable assumptions that V and \rho are
        positive for all elements both the element mass matrices and the global mass matrices
        are symmetric positive definite matrices.

        http://image.diku.dk/svn/OpenTissue/archieve/silcowitz/OpenTissue/dynamics/fem/fem_compute_mass.h
        """
        is_lumped = True
        is_consistent = False
        nnodes = 4
        ndof = 3 * nnodes
        pid = self.property_id[i]
        rho = self.model.elements.properties_solid.psolid.get_density_by_property_id(pid)[0]

        n0, n1, n2, n3 = self.node_ids[i, :]
        V = volume4(positions[self.node_ids[i, 0]],
                    positions[self.node_ids[i, 1]],
                    positions[self.node_ids[i, 2]],
                    positions[self.node_ids[i, 3]])

        mass = rho * V
        if is_lumped:
            mi = mass / 4.
            nnodes = 4
            M = eye(ndof, dtype='float32')
        else:
            mi = mass / 20.
            M = ones((ndof, ndof), dtype='float32')
            for i in range(nnodes):
                j = i * 3
                M[j:j+3, j:j+3] = 2.
        M *= mi
        dofs, nijv = self.get_dofs_nijv(index0s, n0, n1, n2, n3)
        return M, dofs, nijv

    def get_stiffness_matrices(self, model, positions, index0s):
        out = []

        # volume coordinates
        # FEM: Volume I (Zienkiewicz) p.186
        volume6 = volume * 6
        L1 = (a1 + b1 * x + c1 * y + d1 * z) / volume6
        L2 = (a2 + b2 * x + c2 * y + d2 * z) / volume6
        L3 = (a3 + b3 * x + c3 * y + d3 * z) / volume6

        # FEM: Volume I (Zienkiewicz) p.186
        #x = L1*x1 + L2*x2 + L3*x3 + L4*x4
        #y = L1*y1 + L2*y2 + L3*y3 + L4*y4
        #z = L1*z1 + L2*z2 + L3*z3 + L4*z4
        #1 = L1 + L2 + L3 + L4


        for i in range(self.n):
            K, dofs, nijv = self.get_stiffness_matrix(
                i, model, self.positions, index0s)
            out.append(K, dofs, nijv)
        self.add_stiffness(K, dofs, nijv)

    def get_stiffness_matrix(self, i, model, positions, index0s):
        nnodes = 4
        ndof = 3 * nnodes
        pid = self.property_id[i]
        prop = self.model.elements.properties_solid.psolid
        rho = prop.get_density_by_property_id(pid)[0]

        n0, n1, n2, n3 = self.node_ids[i, :]
        xyz1 = positions[self.node_ids[i, 0]]
        xyz2 = positions[self.node_ids[i, 1]]
        xyz3 = positions[self.node_ids[i, 2]]
        xyz4 = positions[self.node_ids[i, 3]]
        vol = volume4(xyz1, xyz2, xyz3, xyz4)

        stiffness = rho * vol
        ki = stiffness / 4.
        nnodes = 4
        K = eye(ndof, dtype='float32')  # not done...


        u = 0.
        v = 0.
        #wts = [-0.57735, 0.57735]
        #for u in wts:
            #for v in wts:
        Ji = array([
            [v - 1.0, -v + 1.0, v + 1.0, -v - 1.0],
            [u - 1.0, -u - 1.0, u + 1.0, -u + 1.0],
        ]) / 4.
        #J = Ji @ xy
        #Jinv = np.linalg.inv(J)
        #det_j = np.linalg.det(J)
        #darea = det_j


        #B1 = Jinv @ Ji
        #print('B1 =\n', B1)
        #N1x, N2x, N3x, N4x = B1[0, :]
        #N1y, N2y, N3y, N4y = B1[1, :]
        #print('Nix =', B1[0, :])


        vol_matrix = np.hstack(
            [1., 1., 1., 1.],
            np.vstack([xyz1, xyz2, xyz3, xyz4]).T,
        )
        ivol_matrix = np.linalg.inv(vol_matrix)
        a1, b1, c1 = ivol_matrix[0, 1:]
        a2, b2, c2 = ivol_matrix[1, 1:]
        a3, b3, c3 = ivol_matrix[2, 1:]
        a4, b4, c4 = ivol_matrix[3, 1:]

        #N1x, N2x, N3x, N4x = v - 1.0, -v + 1.0, v + 1.0, -v - 1.0
        #N1y, N2y, N3y, N4y = u - 1.0, -u - 1.0, u + 1.0, -u + 1.0
        B = array([
            [a1, 0., 0., a2, 0., 0., a3, 0., 0., a4, 0., 0.],
            [0., b1, 0., 0., b2, 0., 0., b3, 0., 0., b4, 0.],
            [0., 0., c1, 0., 0., c2, 0., 0., c3, 0., 0., c4],
            [b1, a1, 0., b2, a2, 0., b3, a3, 0., b4, a4, 0.],
            [0., c1, b1, 0., c2, b2, 0., c3, b3, 0., c4, b4],
            [c1, 0., a1, c2, 0., a2, c3, 0., a3, c4, 0., a4],
        ]) / (6 * vol)

        #N = array([
            #[N1, 0., 0., N2, 0., 0., N3, 0., N4, 0., 0.],
            #[0., N1, 0., 0., N2, 0., 0., N3, 0., N4, 0.],
            #[0., 0., N1, 0., 0., N2, 0., 0., N3, 0., N4],
        #])
        #print('B =\n', B)

        #E = 1.0
        #nu = 0.25

        mid1 = prop.material_id[0]
        mat = self.model.materials.get_solid_material(mid1)
        print(mat)
        E = mat.E[0]
        nu = mat.nu[0]
        G = mat.G[0]

        # [sigma] = [C] * [epsilon]
        #denom = 1 - nu**2
        #C = np.zeros((6, 6), dtype='float64')
        #outside = E / ((1 + nu) * (1 - 2 * nu))
        #C[0, 0] = C[1, 1] = C[2, 2] = (1 - nu) * outside
        #C[3, 3] = C[4, 4] = C[5, 5] = (0.5 - nu) * outside

        if 0:
            ## [stress] = [E] [strain]
            #emat = np.zeros((5, 5), dtype='float64')
            #emat[0, 0] = emat[1, 1] = E / denom
            #emat[1, 0] = emat[0, 1] = (E * nu) / denom
            #emat[2, 2] = emat[3, 3] = emat[4, 4] = G


            ## [M] = [D] * [bending]
            #dmat = np.zeros((5, 5), dtype='float64')
            #D = E * h**3 / (12 * denom)
            #dmat[0, 0] = dmat[1, 1] = D
            #dmat[1, 0] = dmat[0, 1] = D * nu
            #dmat[2, 2] = D * (1. - nu) / 2.
            #dmat[3, 3] = emat[4, 4] = G * h

            # FEM: Volume I (Zienkiewicz) p.132
            dmat2 = np.array(6, 6)
            dmat2[0, 0] = dmat2[1, 1] = dmat2[2, 2] = 1 - nu
            dmat2[0, 1] = dmat2[0, 2] = dmat2[1, 0] = dmat2[2, 0] = nu
            dmat2[3, 3] = dmat2[4, 4] = dmat[5, 5] = (1 - 2 * nu) / 2.
            dmat2 *= E / ((1 + nu) * (1 - 2 * nu))

        #print('C =\n', C)
        #print('thickness =', thickness)
        Ki = B.T @ C @ B
        #print('Ki(%s,%s) =%s\n' % (u, v, Ki))
        #print('Ki(%s,%s) =\n%s\n' % (u, v, list_print(Ki, '%.4e')))
        K += Ki

        #K *= ki
        dofs, nijv = self.get_dofs_nijv(index0s, n0, n1, n2, n3)
        return K, dofs, nijv

    def get_dofs_nijv(self, index0s, n0, n1, n2, n3):
        i0 = index0s[n0]
        i1 = index0s[n1]
        i2 = index0s[n2]
        i3 = index0s[n3]
        dofs = array([
            i0, i0+1, i0+2,
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,
            i3, i3+1, i3+2,
        ], 'int32')
        nijv = [
            # translation
            (n0, 1), (n0, 2), (n0, 3),
            (n1, 1), (n1, 2), (n1, 3),
            (n2, 1), (n2, 2), (n2, 3),
            (n3, 1), (n3, 2), (n3, 3),
        ]
        return dofs, nijv

    def _verify(self, xref=True):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            c = self.centroid()
            v = self.volume()
            assert isinstance(v, float)
            for i in range(3):
                assert isinstance(c[i], float)

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

    def _get_node_locations_by_index(self, i, xyz_cid0):
        """
        :param i:        None or an array of node IDs
        :param xyz_cid0: the node positions as a dictionary
        """
        grid = self.model.grid
        get_node_index_by_node_id = self.model.grid.get_node_index_by_node_id
        node_ids = self.node_ids

        #msg = ', which is required by %s' % self.type
        i1, i2, i3, i4 = self.get_node_indicies(i)
        n1 = xyz_cid0[i1, :]
        n2 = xyz_cid0[i2, :]
        n3 = xyz_cid0[i3, :]
        n4 = xyz_cid0[i4, :]
        return n1, n2, n3, n4

    def get_volume_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more elements.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the volume be summed
        """
        n1, n2, n3, n4 = self._get_node_locations_by_element_id(element_id, xyz_cid0)

        V = zeros(n1.shape[0], self.model.float_fmt)
        for i, n1i, n2i, n3i, n4i in zip(count(), n1, n2, n3, n4):
            V[i] = volume4(n1i, n2i, n3i, n4i)
            i += 1
        return V

    def get_mass_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the mass for one or more CTETRA elements.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the centroid be summed
        """
        if element_id is None:
            element_id = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_node_index()

        V = self.get_volume_by_element_id(element_id, xyz_cid0)
        mid = self.model.properties_solid.get_material_id_by_property_id(self.property_id)
        rho = self.model.materials.get_density_by_material_id(mid)

        mass = V * rho
        if total:
            mass = mass.sum()
        return mass

    def get_centroid_volume(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more elements.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        :param total: should the volume be summed; centroid be averaged (default=False)

        .. seealso:: CTETRA4.volume() and CTETRA4.centroid for more information.
        """
        n1, n2, n3, n4 = self._get_node_locations_by_element_id(element_id, xyz_cid0)
        n = len(element_id)
        volume = zeros(n, self.model.float_fmt)

        i = 0
        for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
            volume[i] = volume4(n1i, n2i, n3i, n4i)
            i += 1

        centroid = (n1 + n2 + n3 + n4) / 4.0
        if total:
            centroid = centroid.mean()
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        assert volume.min() > 0.0, 'volume.min() = %f' % volume.min()
        return centroid, volume

    def get_centroid_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid for one or more elements.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the centroid be averaged
        """
        n1, n2, n3, n4 = self._get_node_locations_by_element_id(element_id, xyz_cid0)
        centroid = (n1 + n2 + n3 + n4) / 4.0
        if total:
            centroid = centroid.mean(axis=0)
        return centroid

    #def get_face_nodes(self, nid, nid_opposite):
        #raise NotImplementedError()
        #nids = self.node_ids[:4]
        #indx = nids.index(nid_opposite)
        #nids.pop(indx)
        #return nids

    def write_card(self, bdf_file, size=8, element_id=None):
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_id)

            if size == 16 or max(self.element_id[i].max(), self.property_id[i].max(),
                                 self.node_ids[i, :].max()) > 1000000000:
                msg = ('CTETRA  %16i%16i%16i%16i\n'
                       '        %16i%16i\n')
                for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                    if eid in self._comments:
                        bdf_file.write(self._comments[eid])
                    data = [eid, pid] + list(n)
                    bdf_file.write(msg)
            else:
                msg = 'CTETRA  %8i%8i%8i%8i%8i%8i\n'
                for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                    if eid in self._comments:
                        bdf_file.write(self._comments[eid])
                    data = [eid, pid] + list(n)
                    bdf_file.write(msg % tuple(data))
