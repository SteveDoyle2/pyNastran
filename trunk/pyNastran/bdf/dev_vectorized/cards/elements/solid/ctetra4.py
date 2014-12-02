from six.moves import zip, range
from itertools import count

from numpy import zeros, arange, dot, cross, searchsorted, array, where, asarray, eye
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (fields, integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)

from pyNastran.bdf.dev_vectorized.cards.elements.solid.solid_element import SolidElement

def volume4(n1, n2, n3, n4):
    r"""
    Gets the volume, :math:`V`, of the tetrahedron.

    .. math:: V = \frac{(a-d) \cdot \left( (b-d) \times (c-d) \right) }{6}
    """
    V = -dot((n1 - n4), cross(n2 - n4, n3 - n4)) / 6.
    return V


class CTETRA4(SolidElement):
    type = 'CTETRA4'
    nnodes = 4
    def __init__(self, model):
        """
        Defines the CTETRA object.

        :param self: the CTETRA object
        :param model: the BDF object
        """
        SolidElement.__init__(self, model)

    def add(self, card, comment=''):
        i = self.i

        #comment = self._comments[i]
        eid = integer(card, 1, 'element_id')
        #if comment:
            #self._comments[eid] = comment

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
        assert len(card) == 7, 'len(CTETRA4 card) = %i' % len(card)
        self.i += 1

    def get_mass_matrix(self, i, model, positions, index0s):
        """
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

        The advantage of lumping is less storage and higher performace. On the downside
        lumping introduces a discontinouty in the displacement field.

        Obrien.shen state that the errors in lumping is negligeble for small-size course
        meshes used in computer graphics. However, for finer meshes the errors becomes
        noticeable.

        There do exist other approaches for computing mass matrices, even mehtods which
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

    def get_stiffness_matrix(self, i, model, positions, index0s):
        nnodes = 4
        ndof = 3 * nnodes
        pid = self.property_id[i]
        rho = self.model.elements.properties_solid.psolid.get_density_by_property_id(pid)[0]

        n0, n1, n2, n3 = self.node_ids[i, :]
        V = volume4(positions[self.node_ids[i, 0]],
                    positions[self.node_ids[i, 1]],
                    positions[self.node_ids[i, 2]],
                    positions[self.node_ids[i, 3]])

        stiffness = rho * V
        ki = stiffness / 4.
        nnodes = 4
        K = eye(ndof, dtype='float32')  # not done...
        K *= ki
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
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()
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

    def _get_node_locations_by_index(self, i, xyz_cid0):
        """
        :param i:        None or an array of node IDs
        :param xyz_cid0: the node positions as a dictionary
        """
        grid = self.model.grid
        get_node_index_by_node_id = self.model.grid.get_node_index_by_node_id
        node_ids = self.node_ids

        msg = ', which is required by %s' % self.type
        n1 = xyz_cid0[get_node_index_by_node_id(node_ids[i, 0], msg), :]
        n2 = xyz_cid0[get_node_index_by_node_id(node_ids[i, 1], msg), :]
        n3 = xyz_cid0[get_node_index_by_node_id(node_ids[i, 2], msg), :]
        n4 = xyz_cid0[get_node_index_by_node_id(node_ids[i, 3], msg), :]
        return n1, n2, n3, n4

    def get_volume_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more elements.

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)
        """
        n1, n2, n3, n4 = self._get_node_locations_by_element_id(element_id, xyz_cid0)

        V = zeros(n1.shape[0], self.model.float)
        for i, n1i, n2i, n3i, n4i in zip(count(), n1, n2, n3, n4):
            V[i] = volume4(n1i, n2i, n3i, n4i)
            i += 1
        return V

    def get_mass_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the mass for one or more CTETRA elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if element_id is None:
            element_id = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_index()

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

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed; centroid be averaged (default=False)

        ..see:: CTETRA4.volume() and CTETRA4.centroid for more information.
        """
        n1, n2, n3, n4 = self._get_node_locations_by_element_id(element_id, xyz_cid0)
        n = len(element_id)
        volume = zeros(n, self.model.float)

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

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be averaged (default=False)
        """
        n1, n2, n3, n4 = self._get_node_locations_by_element_id(element_id, xyz_cid0)
        centroid = (n1 + n2 + n3 + n4) / 4.0
        if total:
            centroid = centroid.mean(axis=0)
        return centroid

    def get_face_nodes(self, nid, nid_opposite):
        raise NotImplementedError()
        nids = self.nodeIDs()[:4]
        indx = nids.index(nid_opposite)
        nids.pop(indx)
        return nids

    def write_bdf(self, f, size=8, element_id=None):
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_id)

            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                card = ['CTETRA', eid, pid, n[0], n[1], n[2], n[3]]
                f.write(print_card_8(card))

