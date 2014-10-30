from six.moves import zip
from itertools import count

from numpy import zeros, arange, dot, cross, searchsorted, array, where, asarray, eye
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card
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
    def __init__(self, model):
        """
        Defines the CTETRA4 object.

        :param self: the CTETRA4 object
        :param model: the BDF object
        """
        SolidElement.__init__(self, model)

    def build(self):
        cards = self._cards
        ncards = len(cards)

        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            self.element_id = zeros(ncards, 'int32')
            self.property_id = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 4), 'int32')

            #comments = {}
            for i, card in enumerate(cards):
                #comment = self._comments[i]
                eid = integer(card, 1, 'eid')
                #if comment:
                    #self._comments[eid] = comment

                #: Element ID
                self.element_id[i] = eid
                #: Property ID
                self.property_id[i] = integer(card, 2, 'pid')
                #: Node IDs
                self.node_ids[i, :] = [
                    integer(card, 3, 'nid1'),
                    integer(card, 4, 'nid2'),
                    integer(card, 5, 'nid3'),
                    integer(card, 6, 'nid4'),
                ]
                assert len(card) == 7, 'len(CTETRA4 card) = %i' % len(card)

            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self._cards = []
            self._comments = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

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
        rho = self.model.elements.properties_solid.psolid.get_density(pid)[0]

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
            for i in xrange(nnodes):
                j = i * 3
                M[j:j+3, j:j+3] = 2.
        M *= mi
        dofs, nijv = self.get_dofs_nijv(index0s, n0, n1, n2, n3)
        return M, dofs, nijv

    def get_stiffness_matrix(self, i, model, positions, index0s):
        nnodes = 4
        ndof = 3 * nnodes
        pid = self.property_id[i]
        rho = self.model.elements.properties_solid.psolid.get_density(pid)[0]

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
        return M, dofs, nijv

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

    def _node_locations_element_id(self, element_id=None, xyz_cid0=None):
        if element_id is None:
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_positions()
        return self._node_locations_i(i, xyz_cid0)

    def _node_locations_i(self, i, xyz_cid0):
        """
        :param i:        None or an array of node IDs
        :param xyz_cid0: the node positions as a dictionary
        """
        index_map = self.model.grid.index_map
        node_ids = self.node_ids
        n1 = xyz_cid0[index_map(node_ids[i, 0]), :]
        n2 = xyz_cid0[index_map(node_ids[i, 1]), :]
        n3 = xyz_cid0[index_map(node_ids[i, 2]), :]
        n4 = xyz_cid0[index_map(node_ids[i, 3]), :]
        return n1, n2, n3, n4

    def get_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CTETRA4 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)
        """
        if element_ids is None:
            element_ids = self.element_id
        n1, n2, n3, n4 = self._node_locations(xyz_cid0)

        n = len(element_ids)
        V = zeros(n, self.model.float)

        for i, n1i, n2i, n3i, n4i in zip(count(), n1, n2, n3, n4):
            V[i] = volume4(n1i, n2i, n3i, n4i)
            i += 1
        return V

    def get_centroid_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CTETRA4 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed; centroid be averaged (default=False)

        ..see:: CTETRA4.volume() and CTETRA4.centroid for more information.
        """
        n1, n2, n3, n4 = self._node_locations_element_id(xyz_cid0)
        n = len(element_ids)
        volume = zeros(n, self.model.float)

        i = 0
        for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
            volume[i] = volume4(n1i, n2i, n3i, n4i)
            i += 1

        centroid = (n1 + n2 + n3 + n4) / 4.0
        if total:
            centroid = centroid.mean()
        assert volume.min() > 0.0, 'volume.min() = %f' % volume.min()
        return centroid, volume

    def get_centroid(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid for one or more CTETRA elements.

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be averaged (default=False)
        """
        n1, n2, n3, n4 = self._node_locations_element_id(element_id, xyz_cid0)
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
                f.write(print_card(card))

    #def slice_by_index(self, i):
        #i = asarray(i)
        #obj = CTETRA4(self.model)
        #obj.n = len(i)
        ##obj._cards = self._cards[i]
        ##obj._comments = obj._comments[i]
        ##obj.comments = obj.comments[i]
        #obj.element_id = self.element_id[i]
        #obj.property_id = self.property_id[i]
        #obj.node_ids = self.node_ids[i, :]
        #return obj
