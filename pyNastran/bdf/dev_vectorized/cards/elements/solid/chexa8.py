from six.moves import zip, range
from numpy import zeros, arange, unique, dot, cross, abs, searchsorted, array, where, asarray
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank, fields)

from pyNastran.bdf.dev_vectorized.cards.elements.solid.solid_element import SolidElement


def quad_area_centroid(n1, n2, n3, n4):
    """
    Gets the area, :math:`A`, and centroid of a quad.::

      1-----2
      |   / |
      | /   |
      4-----3
    """
    a = n1 - n2
    b = n2 - n4
    area1 = 0.5 * norm(cross(a, b), axis=1)
    c1 = (n1 + n2 + n4) / 3.

    a = n2 - n4
    b = n2 - n3
    area2 = 0.5 * norm(cross(a, b), axis=1)
    c2 = (n2 + n3 + n4) / 3.

    area = area1 + area2
    try:
        centroid = (c1 * area1 + c2 * area2) / area
    except FloatingPointError:
        msg = '\nc1=%r\narea1=%r\n' % (c1, area1)
        msg += 'c2=%r\narea2=%r' % (c2, area2)
        raise FloatingPointError(msg)
    return(area, centroid)


class CHEXA8(SolidElement):
    type = 'CHEXA8'
    def __init__(self, model):
        """
        Defines the CHEXA8 object.

        :param self: the CHEXA8 object
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
            self.node_ids = zeros((ncards, 8), 'int32')

            #comments = {}
            for i, card in enumerate(cards):
                #comment = self._comments[i]
                #if comment:
                    #self._comments[eid] = comment

                #: Element ID
                self.element_id[i] = eid = integer(card, 1, 'element_id')
                #: Property ID
                self.property_id[i] = integer(card, 2, 'property_id')
                #: Node IDs
                nids = [
                    integer(card, 3, 'nid1'),
                    integer(card, 4, 'nid2'),
                    integer(card, 5, 'nid3'),
                    integer(card, 6, 'nid4'),
                    integer(card, 7, 'nid5'),
                    integer(card, 8, 'nid6'),
                    integer(card, 9, 'nid7'),
                    integer(card, 10, 'nid8')
                ]
                assert 0 not in nids, '%s\n%s' % (nids, card)
                self.node_ids[i, :] = nids
                assert len(card) == 11, 'len(CHEXA8 card) = %i' % len(card)

            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self._cards = []
            self._comments = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i,nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            c = self.centroid()
            v = self.volume()
            assert isinstance(v, float)
            for i in range(3):
                assert isinstance(c[i], float)

    def get_mass_matrix(self, i, model, positions, index0s):
        nnodes = 8
        ndof = 3 * nnodes
        pid = self.property_id[i]
        rho = self.model.elements.properties_solid.psolid.get_density(pid)[0]

        n0, n1, n2, n3, n4, n5, n6, n7 = self.node_ids[i, :]
        V = volume8(positions[self.node_ids[i, 0]],
                    positions[self.node_ids[i, 1]],
                    positions[self.node_ids[i, 2]],
                    positions[self.node_ids[i, 3]],

                    positions[self.node_ids[i, 4]],
                    positions[self.node_ids[i, 5]],
                    positions[self.node_ids[i, 6]],
                    positions[self.node_ids[i, 7]],
                    )

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
        dofs, nijv = self.get_dofs_nijv(index0s, n0, n1, n2, n3, n4, n5, n6, n7)
        return M, dofs, nijv

    def get_dofs_nijv(self, index0s, n0, n1, n2, n3, n4, n5, n6, n7):
        i0 = index0s[n0]
        i1 = index0s[n1]
        i2 = index0s[n2]
        i3 = index0s[n3]
        i4 = index0s[n4]
        i5 = index0s[n5]
        i6 = index0s[n6]
        i7 = index0s[n7]
        dofs = array([
            i0, i0+1, i0+2,
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,
            i3, i3+1, i3+2,
            i4, i4+1, i4+2,
            i5, i5+1, i5+2,
            i6, i6+1, i6+2,
            i7, i7+1, i7+2,
        ], 'int32')
        nijv = [
            # translation
            (n0, 1), (n0, 2), (n0, 3),
            (n1, 1), (n1, 2), (n1, 3),
            (n2, 1), (n2, 2), (n2, 3),
            (n3, 1), (n3, 2), (n3, 3),
            (n4, 1), (n4, 2), (n4, 3),
            (n5, 1), (n5, 2), (n5, 3),
            (n6, 1), (n6, 2), (n6, 3),
            (n7, 1), (n7, 2), (n7, 3),
        ]
        return dofs, nijv

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
        n5 = xyz_cid0[index_map(node_ids[i, 4]), :]
        n6 = xyz_cid0[index_map(node_ids[i, 5]), :]
        n7 = xyz_cid0[index_map(node_ids[i, 6]), :]
        n8 = xyz_cid0[index_map(node_ids[i, 7]), :]
        return n1, n2, n3, n4, n5, n6, n7, n8

    def get_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CHEXA8 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)

        ..note:: Volume for a CHEXA is the average area of two opposing faces
        times the length between the centroids of those points
        """
        n1, n2, n3, n4, n5, n6, n7, n8 = self._node_locations_element_id(ielement, xyz_cid0)
        (A1, c1) = quad_area_centroid(n1, n2, n3, n4)
        (A2, c2) = quad_area_centroid(n5, n6, n7, n8)
        if total:
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        return volume

    def get_centroid_volume(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CHEXA8 elements.

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed; centroid be averaged (default=False)

        ..see:: CHEXA8.get_volume() and CHEXA8.get_centroid() for more information.
        """
        n1, n2, n3, n4, n5, n6, n7, n8 = self._node_locations_element_id(element_id, xyz_cid0)
        (A1, c1) = quad_area_centroid(n1, n2, n3, n4)
        (A2, c2) = quad_area_centroid(n5, n6, n7, n8)
        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
        volume = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
        if total:
            centroid = centroid.mean()
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        assert volume.min() > 0.0, 'volume.min() = %f' % volume.min()
        return centroid, volume

    def get_centroid(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid for one or more CHEXA8 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be averaged (default=False)
        """
        n1, n2, n3, n4, n5, n6, n7, n8 = self._node_locations_element_id(element_id, xyz_cid0)
        (A1, c1) = quad_area_centroid(n1, n2, n3, n4)
        (A2, c2) = quad_area_centroid(n5, n6, n7, n8)
        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
        if total:
            centroid = centroid.mean(axis=0)
        return centroid

    def get_face_nodes(self, nid, nid_opposite):
        raise NotImplementedError()
        nids = self.nodeIDs()[:8]
        indx = nids.index(nid_opposite)
        nids.pop(indx)
        return nids

    def write_bdf(self, f, size=8, element_id=None):
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_id)

            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i]):
                card = ['CHEXA', eid, pid, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7]]
                f.write(print_card_8(card))


def volume8(n1, n2, n3, n4, n5, n6, n7, n8):
    (A1, c1) = quad_area_centroid(n1, n2, n3, n4)
    (A2, c2) = quad_area_centroid(n5, n6, n7, n8)
    volume = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
    return volume
