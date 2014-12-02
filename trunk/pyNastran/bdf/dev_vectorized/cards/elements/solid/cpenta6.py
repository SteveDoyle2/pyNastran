from six.moves import zip
from numpy import zeros, arange, unique, dot, cross, abs, searchsorted, array, where, asarray
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank, fields)

from pyNastran.bdf.dev_vectorized.cards.elements.solid.solid_element import SolidElement

def tri_area_centroid(n1, n2, n3):
    """
    Gets the area, :math:`A`, and centroid of a triangle.::

      1-----2
      |   /
      | /
      3
    """
    a = n1 - n2
    b = n2 - n3
    area = 0.5 * norm(cross(a, b), axis=1)
    centroid = (n1 + n2 + n3) / 3.

    return(area, centroid)


class CPENTA6(SolidElement):
    type = 'CPENTA6'
    nnodes = 6
    def __init__(self, model):
        """
        Defines the CPENTA object.

        :param self: the CPENTA object
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
            integer(card, 7, 'node_id_5'),
            integer(card, 8, 'node_id_6'),
        ], dtype='int32')
        assert 0 not in nids, '%s\n%s' % (nids, card)
        self.node_ids[i, :] = nids
        assert len(card) == 9, 'len(CPENTA6 card) = %i' % len(card)
        self.i += 1

    def get_mass_matrix(self, i, model, positions, index0s):
        return M, dofs, nijv

    def get_stiffness_matrix(self, i, model, positions, index0s):
        return K, dofs, nijv

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
        n5 = xyz_cid0[get_node_index_by_node_id(node_ids[i, 4], msg), :]
        n6 = xyz_cid0[get_node_index_by_node_id(node_ids[i, 5], msg), :]
        return n1, n2, n3, n4, n5, n6

    def get_volume_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more elements.

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)

        ..note:: Volume for a CPENTA is the average area of two opposing faces
        times the length between the centroids of those points
        """
        n1, n2, n3, n4, n5, n6 = self._get_node_locations_by_element_id(element_id, xyz_cid0)
        (A1, c1) = tri_area_centroid(n1, n2, n3)
        (A2, c2) = tri_area_centroid(n4, n5, n6)
        volume = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
        if total:
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        return volume

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

    def get_centroid_volume_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more elements.

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed; centroid be averaged (default=False)

        ..see:: CPENTA6.get_volume_by_element_id() and CPENTA6.get_centroid_by_element_id
        """
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_from_node_index()

        n1, n2, n3, n4, n5, n6 = self._get_node_locations_by_element_id(element_id, xyz_cid0)
        (A1, c1) = tri_area_centroid(n1, n2, n3)
        (A2, c2) = tri_area_centroid(n4, n5, n6)
        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
        volume = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
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
        if element_id is None:
            element_id = self.element_id
        n1, n2, n3, n4, n5, n6 = self._get_node_locations_by_element_id(element_id, xyz_cid0)
        (A1, c1) = tri_area_centroid(n1, n2, n3)
        (A2, c2) = tri_area_centroid(n4, n5, n6)
        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
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
                card = ['CPENTA', eid, pid, n[0], n[1], n[2], n[3], n[4], n[5]]
                f.write(print_card_8(card))

