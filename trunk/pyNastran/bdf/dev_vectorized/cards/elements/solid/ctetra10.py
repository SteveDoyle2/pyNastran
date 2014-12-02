from six.moves import zip
from numpy import zeros, arange, dot, cross, searchsorted, array, asarray
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16

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


class CTETRA10(SolidElement):
    type = 'CTETRA10'
    nnodes = 10
    def __init__(self, model):
        """
        Defines the CTETRA10 object.

        :param self: the CTETRA10 object
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
        #comment = self._comments[i]
        element_id = integer(card, 1, 'element_id')
        #if comment:
            #self.comments[eid] = comment

        #: Element ID
        self.element_id[i] = element_id
        #: Property ID
        self.property_id[i] = integer(card, 2, 'property_id')
        #: Node IDs
        nids = array([
            integer(card, 3, 'node_id_1'),
            integer(card, 4, 'node_id_2'),
            integer(card, 5, 'node_id_3'),
            integer(card, 6, 'node_id_4'),
            integer_or_blank(card, 7, 'node_id_5', 0),
            integer_or_blank(card, 8, 'node_id_6', 0),
            integer_or_blank(card, 9, 'node_id_7', 0),
            integer_or_blank(card, 10, 'node_id_8', 0),
            integer_or_blank(card, 11, 'node_id_9', 0),
            integer_or_blank(card, 12, 'node_id_10', 0),
        ], dtype='int32')
        self.node_ids[i, :] = nids
        assert len(card) <= 13, 'len(CTETRA10 card) = %i' % len(card)
        self.i += 1

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

    def _node_locations(self, xyz_cid0, n):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_index()

        n1 = xyz_cid0[self.model.grid.get_node_index_by_node_id(self.node_ids[:, 0]), :]
        n2 = xyz_cid0[self.model.grid.get_node_index_by_node_id(self.node_ids[:, 1]), :]
        n3 = xyz_cid0[self.model.grid.get_node_index_by_node_id(self.node_ids[:, 2]), :]
        n4 = xyz_cid0[self.model.grid.get_node_index_by_node_id(self.node_ids[:, 3]), :]
        if n == 4:
            return n1, n2, n3, n4
        assert n == 10, n

        n5 = xyz_cid0[self.model.grid.get_node_index_by_node_id(self.node_ids[:, 4]), :]
        n6 = xyz_cid0[self.model.grid.get_node_index_by_node_id(self.node_ids[:, 5]), :]
        n7 = xyz_cid0[self.model.grid.get_node_index_by_node_id(self.node_ids[:, 6]), :]
        n8 = xyz_cid0[self.model.grid.get_node_index_by_node_id(self.node_ids[:, 7]), :]
        n9 = xyz_cid0[self.model.grid.get_node_index_by_node_id(self.node_ids[:, 8]), :]
        n10 = xyz_cid0[self.model.grid.get_node_index_by_node_id(self.node_ids[:, 9]), :]
        return n1, n2, n3, n4, n5, n6, n7, n8, n9, n10

    def get_volume_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CTETRA10 elements.

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)
        """
        if element_id is None:
            element_id = self.element_id
        n1, n2, n3, n4 = self._node_locations(xyz_cid0, 4)

        n = len(element_id)
        V = zeros(n, self.model.float)

        i = 0
        for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
            V[i] = volume4(n1i, n2i, n3i, n4i)
            i += 1
        return V

    def get_centroid_volume(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CTETRA10 elements.

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)

        ..see:: CTETRA10.volume() and CTETRA10.centroid() for more information.
        """
        if element_id is None:
            element_id = self.element_id

        n1, n2, n3, n4 = self._node_locations(xyz_cid0, 4)
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

    def get_centroid_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid for one or more CTETRA elements.

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if element_id is None:
            element_id = self.element_id

        n1, n2, n3, n4 = self._node_locations(xyz_cid0, 4)
        centroid = (n1 + n2 + n3 + n4) / 4.0
        if total:
            centroid = centroid.mean()
        return centroid

    def get_mass_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the mass for one or more CTETRA10 elements.

        :param element_id: the elements to consider (default=None -> all)
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
                n = [ni if ni != 0 else None for ni in n]
                card = ['CTETRA', eid, pid, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9]]
                f.write(print_card_8(card))

    #def slice_by_index(self, i):
        #i = asarray(i)
        #obj = CTETRA10(self.model)
        #obj.n = len(i)
        ##obj._cards = self._cards[i]
        ##obj._comments = obj._comments[i]
        ##obj.comments = obj.comments[i]
        #obj.element_id = self.element_id[i]
        #obj.property_id = self.property_id[i]
        #obj.node_ids = self.node_ids[i, :]
        #return obj
