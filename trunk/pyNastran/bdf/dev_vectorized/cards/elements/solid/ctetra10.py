from numpy import zeros, arange, dot, cross, searchsorted
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (fields, integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)

def volume4(n1, n2, n3, n4):
    r"""
    Gets the volume, :math:`V`, of the tetrahedron.

    .. math:: V = \frac{(a-d) \cdot \left( (b-d) \times (c-d) \right) }{6}
    """
    V = -dot((n1 - n4), cross(n2 - n4, n3 - n4)) / 6.
    return V


class CTETRA10(object):
    type = 'CTETRA10'
    def __init__(self, model):
        """
        Defines the CTETRA10 object.

        :param self: the CTETRA10 object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        ncards = len(cards)

        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            self.element_id = zeros(ncards, 'int32')
            self.property_id = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 10), 'int32')

            comments = {}
            for i, card in enumerate(cards):
                comment = self._comments[i]
                element_id = integer(card, 1, 'element_id')
                #if comment:
                    #self.comments[eid] = comment

                #: Element ID
                self.element_id[i] = element_id
                #: Property ID
                self.property_id[i] = integer(card, 2, 'property_id')
                #: Node IDs
                self.node_ids[i, :] = fields(integer, card, 'node_ids', i=3, j=13)
                assert len(card) == 13, 'len(CTETRA10 card) = %i' % len(card)

            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self._cards = []
            self._comments = []

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

    def _node_locations(self, xyz_cid0, n=4):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_positions()
        n1 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 0]), :]
        n2 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 1]), :]
        n3 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 2]), :]
        n4 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 3]), :]
        if n4:
            return n1, n2, n3, n4

        n5 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 4]), :]
        n6 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 5]), :]
        n7 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 6]), :]
        n8 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 7]), :]
        n9 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 8]), :]
        n10 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 9]), :]
        return n1, n2, n3, n4, n5, n6, n7, n8, n9, n10

    def get_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CTETRA10 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)
        """
        if element_ids is None:
            element_ids = self.element_id
        n1, n2, n3, n4 = self._node_locations(xyz_cid0)

        n = len(element_ids)
        V = zeros(n, self.model.float)

        i = 0
        for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
            V[i] = volume4(n1i, n2i, n3i, n4i)
            i += 1
        return V

    def get_centroid_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CTETRA10 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)

        ..see:: CTETRA10.volume() and CTETRA10.centroid() for more information.
        """
        if element_ids is None:
            element_ids = self.element_id

        n1, n2, n3, n4 = self._node_locations(xyz_cid0)
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

    def get_centroid(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid for one or more CTETRA elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if element_ids is None:
            element_ids = self.element_id

        n1, n2, n3, n4 = self._node_locations(xyz_cid0)
        centroid = (n1 + n2 + n3 + n4) / 4.0
        if total:
            centroid = centroid.mean()
        return centroid

    def get_mass(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the mass for one or more CTETRA10 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if element_ids is None:
            element_ids = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_positions()

        V = self.get_volume(element_ids, xyz_cid0)
        mid = self.model.properties_solid.get_mid(self.property_id)
        rho = self.model.materials.get_rho(mid)

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

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_ids)
            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                card = ['CTETRA', eid, pid, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9]]
                f.write(print_card(card))
