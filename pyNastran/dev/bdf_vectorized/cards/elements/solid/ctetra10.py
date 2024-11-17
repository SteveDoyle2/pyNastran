from numpy import zeros, arange, cross, searchsorted, array

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import integer, integer_or_blank
from pyNastran.dev.bdf_vectorized.cards.elements.solid.solid_element import SolidElement

def volume4(n1, n2, n3, n4):
    r"""
    Gets the volume, :math:`V`, of the tetrahedron.

    .. math:: V = \frac{(a-d) \cdot \left( (b-d) \times (c-d) \right) }{6}
    """
    V = (n1 - n4) @ cross(n2 - n4, n3 - n4) / -6.
    return V


class CTETRA10(SolidElement):
    type = 'CTETRA10'
    nnodes = 10
    def __init__(self, model):
        """
        Defines the CTETRA10 object.

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
        #comment = self._comments[i]
        element_id = integer(card, 1, 'element_id')
        if comment:
            self.set_comment(eid, comment)

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
        assert len(card) <= 13, 'len(CTETRA10 card) = %i\ncard=%s' % (len(card), card)
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
                self.node_ids[i, 4] = nid_map[nids[4]]
                self.node_ids[i, 5] = nid_map[nids[5]]
                self.node_ids[i, 6] = nid_map[nids[6]]
                self.node_ids[i, 7] = nid_map[nids[7]]
                self.node_ids[i, 8] = nid_map[nids[8]]
                self.node_ids[i, 9] = nid_map[nids[9]]

    def _verify(self, xref=True):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' % (i, nid)
        if xref:
            c = self.centroid()
            v = self.volume()
            assert isinstance(v, float)
            for i in range(3):
                assert isinstance(c[i], float)

    def _node_locations(self, xyz_cid0, n):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_node_index()

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

        Parameters
        ----------
        element_id : (N, ) int ndarray; (default=None -> all)
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the volume be summed
        """
        if element_id is None:
            element_id = self.element_id
        n1, n2, n3, n4 = self._node_locations(xyz_cid0, 4)

        n = len(element_id)
        V = zeros(n, self.model.float_fmt)

        i = 0
        for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
            V[i] = volume4(n1i, n2i, n3i, n4i)
            i += 1
        if total:
            V = V.sum()
        return V

    def get_centroid_volume(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CTETRA10 elements.

        Parameters
        ----------
        element_id : (N, ) int ndarray; (default=None -> all)
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the volume be summed

        .. seealso:: CTETRA10.volume() and CTETRA10.centroid() for more information.
        """
        if element_id is None:
            element_id = self.element_id

        n1, n2, n3, n4 = self._node_locations(xyz_cid0, 4)
        n = len(element_id)
        volume = zeros(n, self.model.float_fmt)

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

        Parameters
        ----------
        element_id : (N, ) int ndarray; (default=None -> all)
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the centroid be summed
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

        Parameters
        ----------
        element_id : (N, ) int ndarray; (default=None -> all)
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

    def get_face_nodes(self, nid, nid_opposite):
        raise NotImplementedError()
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
            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                if eid in self._comments:
                    bdf_file.write(self._comments[eid])
                n = [ni if ni != 0 else None for ni in n]
                card = ['CTETRA', eid, pid, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9]]
                bdf_file.write(print_card_8(card))

    #def slice_by_index(self, i):
        #i = self._validate_slice(i)
        #obj = CTETRA10(self.model)
        #obj.n = len(i)
        ##obj._cards = self._cards[i]
        ##obj._comments = obj._comments[i]
        ##obj.comments = obj.comments[i]
        #obj.element_id = self.element_id[i]
        #obj.property_id = self.property_id[i]
        #obj.node_ids = self.node_ids[i, :]
        #return obj
