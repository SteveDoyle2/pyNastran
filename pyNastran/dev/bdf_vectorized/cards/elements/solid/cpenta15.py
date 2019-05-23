from numpy import arange, searchsorted, array
from numpy.linalg import norm  # type: ignore

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import integer, integer_or_blank
from pyNastran.dev.bdf_vectorized.cards.elements.solid.cpenta6 import tri_area_centroid

from pyNastran.dev.bdf_vectorized.cards.elements.solid.solid_element import SolidElement

class CPENTA15(SolidElement):
    type = 'CPENTA15'
    nnodes = 15
    def __init__(self, model):
        """
        Defines the CPENTA15 object.

        Parameters
        ----------
        model : BDF
            the BDF object

        """
        SolidElement.__init__(self, model)

    def add_card(self, card, comment=''):
        i = self.i

        #comment = self._comments[i]
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
            integer(card, 7, 'node_id_5'),
            integer(card, 8, 'node_id_6'),
            integer_or_blank(card, 9, 'node_id_7', 0),
            integer_or_blank(card, 10, 'node_id_8', 0),
            integer_or_blank(card, 11, 'node_id_9', 0),
            integer_or_blank(card, 12, 'node_id_10', 0),
            integer_or_blank(card, 13, 'node_id_11', 0),
            integer_or_blank(card, 14, 'node_id_12', 0),
            integer_or_blank(card, 15, 'node_id_13', 0),
            integer_or_blank(card, 16, 'node_id_14', 0),
            integer_or_blank(card, 17, 'node_id_15', 0),
        ], dtype='int32')
        self.node_ids[i, :] = nids
        assert len(card) <= 18, 'len(CPENTA15 card) = %i\ncard=%s' % (len(card), card)
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
                self.node_ids[i, 10] = nid_map[nids[10]]
                self.node_ids[i, 11] = nid_map[nids[11]]
                self.node_ids[i, 12] = nid_map[nids[12]]
                self.node_ids[i, 13] = nid_map[nids[13]]
                self.node_ids[i, 14] = nid_map[nids[14]]

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

    def _area_centroid(self, element_id, xyz_cid0):
        (A1, c1) = tri_area_centroid(n1, n2, n3)
        (A2, c2) = tri_area_centroid(n4, n5, n6)
        return (A1, A2, c1, c2)

    def get_volume_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CPENTA15 elements.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the volume be summed

        .. note:: Volume for a CPENTA is the average area of two opposing faces
                  times the length between the centroids of those points
        """
        (A1, A2, c1, c2) = self._area_centroid(self, element_id, xyz_cid0)
        if total:
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        return volume

    def get_centroid_volume_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CPENTA15 elements.

        Parameters
        ----------
        element_id : (N, ) int ndarray; (default=None -> all)
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the volume be summed

        .. seealso:: CPENTA15.get_volume_by_element_id() and CPENTA15.get_centroid_by_element_id() for more information.
        """
        if element_id is None:
            element_id = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_node_index()

        (A1, A2, c1, c2) = self._area_centroid(element_id, xyz_cid0)
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
        Gets the centroid for one or more CPENTA15 elements.

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
        (A1, A2, c1, c2) = self._area_centroid(self, element_ids, xyz_cid0)

        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
        if total:
            centroid = centroid.mean()
        return centroid

    def get_face_nodes(self, nid, nid_opposite):
        raise NotImplementedError()
        #nids = self.node_ids[:4]
        #indx = nids.index(nid_opposite)
        #nids.pop(indx)
        #return nids

    def write_card(self, bdf_file, size=8, element_id=None):
        self.model.log.debug('CPENTA15.write_card; n=%s' % self.n)
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_id)
            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                if eid in self._comments:
                    bdf_file.write(self._comments[eid])
                n = [ni if ni != 0 else None for ni in n]
                card = ['CPENTA', eid, pid, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9],
                        n[10], n[11], n[12], n[13], n[14]]
                bdf_file.write(print_card_8(card).rstrip() + '\n')

    #def slice_by_index(self, i):
        #i = self._validate_slice(i)
        #obj = CPENTA15(self.model)
        #obj.n = len(i)
        ##obj._cards = self._cards[i]
        ##obj._comments = obj._comments[i]
        ##obj.comments = obj.comments[i]
        #obj.element_id = self.element_id[i]
        #obj.property_id = self.property_id[i]
        #obj.node_ids = self.node_ids[i, :]
        #return obj
