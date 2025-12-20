from numpy import arange, searchsorted, array
from numpy.linalg import norm  # type: ignore

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import integer, integer_or_blank
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

from pyNastran.dev.bdf_vectorized.cards.elements.solid.cpyram5 import tri_area_centroid

from pyNastran.dev.bdf_vectorized.cards.elements.solid.solid_element import SolidElement

class CPYRAM13(SolidElement):
    type = 'CPYRAM13'
    nnodes = 13
    def __init__(self, model):
        """
        Defines the CPYRAM13 object.

        Parameters
        ----------
        model : BDF
            the BDF

        """
        SolidElement.__init__(self, model)

    def add_card(self, card, comment: str=''):
        self.model.log.debug('add...')
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
        ], dtype='int32')
        self.node_ids[i, :] = nids
        assert len(card) <= 17, 'len(CPYRAM13 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

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
        return A1, A2, c1, c2

    def get_volume_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CPYRAM13 elements.

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

    def get_centroid_volume_by_element_id(self, element_id=None, xyz_cid0=None,
                                          total: bool=False):
        """
        Gets the centroid and volume for one or more CPYRAM13 elements.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the volume be summed

        .. seealso:: CPYRAM13.get_volume_by_element_id() and CPYRAM13.get_centroid_by_element_id() for more information.
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

    def get_centroid_by_element_id(self, element_id=None, xyz_cid0=None,
                                   total: bool=False):
        """
        Gets the centroid for one or more CPYRAM13 elements.

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
        self.model.log.debug('CPYRAM13.write_card; n=%s' % self.n)
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_id)
            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                if eid in self._comments:
                    bdf_file.write(self._comments[eid])
                n = [ni if ni != 0 else None for ni in n]
                card = ['CPYRAM', eid, pid,
                        n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9],
                        n[10], n[11], n[12]]
                bdf_file.write(print_card_8(card).rstrip() + '\n')

    #def slice_by_index(self, i):
        #i = self._validate_slice(i)
        #obj = CPYRAM13(self.model)
        #obj.n = len(i)
        ##obj._cards = self._cards[i]
        ##obj._comments = obj._comments[i]
        ##obj.comments = obj.comments[i]
        #obj.element_id = self.element_id[i]
        #obj.property_id = self.property_id[i]
        #obj.node_ids = self.node_ids[i, :]
        #return obj
