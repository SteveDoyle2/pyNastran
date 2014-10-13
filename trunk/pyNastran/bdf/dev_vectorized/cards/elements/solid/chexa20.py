from numpy import zeros, arange, dot, cross, searchsorted, array, asarray
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (fields, integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)

from .chexa8 import quad_area_centroid, get_hex_volume
from .solid_element import SolidElement

class CHEXA20(SolidElement):
    type = 'CHEXA20'
    op2_id = 65
    def __init__(self, model):
        """
        Defines the CHEXA20 object.

        :param self: the CHEXA20 object
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
            self.node_ids = zeros((ncards, 20), 'int32')

            comments = {}
            for i, card in enumerate(cards):
                comment = self._comments[i]
                if comment:
                    self.comments[eid] = comment

                #: Element ID
                self.element_id[i] = integer(card, 1, 'element_id')
                #: Property ID
                self.property_id[i] = integer(card, 2, 'property_id')
                #: Node IDs
                self.node_ids[i, :] = fields(integer, card, 'node_ids', i=3, j=23)
                assert len(card) == 23, 'len(CHEXA20 card) = %i' % len(card)

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

    def _get_area_centroid(self, element_ids, xyz_cid0):
        n1, n2, n3, n4, n5, n6, n7, n8 = self._node_locations(xyz_cid0)
        (A1, c1) = quad_area_centroid(n1, n2, n3, n4)
        (A2, c2) = quad_area_centroid(n5, n6, n7, n8)
        return (A1, A2, c1, c2)

    def get_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CHEXA20 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)

        ..note:: Volume for a CHEXA is the average area of two opposing faces
        times the length between the centroids of those points
        """
        volume = get_hex_volume(self, element_ids, xyz_cid0)
        if total:
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        return volume

    def get_centroid_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CHEXA20 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)

        ..see:: CHEXA20.get_volume() and CHEXA20.get_centroid() for more information.
        """
        if element_ids is None:
            element_ids = self.element_id

        (A1, A2, c1, c2) = self._area_centroid(element_ids, xyz_cid0)
        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
        volume = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
        if total:
            centroid = centroid.mean()
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        assert volume.min() > 0.0, 'volume.min() = %f' % volume.min()
        return centroid, volume

    def get_centroid(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid for one or more CHEXA20 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if element_ids is None:
            element_ids = self.element_id
        (A1, A2, c1, c2) = self._area_centroid(self, element_ids, xyz_cid0)
        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
        if total:
            centroid = centroid.mean()
        return centroid

    def get_face_nodes(self, nid, nid_opposite):
        raise NotImplementedError()
        nids = self.nodeIDs()[:4]
        indx = nids.index(nid_opposite)
        nids.pop(indx)
        return nids

    def get_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CHEXA8 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)

        ..note:: Volume for a CHEXA is the average area of two opposing faces
        times the length between the centroids of those points
        """
        volume = get_hex_volume(self, element_ids, xyz_cid0)
        if total:
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        return volume

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_ids)
            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                card = ['CHEXA', eid, pid, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9],
                        n[10], n[11], n[12], n[13], n[14], n[15], n[16], n[17], n[18], n[19]]
                f.write(print_card_8(card))

    #def slice_by_index(self, i):
        #i = asarray(i)
        #obj = CHEXA20(self.model)
        #obj.n = len(i)
        ##obj._cards = self._cards[i]
        ##obj._comments = obj._comments[i]
        ##obj.comments = obj.comments[i]
        #obj.element_id = self.element_id[i]
        #obj.property_id = self.property_id[i]
        #obj.node_ids = self.node_ids[i, :]
        #return obj
