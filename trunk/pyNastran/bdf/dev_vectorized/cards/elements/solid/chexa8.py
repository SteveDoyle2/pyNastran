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
    op2_id = 64
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

    def _node_locations(self, xyz_cid0):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_positions()

        n1 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 0]), :]
        n2 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 1]), :]
        n3 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 2]), :]
        n4 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 3]), :]
        n5 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 4]), :]
        n6 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 5]), :]
        n7 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 6]), :]
        n8 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 7]), :]
        return n1, n2, n3, n4, n5, n6, n7, n8

    def _get_area_centroid(self, element_ids, xyz_cid0):
        n1, n2, n3, n4, n5, n6, n7, n8 = self._node_locations(xyz_cid0)
        (A1, c1) = quad_area_centroid(n1, n2, n3, n4)
        (A2, c2) = quad_area_centroid(n5, n6, n7, n8)
        return (A1, A2, c1, c2)

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

    def get_centroid_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CHEXA8 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed; centroid be averaged (default=False)

        ..see:: CHEXA8.get_volume() and CHEXA8.get_centroid() for more information.
        """
        if element_ids is None:
            element_ids = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_positions()

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
        Gets the centroid for one or more CHEXA8 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be averaged (default=False)
        """
        if element_ids is None:
            element_ids = self.element_id
        (A1, A2, c1, c2) = self._get_area_centroid(element_ids, xyz_cid0)
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

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_ids)

            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i]):
                card = ['CHEXA', eid, pid, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7]]
                f.write(print_card_8(card))

def get_hex_volume(self, element_ids, xyz_cid0):
    (A1, A2, c1, c2) = self._get_area_centroid(element_ids, xyz_cid0)
    volume = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
    return volume
