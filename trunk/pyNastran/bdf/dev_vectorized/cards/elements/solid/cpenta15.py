import cStringIO
from numpy import zeros, arange, dot, cross, searchsorted, array, asarray
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (fields, integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)
from .cpenta6 import area_centroid


class CPENTA15(object):
    type = 'CPENTA15'
    op2_id = 63
    def __init__(self, model):
        """
        Defines the CPENTA15 object.

        :param self: the CPENTA15 object
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
            self.node_ids = zeros((ncards, 15), 'int32')

            comments = {}
            for i, card in enumerate(cards):
                comment = self._comments[i]
                element_id = integer(card, 1, 'element_id')
                #if comment:
                    #self.comments[eid] = comment

                #: Element ID
                self.element_id[i] = eid
                #: Property ID
                self.property_id[i] = integer(card, 2, 'property_id')
                #: Node IDs
                self.node_ids[i, :] = fields(integer, card, 'node_ids', i=3, j=13)
                assert len(card) == 13, 'len(CPENTA15 card) = %i' % len(card)

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

    def _area_centroid(self, element_ids, xyz_cid0):
        (A1, c1) = area_centroid(n1, n2, n3)
        (A2, c2) = area_centroid(n4, n5, n6)
        return (A1, A2, c1, c2)

    def get_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CPENTA15 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)

        ..note:: Volume for a CPENTA is the average area of two opposing faces
        times the length between the centroids of those points
        """
        (A1, A2, c1, c2) = self._area_centroid(self, element_ids, xyz_cid0)
        if total:
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        return volume

    def get_centroid_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CPENTA15 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)

        ..see:: CPENTA15.get_volume() and CPENTA15.get_centroid() for more information.
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
        Gets the centroid for one or more CPENTA15 elements.

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

    def get_mass(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the mass for one or more CPENTA15 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if element_ids is None:
            element_ids = self.element_id
        V = self.get_volume(element_ids, xyz_cid0)

        mid = self.model.properties_solid.get_mid(self.pid)
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
                card = ['CHEXA', eid, pid, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9],
                        n[10], n[11], n[12], n[13], n[14]]
                f.write(print_card_8(card))

    def __getitem__(self, element_ids):
        i = searchsorted(self.element_id, element_ids)
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        i = asarray(i)
        obj = CPENTA15(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.property_id = self.property_id[i]
        obj.node_ids = self.node_ids[i, :]
        return obj

    def __repr__(self):
        f = cStringIO.StringIO()
        f.write('<CPENTA15 object> n=%s\n' % self.n)
        self.write_bdf(f)
        return f.getvalue()