from numpy import zeros, dot, cross, searchsorted
from pyNastran.utils.mathematics import norm_axis as norm

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


class CTETRA4(object):
    type = 'CTETRA4'
    def __init__(self, model):
        """
        Defines the CTETRA4 object.

        :param self: the CTETRA4 object
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
            self.node_ids = zeros((ncards, 4), 'int32')
            
            comments = {}
            for i, card in enumerate(cards):
                comment = self._comments[i]
                eid = integer(card, 1, 'eid')
                if comment:
                    self.comments[eid] = comment
                    
                #: Element ID
                self.element_id[i] = eid
                #: Property ID
                self.property_id[i] = integer(card, 2, 'pid')
                #: Node IDs
                self.node_ids[i, :] = fields(integer, card, 'nid', i=3, j=7)
                assert len(card) == 7, 'len(CTETRA4 card) = %i' % len(card)
                
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

    def volume(self, eids=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CHEXA elements.
        
        :param eids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)
        """
        if eids is None:
            eids = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.positions()

        n1 = xyz_cid0[self.node_ids[:, 0], :]
        n2 = xyz_cid0[self.node_ids[:, 1], :]
        n3 = xyz_cid0[self.node_ids[:, 2], :]
        n4 = xyz_cid0[self.node_ids[:, 3], :]
        
        n = len(eids)
        V = zeros(n, self.model.float)
        
        i = 0
        for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
            V[i] = volume4(n1i, n2i, n3i, n4i)
            i += 1
        return V
        
    def centroid_volume(self, eids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CHEXA elements.
        
        :param eids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)
        
        ..see:: CHEXA.volume() and CHEXA.centroid for more information.
        """
        if eids is None:
            eids = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.positions()
        
        n1 = xyz_cid0[self.node_ids[:, 0], :]
        n2 = xyz_cid0[self.node_ids[:, 1], :]
        n3 = xyz_cid0[self.node_ids[:, 2], :]
        n4 = xyz_cid0[self.node_ids[:, 3], :]
        
        n = len(eids)
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

    def centroid(self, eids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid for one or more CHEXA elements.
        
        :param eids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if eids is None:
            eids = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.positions()

        n1 = xyz_cid0[self.node_ids[:, 0], :]
        n2 = xyz_cid0[self.node_ids[:, 1], :]
        n3 = xyz_cid0[self.node_ids[:, 2], :]
        n4 = xyz_cid0[self.node_ids[:, 3], :]
        centroid = (n1 + n2 + n3 + n4) / 4.0
        if total:
            centroid = centroid.mean()
        return centroid

    def mass(self, eids=None, xyz_cid0=None, total=False):
        """
        Gets the mass for one or more CHEXA elements.
        
        :param eids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if eids is None:
            eids = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.positions()

        V = self.volume(eids, xyz_cid0)

        mid = self.model.properties_solid.get_mid(self.pid)
        rho = self.model.materials.get_rho(mid)
        nsm = self.model.materials.get_nsm(mid)
        rho, nsm = self.model.materials.get_rho_nsm(mid)
        
        mass = V * rho + nsm
        if total:
            mass = mass.sum()
        return mass

    def get_face_nodes(self, nid, nid_opposite):
        raise NotImplementedError()
        nids = self.nodeIDs()[:4]
        indx = nids.index(nid_opposite)
        nids.pop(indx)
        return nids

    def write_bdf(self, f, size=8, eids=None):
        if eids is None:
            for (eid, pid, n) in zip(self.element_id, self.property_id, self.node_ids):
                card = ['CTETRA', eid, pid, n[0], n[1], n[2], n[3]]
                f.write(print_card(card))
        else:
            i = searchsorted(self.element_id, eids)
            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                card = ['CTETRA', eid, pid, n[0], n[1], n[2], n[3]]
                f.write(print_card(card))
            