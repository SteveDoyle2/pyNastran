from numpy import dot, cross

def volume4(n1, n2, n3, n4):
    r"""
    Gets the volume, :math:`V`, of the tetrahedron.
    
    .. math:: V = \frac{(a-d) \cdot \left( (b-d) \times (c-d) \right) }{6}
    """
    V = -dot((n1 - n4), cross(n2 - n4, n3 - n4)) / 6.
    return V


class CTETRA(object):
    type = 'CTETRA'
    def __init__(self, model):
        """
        Defines the CTETRA object.

        :param self: the CTETRA object
        :param model: the BDF object
        """
        self.model = model
        self._ctetra = []
        self._ctetra_comment = []

    def add(self, card, comment):
        self._ctetra.append(card)
        self._ctetra_comment.append(comment)

    def build(self):
        cards = self._ctetra
        ncards = len(cards)

        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            self.element_id = zeros(ncards, 'int32')
            self.property_id = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 4), 'int32')
            
            comments = {}
            for i, card in enumerate(cards):
                comment = self._ctetra_comments[i]
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

    def volume(self, eids=None):
        if eids is None:
            eids = self.eids
        n = len(eids)
        n1 = nodes[self.node_ids[:, 0]]
        n2 = nodes[self.node_ids[:, 1]]
        n3 = nodes[self.node_ids[:, 2]]
        n4 = nodes[self.node_ids[:, 3]]
        
        V = zeros(n, self.model.float)
        
        i = 0
        for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
            V[i] = volume4(n1i, n2i, n3i, n4i)
            i += 1
        return V

    def centroid(self, eids=None, total=False):
        if eids is None:
            eids = self.eids
        n1 = nodes[self.node_ids[:, 0]]
        n2 = nodes[self.node_ids[:, 1]]
        n3 = nodes[self.node_ids[:, 2]]
        n4 = nodes[self.node_ids[:, 3]]
        centroid = (n1 + n2 + n3 + n4) / 4.0
        if total:
            return centroid.mean()
        else:
            return centroid

    def mass(self, eids=None, total=False):
        if eids is None:
            eids = self.eids
        V = self.volume(eids)

        mid = self.model.properties_solid.get_mid(self.pid)
        rho = self.model.materials.get_rho(mid)
        nsm = self.model.materials.get_nsm(mid)
        rho, nsm = self.model.materials.get_rho_nsm(mid)
        
        mass = V * rho + nsm
        if total:
            return mass.sum()
        else:
            return mass

    def get_face_nodes(self, nid, nidOpposite):
        asdf
        nids = self.nodeIDs()[:4]
        indx = nids.index(nidOpposite)
        nids.pop(indx)
        return nids

    def write_bdf(self, f, size=8, eids=None):
        if eids is None:
            for (eid, pid, n1234) in zip(self.element_id, self.property_id, self.node_ids)
                card = ['CTETRA', eid, pid, n1234[0], n1234[1], n1234[2], n1234[3]]
                f.write(print_card(card))
        else:
            i = searchsorted(self.eids, eids)
            for (eid, pid, n1234) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :])
                card = ['CTETRA', eid, pid, n1234[0], n1234[1], n1234[2], n1234[3]]
                f.write(print_card(card))
            