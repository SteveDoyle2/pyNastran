from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)

class CROD(object):
    type = 'CROD'
    def __init__(self, model):
        """
        Defines the CROD object.

        :param self: the CROD object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._crod = []
        self._crod_comments = []

    def add(self, card, comment=None):
        self._crod.append(card)
        self._crod_comment.append(comment)

    def build(self):
        """
        :param self: the CROD object
        :param cards: the list of PCOMP cards
        """
        cards = self._crod
        ncards = len(cards)
        self.n = ncards
        if ncards:
            #: Property ID
            self.eid = zeros(ncards, 'int32')
            self.pid = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 2), 'int32')

            for i, card in enumerate(cards):
                self.eid[i] = integer(card, 1, 'eid')
                self.pid[i] = integer_or_blank(card, 2, 'pid', self.eid)
                self.node_ids = [integer(card, 3, 'n1'),
                                 integer(card, 4, 'n2')]
                assert len(card) == 5, 'len(CROD card) = %i' % len(card)

            i = self.eid.argsort()
            self.eid = self.eid[i]
            self.pid = self.pid[i]
            self.node_ids = self.node_ids[i, :]

            unique_eids = unique(self.eid)
            if len(unique_eids) != len(self.eid):
                raise RuntimeError('There are duplicate CROD IDs...')
            self._crod = []
            self._crod_comments = []
        
    def mass(self, total=False):
        """
        mass = rho * A * L + nsm
        """
        grid_cid0 = self.grid.position()
        p1 = grid_cid0[self.node_ids[:, 0]]
        p2 = grid_cid0[self.node_ids[:, 1]]
        L = p2 - p1
        i = self.model.prod.get_index(pid)
        A = self.model.prod.A[i]
        mid = self.model.prod.mid[i]

        rho, E, J = self.model.Materials.get_rho_E_J(self.mid)
        rho = self.model.Materials.get_rho(self.mid)
        E   = self.model.Materials.get_E(self.mid)
        J   = self.model.Materials.get_J(self.mid)

        mass = norm(L, axis=1) * A * rho + self.nsm
        if total:
            return mass.sum()
        else:
            return mass
        
    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('CROD', self.n))
        return msg

    def write_bdf(self, f, size=8, pids=None):
        for (eid, pid, n12) in (self.eid, self.pid, self.node_ids):
            card = ['CROD', eid, pid, self.n12[0], n12[1] ]
            f.write(print_card(card))