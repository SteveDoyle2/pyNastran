from itertools import izip

from numpy import zeros, array, where, unique, searchsorted

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank)


class POINTAX(object):
    type = 'POINTAX'
    def __init__(self, model):
        """
        Defines the POINTAX object.

        :param self: the POINTAX object
        :param model: the BDF object

        +---------+-----+-----+-----+
        |    1    |  2  |  3  |  4  |
        +=========+=====+=====+=====+
        | POINTAX | NID | RID | PHI |
        +---------+-----+-----+-----+
        """
        self.model = model
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
            self.node_id = zeros(ncards, 'int32')
            self.phi = zeros(ncards, float_fmt)
            self.ring_id = zeros(ncards, 'int32')

            for i, card in enumerate(cards):
                #: Node ID
                self.node_id[i] = integer(card, 1, 'nid')
                self.ring_id[i] = integer(card, 2, 'rid')
                self.phi[i] = double_or_blank(card, 3, 'phi', 0.)

    def positions(self, node_ids=None):
        asdf
        if node_ids is None:
            node_ids = self.node_id
        xyz = xyz.copy()

        n = arange(self.n)
        i = where(self.coord_id != 0)[0]
        if i:
            n = n[i]
            cids = set(list(unique(self.coord_id)))
            for cid in cids:
                i = where(self.coord_id != 0)[0]
                T = self.model.coord.transform(cid)
                xyzi = xyz[n[i], :]
                xyzi = dot(transpose(T), dot(xyzi, T))
        return xyz

    def positions_wrt(self, node_ids=None, cids=None):
        raise NotImplementedError()

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('POINTAX', self.n))
        return msg

    def write_bdf(self, f, size=8):
        if self.n:
            f.write('$POINTAX\n')
            for (nid, rid, phi) in izip(self.point_id, self.ring_id, self.phi):
                card = ['POINTAX', rid, phi]
                f.write(print_card(card, size))

    def __repr__(self):
        msg = "<POINTAX>\n"
        msg += '  nPOINTAX = %i' % self.n