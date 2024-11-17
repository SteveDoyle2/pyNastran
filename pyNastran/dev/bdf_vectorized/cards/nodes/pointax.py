from numpy import zeros, where, unique, transpose, dot, arange

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    double_or_blank)
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard


class POINTAX(VectorizedCard):
    type = 'POINTAX'
    def __init__(self, model):
        """
        Defines the POINTAX object.

        Parameters
        ----------
        model : BDF
           the BDF object

        +---------+-----+-----+-----+
        |    1    |  2  |  3  |  4  |
        +=========+=====+=====+=====+
        | POINTAX | NID | RID | PHI |
        +---------+-----+-----+-----+
        """
        VectorizedCard.__init__(self, model)
        self._cards = []
        self._comments = []

    def add_card(self, card: BDFCard, comment: str=''):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        ncards = len(cards)

        self.n = ncards
        if ncards:
            float_fmt = self.model.float_fmt
            self.node_id = zeros(ncards, 'int32')
            self.phi = zeros(ncards, float_fmt)
            self.ring_id = zeros(ncards, 'int32')

            for i, card in enumerate(cards):
                #: Node ID
                self.node_id[i] = integer(card, 1, 'nid')
                self.ring_id[i] = integer(card, 2, 'rid')
                self.phi[i] = double_or_blank(card, 3, 'phi', 0.)

    def positions(self, node_ids=None):
        raise NotImplementedError('POINTAX.positions')
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

    def write_card(self, bdf_file, size=8, is_double=False):
        if self.n:
            bdf_file.write('$POINTAX\n')
            for (nid, rid, phi) in zip(self.point_id, self.ring_id, self.phi):
                card = ['POINTAX', rid, phi]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def __repr__(self):
        msg = "<POINTAX>\n"
        msg += '  nPOINTAX = %i' % self.n
