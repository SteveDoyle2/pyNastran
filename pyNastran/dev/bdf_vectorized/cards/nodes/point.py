import numpy as np
from pyNastran.bdf.bdf_interface.assign_type import double_or_blank, integer_or_blank, integer

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class POINT:
    type = 'POINT'
    def __init__(self, model):
        """
        Defines the POINT object.

        Parameters
        ----------
        model : BDF
           the BDF object

        +-------+-----+----+----+----+----+
        |   1   |  2  | 3  | 4  | 5  | 6  |
        +=======+=====+====+====+====+====+
        | POINT | NID | CP | X1 | X2 | X3 |
        +-------+-----+----+----+----+----+
        """
        self.model = model
        self.n = 0
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
            self.node_id = np.zeros(ncards, 'int32')
            self.xyz = np.zeros((ncards, 3), float_fmt)
            self.coord_id = np.zeros(ncards, 'int32')

            cp0 = self.model.grdset.cp
            for i, card in enumerate(cards):
                #: Node ID
                self.node_id[i] = integer(card, 1, 'nid')

                #: Grid point coordinate system
                self.coord_id[i] = integer_or_blank(card, 2, 'cp', cp0)

                x = double_or_blank(card, 3, 'x1', 0.)
                y = double_or_blank(card, 4, 'x2', 0.)
                z = double_or_blank(card, 5, 'x3', 0.)
                #: node location in local frame
                self.xyz[i] = [x, y, z]

    def positions(self, node_ids=None):
        if node_ids is None:
            node_ids = self.node_id
        xyz = xyz.copy()

        n = np.arange(self.n)
        i = np.where(self.coord_id != 0)[0]
        if i:
            n = n[i]
            cids = set(list(np.unique(self.coord_id)))
            for cid in cids:
                i = np.where(self.coord_id != 0)[0]
                T = self.model.coord.transform(cid)
                xyzi = xyz[n[i], :]
                xyzi = np.dot(np.transpose(T), np.dot(xyzi, T))
        return xyz

    def positions_wrt(self, node_ids=None, cids=None):
        raise NotImplementedError()

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('POINT', self.n))
        return msg

    def write_card(self, bdf_file, size=8, is_double=False):
        if self.n:
            bdf_file.write('$POINT\n')
            cp0 = self.model.grdset.cp
            Cp = [cpi if cpi != cp0 else '' for cpi in self.cp]
            for (nid, cp, xyz) in zip(self.nid, Cp, self.xyz):
                card = ['POINT', nid, cp, xyz[0], xyz[1], xyz[2]]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def __repr__(self):
        msg = "<POINT>\n"
        msg += '  nPOINT = %i' % self.n
