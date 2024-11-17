from numpy import zeros

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class EPOINT:
    type = 'EPOINT'
    def __init__(self, model):
        self.model = model
        self._cards = []
        self._comments = []
        self.n = 0

    def add_card(self, card: BDFCard, comment: str=''):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        self.n = len(cards)
        if self.n:
            self.epoint = zeros(self.n, 'int32')
            for i, card in enumerate(cards):
                self.epoint[i] = i
            self.epoint.sort()
            self._cards = []
            self._comments = []

    def write_card(self, bdf_file, size=8, is_double=False):
        #.. todo:: collapse the IDs
        if self.n:
            card = ['EPOINT'] + list(self.epoint)
            if size == 8:
                bdf_file.write(print_card_8(card))
            else:
                bdf_file.write(print_card_16(card))

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('EPOINT', self.n))
        return msg
