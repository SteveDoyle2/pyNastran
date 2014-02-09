from numpy import zeros

from pyNastran.bdf.fieldWriter import print_card


class EPOINT(object):
    type = 'EPOINT'
    def __init__(self, model):
        self.model = model
        self._cards = []
        self._comments = []
        self.n = 0

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        self.n = len(cards)
        if self.n:
            self.epoint = zeros(ncards, 'int32')
            for i, card in enumerate(cards):
                self.epoint[i] = i
            self.epoint.sort()
            self._cards = []
            self._comments = []

    def write_bdf(self, f, size=8):
        #..todo:: collapse the IDs
        if self.n:
            card = ['EPOINT'] + list(self.epoint)
            f.write(print_card(card))

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('EPOINT', self.n))
        return msg