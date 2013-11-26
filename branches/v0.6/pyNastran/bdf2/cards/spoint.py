from numpy import zeros

class SPOINT(object):
    type = 'SPOINT'
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
        self.spoint = zeros(ncards, 'int32')
        self.n = ncards
        if ncards:
            for i, card in enumerate(cards):
                self.spoint[i] = i
            self.spoint.sort()
            self._cards = []
            self._comments = []

    def write_bdf(self, f, size=8):
        #..todo:: collapse the IDs
        if self.n:
            card = ['SPOINT'] + list(self.spoint)
        f.write(print_card(card))

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('SPOINT', self.n))
        return msg