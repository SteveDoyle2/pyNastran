class SPOINT(object):
    type = 'SPOINT'
    def __init__(self, model):
        self.model = model
        self._spoint = []
        self._spoint_comments = []

    def add_spoint(self, card, comment):
        self._spoint.append(card)
        self._spoint_comment.append(comment)

    def build(self):
        cards = self._spoint
        self.spoint = zeros(ncards, 'int32')
        self.n = ncards
        for i, card in enumerate(cards):
            pass
        self.spoint.sort()
        self._spoint = []
        self._spoint_comment = []

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