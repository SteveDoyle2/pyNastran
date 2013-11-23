class PSHEAR(object):
    type = 'PSHEAR'
    def __init__(self, model):
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
            self.pid = zeros(ncards, 'int32')
            #: Material ID
            self.mid = zeros(ncards, 'int32')
            self.thickness = zeros(ncards, float_fmt)
            self.nsm = zeros(ncards, float_fmt)
            self.f1 = zeros(ncards, float_fmt)
            self.f2 = zeros(ncards, float_fmt)

            for i, card in enumerate(cards):
                self.pid[i] = integer(card, 1, 'pid')
                self.mid[i] = integer(card, 2, 'mid')
                self.thickness[i] = double(card, 3, 't')
                self.nsm[i] = double_or_blank(card, 4, 'nsm', 0.0)
                self.f1[i] = double_or_blank(card, 5, 'f1', 0.0)
                self.f2[i] = double_or_blank(card, 6, 'f2', 0.0)
                assert self.thickness[i] > 0.0
                #assert self.f1 >= 0.0
                #assert self.f2 >= 0.0
                assert len(card) <= 7, 'len(PSHEAR card) = %i' % len(card)

            i = self.pid.argsort()
            self.pid = self.pid[i]
            self.mid = self.mid[i]
            self.thickness = self.thickness[i]
            self.nsm = self.nsm[i]
            self.f1 = self.f1[i]
            self.f2 = self.f2[i]

            unique_pids = unique(self.pid)
            if len(unique_pids) != len(self.pid):
                raise RuntimeError('There are duplicate PSHEAR IDs...')
            self._cards = []
            self._comments = []

    def write_bdf(self, f, size=8, pids=None):
        pass

