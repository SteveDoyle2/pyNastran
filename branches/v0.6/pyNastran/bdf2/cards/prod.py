class PROD(object):
    type = 'PROD'
    def __init__(self, model):
        """
        Defines the PROD object.

        :param self: the PROD object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._prod = []
        self._prod_comments = []

    def add(self, card, comment):
        print card
        prod
        self._prod.append(card)
        self._prod_comment.append(comment)

    def build(self):
        """
        :param self: the PROD object
        :param cards: the list of PROD cards
        """
        cards = self._prod
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            #: Property ID
            self.pid = zeros(ncards, 'int32')
            self.mid = zeros(ncards, 'int32')
            self.A = zeros(ncards, float_fmt)
            self.J = zeros(ncards, float_fmt)
            self.c = zeros(ncards, float_fmt)
            self.nsm = zeros(ncards, float_fmt)

            for i, card in enumerate(cards):
                self.pid[i] = integer(card, 1, 'pid')
                self.mid[i] = integer(card, 2, 'mid')
                self.A[i] = double(card, 3, 'A')
                self.J[i] = double_or_blank(card, 4, 'J', 0.0)
                self.c[i] = double_or_blank(card, 5, 'c', 0.0)
                self.nsm[i] = double_or_blank(card, 6, 'nsm', 0.0)
                assert len(card) <= 7, 'len(PROD card) = %i' % len(card)

            i = self.pid.argsort()
            self.pid = self.pid[i]
            self.mid = self.mid[i]
            self.A = self.A[i]
            self.J = self.J[i]
            self.c = self.c[i]
            self.nsm = self.nsm[i]

            unique_pids = unique(self.pid)
            if len(unique_pids) != len(self.pid):
                raise RuntimeError('There are duplicate PROD IDs...')
            self._prods = []
            self._prod_comments = []
        
    def write_bdf(self, f, size=8, pids=None):
        for (pid, mid, A, J, c, nsm) in zip(
             self.pid, self.A, self.J, self.c, self.nsm):

            #self.mid = integer(card, 4, 'mid')
            #self.A = double(card, 5, 'A')
            #self.j = double_or_blank(card, 6, 'j', 0.0)
            #self.c = double_or_blank(card, 7, 'c', 0.0)
            #self.nsm = double_or_blank(card, 8, 'nsm', 0.0)

            card = ['PROD', pid, mid, A, J, c, nsm]
            f.write(print_card(card))