from numpy import zeros, unique

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, integer_double_or_blank, integer_string_or_blank,
    string_or_blank, blank)


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
        self._cards = []
        self._comments = []

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        """
        :param self: the PROD object
        :param cards: the list of PROD cards
        """
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            self.material_id = zeros(ncards, 'int32')
            self.A = zeros(ncards, float_fmt)
            self.J = zeros(ncards, float_fmt)
            self.c = zeros(ncards, float_fmt)
            self.nsm = zeros(ncards, float_fmt)

            for i, card in enumerate(cards):
                self.property_id[i] = integer(card, 1, 'pid')
                self.material_id[i] = integer(card, 2, 'mid')
                self.A[i] = double(card, 3, 'A')
                self.J[i] = double_or_blank(card, 4, 'J', 0.0)
                self.c[i] = double_or_blank(card, 5, 'c', 0.0)
                self.nsm[i] = double_or_blank(card, 6, 'nsm', 0.0)
                assert len(card) <= 7, 'len(PROD card) = %i' % len(card)

            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]
            self.A = self.A[i]
            self.J = self.J[i]
            self.c = self.c[i]
            self.nsm = self.nsm[i]

            unique_pids = unique(self.property_id)
            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PROD IDs...')
            self._cards = []
            self._comments = []
        
    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('PROD', self.n))
        return msg

    def write_bdf(self, f, size=8, pids=None):
        for (pid, mid, A, J, c, nsm) in zip(
             self.property_id, self.material_id, self.A, self.J, self.c, self.nsm):

            #self.mid = integer(card, 4, 'mid')
            #self.A = double(card, 5, 'A')
            #self.j = double_or_blank(card, 6, 'j', 0.0)
            #self.c = double_or_blank(card, 7, 'c', 0.0)
            #self.nsm = double_or_blank(card, 8, 'nsm', 0.0)

            card = ['PROD', pid, mid, A, J, c, nsm]
            f.write(print_card(card))