from numpy import array, zeros, arange, concatenate, searchsorted, where, unique

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)


class PCOMP(object):
    type = 'PCOMP'
    def __init__(self, model):
        """
        Defines the PCOMP object.

        :param self: the PCOMP object
        :param model: the BDF object
        :param cards: the list of PCOMP cards
        """
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
        #return
        if ncards:
            #: Property ID
            self.pid = zeros(ncards, 'int32')

            ncards = len(cards)
            for i, card in enumerate(cards):
                self.pid[i] = integer(card, 1, 'pid')
            i = self.pid.argsort()
            self.pid = self.pid[i]
            unique_pids = unique(self.pid)

            if len(unique_pids) != len(self.pid):
                raise RuntimeError('There are duplicate PCOMP IDs...')
            self._cards = []
            self._comments = []
        
    def write_bdf(self, f, size=8, property_ids=None):
        pass