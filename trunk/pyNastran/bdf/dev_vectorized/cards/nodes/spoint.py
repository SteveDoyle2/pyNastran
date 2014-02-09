from numpy import zeros

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.cards.baseCard import expand_thru, collapse_thru
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank, integer_or_string)

class SPOINT(object):
    type = 'SPOINT'
    def __init__(self, model):
        self.model = model
        self._comments = []
        self.spoint = set([])
        self.n = 0

    def add(self, card, comment):
        fields = []
        for i in range(1, len(card)):
            field = integer_or_string(card, i, 'ID%i' % i)
            fields.append(field)
        self.spoint.union( set(expand_thru(fields)) )
        self._comments.append(comment)

    def build(self):
        self.n = len(self.spoint)
        if self.n:
            self.spoint = array(self.spoint, 'int32')
            self.spoint.sort()
            self._comments = []

    def write_bdf(self, f, size=8):
        #..todo:: collapse the IDs
        if self.n:
            card = ['SPOINT'] + collapse_thru(self.spoint)
            f.write(print_card(card))

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('SPOINT', self.n))
        return msg