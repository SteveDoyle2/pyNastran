from numpy import zeros, array

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.cards.baseCard import expand_thru, collapse_thru
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank, integer_or_string)

class SPOINT(object):
    type = 'SPOINT'
    def __init__(self, model):
        self.model = model
        self._comments = []
        self.spoint = set([])
        self.n = 0

    def __len__(self):
        return self.n

    def add(self, card, comment=''):
        fields = []
        for i in range(1, len(card)):
            field = integer_or_string(card, i, 'ID%i' % i)
            fields.append(field)

        ex = expand_thru(fields)
        ex = set(ex)

        self.spoint.update(ex)
        self._comments.append(comment)
        self.n = len(self.spoint)

    def remove(self, values):
        s = set(values)
        for si in s:
            self.spoint.remove(si)
            self.n -= 1

    def discard(self, values):
        s = set(values)
        for si in s:
            self.spoint.discard(si)
        self.n = len(self.spoint)

    def max(self):
        return max(self.spoint)

    def min(self):
        return min(self.spoint)

    def build(self):
        if self.n:
            self._comments = []

    def write_bdf(self, f, size=8, is_double=False):
        #..todo:: collapse the IDs
        if self.n:
            spoint = list(self.spoint)
            spoint.sort()
            card = ['SPOINT'] + collapse_thru(spoint)
            f.write(print_card_8(card))

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('SPOINT', self.n))
        return msg

    def __iter__(self):
        for spoint in self.spoint:
            yield spoint

    #def values(self):
        #for spoint in self.spoint:
            #yield spoint

    #def items(self):
        #mids = self.material_id
        ##self.model.log.debug('mids = %s' % mids)
        #for mid in mids:
            #yield mid, self.__getitem__(mid)

    #def iteritems(self):
        #return self.items()

    #def iterkeys(self):
        #return self.keys()

