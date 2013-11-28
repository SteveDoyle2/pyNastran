from numpy import arange, array, zeros, searchsorted, unique

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string, string_or_blank)


class PLOAD2(object):
    type = 'PLOAD2'

    def __init__(self, model):
        """
        Defines the PLOAD2 object.

        :param self: the PLOAD2 object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        #self._cards = []
        self._comments = []

        self._load_id = []
        self._element_ids = []
        self._p = []
        
    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        #print "force", unique_lid, i
        if len(i):
            f = PLOAD2(self.model)
            f.load_id = self.load_id[i]
            f.element_id = self.element_id[i]
            f.p = self.p[i]
            f.n = len(i)
            return f
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        f = PLOAD2(self.model)
        f.load_id = self.load_id
        f.element_id = self.element_id
        f.p = self.p[i]
        f.n = self.n
        return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def add(self, card, comment=None):
        #self._comments.append(comment)
        self._load_id.append(integer(card, 1, 'load_id'))
        self._p.append(double(card, 2, 'p'))

        if integer_string_or_blank(card, 4, 'THRU') == 'THRU':
            e1 = integer(card, 3, 'Element1')
            e2 = integer(card, 5, 'Element1')
            eids = [i for i in xrange(e1, e2 + 1)]
        else:
            eids = fields(integer, card, 'eid', i=3, j=len(card))
        assert len(card) == 6, 'len(PLOAD2 card) = %i' % len(card)
        self.n += len(eids)
        self._eids.append(eids)

    def build(self):
        #if comment:
            #self._comment = comment
        #cards = self._cards
        ncards = self.n
        
        float_fmt = self.model.float
        self.load_id = zeros(ncards, 'int32')
        self.element_id = zeros(ncards, 'int32')
        self.p = zeros(ncards, float_fmt)

        self.n = ncards
        istart = 0
        iend = None
        for i in xrange(self.n):
            element_ids = self._element_ids[i]
            n = len(eids)
            iend = istart + n

            self.load_id[istart:iend] = self._load_id[i]
            self.element_id[istart:iend] = element_ids
            self.p[istart : iend] = self._p[i]
            assert len(card) <= 9, 'len(PLOAD2 card) = %i' % len(card)

            self._load_id = []
            self._element_ids = []
            self._p = []
            self._comments = []

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('PLOAD2', self.n))
        return msg

    def write_bdf(self, f, size=8, load_ids=None):
        aaa
        if self.n:
            if load_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(load_ids, self.load_id)

            for (load_id, element_id, p) in zip(self.load_id[i], self.element_id[i], self.p[i]):
                card = ['PLOAD2', load_id, element_id, p]
                f.write(print_card(card))

