from six.moves import zip, range
from numpy import arange, array, zeros, searchsorted, unique, full, nan, where

from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_string_or_blank, string_or_blank)


class PLOAD4(object):
    type = 'PLOAD4'

    def __init__(self, model):
        """
        Defines the PLOAD4 object.

        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

        #self._load_id = []
        #self._element_ids = []
        #self._p = []

    #def __getitem__(self, i):
        #unique_lid = unique(self.load_id)
        ##print "force", unique_lid, i
        ##if len(i):
        #f = PLOAD4(self.model)
        #f.load_id = self.load_id[i]
        #f.element_id = self.element_id[i]
        #f.p = self.p[i]
        #f.n = len(i)
        #return f
        ##raise RuntimeError('len(i) = 0')

    def __getitem__(self, load_id):
        i = where(load_id == self.load_id)[0]
        return self.slice_by_index(i)

    def __mul__(self, value):
        raise NotImplementedError()
        #f = PLOAD4(self.model)
        #f.load_id = self.load_id
        #f.element_id = self.element_id
        #f.p = self.p[i]
        #f.n = self.n
        #return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def add(self, card, comment=None):
        self._cards.append(card)
        self._comments.append(comment)

    def allocate(self, ncards):
        float_fmt = self.model.float
        self.load_id = zeros(ncards, 'int32')
        #self.element_id = zeros(ncards, 'int32')
        self.pressures = zeros((ncards, 4), 'int32')

        self.element_ids = {}
        for i in range(ncards):
            self.element_ids[i] = []

        self.g1 = full(ncards, nan, 'int32')
        self.g34 = full(ncards, nan, 'int32')
        self.ldir = full(ncards, nan, '|S4')
        self.sorl = full(ncards, nan, '|S4')
        self.cid = zeros(ncards, dtype='int32')
        self.NVector = zeros((ncards, 3), dtype=float_fmt)

    def build(self):
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float

            self.load_id = zeros(ncards, 'int32')
            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.pressures = zeros((ncards, 4), 'int32')

            element_ids = {}
            for i in range(ncards):
                element_ids[i] = []

            self.g1 = full(ncards, nan, 'int32')
            self.g34 = full(ncards, nan, 'int32')

            self.ldir = full(ncards, nan, '|S4')
            self.sorl = full(ncards, nan, '|S4')

            self.cid = zeros(ncards, dtype='int32')
            self.NVector = zeros((ncards, 3), dtype=float_fmt)

            for i, card in enumerate(cards):
                self.load_id[i] = integer(card, 1, 'load_id')
                eid = integer(card, 2, 'element_id')
                self.element_id[i] = eid
                p1 = double_or_blank(card, 3, 'p1', 0.0)
                p = [p1,
                     double_or_blank(card, 4, 'p2', p1),
                     double_or_blank(card, 5, 'p3', p1),
                     double_or_blank(card, 6, 'p4', p1)]
                self.pressures[i, :] = p

                self.element_ids[i] = [eid]
                if(integer_string_or_blank(card, 7, 'g1/THRU') == 'THRU' and
                   integer_or_blank(card, 8, 'eid2')):  # plates
                    eid2 = integer(card, 8, 'eid2')
                    if eid2:
                        self.element_ids[i] = list(unique(expand_thru([self.eid, 'THRU', eid2],
                                                   set_fields=False, sort_fields=False)))
                    #self.g1 = None
                    #self.g34 = None
                else:
                    #: used for CPENTA, CHEXA
                    self.element_ids[i] = [self.eid]
                    #: used for solid element only
                    self.g1[i] = integer_or_blank(card, 7, 'g1')
                    #: g3/g4 - different depending on CHEXA/CPENTA or CTETRA
                    self.g34[i] = integer_or_blank(card, 8, 'g34')

                #: Coordinate system identification number. See Remark 2.
                #: (Integer >= 0;Default=0)
                self.cid[i] = integer_or_blank(card, 9, 'cid', 0)
                self.NVector[i, :] = [
                    double_or_blank(card, 10, 'N1', 0.0),
                    double_or_blank(card, 11, 'N2', 0.0),
                    double_or_blank(card, 12, 'N3', 0.0), ]
                self.sorl[i] = string_or_blank(card, 13, 'sorl', 'SURF')
                self.ldir[i] = string_or_blank(card, 14, 'ldir', 'NORM')
                assert len(card) <= 15, 'len(PLOAD4 card) = %i\ncard=%s' % (len(card), card)

            i = self.load_id.argsort()
            #self.element_id = self.element_id[i]
            self.pressures = self.pressures[i, :]
            #self.node_ids = self.node_ids[i, :]

            #element_ids = {}
            #for j in range(ncards):
                #element_ids[j] = element_ids[i[j]]

            self.g1 = self.g1[i]
            self.g34 = self.g34[i]
            self.ldir = self.ldir[i]
            self.sorl = self.sorl[i]
            self.cid = self.cid[i]
            self.NVector = self.NVector[i, :]
            self.n += len(eids)

            self.load_cases = {}
            self.load_ids = unique(self.load_id)
        else:
            self.load_id = array([], dtype='int32')

    #def build(self):
        # PLOAD2
        #ncards = self.n

        #float_fmt = self.model.float
        #self.load_id = zeros(ncards, 'int32')
        #self.element_id = zeros(ncards, 'int32')
        #self.p = zeros(ncards, float_fmt)

        #self.n = ncards
        #istart = 0
        #iend = None
        #for i in range(self.n):
            #element_ids = self._element_ids[i]
            #n = len(element_ids)
            #iend = istart + n

            #self.load_id[istart:iend] = self._load_id[i]
            #self.element_id[istart:iend] = element_ids
            #self.p[istart : iend] = self._p[i]

            #self._load_id = []
            #self._element_ids = []
            #self._p = []
            #self._comments = []

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('PLOAD4', self.n))
        return msg

    def write_card(self, f, size=8, load_ids=None):
        if self.n:
            if load_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(load_ids, self.load_id)

            for (load_id, element_id, p) in zip(self.load_id[i], self.element_id[i], self.pressures[i]):
                card = ['PLOAD4', load_id, element_id, p]
                f.write(print_card_8(card))

