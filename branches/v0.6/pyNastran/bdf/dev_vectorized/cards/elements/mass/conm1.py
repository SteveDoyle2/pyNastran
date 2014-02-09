import StringIO
from itertools import izip

from numpy import zeros, array, arange, unique, searchsorted, where

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank, double_or_blank)


class CONM1(object):
    type = 'CONM1'
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
        ncards = len(cards)

        self.n = ncards
        if ncards:
            self.element_id = zeros(ncards, 'int32')
            self.node_id = zeros(ncards, 'int32')
            self.coord_id = zeros(ncards, 'int32')
            
            float_fmt = self.model.float
            m = zeros((ncards, 6, 6), float_fmt)
            self.mass_matrix = mass

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'eid')
                self.node_id[i] = integer(card, 2, 'nid')
                self.coord_id[i] = integer_or_blank(card, 3, 'cid', 0)

                m[i, 0, 0] = double_or_blank(card, 4, 'M11', 0.)
                m[i, 1, 0] = double_or_blank(card, 5, 'M21', 0.)
                m[i, 1, 1] = double_or_blank(card, 6, 'M22', 0.)
                m[i, 2, 0] = double_or_blank(card, 7, 'M31', 0.)
                m[i, 2, 1] = double_or_blank(card, 8, 'M32', 0.)
                m[i, 2, 2] = double_or_blank(card, 9, 'M33', 0.)
                m[i, 3, 0] = double_or_blank(card, 10, 'M41', 0.)
                m[i, 3, 1] = double_or_blank(card, 11, 'M42', 0.)
                m[i, 3, 2] = double_or_blank(card, 12, 'M43', 0.)
                m[i, 3, 3] = double_or_blank(card, 13, 'M44', 0.)
                m[i, 4, 0] = double_or_blank(card, 14, 'M51', 0.)
                m[i, 4, 1] = double_or_blank(card, 15, 'M52', 0.)
                m[i, 4, 2] = double_or_blank(card, 16, 'M53', 0.)
                m[i, 4, 3] = double_or_blank(card, 17, 'M54', 0.)
                m[i, 4, 4] = double_or_blank(card, 18, 'M55', 0.)
                m[i, 5, 0] = double_or_blank(card, 19, 'M61', 0.)
                m[i, 5, 1] = double_or_blank(card, 20, 'M62', 0.)
                m[i, 5, 2] = double_or_blank(card, 21, 'M63', 0.)
                m[i, 5, 3] = double_or_blank(card, 22, 'M64', 0.)
                m[i, 5, 4] = double_or_blank(card, 23, 'M65', 0.)
                m[i, 5, 5] = double_or_blank(card, 24, 'M66', 0.)
                assert len(card) <= 25, 'len(CONM1 card) = %i' % len(card)
            self._cards = []
            self._comments = []

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_ids)

            Cid = [cid if cid != 0 else '' for cid in self.coord_id[i]]
            for (eid, nid, cid, m) in izip(self.element_id[i], self.node_id[i], Cid, self.mass_matrix):
                    card = ['CONM1', eid, nid, cid, m[0, 0], m[1, 0], m[1, 1],
                              m[2, 0], m[2, 1], m[2, 2], m[3, 0], m[3, 1], m[3, 2],
                              m[3, 3], m[4, 0], m[4, 1], m[4, 2], m[4, 3], m[4, 4],
                              m[5, 0], m[5, 1], m[5, 2], m[5, 3], m[5, 4], m[5, 5]]
                    f.write(print_card(card, size=size))

    def __repr__(self):
        f = StringIO.StringIO()
        self.write_bdf(f)
        return f.getvalue()