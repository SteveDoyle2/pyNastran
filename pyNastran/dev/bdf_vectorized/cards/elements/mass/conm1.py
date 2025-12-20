from numpy import zeros, array, arange, searchsorted

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank, double_or_blank)
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard

class CONM1(VectorizedCard):
    type = 'CONM1'
    def __init__(self, model):
        VectorizedCard.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.element_id = zeros(ncards, 'int32')
            self.node_id = zeros(ncards, 'int32')
            self.coord_id = zeros(ncards, 'int32')

            float_fmt = self.model.float_fmt
            m = zeros((ncards, 6, 6), float_fmt)
            self.mass_matrix = m

    def add_card(self, card: BDFCard, comment: str=''):
        eid = integer(card, 1, 'eid')
        if comment:
            self.set_comment(eid, comment)
        i = self.i
        self.element_id[i] = eid
        self.node_id[i] = integer(card, 2, 'nid')
        self.coord_id[i] = integer_or_blank(card, 3, 'cid', 0)

        m = self.mass_matrix
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
        assert len(card) <= 25, 'len(CONM1 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        if self.n:
            self._cards = []
            self._comments = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def get_mass_matrix(self, i):
        return self.mass_matrix[i, :, :]

    def get_mass_by_element_id(self, element_id=None, total=False):
        m = self.mass_matrix
        if element_id is None:
            mm = m[:, 0, 0] + m[:, 1, 1] + m[:, 2, 2]
        else:
            i = searchsorted(self.element_id, element_id)
            mm = m[i, 0, 0] + m[i, 1, 1] + m[i, 2, 2]
        if total:
            return mm.sum()
        else:
            return mm

    def write_card(self, bdf_file, size=8, element_id=None):
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_id)

            Cid = [cid if cid != 0 else '' for cid in self.coord_id[i]]
            for (eid, nid, cid, m) in zip(self.element_id[i], self.node_id[i], Cid, self.mass_matrix):
                if eid in self._comments:
                    bdf_file.write(self._comments[eid])
                card = [
                    'CONM1', eid, nid, cid, m[0, 0], m[1, 0], m[1, 1],
                    m[2, 0], m[2, 1], m[2, 2], m[3, 0], m[3, 1], m[3, 2],
                    m[3, 3], m[4, 0], m[4, 1], m[4, 2], m[4, 3], m[4, 4],
                    m[5, 0], m[5, 1], m[5, 2], m[5, 3], m[5, 4], m[5, 5]
                ]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))
