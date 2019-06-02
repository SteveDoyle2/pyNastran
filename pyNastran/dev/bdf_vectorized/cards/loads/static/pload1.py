from numpy import arange, array, zeros, searchsorted, unique

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    double, double_or_blank, string)

from pyNastran.dev.bdf_vectorized.cards.loads.vectorized_load import VectorizedLoad

class PLOAD1:
    type = 'PLOAD1'
    valid_types = ['FX', 'FY', 'FZ', 'FXE', 'FYE', 'FZE',
                   'MX', 'MY', 'MZ', 'MXE', 'MYE', 'MZE']
    valid_scales = ['LE', 'FR', 'LEPR', 'FRPR'] # LE: length-based; FR: fractional; PR:projected

    def __init__(self, model):
        """
        Defines the PLOAD object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        #print("force", unique_lid, i)
        if len(i):
            f = PLOAD1(self.model)
            f.load_id = self.load_id[i]
            f.element_id = self.element_id[i]
            f.Type = self.Type[i]
            f.scale = self.scale[i]
            f.x1 = self.x1[i]
            f.x2 = self.x2[i]
            f.p1 = self.p1[i]
            f.p2 = self.p2[i]
            f.n = len(i)
            return f
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        f = PLOAD1(self.model)
        f.load_id = self.load_id
        f.element_id = self.element_id
        f.Type = self.Type
        f.scale = self.scale * value
        f.x1 = self.x1 * value
        f.x2 = self.x2 * value
        f.p1 = self.p1
        f.p2 = self.p2
        f.n = self.n
        return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def add_card(self, card, comment=None):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        #if comment:
            # self.comment = comment
        cards = self._cards
        ncards = len(cards)

        float_fmt = self.model.float_fmt
        self.load_id = zeros(ncards, 'int32')
        self.element_id = zeros(ncards, 'int32')
        self.Type = array([''] * ncards, '|U4')
        self.scale = array([''] * ncards, '|U4')
        self.x1 = zeros(ncards, float_fmt)
        self.x2 = zeros(ncards, float_fmt)
        self.p1 = zeros(ncards, float_fmt)
        self.p2 = zeros(ncards, float_fmt)

        self.n = ncards
        for i, card in enumerate(cards):
            self.load_id[i] = integer(card, 1, 'load_id')
            self.element_id[i] = integer(card, 2, 'eid')
            Type = string(card, 3, 'Type ("%s")' % '",  "'.join(self.valid_types))
            scale = string(card, 4, 'scale ("%s")' % '", "'.join(self.valid_scales))
            self.x1[i] = double(card, 5, 'x1')
            self.p1[i] = double(card, 6, 'p1')
            self.x2[i] = double_or_blank(card, 7, 'x2', self.x1)
            self.p2[i] = double_or_blank(card, 8, 'p2', self.p1)
            assert len(card) <= 9, 'len(PLOAD1 card) = %i\ncard=%s' % (len(card), card)

            if Type not in self.valid_types:
                msg = '%r is an invalid type on the PLOAD1 card' % Type
                raise RuntimeError(msg)
            self.Type[i] = Type

            assert 0.0 <= self.x1[i] <= self.x2[i]
            if scale in ['FR', 'FRPR']:
                assert self.x1[i] <= 1.0, 'x1=%r' % self.x1[i]
                assert self.x2[i] <= 1.0, 'x2=%r' % self.x2[i]
            assert scale in self.valid_scales, '%s is an invalid scale on the PLOAD1 card' % scale
            self.scale[i] = scale
            self._cards = []
            self._comments = []

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('PLOAD1', self.n))
        return msg

    def write_card(self, bdf_file, size=8, is_double=False, load_id=None):
        if self.n:
            if load_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(load_ids, self.load_id)

            for (load_id, element_id, Type, scale, x1, p1, x2, p2) in zip(self.load_id[i], self.element_id[i],
                    self.Type[i], self.scale[i], self.x1[i], self.p1[i], self.x2[i], self.p2[i]):
                card = ['PLOAD1', load_id, element_id, Type, scale, x1, p1, x2, p2]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))
