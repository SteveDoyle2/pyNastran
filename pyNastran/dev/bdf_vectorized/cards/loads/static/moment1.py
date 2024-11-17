from io import StringIO

from numpy import zeros, unique

#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double, double_or_blank)

from pyNastran.dev.bdf_vectorized.cards.loads.vectorized_load import VectorizedLoad

class MOMENT1(VectorizedLoad):
    type = 'MOMENT1'
    def __init__(self, model):
        """
        Defines the MOMENT1 object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        VectorizedLoad.__init__(self, model)

    def get_load_ids(self):
        return unique(self.load_id)

    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        if len(i):
            f = MOMENT1(self.model)
            f.load_id = self.load_id[i]
            f.node_id = self.node_id[i]
            f.coord_id = self.coord_id[i]
            f.mag = self.mag[i]
            f.xyz = self.xyz[i]
            f.n = len(i)
            return f
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        f = MOMENT1(self.model)
        f.load_id = self.load_id
        f.node_id = self.node_id
        f.coord_id = self.coord_id
        f.mag = self.mag * value
        f.xyz = self.xyz * value
        f.n = self.n
        return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i
        self.load_id[i] = integer(card, 1, 'sid')
        self.node_id[i] = integer(card, 2, 'node')
        self.coord_id[i] = integer_or_blank(card, 3, 'cid', 0)
        self.mag[i] = double(card, 4, 'mag')
        xyz = [double_or_blank(card, 5, 'X1', 0.0),
               double_or_blank(card, 6, 'X2', 0.0),
               double_or_blank(card, 7, 'X3', 0.0)]
        self.xyz[i] = xyz
        assert len(card) <= 8, 'len(MOMENT card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
        float_fmt = self.model.float_fmt
        self.load_id = zeros(ncards, 'int32')
        self.node_id = zeros(ncards, 'int32')
        self.coord_id = zeros(ncards, 'int32')
        self.mag = zeros(ncards, float_fmt)
        self.xyz = zeros((ncards, 3), float_fmt)

    def build(self):
        if self.n:
            i = self.load_id.argsort()
            self.load_id = self.load_id[i]
            self.node_id = self.node_id[i]
            self.coord_id = self.coord_id[i]
            self.mag = self.mag[i]
            self.xyz = self.xyz[i]
            self._comments = []

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('MOMENT', self.n))
        return msg

    def write_card_by_index(self, bdf_file, size=8, is_double=False, i=None):
        for (lid, nid, cid, mag, xyz) in zip(
             self.load_id, self.node_id, self.coord_id, self.mag, self.xyz):

            card = ['MOMENT1', lid, nid, cid, mag, xyz[0], xyz[1], xyz[2]]
            if size == 8:
                bdf_file.write(print_card_8(card))
            else:
                bdf_file.write(print_card_16(card))

    def __repr__(self):
        string_io = StringIO()
        self.write_card(string_io)
        return string_io.getvalue().rstrip()
