from numpy import searchsorted, arange, zeros
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard
from pyNastran.bdf.field_writer_8 import print_card_8
#from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import integer, integer_or_blank


class PAERO1(VectorizedCard):
    """
    Defines associated bodies for the panels in the Doublet-Lattice method.::

      PAERO1 PID B1 B2 B3 B4 B5 B6
    """
    type = 'PAERO1'

    def __init__(self, model):
        VectorizedCard.__init__(self, model)

    def allocate(self, ncards):
        self.n = ncards
        if self.n:
            #float_fmt = self.model.float_fmt
            self.property_id = zeros(ncards, dtype='int32')
            self.b = zeros((ncards, 6), dtype='int32')

    def add_card(self, card, comment=''):
        i = self.i
        if comment:
            self.comment = comment
        self.property_id[i] = integer(card, 1, 'pid')
        for j in range(2, len(card)):
            self.b[i, j-2] = integer_or_blank(card, j, 'b%i' % j-2, 0)
        self.i += 1

    def build(self):
        if self.n:
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.b = self.b[i, :]

    #def get_bodies(self):
        #return self.Bi

    def write_card(self, bdf_file, size, is_double, property_id=None):
        if self.n:
            assert size in [8, 16], size
            assert is_double in [True, False], is_double

            if property_id is None:
                i = arange(self.n)
            else:
                #assert len(unique(property_id))==len(property_id), unique(property_id)
                i = searchsorted(self.property_id, property_id)
            #print(self.property_id)
            for j, pid in enumerate(self.property_id[i]):
                b = self.b[j, :]
                b2 = [bi for bi in b if bi > 0]
                list_fields = ['PAERO1', pid] + b2
                bdf_file.write(print_card_8(list_fields))
        #return self.comment + print_card_8(card)
