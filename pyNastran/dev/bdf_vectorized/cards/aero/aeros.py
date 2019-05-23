from numpy import zeros
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard

from pyNastran.dev.bdf_vectorized.bdf_interface.assign_type import integer_or_blank, double
#from pyNastran.bdf.field_writer_8 import set_blank_if_default
#from pyNastran.bdf.cards.utils import wipe_empty_fields
from pyNastran.bdf.field_writer_8 import print_card_8

class AEROS(VectorizedCard):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.

    +-------+-------+-------+------+------+-------+------+-------+
    | 1     | 2     | 3     | 4    | 5    | 6     | 7    |   8   |
    +-------+-------+-------+------+------+-------+------+-------+
    | AEROS | ACSID | RCSID | REFC | REFB | REFS  |SYMXZ | SYMXY |
    +-------+-------+-------+------+------+-------+------+-------+

    +-------+-------+-------+------+------+-------+------+-------+
    | AEROS | 10    | 20    | 10.  | 100. | 1000. | 1    |       |
    +-------+-------+-------+------+------+-------+------+-------+

    """
    type = 'AEROS'
    _field_map = {
        1: 'acsid', 2:'rcsid', 3:'cRef', 4:'bRef', 5:'Sref',
        6:'symXZ', 7:'symXY',
    }

    def __init__(self, model):
        VectorizedCard.__init__(self, model)

    def __iter__(self):
        if len(self) == 0:
            raise StopIteration
        for acsid in self.acsid:
            yield acsid

    def allocate(self, ncards):
        self.n = ncards
        self.model.log.debug('AEROS.n=%s' % self.n)
        if self.n:
            float_fmt = self.model.float_fmt
            self.acsid = zeros(ncards, dtype='int32')
            self.rcsid = zeros(ncards, dtype='int32')
            self.cRef = zeros(ncards, dtype=float_fmt)
            self.bRef = zeros(ncards, dtype=float_fmt)
            self.SRef = zeros(ncards, dtype=float_fmt)
            self.symXZ = zeros(ncards, dtype='int32')
            self.symXY = zeros(ncards, dtype='int32')

    def add_card(self, card, comment=''):
        if comment:
            self.comment = comment
        i = self.i
        self.acsid[i] = integer_or_blank(card, 1, 'acsid', 0)
        self.rcsid[i] = integer_or_blank(card, 2, 'rcsid', 0)
        self.cRef[i] = double(card, 3, 'cRef')
        self.bRef[i] = double(card, 4, 'bRef')
        self.Sref[i] = double(card, 5, 'Sref')
        self.symXZ[i] = integer_or_blank(card, 6, 'symXZ', 0)
        self.symXY[i] = integer_or_blank(card, 7, 'symXY', 0)
        assert len(card) <= 8, 'len(AEROS card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def add_op2(self, data):
        self.acsid = data[0]
        self.rcsid = data[1]
        self.cRef = data[2]
        self.bRef = data[3]
        self.SRef = data[4]
        self.symXZ = data[5]
        self.symXY = data[6]
        assert len(data) == 7, 'data = %s' % data

    def build(self):
        if self.n:
            i = self.acsid.argsort()
            self.acsid = self.acsid[i]
            self.rcsid = self.rcsid[i]
            self.cRef = self.cRef[i]
            self.bRef = self.bRef[i]
            self.SRef = self.SRef[i]
            self.symXZ = self.symXZ[i]
            self.symXY = self.symXY[i]


    def write_card(self, bdf_file, size, is_double):
        if self.n == 0:
            return
        for acsid, rcsid, c, b, S, symXZ, symXY in zip(
                self.acsid, self.rcsid, self.cRef, self.bRef, self.SRef, self.symXZ, self.symXY):

            #symXZ = set_blank_if_default(self.symXZ, 0)
            #symXY = set_blank_if_default(self.symXY, 0)
            list_fields = ['AEROS', acsid, rcsid, c,
                           b, S, symXZ, symXY]
            bdf_file.write(print_card_8(list_fields))
