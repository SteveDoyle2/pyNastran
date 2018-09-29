from numpy import zeros
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard
from pyNastran.dev.bdf_vectorized.bdf_interface.assign_type import integer_or_blank, double_or_blank, double
from pyNastran.bdf.field_writer_8 import print_card_8


class AERO(VectorizedCard):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.

    +------+-------+----------+------+--------+-------+-------+
    | 1    | 2     | 3        | 4    | 5      | 6     | 7     |
    +------+-------+----------+------+--------+-------+-------+
    | AERO | ACSID | VELOCITY | REFC | RHOREF | SYMXZ | SYMXY |
    +------+-------+----------+------+--------+-------+-------+
    | AERO | 3     | 1.3+4    | 100. |  1.-5  | 1     | -1    |
    +------+-------+----------+------+--------+-------+-------+
    """
    type = 'AERO'
    _field_map = {
        1: 'acsid', 2:'velocity', 3:'cRef', 4:'rhoRef', 5:'symXZ',
        6:'symXY',
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
        self.model.log.debug('AERO.n=%s' % self.n)
        if self.n:
            float_fmt = self.model.float_fmt
            self.acsid = zeros(ncards, dtype='int32')
            self.velocity = zeros(ncards, dtype=float_fmt)
            self.cRef = zeros(ncards, dtype=float_fmt)
            self.rhoRef = zeros(ncards, dtype=float_fmt)
            self.symXZ = zeros(ncards, dtype='int32')
            self.symXY = zeros(ncards, dtype='int32')

    def add_card(self, card, comment=''):
        i = self.i
        if comment:
            self.comment = comment
        self.acsid[i] = integer_or_blank(card, 1, 'acsid', 0)
        self.velocity[i] = double_or_blank(card, 2, 'velocity')
        self.cRef[i] = double(card, 3, 'cRef')
        self.rhoRef[i] = double(card, 4, 'rhoRef')
        self.symXZ[i] = integer_or_blank(card, 5, 'symXZ', 0)
        self.symXY[i] = integer_or_blank(card, 6, 'symXY', 0)
        assert len(card) <= 7, 'len(AERO card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def add_op2(self, data):
        self.acsid = data[0]
        self.velocity = data[1]
        self.cRef = data[2]
        self.rhoRef = data[3]
        self.symXZ = data[4]
        self.symXY = data[5]
        assert len(data) == 6, 'data = %s' % data

        # T is the tabular function
        #angle = self.wg*self.t*(t-(x-self.x0)/self.V)

    def build(self):
        if self.n:
            i = self.acsid.argsort()
            self.acsid = self.acsid[i]
            self.velocity = self.velocity[i]
            self.cRef = self.cRef[i]
            self.rhoRef = self.rhoRef[i]
            self.symXZ = self.symXZ[i]
            self.symXY = self.symXY[i]

    def write_card(self, bdf_file, size, is_double):
        #card = self.repr_fields()
        #symXZ = set_blank_if_default(self.symXZ, 0)
        #symXY = set_blank_if_default(self.symXY, 0)
        if self.n:
            for acsid, V, c, rho, xz, xy in zip(self.acsid, self.velocity, self.cRef,
                                                self.rhoRef, self.symXZ, self.symXY):

                list_fields = ['AERO', acsid, V, c, rho, xz, xy]
                bdf_file.write(print_card_8(list_fields))
