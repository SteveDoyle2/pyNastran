from numpy import zeros, arange, searchsorted

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank, string_or_blank)
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.utils import wipe_empty_fields


class SPLINE1(VectorizedCard):
    """
    Surface Spline Methods
    Defines a surface spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by regular arrays of
    aerodynamic points.::

      +---------+-------+-------+------+------+------+----+------+-------+
      | SPLINE1 | EID   | CAERO | BOX1 | BOX2 | SETG | DZ | METH | USAGE |
      +---------+-------+-------+------+------+------+----+------+-------+
      | NELEM   | MELEM |
      +---------+-------+-------+------+------+------+----+
      | SPLINE1 |   3   |  111  | 115  | 122  |  14  | 0. |
      +---------+-------+-------+------+------+------+----+
    """
    type = 'SPLINE1'
    _field_map = {
        1: 'eid', 2:'caero', 3:'box1', 4:'box2', 5:'setg', 6:'dz',
        7: 'method', 8:'usage', 9:'nelements', 10:'melements',
    }
    def __init__(self, model):
        """
        ::

          1
          | \
          |   \
          |     \
          |      4
          |      |
          |      |
          2------3
        """
        VectorizedCard.__init__(self, model)

    def allocate(self, ncards):
        self.n = ncards
        if self.n:
            float_fmt = self.model.float_fmt
            self.element_id = zeros(ncards, dtype='int32')
            self.caero = zeros(ncards, dtype='int32')
            self.box1 = zeros(ncards, dtype='int32')
            self.box2 = zeros(ncards, dtype='int32')
            self.setg = zeros(ncards, dtype='int32')
            self.dz = zeros(ncards, dtype=float_fmt)
            self.method = zeros(ncards, dtype='|U8')
            self.usage = zeros(ncards, dtype='|U8')
            self.nelements = zeros(ncards, dtype='int32')
            self.melements = zeros(ncards, dtype='int32')

    def add_card(self, card, comment=''):
        i = self.i
        self.element_id[i] = integer(card, 1, 'element_id')
        self.caero[i] = integer(card, 2, 'caero')
        self.box1[i] = integer(card, 3, 'box1')
        self.box2[i] = integer(card, 4, 'box2')
        self.setg[i] = integer(card, 5, 'setg')
        self.dz[i] = double_or_blank(card, 6, 'dz', 0.0)
        self.method[i] = string_or_blank(card, 7, 'method', 'IPS')
        self.usage[i] = string_or_blank(card, 8, 'usage', 'BOTH')
        self.nelements[i] = integer_or_blank(card, 9, 'nelements', 10)
        self.melements[i] = integer_or_blank(card, 10, 'melements', 10)
        assert self.nelements[i] > 0, 'nelements = %s' % self.nelements[i]
        assert self.melements[i] > 0, 'melements = %s' % self.melements[i]
        assert len(card) <= 11, 'len(SPLINE1 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def add_op2(self, data, comment=''):
        i = self.i
        self.element_id[i] = data[0]
        self.caero[i] = data[1]
        self.box1[i] = data[2]
        self.box2[i] = data[3]
        self.setg[i] = data[4]
        self.dz[i] = data[5]
        self.method[i] = data[6]
        self.usage[i] = data[7]
        self.nelements[i] = data[8]
        self.melements[i] = data[9]
        assert len(data) == 10, 'data = %s' % data
        self.i += 1

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.caero = self.caero[i]
            self.box1 = self.box1[i]
            self.box2 = self.box2[i]
            self.setg = self.setg[i]
            self.dz = self.dz[i]
            self.method = self.method[i]
            self.usage = self.usage[i]
            self.nelements = self.nelements[i]
            self.melements = self.melements[i]

    def write_card(self, bdf_file, size=8, is_double=True, element_id=None):
        assert size in [8, 16], size
        assert is_double in [True, False], is_double
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                #assert len(unique(element_id))==len(element_id), unique(element_id)
                i = searchsorted(self.element_id, element_id)

            for (eid, caero, box1, box2, setg, dz, method, usage, nelements, melements,
                    ) in zip(self.element_id[i], self.caero[i], self.box1[i],
                    self.box2[i], self.setg[i], self.dz[i],
                    self.method[i], self.usage[i], self.nelements[i], self.melements[i]):

                sdz = set_blank_if_default(dz, 0.)
                smethod = set_blank_if_default(method, 'IPS')
                susage = set_blank_if_default(usage, 'BOTH')
                snelements = set_blank_if_default(nelements, 10)
                smelements = set_blank_if_default(melements, 10)
                list_fields = ['SPLINE1', eid, caero, box1, box2,
                  setg, sdz, smethod, susage, snelements, smelements]
                list_fields = wipe_empty_fields(list_fields)
                bdf_file.write(print_card_8(list_fields))
