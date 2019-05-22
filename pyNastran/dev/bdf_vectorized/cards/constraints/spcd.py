import numpy as np
from numpy import array

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class SPCD:
    """
    Defines an enforced displacement value for static analysis and an enforced
    motion value (displacement, velocity or acceleration) in dynamic analysis.

     +------+-----+-----+-----+------+----+----+----+
     |   1  |  2  |  3  |  4  |   5  |  6 | 7  |  8 |
     +======+=====+=====+=====+======+====+====+====+
     | SPCD | SID |  G1 | C1  |  D1  | G2 | C2 | D2 |
     +------+-----+-----+-----+------+----+----+----+
     | SPCD | 100 | 32  | 436 | -2.6 | 5  | 2  | .9 |
     +------+-----+-----+-----+------+----+----+----+
    """
    type = 'SPCD'

    def __init__(self, model):
        self.model = model
        self.n = 0
        self.i = 0

        self._comments = []
        self.constraint_id = None
        #self.grid_id = []
        #self.components = []

    def allocate(self, card_count):
        ncards = card_count['SPCD']
        if ncards:
            self.n = ncards
            #print('ngrid=%s' % self.n)
            float_fmt = self.model.float_fmt
            self.node_id = np.zeros(ncards, 'int32')
            self.components = np.zeros(ncards, 'int32')
            self.enforced_motion = np.zeros(ncards, float_fmt)

    def add(self, constraint_id, node_id, dofs, enforced_motion, comment):
        i = self.i
        assert enforced_motion != 0.0, enforced_motion

        self._comments.append(comment)
        self.node_id[i] = node_id
        self.components[i] = dofs
        self.enforced_motion[i] = enforced_motion

        if self.constraint_id == constraint_id:
            pass
        elif self.constraint_id is None:
            self.constraint_id = constraint_id
        else:
            msg = ('self.constraint_id == constraint_id; constraint_id=%r '
                   'expected; found=%r' % (self.constraint_id, constraint_id))
            raise RuntimeError(msg)
            #assert self.constraint_id == constraint_id, ''
        self.i += 1

    def build(self):
        #float_fmt = self.model.float_fmt
        self.n = len(self.self.constraint_id)

        self.node_id = np.array(self.node_id)
        self.components = np.array(self.components)
        self.enforced_motion = np.array(self.enforced_motion)

    def write_card(self, bdf_file, size=8):
        if self.n:
            for (nid, constraint, enforced) in zip(self.node_id, self.components, self.enforced_motion):
                fields = ['SPCD', self.constraint_id, nid, constraint, enforced]
                if size == 8:
                    bdf_file.write(print_card_8(fields))
                else:
                    bdf_file.write(print_card_16(fields))
