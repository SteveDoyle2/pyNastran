from itertools import izip, count

from numpy import array, zeros, unique, searchsorted, arange

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank,
    components, components_or_blank)


class SPCD(object):
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis).::

      SPCD SID G1 C1 D1   G2 C2 D2
      SPCD 2   32 3  -2.6  5
    """
    type = 'SPC'

    def __init__(self, model):
        self.model = model
        self.n = 0

        self._comments = []
        self.constraint_id = None
        self.grid_id = []
        self.components = []

    def add(self, i, constraint_id, node_id, dofs, enforced_motion, comment):
        if i == 0:
            n = 2
        elif i == 1:
            n = 5
        else:
            raise RuntimeError('i =', i)

        assert enforced_motion != 0.0

        self._comments.append(comment)
        self.grid_id.append(node_id)
        self.components.append(dofs)
        self.enforced_motion.append(enforced_motion)

        if self.constraint_id is None:
            self.constraint_id = constraint_id
        else:
            raise RuntimeError('self.constraint_id == constraint_id; constraint_id=%r expected; found=%r' % (self.constraint_id. constraint_id))
            assert self.constraint_id == constraint_id, ''

    def build(self):
        float_fmt = self.model.float
        self.n = len(self.self.constraint_id)

        self.grid_id = array(self.grid_id)
        self.components = array(self.components)
        self.enforced_motion = array(self.enforced_motion)

    def write_bdf(self, f, size=8):
        if self.n:
            fields = ['SPCD', self.constraint_id]
            for (nid, constraint, enforced) in izip(self.gids, self.components, self.enforced_motion):
                fields += [nid, constraint, enforced]
            f.write(card)