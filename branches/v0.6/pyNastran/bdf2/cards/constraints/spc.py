from collections import defaultdict
from itertools import izip, count

from numpy import array, zeros, unique, searchsorted, arange

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank,
    components, components_or_blank)


def get_spc_constraint(card, i):
    if i == 0:
        constraint_id = integer(card, 1, 'sid')
        node_id = integer(card, 2, 'G1')
        dofs = components_or_blank(card, 3, 'C1', 0)
        enforced_motion = double_or_blank(card, 4, 'D1', 0.0)
    elif i == 1:
        constraint_id = integer(card, 1, 'sid')
        node_id = integer_or_blank(card, 5, 'G2')
        dofs = components_or_blank(card, 6, 'C2', 0)
        enforced_motion = double_or_blank(card, 7, 'D2', 0.0)
    else:
        raise RuntimeError('i =', i)

    return constraint_id, node_id, dofs, enforced_motion
    
class SPC(object):
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis).::

      SPC SID G1 C1 D1   G2 C2 D2
      SPC 2   32 3  -2.6  5
    """
    type = 'SPC'

    def __init__(self, model):
        self.model = model
        self.n = 0

        self._comments = []
        self.constraint_id = None
        self.grid_id = None
        self.components = defaultdict(list)

    def add(self, constraint_id, node_id, dofs, enforced_motion, comment):
        if node_id is None:
            return
        assert enforced_motion == 0.0

        self._comments.append(comment)
        self.components[dofs].append(node_id)

        if self.constraint_id is None:
            self.constraint_id = constraint_id
        elif self.constraint_id != constraint_id:
            raise RuntimeError('self.constraint_id == constraint_id; constraint_id=%r expected; found=%r' % (self.constraint_id. constraint_id))

    def build(self):
        float_fmt = self.model.float
        self.n = len(self.components)
        if self.n:
            self.grid_id = array(self.grid_id)
            for dof, nodes in self.components.iteritems():
                self.components[dof] = array(nodes)

    def write_bdf(self, f, size=8):
        if self.n:
            for dof, node_ids in sorted(self.components.iteritems()):
                card = ['SPC', self.constraint_id]
                for node_id in node_ids:
                    card += [node_id, dof, 0.0]
                f.write(print_card(card))