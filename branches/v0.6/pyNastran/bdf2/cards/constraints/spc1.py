from collections import defaultdict
from itertools import izip, count

from numpy import array

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.cards.baseCard import expand_thru
from pyNastran.bdf2.bdf_interface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, components, components_or_blank)

def get_spc1_constraint(card):
    constraint_id = integer(card, 1, 'constraint_id')
    dofs = components(card, 2, 'constraints')  # 246 = y; dx, dz dir

    node_ids = card.fields(3)
    node_ids = expand_thru(node_ids)

    assert isinstance(constraint_id, int), constraint_id
    return constraint_id, dofs, node_ids


class SPC1(object):
    """
    ::

      SPC1 SID C G1 G2 G3 G4 G5 G6
      G7 G8 G9 -etc.-

      SPC1     3       246     209075  209096  209512  209513  209516
      SPC1 3 2 1 3 10 9 6 5
      2 8

      SPC1 SID C    G1 THRU G2
      SPC1 313 12456 6 THRU 32
    """
    def __init__(self, model):
        self.model = model
        self.components = defaultdict(list)

    def add(self, constraint_id, dofs, node_ids, comment):
        #if comment:
            #self._comment = comment
        assert isinstance(constraint_id, int), constraint_id
        self.constraint_id = constraint_id
        self.components[dofs].append(node_ids)

    def build(self):
        for comp, nodes_lists in self.components.iteritems():
            nodes2 = []
            for nodes in nodes_lists:
                nodes2 += nodes
            nodes2.sort()
            self.components[comp] = nodes2

    def write_bdf(self, f, size=8):
        for comp, nodes in self.components.iteritems():
            card = ['SPC1', self.constraint_id, comp] + list(nodes)
            f.write(print_card(card))