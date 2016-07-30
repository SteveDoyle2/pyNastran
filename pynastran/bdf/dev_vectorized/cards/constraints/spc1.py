from six import iteritems
from six.moves import StringIO
from collections import defaultdict
#from itertools import count

#from numpy import array

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.dev_vectorized.bdf_interface.assign_type import (integer,
    components)

def get_spc1_constraint(card):
    constraint_id = integer(card, 1, 'constraint_id')
    dofs = components(card, 2, 'constraints')  # 246 = y; dx, dz dir

    node_ids = card.fields(3)
    node_ids = expand_thru(node_ids)

    assert isinstance(constraint_id, int), constraint_id
    return constraint_id, dofs, node_ids


class SPC1(object):
    """
    +------+-----+------+--------+--------+--------+--------+--------+-----+
    | SPC1 | SID |  C   |   G1   | G2     |   G3   |   G4   |   G5   | G6  |
    +------+-----+------+--------+--------+--------+--------+--------+-----+
    |      | G7  |  G8  |   G9   | -etc.- |        |        |        |     |
    +------+-----+------+--------+--------+--------+--------+--------+-----+

    +------+-----+------+--------+--------+--------+--------+--------+-----+
    | SPC1 |  3  | 246  | 209075 | 209096 | 209512 | 209513 | 209516 |     |
    +------+-----+------+--------+--------+--------+--------+--------+-----+
    | SPC1 |  3  |   2  |   1    |    3   |   10   |   9    |   6    |  5  |
    +------+-----+------+--------+--------+--------+--------+--------+-----+
    |      |  2  |   8  |        |        |        |        |        |     |
    +------+-----+------+--------+--------+--------+--------+--------+-----+

    +------+-----+-------+-------+--------+--------+--------+--------+-----+
    | SPC1 | SID | C     |   G1  |  THRU  |   G2   |        |        |     |
    +------+-----+-------+-------+--------+--------+--------+--------+-----+
    | SPC1 | 313 | 12456 |   6   |  THRU  |   32   |        |        |     |
    +------+-----+-------+-------+--------+--------+--------+--------+-----+
    """
    type = 'SPC1'
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
        for comp, nodes_lists in iteritems(self.components):
            nodes2 = []
            for nodes in nodes_lists:
                nodes2 += nodes
            nodes2.sort()
            self.components[comp] = nodes2

    def write_card(self, f, size=8):
        for comp, nodes in iteritems(self.components):
            card = ['SPC1', self.constraint_id, comp] + list(nodes)
            if size == 8:
                f.write(print_card_8(card))
            else:
                f.write(print_card_16(card))

    def __repr__(self):
        f = StringIO()
        self.write_card(f)
        return f.getvalue().rstrip()
