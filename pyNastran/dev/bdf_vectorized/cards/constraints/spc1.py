from io import StringIO
from collections import defaultdict
#from itertools import count

import numpy as np
#from numpy import array

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    parse_components)

def get_spc1_constraint(card):
    constraint_id = integer(card, 1, 'constraint_id')
    dofs = parse_components(card, 2, 'constraints')  # 246 = y; dx, dz dir

    node_ids = card.fields(3)
    node_ids = expand_thru(node_ids)

    assert isinstance(constraint_id, int), constraint_id
    return constraint_id, dofs, node_ids


class SPC1:
    """
    +------+-----+------+--------+--------+--------+--------+--------+-----+
    | SPC1 | SID |  C   |   G1   | G2     |   G3   |   G4   |   G5   | G6  |
    +------+-----+------+--------+--------+--------+--------+--------+-----+
    |      | G7  |  G8  |   G9   |  etc.  |        |        |        |     |
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
        self.n = 0

    def add(self, constraint_id, dofs, node_ids, comment):
        #if comment:
            # self.comment = comment
        assert isinstance(constraint_id, int), constraint_id
        self.constraint_id = constraint_id
        self.components[dofs] += node_ids

    def add_card(self, card: BDFCard, comment: str=''):
        #if comment:
            # self.comment = comment
        constraint_id = integer(card, 1, 'conid')
        dofs = parse_components(card, 2, 'constraints')  # 246 = y; dx, dz dir
        node_ids = card.fields(3)

        assert isinstance(constraint_id, int), constraint_id
        self.constraint_id = constraint_id
        self.components[dofs] += node_ids
        self.n += 1

    #def allocate(self, card_count):
        #pass

    def build(self):
        for comp, nodes_lists in self.components.items():
            nodes2 = []
            for nodes in nodes_lists:
                nodes2 += nodes
            nodes2.sort()
            self.components[comp] = np.array(nodes2, dtype='int32')

    def update(self, maps):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        nid_map = maps['node']
        components = {}
        for dof, nids in self.components.items():
            components[dof] = [nid_map[nid] for nid in nids]
        self.components = components
        # TODO: constraint_map...

    def write_card(self, bdf_file, size=8):
        for comp, nodes in self.components.items():
            nodes = np.array(nodes, dtype='int32')
            nodes = np.unique(nodes)
            dnid = nodes.max() - nodes.min() + 1
            nnodes = len(nodes)
            #print('dnid=%s nnodes=%s' % (dnid, nnodes))
            if dnid == len(nodes):
                card = ['SPC1', self.constraint_id, comp, nodes.min(), 'THRU', nodes.max()]
            else:
                card = ['SPC1', self.constraint_id, comp] + list(nodes)

            if size == 8:
                bdf_file.write(print_card_8(card))
            else:
                bdf_file.write(print_card_16(card))

    def __repr__(self):
        f = StringIO()
        self.write_card(f)
        return f.getvalue().rstrip()
