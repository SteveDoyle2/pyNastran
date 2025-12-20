from collections import defaultdict

#import numpy as np
from numpy import array

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
#from pyNastran.bdf.cards.base_card import BaseCard, expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, components_or_blank)


def get_spc_constraint(card, i):
    """parses an SPC/SPCD card"""
    if i == 0:
        constraint_id = integer(card, 1, 'sid')
        node_id = integer(card, 2, 'G1')
        dofs = components_or_blank(card, 3, 'C1', '0')
        enforced_motion = double_or_blank(card, 4, 'D1', 0.0)
    elif i == 1:
        constraint_id = integer(card, 1, 'sid')
        node_id = integer_or_blank(card, 5, 'G2')
        dofs = components_or_blank(card, 6, 'C2', '0')
        enforced_motion = double_or_blank(card, 7, 'D2', 0.0)
    else:
        raise RuntimeError('i =', i)

    return constraint_id, node_id, dofs, enforced_motion

class SPC:
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis).

     +-----+-----+----+----+------+----+----+----+
     | SPC | SID | G1 | C1 |  D1  | G2 | C2 | D2 |
     +-----+-----+----+----+------+----+----+----+
     | SPC |  2  | 32 | 3  | -2.6 |  5 |    |    |
     +-----+-----+----+----+------+----+----+----+
    """
    type = 'SPC'

    def __init__(self, model):
        self.model = model
        self.n = 0

        self._comments = []
        self.constraint_id = None
        self.grid_id = []
        self.components = defaultdict(list)
        self.value = []

    def allocate(self, card_count):
        return
        #ncards = card_count['SPC']
        #if ncards:
            #self.n = ncards
            #print('ngrid=%s' % self.n)
            #float_fmt = self.model.float_fmt
            #self.node_id = zeros(ncards, 'int32')
            #self.components = zeros(ncards, 'int32')
            #self.enforced_motion = zeros(ncards, float_fmt)

    def add(self, constraint_id, node_id, dofs, enforced_motion, comment):
        assert enforced_motion == 0.0

        self._comments.append(comment)
        #self.model.log.debug('dofs=%r node_id=%r' % (dofs, node_id))
        self.components[dofs].append(node_id)

        if self.constraint_id == constraint_id:
            pass
        elif self.constraint_id is None:
            self.constraint_id = constraint_id
        elif self.constraint_id != constraint_id:
            msg = 'self.constraint_id == constraint_id; constraint_id=%r expected; found=%r' % (
                self.constraint_id. constraint_id)
            raise RuntimeError(msg)
        self.n += 1

    def build(self):
        self.n = len(self.components)
        if self.n:
            self.grid_id = array(self.grid_id)
            for dof, nodes in self.components.items():
                self.components[dof] = array(nodes)

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
        if self.n:
            for dof, node_ids in sorted(self.components.items()):
                card = ['SPC', self.constraint_id]
                for node_id in node_ids:
                    card += [node_id, dof, 0.0]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))
