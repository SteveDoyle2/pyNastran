from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef, TableData


class Node(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.grid = GRID(self._h5n, self)

    def path(self):
        return self._input.path() + ['NODE']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass

########################################################################################################################


class GRID(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/NODE/GRID')

    @staticmethod
    def from_bdf(card):
        data = [card.nid, card.cp, card.xyz, card.cd, card.ps, card.seid]
        return TableData([data])

########################################################################################################################

