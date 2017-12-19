from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef, TableData


class Load(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.force = FORCE(self._h5n, self)
        self.moment = MOMENT(self._h5n, self)

    def path(self):
        return self._input.path() + ['LOAD']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass

########################################################################################################################


class FORCE(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/FORCE')

    @staticmethod
    def from_bdf(card):
        data = []
        for _ in card:
            data.append([_.sid, _.node, _.cid, _.mag, _.xyz])
        return TableData(list(sorted(data, key=lambda x: x[1])))

########################################################################################################################


class MOMENT(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/MOMENT')

    @staticmethod
    def from_bdf(card):
        data = []
        for _ in card:
            data.append([_.sid, _.node, _.cid, _.mag, _.xyz])
        return TableData(list(sorted(data, key=lambda x: x[1])))
