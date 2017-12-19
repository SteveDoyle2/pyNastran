from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef, TableData


class CoordinateSystem(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.cord1c = CORD1C(self._h5n, self)
        self.cord1r = CORD1R(self._h5n, self)
        self.cord1s = CORD1S(self._h5n, self)
        self.cord2c = CORD2C(self._h5n, self)
        self.cord2r = CORD2R(self._h5n, self)
        self.cord2s = CORD2S(self._h5n, self)
        self.cord3g = CORD3G(self._h5n, self)
        self.cord3r = CORD3R(self._h5n, self)
        # self.transformation = TRANSFORMATION(self._h5n, self)

    def path(self):
        return self._input.path() + ['COORDINATE_SYSTEM']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass


########################################################################################################################


class CORD1C(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD1C')

########################################################################################################################


class CORD1R(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD1R')

########################################################################################################################

class CORD1S(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD1S')

########################################################################################################################


def _cord2c_from_bdf(card):
    data = TableData()
    e1 = card.e1
    e2 = card.e2
    e3 = card.e3
    data.data = [[card.cid, card.rid, e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], e3[0], e3[1], e3[2]]]
    return data


class CORD2C(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD2C')

    from_bdf = staticmethod(_cord2c_from_bdf)

########################################################################################################################


class CORD2R(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD2R')

    from_bdf = staticmethod(_cord2c_from_bdf)

########################################################################################################################


class CORD2S(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD2S')

    from_bdf = staticmethod(_cord2c_from_bdf)

########################################################################################################################


class CORD3G(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD3G')

########################################################################################################################


class CORD3R(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD3R')

########################################################################################################################


# class TRANSFORMATION(CardTable):
#     table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/TRANSFORMATION/IDENTITY')

########################################################################################################################
