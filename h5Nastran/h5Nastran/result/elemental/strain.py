from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from ..result_table import ResultTable, TableDef


class Strain(object):
    def __init__(self, h5n, elemental):
        self._h5n = h5n
        self._elemental = elemental

        self.bar = BAR(self._h5n, self)
        self.beam = BEAM(self._h5n, self)
        self.quad4 = QUAD4(self._h5n, self)
        self.quad_cn = QUAD4_CN(self._h5n, self)
        self.tria3 = TRIA3(self._h5n, self)

    def path(self):
        return self._elemental.path() + ['STRAIN']

########################################################################################################################


class BAR(ResultTable):
    result_type = 'ELEMENT STRAINS 34 BAR REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BAR', result_type)


########################################################################################################################


class BEAM(ResultTable):
    result_type = 'ELEMENT STRAINS 2 BEAM REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BEAM', result_type)


########################################################################################################################


class QUAD4(ResultTable):
    result_type = 'ELEMENT STRAINS 33 QUAD4 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD4', result_type)


########################################################################################################################

# TODO: QUAD4_CN strain/stress not in spec, format from QRG, guessing column id's for now

class QUAD4_CN_DEF(object):
    name = 'QUAD4_CN'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = [
        ('EID', '<i8', ()), ('TYPE', 'S4', ()),
        ('STRESS',
            [
                ('GRID', '<i8', ()), ('FD1', '<f8', ()), ('X1', '<f8', ()), ('Y1', '<f8', ()), ('XY1', '<f8', ()),
                ('A1', '<f8', ()), ('FP1', '<f8', ()), ('FM1', '<f8', ()), ('FVM1', '<f8', ()),
                ('FD2', '<f8', ()), ('X2', '<f8', ()), ('Y2', '<f8', ()), ('XY2', '<f8', ()),
                ('A2', '<f8', ()), ('FP2', '<f8', ()), ('FM2', '<f8', ()), ('FVM2', '<f8', ()),
            ],
         (5,)
         ),
        ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    subtables = []
    same_as = None


class QUAD4_CN(ResultTable):
    result_type = 'ELEMENT STRAINS 144 QUAD4C REAL'
    table_def = TableDef.create(QUAD4_CN_DEF, result_type)


########################################################################################################################


class TRIA3(ResultTable):
    result_type = 'ELEMENT STRAINS 74 TRIA3 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA3', result_type)

########################################################################################################################
