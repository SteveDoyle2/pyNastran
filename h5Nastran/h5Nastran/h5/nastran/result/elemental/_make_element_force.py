

_data = 'BAR 34;BARS 100;BARS_CPLX 100;BAR_CPLX 34;BEAM 2;BEAM3 184;BEAM3_CPLX 184;BEAM_CPLX 2;BEND 69;BEND_CPLX 69;BUSH 102;BUSH1D 40;BUSH_CPLX 102;CONE 35;CONROD 10;CONROD_CPLX 10;CONV 110;CONVM 111;DAMP1 20;DAMP1_CPLX 20;DAMP2 21;DAMP2_CPLX 21;DAMP3 22;DAMP3_CPLX 22;DAMP4 23;DAMP4_CPLX 23;ELAS1 11;ELAS1_CPLX 11;ELAS2 12;ELAS2_CPLX 12;ELAS3 13;ELAS3_CPLX 13;ELAS4 14;ELAS4_CPLX 14;FAST 126;FAST_CPLX 126;GAP 38;GRAD_FLUX -1;HBDYE 107;HBDYG 108;HBDYP 109;QUAD4 33;QUAD4_CN 144;QUAD4_CN_CPLX 144;QUAD4_COMP 95;QUAD4_CPLX 33;QUAD8 64;QUAD8_COMP 96;QUAD8_CPLX 64;QUADR 82;QUADR_CPLX 82;RAC2D 60;RAC2D_CPLX 60;RAC3D 61;RAC3D_CPLX 61;RADBC 115;RADINT 155;ROD 1;ROD_CPLX 1;SEAM 119;SEAM_CPLX 119;SHEAR 4;SHEAR_CPLX 4;TRIA3 74;TRIA3_COMP 97;TRIA3_CPLX 74;TRIA6 75;TRIA6_COMP 98;TRIA6_CPLX 75;TRIAR 70;TRIAR_CPLX 70;TUBE 3;TUBE_CPLX 3;VISC 24;VISC_CPLX 24;WELD 200;WELDC 117;WELDC_CPLX 117;WELDP 118;WELDP_CPLX 118;WELD_CPLX 200'

_data = [_.split(' ') for _ in _data.split(';')]

lines = """from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from ..result_table import ResultTable, TableDef


class ElementForce(object):
    def __init__(self, h5n, elemental):
        self._h5n = h5n
        self._elemental = elemental

__data__

    def path(self):
        return self._elemental.path() + ['ELEMENT_FORCE']
        

"""

_data_ = []

lines = [lines]


def make_class(name, num, real_complex):

    result_name = name.upper().replace('_CPLX', '').replace('_CN', 'C').replace('_COMP', 'LC')

    _lines = """########################################################################################################################
    
    
class %s(ResultTable):
    result_type = 'ELEMENT FORCES %d %s %s'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/%s', result_type)
        

""" % \
             (name.upper(), num, result_name, real_complex, name.upper())

    lines[0] += _lines
    _data_.append('        self.%s = %s(self._h5n, self)' % (name.lower(), name.upper()))


for i in range(len(_data)):
    name, num = _data[i]
    num = int(num)
    if num == -1:
        _data_.append('        self.%s = None  # skipping for now' % name.lower())
        continue

    if '_CPLX' in name:
        result_type = 'COMPLEX'
    elif '_RR' in name:
        result_type = 'RANDOM'
    else:
        result_type = 'REAL'

    make_class(name, num, result_type)


lines = lines[0].replace('__data__', '\n'.join(_data_))

with open('_element_forces.txt', 'w') as f:
    f.write(lines)
