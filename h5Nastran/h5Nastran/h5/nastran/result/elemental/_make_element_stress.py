
_data = 'AXIF2 47;AXIF2_CPLX 47;AXIF3 48;AXIF3_CPLX 48;AXIF4 49;AXIF4_CPLX 49;AXISYM 241;BAR 34;BARS 100;BARS_CPLX 100;BAR_CPLX 34;BAR_NL 240;BAR_RR 34;BEAM 2;BEAM3 184;BEAM3_CPLX 184;BEAM_CPLX 2;BEAM_NL 94;BEAM_RR 2;BUSH 102;BUSH1D 40;BUSH1D_CPLX 40;BUSH1D_RR 40;BUSH_CPLX 102;BUSH_NL 226;CONE 35;CONROD 10;CONROD_CPLX 10;CONROD_NL 92;CONROD_RR 10;ELAS1 11;ELAS1_CPLX 11;ELAS1_NL 224;ELAS2 12;ELAS2_CPLX 12;ELAS3 13;ELAS3_CPLX 13;ELAS3_NL 225;ELEM_COMP -1;EXTREME_FIBRE -1;EXTREME_FIBRE_CPLX -1;FAST 126;FAST_CPLX 126;GAP 38;GAP_NL 86;HEXA 67;HEXA20_27FDNL 207;HEXA20_8FDNL 202;HEXA20_FD 163;HEXA_CPLX 67;HEXA_FD 140;HEXA_FDNL 218;HEXA_NL 93;IFHEXA 65;IFPENTA 66;PENTA 68;PENTA15_21FDNL 209;PENTA15_6FDNL 204;PENTA15_FD 165;PENTA_FDNL 220;PENTA_NL 91;QUAD4 33;QUAD4_COMP 95;QUAD4_COMP_CPLX 95;QUAD4_CPLX 33;QUAD4_FD 139;QUAD4_FDNL 201;QUAD4_FD_CPLX 139;QUAD4_NL 90;QUAD8 64;QUAD8_4FDNL 219;QUAD8_9FDNL 208;QUAD8_COMP 96;QUAD8_COMP_CPLX 96;QUAD8_CPLX 64;QUAD8_FD 164;QUAD8_FD_CPLX 164;QUADR 82;QUADR_COMP 232;QUADR_COMP_CPLX 232;QUADR_CPLX 82;QUADR_FD -1;QUADR_FD_CPLX -1;QUADR_NL 172;QUADX4_FD 170;QUADX4_FDNL 214;QUADX4_FD_CPLX 170;QUADX8_4FDNL 223;QUADX8_9FDNL 215;QUADX8_FD 171;QUADX8_FD_CPLX 171;QUAD_CN 144;QUAD_CN_CPLX 144;RAC2D 60;RAC3D 61;ROD 1;ROD_CPLX 1;ROD_NL 89;ROD_RR 1;SEAM 119;SEAMP 159;SEAMP_CPLX 159;SHEAR 4;SHEAR_CPLX 4;SHEAR_RR 4;SLOT3 50;SLOT3_CPLX 50;SLOT4 51;SLOT4_CPLX 51;TETRA 39;TETRA10_4FDNL 221;TETRA10_5FDNL 210;TETRA10_FD 166;TETRA4_FDNL 216;TETRA_CPLX 39;TETRA_FD 161;TETRA_FDNL 205;TETRA_NL 85;TRIA3 74;TRIA3_1FDNL 206;TRIA3_3FDNL 217;TRIA3_COMP 97;TRIA3_COMP_CPLX 97;TRIA3_CPLX 74;TRIA3_FD 162;TRIA3_FD_CPLX 162;TRIA3_NL 88;TRIA6 75;TRIA6_COMP 98;TRIA6_COMP_CPLX 98;TRIA6_CPLX 75;TRIA6_FD 167;TRIA6_FDNL 211;TRIA6_FD_CPLX 167;TRIAR 70;TRIAR_1FD -1;TRIAR_1FD_CPLX -1;TRIAR_4FD -1;TRIAR_4FD_CPLX -1;TRIAR_COMP 233;TRIAR_COMP_CPLX 233;TRIAR_CPLX 70;TRIAR_NL 173;TRIAX3_1FDNL 212;TRIAX3_3FDNL 222;TRIAX3_FD 168;TRIAX3_FD_CPLX 168;TRIAX6 53;TRIAX6_CPLX 53;TRIAX6_FD 169;TRIAX6_FDNL 213;TRIAX6_FD_CPLX 169;TUBE 3;TUBE_CPLX 3;TUBE_NL 87;TUBE_RR 3;VISC_CPLX 24;VISC_RR 24;WELD 200;WELD_CPLX 200;WELDC 117;WELDC_CPLX 117;WELDP 118;WELDP_CPLX 118'

_data = [_.split(' ') for _ in _data.split(';')]

lines = """from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from ..result_table import ResultTable, TableDef


class Stress(object):
    def __init__(self, h5n, elemental):
        self._h5n = h5n
        self._elemental = elemental

__data__

    def path(self):
        return self._elemental.path() + ['STRESS']


"""

_data_ = []

lines = [lines]


def make_class(name, num, real_complex):

    result_name = name.upper().replace('_CPLX', '').replace('_CN', 'C').replace('_COMP', 'LC')

    _lines = """########################################################################################################################


class %s(ResultTable):
    result_type = 'ELEMENT STRESSES %d %s %s'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRESS/%s', result_type)


""" % \
             (name.upper(), num, result_name, real_complex, name.upper())

    lines[0] += _lines
    _data_.append('        self.%s = %s(self._h5n, self)' % (name.lower(), name.upper()))


for i in range(len(_data)):
    name, num = _data[i]
    num = int(num)
    if num == -1:
        _data_.append('        # self.%s = None  # skipping for now' % name.lower())
        continue

    if '_CPLX' in name:
        result_type = 'COMPLEX'
    elif '_RR' in name:
        result_type = 'RANDOM'
    else:
        result_type = 'REAL'

    make_class(name, num, result_type)


lines = lines[0].replace('__data__', '\n'.join(_data_))

with open('_element_stresses.txt', 'w') as f:
    f.write(lines)
