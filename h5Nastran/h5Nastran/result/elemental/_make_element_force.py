

_data = """34 CBAR TRUE
100 CBARS TRUE
2 CBEAM TRUE
184 CBEAM3 TRUE
69 _BEND TRUE
102 CBUSH TRUE
35 CCONEAX FALSE
20 CDAMP1 TRUE
21 CDAMP2 TRUE
22 CDAMP3 TRUE
23 CDAMP4 TRUE
55 CDUM3 TRUE
56 CDUM4 TRUE
57 CDUM5 TRUE
58 CDUM6 TRUE
59 CDUM7 TRUE
60 CDUM8 TRUE
61 CDUM9 TRUE
11 CELAS1 TRUE
12 CELAS2 TRUE
13 CELAS3 TRUE
14 CELAS4 TRUE
38 CGAP FALSE
10 _CONROD TRUE
33 CQUAD4 TRUE
95 CQUAD4_COMP FALSE
144 CQUAD4_CN TRUE
64 CQUAD8 TRUE
96 CQUAD8_COMP TRUE
82 CQUADR TRUE
235 CQUADR_NL TRUE
1 CROD TRUE
4 CSHEAR TRUE
74 CTRIA3 TRUE
97 CTRIA3_COMP FALSE
75 CTRIA6 TRUE
98 CTRIA6_COMP FALSE
70 CTRIAR TRUE
236 CTRIAR_NL TRUE
3 CTUBE TRUE
24 CVISC TRUE
118 CWELDP TRUE
117 CWELDC TRUE
200 CWELD TRUE
189 _VUQUAD TRUE
190 _VUTRIA TRUE
191 _VUBEAM TRUE"""

_data = _data.split('\n')

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
    _lines = """########################################################################################################################
    
    
class %s(ResultTable):
    result_type = 'ELEMENT FORCES %d %s %s'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/%s', result_type)
        

""" % \
             (name.upper(), num, name.upper(), real_complex, name.upper())

    lines[0] += _lines
    _data_.append('        self.%s = %s(self._h5n, self)' % (name.lower(), name.upper()))


for i in range(len(_data)):
    num, name, has_complex = _data[i].split()
    num = int(num)
    has_complex = bool(has_complex)

    make_class(name[1:], num, 'REAL')

    if has_complex is True:
        make_class(name[1:] + '_CPLX', num, 'COMPLEX')


lines = lines[0].replace('__data__', '\n'.join(_data_))

with open('_element_forces.txt', 'w') as f:
    f.write(lines)
