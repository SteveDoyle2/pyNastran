

_data = """47 CAXIF2 TRUE
48 CAXIF3 TRUE
49 CAXIF4 TRUE
241 CAXISYM FALSE
34 CBAR TRUE
100 CBARS TRUE
238 CBAR_NL FALSE
2 CBEAM TRUE
94 CBEAM FALSE
239 CBEAM FALSE
184 CBEAM3 TRUE
69 CBEND TRUE
102 CBUSH TRUE
40 CBUSH1D FALSE
35 CCONEAX FALSE
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
86 CGAP FALSE
67 CHEXA TRUE
93 CHEXA FALSE
202 CHEXAFD FALSE
207 CHEXAFD 207 FALSE
65 CIFHEX FALSE
66 CIFPENT FALSE
73 CIFQDX FALSE
63 CIFQUAD FALSE
10 _CONROD TRUE
92 _CONROD FALSE
68 CPENTA TRUE
91 CPENTA FALSE
204 CPENTAFD FALSE
209 CPENTAFD FALSE
33 CQUAD4 TRUE
90 CQUAD4 FALSE
95 CQUAD4 TRUE
144 CQUAD4 TRUE
64 CQUAD8 TRUE
96 CQUAD8 TRUE
201 CQUADFD FALSE
208 CQUADFD FALSE
82 CQUADR TRUE
172 CQUADR FALSE
232 CQUADR TRUE
18 CQUADX TRUE
214 CQUADXFD FALSE
215 CQUADXFD FALSE
1 CROD TRUE
89 CROD FALSE
4 CSHEAR TRUE
50 CSLOT3 TRUE
51 CSLOT4 TRUE
39 CTETRA TRUE
85 CTETRA FALSE
205 CTETRAFD FALSE
210 CTETRAFD FALSE
74 CTRIA3 TRUE
97 CTRIA3 TRUE
88 CTRIA3 FALSE
75 CTRIA6 TRUE
98 CTRIA6 TRUE
206 CTRIAFD FALSE
211 CTRIAFD FALSE
70 CTRIAR TRUE
173 CTRIAR FALSE
233 CTRIAR TRUE
17 CTRIAX TRUE
53 CTRIAX6 TRUE
212 CTRIAXFD FALSE
213 CTRIAXFD FALSE
3 CTUBE TRUE
87 CTUBE FALSE
118 CWELDP TRUE
117 CWELDC TRUE
200 CWELD TRUE"""

_data = _data.split('\n')

lines = """from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from ..result_table import ResultTable, TableDef


class ElementStress(object):
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
    _lines = """########################################################################################################################
    
    
class %s(ResultTable):
    result_type = 'ELEMENT STRESSES %d %s %s'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRESS/%s', result_type)
        

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

with open('_element_stresses.txt', 'w') as f:
    f.write(lines)
