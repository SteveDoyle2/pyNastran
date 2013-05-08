## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import os

#print "f = ",op4.__file__
from numpy import ndarray
import unittest
from pyNastran.op4.op4 import OP4
#from pyNastran.op4.cop4 import OP4 as cOP4

import pyNastran.op4.test
op4Path = pyNastran.op4.test.__path__[0]
#print(op4Path)


def pass_test1():
    fh = cOP4(os.path.abspath('mat_b_dn.op4'), 'r')
    fh.print_header()
    #print fh.nmat = 9

    # crash
    a, b, c = fh.Load(nmat=3, skip=0)
    print(a)
    print(b)
    print(c)


def failed_test1():
    fh = cOP4('mat_b_dn.op4', 'r')
    fh.print_header()
    #print fh.nmat = 9

    a, b, c = fh.Load(nmat=3, skip=0)
    print(a)
    print(b)
    print(c)


def pass_test2():
    fh = cOP4(os.path.abspath('mat_b_dn.op4'), 'r')
    #print fh.nmat = 9

    # crash with "unnamed is sparse, skipping for now"
    (a, b, c, d, f, g, h, i) = fh.Load(nmat=9, skip=0)
    print(a)
    print(b)
    print(c)


def failed_test2():
    fh = cOP4('mat_b_dn.op4', 'r')
    #print fh.nmat = 9

    # ValueError:  need more than 8 values to unpack
    (a, b, c, d, f, g, h, i) = fh.Load(nmat=9, skip=0)
    print(a)
    print(b)
    print(c)


class TestOP4(unittest.TestCase):
    def test_op4_binary(self):
        for fname in ['mat_b_dn.op4',
                      'mat_b_s1.op4',
                      'mat_b_s2.op4',
                      ]:
            op4 = OP4()

            matrices = op4.read_op4(os.path.join(op4Path, fname))
            for name, (form, matrix) in sorted(matrices.items()):
                #print("name = %s" % (name))
                if isinstance(matrix, ndarray):
                    print(matrix)
                else:
                    #print(matrix.todense())
                    print(matrix)

    def test_op4_ascii(self):
        for fname in ['mat_t_dn.op4',
                      'mat_t_s1.op4',
                      'mat_t_s2.op4',
                      ]:
            op4 = OP4()
            matrices = op4.read_op4(os.path.join(op4Path, fname))
            for name, (form, matrix) in sorted(matrices.items()):
                #print("name = %s" % name)
                if isinstance(matrix, ndarray):
                    print(matrix)
                else:
                    #print(matrix.todense())
                    print(matrix)

if __name__ == '__main__':
    #failed_test1()
    #print "*********"
    #failed_test2()

    #pass_test1()
    #pass_test2()
    unittest.main()
