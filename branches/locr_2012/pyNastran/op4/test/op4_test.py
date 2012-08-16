import os
import sys

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


class OP4_Test(unittest.TestCase):
    def test_op4_binary(self):
        for fname in ['mat_b_dn.op4',
                      'mat_b_s1.op4',
                      'mat_b_s2.op4',
                      ]:
            op4 = OP4()

            matrices = op4.readOP4(os.path.join(op4Path, fname))
            for name, (form, matrix) in sorted(matrices.items()):
                print("name = %s" % (name))
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
            matrices = op4.readOP4(os.path.join(op4Path, fname))
            for name, (form, matrix) in sorted(matrices.items()):
                print("name = %s" % (name))
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
