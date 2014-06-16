import os

from numpy import ndarray, eye, array_equal, complex64, complex128, zeros
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
            for name, (form, matrix) in sorted(matrices.iteritems()):
                #print("name = %s" % (name))
                if isinstance(matrix, ndarray):
                    pass
                    #print(matrix)
                else:
                    #print(matrix.todense())
                    pass
                    #print(matrix)

    def test_op4_ascii(self):
        for fname in ['mat_t_dn.op4',
                      'mat_t_s1.op4',
                      'mat_t_s2.op4',
                      ]:
            op4 = OP4()
            matrices = op4.read_op4(os.path.join(op4Path, fname))
            for name, (form, matrix) in sorted(matrices.iteritems()):
                #print("name = %s" % name)
                if isinstance(matrix, ndarray):
                    #print(matrix)
                    pass
                else:
                    pass
                    #print(matrix.todense())
                    #print(matrix)

    def test_EYE10(self):
        for fname in ['mat_t_dn.op4',
                      'mat_t_s1.op4',
                      'mat_t_s2.op4',
                      ]:
            op4 = OP4()
            matrices = op4.read_op4(os.path.join(op4Path, fname))
            (form, A) = matrices['EYE10']
            self.assertEquals(form, 6)  # form=6 -> Symmetric
            if 's' in fname:  # sparse
                self.assertTrue(array_equal(A.row, range(10)))
                self.assertTrue(array_equal(A.col, range(10)))
                self.assertTrue(array_equal(A.data, [1] * 10))
            else: # real
                E = eye(10)
                self.assertTrue(array_equal(A, E))

    def test_EYE5CD(self):
        for fname in ['mat_t_dn.op4',
                      'mat_t_s1.op4',
                      'mat_t_s2.op4',
                      ]:
            op4 = OP4()
            matrices = op4.read_op4(os.path.join(op4Path, fname))
            (form, A) = matrices['EYE5CD']
            self.assertEquals(form, 6)  # form=6 -> Symmetric
            if 's' in fname:  # sparse
                self.assertTrue(array_equal(A.row, range(5)))
                self.assertTrue(array_equal(A.col, range(5)))
                self.assertTrue(array_equal(A.data, [-1+1j] * 5))
            else: # real
                E = -eye(5) + 1j*eye(5)
                self.assertTrue(array_equal(A, E))

    def test_NULL(self):
        for fname in ['mat_t_dn.op4',
                      'mat_t_s1.op4',
                      'mat_t_s2.op4',
                      ]:
            #print('-------%s-------' % fname)
            op4 = OP4()
            matrices = op4.read_op4(os.path.join(op4Path, fname))
            (form, A) = matrices['NULL']
            self.assertEquals(form, 6)  # form=6 -> Symmetric
            #print A.shape

            # kind of strange that the NULL matrix is dense...
            #if 's' in fname:  # sparse
            if 0:
                self.assertTrue(array_equal(A.row, range(3)))
                self.assertTrue(array_equal(A.col, range(3)))
                msg = 'fname=%s NULL sparse matrix error' % fname
                self.assertTrue(array_equal(A.data, [0] * 3), msg)
            else: # real
                E = zeros((3,3))
                msg = 'fname=%s NULL dense matrix error' % fname
                self.assertTrue(array_equal(A, E))
            del A

if __name__ == '__main__':
    #failed_test1()
    #print "*********"
    #failed_test2()

    #pass_test1()
    #pass_test2()
    unittest.main()
