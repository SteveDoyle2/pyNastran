from __future__ import print_function
import os
import unittest
from six import iteritems, PY2

import numpy as np
from numpy import ones, reshape, arange
from numpy import ndarray, eye, array_equal, zeros
from pyNastran.op4.op4 import OP4, read_op4

import pyNastran.op4.test
op4_path = pyNastran.op4.test.__path__[0]


class TestOP4(unittest.TestCase):

    @staticmethod
    def test_op4_binary():
        fnames = [
            'mat_b_dn.op4',
            'mat_b_s1.op4',
            'mat_b_s2.op4',
        ]
        for fname in fnames:
            matrices = read_op4(os.path.join(op4_path, fname))
            for name, (form, matrix) in sorted(iteritems(matrices)):
                #print("name = %s" % (name))
                if isinstance(matrix, ndarray):
                    pass
                    #print(matrix)
                else:
                    #print(matrix.todense())
                    pass
                    #print(matrix)

    @staticmethod
    def test_op4_ascii():
        fnames = [
            'mat_t_dn.op4',
            'mat_t_s1.op4',
            'mat_t_s2.op4',
        ]
        for fname in fnames:
            op4 = OP4()
            matrices = op4.read_op4(os.path.join(op4_path, fname))
            for name, (form, matrix) in sorted(iteritems(matrices)):
                #print("name = %s" % name)
                if isinstance(matrix, ndarray):
                    #print(matrix)
                    pass
                else:
                    pass
                    #print(matrix.todense())
                    #print(matrix)

    def test_eye10(self):
        fnames = [
            'mat_t_dn.op4',
            'mat_t_s1.op4',
            'mat_t_s2.op4',
        ]
        for fname in fnames:
            op4 = OP4()
            matrices = op4.read_op4(os.path.join(op4_path, fname))
            (form, A) = matrices['EYE10']
            self.assertEqual(form, 6)  # form=6 -> Symmetric
            if 's' in fname:  # sparse
                self.assertTrue(array_equal(A.row, range(10)))
                self.assertTrue(array_equal(A.col, range(10)))
                self.assertTrue(array_equal(A.data, [1] * 10))
            else: # real
                E = eye(10)
                self.assertTrue(array_equal(A, E))

    def test_eye5cd(self):
        fnames = [
            'mat_t_dn.op4',
            'mat_t_s1.op4',
            'mat_t_s2.op4',
        ]
        for fname in fnames:
            op4 = OP4()
            matrices = op4.read_op4(os.path.join(op4_path, fname))
            (form, A) = matrices['EYE5CD']
            self.assertEqual(form, 6)  # form=6 -> Symmetric
            if 's' in fname:  # sparse
                self.assertTrue(array_equal(A.row, range(5)))
                self.assertTrue(array_equal(A.col, range(5)))
                self.assertTrue(array_equal(A.data, [-1+1j] * 5))
            else: # real
                E = -eye(5) + 1j*eye(5)
                self.assertTrue(array_equal(A, E))

    def test_null(self):
        fnames = [
            'mat_t_dn.op4',
            'mat_t_s1.op4',
            'mat_t_s2.op4',
        ]
        for fname in fnames:
            #print('-------%s-------' % fname)
            op4 = OP4()
            matrices = op4.read_op4(os.path.join(op4_path, fname))
            (form, A) = matrices['NULL']
            self.assertEqual(form, 6)  # form=6 -> Symmetric
            #print A.shape

            # kind of strange that the NULL matrix is dense...
            #if 's' in fname:  # sparse
            if 0:
                self.assertTrue(array_equal(A.row, range(3)))
                self.assertTrue(array_equal(A.col, range(3)))
                msg = 'fname=%s NULL sparse matrix error' % fname
                self.assertTrue(array_equal(A.data, [0] * 3), msg)
            else: # real
                E = zeros((3, 3))
                msg = 'fname=%s NULL dense matrix error' % fname
                self.assertTrue(array_equal(A, E))
            del A

    def test_bad_inputs_1(self):
        op4 = OP4()
        op4_filename = os.path.join(op4_path, 'bad_inputs.op4')
        form1 = 1
        A1 = ones((3, 3), dtype='float64')
        matrices = {
            'A1': (form1, A1),
            'A2': ('bad', A1),
            'A3': (form1, 'bad'),
        }
        with self.assertRaises(ValueError):
            op4.write_op4(op4_filename, matrices, name_order=None, precision='default_bad',
                          is_binary=False)
        with self.assertRaises(ValueError):
            op4.write_op4(op4_filename, matrices, name_order='A1', precision='default',
                          is_binary='bad')
        with self.assertRaises(ValueError):
            op4.write_op4(op4_filename, matrices, name_order='A2', precision='default',
                          is_binary=True)
        with self.assertRaises(NotImplementedError):
            op4.write_op4(op4_filename, matrices, name_order='A3', precision='default',
                          is_binary=True)

        # now lets write the op4, so we can test bad reading
        op4.write_op4(op4_filename, matrices, name_order='A1', precision='default',
                      is_binary=False)
        with self.assertRaises(ValueError):
            matrices2 = op4.read_op4(op4_filename, precision='default_bad')
        with self.assertRaises(IOError):
            matrices2 = op4.read_op4('op4_filename', precision='default')

        # now the inputs are valid, so this works
        matrices2 = op4.read_op4(op4_filename, precision='default')

    @staticmethod
    def test_file_obj_ascii():
        op4 = OP4()
        form1 = 1
        A1 = ones((3, 3), dtype='float64')
        matrices = {
            'A1': (form1, A1),
        }
        if PY2:
            wb = 'wb'
        else:
            wb = 'w'
        with open(os.path.join(op4_path, 'file_ascii.op4'), wb) as op4_file:
            op4.write_op4(op4_file, matrices, name_order='A1', precision='default',
                          is_binary=False)

    @staticmethod
    def test_file_obj_binary():
        op4 = OP4()
        form1 = 1
        A1 = ones((3, 3), dtype='float64')
        matrices = {
            'A1': (form1, A1),
        }
        with open(os.path.join(op4_path, 'file_binary.op4'), 'wb') as op4_file:
            op4.write_op4(op4_file, matrices, name_order='A1', precision='default',
                          is_binary=True)

    def test_square_matrices_1(self):
        op4 = OP4()
        #matrices = op4.read_op4(os.path.join(op4_path, fname))
        form1 = 1
        form2 = 2
        form3 = 2
        A1 = np.matrix(ones((3, 3), dtype='float64'))
        A2 = reshape(arange(9, dtype='float64'), (3, 3))
        A3 = np.matrix(ones((1, 1), dtype='float32'))
        matrices = {
            'A1': (form1, A1),
            'A2': (form2, A2),
            'A3': (form3, A3),
        }

        for (is_binary, fname) in [(False, 'small_ascii.op4'), (True, 'small_binary.op4')]:
            op4_filename = os.path.join(op4_path, fname)
            op4.write_op4(op4_filename, matrices, name_order=None, precision='default',
                          is_binary=False)
            matrices2 = op4.read_op4(op4_filename, precision='default')
            (form1b, A1b) = matrices2['A1']
            (form2b, A2b) = matrices2['A2']
            self.assertEqual(form1, form1b)
            self.assertEqual(form2, form2b)

            (form1b, A1b) = matrices2['A1']
            (form2b, A2b) = matrices2['A2']
            (form3b, A3b) = matrices2['A3']
            self.assertEqual(form1, form1b)
            self.assertEqual(form2, form2b)
            self.assertEqual(form3, form3b)

            self.assertTrue(array_equal(A1, A1b))
            self.assertTrue(array_equal(A2, A2b))
            self.assertTrue(array_equal(A3, A3b))
            del A1b, A2b, A3b
            del form1b, form2b, form3b


if __name__ == '__main__':  # pragma: no cover
    on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
    if not on_rtd:
        unittest.main()
