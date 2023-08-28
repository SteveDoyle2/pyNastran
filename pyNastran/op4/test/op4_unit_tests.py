"""runs various OP4 tests"""
import os
import unittest

import numpy as np
from numpy import ones, reshape, arange
from numpy import ndarray, eye, array_equal, zeros
import scipy
import scipy.sparse
from scipy.sparse import coo_matrix  # type: ignore

from pyNastran.op4.op4 import OP4, read_op4, Matrix
import pyNastran.op4.test

OP4_PATH = pyNastran.op4.test.__path__[0]
PKG_PATH = pyNastran.__path__[0]

#coo_matrix = scipy.sparse._coo.coo_matrix

class TestOP4(unittest.TestCase):
    """runs various OP4 tests"""
    def test_op4_binary(self):
        fnames = [
            'mat_b_dn.op4',
            'mat_b_s1.op4',
            'mat_b_s2.op4',
        ]
        for fname in fnames:
            op4_filename = os.path.join(OP4_PATH, fname)
            matrices = read_op4(op4_filename)

            for unused_name, matrix in sorted(matrices.items()):
                assert isinstance(matrix, Matrix), type(matrix)
                form = matrix.form
                data = matrix.data
                assert form in {1, 2, 6}, form
                #print("name = %s" % (name))
                if isinstance(data, ndarray):
                    pass
                    #print(data)
                else:
                    #print(data.toarray())
                    assert isinstance(data, coo_matrix), type(data)
                    #print(data)
                matrix.write_dmi()

    def test_op4_ascii(self):
        fnames = [
            'mat_t_dn.op4',
            'mat_t_s1.op4',
            'mat_t_s2.op4',
        ]
        for fname in fnames:
            op4_filename = os.path.join(OP4_PATH, fname)
            matrices = read_op4(op4_filename)
            for unused_name, matrix in sorted(matrices.items()):
                #print("name = %s" % name)
                data = matrix.data
                if isinstance(data, ndarray):
                    #print(data)
                    pass
                else:
                    pass
                    #print(data.toarray())
                    #print(data)
                matrix.write_dmi()

    def test_eye10(self):
        """tests the EYE10 matrices"""
        fnames = [
            'mat_t_dn.op4',
            'mat_t_s1.op4',
            'mat_t_s2.op4',
        ]
        for fname in fnames:
            op4_filename = os.path.join(OP4_PATH, fname)
            matrices = read_op4(op4_filename, debug=False)
            A = matrices['EYE10']
            mat = A.data
            form = A.form
            self.assertEqual(form, 6)  # form=6 -> Symmetric
            if 's' in fname:  # sparse
                self.assertTrue(array_equal(mat.row, range(10)))
                self.assertTrue(array_equal(mat.col, range(10)))
                self.assertTrue(array_equal(mat.data, [1] * 10))
            else: # real
                eye_matrix = eye(10)
                self.assertTrue(array_equal(mat, eye_matrix))
            A.write_dmi()

    def test_eye5cd(self):
        """tests the EYE5CD matrices"""
        fnames = [
            'mat_t_dn.op4',
            'mat_t_s1.op4',
            'mat_t_s2.op4',
        ]
        for fname in fnames:
            op4_filename = os.path.join(OP4_PATH, fname)
            matrices = read_op4(op4_filename, debug=False)
            A = matrices['EYE5CD']
            form = A.form
            mat = A.data
            self.assertEqual(form, 6)  # form=6 -> Symmetric
            if 's' in fname:  # sparse
                self.assertTrue(array_equal(mat.row, range(5)))
                self.assertTrue(array_equal(mat.col, range(5)))
                self.assertTrue(array_equal(mat.data, [-1+1j] * 5))
            else: # real
                eye_matrix = -eye(5) + 1j * eye(5)
                self.assertTrue(array_equal(mat, eye_matrix))
            A.write_dmi()

    def test_null(self):
        """tests the NULL matrices"""
        fnames = [
            'mat_t_dn.op4',
            'mat_t_s1.op4',
            'mat_t_s2.op4',
        ]
        for fname in fnames:
            #print('-------%s-------' % fname)
            op4_filename = os.path.join(OP4_PATH, fname)
            matrices = read_op4(op4_filename, debug=False)
            amat = matrices['NULL']
            amat.write_dmi()
            form = amat.form
            mat = amat.data
            self.assertEqual(form, 6)  # form=6 -> Symmetric
            #print A.shape

            # kind of strange that the NULL matrix is dense...
            #if 's' in fname:  # sparse
            if 0:
                self.assertTrue(array_equal(mat.row, range(3)))
                self.assertTrue(array_equal(mat.col, range(3)))
                msg = 'fname=%s NULL sparse matrix error' % fname
                self.assertTrue(array_equal(mat.data, [0] * 3), msg)
            else: # real
                zero_matrix = zeros((3, 3))
                msg = 'fname=%s NULL dense matrix error' % fname
                self.assertTrue(array_equal(mat, zero_matrix))
            del amat

    def test_bad_inputs_1(self):
        op4 = OP4(debug=False)
        op4_filename = os.path.join(OP4_PATH, 'bad_inputs.op4')
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
            op4.read_op4(op4_filename, precision='default_bad')
        with self.assertRaises(IOError):
            op4.read_op4('op4_filename', precision='default')

        # now the inputs are valid, so this works
        unused_matrices2 = op4.read_op4(op4_filename, precision='default')

    def test_file_obj_ascii(self):
        """tests ascii writing"""
        op4 = OP4(debug=False)
        form1 = 1
        amat = ones((3, 3), dtype='float64')
        matrices = {
            'A1': (form1, amat),
        }
        op4_filename = os.path.join(OP4_PATH, 'file_ascii.op4')
        with open(op4_filename, 'w') as op4_file:
            op4.write_op4(op4_file, matrices, name_order='A1', precision='default',
                          is_binary=False)
        os.remove(op4_filename)

        matrices = {
            'A1': Matrix('A1', form1, data=amat),
        }
        op4_filename = os.path.join(OP4_PATH, 'file_ascii.op4')
        with open(op4_filename, 'w') as op4_file:
            op4.write_op4(op4_file, matrices, name_order='A1', precision='default',
                          is_binary=False)
        os.remove(op4_filename)

    @staticmethod
    def test_file_obj_binary():
        op4 = OP4(debug=False)
        form1 = 1
        A1 = ones((3, 3), dtype='float64')
        matrices = {
            'A1': (form1, A1),
        }
        op4_filename = os.path.join(OP4_PATH, 'file_binary.op4')
        with open(op4_filename, 'wb') as op4_file:
            op4.write_op4(op4_file, matrices, name_order='A1', precision='default',
                          is_binary=True)
        os.remove(op4_filename)

    def test_square_matrices_1(self):
        """tests reading/writing square matrices (A1, A2, A3)"""
        op4 = OP4(debug=False)
        #matrices = op4.read_op4(os.path.join(OP4_PATH, fname))
        form1 = 1
        form2 = 2
        form3 = 2
        A1 = np.ones((3, 3), dtype='float64')
        A2 = reshape(arange(9, dtype='float64'), (3, 3))
        A3 = np.ones((1, 1), dtype='float32')
        matrices = {
            'A1': Matrix('A1', form1, data=A1),
            'A2': Matrix('A2', form2, data=A2),
            'A3': Matrix('A3', form3, data=A3),
        }

        for (unused_is_binary, fname) in [(False, 'small_ascii.op4'), (True, 'small_binary.op4')]:
            op4_filename = os.path.join(OP4_PATH, fname)
            op4.write_op4(op4_filename, matrices, name_order=None, precision='default',
                          is_binary=False)
            matrices2 = op4.read_op4(op4_filename, precision='default')
            A1b = matrices2['A1']
            A2b = matrices2['A2']
            self.assertEqual(form1, A1b.form)
            self.assertEqual(form2, A2b.form)

            A1b = matrices2['A1']
            A2b = matrices2['A2']
            A3b = matrices2['A3']
            self.assertEqual(form1, A1b.form)
            self.assertEqual(form2, A2b.form)
            self.assertEqual(form3, A3b.form)

            self.assertTrue(array_equal(A1.data, A1b.data))
            self.assertTrue(array_equal(A2.data, A2b.data))
            self.assertTrue(array_equal(A3.data, A3b.data))
            del A1b, A2b, A3b

    #def test_compress_column(self):
        #compress_column([14, 15, 16, 20, 21, 22, 26, 27, 28])

    def test_qhh_reading(self):
        """tests QHH reading"""
        op4_filename = os.path.join(PKG_PATH, '..', 'models', 'aero',
                                    'bah_plane', 'bah_plane_qhh.op4')
        assert os.path.exists(op4_filename), op4_filename
        matrices = read_op4(op4_filename=op4_filename)
        qhh = matrices['QHH']
        forms = qhh.form
        assert len(forms) == 30, forms
        #print('forms=', forms)
        #print('qhh=', np.dstack(qhh).shape)

    def test_main(self):
        """tests various matrices"""
        #from pyNastran.op4.utils import write_dmig

        op4_filenames = [
            # name, write_binary
            ('mat_t_dn.op4', False),
            ('mat_t_s1.op4', False),
            ('mat_t_s2.op4', False),
            ('mat_b_dn.op4', False),
            ('mat_b_s1.op4', False),
            ('mat_b_s2.op4', False),
            #'b_sample.op4',
            #'binary.op4',
        ]

        #matrix_names = 'EYE10' # identity
        #matrix_names = 'LOW'
        #matrix_names = 'RND1RS' # real,single
        #matrix_names = 'RND1RD' # real,double
        #matrix_names = 'RND1CS' # complex,single
        #matrix_names = 'RND1CD' # complex,double
        #matrix_names = 'STRINGS'
        #matrix_names = 'EYE5CD' # complex identity
        matrix_names = None
        strings = get_matrices()

        unused_is_big_mat = True
        with open('ascii.op4', 'w') as op4_filea, open('binary.op4', 'wb') as op4_fileb:

            op4 = OP4(debug=False)
            matrices = {'strings' : (2, strings)}
            name = 'strings'
            op4.write_op4(op4_filea, matrices, name_order=name,
                          is_binary=False)
            #op4.write_op4(op4_fileb, matrices, name_order=name,
                          #is_binary=True)

            for op4_filename, write_binary in op4_filenames:
                op4 = OP4()
                op4.endian = '>'
                #if 't' in fname:
                #else:
                    #f = open('binary.op4', 'wb')
                op4_filename = os.path.join(OP4_PATH, op4_filename)
                matrices = op4.read_op4(op4_filename, matrix_names=matrix_names,
                                        precision='default')
                #print("keys = %s" % matrices.keys())
                #print("fname=%s" % fname)
                for name, unused_matrix in sorted(matrices.items()):
                    op4.write_op4(op4_filea, matrices, name_order=name,
                                  is_binary=False)
                    if write_binary:
                        op4.write_op4(op4_fileb, matrices, name_order=name,
                                      is_binary=True)

        os.remove('ascii.op4')
        os.remove('binary.op4')

    def test_op4_plate(self):
        """tests sparse binary example"""
        op4_filename = os.path.join(OP4_PATH, 'testplate_kgg.op4')
        matrices = read_op4(op4_filename=op4_filename,
                            matrix_names=None, precision='default', debug=False, log=None)
        Kgg = matrices['KGG'].data.toarray()

        op4_filename = os.path.join(OP4_PATH, 'testplate_kgg_ascii.op4')
        matrices = read_op4(op4_filename=op4_filename,
                            matrix_names=None, precision='default', debug=False, log=None)
        Kgga = matrices['KGG'].data.toarray()
        assert np.allclose(Kgg, Kgga)
        #for line in Kgg:
            #print(line)

def get_matrices():
    """creates dummy matrices"""
    strings = np.array([
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
        [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
        [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
        [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
        [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype='float32') # f?
    return strings


if __name__ == '__main__':  # pragma: no cover
    ON_RTD = os.environ.get('READTHEDOCS', None) == 'True'
    if not ON_RTD:
        unittest.main()
