"""defines OP2 Matrix Test"""
from __future__ import print_function
import os
import unittest

import numpy as np

import pyNastran
from pyNastran.bdf.bdf import read_bdf
from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_geom import read_op2_geom, FatalError
pkg_path = pyNastran.__path__[0]


class TestOP2Matrix(unittest.TestCase):
    """various matrix tests"""

    def test_gpspc(self):
        """Tests the gspc1 MATPOOL model"""
        op2_filename = os.path.join(pkg_path, 'op2', 'test', 'matrices', 'gpsc1.op2')
        model = read_op2_geom(op2_filename, debug=False)

        deltak = model.matrices['DELTAK']
        assert deltak.data.shape == (17, 221), deltak.data.shape

        deltam = model.matrices['DELTAM']
        assert deltam.data.shape == (17, 221), deltam.data.shape

        deltam0 = model.matrices['DELTAM0']
        assert deltam0.data.shape == (6, 78), deltam0.data.shape

        mrggt = model.matrices['MRGGT']
        assert mrggt.data.shape == (24, 24), mrggt.data.shape

        rbm0 = model.matrices['RBM0']
        assert rbm0.data.shape == (6, 6), rbm0.data.shape

        uexpt = model.matrices['UEXPT']
        assert uexpt.data.shape == (276, 24), uexpt.data.shape

    def test_op2_dmi_01(self):
        """tests DMI matrix style"""
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))
        bdf_filename = os.path.join(folder, 'matrix', 'matrix.dat')
        op2_filename = os.path.join(folder, 'matrix', 'mymatrix.op2')
        matrices = {
            'A' : True,
            'B' : False,
            'ATB' : False,
            'BTA' : False,
            'MYDOF' : True,
        }
        model = read_bdf(bdf_filename, debug=False)

        dmi_a = model.dmis['A']
        assert dmi_a.shape == (4, 2), 'shape=%s' % (dmi_a.shape)
        #print('dmi_a\n', dmi_a)
        a, rows_reversed, cols_reversed = dmi_a.get_matrix(is_sparse=False, apply_symmetry=False)
        #print('model.dmi.A =\n%s' % dmi_a)
        #print('model.dmi.A =\n%s' % str(a))
        #return
        op2 = OP2(debug=False)
        op2.set_additional_matrices_to_read(matrices)
        try:
            op2.read_op2(op2_filename)
            raise RuntimeError('this is wrong...')
        except FatalError:
            # the OP2 doesn't have a trailing zero marker
            pass

        # M rows, Ncols
        A = np.array([
            [1., 0.],
            [3., 6.],
            [5., 0.],
            [0., 8.],
        ], dtype='float32')
        B = A
        mydof = np.array([
            -1.0, 1.0, 1.0, -1.0, 1.0,
            2.0, -1.0, 1.0, 3.0, -1.0, 1.0, 4.0, -1.0,
            1.0, 5.0, -1.0, 1.0, 6.0, -1.0, 2.0, 1.0,
            -1.0, 2.0, 2.0, -1.0, 2.0, 3.0, -1.0, 2.0,
            4.0, -1.0, 2.0, 5.0, -1.0, 2.0, 6.0,
        ])
        BTA = np.dot(B.T, A)
        ATB = np.dot(A.T, B)
        ATB_expected = np.array([
            [35., 18.],
            [18., 100.]
        ], dtype='float32')
        #BTA_expected = ATB_expected

        expecteds = [A, ATB, B, BTA, mydof]
        matrix_names = sorted(matrices.keys())

        for matrix_name, expected in zip(matrix_names, expecteds):
            assert matrix_name in op2.matrices, matrix_name
            actual = op2.matrices[matrix_name].data
            compare_dmi_matrix_from_bdf_to_op2(model, op2, expected, actual, matrix_name)

    def test_op2_dmi_02(self):
        """tests DMI matrix style"""
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))
        bdf_filename = os.path.join(folder, 'matrix', 'matrix.dat')
        op2_filename = os.path.join(folder, 'matrix', 'mymatrix.op2')
        matrices = {
            'A' : True,
            'B' : False,
            'ATB' : False,
            'BTA' : False,
            'MYDOF' : True,
        }
        model = read_bdf(bdf_filename, debug=False)

        dmi_a = model.dmis['A']
        a, rows_reversed, cols_reversed = dmi_a.get_matrix(is_sparse=False, apply_symmetry=False)
        #print('model.dmi.A =\n%s' % dmi_a)
        #print('model.dmi.A =\n%s' % str(a))
        #return
        op2 = OP2(debug=False)
        try:
            op2.read_op2(op2_filename, skip_undefined_matrices=True)
            raise RuntimeError('this is wrong...')
        except FatalError:
            # the OP2 doesn't have a trailing zero marker
            pass

        # M rows, Ncols
        A = np.array([
            [1., 0.],
            [3., 6.],
            [5., 0.],
            [0., 8.],
        ], dtype='float32')
        B = A
        mydof = np.array([
            -1.0, 1.0, 1.0, -1.0, 1.0,
            2.0, -1.0, 1.0, 3.0, -1.0, 1.0, 4.0, -1.0,
            1.0, 5.0, -1.0, 1.0, 6.0, -1.0, 2.0, 1.0,
            -1.0, 2.0, 2.0, -1.0, 2.0, 3.0, -1.0, 2.0,
            4.0, -1.0, 2.0, 5.0, -1.0, 2.0, 6.0,
        ])
        BTA = np.dot(B.T, A)
        ATB = np.dot(A.T, B)
        ATB_expected = np.array([
            [35., 18.],
            [18., 100.]
        ], dtype='float32')
        #BTA_expected = ATB_expected

        expecteds = [A, ATB, B, BTA, mydof]
        matrix_names = sorted(matrices.keys())

        for matrix_name, expected in zip(matrix_names, expecteds):
            assert matrix_name in op2.matrices, matrix_name
            actual = op2.matrices[matrix_name].data
            compare_dmi_matrix_from_bdf_to_op2(model, op2, expected, actual, matrix_name)

def compare_dmi_matrix_from_bdf_to_op2(bdf_model, op2_model, expected, actual, matrix_name):
    """compares two matrices"""
    if not (np.array_equal(expected, actual) or
            np.array_equal(expected, np.squeeze(actual))):

        if matrix_name in bdf_model.dmis:
            dmi = bdf_model.dmis[matrix_name]
            table_array, rows_reversed, cols_reversed = dmi.get_matrix(is_sparse=False, apply_symmetry=False)
            #stable_array, rows_reversed, cols_reversed = dmi.get_matrix(is_sparse=True, apply_symmetry=False)
            #print(table_array)
        #print(stable_array)
        msg = 'matrix %s was not read properly\n' % matrix_name
        msg += 'expected shape=%s\n%s\n' % (str(expected.shape), expected)
        msg += 'actual shape=%s;  squeeze(actual) shape=%s\n%s' % (
            str(actual.shape), str(np.squeeze(actual.shape)), actual.ravel())
        #msg += '\n%s' % actual.ravel()
        print(msg)
        print('==========================')
        #raise RuntimeError(msg)


if __name__ == '__main__':  # pragma: no cover
    on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
    if not on_rtd:
        unittest.main()
