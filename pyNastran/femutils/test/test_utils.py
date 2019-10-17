"""tests general femutils"""
# -*- coding: utf-8 -*-
# pylint:  disable=R0201,C0103
__all__ = ['TestMatrix3d', 'TestNumpyUtils', 'TestFemIO']

import os
import unittest
from io import StringIO, BytesIO

import numpy as np
from numpy.testing import assert_equal, assert_array_equal

import pyNastran
from pyNastran.femutils.io import loadtxt_nice, savetxt_nice
from pyNastran.femutils.matrix3d import dot_n33_n33, transpose3d, triple_n33_n33, triple_n33_33
from pyNastran.femutils.utils import augmented_identity, perpendicular_vector, perpendicular_vector2d
from pyNastran.femutils.coord_transforms import cylindrical_rotation_matrix

from pyNastran.femutils.test.utils import is_array_close

PKG_PATH = pyNastran.__path__[0]


#class TestNan(unittest.TestCase):
class TestMatrix3d(unittest.TestCase):
    """tests functions in femutils.matrix3d"""
    def test_dot_n33_n33(self):
        """tests dot_n33_n33"""
        A = np.array([
            [[1., 0., 0.],
             [0., 1., 0.],
             [0., 0., 1.],],

            [[1., 0., 0.],
             [0., 1., 0.],
             [0., 0., 1.],],

            [[1., 0., 0.],
             [0., 1., 0.],
             [0., 0., 1.],],
        ])
        theta = np.radians([0., 45., 90])
        B = cylindrical_rotation_matrix(theta, dtype='float64')
        C = dot_n33_n33(A, B)
        #print('-------')
        #for Ci in C:
            #print(Ci.shape)
            #print(Ci)
            #print('-------')
        ## TODO: not compared

    def test_trapose3d(self):
        """tests transpose3d"""
        A = np.array([
            [[1., 2., 3.],
             [4., 5., 6.],
             [7., 8., 9.],],

            [[1., 2., 3.],
             [4., 5., 6.],
             [7., 8., 9.],],

            [[1., 2., 3.],
             [4., 5., 6.],
             [7., 8., 9.],],
        ])

        B_expected = np.array([
            [[1., 4., 7.],
             [2., 5., 8.],
             [3., 6., 9.],],

            [[1., 4., 7.],
             [2., 5., 8.],
             [3., 6., 9.],],

            [[1., 4., 7.],
             [2., 5., 8.],
             [3., 6., 9.],],
        ])
        B_actual = transpose3d(A)
        #print('transpose B_actual:')
        #print(B_actual)
        assert is_array_close(B_expected, B_actual)

    def test_triple_n33_n33(self):
        """tests triple_n33_n33"""
        A = np.array([
            [[1., 2., 3.],
             [4., 5., 6.],
             [7., 8., 9.],],

            [[1., 2., 3.],
             [4., 5., 6.],
             [7., 8., 9.],],
        ])
        T = np.array([
            [[0., 1., 0.],
             [0., 0., 1.],
             [1., 0., 9.],],

            [[0., 1., 0.],
             [0., 0., 1.],
             [1., 0., 9.],],
        ])
        TtAT_actual = triple_n33_n33(A, T, tranpose=False)
        TtAT_expected = [
            [[5., 6., 58.],
             [8., 9., 88.],
             [74., 84., 820.],],

            [[5., 6., 58.],
             [8., 9., 88.],
             [74., 84., 820.],],
        ]
        assert is_array_close(TtAT_expected, TtAT_actual)
        TATt_actual = triple_n33_n33(A, T, tranpose=True)
        TATt_expected = [
            [[9., 7., 89.],
             [3., 1., 29.],
             [87., 67., 860.],],

            [[9., 7., 89.],
             [3., 1., 29.],
             [87., 67., 860.],],
        ]
        assert is_array_close(TATt_expected, TATt_actual)

    def test_triple_n33_33(self):
        """tests test_triple_n33_33"""
        A = np.array([
            [[1., 2., 3.],
             [4., 5., 6.],
             [7., 8., 9.],],

            [[1., 2., 3.],
             [4., 5., 6.],
             [7., 8., 9.],],
        ])
        T = np.array([
            [0., 1., 0.],
            [0., 0., 1.],
            [1., 0., 9.],
        ])
        TtAT_actual = triple_n33_33(A, T, tranpose=False)
        TtAT_expected = [
            [[5., 6., 58.],
             [8., 9., 88.],
             [74., 84., 820.],],

            [[5., 6., 58.],
             [8., 9., 88.],
             [74., 84., 820.],],
        ]
        assert is_array_close(TtAT_expected, TtAT_actual)
        TATt_actual = triple_n33_33(A, T, tranpose=True)
        TATt_expected = [
            [[9., 7., 89.],
             [3., 1., 29.],
             [87., 67., 860.],],

            [[9., 7., 89.],
             [3., 1., 29.],
             [87., 67., 860.],],
        ]
        assert is_array_close(TATt_expected, TATt_actual)


class TestNumpyUtils(unittest.TestCase):
    """tests functions in femutils.utils"""
    def test_perpendicular_vector(self):
        """tests perpendicular_vector"""
        with self.assertRaises(ValueError):
            perpendicular_vector([0., 0., 0.])

        a1 = perpendicular_vector([1., 0., 0.])
        self.assertTrue(is_array_close(a1, [0., 1., 0.]), msg=str(a1))

        a2 = perpendicular_vector([1., 1., 0.])
        self.assertTrue(is_array_close(a2, [0., 0., 1.]), msg=str(a2))

        a3 = perpendicular_vector([1., 1., 1.])
        self.assertTrue(is_array_close(a3, [1., 1., -2.]), msg=str(a3))

        a1 = perpendicular_vector2d([1., 0., 0.])
        a1 = perpendicular_vector2d((1., 0., 0.))

        #-----------------------------
        expected = np.array([
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [1., 1., 0.],
            [1., 1., 1.],
        ])
        out = perpendicular_vector2d(expected)
        expected2 = np.array([     # input
            [0., 0., 1.],   # [1., 0., 0.],
            [0., 0., 1.],   # [0., 1., 0.],
            [0., 1., 0.],   # [0., 0., 1.],
            [0., 0., 1.],   # [1., 1., 0.],
            [1., 1., -2.],  # [1., 1., 1.],
        ])
        assert np.allclose(out, expected2)
        #print('out')
        #print(out)
        #print('-----------')
        #print('expected')
        #print(v2)
        #print('-----------')
        #print('diff')
        #print(v2 - out)

    def test_augmented_identity(self):
        """tests augmented_identity"""
        expected_array = np.array([
            [1., 0., 0., 0.],
            [0., 1., 0., 0.],
            [0., 0., 1., 0.],
        ])
        actual_array = augmented_identity(3, 4)
        msg = 'expected:\n%s\nactual:\n%s' % (expected_array, actual_array)
        assert np.array_equal(expected_array, actual_array), msg

class TestFemIO(unittest.TestCase):
    """tests functions in femutils.io"""

    def test_file_roundtrip(self):
        """per numpy"""
        a = np.array([(1, 2), (3, 4)])
        np.savetxt('temp.txt', a)
        b = loadtxt_nice('temp.txt')
        assert_array_equal(a, b)
        os.remove('temp.txt')

    def test_record(self):
        """per numpy"""
        c = StringIO()
        c.write('1 2\n3 4')
        c.seek(0)

        x = loadtxt_nice(c, dtype=[('x', np.int32), ('y', np.float32)])
        unused_x2 = np.loadtxt(c, dtype=[('x', np.int32), ('y', np.float32)])

        #print('x =', x, type(x2))
        #print('x2 =', x2, type(x2))
        #a = np.array([(1, 2), (3, 4)], dtype=[('x', 'i4'), ('y', 'i4')])
        #assert_array_equal(x, a)

        d = StringIO()
        d.write('M 64 75.0\nF 25 60.0')
        d.seek(0)
        mydescriptor = {'names': ('gender', 'age', 'weight'),
                        'formats': ('S1', 'i4', 'f4')}
        b = np.array([('M', 64.0, 75.0),
                      ('F', 25.0, 60.0)], dtype=mydescriptor)
        y = loadtxt_nice(d, dtype=mydescriptor)
        #assert_array_equal(y, b)

    #def test_loadtxt_fields_subarrays(self):
        #"""per numpy"""
        #from numpy.testing import assert_, assert_equal, temppath, assert_array_equal
        ## For ticket #1936
        #dt = [("a", 'u1', 2), ("b", 'u1', 2)]
        #x = np.loadtxt(StringIO("0 1 2 3"), dtype=dt)
        #assert_equal(x, np.array([((0, 1), (2, 3))], dtype=dt))

        #dt = [("a", [("a", 'u1', (1, 3)), ("b", 'u1')])]
        #x = np.loadtxt(StringIO("0 1 2 3"), dtype=dt)
        #assert_equal(x, np.array([(((0, 1, 2), 3),)], dtype=dt))

        #dt = [("a", 'u1', (2, 2))]
        #x = np.loadtxt(StringIO("0 1 2 3"), dtype=dt)
        #assert_equal(x, np.array([(((0, 1), (2, 3)),)], dtype=dt))

        #dt = [("a", 'u1', (2, 3, 2))]
        #x = np.loadtxt(StringIO("0 1 2 3 4 5 6 7 8 9 10 11"), dtype=dt)
        #data = [((((0, 1), (2, 3), (4, 5)), ((6, 7), (8, 9), (10, 11))),)]
        #assert_equal(x, np.array(data, dtype=dt))


    #def test_complex_arrays(self):
        #"""per numpy"""
        #from io import BytesIO, StringIO
        #from numpy.testing import assert_, assert_equal, temppath, assert_array_equal
        #ncols = 2
        #nrows = 2
        #a = np.zeros((ncols, nrows), dtype=np.complex128)
        #re = np.pi
        #im = np.e
        #a[:] = re + 1.0j * im

        ## One format only
        #c = BytesIO()
        #np.savetxt(c, a, fmt=' %+.3e')
        #c.seek(0)
        #lines = c.readlines()
        #assert_equal(
            #lines,
            #[b' ( +3.142e+00+ +2.718e+00j)  ( +3.142e+00+ +2.718e+00j)\n',
             #b' ( +3.142e+00+ +2.718e+00j)  ( +3.142e+00+ +2.718e+00j)\n'])

        ## One format for each real and imaginary part
        #c = BytesIO()
        #np.savetxt(c, a, fmt='  %+.3e' * 2 * ncols)
        #c.seek(0)
        #lines = c.readlines()
        #assert_equal(
            #lines,
            #[b'  +3.142e+00  +2.718e+00  +3.142e+00  +2.718e+00\n',
             #b'  +3.142e+00  +2.718e+00  +3.142e+00  +2.718e+00\n'])

        ## One format for each complex number
        #c = BytesIO()
        #np.savetxt(c, a, fmt=['(%.3e%+.3ej)'] * ncols)
        #c.seek(0)
        #lines = c.readlines()
        #assert_equal(
            #lines,
            #[b'(3.142e+00+2.718e+00j) (3.142e+00+2.718e+00j)\n',
             #b'(3.142e+00+2.718e+00j) (3.142e+00+2.718e+00j)\n'])

    def test_complex_negative_exponent(self):
        """per numpy"""
        # Previous to 1.15, some formats generated x+-yj, gh 7895
        ncols = 2
        nrows = 2
        a = np.zeros((ncols, nrows), dtype=np.complex128)
        re = np.pi
        im = np.e
        a[:] = re - 1.0j * im
        c = BytesIO()
        savetxt_nice(c, a, fmt='%.3e')
        c.seek(0)
        lines = c.readlines()
        assert_equal(
            lines,
            [b' (3.142e+00-2.718e+00j)  (3.142e+00-2.718e+00j)\n',
             b' (3.142e+00-2.718e+00j)  (3.142e+00-2.718e+00j)\n'])

    def test_loadtxt_nice(self):
        """tests that we can reimplement loadtxt so it has good error messages"""
        str_data = StringIO("1,0,2\n3,0,4")
        x1, y1 = np.loadtxt(str_data, delimiter=',', usecols=(0, 2), unpack=True)
        x2, y2 = loadtxt_nice(str_data, delimiter=',', usecols=(0, 2), unpack=True)
        #print('x1=%s y1=%s' % (x1, y1))
        #print('x2=%s y2=%s\n+' % (x2, y2))
        assert np.array_equal(x1, x2), 'x1=%s x2=%s' % (x1, x2)
        assert np.array_equal(y1, y2), 'y1=%s y2=%s' % (y1, y2)

        str_data = StringIO("#1,0,2\n3,0,4")
        x1, y1 = np.loadtxt(str_data, delimiter=',', usecols=(0, 2), unpack=True)
        x2, y2 = loadtxt_nice(str_data, delimiter=',', usecols=(0, 2), unpack=True)
        #print('x1=%s y1=%s' % (x1, y1))
        #print('x2=%s y2=%s' % (x2, y2))
        assert np.array_equal(x1, x2), 'x1=%s x2=%s' % (x1, x2)
        assert np.array_equal(y1, y2), 'y1=%s y2=%s' % (y1, y2)

        str_data = StringIO("#1,0,2\n3,0,4")
        x1, y1 = np.loadtxt(str_data, delimiter=',', usecols=(0, 2), unpack=True, ndmin=1)
        x2, y2 = loadtxt_nice(str_data, delimiter=',', usecols=(0, 2), unpack=True, ndmin=1)
        #print('x1=%s y1=%s' % (x1, y1))
        #print('x2=%s y2=%s' % (x2, y2))
        assert np.array_equal(x1, x2), 'x1=%s x2=%s' % (x1, x2)
        assert np.array_equal(y1, y2), 'y1=%s y2=%s' % (y1, y2)

        str_data = StringIO("#1,0,2\n3,0,4")
        x1, y1 = np.loadtxt(str_data, delimiter=',', usecols=(0, 2), unpack=True, ndmin=2)
        x2, y2 = loadtxt_nice(str_data, delimiter=',', usecols=(0, 2), unpack=True, ndmin=2)
        #print('x1=%s y1=%s' % (x1, y1))
        #print('x2=%s y2=%s' % (x2, y2))
        assert np.array_equal(x1, x2), 'x1=%s x2=%s' % (x1, x2)
        assert np.array_equal(y1, y2), 'y1=%s y2=%s' % (y1, y2)

    def test_savetxt_nice(self):
        """tests that we can reimplement savetxt so it works on unicode for unicode file handlers"""
        A = np.eye(10)
        csv_filename = 'savetxt_real.csv'
        savetxt_nice(csv_filename, A, fmt='%.18e', delimiter=',', newline='\n',
                     header='', footer='', comments='# ')

        with self.assertRaises(ValueError):
            loadtxt_nice(csv_filename, delimiter=' ', skiprows=0, comments='#',
                         dtype=np.float64, converters=None,
                         usecols=None, unpack=False, ndmin=0)

        A2 = loadtxt_nice(csv_filename, delimiter=',', skiprows=0, comments='#',
                          dtype=np.float64, converters=None,
                          usecols=None, unpack=False, ndmin=0)
        assert np.array_equal(A, A2), 'expected:\n%s\nactual:\n%s' % (A, A2)
        os.remove(csv_filename)

        csv_filename = 'savetxt_complex.csv'
        B = np.eye(10, dtype='complex128') - 2 * A*1j
        savetxt_nice(csv_filename, B, fmt='%.18e', delimiter=',', newline='\n',
                     header='', footer='', comments='# ')
        with self.assertRaises(ValueError):  ## TODO: mistake
            unused_B2 = loadtxt_nice(csv_filename, delimiter=',', skiprows=0, comments='#',
                                     dtype=np.float64, converters=None,
                                     usecols=None, unpack=False, ndmin=0)
            #assert np.array_equal(B, B2), 'expected:\n%s\nactual:\n%s' % (B, B2)
        os.remove(csv_filename)

        if 0:  ## TODO: not done with filehandle test
            with open(csv_filename, 'w') as csv_file:
                savetxt_nice(csv_file, B, fmt='%.18e', delimiter=',', newline='\n',
                             header='', footer='', comments='# ')
            os.remove(csv_filename)

        with self.assertRaises(FileNotFoundError):
            B2 = loadtxt_nice('missing.txt', delimiter=',', skiprows=0, comments='#',
                              dtype=np.float64, converters=None,
                              usecols=None, unpack=False, ndmin=0)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
