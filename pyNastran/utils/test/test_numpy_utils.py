# -*- coding: utf-8 -*-

import os
import unittest

from six import PY2, StringIO
import numpy as np

import pyNastran
from pyNastran.utils.numpy_utils import (
    loadtxt_nice, augmented_identity, savetxt_nice,
    perpendicular_vector, perpendicular_vector2d,
    dot3d,
)
from pyNastran.utils.numpy_functions.coord_utils import (
    cylindrical_rotation_matrix,
)

PKG_PATH = pyNastran.__path__[0]


def is_array_close(v1, v2):
    """are two arrays close"""
    return np.all(np.isclose(v1, v2))

class TestNumpyUtils(unittest.TestCase):
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

    def coords_from_vector_1d(self):
        """tests coords_from_vector_1d"""
        v = [ # duplicate
            [0, 0., 1.],
            [0., 0., 1.],
        ]
        expected = np.array([
            [[0., 0., 1.],
             [ 0., 1., 0.],
             [-1., 0., 0.],],

            [[ 0., 0., 1.],
             [ 0., 1., 0.],
             [-1., 0., 0.],],
        ])
        #out = perpendicular_vector2d(v)
        out = coords_from_vector_1d(v)
        assert np.allclose(out, expected)

    def test_cylindrical_rotation_matrix(self):
        """tests cylindrical_rotation_matrix"""
        theta = [0., 0.]
        coords = cylindrical_rotation_matrix(theta, dtype='float64')

        theta = np.radians([0., 45., 90.])
        coords = cylindrical_rotation_matrix(theta, dtype='float64')
        #print(coords)

    def test_dot3d(self):
        """tests dot3d"""
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
        C = dot3d(A, B)
        print('-------')
        for Ci in C:
            print(Ci.shape)
            print(Ci)
            print('-------')

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

    def test_loadtxt_01(self):
        """tests that we can reimplement loadtxt so it doesn't suck"""
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
        """tests that we can reimplement savetxt so it doesn't suck"""
        A = np.eye(10)
        csv_filename = 'savetxt_real.csv'
        savetxt_nice(csv_filename, A, fmt='%.18e', delimiter=',', newline='\n',
                     header='', footer='', comments='# ')

        with self.assertRaises(ValueError):
            loadtxt_nice(csv_filename, delimiter=' ', skiprows=0, comment='#',
                         dtype=np.float64, converters=None,
                         usecols=None, unpack=False, ndmin=0)

        A2 = loadtxt_nice(csv_filename, delimiter=',', skiprows=0, comment='#',
                          dtype=np.float64, converters=None,
                          usecols=None, unpack=False, ndmin=0)
        assert np.array_equal(A, A2), 'expected:\n%s\nactual:\n%s' % (A, A2)
        os.remove(csv_filename)

        csv_filename = 'savetxt_complex.csv'
        B = np.eye(10, dtype='complex128') - 2 * A*1j
        savetxt_nice(csv_filename, B, fmt='%.18e', delimiter=',', newline='\n',
                     header='', footer='', comments='# ')
        with self.assertRaises(ValueError):  ## TODO: mistake
            B2 = loadtxt_nice(csv_filename, delimiter=',', skiprows=0, comment='#',
                              dtype=np.float64, converters=None,
                              usecols=None, unpack=False, ndmin=0)
            #assert np.array_equal(B, B2), 'expected:\n%s\nactual:\n%s' % (B, B2)
        os.remove(csv_filename)

        if 0:  ## TODO: not done with filehandle test
            from codecs import open
            with open(csv_filename, 'w') as csv_file:
                savetxt_nice(csv_file, B, fmt='%.18e', delimiter=',', newline='\n',
                             header='', footer='', comments='# ')
            os.remove(csv_filename)

        if PY2:
            with self.assertRaises(IOError):
                B2 = loadtxt_nice('missing.txt', delimiter=',', skiprows=0, comment='#',
                                  dtype=np.float64, converters=None,
                                  usecols=None, unpack=False, ndmin=0)
        else:
            with self.assertRaises(FileNotFoundError):
                B2 = loadtxt_nice('missing.txt', delimiter=',', skiprows=0, comment='#',
                                  dtype=np.float64, converters=None,
                                  usecols=None, unpack=False, ndmin=0)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
