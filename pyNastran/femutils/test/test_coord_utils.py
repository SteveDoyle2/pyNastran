"""tests coordinate system related femutils"""
# -*- coding: utf-8 -*-
# pylint: disable=R0201, C0103
import unittest

import numpy as np

from pyNastran.femutils.coord_transforms import (
    xyz_to_rtz_array, xyz_to_rtp_array,
    rtz_to_xyz_array, rtp_to_xyz_array,
    rtz_to_rtp_array, rtp_to_rtz_array,
    cylindrical_rotation_matrix,
    spherical_rotation_matrix,
)
from pyNastran.femutils.coord_utils import (
    coords_from_vector_1d,
    coordinate_system_from_vector_2d_tri,
    coordinate_system_from_vector_2d_quad,
    shape4, shape4_to_xyz,
)
#from .utils import is_array_close

__all__ = ['TestShapeFunction', 'TestCoordUtils']

class TestShapeFunction(unittest.TestCase):
    """various shape function tests"""
    def test_shape4_to_xyz(self):
        """
        tests:
         - shape4
         - shape4_to_xyz

        """
        shape = [0.0, 0.0]
        p1 = [0., 0., 0.]
        p2 = [0., 1., 0.]
        p3 = [1., 1., 0.]
        p4 = [1., 0., 0.]
        nquads = 1
        p1234 = [[p1, p2, p3, p4]]
        #nquads, npoints, nxyz = p1234.shape
        #assert npoints == 4 and nxyz == 3, 'shape=(nquads, npoints, nxyz) = %s' % str(p1234.shape)
        n4 = shape4(shape)
        assert n4.shape == (nquads, 4), 'shape=%s' % str(n4.shape)

        pcentroid = shape4_to_xyz(p1234, n4)
        assert np.array_equal(pcentroid, [[0.5, 0.5, 0.]])
        #p1234v = np.vstack([p1234, p1234])
        #n4v = np.vstack([n4, n4])
        #print('n4v.shape =', n4v.shape)
        #pcentroid = shape4_to_xyz(p1234v, n4v)
        #print(pcentroid)

    #def test_shape8_to_xyz(self):
        #shape = [0.5, 0.5, 0.5]
        #n = shape8(shape)

class TestCoordUtils(unittest.TestCase):
    """tests coordinate system related femutils"""
    def test_xyz_to_rtz_array(self):
        """
        tests:
         - xyz_to_rtz_array
         - rtz_to_xyz_array
        """
        xyz1 = [
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 2., 0.],
            [1., 2., 3.],
        ]
        rtz = xyz_to_rtz_array(xyz1)
        xyz2 = rtz_to_xyz_array(rtz)
        assert np.allclose(xyz1, xyz2), 'xyz1:\n%s\nxyz2:\n%s' % (xyz1, xyz2)

    def test_xyz_to_rtp_array(self):
        """
        tests:
         - xyz_to_rtp_array
         - rtp_to_xyz_array
        """
        xyz1 = [
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 2., 0.],
            [1., 2., 3.],
        ]
        rtp = xyz_to_rtp_array(xyz1)
        xyz2 = rtp_to_xyz_array(rtp)
        assert np.allclose(xyz1, xyz2), 'xyz1:\n%s\nxyz2:\n%s' % (xyz1, xyz2)

    def test_rtz_to_rtp_array(self):
        """
        tests:
         - rtz_to_rtp_array
         - rtp_to_rtz_array
        """
        rtz1 = [
            [0., 0., 0.],

            [1., 0., 0.],
            [1., 45., 0.],
            [1., 90., 0.],

            [1., 0., 2.],
            [1., 45., 2.],
            [1., 90., 2.],
        ]
        rtp = rtz_to_rtp_array(rtz1)
        rtz2 = rtp_to_rtz_array(rtp)
        assert np.allclose(rtz1, rtz2), 'rtz1:\n%s\nrtz2:\n%s' % (rtz1, rtz2)

    #---------------------------------------------------------------------------

    def test_cylindrical_rotation_matrix(self):
        """tests cylindrical_rotation_matrix

        Verifies:
        - Orthogonality (R @ R.T = I)
        - det = 1 (proper rotation)
        - Known angle values:
          theta=0: [[1,0,0],[0,1,0],[0,0,1]] (identity)
          theta=90: [[0,-1,0],[1,0,0],[0,0,1]]
        """
        theta = np.radians([0., 45., 90., 180., 270.])
        R = cylindrical_rotation_matrix(theta, dtype='float64')
        assert R.shape == (5, 3, 3)

        for i in range(5):
            assert np.allclose(R[i] @ R[i].T, np.eye(3), atol=1e-15)
            assert np.isclose(np.linalg.det(R[i]), 1.0)

        # theta=0: identity
        assert np.allclose(R[0], np.eye(3), atol=1e-15)
        # theta=90: [[0,-1,0],[1,0,0],[0,0,1]]
        expected_90 = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        assert np.allclose(R[2], expected_90, atol=1e-15)
        # theta=180: [[-1,0,0],[0,-1,0],[0,0,1]]
        expected_180 = np.array([[-1., 0., 0.], [0., -1., 0.], [0., 0., 1.]])
        assert np.allclose(R[3], expected_180, atol=1e-15)

    def test_spherical_rotation_matrix(self):
        """tests spherical_rotation_matrix

        Verifies:
        - Orthogonality (R @ R.T = I) at multiple angles
        - det = 1 (proper rotation, not reflection)
        - Known angle values:
          (theta=0, phi=0): e_r=[0,0,1], e_theta=[1,0,0], e_phi=[0,1,0]
          (theta=90, phi=0): e_r=[1,0,0], e_theta=[0,0,-1], e_phi=[0,1,0]
          (theta=90, phi=90): e_r=[0,1,0], e_theta=[0,0,-1], e_phi=[-1,0,0]
        """
        thetar = np.radians([0., 90., 90., 45., 135., 60.])
        phir = np.radians([0., 0., 90., 45., 270., 30.])
        R = spherical_rotation_matrix(thetar, phir)
        assert R.shape == (6, 3, 3)

        for i in range(6):
            assert np.allclose(R[i] @ R[i].T, np.eye(3), atol=1e-14), (
                f'not orthogonal at i={i}')
            assert np.isclose(np.linalg.det(R[i]), 1.0), (
                f'det != 1 at i={i}: {np.linalg.det(R[i])}')

        # (theta=0, phi=0): pole
        assert np.allclose(R[0, :, 0], [0., 0., 1.], atol=1e-15)  # e_r
        assert np.allclose(R[0, :, 1], [1., 0., 0.], atol=1e-15)  # e_theta
        assert np.allclose(R[0, :, 2], [0., 1., 0.], atol=1e-15)  # e_phi

        # (theta=90, phi=0): +x axis
        assert np.allclose(R[1, :, 0], [1., 0., 0.], atol=1e-15)
        assert np.allclose(R[1, :, 1], [0., 0., -1.], atol=1e-15)
        assert np.allclose(R[1, :, 2], [0., 1., 0.], atol=1e-15)

        # (theta=90, phi=90): +y axis
        assert np.allclose(R[2, :, 0], [0., 1., 0.], atol=1e-15)
        assert np.allclose(R[2, :, 1], [0., 0., -1.], atol=1e-15)
        assert np.allclose(R[2, :, 2], [-1., 0., 0.], atol=1e-15)

    def test_coords_from_vector_1d(self):
        """tests coords_from_vector_1d"""
        v = [ # duplicate
            [0, 0., 1.],
            [0., 0., 1.],
        ]
        expected = np.array([
            [[0., 0., 1.],
             [0., 1., 0.],
             [-1., 0., 0.],],

            [[0., 0., 1.],
             [0., 1., 0.],
             [-1., 0., 0.],],
        ])
        #out = perpendicular_vector2d(v)
        out = coords_from_vector_1d(v)
        assert np.allclose(out, expected)

    def test_coordinate_system_from_vector_2d_tri(self):
        """tests coordinate_system_from_vector_2d_tri"""
        xyz1 = [0., 0., 0.]
        xyz2 = [1., 0., 0.]
        xyz3 = [0., 1., 0.]
        coords1 = coordinate_system_from_vector_2d_tri(xyz1, xyz2, xyz3)
        assert np.allclose(coords1[0, 2, 2], 1.0), '\n' + str(coords1[0, :, :])

        xyz1 = [
            [0., 0., 0.],
            [0., 0., 0.],
        ]
        xyz2 = [
            [1., 0., 0.],
            [1., 0., 0.],
        ]
        xyz3 = [
            [0., 1., 0.],
            [0., 1., 0.],
        ]
        coords2 = coordinate_system_from_vector_2d_tri(xyz1, xyz2, xyz3)
        two_coords = np.vstack([coords1, coords1])
        assert np.array_equal(two_coords, coords2)

    def test_coordinate_system_from_vector_2d_quad(self):
        """tests coordinate_system_from_vector_2d_quad"""
        xyz1 = [0., 0., 0.]
        xyz2 = [1., 0., 0.]
        xyz3 = [1., 1., 0.]
        xyz4 = [0., 1., 0.]
        coords1 = coordinate_system_from_vector_2d_quad(xyz1, xyz2, xyz3, xyz4)
        assert np.allclose(coords1[0, 2, 2], 1.0), '\n' + str(coords1[0, :, :])

        xyz1 = [
            [0., 0., 0.],
            [0., 0., 0.],
        ]
        xyz2 = [
            [1., 0., 0.],
            [1., 0., 0.],
        ]
        xyz3 = [
            [1., 1., 0.],
            [1., 1., 0.],
        ]
        xyz4 = [
            [0., 1., 0.],
            [0., 1., 0.],
        ]
        coords2 = coordinate_system_from_vector_2d_quad(xyz1, xyz2, xyz3, xyz4)
        two_coords = np.vstack([coords1, coords1])
        assert np.array_equal(two_coords, coords2)
