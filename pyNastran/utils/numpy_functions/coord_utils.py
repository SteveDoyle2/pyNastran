"""
Defines general coordinate system related functions including:
 - coords = coords_from_vector_1d(v_array)
 - coords = coordinate_system_from_vector_2d_tri(xyz1, xyz2, xyz3)
 - coords = coordinate_system_from_vector_2d_quad(xyz1, xyz2, xyz3, xyz4)
 - coords = coordinate_system_from_vector_2d_tri_theta(xyz1, xyz2, xyz3, theta, dtype=None)
 - coords = cylindrical_rotation_matrix(thetar, dtype='float64')
"""
from __future__ import print_function, absolute_import, division

import numpy as np
from .numpy_utils import cross2d, perpendicular_vector2d
from .matrix3d import axes_stack, normalize_vector2d, dot3d


def coords_from_vector_1d(v_array):
    """
    Gets the coordinate systems for a series of 1D vectors.
    Fakes the j and k axes.

    Parameters
    ----------
    v_array : (n, 3)
        the i vectors

    Returns
    -------
    coords : (n, 3, 3)
        the coordinate systems
    """
    v = np.atleast_2d(v_array)

    i, nmag = normalize_vector2d(v)
    j = perpendicular_vector2d(i)
    k = cross2d(i, j)
    coords = axes_stack(i, j, k, nmag)
    return coords

def coordinate_system_from_vector_2d_tri(xyz1, xyz2, xyz3):
    v21 = xyz2 - xyz1
    v31 = xyz3 - xyz2
    normal = cross2d(v21, v31)
    coords = _coordinate_system_from_vector_2d(v21, normal)
    return coords

def coordinate_system_from_vector_2d_quad(xyz1, xyz2, xyz3, xyz4):
    v21 = xyz2 - xyz1
    v31 = xyz3 - xyz1
    v42 = xyz4 - xyz2
    normal = cross2d(v31, v42)
    coords = _coordinate_system_from_vector_2d(v21, normal)
    return coords

def _coordinate_system_from_vector_2d(v21, normal):
    """
         3
        / \      4----------3
       /   \     |          |
      /     \    |          |
    1--------2   1----------2
    Parameters
    ----------
    theta is defined relative to 1-2, pivoting about 1
    """
    i, nmag = normalize_vector2d(v21)
    k, unused_nmag = normalize_vector2d(normal)
    j = cross2d(k, i)
    coords = axes_stack(i, j, k, nmag)
    return coords

def coordinate_system_from_vector_2d_tri_theta(xyz1, xyz2, xyz3, theta, dtype=None):
    """
         3
        / \
       /   \
      /     \
    1--------2
    theta is defined relative to 1-2, pivoting about 1
    """
    if dtype is None:
        dtype = xyz1.dtype
    dtype = xyz1.dtype
    coords = coordinate_system_from_vector_2d_tri(xyz1, xyz2, xyz3)
    rotation = cylindrical_rotation_matrix(theta, dtype=dtype)
    coords_transformed = dot3d(rotation, coords)
    return coords_transformed

def cylindrical_rotation_matrix(thetar, dtype='float64'):
    """
    Craates a series transformation matrices to rotate by some angle theta

    Parameters
    ----------
    thetar : (n, ) float ndarray
        the theta in radians
    dtype : dtype/str
        the type of the output matrix

    Returns
    -------
    rotation : (ntheta, 3, 3)
        the rotation matrices
    """
    theta = np.asarray(thetar, dtype=dtype)
    ntheta = len(theta)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    #rotation = np.array([
        #[cos(theta), -sin(theta), zero],
        #[sin(theta),  cos(theta), zero],
        #[0.zero,       zero,       one],
    #], dtype=dtype)

    # can this be faster?
    # it definitely could look nicer
    rotation = np.zeros((ntheta, 3, 3), dtype=dtype)
    rotation[:, 0, 0] = cos_theta
    rotation[:, 0, 1] = -sin_theta
    rotation[:, 1, 1] = cos_theta
    rotation[:, 1, 0] = sin_theta
    rotation[:, 2, 2] = 1.

    #print('---------')
    #for rot in rotation:
        #print(np.squeeze(rot))
        #print('---------')
    return rotation
