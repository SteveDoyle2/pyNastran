"""
Various numpy dependent 3D matrix functionsare defined in this file.
This includes:
 - mag = norm2d(v)
 - normalize_vector2d(v)
 - ijk = axes_stack(i, j, k, nmag)
 - C = dot3d(A, B)
 - Tt = transpose3d(T)
 - C = triple(A, T)
 - C = triple_transpose(A, T)
"""
from __future__ import print_function
from itertools import count
import numpy as np

def norm2d(v):
    mag = np.linalg.norm(v, axis=1)
    assert v.shape[0] == len(mag)
    return mag

def normalize_vector2d(v):
    mag = norm2d(v)
    nmag = len(mag)
    i = v / mag[:, np.newaxis]
    return i, nmag

def axes_stack(i, j, k, nmag):
    i.shape = (nmag, 1, 3)
    j.shape = (nmag, 1, 3)
    k.shape = (nmag, 1, 3)
    ijk = np.hstack([i, j, k])
    return ijk

def dot3d(A, B, debug=True):
    """
    Multiplies two N x 3 x 3 matrices together

    Parameters
    ----------
    A, B : (n, 3, 3) float ndarray
        the set of matrices to multiply

    Returns
    -------
    C : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix multiplies
    """
    C = np.matmul(A, B)
    if debug:
        dtype = A.dtype
        #print('------------------------')
        D = np.zeros(A.shape, dtype=dtype)
        for i, Ai, Bi in zip(count(), A, B):
            D[i, :, :] = Ai.dot(Bi)
            #print(D[i, :, :])
            #print('------------------------')
    assert np.all(np.allclose(C, D))
    return C

def transpose3d(T):
    """
    Returns the transpose in 3d

    Parameters
    ----------
    A : (n, a, b)
        the set of matrices to transpose

    Returns
    -------
    transpose : (n, b, a) float ndarray
        the transposed matrix
    """
    return np.transpose(T, axes=(1, 2))

def triple(A, T):
    """
    Calculates the matrix triple product  for a series of::

        triple[n, :, :] = T.T @ A @ T

    Parameters
    ----------
    A, T : (n, 3, 3)
        the set of matrices to multiply

    Returns
    -------
    triple : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix triples

    """
    return np.matmul(transpose3d(T), np.matmul(A, T))

def triple_transpose(A, T):
    """
    Calculates the matrix triple product for a series of::

        triple[n, :, :] = T @ A @ T.T

    Parameters
    ----------
    A, T : (n, 3, 3)
        the set of matrices to multiply

    Returns
    -------
    triple : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix triples

    """
    return np.matmul(T, np.matmul(A, transpose3d(T)))
