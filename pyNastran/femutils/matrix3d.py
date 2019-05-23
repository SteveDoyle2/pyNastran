"""
Various numpy dependent 3D matrix functionsare defined in this file.
This includes:
 - mag = norm2d(v)
 - normalize_vector2d(v)
 - ijk = axes_stack(i, j, k, nmag)
 - C = dot_n33_n33(A, B)
 - Tt = transpose3d(T)
 - C = triple(A, T, tranpose=False)

"""
# pylint: disable=C0103
from itertools import count
import numpy as np

def norm2d(v):
    """takes N norms of a (N,3) set of vectors"""
    mag = np.linalg.norm(v, axis=1)
    assert v.shape[0] == len(mag)
    return mag

def normalize_vector2d(v):
    """normalzes a series of (N,3) vectors"""
    mag = norm2d(v)
    nmag = len(mag)
    i = v / mag[:, np.newaxis]
    return i, nmag

def axes_stack(i, j, k, nmag):
    """stack coordinate axes in 3d"""
    i.shape = (nmag, 1, 3)
    j.shape = (nmag, 1, 3)
    k.shape = (nmag, 1, 3)
    ijk = np.hstack([i, j, k])
    return ijk

def dot_33_n33(A, B, debug=True):
    """
    Multiplies a (3x3) matrix by a Nx3x3 matrix

    Parameters
    ----------
    A : (3, 3) float ndarray
        the transformation matrix
    B : (n, 3, 3) float ndarray
        the set of matrices to multiply

    Returns
    -------
    C : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix multiplies
    """
    C = np.matmul(A, B)
    assert A.shape == (3, 3), A.shape
    assert len(B.shape) == 3, B.shape
    assert B.shape[1:] == (3, 3), B.shape
    if debug:
        dtype = A.dtype
        #print('------------------------')
        D = np.zeros(B.shape, dtype=dtype)
        print('A.shape =', A.shape)
        for i, Bi in zip(count(), B):
            print('Bi.shape =', Bi.shape)
            ABi = A.dot(Bi)
            print('A @ Bi.shape =', ABi.shape)
            D[i, :, :] = ABi
            #print(D[i, :, :])
            #print('------------------------')
        assert np.all(np.allclose(C, D))
    return D

def dot_n33_33(A, B, debug=True):
    """
    Multiplies a (3x3) matrix by a Nx3x3 matrix

    Parameters
    ----------
    A : (n, 3, 3) float ndarray
        the set of matrices to multiply
    B : (3, 3) float ndarray
        the transformation matrix

    Returns
    -------
    C : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix multiplies
    """
    assert len(A.shape) == 3, A.shape
    assert A.shape[1:] == (3, 3), A.shape
    assert B.shape == (3, 3), B.shape
    #C = np.matmul(B, A.T) # 3x3x4
    #C = B.dot(A)  # 3x4x3
    C = np.tensordot(A, B, axes=(1))

    #print('C (nx3x3 @ 33).shape = ', C.shape)
    if debug:
        dtype = A.dtype
        #print('------------------------')
        #print('dot_n33_33: A.shape=%s; B.shape=%s' % (str(A.shape), str(B.shape)))
        D = np.zeros(A.shape, dtype=dtype)
        for i, Ai in zip(count(), A):
            AiB = Ai.dot(B)
            #print(AiB.shape)
            D[i, :, :] = AiB
            #print(D[i, :, :])
            #print('------------------------')
        assert np.all(np.allclose(C, D))
    return C

def dot_n33_n33(A, B, debug=True):
    """
    Multiplies two matrices together

    Parameters
    ----------
    A, B : (n, 3, 3) float ndarray
        the set of matrices to multiply

    Returns
    -------
    C : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix multiplies
    """
    assert len(A.shape) == 3, A.shape
    assert A.shape[1:] == (3, 3), A.shape
    assert len(B.shape) == 3, B.shape
    assert B.shape[1:] == (3, 3), B.shape
    #C = np.matmul(A, B)
    if debug:
        dtype = A.dtype
        #print('------------------------')
        D = np.zeros(A.shape, dtype=dtype)
        for i, Ai, Bi in zip(count(), A, B):
            D[i, :, :] = Ai.dot(Bi)
            #print(D[i, :, :])
            #print('------------------------')
    #if not np.all(np.allclose(C, D)):
        #print('C:\n%s'% C)
        #print('D:\n%s'% D)
    return D

def dot_n33_n3(A, B, debug=True):
    """
    Multiplies two N x 3 x 3 matrices together

    Parameters
    ----------
    A : (n, 3, 3) float ndarray
        the first matrix
    B : (n, 3) float ndarray
        the second matrix

    Returns
    -------
    C : (n, 3) float ndarray
        the set of 3 x 3 matrix multiplies
    """
    dtype = A.dtype
    D = np.zeros(B.shape, dtype=dtype)
    for i, Ai, Bi in zip(count(), A, B):
        #print('------------')
        #print('dot_n33_n3')
        #print(Ai)
        #print(Bi)
        D[i, :] = Ai.dot(Bi)
    return D

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
    return np.transpose(T, axes=(0, 2, 1))

def triple_n33_n33(A, T, tranpose=False, debug=True):
    """
    Calculates the matrix triple product  for a series of::

        triple[n, :, :] = T.T @ A @ T  # tranpose=False
        triple[n, :, :] = T @ A @ T.T  # transpose=True

    Parameters
    ----------
    A, T : (n, 3, 3)
        the set of matrices to multiply
    tranpose : bool; default=False
        transposes T

    Returns
    -------
    triple : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix triples

    """
    assert A.shape == T.shape, 'A.shape=%s T.shape=%s' % (str(A.shape), str(T.shape))
    if tranpose:
        C = np.matmul(transpose3d(T), np.matmul(A, T))
    else:
        C = np.matmul(T, np.matmul(A, transpose3d(T)))

    if debug:
        D = np.full(A.shape, np.nan)
        if tranpose:
            for i, Ai, Ti in zip(count(), A, T):
                Dia = Ti.T.dot(Ai.dot(Ti))
                D[i, :, :] = Dia
        else:
            for i, Ai, Ti in zip(count(), A, T):
                Dib = Ti.dot(Ai.dot(Ti.T))
                D[i, :, :] = Dib
        assert np.all(np.allclose(C, D)), 'tranpose=%s' % tranpose
    return C

def triple_n33_33(A, T, tranpose=False, debug=True):
    """
    Calculates the matrix triple product  for a series of::

        triple[n, :, :] = T.T @ A @ T  # tranpose=False
        triple[n, :, :] = T @ A @ T.T  # transpose=True

    Parameters
    ----------
    A, T : (n, 3, 3)
        the set of matrices to multiply
    tranpose : bool; default=False
        transposes T

    Returns
    -------
    triple : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix triples

    """
    assert list(A.shape)[1:] == list(T.shape), 'A.shape=%s T.shape=%s' % (str(A.shape), str(T.shape))
    if tranpose:
        C = np.matmul(T.T, np.matmul(A, T))
    else:
        C = np.matmul(T, np.matmul(A, T.T))

    if debug:
        D = np.full(A.shape, np.nan)
        if tranpose:
            for i, Ai in zip(count(), A):
                Dia = T.T.dot(Ai.dot(T))
                D[i, :, :] = Dia
        else:
            for i, Ai, Ti in zip(count(), A, T):
                Dib = T.dot(Ai.dot(T.T))
                D[i, :, :] = Dib
        assert np.all(np.allclose(C, D)), 'tranpose=%s' % tranpose
    return C
