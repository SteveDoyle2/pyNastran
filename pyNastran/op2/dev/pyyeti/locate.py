# -*- coding: utf-8 -*-
"""
Tools for locating data or subarrays inside other arrays.

@author: Tim Widrick
"""

import numpy as np


def findvals(m, v):
    """
    Get partition vector for all occurrences of all values in `v` in
    `m`.

    Parameters
    ----------
    m : array
        Array to be searched.
    v : array
        Array of values to find in m.

    Returns
    -------
    pv : 1d ndarray
        Values are indexes into `m` of any value in `v`.  Will be
        empty if `m` has none of the values in `v`.

    `m` is flattened to 1d before searching (using column-major
    ordering 'F').  The values in `pv` correspond to::
          [  0      r  ...
             1    r+1
           ...    ...
           r-1   2r-1  ... r*c-1 ]  where m is r x c

    Examples
    --------
    >>> import numpy as np
    >>> import locate
    >>> m = np.array([[10, 20], [30, 20]])
    >>> locate.findvals(m, 20)
    array([2, 3])
    >>> locate.findvals(m, 30)
    array([1])
    >>> locate.findvals(m, 100)
    array([], dtype=int64)

    """
    m = np.atleast_1d(m)
    v = np.atleast_1d(v)
    m = m.flatten(order='F')
    v = v.flatten()
    pv = np.zeros(len(m), dtype=bool)
    for i in range(len(v)):
        pv |= m == v[i]
    return pv.nonzero()[0]


def find_subsequence(seq, subseq):
    """
    Returns indices of where subseq occurs in seq.  Both are 1d numpy
    arrays.

    Parameters
    ----------
    seq : array
        1D array to search in.
    subseq : array
        1D array to search for.

    Returns
    -------
    pv : array
        1D numpy array of indices:
         - length will be equal to the number of occurrences of subseq
         - the indices are to the start of each subseq in seq
        Will be empty if subseq is not found in seq.

    Examples
    --------
    >>> import locate
    >>> a = [1, 2, 3, 4, 5, 6, 2, 3]
    >>> sub = [2, 3]
    >>> locate.find_subsequence(a,sub)
    array([1, 6])
    >>> locate.find_subsequence(a,[6, 5])
    array([], dtype=int64)

    """
    seq = np.asarray(seq).reshape(-1)
    subseq = np.asarray(subseq).reshape(-1)
    target = subseq @ subseq
    candidates = np.where(np.correlate(seq, subseq,
                                       mode='valid') == target)[0]
    # some of the candidates entries may be false positives; check:
    check = candidates[:, np.newaxis] + np.arange(len(subseq))
    mask = np.all((np.take(seq, check) == subseq), axis=-1)
    return candidates[mask]


def find_rows(matrix, row):
    """
    Returns indices of where row occurs in matrix.

    Parameters
    ----------
    matrix : array
        2d numpy array.
    row : array
        1d numpy array.

    Returns
    -------
    pv : array
        A 1d numpy array of row indices.  Will be empty if row is not
        found or if length(row) != cols(matrix).

    Examples
    --------
    >>> import numpy as np
    >>> import locate
    >>> mat = np.array([[7, 3], [6, 8], [4, 0],
    ...                 [9, 2], [1, 5], [6, 8]])
    >>> locate.find_rows(mat,np.array([1, 2]))
    array([], dtype=int64)
    >>> pv = locate.find_rows(mat,np.array([6, 8]))
    >>> pv
    array([1, 5])
    >>> mat[pv,:]
    array([[6, 8],
           [6, 8]])

    """
    (r1, c1) = np.shape(matrix)
    c2 = len(row)
    if c1 != c2:
        return np.array([], dtype=int)
    i = find_subsequence(matrix.flatten('C'), row)
    pv = np.mod(i, c1) == 0
    return i[pv] // c1


def get_intersection(D1, D2, keep=0):
    """
    Get row intersection partition vectors between two matrices or
    vectors.

    Parameters
    ----------
    D1 : array
        1d or 2d array.
    D2 : array
        1d or 2d array.
    keep : integer
        0, 1 or 2:
           - if 0, rows are only swapped in the larger matrix
           - if 1, rows are only swapped in D2
           - if 2, rows are only swapped in D1

    Returns
    -------
    tuple: (pv1, pv2)
        pv1 : array
            Row index vector into D1.
        pv2 : array
            Row index vector into D2.

    `pv1` and `pv2` are found such that:
    ::
        D1[pv1] == D2[pv2]
        (Note for matrices:  M[i] == M[i, :])

    For matrices, the number of columns in D1 and D2 must be equal to
    get non-empty results.

    Examples
    --------
    >>> import numpy as np
    >>> import locate
    >>> mat1 = np.array([[7, 3], [6, 8], [4, 0], [9, 2], [1, 5]])
    >>> mat2 = np.array([[9, 2], [1, 5], [7, 3]])
    >>> pv1, pv2 = locate.get_intersection(mat1, mat2)
    >>> pv1
    array([3, 4, 0])
    >>> pv2
    array([0, 1, 2])
    >>> np.all(mat1[pv1] == mat2[pv2])
    True
    >>> locate.get_intersection(mat1, mat2, 1)
    (array([0, 3, 4]), array([2, 0, 1]))
    >>> locate.get_intersection(mat2, mat1, 2)
    (array([2, 0, 1]), array([0, 3, 4]))
    >>> locate.get_intersection(mat2, mat1)
    (array([0, 1, 2]), array([3, 4, 0]))
    >>> mat3 = np.array([[1,2,3]])
    >>> locate.get_intersection(mat1, mat3)
    (array([], dtype=int64), array([], dtype=int64))

    """
    D1 = np.array(D1)
    D2 = np.array(D2)
    if D1.ndim == D2.ndim == 1:
        c1 = c2 = 1
        r1 = len(D1)
        r2 = len(D2)
        D1 = np.atleast_2d(D1)
        D2 = np.atleast_2d(D2)
        D1 = D1.T
        D2 = D2.T
    else:
        D1 = np.atleast_2d(D1)
        D2 = np.atleast_2d(D2)
        (r1, c1) = np.shape(D1)
        (r2, c2) = np.shape(D2)
    if c1 != c2:
        return np.array([], dtype=int), np.array([], dtype=int)
    # loop over the smaller one; index into 'd2' can be in any order
    if r1 <= r2:
        r = r1
        d1 = D1
        d2 = D2
        switch = False
    else:
        r = r2
        d1 = D2
        d2 = D1
        switch = True
    pv1 = np.zeros(r, dtype=int)
    pv2 = np.zeros(r, dtype=int)
    j = 0
    for i in range(r):
        l = find_rows(d2, d1[i])
        if l.size > 0:
            pv1[j] = i
            pv2[j] = l[0]
            j += 1
    if j == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    if switch:
        t = pv1[:j]
        pv1 = pv2[:j]
        pv2 = t
    else:
        pv1 = pv1[:j]
        pv2 = pv2[:j]
    if switch and keep == 1:
        si = pv1.argsort()
        return pv1.take(si), pv2.take(si)
    elif not switch and keep == 2:
        si = pv2.argsort()
        return pv1.take(si), pv2.take(si)
    return pv1, pv2


def find2zo(pv, n):
    """
    Return a True/False vector of length n where the True values are
    located according to pv.

    Example:
    >>> import numpy as np
    >>> import locate
    >>> pv = np.array([0,3,5])
    >>> locate.find2zo(pv,8)
    array([ True, False, False,  True, False,  True, False, False], dtype=bool)

    """
    tf = np.zeros(n, dtype=bool)
    tf[pv] = True
    return tf


def flippv(pv, n):
    """Flips the meaning of an index partition vector.

    Parameters
    ----------
    pv : ndarray
        The index partition to flip.
    n : integer
        The length of the dimension to partition.

    Returns
    -------
    notpv : ndarray
        The complement of pv.

    Example:
    >>> import numpy as np
    >>> import locate
    >>> pv = np.array([0,3,5])
    >>> locate.flippv(pv,8)
    array([1, 2, 4, 6, 7])

    """
    tf = np.ones(n, dtype=bool)
    tf[pv] = False
    return tf.nonzero()[0]


def list_intersection(L1, L2):
    """
    Get list intersection partition vectors between two lists

    Parameters
    ----------
    L1 : list
       List 1; the output vectors maintain the order of `L1`.
    L1 : list
       List 2.

    Returns
    -------
    tuple: (pv1, pv2)
       pv1 : ndarray.  Index vector into L1.

       pv2 : ndarray.  Index vector into L2.

    `pv1` and `pv2` are found such that:
      [L1[i] for i in pv1] == [L2[i] for i in pv2]

    Examples
    --------
    >>> import locate
    >>> pv1, pv2 = locate.list_intersection(['a', 3, 'z', 0],
    ...                                     [0, 'z', 1, 'a'])
    >>> pv1
    array([0, 2, 3])
    >>> pv2
    array([3, 1, 0])
    >>> locate.list_intersection(['a', 'b'], [1, 2])
    (array([], dtype=int64), array([], dtype=int64))

    """
    inters = set(L1) & set(L2)
    r = len(inters)
    if r == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    pv1 = np.zeros(r, dtype=int)
    pv2 = np.zeros(r, dtype=int)
    j = 0
    for j, item in enumerate(inters):
        pv1[j] = L1.index(item)
        pv2[j] = L2.index(item)
    si = pv1.argsort()
    return pv1[si], pv2[si]


if __name__ == '__main__':
    mat1 = np.array([[7, 3], [6, 8], [4, 0], [9, 2], [1, 5]])
    mat2 = np.array([[9, 2], [1, 5], [7, 3]])
    pv1, pv2 = get_intersection(mat1, mat2)
    assert np.all(np.array([3, 4, 0]) == pv1)
    assert np.all(np.array([0, 1, 2]) == pv2)
    pv1, pv2 = get_intersection(mat1, mat2, 1)
    assert np.all(np.array([0, 3, 4]) == pv1)
    assert np.all(np.array([2, 0, 1]) == pv2)
    pv1, pv2 = get_intersection(mat1, mat2, 2)
    assert np.all(np.array([3, 4, 0]) == pv1)
    assert np.all(np.array([0, 1, 2]) == pv2)
    pv = np.array([0, 3, 5])
    tf = find2zo(pv, 8)
    assert np.all(np.array([True, False, False, True, False, True,
                            False, False]) == tf)

    import doctest
    doctest.testmod()
