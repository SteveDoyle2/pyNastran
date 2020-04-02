"""
Various numpy dependent mathematical functions are defined in this file.
This includes:
 - c = cross2d(a, b)
 - row_col_pairs = unique2d(array)
 - row_col_pairs, optional_index, optional_inverse = unique_rows(
       return_index=False, return_inverse=False):
 - augmented_identity(A)

"""
import numpy as np

#ver = np.lib.NumpyVersion(np.__version__)
#if ver < '1.13.0':

def pivot_table(data, rows, cols):
    """PCOMP: rows=element_ids, cols=layer"""
    ncount = len(rows)
    icount = np.arange(ncount)
    assert len(data.shape) in [2, 3], data.shape
    nresults = data.shape[-1]

    rows_new, row_pos_new = np.unique(rows, return_inverse=True)
    cols_new, col_pos_new = np.unique(cols, return_inverse=True)
    nrows = len(rows_new)
    ncols = len(cols_new)

    nshape = len(data.shape)
    if nshape == 3:
        ntimes = data.shape[0]
        shape2 = (ntimes, nrows, ncols, nresults)
    elif nshape == 2:
        shape2 = (nrows, ncols, nresults)
    else:  # pragma: no cover
        raise RuntimeError(nshape)

    pivot_table = np.full((nrows, ncols), -1, dtype='int32')
    pivot_table[row_pos_new, col_pos_new] = icount
    #print(pivot_table)

    ipivot_row, ipivot_col = np.where(pivot_table != -1)
    data2 = np.full(shape2, np.nan, dtype=data.dtype)

    if nshape == 3:
        data2[:, ipivot_row, ipivot_col, :] = data[:, icount, :]
    elif  nshape == 2:
        data2[ipivot_row, ipivot_col, :] = data[icount, :]
    else:  # pragma: no cover
        raise NotImplementedError(nshape)
    return data2, rows_new


def unique2d(a):
    """
    Gets the unique pairs in a 2D vector where the pairs are defined:
    (column 0, column 1).

    Parameters
    ----------
    a : (n,2) ndarray
        the input data

    Returns
    -------
    u : (m,2)
        the unique values in a

    .. note:: this is intended to be used to find unique rows of
              element-id/property-id or property-id/material-id pairs
    .. note:: it works by finding the unique complex numbers and doesn't
              extend well to a 3 column pair
    """
    x, y = a.T
    b = x + y*1.0j
    idx = np.unique(b, return_index=True)[1]
    return a[idx]

#def unique_rows(data):
    #"""
    #finds the unique rows of a numpy array
    #"""
    #uniq = unique(data.view(data.dtype.descr * data.shape[1]))
    #return uniq.view(data.dtype).reshape(-1, data.shape[1])

def duplicates(ids):
    """finds the duplicate ids"""
    counts = np.bincount(ids)
    return np.where(counts > 1)[0]

def is_monotonic(int_array):
    """is the array monotonic?"""
    return np.all(int_array[1:] >= int_array[:-1])

def unique_rows(A, return_index=False, return_inverse=False):
    """
    Similar to MATLAB's unique(A, 'rows'), this returns B, I, J
    where B is the unique rows of A and I and J satisfy
    A = B[J,:] and B = A[I,:]

    Returns
    -------
    I : ndarray
        the index array;
        returns if return_index=True
    J : ndarray
        the inverse array;
        returns if return_inverse=True

    Example (not tested)
    --------------------
    >>> B       = unique_rows(A, return_index=False, return_inverse=False)
    >>> B, I    = unique_rows(A, return_index=True,  return_inverse=False)
    >>> B, J    = unique_rows(A, return_index=False, return_inverse=True)
    >>> B, I, J = unique_rows(A, return_index=True,  return_inverse=True)

    per https://github.com/numpy/numpy/issues/2871
    """
    A = np.require(A, requirements='C')
    assert A.ndim == 2, 'array must be 2D; shape=%s' % str(A.shape)

    B = np.unique(A.view([('', A.dtype)] * A.shape[1]),
                  return_index=return_index,
                  return_inverse=return_inverse)

    if return_index or return_inverse:
        return (B[0].view(A.dtype).reshape((-1, A.shape[1]), order='C'),) \
            + B[1:]
    else:
        return B.view(A.dtype).reshape((-1, A.shape[1]), order='C')

def cross2d(a, b):
    """
    Interface to np.cross for 2d matrices

    [cx0, cy0, cz0] = [ax0, ay0, az0]   [bx0, by0, bz0]
    |cx1, cy1, cz1| = |ax1, ay1, az1| x |bx1, by1, bz1|
    |cx2, cy2, cz2| = |ax2, ay2, az2|   |bx2, by2, bz2|
    [cx3, cy3, cz3] = [ax3, ay3, az3]   [bx3, by3, bz3]
    """
    # axisa=-1, axisb=-1, axisc=-1,
    return np.cross(a, b, axis=1)

def augmented_identity(nx, ny):
    """
    Creates an Identity Matrix augmented with zeros.
    The location of the extra zeros depends on nx/ny.

    .. code-block:: python

       [ 1, 0, 0, 0 ]
       [ 0, 1, 0, 0 ]
       [ 0, 0, 1, 0 ]

    """
    eye = np.eye(max(nx, ny), dtype='float64')
    return eye[:nx, :ny]

def perpendicular_vector(v):
    """
    Finds an arbitrary perpendicular vector to *v*.

    Parameters
    -----------
    v : (3, ) float ndarray
        Vector

    Returns
    --------
    u : (3, ) float ndarray
        Perpendicular vector

    """
    # for two vectors (x, y, z) and (a, b, c) to be perpendicular,
    # the following equation has to be fulfilled
    #     xyz dot abc = 0
    #     0 = ax + by + cz = 0

    # x = y = z = 0 is not an acceptable solution
    if v[0] == v[1] == v[2] == 0.:
        raise ValueError('zero-vector')

    # If one dimension is zero, this can be solved by setting that to
    # non-zero and the others to zero. Example: (4, 2, 0) lies in the
    # x-y-Plane, so (0, 0, 1) is orthogonal to the plane.
    if v[0] == 0:
        return np.array([1., 0., 0.])
    if v[1] == 0:
        return np.array([0., 1., 0.])
    if v[2] == 0:
        return np.array([0., 0., 1.])

    # arbitrarily set a = b = 1
    # then the equation simplifies to
    #     c = -(x + y)/z
    return np.array([1., 1., -1.0 * (v[0] + v[1]) / v[2]])

def perpendicular_vector2d(v_array):
    """
    Finds an array of arbitrary perpendicular vector to *v_array*.

    Parameters
    -----------
    v_array : (3, ) or (n, 3) float ndarray
        Vector

    Returns
    --------
    v_array : (n, 3) float ndarray
        Perpendicular vector

    """
    v = np.atleast_2d(v_array)
    # for two vectors (x, y, z) and (a, b, c) to be perpendicular,
    # the following equation has to be fulfilled
    #     xyz dot abc = 0
    #     0 = ax + by + cz = 0

    # x = y = z = 0 is not an acceptable solution
    v1 = v[:, 0]
    v2 = v[:, 1]
    v3 = v[:, 2]

    v1_equals_0 = v1 == 0.
    v2_equals_0 = v2 == 0.
    v3_equals_0 = v3 == 0.
    ibad = np.where(v1_equals_0 & v2_equals_0 & v3_equals_0)[0]
    if ibad.size:
        raise ValueError('zero-vector')

    vout = np.full(v.shape, np.nan, dtype=None, order='C')
    # If one dimension is zero, this can be solved by setting that to
    # non-zero and the others to zero. Example: (4, 2, 0) lies in the
    # x-y-Plane, so (0, 0, 1) is orthogonal to the plane.
    iv1 = np.where(v1_equals_0)[0]
    iv2 = np.where(v2_equals_0)[0]
    iv3 = np.where(v3_equals_0)[0]
    if iv1.size:
        vout[iv1] = [1., 0., 0.]
    if iv2.size:
        vout[iv2] = [0., 1., 0.]
    if iv3.size:
        vout[iv3] = [0., 0., 1.]

    all_true = v1_equals_0 | v2_equals_0 | v3_equals_0
    #print('iv1 =', iv1)
    #print('iv2 =', iv2)
    #print('iv3 =', iv3)
    #print('all_true =', all_true)
    if np.all(all_true):
        return vout
    is_3d = np.where(~all_true)[0]
    #print('----------')
    #print('is_3d =', is_3d)
    #print(vout)

    # arbitrarily set a = b = 1
    # then the equation simplifies to
    #     c = -(x + y)/z
    vout[is_3d, :2] = 1.
    vout[is_3d, 2] = -1. * (v[is_3d, 0] + v[is_3d, 1]) / v[is_3d, 2]
    return vout

def underflow_norm(x, ord=None, axis=None, keepdims=False):
    """see numpy.linalg.norm"""
    try:
        normi = np.linalg.norm(x, axis=axis)
    except FloatingPointError:
        # the dreaded underflow
        if x.dtype == np.float32:
            x = x.astype('float64')
            normi = np.linalg.norm(x, axis=axis)
        else:  # pragma: no cover
            # the next step would be to nan the min=max depending on the axis
            raise
    return normi
