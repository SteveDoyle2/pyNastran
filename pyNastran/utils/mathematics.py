# -*- coding: utf-8 -*-
# pylint: disable=C0103
"""
Various mathematical functions are defined in this file.  This includes:
 - gauss(n)
 - get_abs_index(data, axis=1)
 - get_abs_max(min_values, max_values)
 - get_max_index(data, axis=1)
 - get_min_index(data, axis=1)
 - integrate_positive_unit_line(x, y, min_value=0.)
 - integrate_unit_line(x, y)
 - is_float_ranged(a, x, b)
 - is_list_ranged(a, List, b)
 - list_print(list_a, tol=1e-8, float_fmt='%-3.2g', zero_fmt='    0')
 - print_annotated_matrix(A, row_names=None, col_names=None, tol=1e-8)
 - print_matrix(A, tol=1e-8)
 - reduce_matrix(matrix_a, nids)
 - roundup(value, round_increment=100)
 - solve_tridag(A, D)
 - unique2d(a)

All beams are LineProperty objects.
Multi-segment beams are IntegratedLineProperty objects.

"""
from math import sqrt, ceil

from numpy import (float32, float64, complex64, complex128, array, cross,
                   allclose, argmax, argmin, arange)
import numpy as np
from numpy.linalg import norm  # type: ignore

#from scipy.linalg import solve_banded  # type: ignore
from scipy.integrate import quad  # type: ignore

# should future proof this as it handles 1.9.0.dev-d1dbf8e, 1.10.2, and 1.6.2
#_numpy_version = [int(i) for i in numpy.__version__.split('.') if i.isdigit()]

# def vectorized_searchsorted_eq(eids_all, eids):
    # """
    # Vectorizes where to find all locations for values in array

    # TODO: there has to be a better function to do this...
    # """
    # #i = zeros(eids.shape, dtype='int32')
    # i = []
    # for eid in eids:
        # i.append(where(eids_all == eid)[0])
    # return hstack(i)

def get_abs_max(min_values, max_values, dtype='float32'):
    """Get return the value with the greatest magnitude, preserving sign."""
    min_values = np.asarray(min_values)
    max_values = np.asarray(max_values)
    imin = np.abs(min_values) > np.abs(max_values)
    out = np.zeros(min_values.shape, dtype=dtype)
    out[imin] = min_values[imin]
    out[~imin] = max_values[~imin]
    return out

    #nvalues = len(min_values)
    #data = array([min_values, max_values], dtype=dtype)
    #i = argmax(abs(data), axis=0)
    #assert len(i) == nvalues
    #k = arange(nvalues, dtype='int32')
    #return data[i[:], k]


def get_abs_index(data, axis=1):
    """
    Gets the maximum absolute value of a 2D matrix along an axis

    Examples
    --------
    >>> data = [
            [4.0, 2.2, 3.0, 5.0, 2.2]  # subcase 1
            [4.1, 2.1, 3.1, 5.1, 2.1], # subcase 2
        ]
    >>> max_values, index = get_min_index(data, axis=1)
    >>> out
    [4.1, 2.2, 3.1, 5.1, 2.2]

    >>> index
    [1, 0, 1, 1, 0]
    """
    nvalues = data.shape[axis]
    # isubcase, nelements
    axis2 = abs(axis - 1)
    i = argmax(abs(data), axis=axis2)
    assert len(i) == nvalues, 'data.shape=%s len(i)=%s nvalues=%s' % (str(data.shape), len(i), nvalues)
    k = arange(nvalues, dtype='int32')
    return data[i[:], k], i

def get_max_index(data, axis=1):
    """
    Gets the maximum values of a 2D matrix along an axis

    Examples
    --------
    >>> data = [
            [4.0, 2.2, 3.0, 5.0, 2.2]  # subcase 1
            [4.1, 2.1, 3.1, 5.1, 2.1], # subcase 2
        ]
    >>> max_values, index = get_max_index(data, axis=1)
    >>> out
    [4.1, 2.2, 3.1, 5.1, 2.2]

    >>> index
    [1, 0, 1, 1, 0]
    """
    nvalues = data.shape[axis]
    # isubcase, nelements
    axis2 = abs(axis - 1)
    i = argmax(data, axis=axis2)
    assert len(i) == nvalues, 'data.shape=%s len(i)=%s nvalues=%s' % (str(data.shape), len(i), nvalues)
    k = arange(nvalues, dtype='int32')
    return data[i[:], k], i

def get_min_index(data, axis=1):
    """
    Gets the minimum values of a 2D matrix along an axis

    Examples
    --------
    >>> data = [
            [4.0, 2.2, 3.0, 5.0, 2.2]  # subcase 1
            [4.1, 2.1, 3.1, 5.1, 2.1], # subcase 2
        ]
    >>> min_values, index = get_min_index(data, axis=1)
    >>> out
    [4.0, 2.1, 3.0, 5.0, 2.1]

    >>> index
    [0, 1, 0, 0, 1]
    """
    nvalues = data.shape[axis]
    axis2 = abs(axis - 1)
    i = argmin(data, axis=axis2)
    assert len(i) == nvalues, 'data.shape=%s len(i)=%s nvalues=%s' % (str(data.shape), len(i), nvalues)
    k = arange(nvalues, dtype='int32')
    return data[i[:], k], i


def integrate_unit_line(x, y):
    """
    Integrates a line of length 1.0 by linear interpolation

    Parameters
    ----------
    x : List[float]
        the independent variable
    y : List[float]
        the dependent variable

    Returns
    -------
    integrated_value : float
        the area under the curve
    """
    if len(set(y)) == 1:
        return y[0]  # (x1-x0 = 1., so yBar*1 = yBar)
    try:
        assert len(x) == len(y), 'x=%s y=%s' % (x, y)
        # integrate the area; y=f(x); A=integral(y*dx,x)

        #f = np.interp(_xi, x, y, left=y[0], right=y[-1])
        out = quad(np.interp, 0., 1., args=(x, y, y[0], y[-1]))
    except:
        # print('spline Error x=%s y=%s' % (x, y))
        raise
    return out[0]




def integrate_positive_unit_line(x, y, min_value=0.):
    """
    Integrates a line of length 1.0 by linear interpolation

    Parameters
    ----------
    x : List[float]
        the independent variable
    y : List[float]
        the dependent variable
    min_value : float; default=0.0
        ???

    Returns
    -------
    integrated_value : float
        the area under the curve
    """

    for i, yi in enumerate(y):
        if yi < min_value:
            raise ValueError('y%i=%s and must be greater than %s' % (i+1, yi, min_value))
    return integrate_unit_line(x, y)


def reduce_matrix(matrix_a, nids):
    """
    takes a list of ids and removes those rows and cols
    """
    nrows = len(nids)
    matrix_b = np.zeros((nrows, nrows), dtype='float64')
    for i, irow in enumerate(nids):
        for j, jcol in enumerate(nids):
            matrix_b[i, j] = matrix_a[irow, jcol]
    return matrix_b


def is_list_ranged(a, List, b):
    """
    Returns true if a<= x <= b or a-x < 0 < b-x

    Parameters
    ----------
    a : float
        the lower bound value (inclusive)
    x : List[float, ...]
        the search values
    b: float
        the upper bound value (inclusive)

    Returns
    -------
    is_ranged : bool
        True/False
    """
    for x in List:
        if not is_float_ranged(a, x, b):
            return False
    return True


def is_float_ranged(a, x, b):
    """
    Returns true if a<= x <= b or a-x < 0 < b-x.

    Parameters
    ----------
    a : float
        the lower bound value (inclusive)
    x : List[float, ...]
        the search values
    b: float
        the upper bound value (inclusive)

    Returns
    -------
    is_ranged : bool
        True/False
    """
    if (not a < x) and (not allclose(x, a)):
        return False

    if (not x < b) and (not allclose(x, b)):
        return False
    return True


def print_annotated_matrix(A, row_names=None, col_names=None, tol=1e-8):
    """
    Takes a list/dictionary and annotates the row number with that value
    indicies go from 0 to N
    """
    B = array(A)
    if row_names is None:
        row_names = [i for i in range(B.shape[0])]

    rwidth = max([len(str(row_names[i])) for i in range(len(row_names))])
    row_fmt = '%%-%ss' % rwidth

    header = ''
    if col_names is not None:
        col_name_list = [str(col_names[i]) for i in col_names]
        cwidth = max([len(name) for name in col_name_list])

        cwidth = 5
        col_fmt = '%%-%ss ' % cwidth
        #print("col_fmt = ", col_fmt)
        header = row_fmt % '' + '   ' + col_fmt * len(col_names) % tuple(col_name_list) + '\n'
        float_fmt = '%%-%i.2f' % cwidth

    c = header + ''.join([row_fmt % (str(row_names[i])) + ' ' +
                          list_print(B[i, :], tol, float_fmt=float_fmt)
                          + '\n' for i in range(B.shape[0])])
    return c


def print_matrix(A, tol=1e-8):
    """prints a 2d matrix in a readable format"""
    B = array(A)
    return ''.join([list_print(B[i, :], tol) + '\n' for i in range(B.shape[0])])


def list_print(list_a, tol=1e-8, float_fmt='%-3.2g', zero_fmt='    0'):
    """prints a list / numpy array in a readable format"""
    if len(list_a) == 0:
        return '[]'

    def _print(a):
        if isinstance(a, str):
            return a
        for i, j in ((float, float_fmt), (float32, float_fmt),
                     (float64, float_fmt), (int, '%3i')):
            if isinstance(a, i):
                if abs(a) < tol:
                    return zero_fmt
                return j % (0. if abs(a) < tol else a)

        if isinstance(a, (complex, complex64, complex128)):
            return '%4s%4s' % ('0' if abs(a.real) < 1e-8 else '%.4g' % (a.real),
                               '' if abs(a.imag) < 1e-8 else '%+.4gj' % (a.imag))
        try:
            print("list_print: type(a) is not supported... %s" % (type(a)))
            return ' %g' % a
        except TypeError:
            print("a = %r" % a)
            raise

    return '[ '+ ', '.join([_print(a) for a in list_a])+ ']'


#def solve_tridag(A, D):
    #"""
    #Solves a tridagonal matrix [A]{x}={b} for {x}

    #Parameters
    #----------
    #A : (N,) float ndarray
        #main diagonal
    #D : (N-1,) float ndarray)
        #off diagonal

    #Returns
    #-------
    #x : (N, )
        #the result
    #"""
    ## Find the diagonals
    #ud = insert(diag(A, 1), 0, 0)  # upper diagonal
    #d = diag(A)  # main diagonal
    #ld = insert(diag(A, -1), len(d) - 1, 0)  # lower diagonal
    ## simplified matrix
    #ab = np.matrix([ud, d, ld])
    #return solve_banded((1, 1), ab, D, overwrite_ab=True, overwrite_b=True)


Area = lambda a, b: 0.5 * norm(cross(a, b))

def gauss(n):  # pragma: no cover
    r"""
    A quadrature rule: an approximation of the definite integral of a function.
    Currently implementation supports up to 5 quadrature points.

    Function returns following values depending on n (number of points):

    * n = 1:

     * \f$ 0 \f$ --> \f$ 2 \f$

    * n = 2:

     * \f$ \pm 1/\sqrt{3} \f$ --> \f$ 1 \f$

    * n = 3

     * \f$ 0 \f$   --> \f$ 8/9 \f$
     * \f$ \pm\sqrt{3/5} \f$ --> \f$ 5/9 \f$

    * n = 4:

     * \f$ \pm\sqrt{\left( 3 - 2\sqrt{6/5} \right)/7} \f$ --> \f$ (18+\sqrt{30})/36 \f$
     * \f$ \pm\sqrt{\left( 3 + 2\sqrt{6/5} \right)/7} \f$ --> \f$ (18-\sqrt{30})/36 \f$

    * n = 5:

     - \f$ 0 \f$ --> \f$ 128/225 \f$
     - \f$ \pm\frac{1}{3}\sqrt{5-2\sqrt{10/7}} \f$ --> \f$ (322+13\sqrt{70})/900 \f$
     - \f$ \pm\frac{1}{3}\sqrt{5+2\sqrt{10/7}} \f$ --> \f$ (322-13\sqrt{70})/900 \f$


    :param n:       Number of quadrature points
    :returns lists: points and corresponding weights, sorted by points value

    .. seealso:: http://en.wikipedia.org/wiki/Gaussian_quadrature"""
    if n == 1:
        points = [0.]
        weights = [2.]
    elif n == 2:
        p = 1. / sqrt(3)
        points = [-p, p]
        weights = [1., 1.]
    elif n == 3:
        p = sqrt(3 / 5.)
        points = [-p, 0., p]
        weights = [5 / 9., 8 / 9., 5 / 9.]
    elif n == 4:
        p1 = (3 - 2. * sqrt(6 / 5)) / 7.
        p2 = (3 + 2. * sqrt(6 / 5)) / 7.
        w1 = (18 + sqrt(30)) / 36.
        w2 = (18 - sqrt(30)) / 36.
        points = [-p2, -p1, p1, p2]
        weights = [w2, w1, w1, w2]
    elif n == 5:
        p1 = 1 / 3. * sqrt(5 - 2 * sqrt(10. / 7.))
        p2 = 1 / 3. * sqrt(5 + 2 * sqrt(10. / 7.))
        w1 = (322 + 13 * sqrt(70)) / 900.
        w2 = (322 - 13 * sqrt(70)) / 900.
        points = [-p2, -p1, 0, p1, p2]
        weights = [w2, w1, 128 / 225., w1, w2]
    else:
        raise NotImplementedError('The current implementation only supports up to '
                                  '5 quadrature points')
    return points, weights

def roundup(value, round_increment=100):
    """
    Rounds up to the next N.

    Parameters
    ----------
    value : int
        the value to round up
    round_increment : int
        the increment to round by

    .. python

      >>> 100 = roundup(10)
      >>> 200 = roundup(105)
      >>> 300 = roundup(200)
      >>> 1000 = roundup(200, 1000)
      >>> 2000 = roundup(1000, 1000)
      >>> 2000 = roundup(1001, 1000)

    .. note :: this function is used to ensure that renumbering is more
               obvious when testing
    """
    return int(ceil(value / float(round_increment))) * round_increment
