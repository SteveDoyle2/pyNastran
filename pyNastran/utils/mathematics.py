# -*- coding: utf-8 -*-
# pylint: disable=C0103
from __future__ import print_function
from six import string_types
from six.moves import range

from math import sqrt
from numpy import (float32, float64, complex64, complex128, array, cross,
                   allclose, zeros, matrix, insert, diag, eye, argmax, argmin, arange)
from numpy.linalg import norm

from scipy.linalg import solve_banded
from scipy.interpolate import splrep, splev
from scipy.integrate import quad


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

def get_abs_max(min_values, max_values):
    """Get return the value with the greatest magnitude, preserving sign."""
    nvalues = len(min_values)
    data = array([min_values, max_values], dtype='float32')
    i = argmax(abs(data), axis=0)
    assert len(i) == nvalues
    # return data[i, :]
    k = arange(nvalues, dtype='int32')
    return data[i[:], k]


def get_abs_index(data, axis=1):
    """
    Gets the maximum absolute value of a 2D matrix along an axis

    Example
    -------
    >>> data = [
            [4.0, 2.2, 3.0, 5.0, 2.2]  # subcase 1
            [4.1, 2.1, 3.1, 5.1, 2.1], # subcase 2
        ]
    >>> max_values, index = get_min_index(data, axis=1)
    >>> out   = [4.1, 2.2, 3.1, 5.1, 2.2]
    >>> index = [1, 0, 1, 1, 0]
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

    Example
    -------
    >>> data = [
            [4.0, 2.2, 3.0, 5.0, 2.2]  # subcase 1
            [4.1, 2.1, 3.1, 5.1, 2.1], # subcase 2
        ]
    >>> max_values, index = get_min_index(data, axis=1)
    >>> out   = [4.1, 2.2, 3.1, 5.1, 2.2]
    >>> index = [1, 0, 1, 1, 0]
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

    Example
    -------
    >>> data = [
            [4.0, 2.2, 3.0, 5.0, 2.2]  # subcase 1
            [4.1, 2.1, 3.1, 5.1, 2.1], # subcase 2
        ]
    >>> min_values, index = get_min_index(data, axis=1)
    >>> out   = [4.0, 2.1, 3.0, 5.0, 2.1]
    >>> index = [0, 1, 0, 0, 1]
    """
    nvalues = data.shape[axis]
    axis2 = abs(axis - 1)
    i = argmin(data, axis=axis2)
    assert len(i) == nvalues, 'data.shape=%s len(i)=%s nvalues=%s' % (str(data.shape), len(i), nvalues)
    k = arange(nvalues, dtype='int32')
    return data[i[:], k], i


def integrate_line(x, y):
    """
    Integrates a line of length 1.0

    Parameters
    ----------
    x : List[float]
        the independent variable
    y : List[float]
        the dependent variable

    :returns integrated_value: the area under the curve
    """
    if len(set(y)) == 1:
        return y[0]  # (x1-x0 = 1., so yBar*1 = yBar)
    try:
        assert len(x) == len(y), 'x=%s y=%s' % (x, y)
        # integrate the area; y=f(x); A=integral(y*dx,x)
        out = quad(splev, 0., 1., args=(build_spline(x, y)))
    except:
        # print('spline Error x=%s y=%s' % (x, y))
        raise
    return out[0]


def build_spline(x, y):
    """
    Builds a cubic spline or 1st order spline if there are less than 3 terms

    Parameters
    ----------
    x : List[float]
        the independent variable
    y : List[float]
        the dependent variable

    :returns splrep: a splrep object (linear or cubic spline depending
                     on the length of x)

    .. note:: a 1st order spline is the same as linear interpolation
    """
    # build a linearly interpolated representation or cubic one
    return splrep(x, y, k=1) if len(x) < 3 else splrep(x, y)


def integrate_positive_line(x, y, minValue=0.):
    """
    Integrates a line of length 1.0

    Parameters
    ----------
    x : List[float]
        the independent variable
    y : List[float]
        the dependent variable

    :returns integrated_value: the area under the curve
    """
    if len(set(y)) == 1:
        return y[0]  # (x1-x0 = 1., so yBar*1 = yBar)
    try:
        assert len(x) == len(y), 'x=%s y=%s' % (x, y)
        # now integrate the area
        eval_posit_spline = lambda x, spl, min_val: max(splev([x], spl), min_val)
        out = quad(eval_posit_spline, 0., 1., args=(build_spline(x, y), minValue))
    except:
        raise RuntimeError('spline Error x=%s y=%s' % (x, y))
    return out[0]


def reduce_matrix(matA, nids):
    """
    takes a list of ids and removes those rows and cols
    """
    nRows = len(nids)
    matB = matrix(zeros((nRows, nRows), dtype='float64'))
    for i, irow in enumerate(nids):
        for j, jcol in enumerate(nids):
            matB[i, j] = matA[irow, jcol]
    return matB


def is_list_ranged(a, List, b):
    """
    Returns true if a<= x <= b or a-x < 0 < b-x

    :param a: the lower bound value (inclusive)
    :param x: the search values
    :param b: the upper bound value (inclusive)

    :returns is_ranged: True/False
    """
    for x in List:
        if not is_float_ranged(a, x, b):
            return False
    return True


def is_float_ranged(a, x, b):
    """
    Returns true if a<= x <= b or a-x < 0 < b-x.

    :param a: the lower bound value (inclusive)
    :param x: the search value
    :param b: the upper bound value (inclusive)

    :returns is_ranged: True/False
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
    B = array(A)
    return ''.join([list_print(B[i, :], tol) + '\n' for i in range(B.shape[0])])


def list_print(listA, tol=1e-8, float_fmt='%-3.2g', zero_fmt='    0'):
    if len(listA) == 0:
        return '[]'

    def _print(a):
        if isinstance(a, string_types):
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

    return '[ '+ ', '.join([_print(a) for a in listA])+ ']'


def augmented_identity(A):
    """
    Creates an Identity Matrix augmented with zeros.
    The location of the extra zeros depends on A.

    .. code-block:: python

      [ 1, 0, 0, 0 ]
      [ 0, 1, 0, 0 ]
      [ 0, 0, 1, 0 ]
    """
    (nx, ny) = A.shape
    I = eye(max(nx, ny), 'float64')
    return I[:nx, :ny]


def solve_tridag(A, D):
    """
    Solves a tridagonal matrix [A]{x}={b} for {x}

    :param A: main diagonal (length=N)
    :param D: off diagonal (length=N-1)
    :returns: x
    """
    # Find the diagonals
    ud = insert(diag(A, 1), 0, 0)  # upper diagonal
    d = diag(A)  # main diagonal
    ld = insert(diag(A, -1), len(d) - 1, 0)  # lower diagonal
    # simplified matrix
    ab = matrix([ud, d, ld])
    return solve_banded((1, 1), ab, D, overwrite_ab=True, overwrite_b=True)


Area = lambda a, b: 0.5 * norm(cross(a, b))

def gauss(n):
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
        return [0.], [2.]
    elif n == 2:
        p = 1. / sqrt(3)
        return [-p, p], [1., 1.]
    elif n == 3:
        p = sqrt(3 / 5.)
        return [-p, 0., p], [5 / 9., 8 / 9., 5 / 9.]
    elif n == 4:
        p1 = (3 - 2. * sqrt(6 / 5)) / 7.
        p2 = (3 + 2. * sqrt(6 / 5)) / 7.
        w1 = (18 + sqrt(30)) / 36.
        w2 = (18 - sqrt(30)) / 36.
        return [-p2, -p1, p1, p2], [w2, w1, w1, w2]
    elif n == 5:
        p1 = 1 / 3. * sqrt(5 - 2 * sqrt(10. / 7.))
        p2 = 1 / 3. * sqrt(5 + 2 * sqrt(10. / 7.))
        w1 = (322 + 13 * sqrt(70)) / 900.
        w2 = (322 - 13 * sqrt(70)) / 900.
        return [-p2, -p1, 0, p1, p2], [w2, w1, 128 / 225., w1, w2]

    raise NotImplementedError('The current implementation only supports up to '
                              '5 quadrature points')
