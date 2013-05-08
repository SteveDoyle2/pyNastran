## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
# -*- coding: utf-8 -*-
# pylint: disable=C0103
from __future__ import print_function

import sys
from numpy import (float32, float64, complex64, complex128, array, cross,
                   allclose, zeros, matrix, insert, diag)
from numpy.linalg import norm  # , solve

from scipy.linalg import solve_banded
from scipy.interpolate import splrep, splev
from scipy.integrate import quad

from math import sqrt

if sys.version_info < (3, 0):
    """
    "fixes" bug where scipy screws up return code handling
    for example:
        >>> import sys
        >>> import scipy.sparse
        >>> sys.exit(1)

    The program's return code is 0
    .. note:: Python v3.0+ doesn't have scipy.weave
    """
    import scipy.weave


def integrate_line(x, y):
    """
    Integrates a line of length 1.0
    :x: the independent variable
    :y: the dependent variable
    """
    if len(set(y)) == 1:
        return y[0]  # (x1-x0 = 1., so yBar*1 = yBar)
    try:
        assert len(x) == len(y), 'x=%s y=%s' % (x, y)
        # integrate the area; y=f(x); A=integral(y*dx,x)
        out = quad(splev, 0., 1., args=(build_spline(x, y)))
    except:
        print('spline Error x=%s y=%s' % (x, y))
        raise
    return out[0]


def build_spline(x, y):
    """
    Builds a cubic spline or 1st order spline if there are less than 3 terms
    :param x: the independent variable
    :param y: the dependent variable
    .. note:: a 1st order spline is the same as linear interpolation
    """
    # build a linearly interpolated representation or cubic one
    return splrep(x, y, k=1) if len(x) < 3 else splrep(x, y)


def integrate_positive_line(x, y, minValue=0.):
    """
    Integrates a line of length 1.0
    :param x: the independent variable
    :param y: the dependent variable
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
    matB = matrix(zeros((nRows, nRows), 'd'))

    for i, irow in enumerate(nids):
        for j, jcol in enumerate(nids):
            matB[i, j] = matA[irow, jcol]
    return matB


def is_list_ranged(a, List, b):
    """
    Returns true if a<= x <= b
    or a-x < 0 < b-x
    """
    for x in List:
        if not is_float_ranged(a, x, b):
            return False
    return True


def is_float_ranged(a, x, b):
    """
    Returns true if a<= x <= b
    or a-x < 0 < b-x
    """
    if (not a < x) and (not allclose(x, a)):
        return False

    if (not x < b) and (not allclose(x, b)):
        return False

    return True


def print_annotated_matrix(A, rowNames=None, tol=1e-8):
    """
    takes a list/dictionary and annotates the row number with that value
    indicies go from 0 to N
    """
    if rowNames is None:
        rowNames = [i for i in xrange(A.shape[0])]
    B = array(A)
    return ''.join([ '%-2s' % (str(rowNames[i])) + ' ' + list_print(B[i, :], tol)
                     + '\n' for i in xrange(B.shape[0])])


def print_matrix(A, tol=1e-8):
    B = array(A)
    return ''.join([list_print(B[i, :], tol) + '\n' for i in xrange(B.shape[0])])


def list_print(listA, tol=1e-8):
    if len(listA) == 0:
        return '[]'

    def _print(a):
        if isinstance(a, str):
            return a
        for i, j in ((float, '%-3.2g'), (float32, '%-3.2g'),
                     (float64, '%-3.2g'), (int, '%3i') ):
            if isinstance(a, i):
                return j % (0. if abs(a) < tol else a)

        if isinstance(a, complex) or isinstance(a, complex64) or isinstance(a, complex128):
            return '%4s%4s' % ('0' if abs(a.real) < 1e-8 else '%.4g' % (a.real),
                                '' if abs(a.imag) < 1e-8 else '%+.4gj' % (a.imag))
        try:
            print("list_print: type(a) is not supported... %s" % (type(a)))
            return ' %g' % a
        except TypeError:
            print("a = |%s|" % a)
            raise

    return '[ '+ ', '.join([_print(a) for a in listA])+ ']'


def augmented_identity(A):
    """
    Creates an Identity Matrix augmented with zeros.
    The location of the extra zeros depends on A.
    """
    (nx, ny) = A.shape
    I = zeros([nx, ny], 'd')

    for i in xrange(nx):
        if i == nx or i == ny:
            break
        I[i][i] = 1.
    return I


def solve_tridag(A, D):
    # Find the diagonals
    ud = insert(diag(A, 1), 0, 0)  # upper diagonal
    d = diag(A)  # main diagonal
    ld = insert(diag(A, -1), len(d) - 1, 0)  # lower diagonal
    # simplified matrix
    ab = matrix([ud, d, ld])
    return solve_banded((1, 1), ab, D, overwrite_ab=True, overwrite_b=True)


Area = lambda a, b: 0.5 * norm(cross(a, b))

def centroid_triangle(n1, n2, n3):
    centroid = (n1 + n2 + n3) / 3.
    return centroid

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


    :param n: Number of quadrature points
    @retval Two lists: points and corresponding weights, sorted by points value
    .. seealso:: http://en.wikipedia.org/wiki/Gaussian_quadrature"""
    if n == 1:
        return ([0.], [2.])
    if n == 2:
        p = 1. / sqrt(3)
        return ([-p, p], [1., 1.])
    if n == 3:
        p = sqrt(3 / 5.)
        return ([-p, 0., p], [5 / 9., 8 / 9., 5 / 9.])
    if n == 4:
        p1 = (3 - 2. * sqrt(6 / 5)) / 7.
        p2 = (3 + 2. * sqrt(6 / 5)) / 7.
        w1 = (18 + sqrt(30)) / 36.
        w2 = (18 - sqrt(30)) / 36.
        return ([-p2, -p1, p1, p2], [w2, w1, w1, w2])
    if n == 5:
        p1 = 1 / 3. * sqrt(5 - 2 * sqrt(10. / 7.))
        p2 = 1 / 3. * sqrt(5 + 2 * sqrt(10. / 7.))
        w1 = (322 + 13 * sqrt(70)) / 900.
        w2 = (322 - 13 * sqrt(70)) / 900.
        return ([-p2, -p1, 0, p1, p2], [w2, w1, 128 / 225., w1, w2])

    raise NotImplementedError('The current implementation only supports up to '
                              '5 quadrature points')