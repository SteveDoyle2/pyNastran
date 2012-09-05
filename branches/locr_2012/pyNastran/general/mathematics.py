# -*- coding: utf-8 -*-
from __future__ import print_function
#import sys

from numpy import (float32, float64, complex64, complex128, array, cross, 
                   allclose, zeros, matrix, insert, diag)
from numpy.linalg import norm  # , solve

from scipy.linalg import solve_banded
from scipy.interpolate import splrep, splev
from scipy.integrate import quad

from math import sqrt

#if sys.version_info < (3, 0):
#    # "fixes" bug where scipy screws up return code handling
#    import scipy.weave


def integrateLine(x, y):
    """
    Integrates a line of length 1.0
    @param x the independent variable
    @param y the   dependent variable
    """
    if len(set(y)) == 1:
        return y[0]  # (x1-x0 = 1., so yBar*1 = yBar)
    try:
        assert len(x) == len(y), 'x=%s y=%s' % (x, y)
        spline = buildSpline(x, y)
        #y = splev(xi,spline)
        # integrate the area; y=f(x); A=integral(y*dx,x)
        out = quad(splev, 0., 1., args=(spline))
        
        A = out[0]
    except:
        print('spline Error x=%s y=%s' % (x, y))
        raise
    return A


def evaluatePositiveSpline(x, *args):
    spline, minValue = args
    y = splev([x], spline)
    return max(y, minValue)


def buildSpline(x, y):
    """
    builds a cubic spline or 1st order spline if there are less than 3 terms
    @param x the independent variable
    @param y the   dependent variable
    @note a 1st order spline is the same as linear interpolation
    """
    if len(x) < 3:
        # build a linearly interpolated representation
        spline = splrep(x, y, k=1)
    else:
        spline = splrep(x, y)  # build a cubic spline representation
    return spline


def integratePositiveLine(x, y, minValue=0.):
    """
    Integrates a line of length 1.0
    @param x the independent variable
    @param y the   dependent variable
    """
    if len(set(y)) == 1:
        return y[0]  # (x1-x0 = 1., so yBar*1 = yBar)
    try:
        assert len(x) == len(y), 'x=%s y=%s' % (x, y)
        spline = buildSpline(x, y)
        out = quad(evaluatePositiveSpline, 0., 1., args=(spline,
                                                         minValue))  # now integrate the area
        A = out[0]
    except:
        raise RuntimeError('spline Error x=%s y=%s' % (x, y))
    return A


def reduceMatrix(matA, nids):
    """
    takes a list of ids and removes those rows and cols
    """
    nRows = len(nids)
    matB = matrix(zeros((nRows, nRows), 'd'))

    for i, irow in enumerate(nids):
        for j, jcol in enumerate(nids):
            matB[i, j] = matA[irow, jcol]
    return matB


def isListRanged(a, List, b):
    """
    Returns true if a<= x <= b
    or a-x < 0 < b-x
    """
    for x in List:
        if not isFloatRanged(a, x, b):
            return False
    return True


def isFloatRanged(a, x, b):
    """
    Returns true if a<= x <= b
    or a-x < 0 < b-x
    """
    if not a < x:
        if not allclose(x, a):
            return False

    if not x < b:
        if not allclose(x, b):
            return False

    return True


def printAnnotatedMatrix(A, rowNames=None, tol=1e-8):
    """
    takes a list/dictionary and annotates the row number with that value
    indicies go from 0 to N
    """
    if rowNames is None:
        rowNames = [i for i in xrange(A.shape[0])]
    B = array(A)
    msg = ''
    for i in xrange(B.shape[0]):
        rowName = str(rowNames[i])
        msg += '%-2s' % (rowName) + ' ' + ListPrint(B[i, :], tol) + '\n'
    return msg


def printMatrix(A, tol=1e-8):
    B = array(A)
    msg = ''
    for i in xrange(B.shape[0]):
        msg += ListPrint(B[i, :], tol) + '\n'
    #for row in A:
    #    msg += ListPrint(row)+'\n'
    ###
    return msg
###


def ListPrint(listA, tol=1e-8):
    if len(listA) == 0:
        return '[]'
    ###

    msg = '['
    for a in listA:
        if isinstance(a, str):
            msg += ' %s,' % (a)
        elif isinstance(a, float) or isinstance(a, float32) or isinstance(a, float64):
            if abs(a) < tol:
                a = 0.
            msg += ' %-3.2g,' % (a)
        elif isinstance(a, int):
            if abs(a) < tol:
                a = 0.
            msg += ' %3i,' % (a)
        elif isinstance(a, complex64) or isinstance(a, complex128):
            r = '%.4g' % (a.real)
            i = '%+.4gj' % (a.imag)

            if abs(a.real) < tol:
                r = '0'
            if abs(a.imag) < tol:
                i = ''
            temp = '%4s%4s' % (r, i)
            msg += ' %8s,' % (temp.strip())
        else:
            try:
                print("type(a) is not supported...%s" % (type(a)))
                msg += ' %g,' % (a)
            except TypeError:
                print("a = |%s|" % (a))
                raise

    msg = msg[:-1]
    msg += ' ]'
    return msg


def augmentedIdentity(A):
    """
    Creates an Identity Matrix augmented with zeros.
    The location of the extra zeros depends on A.
    """
    (nx, ny) = A.shape
    #(ny,nx) = A.shape
    I = zeros([nx, ny], 'd')

    for i in xrange(nx):
        if i == nx or i == ny:
            break
        I[i][i] = 1.
    return I


def solveTridag(A, D):
    # Find the diagonals
    ud = insert(diag(A, 1), 0, 0)  # upper diagonal
    d = diag(A)  # main diagonal
    ld = insert(diag(A, -1), len(d) - 1, 0)  # lower diagonal

    # simplified matrix
    ab = matrix([ud, d, ld])
    #print "ab = ",ab
    return solve_banded((1, 1), ab, D, overwrite_ab=True, overwrite_b=True)


Area = lambda a, b: 0.5 * norm(cross(a, b))


def AreaNormal(nodes):
    """
    @code
    Returns area,unitNormal
    n = Normal = a x b
    Area   = 1/2 * |a x b|
    V = <v1,v2,v3>
    |V| = sqrt(v1^0.5+v2^0.5+v3^0.5) = norm(V)

    Area = 0.5 * |n|
    unitNormal = n/|n|
    @endcode
    """
    (n0, n1, n2) = nodes
    a = n0 - n1
    b = n0 - n2
    vector = cross(a, b)
    length = norm(vector)
    normal = vector / length
    area = 0.5 * length
    if allclose(norm(normal), 1.) == False:
        print("a = ", a)
        print("b = ", b)
        print("normal = ", normal)
        print("length = ", length)
        raise RuntimeError('check...')
    return (area, normal)


def Triangle_AreaCentroidNormal(nodes): 
    '''
    Returns area,centroid,unitNormal
    
    @param nodes:
      list of three triangle vertices
    '''
    
    (area, normal) = AreaNormal(nodes)
    return (area, centroid(*nodes), normal) #pylint: disable=W0142


def Normal(a, b):
    """Finds the unit normal vector of 2 vectors"""
    vector = cross(a, b)
    normal = vector / norm(vector)
    assert allclose(norm(normal), 1.)
    return normal


def centroid(A, B, C):
    """Returns the centroid of a triangle A, B, C"""
    return (A + B + C) / 3.



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
    
    
    @param n
      Number of quadrature points
    @retval
      Two lists: points and corresponding weights, sorted by points value
      
    @see http://en.wikipedia.org/wiki/Gaussian_quadrature
    """
 
    if n == 1:
        return ([0.], [2.])
    if n == 2:
        p = 1. / sqrt(3)
        return ([-p, p], [1.,1.])
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
