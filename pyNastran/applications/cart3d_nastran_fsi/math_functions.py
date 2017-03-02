from __future__ import print_function
import sys
from six import string_types
from numpy import array, cross, allclose
from numpy.linalg import norm, solve

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
    if not a < x:
        if not allclose(x, a):
            return False

    if not x < b:
        if not allclose(x, b):
            return False
    return True


def print_matrix(A):
    msg = ''
    for row in A:
        msg += list_print(row) + '\n'
    return msg

def list_print(listA):
    if len(listA) == 0:
        return '[]'

    msg = '['
    for a in listA:
        if isinstance(a, string_types):
            msg += ' %s,' % (a)
        elif isinstance(a, float):
            msg += ' %-4.2f,' % (a)
        elif isinstance(a, int):
            msg += ' %g,' % (a)
        else:
            try:
                msg += ' %g,' % (a)
            except TypeError:
                print("a = %r" % (a))
                raise

    msg = msg[:-1]
    msg += ' ]'
    return msg


def pierce_plane_vector(p0, p1, p2, pA, pB, pierced_elements):
    """
    http://en.wikipedia.org/wiki/Line-plane_intersection

    [t] = [xa-xb, x1-x0, x2-x0]-1   [xa-x0]
    [u] = [ya-yb, y1-y0, y2-y0]   * [ya-y0]
    [v] = [za-xb, x1-z0, z2-z0]     [za-z0]

    Ax=b          A^-1*b = x
    xa0=mat*tuv   mat*xa0=tuv

    p0,p1,p2 on plane
    a,b are points on the vector
    """
    mat = array([pA-pB, p1-p0, p2-p0]) # numpy array-matrix
    #mat = array([[xa-xb, x1-x0, x2-x0],
    #             [ya-yb, y1-y0, y2-y0],
    #             [za-xb, x1-z0, z2-z0]])
    #xa0 = array([[xa-x0],
    #             [ya-y0],
    #             [za-z0]])
    xa0 = pA-p0  # numpy array
    tuv = solve(mat, xa0)
    return tuv

def distance(x, y):
    """Finds the euclidian/spherical (2D/3D) distance between 2 points"""
    return norm(x - y)


def shepard_SMS_weight(n, nodes):
    """
    a = (R-hi/R*hi)^2
    wi = a/sum(a)

    hi = distance from interpolation point to scatter point
    R  = distance from interpolation point to most distant scatter point

    http://www.ems-i.com/smshelp/Data_Module/Interpolation/Inverse_Distance_Weighted.htm
    Franke & Nielson, 1980
    """
    dists = []
    for nodei in nodes:
        hi = distance(n, nodei)
        dists.append(hi)
    R = min(dists)

    aWeights = []
    for hi in dists:
        a = ((R-hi) / R*hi)**2.
        aWeights.append(a)
    aWeights = array(aWeights)
    weights = aWeights / sum(aWeights)
    return weights


def shepard_weight(n, nodes):
    """
    Finds the weightings based on the distance weighted average method
    In:
        n = node
        nodes = alternate nodes
    Out:
        weights on all nodes (percentage on each node)
    Formula:
        weightI = (1/distanceI)  / sum(1/distanceI) - correct
    http://www.ems-i.com/smshelp/Data_Module/Interpolation/Inverse_Distance_Weighted.htm
    """
    invDists = []
    invDistSum = 0.
    d = 0.
    for nodei in nodes:
        hi = distance(n, nodei)
        d += hi
        if allclose(hi, 0.):
            hi = 0.000001  # make this the dominating node

        invDist = 1.0 / hi
        invDists.append(invDist)
        invDistSum += invDist

    invDists = array(invDists)
    weights = invDists / invDistSum
    return weights


def area_weight(n, n1, n2, n3):
    """
    Finds the weightings based on the barycentric coordinates weighted average method
    http://www.ems-i.com/smshelp/Data_Module/Interpolation/Inverse_Distance_Weighted.htm
    """
    #a = area(n1, n2, n3)
    a1 = area(a, a2, a3)
    a2 = area(a, a1, a3)
    a3 = area(a, a1, a2)
    a = a1 + a2 + a3
    b1 = a1 / a
    b2 = a2 / a
    b3 = a3 / a

    x = b1 * x1 + b2 * x2 + b3 * x3
    #x = (a1 * x1 + a2 * x2 + a3 * x3) / a
    return w1, w2, w3


def get_triangle_weights(n, n1, n2, n3):
    """
    Returns the weighting for each node such that:
        <F1,F2,F3> = F*<w1,w2,w3>
    on the 3 nodes of a triangle
    """
    weights = shepard_weight(n, [n1, n2, n3])
    w1, w2, w3 = weights
    #w1, w2, w3 = areaWeight(n, n1, n2, n3)
    return w1, w2, w3


def Area(a, b):
    """gets the quad/tri area"""
    return 0.5 * norm(cross(a, b))


def AreaNormal(nodes):
    """
    Returns area,unitNormal
    n = Normal = a x b
    Area   = 1/2 * |a x b|
    V = <v1,v2,v3>
    |V| = sqrt(v1^0.5+v2^0.5+v3^0.5) = norm(V)

    Area = 0.5 * |n|
    unitNormal = n/|n|
    """
    (n0, n1, n2) = nodes
    a = n0 - n1
    b = n0 - n2
    vector = cross(a, b)
    length = norm(vector)
    normal = vector / length
    area = 0.5 * length
    if not allclose(norm(normal), 1.):
        print("a = ", a)
        print("b = ", b)
        print("normal = ", normal)
        print("length = ", length)
        sys.exit('check...')
    return area, normal

def triangle_area_centroid_normal(nodes):
    """Returns area,centroid,unitNormal"""
    (area, normal) = AreaNormal(nodes)
    n1, n2, n3 = nodes[0], nodes[1], nodes[2]
    centroid = Centroid(n1, n2, n3)
    return area, centroid, normal


def Normal(a, b):
    """finds the unit normal vector of 2 vectors"""
    vector = cross(a, b)
    length = norm(vector)
    normal = vector / length
    assert allclose(norm(normal), 1.)
    return normal


def Centroid(A, B, C):
    """returns the centroid of a triangle"""
    centroid = (A + B + C) / 3.
    return centroid


def main(): # pragma: no cover
    n1 = array([0., 0., 0.])
    n2 = array([1., 1., 1.])
    n3 = array([1., 0., 0.])
    n4 = array([5., 3., 0.])
    n5 = array([2., 0., 4.])

    n2 = array([0., 1., 0.])
    c1 = Centroid(n1, n2, n3)
    n = Normal(n5, n2)
    print("norm = ", n, norm(n))

    area, centroid, normal = triangle_area_centroid_normal([n1, n2, n3])
    print("area=%s centroid=%s normal=%s" % (area, centroid, normal))


if __name__ == '__main__': # pragma: no cover
    main()
