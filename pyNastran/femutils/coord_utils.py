"""
Defines general coordinate system related functions including:
 - coords = coords_from_vector_1d(v_array)
 - coords = coordinate_system_from_vector_2d_tri(xyz1, xyz2, xyz3)
 - coords = coordinate_system_from_vector_2d_quad(xyz1, xyz2, xyz3, xyz4)
 - coords = coordinate_system_from_vector_2d_tri_theta(xyz1, xyz2, xyz3, theta, dtype=None)

"""
import numpy as np
from .utils import cross2d, perpendicular_vector2d
from .matrix3d import axes_stack, normalize_vector2d


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
    """computes the elemental for a tri element"""
    xyz1 = np.atleast_2d(xyz1)
    xyz2 = np.atleast_2d(xyz2)
    xyz3 = np.atleast_2d(xyz3)

    v21 = xyz2 - xyz1
    v31 = xyz3 - xyz2
    normal = cross2d(v21, v31)
    coords = _coordinate_system_from_vector_2d(v21, normal)
    return coords

def coordinate_system_from_vector_2d_quad(xyz1, xyz2, xyz3, xyz4):
    """computes the elemental for a quad element"""
    xyz1 = np.atleast_2d(xyz1)
    xyz2 = np.atleast_2d(xyz2)
    xyz3 = np.atleast_2d(xyz3)
    xyz4 = np.atleast_2d(xyz4)

    v21 = xyz2 - xyz1
    v31 = xyz3 - xyz1
    v42 = xyz4 - xyz2
    normal = cross2d(v31, v42)
    coords = _coordinate_system_from_vector_2d(v21, normal)
    return coords

def shape4(xy):
    """
    xy : (n, 2)
    """
    xy = np.atleast_2d(xy)
    x = xy[:, 0]
    y = xy[:, 1]
    #xm1 = x - 1
    #ym1 = y - 1
    #xp1 = x - 1
    #yp1 = y - 1
    n1 = (1-x) * (1-y)
    n2 = (1+x) * (1-y)
    n3 = (1+x) * (1+y)
    n4 = (1-x) * (1+y)
    n = np.hstack([n1, n2, n3, n4]) / 4
    return np.atleast_2d(n)

def tet_coord(p1, p2, p3, p4):  # pragma: no cover
    """doesn't compute the elemental for a tet element"""
    ## R vector joins midpoints of edges G1-G2 and G3-G4.
    p12 = (p1 + p2) / 2.
    p34 = (p3 + p4) / 2.
    r = p34 - p12

    ## S vector joins midpoints of edges G1-G3 and G2-G4.
    p13 = (p1 + p3) / 2.
    p24 = (p2 + p4) / 2.
    s = p24 - p13

    ## T vector joins midpoints of edges G1-G4 and G2-G3.
    p14 = (p1 + p4) / 2.
    p23 = (p2 + p3) / 2.
    t = p23 - p14
    origin = p1

    # The element coordinate system is chosen as close as possible to the
    # R, S, and T vectors and points in the same general direction.
    #
    # Mathematically speaking, the coordinate system is computed in such
    # a way that, if the R, S, and T vectors are described in the element
    # coordinate system, a 3x3 positive definite symmetric matrix would
    # be produced.
    #
    # ???
    return r, s, t, origin

def shape8(xyz):  # pragma: no cover
    """
    xyz : (n, 3)
    """
    xyz = np.atleast_2d(xyz)
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    #xm1 = x - 1
    #ym1 = y - 1
    #zm1 = z - 1
    #xp1 = x - 1
    #yp1 = y - 1
    #zp1 = z - 1
    n1 = (1-x) * (1-y) * (1-z)
    n2 = (1+x) * (1-y) * (1-z)
    n3 = (1+x) * (1+y) * (1-z)
    n4 = (1-x) * (1+y) * (1-z)

    n5 = (1-x) * (1-y) * (1+z)
    n6 = (1+x) * (1-y) * (1+z)
    n7 = (1+x) * (1+y) * (1+z)
    n8 = (1-x) * (1+y) * (1+z)
    n = np.hstack([n1, n2, n3, n4, n5, n6, n7, n8]) / 8
    return n

def shape4_to_xyz(p1234, n4):
    """
    p1234 : (n, 4, 3)
        the points of the quad
    shape : (n, 2)
        the values for the she function
    """
    #print('n4 = ', n4)
    p1234 = np.asarray(p1234)
    nquads, npoints, nxyz = p1234.shape
    assert npoints == 4 and nxyz == 3, 'shape=(nquads, npoints, nxyz) = %s' % str(p1234.shape)
    #n = shape4(shape)
    assert n4.shape == (nquads, 4), 'shape=%s' % str(n4.shape)

    # can't I do this more efficiently?
    p = (
        p1234[:, 0] * n4[:, 0] +
        p1234[:, 1] * n4[:, 1] +
        p1234[:, 2] * n4[:, 2] +
        p1234[:, 3] * n4[:, 3]
    )
    #xyz = np.matmul(p1234, n4)
    #xyz = p1234 @ n4
    #xyz = np.tensordot(p1234, n4, axes=0)
    #xyz = np.tensordot(p1234, n4, axes=1)
    #xyz = (p1234 * n4)#.sum(axis=3)
    #np.multiply()
    #print(xyz)
    #assert xyz.shape == (nquads, 3), 'shape=%s' % str(p1234.shape)
    return p

def hexa_coord(xyz1, xyz2, xyz3, xyz4,
               xyz5, xyz6, xyz7, xyz8):  # pragma: no cover
    """doesn't computes the elemental for a hexa element"""
    xyz1 = np.atleast_2d(xyz1)
    xyz2 = np.atleast_2d(xyz2)
    xyz3 = np.atleast_2d(xyz3)
    xyz4 = np.atleast_2d(xyz4)

    xyz5 = np.atleast_2d(xyz5)
    xyz6 = np.atleast_2d(xyz6)
    xyz7 = np.atleast_2d(xyz7)
    xyz8 = np.atleast_2d(xyz8)

    ## R vector joins the centroids of faces G4-G1-G5-G8 and G3-G2-G6-G7.
    c4158 = (xyz4 + xyz1 + xyz5 + xyz8) / 4.
    c3267 = (xyz3 + xyz2 + xyz7 + xyz7) / 4.
    r = c3267 - c4158

    ## S vector joins the centroids of faces G1-G2-G6-G5 and G4-G3-G7-G8.
    c1265 = (xyz1 + xyz2 + xyz6 + xyz5) / 4.
    c4378 = (xyz4 + xyz3 + xyz7 + xyz8) / 4.
    s = c4378 - c1265

    ## T vector joins the centroids of faces G1-G2-G3-G4 and G5-G6-G7-G8.
    c1234 = (xyz1 + xyz2 + xyz3 + xyz4) / 4.
    c5678 = (xyz5 + xyz6 + xyz7 + xyz8) / 4.
    t = c5678 - c1234

    return r, s, t
    #v21 = xyz2 - xyz1
    #v31 = xyz3 - xyz1
    #v42 = xyz4 - xyz2
    #normal = cross2d(v31, v42)
    #coords = _coordinate_system_from_vector_2d(v21, normal)
    #return coords

def _coordinate_system_from_vector_2d(v21, normal):
    r"""
    Helper method::

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

#def coordinate_system_from_vector_2d_tri_theta(xyz1, xyz2, xyz3, theta, dtype=None):
    #"""
         #3
        #/ \
       #/   \
      #/     \
    #1--------2
    #theta is defined relative to 1-2, pivoting about 1
    #"""
    #if dtype is None:
        #dtype = xyz1.dtype
    #dtype = xyz1.dtype
    #coords = coordinate_system_from_vector_2d_tri(xyz1, xyz2, xyz3)
    #rotation = cylindrical_rotation_matrix(theta, dtype=dtype)
    #coords_transformed = dot3d(rotation, coords)
    #return coords_transformed
