import numpy as np
from pyNastran.bdf.cards.aero.utils import (
    points_elements_from_quad_points, create_axisymmetric_body)

def create_structured_cquad4s(model, pid,
                              p1, p2, p3, p4, nx, ny, nid=1, eid=1):
    """
    Parameters
    ----------
    p1 : (3, ) float ndarray
        leading edge root
    p2 : (3, ) float ndarray
        trailing edge root
    p3 : (3, ) float ndarray
        trailing edge tip
    p4 : (3, ) float ndarray
        leading edge tip
    nx / ny : int
        nx : points in the chordwise direction in percentage of the chord
        ny : points in the spanwise direction in percentage of the span
    nid / eid : int
        node / element id offset

    Returns
    -------
    points (nchord, nspan) float ndarray; might be backwards???
        the points
    elements (nquads, 4) int ndarray
        series of quad elements
        nquads = (nchord-1) * (nspan-1)
    """
    nid0 = nid
    x = np.linspace(0., 1., nx + 1)
    y = np.linspace(0., 1., ny + 1)
    points, elements = points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32')
    for point in points:
        model.add_grid(nid, point)
        nid += 1

    for node_ids in elements + nid0:
        model.add_cquad4(eid, pid, node_ids)
        eid += 1
    return nid, eid

def cone3d(model, pid,
           xs, radius, nx=5, ntheta=10, endpoint=True):
    """
    create a cone by segments in 3d

    xs     = [0., 1., 2.]
    radius = [1., 2., 4.]
    >>> cone(model, pid, xs, radius1, radius2)
    """
    xstation = np.asarray(xstation)
    ystation = np.zeros(xstation.shape)
    zstation = np.zeros(xstation.shape)
    aspect_ratio = 1.0
    xyz_elems = create_axisymmetric_body(
        nx, aspect_ratio,
        xstation, ystation, zstation, radii,
        p1, dy, dz)
    return

