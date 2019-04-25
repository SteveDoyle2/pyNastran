import numpy as np
from pyNastran.bdf.cards.aero.utils import (
    points_elements_from_quad_points, create_axisymmetric_body)

def create_structured_cquad4s(model, pid,
                              p1, p2, p3, p4, nx, ny, nid=1, eid=1, theta_mcid=0.):
    """
    Parameters
    ----------
    p1 / p2 / p3 / p4 : (3, ) float ndarray
        points defining the quad
    nx : int
        points in the p1-p2 direction
    ny : int
        points in the p1-p4 direction
    nid / eid : int
        node / element id offset

    Returns
    -------
    nid : int
        ???
    eid : int
        ???
    """
    nid0 = nid
    x = np.linspace(0., 1., nx + 1)
    y = np.linspace(0., 1., ny + 1)
    points, elements = points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32')
    for point in points:
        model.add_grid(nid, point)
        nid += 1

    for node_ids in elements + nid0:
        model.add_cquad4(eid, pid, node_ids, theta_mcid=theta_mcid)
        eid += 1
    return nid, eid

#def cone3d(model, pid,
           #xs, radius, nx=5, ntheta=10, endpoint=True):
    #"""
    #create a cone by segments in 3d

    #xs     = [0., 1., 2.]
    #radius = [1., 2., 4.]
    #>>> cone(model, pid, xs, radius1, radius2)
    #"""
    #xstation = np.asarray(xstation)
    #ystation = np.zeros(xstation.shape)
    #zstation = np.zeros(xstation.shape)
    #aspect_ratio = 1.0
    #xyz_elems = create_axisymmetric_body(
        #nx, aspect_ratio,
        #xstation, ystation, zstation, radii,
        #p1, dy, dz)
    #return
