import numpy as np
from pyNastran.bdf.cards.aero.utils import (
    points_elements_from_quad_points, create_axisymmetric_body)

def create_structured_chexas(model, pid,
                             x, y, z, nx, ny, nz, eid=1, nid=1):
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
    nnodes = nx * ny * nz
    nelements = (nx - 1) * (ny - 1) * (nz - 1)
    assert nelements > 0, f'nx={nx} ny={ny} nz={nz} nelements={nelements}'
    #points, elements = points_elements_from_cube_points(
        #p1, p2, p3, p4,
        #p5, p6, p7, p8,
        #x, y, z, dtype='int32')
    xv, yv, zv = np.meshgrid(x, y, z)

    for point in zip(xv.ravel(), yv.ravel(), zv.ravel()):
        model.add_grid(nid, point)
        nid += 1

    node_ids = np.arange(nnodes, dtype='int32').reshape(nx, ny, nz)
    #print(node_ids)
    p1 = node_ids[:nx-1, :ny-1, :nz-1].ravel()
    p2 = node_ids[1:,    :ny-1, :nz-1].ravel()
    p3 = node_ids[1:,       1:, :nz-1].ravel()
    p4 = node_ids[:nx-1,    1:, :nz-1].ravel()

    p5 = node_ids[:nx-1, :ny-1, 1:].ravel()
    p6 = node_ids[1:,    :ny-1, 1:].ravel()
    p7 = node_ids[1:,       1:, 1:].ravel()
    p8 = node_ids[:nx-1,    1:, 1:].ravel()

    elements = np.vstack([p1, p2, p3, p4, p5, p6, p7, p8]).T
    #print(elements)
    #print(elements.shape)

    for node_ids in (elements + nid0).tolist():
        #print(node_ids, type(node_ids))
        model.add_chexa(eid, pid, node_ids)
        eid += 1
    return nid, eid

#def points_elements_from_cube_points(p1, p2, p3, p4,
                                     #p5, p6, p7, p8,
                                     #x, y, z, dtype='int32'):
    #"""
    #Creates nodes and elements in a structured grid given 8 points.

    #Parameters
    #----------
    #p1-p8 : (3, ) float ndarray
        #corner point locations
    #x, y, z : (n, ) float ndarray
        #percentage in x, y, and z directions
    #dtype : str; default='int32'
        #the type of elements

    #Returns
    #-------
    #points (nx, ny, nz, 3) float ndarray
        #the points
    #elements (nquads, 8) int ndarray
        #series of hexa elements
        #nhexas = (nx-1) * (ny-1) * (nz-1)
    #"""
    #nx = x.shape[0]
    #ny = y.shape[0]
    #nz = z.shape[0]

    #elements = elements_from_cube(nx, ny, nz, dtype=dtype)
    #npoints = nx * ny

    ## shape the vectors so we can multiply them
    #x = x.reshape((1, nx))
    #y = y.reshape((1, ny))
    #z = y.reshape((1, nz))
    #p1 = np.asarray(p1).reshape(1, 3)
    #p2 = np.asarray(p2).reshape(1, 3)
    #p3 = np.asarray(p3).reshape(1, 3)
    #p4 = np.asarray(p4).reshape(1, 3)

    #p5 = np.asarray(p5).reshape(1, 3)
    #p6 = np.asarray(p6).reshape(1, 3)
    #p7 = np.asarray(p7).reshape(1, 3)
    #p8 = np.asarray(p8).reshape(1, 3)

    ## x repeats ny times and varies slowly
    ## y repeats nx times and varies quickly
    #xv = np.repeat(x, ny, axis=1).reshape(npoints, 1)
    #yv = np.repeat(y, nx, axis=0).reshape(npoints, 1)

    ## calculate the points a and b xv% along the chord
    #a = xv * p2 + (1 - xv) * p1
    #b = xv * p3 + (1 - xv) * p4

    ## calculate the point yv% along the span
    #points = yv * b + (1 - yv) * a
    #assert points.shape == (npoints, 3), 'npoints=%s shape=%s' % (npoints, str(points.shape))

    ## create a matrix with the point counter
    ##ipoints = np.arange(npoints, dtype='int32').reshape((nx, ny))

    #return points, elements

def create_structured_cquad4s(model, pid,
                              p1, p2, p3, p4,
                              nx, ny, nid=1, eid=1, theta_mcid=0.):
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
    points, elements = points_elements_from_quad_points(
        p1, p2, p3, p4, x, y, dtype='int32')
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
