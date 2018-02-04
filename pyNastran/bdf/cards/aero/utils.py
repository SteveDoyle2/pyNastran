"""
defines:
 - elements = elements_from_quad(nx, ny)
 - points, elements = points_elements_from_quad_points(p1, p2, p3, p4, x, y)
"""
import numpy as np

def elements_from_quad(nx, ny):
    """
    Creates an array of rectilinear mesh of nodes and then
    grabs indexs it to get the elements
    """
    assert nx > 1
    assert ny > 1

    nelements = (nx - 1) * (ny - 1)
    npoints = nx * ny

    # create a matrix with the point counter
    ipoints = np.arange(npoints, dtype='int32').reshape((nx, ny))

    # move around the CAERO quad and apply ipoints
    elements = np.zeros((nelements, 4), dtype='int32')
    elements[:, 0] = ipoints[:-1, :-1].ravel()  # (i,  j  )
    elements[:, 1] = ipoints[1:, :-1].ravel()   # (i+1,j  )
    elements[:, 2] = ipoints[1:, 1:].ravel()    # (i+1,j+1)
    elements[:, 3] = ipoints[:-1, 1:].ravel()   # (i,j+1  )
    return elements

def points_elements_from_quad_points(p1, p2, p3, p4, x, y):
    """
    Creates nodes and elements in a structured grid given 4 points.
    Used to make an CAERO1 panel.

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
    x : (nchord, ) float ndarray
        points in the chordwise direction in percentage of the chord
    y : (nspan, ) float ndarray
        points in the spanwise direction in percentage of the span

    Returns
    -------
    points (nchord, nspan) float ndarray; might be backwards???
        the points
    elements (nquads, 4) int ndarray
        series of quad elements
        nquads = (nchord-1) * (nspan-1)
    """
    nx = x.shape[0]
    ny = y.shape[0]

    elements = elements_from_quad(nx, ny)
    npoints = nx * ny

    # shape the vectors so we can multiply them
    x = x.reshape((1, nx))
    y = y.reshape((1, ny))
    p1 = np.asarray(p1).reshape(1, 3)
    p2 = np.asarray(p2).reshape(1, 3)
    p3 = np.asarray(p3).reshape(1, 3)
    p4 = np.asarray(p4).reshape(1, 3)

    # x repeats ny times and varies slowly
    # y repeats nx times and varies quickly
    xv = np.repeat(x, ny, axis=1).reshape(npoints, 1)
    yv = np.repeat(y, nx, axis=0).reshape(npoints, 1)

    # calculate the points a and b xv% along the chord
    a = xv * p2 + (1 - xv) * p1
    b = xv * p3 + (1 - xv) * p4

    # calculate the point yv% along the span
    points = yv * b + (1 - yv) * a
    assert points.shape == (npoints, 3), 'npoints=%s shape=%s' % (npoints, str(points.shape))

    # create a matrix with the point counter
    #ipoints = np.arange(npoints, dtype='int32').reshape((nx, ny))

    return points, elements

