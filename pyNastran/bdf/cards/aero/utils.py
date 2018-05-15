"""
defines:
 - elements = elements_from_quad(nx, ny, dtype='int32')
 - points, elements = points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32')

"""
import numpy as np

def elements_from_quad(nx, ny, dtype='int32'):
    """
    Creates an array of rectilinear mesh of nodes and then
    grabs indexs it to get the elements

    """
    assert nx > 1
    assert ny > 1

    nelements = (nx - 1) * (ny - 1)
    npoints = nx * ny

    # create a matrix with the point counter
    ipoints = np.arange(npoints, dtype=dtype).reshape((nx, ny))

    # move around the CAERO quad and apply ipoints
    elements = np.zeros((nelements, 4), dtype=dtype)
    elements[:, 0] = ipoints[:-1, :-1].ravel()  # (i,  j  )
    elements[:, 1] = ipoints[1:, :-1].ravel()   # (i+1,j  )
    elements[:, 2] = ipoints[1:, 1:].ravel()    # (i+1,j+1)
    elements[:, 3] = ipoints[:-1, 1:].ravel()   # (i,j+1  )
    return elements

def tri_cap(nelements):
    """
    ::

      1
      *    *2
      |   /      *3
      |  /     /
      | /   /
      0  /---------*4

    create a matrix with the point counter

    """
    npoints = nelements
    ipoints = np.arange(npoints, dtype='int32')

    # move around a circle and apply ipoints
    elements = np.zeros((nelements, 3), dtype='int32')

    # the 0 index defines the center point
    elements[:len(ipoints)-1, 1] = ipoints[1:]
    elements[:len(ipoints)-1, 2] = np.hstack([ipoints[2:], 1])
    return elements

def points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32'):
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
    dtype : str; default='int32'
        the type of elements

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

    elements = elements_from_quad(nx, ny, dtype=dtype)
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


def make_monpnt1s_from_cids(model, nids, cids, cid_to_inids,
                            delete_unused_coords=True):
    """
    Creates MONPNT1s, AECOMPs, and SET1s by a series of coordinate systems

    Parameters
    ----------
    model :: BDF()
        the model object
    nids : (nnodes, ) int ndarray
        all the nodes in the model
    cid_to_inids : Dict[cid]=inids
        cid : int
            coord id
        inids : (nnodes_in_cid, ) int ndarray
            the indicies in nids in sorted order
    cids : List[int]
        cids to create
    delete_unused_coords : bool; default=True
        delete coordinate systems from the model that aren't used

    Note
    ----
    Doesn't write duplicate sets
    """
    nnodes_old = 0
    origin_old = np.zeros(3)
    for cid in cids:
        origin = model.coords[cid].origin
        xyz = origin
        inids = cid_to_inids[cid]
        nidsi = nids[inids]
        nnodes = len(nidsi)
        if nnodes == 0:
            del model.coords[cid]
            continue

        # a coordinate system is the same and can be skipped if all the
        # data about is is the same
        #
        # TODO: this is not good enough in the general case, but is probably
        #       ok practically
        if nnodes == nnodes_old and np.allclose(origin_old, origin):
            if delete_unused_coords:
                del model.coords[cid]
            continue

        label = 'xyz %s %s %s' % (xyz[0], xyz[1], xyz[2])
        #aecomp = cid

        ids = nidsi
        model.add_set1(cid, ids, is_skin=False, comment='')
        aecomp_name = 'ae%i' % cid
        list_type = 'SET1'
        lists = [cid]
        model.add_aecomp(aecomp_name, list_type, lists, comment='')

        name = 'c%i' % cid
        axes = '123456'
        model.add_monpnt1(name, label, axes, aecomp_name, [0., 0., 0.], cp=cid, cd=None, comment='')
        nnodes_old = nnodes
        origin_old = origin
