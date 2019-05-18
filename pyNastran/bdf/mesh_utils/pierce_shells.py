"""
Defines:
 - pierce_shell_model(bdf_filename, xyz_points, tol=1.0)
"""
from itertools import count
from typing import List, Optional
import numpy as np
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_equivalence import _get_tree


def quad_intersection(orig, direction, v0, v1, v2, v3):
    # type: (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray) -> Optional[np.ndarray]
    """
    Pierces a quad

    Parameters
    ----------
    orig : (3, ) float ndarray
        the point to pierce
    direction : (3, ) float ndarray
        the pierce vector
    v0, v1, v2, v3 : (3, ) float ndarray
        the xyz points of the quad

    Returns
    -------
    p : (3, ) float ndarray or None
        ndarray : the pierce point
        None : failed pierce

    """
    ans = triangle_intersection(orig, direction, v0, v1, v2)
    if ans is not None:
        return ans
    return triangle_intersection(orig, direction, v0, v2, v3)

def triangle_intersection(orig, direction, v0, v1, v2):
    # type: (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray) -> Optional[np.ndarray]
    """
    Pierces a triangle

    Parameters
    ----------
    orig : (3, ) float ndarray
        the point to pierce
    direction : (3, ) float ndarray
        the pierce vector
    v0, v1, v2 : (3, ) float ndarray
        the xyz points of the triangle

    Returns
    -------
    p : (3, ) float ndarray or None
        ndarray : the pierce point
        None : failed pierce

    """
    e1 = v1 - v0
    e2 = v2 - v0
    # Calculate planes normal vector
    pvec = np.cross(direction, e2)
    det = e1.dot(pvec)

    # Ray is parallel to plane
    if abs(det) < 1e-8:
        return None

    inv_det = 1 / det
    tvec = orig - v0
    u = tvec.dot(pvec) * inv_det
    if u < 0 or u > 1:
        return None
    qvec = np.cross(tvec, e1)
    v = direction.dot(qvec) * inv_det
    if v < 0 or u + v > 1:
        return None
    return orig + direction * (e2.dot(qvec) * inv_det)


def pierce_shell_model(bdf_filename, xyz_points, tol=1.0):
    # type: (Union[BDF, str], Any, float) -> List[int], np.ndarray, List[List[int]]
    """
    Pierces a shell model with a <0., 0., 1.> vector.  In other words,
    models are pierced in the xy plane.

    Parameters
    ----------
    bdf_filename : str / BDF()
        the model to run
    xyz_points : (npoints, 3) float ndarray
        the xyz_points to pierce
    tol : float; default=1.0
        the pierce tolerance
        pick a value that is ~3x the max local element edge length

    Returns
    -------
    eids_pierce : List[int]
        int : The element ids that were pierced.
              If multiple elements are pierced, the one with the largest
              pierced z value will be returned.
        None : invalid pierce
    xyz_pierces_max : List[float ndarray, None]
        ndarray : pierce location
        None : invalid pierce
    node_ids : List[int ndarray, None]
        ndarray : pierced element's nodes
        None : invalid pierce

    """
    xyz_points = np.asarray(xyz_points)
    assert xyz_points.shape[1] == 3, xyz_points.shape
    xy_points = xyz_points[:, :2].copy()
    assert xy_points.shape[1] == 2, xy_points.shape

    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename)

    eids = []
    centroids = []
    etypes_to_skip = [
        'CHEXA', 'CPENTA', 'CTETRA', 'CPYRAM',
        'CBUSH', 'CBUSH1D', 'CBUSH2D',
        'CBEAM', 'CROD', 'CONROD',
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    ]
    for eid, elem in model.elements.items():
        if elem.type in ['CQUAD4', 'CTRIA3']:
            centroid = elem.Centroid()
        elif elem.type in etypes_to_skip:
            continue
        else:
            print(elem.type)
            #raise NotImplementedError(elem)
            continue
        eids.append(eid)
        centroids.append(centroid)
    eids = np.array(eids, dtype='int32')
    assert len(eids) > 0, 'eids=%s\n' % eids
    centroids_xy = np.array(centroids, dtype='float64')[:, :2]

    #print('centroids:\n%s\n' % str(centroids_xy))
    assert centroids_xy.shape[1] == 2, centroids_xy.shape


    msg = 'which is required by pierce_shell_model'
    kdt = _get_tree(centroids_xy, msg=msg)
    results = kdt.query_ball_point(xy_points, tol)
    nresults = len(results)
    iresults = np.arange(nresults)
    #print(results)
    #print(iresults)

    direction = np.array([0., 0., 1.])
    ipoints = []
    eids_pierce = []
    xyz_pierces_max = []
    node_ids = []
    for i, xyz_point, eidsi in zip(count(), xyz_points, results):
        model.log.debug('eidsi = [%s]' % ', '.join([str(eidii).rstrip('L') for eidii in eidsi]))
        if not eidsi:
            ipoints.append(i)
            eids_pierce.append(None)
            xyz_pierces_max.append(None)
            node_ids.append(None)
            model.log.warning('skipping %s because it failed tolerancing (tol=%s)' % (xyz_point, tol))
            continue

        model.log.debug('xyz_point = %s' % xyz_point)
        eidsi2 = eids[eidsi]
        eids_newi = []
        piercesi = []
        zpiercesi = []
        for eid in eidsi2:
            elem = model.elements[eid]
            if elem.type == 'CQUAD4':
                v0 = elem.nodes_ref[0].get_position()
                v1 = elem.nodes_ref[1].get_position()
                v2 = elem.nodes_ref[2].get_position()
                v3 = elem.nodes_ref[3].get_position()
                xyz_pierce = quad_intersection(xyz_point, direction, v0, v1, v2, v3)
            elif elem.type == 'CTRIA3':
                v0 = elem.nodes_ref[0].get_position()
                v1 = elem.nodes_ref[1].get_position()
                v2 = elem.nodes_ref[2].get_position()
                xyz_pierce = triangle_intersection(xyz_point, direction, v0, v1, v2)
            else:
                raise NotImplementedError(elem)

            if xyz_pierce is None:
                #model.log.warning('failed to pierce eid=%s' % eid)
                continue

            zpierce = xyz_pierce[2]
            model.log.debug('  eid=%s xyz_pierce=%s zpierce=%s' % (eid, xyz_pierce, zpierce))
            eids_newi.append(eid)
            piercesi.append(xyz_pierce)
            zpiercesi.append(zpierce)

        if len(zpiercesi) == 0:
            ipoints.append(i)
            eids_pierce.append(None)
            xyz_pierces_max.append(None)
            node_ids.append(None)
            model.log.warning('skipping %s because no pierces found (tol=%s)' % (xyz_point, tol))
            continue

        #print('zpiercesi = %s' % zpiercesi)
        zpiercesi_max = max(zpiercesi)
        ipierce_max = zpiercesi.index(zpiercesi_max)
        eid_max = eids_newi[ipierce_max]
        xyz_piercei_max = piercesi[ipierce_max]
        model.log.debug('i=%s xyz_point=%s xyz_piercei_max=%s zpiercesi=%s zmax=%s eid_max=%s' % (
            i, xyz_point, xyz_piercei_max,
            zpiercesi, zpiercesi_max, eid_max))

        ipoints.append(i)
        eids_pierce.append(eid_max)
        xyz_pierces_max.append(xyz_piercei_max)
        #print(model.elements[eid_max])
        node_ids.append(model.elements[eid_max].node_ids)

    xyz_pierces_max = np.array(xyz_pierces_max)
    model.log.debug('ipoints=%s' % ipoints)
    model.log.info('eids_pierce=%s' % eids_pierce)
    model.log.info('xyz_pierces_max:\n%s' % xyz_pierces_max)
    model.log.info('node_ids=%s' % node_ids)

    return eids_pierce, xyz_pierces_max, node_ids
