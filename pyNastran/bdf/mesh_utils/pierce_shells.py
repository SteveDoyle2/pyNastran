from itertools import count
from six import iteritems
from six.moves import zip
import numpy as np
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_equivalence import _get_tree


def cross(a,b):
    return np.array([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]], 'd')

def dot(a,b):
    return a.dot(b)

def quad_intersection(orig, direction, v0, v1, v2, v3):
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
    if not ans is None:
        return ans
    return triangle_intersection(orig, direction, v0, v2, v3)

def triangle_intersection(orig, direction, v0, v1, v2):
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
    pvec = cross(direction, e2);
    det = dot(e1, pvec);
    # Ray is parallel to plane
    if abs(det) < 1e-8:
        return None

    inv_det = 1 / det
    tvec = orig - v0
    u = dot(tvec, pvec) * inv_det;
    if (u < 0 or u > 1):
        return None
    qvec = cross(tvec, e1)
    v = dot(direction, qvec) * inv_det
    if v < 0 or u + v > 1:
        return None
    return orig + direction*(dot(e2, qvec) * inv_det)

def test_intersect():
    p0 = np.array([0,0,0], 'd')
    p1 = np.array([1,0,0], 'd')
    p2 = np.array([0,1,0], 'd')
    p3 = np.array([1,1,0], 'd')

    v = np.array([0,0,-1], 'd')
    for i in range(10):
        for j in range(10):
            p = np.array([i*.2-.5, j*.2-.5, 1.], 'd')
            print(i, j, p,
                  triangle_intersection(p, v, p0, p1, p2),
                  quad_intersection(p, v, p0, p1, p3, p2))

def pierce_shell_model(bdf_filename, xyz_points, tol=1.0):
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
        The element ids that were pierced.
        If multiple elements are pierced, the one with the largest
        pierced z value will be returned.
    """
    xyz_points = np.asarray(xyz_points)
    assert xyz_points.shape[1] == 3, xyz_points.shape
    xy_points = xyz_points[:, :2]
    assert xy_points.shape[1] == 2, xy_points.shape

    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename)

    eids = []
    centroids = []
    etypes_to_skip = [
        'CHEXA', 'CPENTA', 'CTETRA',
        'CBUSH', 'CBEAM'
    ]
    for eid, elem in iteritems(model.elements):
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

    print('centroids:\n%s\n' % str(centroids_xy))
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
                v0 = elem.nodes[0].get_position()
                v1 = elem.nodes[1].get_position()
                v2 = elem.nodes[2].get_position()
                v3 = elem.nodes[3].get_position()
                xyz_pierce = quad_intersection(xyz_point, direction, v0, v1, v2, v3)
            elif elem.type == 'CTRIA3':
                v0 = elem.nodes[0].get_position()
                v1 = elem.nodes[1].get_position()
                v2 = elem.nodes[2].get_position()
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
        #zpiercesi = []

        if len(zpiercesi) == 0:
            ipoints.append(i)
            eids_pierce.append(None)
            xyz_pierces_max.append(None)
            node_ids.append(None)
            model.log.warning('skipping %s because no pierces found (tol=%s)' % (xyz_point, tol))
            continue

        print('zpiercesi = %s' % zpiercesi)
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
    model.log.info('ipoints=%s' % ipoints)
    model.log.info('eids_pierce=%s' % eids_pierce)
    model.log.info('xyz_pierces_max:\n%s' % xyz_pierces_max)
    model.log.info('node_ids=%s' % node_ids)

    return eids_pierce, xyz_pierces_max, node_ids

def test_pierce_model():
    model = BDF()
    pid = 10
    mid1 = 100
    model = BDF(debug=False)

    # intersects (min)
    model.add_grid(1, xyz=[0., 0., 0.])
    model.add_grid(2, xyz=[1., 0., 0.])
    model.add_grid(3, xyz=[1., 1., 0.])
    model.add_grid(4, xyz=[0., 1., 0.])
    model.add_cquad4(1, pid, [1, 2, 3, 4])

    # intersects (max)
    model.add_grid(5, xyz=[0., 0., 1.])
    model.add_grid(6, xyz=[1., 0., 1.])
    model.add_grid(7, xyz=[1., 1., 1.])
    model.add_grid(8, xyz=[0., 1., 1.])
    model.add_cquad4(2, pid, [5, 6, 7, 8])

    # intersects (mid)
    model.add_grid(9, xyz=[0., 0., 0.5])
    model.add_grid(10, xyz=[1., 0., 0.5])
    model.add_grid(11, xyz=[1., 1., 0.5])
    model.add_grid(12, xyz=[0., 1., 0.5])
    model.add_cquad4(3, pid, [9, 10, 11, 12])

    # doesn't intersect
    model.add_grid(13, xyz=[10., 0., 0.])
    model.add_grid(14, xyz=[11., 0., 0.])
    model.add_grid(15, xyz=[11., 1., 0.])
    model.add_grid(16, xyz=[10., 1., 0.])
    model.add_cquad4(4, pid, [13, 14, 15, 16])

    model.add_pshell(pid, mid1=mid1, t=2.)

    E = 1.0
    G = None
    nu = 0.3
    model.add_mat1(mid1, E, G, nu, rho=1.0)
    model.validate()

    model.cross_reference()

    xyz_points = [
        [0.4, 0.6, 0.],
        [-1., -1, 0.],
    ]
    pierce_shell_model(model, xyz_points)


def test_pierce_model2(bdf_filename, points_filename):
    xyz_points = np.loadtxt(points_filename)
    model = read_bdf(bdf_filename, xref=False)
    model.cross_reference(xref=True, xref_nodes=True, xref_elements=True,
                         xref_nodes_with_elements=True,
                         xref_properties=True,
                         xref_masses=True,
                         xref_materials=True,
                         xref_loads=True,
                         xref_constraints=False,
                         xref_aero=True,
                         xref_sets=True,
                         xref_optimization=True)
    pierce_shell_model(model, xyz_points)


if __name__ == '__main__':
    #test_pierce_model()
    test_pierce_model2(
        r'C:\NASA\Baseline_05_24_17\BDF\Wing.bdf',
        r'C:\NASA\Baseline_05_24_17\BDF\pierce_xy.txt')
