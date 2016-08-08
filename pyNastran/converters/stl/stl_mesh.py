"""
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.479.8237&rep=rep1&type=pdf
http://math.stackexchange.com/questions/4322/check-whether-a-point-is-within-a-3d-triangle
https://en.wikipedia.org/wiki/Barycentric_coordinate_system
http://en.wikipedia.org/wiki/Line-plane_intersection
http://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
"""
from __future__ import print_function
import numpy as np

from pyNastran.converters.stl.stl import read_stl
from scipy.spatial import cKDTree
import scipy.interpolate

def projected_barycentric_coord(p, q, u, v):
    """
    points p, q
    vector u, v
        3
       *v
      /  \    <----p
    q*----*u
    1       2
    u = p2 - p1
    v = p3 - p1
    """
    n = cross(u, v)
    one_over_4_area_squared = 1.0 / dot(n, n)
    w = p - q
    b[2] = dot(cross(u, w), n) * one_over_4_area_squared
    b[1] = dot(cross(w, v), n) * one_over_4_area_squared
    b[0] = 1.0 - b[1] - b[2]
    return b

def project_points_onto_stl(stl, points):
    """
    Parameters
    ----------
    nodes : (n, 3) ndarray floats
        The nodes on the surface.
    elements : (n, 3) ndarray ints
        The elements on the surface.
    """
    nodes = stl.nodes

    elements = stl.elements

    if not hasattr(stl, 'centroids'):
        n1 = elements[:, 0]
        n2 = elements[:, 1]
        n3 = elements[:, 2]
        p1 = nodes[n1, :]
        p2 = nodes[n2, :]
        p3 = nodes[n3, :]
        centroids = (p1 + p2 + p3) / 3.
        stl.centroids = centroids

        tree = cKDTree(centroids, leafsize=16, compact_nodes=True,
                       copy_data=False, balanced_tree=True)
        stl.tree = tree
    #tree = scipy.spatial.KDTree(data, leafsize=10)
    #tree.query_ball_point(x, r, p=2., eps=0)
    #tree.query_ball_tree(other, r, p=2., eps=0)
    #tree.query_pairs(r, p=2., eps=0)
    #tree.sparse_distance_matrix(other, max_distance, p=2.)

    tree = stl.tree
    #d : array of floats
        #The distances to the nearest neighbors.
        #If x has shape tuple+(self.m,), then d has shape tuple+(k,).
        #Missing neighbors are indicated with infinite distances.
    #i : ndarray of ints
        #The locations of the neighbors in self.data.
        #If `x` has shape tuple+(self.m,), then `i` has shape tuple+(k,).
        #Missing neighbors are indicated with self.n.
    dist, i = tree.query(points, k=1, eps=0, p=2,
                          distance_upper_bound=np.inf, n_jobs=1)

    # distance from centroid to point, such that we get the element id directly
    print('i =', i)
    n1 = elements[i, 0]
    n2 = elements[i, 1]
    n3 = elements[i, 2]
    p1 = nodes[n1, :]
    p2 = nodes[n2, :]
    p3 = nodes[n3, :]
    u = p2 - p1
    v = p3 - p1
    #w = points_rotated - p1
    n = np.cross(u, v)
    #n2 = 1 / n**2
    #gamma_a = np.dot(np.cross(u, w), n) / n2
    #gamma_b = np.dot(np.cross(u, w), n) / n2
    try:
        nmag = np.linalg.norm(n, axis=1)
    except ValueError:
        print('n.shape =', n.shape)
        raise
    #area = nmag / 2.
    assert nmag.size == n1.size, 'nmag.size=%s n1.size=%s' % (nmag.size, n1.size)

    print('n1 =', n1)
    print('n2 =', n2)
    print('n3 =', n3)

    p = points
    a = p1
    b = p2
    c = p3
    # http://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
    pc = p - c
    alpha = np.linalg.norm(np.cross(p - b, pc)) / nmag
    beta = np.linalg.norm(np.cross(pc, p - a)) / nmag
    gamma = 1 - alpha - beta
    print('alpha =', alpha)
    print('beta =', beta)
    print('gamma =', gamma)
    #print('a*p =', alpha[:, np.newaxis] * p1)
    p_prime = alpha[:, np.newaxis] * p1 + beta[:, np.newaxis] * p2 + gamma[:, np.newaxis] * p3
    print('p_prime =\n', p_prime)
    #tree.query_ball_point(x, r, p=2., eps=0)
    #tree.query_ball_tree(other, r, p=2., eps=0)
    #tree.query_pairs(r, p=2., eps=0)
    #tree.sparse_distance_matrix(other, max_distance, p=2.)
    return p_prime

def project_line_onto_stl(stl, pa, pb, npoints=11):
    # top down projection
    normal = np.array([0., 0., -1.], dtype='float32')

    #max_z = nodes[:, 2].max()
    #min_z = nodes[:, 2].min()
    # TODO: rotate if want a new normal
    #dz = max_z - min_z
    #dzi = dz / 20.
    #points_rotated = points

    #out_points = project_points_onto_stl(stl, points)
    # TODO: rotate if want a new normal

    p = np.linspace(0., 1., num=npoints, endpoint=True)
    p21 = pb - pa
    ratio = p21 / np.linalg.norm(p21)
    print('p =', p)
    print('ratio =', ratio)
    points = pa + p[:, np.newaxis] * ratio
    print('points =', points)
    out_points = project_points_onto_stl(stl, points)
    return out_points

def project_curve_onto_stl(stl, points, npoints=11):
    # top down projection
    normal = np.array([0., 0., -1.], dtype='float32')

    #max_z = nodes[:, 2].max()
    #min_z = nodes[:, 2].min()
    # TODO: rotate if want a new normal
    #dz = max_z - min_z
    #dzi = dz / 20.
    #points_rotated = points

    #out_points = project_points_onto_stl(stl, points)
    # TODO: rotate if want a new normal

    # create interpolation curve from points
    p2 = points[1:, :]
    p1 = points[:-1, :]
    dx = np.linalg.norm(p2 - p1, axis=1)
    assert dx.size == p1.shape[0]
    t = dx.sum()

    pa = points[0, :]
    dp = points - pa
    dx2 = np.linalg.norm(dp, axis=1)
    t = dx2.sum()

    # http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.interpolate.interp1d.html
    func = scipy.interpolate.interp1d(t, dx2, kind='cubic', axis=-1,
                                      copy=True,
                                      bounds_error=None,
                                      fill_value=np.nan,
                                      assume_sorted=False) # cubic spline
    p = np.linspace(0., t, num=npoints, endpoint=True)
    t2 = func(p)
    dx = func(t2) + pa
    #p21 = pb - pa
    #ratio = p21 / np.linalg.norm(p21)
    #print('p =', p)
    #print('ratio =', ratio)
    points = pa + dx
    print('points =', points)
    out_points = project_points_onto_stl(stl, points)
    return out_points

def main():
    import os
    import pyNastran

    pkg_path = pyNastran.__path__[0]
    stl_filename = os.path.join(pkg_path, 'converters', 'stl', 'sphere.stl')

    stl = read_stl(stl_filename)
    #XYZ Global = (2.0035907914418716, 1.3287668328026303, 2.873731014735773)
    #NodeID = 142; xyz=(1.88823, 1.5, 2.94889)
    #lineNo=2110 annotate_cell_picker()
    #XYZ Global = (1.9419959964242275, 1.141259948469464, 2.869267723165781)
    #NodeID = 141; xyz=(1.93018, 1.02165, 2.85504)
    #lineNo=2110 annotate_cell_picker()
    #XYZ Global = (2.1320656653448338, 1.4367816967143772, 2.83778333777658)
    #NodeID = 137; xyz=(2.25, 1.5, 2.79904)
    # nids = [142, 137, 141]
    # 2.0035907914418716, 1.3287668328026303, 2.873731014735773
    points = np.array([
        [2.0035907914418716, 1.3287668328026303, 2.873731014735773],
        [2.25, 1.5, 2.79904],
        [2.25, 1.5, 2.79903],
    ], dtype='float32')
    pa = points[0, :]
    pb = points[1, :]
    out_points = project_points_onto_stl(stl, points)
    out_points2 = project_line_onto_stl(stl, pa, pb, npoints=11)
    #out_points3 = project_curve_onto_stl(stl, points, npoints=11)

def build():
    from pyNastran.bdf.bdf import BDF, GRID, CQUAD4, PSHELL, PROD, MAT1
    model = BDF(debug=False)
    xyz1 = [0., 0., 0.]
    xyz2 = [0., 1., 0.]
    xyz3 = [1., 1., 0.]
    xyz4 = [1., 0., 0.]
    model.nodes[1] = GRID(nid=1, cp=0, xyz=xyz1)
    model.nodes[2] = GRID(nid=2, cp=0, xyz=xyz2)
    model.nodes[3] = GRID(nid=3, cp=0, xyz=xyz3)
    model.nodes[4] = GRID(nid=4, cp=0, xyz=xyz4)
    model.elements[1] = CQUAD4(eid=1, pid=100, nids=[1,2,3,4])
    model.properties[100] = PSHELL(pid=100, mid1=1000, t=0.1)
    model.materials[1000] = MAT1(mid=1000, E=1e7, G=None, nu=0.3)

if __name__ == '__main__':
    build()
    main()
