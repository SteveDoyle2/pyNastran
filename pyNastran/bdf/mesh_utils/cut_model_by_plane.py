from __future__ import print_function

from six import iterkeys
import numpy as np
from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model


def cut_model_by_axes(bdf_filename, view_up, p1, p2, tol,
                      nodal_result, plane_atol=1e-5):
    assert isinstance(tol, float), tol
    p1 = np.asarray(p1)
    p2 = np.asarray(p2)
    view_up = np.asarray(view_up)

    nids, xyz_cid0, edges = _setup(bdf_filename)

    local_points_array, global_points_array, result_array = cut_model(
        nids, xyz_cid0, edges, view_up, p1, p2, tol,
        nodal_result, plane_atol=plane_atol)
    return local_points_array, global_points_array, result_array

def _setup(bdf_filename):
    model = get_bdf_model(bdf_filename, xref=False, log=None, debug=False)
    out = model.get_xyz_in_coord_array(cid=0, fdtype='float64', idtype='int32')
    nid_cp_cd, xyz_cid0, unused_xyz_cp, unused_icd_transform, unused_icp_transform = out
    nids = nid_cp_cd[:, 0]
    #eid_to_edge_map, nid_to_edge_map, edge_to_eid_map = create_maps(model)
    #model = BDF()
    out = model._get_maps(eids=None, map_names=None,
                          consider_0d=False, consider_0d_rigid=False,
                          consider_1d=False, consider_2d=True, consider_3d=False)
    edge_to_eid_map = out['edge_to_eid_map']
    edges = iterkeys(edge_to_eid_map)
    return nids, xyz_cid0, edges

def cut_model_by_coord(bdf_filename, coord, tol,
                       nodal_result, plane_atol=1e-5):
    """
    Cuts a Nastran model with a cutting plane

    Parameters
    ----------
    bdf_filename : str / BDF
        str : the bdf filename
        model : a properly configurated BDF object
    coord : Coord
        the coordinate system to cut the model with
    nodal_result : (nelements, ) float np.ndarray
        the result to cut the model with
    plane_atol : float; default=1e-5
        the tolerance for a line that's located on the y=0 local plane
    """
    assert isinstance(tol, float), tol
    nids, xyz_cid0, edges = _setup(bdf_filename)
    local_points_array, global_points_array, result_array = _cut_model_by_coord(
        nids, xyz_cid0, edges, coord, tol,
        nodal_result, plane_atol=plane_atol)
    return local_points_array, global_points_array, result_array

def cut_model(nids, xyz_cid0, edges, view_up, p1, p2, tol,
              nodal_result, plane_atol=1e-5):
    """
    Parameters
    ----------
    edges : (nedges, 2) int ndarray
        the edges
    view_up : (3,) float ndarray
        the up/z vector
    p1 / p2 : (3,) float ndarray
        the start/end points
    tol : float
        the plane tolerance
    """
    #view_up = camera.GetViewUp()

    #print('p1 =', p1)
    #print('p2 =', p2)
    z = view_up
    x = p2 - p1
    origin = p1
    i = x / np.linalg.norm(x)
    k = z / np.linalg.norm(z)

    # j axis (y direction) is normal to the plane
    j = np.cross(k, i)
    #print("i = ", i)
    #print("j = ", j)
    #print("k = ", k)

    cid = 1
    zaxis = origin + z
    xzplane = origin + x
    coord = CORD2R(cid, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane,
                   comment='')
    local_points_array, global_points_array, result_array = _cut_model_by_coord(
        nids, xyz_cid0, edges, coord, tol, nodal_result,
        plane_atol=plane_atol)
    return local_points_array, global_points_array, result_array

def _merge_bodies(local_points_array, global_points_array, result_array):
    local_points_dict = {}
    global_points_dict = {}
    result_dict = {}
    #return local_points_dict, global_points_dict, result_dict
    return NotImplementedError()

def _cut_model_by_coord(nids, xyz_cid0, edges, coord, tol,
                        nodal_result, plane_atol=1e-5):
    """
    Cuts a Nastran model with a cutting plane

    Parameters
    ----------
    nids : (nnodes, ) int ndarray
        the node ids in the model
    xyz_cid0 : (nnodes, 3) float ndarray
        the node xyzs in the model
    edges : ???
        the edges of the model
    coord : Coord
        the coordinate system to cut the model with
    tol : float
        the tolerance to filter edges (using some large value) to prevent
        excessive computations
    nodal_result : (nelements, ) float np.ndarray
        the result to cut the model with
    plane_atol : float; default=1e-5
        the tolerance for a line that's located on the y=0 local plane
    """
    xyz_cid = coord.transform_node_to_local_array(xyz_cid0)

    # y direction is normal to the plane
    y = xyz_cid[:, 1]
    abs_y = np.abs(y)
    #print('tol =', tol)
    iclose = np.where(abs_y <= tol)
    nids_close = nids[iclose]
    #print('nids_close =', nids_close)

    close_edges = get_close_edges(edges, nids_close)
    close_edges_array = np.array(close_edges)
    #print('close_edges_array:')
    #for edge in close_edges_array:
        #print(edge)
    #print(close_edges_array)
    #aaa

    iclose_edges_array = np.searchsorted(nids, close_edges_array.ravel()).reshape(
        close_edges_array.shape)

    #print('iclose_edges_array:')
    #print(iclose_edges_array)

    local_points_array, global_points_array, result_array = slice_shell_elements(
        xyz_cid0, xyz_cid, iclose_edges_array, nodal_result, plane_atol=plane_atol)

    #print(coord)
    return local_points_array, global_points_array, result_array

def get_close_edges(edges, nids_close):
    """this seems like it could be a lot faster"""
    #n1 = edges[:, 0]
    #n2 = edges[:, 1]
    close_edges = []
    for edge in edges:
        (n1, n2) = edge
        if n1 not in nids_close and n2 not in nids_close:
            continue
        close_edges.append(edge)
    return close_edges

def slice_shell_elements(xyz_cid0, xyz_cid, edges, nodal_result, plane_atol=1e-5):
    """
    Slices the shell elements

    Parameters
    ----------
    xyz_cid0 : (nnodes, 3) float ndarray
        the node xyzs in the model
    xyz_cid : (nnodes, 3) float ndarray
        the node xyzs in the model in the local plane
    edges : ???
        the edges of the model
    nodal_result : (nelements, ) float np.ndarray
        the result to cut the model with
    plane_atol : float; default=1e-5
        the tolerance for a line that's located on the y=0 local plane
    """
    plane_bdf_filename = 'plane.bdf'
    fbdf = open(plane_bdf_filename, 'w')
    fbdf.write('$pyNastran: punch=True\n')
    fbdf.write('MAT1,1,3.0e7,,0.3\n')
    cid = 0
    local_points = []
    global_points = []
    nid_new = 1
    eid_new = 1
    mid = 1
    area = 1.0
    J = 1.0

    result = []
    for edge in edges:
        (nid1, nid2) = edge
        xyz1_local = xyz_cid[nid1]
        xyz2_local = xyz_cid[nid2]
        xyz1_global = xyz_cid0[nid1]
        xyz2_global = xyz_cid0[nid2]
        y1_local = xyz1_local[1]
        y2_local = xyz2_local[1]
        abs_y1_local = np.abs(y1_local)
        abs_y2_local = np.abs(y2_local)
        is_same_sign = np.sign(y1_local) == np.sign(y2_local)
        is_far_from_plane = abs_y1_local > plane_atol and abs_y2_local > plane_atol

        #print('edge=%s' % (edge))
        #print('  xyz1-local=%s xyz2-local=%s' % (xyz1_local, xyz2_local))
        #print('  xyz1-global=%s xyz2-global=%s' % (xyz1_global, xyz2_global))
        #print('  is_same_sign=%s is_far_from_plane=%s' % (is_same_sign, is_far_from_plane))
        if is_far_from_plane and is_same_sign:
            #print('skip y1_local=%.3f y2_local=%.3f plane_atol=%.e' % (y1_local, y2_local, plane_atol))
            continue
        elif np.allclose(y1_local, y2_local, atol=plane_atol):
            #print('  y-sym; nid1=%s nid2=%s edge=%s' % (nid1, nid2, str(edge)))
            #print('     xyz1=%s xyz2=%s' % (xyz1_global, xyz2_global))
            out_grid1 = ['GRID', nid_new, cid, ] + list(xyz1_local)
            out_grid2 = ['GRID', nid_new + 1, cid, ] + list(xyz2_local)
            conrod = ['CONROD', eid_new, nid_new, nid_new + 1, mid, area, J]
            conm2 = ['CONM2', eid_new+1, nid_new, 0, 100.]
            fbdf.write(print_card_8(out_grid1))
            fbdf.write(print_card_8(out_grid2))
            fbdf.write(print_card_8(conrod))
            fbdf.write(print_card_8(conm2))

            local_points.append(xyz1_local)
            local_points.append(xyz2_local)

            global_points.append(xyz1_global)
            global_points.append(xyz2_global)

            x1, y1, z1 = xyz1_local
            x2, y2, z2 = xyz2_local
            x1g, y1g, z1g = xyz1_global
            x2g, y2g, z2g = xyz2_global
            result1 = nodal_result[nid1]
            result2 = nodal_result[nid2]
            result.append([nid1, x1, y1, z1, x1g, y1g, z1g, result1])
            result.append([nid2, x2, y2, z2, x1g, y1g, z1g, result2])

            nid_new += 2
            eid_new += 2

        elif is_same_sign:  # Labs == Lpos
            # same sign, so no crossing
            #print('*edge =', edge)
            #print("  xyz1_global=%s xyz2_global=%s" % (xyz1_global, xyz2_global))
            #print("  xyz1_local =%s xyz2_local =%s" % (xyz1_local, xyz2_local))
            continue
        else:
            #print('edge =', edge)
            # these edges have crossings
            # reworking:
            #  y = m*x + b
            #
            # with the long form:
            #  y = (y2-y1) / (x2-x1) * (x-x1) + y1
            #
            # to get:
            #   y = y2 * (x-x1)/(x2-x1) + y1 * (1 - (x-x1)/(x2-x1))
            #
            # or:
            #   p = (x-x1)/(x2-x1)  # percent
            #   y = y2 * p + y1 * (1 - p)
            #
            # then we sub the entire point for y
            # and sub out x for the y-coordinate
            #
            # then we just crank the formula where we set the value of "x" to 0.0
            # where "x" is the y-coordinate
            percent = (0. - y1_local) / (y2_local - y1_local)
            avg_local = xyz2_local  * percent + xyz1_local  * (1 - percent)
            avg_global = xyz2_global * percent + xyz1_global * (1 - percent)

            #print('  nid1=%s nid2=%s edge=%s' % (nid1, nid2, str(edge)))
            #print('    xyz1_local=%s xyz2_local=%s' % (xyz1_local, xyz2_local))
            #print('    avg_local=%s' % avg_local)
            #print('    avg_global=%s' % avg_global)
            out_grid1 = ['GRID', nid_new, None, ] + list(xyz1_local)
            out_grid2 = ['GRID', nid_new+1, None, ] + list(xyz2_local)
            conrod = ['CONROD', eid_new, nid_new, nid_new + 1, mid, area, J]
            #conm2 = ['CONM2', eid_new, nid_new, 0, 100.]

            fbdf.write(print_card_8(out_grid1))
            fbdf.write(print_card_8(out_grid2))
            fbdf.write(print_card_8(conrod))
            #fbdf.write(print_card_8(conm2))
            local_points.append(avg_local)
            global_points.append(avg_global)

            xl, yl, zl = avg_local
            xg, yg, zg = avg_global
            result1 = nodal_result[nid1]
            result2 = nodal_result[nid2]
            resulti = result2  * percent + result1  * (1 - percent)
            result.append([0, xl, yl, zl, xg, yg, zg, resulti])

            nid_new += 2
            eid_new += 1
        assert len(local_points) == len(result)
        fbdf.write('$------\n')
    fbdf.close()
    local_points_array = np.array(local_points)
    global_points_array = np.array(global_points)
    result_array = np.array(result)
    return local_points_array, global_points_array, result_array


def test_cut_model_by_axes():
    bdf_filename = r'C:\NASA\m4\formats\git\pyNastran\models\bwb\bwb_saero.bdf'
    p1 = [0., 0., 0.]
    p2 = [20., 0., 0.]
    tol = 10.
    view_up = [0., 0., .2]
    theta = np.radians(np.linspace(0., 360., num=10135))
    nodal_result = np.sin(theta)
    unused_local_points_array, unused_global_points_array, result_array = cut_model_by_axes(
        bdf_filename, view_up, p1, p2, tol, nodal_result)

if __name__ == '__main__':
    test_cut_model_by_axes()
