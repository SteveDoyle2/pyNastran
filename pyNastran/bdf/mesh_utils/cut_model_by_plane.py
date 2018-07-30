"""
defines:
 - local_points_array, global_points_array, result_array = cut_edge_model_by_axes(
        bdf_filename, view_up, p1, p2, tol,
        nodal_result, plane_atol=1e-5)
 - local_points_array, global_points_array, result_array = cut_edge_model_by_coord(
        bdf_filename, coord, tol,
        nodal_result, plane_atol=1e-5)
 - slice_edges(xyz_cid0, xyz_cid, edges, nodal_result, plane_atol=1e-5)
"""
from __future__ import print_function

from six import iterkeys, iteritems
import numpy as np
from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model
NUMPY_VERSION = np.lib.NumpyVersion(np.__version__)
if NUMPY_VERSION < '1.13.0':
    raise RuntimeError('numpy version=%s required>=1.13.0' % NUMPY_VERSION)


#def cut_edge_model_by_axes(bdf_filename, view_up, p1, p2, tol,
                           #nodal_result, plane_atol=1e-5):
    #assert isinstance(tol, float), tol
    #p1 = np.asarray(p1)
    #p2 = np.asarray(p2)
    #view_up = np.asarray(view_up)

    #nids, xyz_cid0, edges = _setup_edges(bdf_filename)

    #local_points_array, global_points_array, result_array = _cut_model(
        #nids, xyz_cid0, edges, view_up, p1, p2, tol,
        #nodal_result, plane_atol=plane_atol)
    #return local_points_array, global_points_array, result_array

def _setup_edges(bdf_filename):
    """helper method"""
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

def _setup_faces(bdf_filename):
    """helper method"""
    model = get_bdf_model(bdf_filename, xref=False, log=None, debug=False)
    out = model.get_xyz_in_coord_array(cid=0, fdtype='float64', idtype='int32')
    nid_cp_cd, xyz_cid0, unused_xyz_cp, unused_icd_transform, unused_icp_transform = out
    nids = nid_cp_cd[:, 0]
    #eid_to_edge_map, nid_to_edge_map, edge_to_eid_map = create_maps(model)
    #model = BDF()
    face_eids = []
    faces = []
    shells = set([
        'CTRIA3', 'CTRIAX', 'CTRIA6', 'CTRIAX6',
        'CQUAD4', 'CQUAD', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUADX8',
        'CSHEAR'])
    for eid, elem in iteritems(model.elements):
        if elem.type in shells:
            #if elem.type == 'CQUAD4':
                # split to 2 faces...not done
            #elif elem.type == 'CTRIA3':
            face_eids.append(eid)
            faces.append(elem.node_ids)

    #out = model._get_maps(eids=None, map_names=None,
                          #consider_0d=False, consider_0d_rigid=False,
                          #consider_1d=False, consider_2d=True, consider_3d=False)
    #edge_to_eid_map = out['edge_to_eid_map']
    #faces = iterkeys(edge_to_eid_map)
    return nids, xyz_cid0, faces, face_eids

def cut_edge_model_by_coord(bdf_filename, coord, tol,
                            nodal_result, plane_atol=1e-5, csv_filename=None):
    """
    Cuts a Nastran model with a cutting plane

    Parameters
    ----------
    bdf_filename : str / BDF
        str : the bdf filename
        model : a properly configurated BDF object
    coord : Coord
        the coordinate system to cut the model with
    tol : float
        the tolerance to filter edges (using some large value) to prevent
        excessive computations
    nodal_result : (nelements, ) float np.ndarray
        the result to cut the model with
    plane_atol : float; default=1e-5
        the tolerance for a line that's located on the y=0 local plane
    csv_filename : str; default=None
        None : don't write a csv
        str : write a csv
    """
    assert isinstance(tol, float), tol
    nids, xyz_cid0, edges = _setup_edges(bdf_filename)
    local_points_array, global_points_array, result_array = _cut_edge_model_by_coord(
        nids, xyz_cid0, edges, coord, tol,
        nodal_result, plane_atol=plane_atol)
    if csv_filename:
        export_edge_cut(csv_filename, result_array)
    return local_points_array, global_points_array, result_array

def cut_face_model_by_coord(bdf_filename, coord, tol,
                            nodal_result, plane_atol=1e-5, csv_filename=None):
    """
    Cuts a Nastran model with a cutting plane

    Parameters
    ----------
    bdf_filename : str / BDF
        str : the bdf filename
        model : a properly configurated BDF object
    coord : Coord
        the coordinate system to cut the model with
    tol : float
        the tolerance to filter faces (using some large value) to prevent
        excessive computations
    nodal_result : (nelements, ) float np.ndarray
        the result to cut the model with
    plane_atol : float; default=1e-5
        the tolerance for a line that's located on the y=0 local plane
    csv_filename : str; default=None
        None : don't write a csv
        str : write a csv
    """
    assert isinstance(tol, float), tol
    nids, xyz_cid0, faces, face_eids = _setup_faces(bdf_filename)
    unique_geometry_array, unique_results_array = _cut_face_model_by_coord(
        nids, xyz_cid0, faces, face_eids, coord, tol,
        nodal_result, plane_atol=plane_atol)
    if csv_filename:
        export_face_cut(csv_filename, unique_geometry_array, unique_results_array)
    print('unique_geometry_array, unique_results_array =', unique_geometry_array, unique_results_array)
    return unique_geometry_array, unique_results_array

def export_edge_cut(csv_filename, result_array):
    #nints = geometry_array.shape[1]
    #nfloats = results_array.shape[1]
    #max_int = geometry_array.max()
    #len_max_int = len(str(max_int))
    #fmt = ('%%%ii,' % (len_max_int)) * nints + '%19.18e,' * nfloats
    #fmt = fmt.rstrip(',')
    #print(fmt)
    #X = np.concatenate((geometry_array, results_array), axis=1)
    header = 'x, y, z, Cp'
    np.savetxt(csv_filename, result_array, delimiter=',', newline='\n', header=header,
               footer='', comments='# ', encoding=None)

def export_face_cut(csv_filename, geometry_array, results_array):
    nints = geometry_array.shape[1]
    nfloats = results_array.shape[1]
    max_int = geometry_array.max()
    len_max_int = len(str(max_int))
    fmt = ('%%%ii,' % (len_max_int)) * nints + '%19.18e,' * nfloats
    fmt = fmt.rstrip(',')
    X = np.concatenate((geometry_array, results_array), axis=1)
    header = 'eid, nid1, nid2, x, y, z, Cp'
    np.savetxt(csv_filename, X, fmt=fmt, newline='\n', header=header,
               footer='', comments='# ', encoding=None)

def _determine_cord2r(origin, zaxis, xzplane):
    k = zaxis / np.linalg.norm(zaxis)
    iprime = xzplane / np.linalg.norm(xzplane)
    j = np.cross(k, iprime)
    j /= np.linalg.norm(j)
    i = np.cross(j, k)
    return i, k, origin, zaxis, xzplane

def _project_z_axis(p1, p2, z_global):
    x = p2 - p1
    iprime = x / np.linalg.norm(x)
    k = z_global / np.linalg.norm(z_global)
    j = np.cross(k, iprime)
    jhat = j / np.linalg.norm(j)
    i = np.cross(jhat, k)

    origin = p1
    zaxis = p1 + k
    xzplane = p1 + i
    return i, k, origin, zaxis, xzplane


def _p1_p2_zaxis_to_cord2r(model, p1, p2, zaxis, method='Z-Axis Projection',
                           cid_p1=0, cid_p2=0, cid_zaxis=0):
    """Creates the coordinate system that will define the cutting plane"""
    p1 = np.asarray(p1)
    p2 = np.asarray(p2)
    zaxis = np.asarray(zaxis)
    print("coord:")
    print('  p1 =', p1)
    print('  p2 =', p2)
    print('  zaxis =', zaxis)

    xyz1 = model.coords[cid_p1].transform_node_to_global(p1)
    xyz2 = model.coords[cid_p2].transform_node_to_global(p2)
    z_global = model.coords[cid_zaxis].transform_node_to_global(zaxis)
    if method == 'CORD2R':
        i, k, origin, zaxis, xzplane = _determine_cord2r(xyz1, xyz2, z_global)
    elif method == 'Z-Axis Projection':
        i, k, origin, zaxis, xzplane = _project_z_axis(xyz1, xyz2, z_global)
    else:
        raise NotImplementedError("method=%r; valid_methods=['CORD2R', 'Z-Axis Projection']")
    return xyz1, xyz2, z_global, i, k, origin, zaxis, xzplane

def _cut_model(nids, xyz_cp, edges, view_up, p1, p2, tol,
               nodal_result, plane_atol=1e-5):
    """
    Helper method for cut_edge_model_by_axes

    Parameters
    ----------
    edges : (nedges, 2) int ndarray
        the edges
    view_up : (3,) float ndarray
        the up/z vector
    p1 / p2 : (3,) float ndarray
        the start/end points
    tol : float
        the tolerance to filter edges (using some large value) to prevent
        excessive computations
    plane_atol : float; default=1e-5
        the tolerance for a line that's located on the y=0 local plane
    """
    #view_up = camera.GetViewUp()

    #print('p1 =', p1)
    #print('p2 =', p2)
    z = view_up
    x = p2 - p1
    origin = p1
    #i = x / np.linalg.norm(x)
    #k = z / np.linalg.norm(z)

    # j axis (y direction) is normal to the plane
    #j = np.cross(k, i)
    #print("i = ", i)
    #print("j = ", j)
    #print("k = ", k)

    cid = 1
    zaxis = origin + z
    xzplane = origin + x
    coord = CORD2R(cid, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane,
                   comment='')
    local_points_array, global_points_array, result_array = _cut_edge_model_by_coord(
        nids, xyz_cp, edges, coord, tol, nodal_result,
        plane_atol=plane_atol)
    return local_points_array, global_points_array, result_array

#def _merge_bodies(local_points_array, global_points_array, result_array):
    #local_points_dict = {}
    #global_points_dict = {}
    #result_dict = {}
    #return local_points_dict, global_points_dict, result_dict
    #return NotImplementedError()

def _cut_edge_model_by_coord(nids, xyz_cid0, edges, coord, tol,
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
    #print('nids_close =', nids_close.tolist())

    # speedup
    close_edges = get_close_edges(edges, nids_close)
    close_edges_array = np.array(close_edges)

    #print('close_edges_array:')
    #for edge in close_edges_array:
        #print(edge)
    #print(close_edges_array)

    iclose_edges_array = np.searchsorted(nids, close_edges_array.ravel()).reshape(
        close_edges_array.shape)

    #print('iclose_edges_array:')
    #print(iclose_edges_array)

    local_points_array, global_points_array, result_array = slice_edges(
        xyz_cid0, xyz_cid, iclose_edges_array, nodal_result, plane_atol=plane_atol)

    #print(coord)
    return local_points_array, global_points_array, result_array

def _cut_face_model_by_coord(nids, xyz_cid0, faces, face_eids, coord, tol,
                             nodal_result, plane_atol=1e-5):
    """
    Cuts a Nastran model with a cutting plane

    Parameters
    ----------
    nids : (nnodes, ) int ndarray
        the node ids in the model
    xyz_cid0 : (nnodes, 3) float ndarray
        the node xyzs in the model
    faces : ???
        the faces of the model
    face_eids : List[int]
        the parent element
    coord : Coord
        the coordinate system to cut the model with
    tol : float
        the tolerance to filter faces (using some large value) to prevent
        excessive computations
    nodal_result : (nelements, ) float np.ndarray
        the result to cut the model with
    plane_atol : float; default=1e-5
        the tolerance for a line that's located on the y=0 local plane
    """
    xyz_cid = coord.transform_node_to_local_array(xyz_cid0)
    #face_eids = np.asarray(face_eids)

    # y direction is normal to the plane
    y = xyz_cid[:, 1]
    abs_y = np.abs(y)
    #print('tol =', tol)
    iclose = np.where(abs_y <= tol)
    nids_close = nids[iclose]
    #print('nids_close =', nids_close.tolist())

    close_faces, close_face_eids = get_close_faces(faces, face_eids, nids_close)
    close_faces_array = np.array(close_faces)
    close_face_eids_array = np.array(close_face_eids)

    #print('close_faces_array:')
    #for edge in close_faces_array:
        #print(face)
    #print(close_faces_array)

    assert np.array_equal(nids, np.unique(nids)), 'not sorted or unique'
    iclose_faces_array = np.searchsorted(nids, close_faces_array.ravel()).reshape(
        close_faces_array.shape)

    #print('iclose_edges_array:')
    #print(iclose_edges_array)

    unique_geometry_array, unique_results_array = slice_faces(
        nids, xyz_cid0, xyz_cid, iclose_faces_array, close_face_eids_array,
        nodal_result, plane_atol=plane_atol)

    #print(coord)
    return unique_geometry_array, unique_results_array

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

def get_close_faces(faces, face_eids, unused_nids_close):
    """this seems like it could be a lot faster"""
    # who cares?
    return faces, face_eids

    #n1 = edges[:, 0]
    #n2 = edges[:, 1]
    #close_faces = []
    #close_face_eids = []
    #for eid, face in zip(face_eids, faces):
        #if not all((nid in nids_close for nid in face)):
            #continue
        #close_faces.append(face)
        #close_face_eids.append(eid)
    #return close_faces, close_face_eids

def slice_edges(xyz_cid0, xyz_cid, edges, nodal_result, plane_atol=1e-5):
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
    plane_bdf_filename = 'plane_edge.bdf'
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
        (inid1, inid2) = edge
        xyz1_local = xyz_cid[inid1]
        xyz2_local = xyz_cid[inid2]
        xyz1_global = xyz_cid0[inid1]
        xyz2_global = xyz_cid0[inid2]
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
            #print('skip y1_local=%.3f y2_local=%.3f plane_atol=%.e' % (
                #y1_local, y2_local, plane_atol))
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
            result1 = nodal_result[inid1]
            result2 = nodal_result[inid2]
            result.append([inid1, x1, y1, z1, x1g, y1g, z1g, result1])
            result.append([inid2, x2, y2, z2, x2g, y2g, z2g, result2])

            nid_new += 2
            eid_new += 2

        elif is_same_sign:  # Labs == Lpos
            # same sign, so no crossing
            #print('*edge =', edge)
            #print("  xyz1_global=%s xyz2_global=%s" % (xyz1_global, xyz2_global))
            #print("  xyz1_local =%s xyz2_local =%s" % (xyz1_local, xyz2_local))
            continue
        else:
            eid_new, nid_new = _interpolate_bar_to_node(
                eid_new, nid_new, mid, area, J, fbdf,
                inid1, inid2,
                xyz1_local, xyz2_local,
                xyz1_global, xyz2_global,
                nodal_result,
                local_points, global_points, result)
        assert len(local_points) == len(result)
        fbdf.write('$------\n')
    fbdf.close()
    local_points_array = np.array(local_points)
    global_points_array = np.array(global_points)
    result_array = np.array(result)
    return local_points_array, global_points_array, result_array

def _interpolate_bar_to_node(eid_new, nid_new, mid, area, J, fbdf,
                             inid1, inid2,
                             xyz1_local, xyz2_local,
                             xyz1_global, xyz2_global,
                             nodal_result,
                             local_points, global_points, result):
    """
    These edges have crossings.  We rework:
     y = m*x + b

    into the long form:
     y = (y2-y1) / (x2-x1) * (x-x1) + y1

    to get:
      y = y2 * (x-x1)/(x2-x1) + y1 * (1 - (x-x1)/(x2-x1))

    or:
      p = (x-x1)/(x2-x1)  # percent
      y = y2 * p + y1 * (1 - p)

    Then we sub the y for the point (3 floats) and sub out x for the
    y-coordinate:
        percent = (y - y1_local) / (y2_local - y1_local)
        avg_xyz = xyz2  * percent + xyz1  * (1 - percent)

    Then we just crank the formula where we set the value of "y" to 0.0:
        percent = (0. - y1_local) / (y2_local - y1_local)
    """
    #print('edge =', edge)
    y1_local = xyz1_local[1]
    y2_local = xyz2_local[1]
    percent = (0. - y1_local) / (y2_local - y1_local)
    avg_local = xyz2_local  * percent + xyz1_local  * (1 - percent)
    avg_global = xyz2_global * percent + xyz1_global * (1 - percent)

    #print('  nid1=%s nid2=%s edge=%s' % (inid1, inid2, str(edge)))
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
    result1 = nodal_result[inid1]
    result2 = nodal_result[inid2]
    resulti = result2  * percent + result1  * (1 - percent)
    result.append([0, xl, yl, zl, xg, yg, zg, resulti])

    nid_new += 2
    eid_new += 1
    return eid_new, nid_new

def slice_faces(nodes, xyz_cid0, xyz_cid, faces, face_eids, nodal_result, plane_atol=1e-5):
    """
    Slices the shell elements

    Parameters
    ----------
    xyz_cid0 : (nnodes, 3) float ndarray
        the node xyzs in the model
    xyz_cid : (nnodes, 3) float ndarray
        the node xyzs in the model in the local plane
    faces : ???
        the faces of the model
    nodal_result : (nelements, ) float np.ndarray
        the result to cut the model with
    plane_atol : float; default=1e-5
        the tolerance for a line that's located on the y=0 local plane
    """
    plane_bdf_filename = 'plane_face.bdf'
    fbdf = open(plane_bdf_filename, 'w')
    fbdf.write('$pyNastran: punch=True\n')
    fbdf.write('MAT1,1,3.0e7,,0.3\n')
    fbdf.write('MAT1,2,3.0e7,,0.3\n')
    #cid = 0
    local_points = []
    global_points = []
    nid_new = 1
    eid_new = 1
    mid = 1
    area = 1.0
    J = 1.0

    result = []
    geometry = []
    tri_faces = faces
    for eid, face in zip(face_eids, tri_faces):
        if len(face) == 4:
            print('skipping face=%s' % face)
            continue
        (inid1, inid2, inid3) = face
        xyz1_local = xyz_cid[inid1]
        xyz2_local = xyz_cid[inid2]
        xyz3_local = xyz_cid[inid3]
        xyz1_global = xyz_cid0[inid1]
        xyz2_global = xyz_cid0[inid2]
        xyz3_global = xyz_cid0[inid3]
        y1_local = xyz1_local[1]
        y2_local = xyz2_local[1]
        y3_local = xyz3_local[1]
        y_local = xyz_cid0[face, 1]
        abs_y_local = np.abs(y_local)
        #abs_y1_local = np.abs(y1_local)
        #abs_y2_local = np.abs(y2_local)
        #abs_y3_local = np.abs(y3_local)
        is_same_sign = np.sign(y1_local) == np.sign(y2_local) == np.sign(y3_local)
        is_far_from_plane = all(np.greater(abs_y_local, plane_atol))
        #print("  y_local = %s" % y_local)
        #print('  xyz1-local=%s xyz2-local=%s' % (xyz1_local, xyz2_local))
        #print('  xyz1-global=%s xyz2-global=%s' % (xyz1_global, xyz2_global))
        #print('  is_same_sign=%s is_far_from_plane=%s' % (is_same_sign, is_far_from_plane))
        if is_far_from_plane and is_same_sign:
            print('  far-face=%s' % (face))
            #print('skip y1_local=%.3f y2_local=%.3f plane_atol=%.e' % (
                #y1_local, y2_local, plane_atol))
            continue
        elif np.allclose(y_local, y1_local, atol=plane_atol): # dot
            continue
            #print('  on_edge y12-face=%s' % (face))
            #print('  xyz1-local=%s' % (xyz1_local,))
            #print('  xyz2-local=%s' % (xyz2_local,))
            #print('  xyz3-local=%s' % (xyz3_local,))
            #eid_new, nid_new = _face_on_edge(eid_new, nid_new, face,
                                             #xyz1_local, xyz2_local,
                                             #cid, mid, area, J, bdf, geometry, result)
        #elif np.allclose(y2_local, y3_local, atol=plane_atol):
        #elif np.allclose(y1_local, y3_local, atol=plane_atol):

        elif is_same_sign:  # Labs == Lpos
            print('  same sign-face=%s' % (face))
            # same sign, so no crossing
            #print('*edge =', edge)
            #print("  xyz1_global=%s xyz2_global=%s" % (xyz1_global, xyz2_global))
            #print("  xyz1_local =%s xyz2_local =%s" % (xyz1_local, xyz2_local))
            continue
        else:
            # a crossing
            print('  intersection-eid=%s face=%s' % (eid, face))
            #print('  is_same_sign=%s is_far_from_plane=%s' % (is_same_sign, is_far_from_plane))
            eid_new, nid_new = _interpolate_face_to_bar(
                nodes, eid, eid_new, nid_new, mid, area, J, fbdf,
                inid1, inid2, inid3,
                xyz1_local, xyz2_local, xyz3_local,
                xyz1_global, xyz2_global, xyz3_global,
                nodal_result,
                local_points, global_points,
                geometry, result, plane_atol)

        if len(local_points) != len(result):
            msg = 'lengths are not equal; local_points=%s result=%s' % (
                len(local_points), len(result))
            raise RuntimeError(msg)
        fbdf.write('$------\n')
        print('----------------------')
    fbdf.close()
    if len(geometry) == 0:
        print('bad')
        return None, None
    unused_local_points_array = np.array(local_points)
    unused_global_points_array = np.array(global_points)
    results_array = np.array(result)

    #for g in geometry:
        #print(g)
    geometry_array = np.array(geometry, dtype='int32')
    unique_geometry_array, unique_results_array = _unique_face_rows(
        geometry_array, results_array, nodes)
    return unique_geometry_array, unique_results_array

def _unique_face_rows(geometry_array, results_array, nodes):
    print(geometry_array)
    #iedges = geometry_array[:, 1:]
    #geometry_array[:, 1:] = nodes[iedges.flatten()].reshape(iedges.shape)
    #aaa
    #myrow = None
    for irow, row in enumerate(geometry_array):
        #if row[0] == 11029:
            #myrow = irow
            #print('ids=%s' % row[1:])
            #print('nids=%s' % nodes[row[1:]])
        out = nodes[row[1:]]
        geometry_array[irow, 1:] = out

    print('geometry_array:')
    print(geometry_array)
    unused_edges = geometry_array[:, 1:]
    #print(','.join([str(val) for val in np.unique(edges)]))
    #print('geom =', geometry_array[myrow, :])

    unused_unique_edges, unique_edge_index = np.unique(
        geometry_array[:, 1:],
        return_index=True, return_inverse=False, return_counts=False,
        axis=0)
    unique_geometry_array = geometry_array[unique_edge_index, :]
    unique_results_array = results_array[unique_edge_index, :]
    #for i in unique_edge_index:
        #(eid, nid1, nid2) = geometry_array[i, :]
        #if nid1 == 192 or nid2 == 192:
            #print('i=%s eid=%s nid1=%s nid2=%s' % (i, eid, nid1, nid2))
        #print(geometry_array[i, :])

    unique_geometry_array2, unique_results_array2 = _connect_face_rows(
        unique_geometry_array, unique_results_array)
    return unique_geometry_array2, unique_results_array2

def _connect_face_rows(geometry_array, results_array):
    return geometry_array, results_array

    # TODO: need to handle the dot case first
    #edges = geometry_array[:, 1:]
    #nodes, counts = np.unique(edges, return_counts=True)
    #print('nodes =', nodes)
    #print('counts =', counts)
    #ibad = np.where(counts != 2)
    #print('unique = ', nodes[ibad])
    #print(','.join([str(val) for val in nodes[ibad]]))
    #return geometry_array2, results_array2

def _face_on_edge(eid, eid_new, nid_new, mid, area, J, fbdf,
                  inid1, inid2, unused_inid3,
                  xyz1_local, xyz2_local, unused_xyz3_local,
                  xyz1_global, xyz2_global, unused_xyz3_global,
                  nodal_result,
                  local_points, global_points,
                  geometry, result, unused_plane_atol):
    #raise NotImplementedError('on edge-y1')
    #print('  y-sym; nid1=%s nid2=%s edge=%s' % (nid1, nid2, str(edge)))
    #print('     xyz1=%s xyz2=%s' % (xyz1_global, xyz2_global))

    cid = 0

    # TODO: hardcoded and wrong
    nid_a = inid1
    nid_b = inid2

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

    unused_x1, unused_y1, unused_z1 = xyz1_local
    unused_x2, unused_y2, unused_z2 = xyz2_local
    x1g, y1g, z1g = xyz1_global
    x2g, y2g, z2g = xyz2_global
    result1 = nodal_result[inid1]
    result2 = nodal_result[inid2]
    #nid1, x1, y1, z1,
    #nid2, x2, y2, z2,

    geometry.append([eid, nid_new, nid_a, nid_b])
    result.append([x1g, y1g, z1g, result1])
    result.append([x2g, y2g, z2g, result2])

    eid_new += 2
    nid_new += 2
    return eid_new, nid_new

def _interpolate_face_to_bar(nodes, eid, eid_new, nid_new, mid, area, J, fbdf,
                             inid1, inid2, inid3,
                             xyz1_local, xyz2_local, xyz3_local,
                             xyz1_global, xyz2_global, xyz3_global,
                             nodal_result,
                             local_points, global_points,
                             geometry, result, plane_atol):
    """
    These edges have crossings.  We rework:
     y = m*x + b

    into the long form:
     y = (y2-y1) / (x2-x1) * (x-x1) + y1

    to get:
      y = y2 * (x-x1)/(x2-x1) + y1 * (1 - (x-x1)/(x2-x1))

    or:
      p = (x-x1)/(x2-x1)  # percent
      y = y2 * p + y1 * (1 - p)

    Then we sub the y for the point (3 floats) and sub out x for the
    y-coordinate:
        percent = (y - y1_local) / (y2_local - y1_local)
        avg_xyz = xyz2  * percent + xyz1  * (1 - percent)

    Then we just crank the formula where we set the value of "y" to 0.0:
        percent = (0. - y1_local) / (y2_local - y1_local)

    That's how you do 1 edge, so we do this 3 times.  One of the edges
    won't be a crossing (the percent is not between 0 and 1.), but 2
    edges are.  Thus, two points create a line.

    We also need to handle the dot case.  We're using a triangle
    (nodes 1, 2, and 3), so we have 3 vectors:
        e0 = e12 = p2 - p1
        e1 = e13 = p3 - p1
        e2 = e23 = p3 - p2

    As metioned previously, only two vectors are used (e.g., e12 and e13).
    When combined with the percentage, we find that for a dot, using e12
    and e13, node 1 must be a source (both vectors originate from node 1).
    Thus the percentages for e12=0. and e13=0.  Similarly, node 3 is a
    sink (both vectors end at node 3) and node 2 is a corner/mixed (one
    vector ends at node 2).  In summary:

    Node Combination Percentages for Dot
    ==== =========== ===================
    1    e12, e13    0., 0.
    2    e12, e23    1., 0.
    3    e13, e23    1., 1.
    """
    #print('edge =', edge)
    if eid == 11029:
        print('eid=%s inid1=%s, inid2=%s, inid3=%s' % (eid, inid1, inid2, inid3))
        print('nid1=%s, nid2=%s, nid3=%s' % (nodes[inid1], nodes[inid2], nodes[inid3]))

    edgesi = (
        # (nid_index, xyz in local frame, xyz in global frame
        ((inid1, xyz1_local, xyz1_global), (inid2, xyz2_local, xyz2_global)), # edge 1-2
        ((inid2, xyz2_local, xyz2_global), (inid3, xyz3_local, xyz3_global)), # edge 2-3
        ((inid1, xyz1_local, xyz1_global), (inid3, xyz3_local, xyz3_global)), # edge 1-3
    )
    nid_a_prime = nid_new
    nid_b_prime = nid_new + 1

    #projected_points = []
    #lengths = []

    # we need to prevent dots
    msg = ''
    results_temp = []
    geometry_temp = []
    i_values = []
    percent_values = []
    local_points_temp = []
    global_points_temp = []
    for i, (edge1, edge2) in enumerate(edgesi):
        (inid_a, p1_local, p1_global) = edge1
        (inid_b, p2_local, p2_global) = edge2
        print('  inid_a=%s, p1_local=%s, p1_global=%s' % (inid_a, p1_local, p1_global))
        print('  inid_b=%s, p2_local=%s, p2_global=%s' % (inid_b, p2_local, p2_global))
        py1_local = p1_local[1]
        py2_local = p2_local[1]
        #length = np.linalg.norm(p2_global - p1_global)
        #lengths.append(length)
        dy = py2_local - py1_local
        if np.allclose(dy, 0.0, atol=plane_atol):
            continue

        #dy21 = y2_local - y1_local
        #dy31 = y3_local - y1_local
        #dy23 = y2_local - y3_local

        # the second number is on the top
        percent = (0. - py1_local) / dy
        abs_percent_shifted = abs(percent - 0.5)
        #print('  percent = %s' % percent)
        #print('  abs_percent_shifted = %s' % abs_percent_shifted)

        # catching the case where all edges will intersect with the plane
        # if the edges are extended to infinity
        #
        # a "valid" percent is ranged from [0.-tol, 1.+tol], so:
        # b = [0.-tol, 1.+tol] - 0.5 = [-0.5-tol, 0.5+tol] # is the same thing
        # in_range = abs(b) < 0.5+tol
        #
        in_range = abs_percent_shifted < 0.5 + plane_atol
        if not in_range:
            #print('  **too big...\n')
            continue
        cut_edgei = [inid_a, inid_b]
        cut_edgei.sort()
        avg_local = p2_local  * percent + p1_local  * (1 - percent)
        avg_global = p2_global * percent + p1_global * (1 - percent)
        #projected_points.append(avg_global)

        xl, yl, zl = avg_local
        xg, yg, zg = avg_global

        local_points_temp.append(avg_local)
        global_points_temp.append(avg_global)

        result1 = nodal_result[inid_a]
        result2 = nodal_result[inid_b]
        resulti = result2  * percent + result1  * (1 - percent)

        print('  inid1=%s inid2=%s edge1=%s' % (inid1, inid2, str(edge1)))
        print('    xyz1_local=%s xyz2_local=%s' % (xyz1_local, xyz2_local))
        print('    avg_local=%s' % avg_local)
        print('    avg_global=%s' % avg_global)
        sid = 1
        out_grid = ['GRID', nid_new, None, ] + list(avg_local)
        out_temp = ['TEMP', sid, nid_new, resulti] #+ resulti.tolist()
        msg += print_card_8(out_grid)
        msg += print_card_8(out_temp)
        #print('  ', out_grid)
        #print('  plane_atol=%s dy=%s\n' % (plane_atol, dy))

        geometry_temp.append([eid] + cut_edgei)
        results_temp.append([xl, yl, zl, xg, yg, zg, resulti])  # TODO: doesn't handle results of length 2+
        i_values.append(i)
        percent_values.append(percent)
        nid_new += 1

    #p1 = global_points[-2]
    #p2 = global_points[-1]
    #dxyz = np.linalg.norm(p2 - p1)
    if _is_dot(i_values, percent_values, plane_atol):
        print('dot!!!')
        mid = 2
        return eid_new, nid_new

    fbdf.write(msg)
    local_points.extend(local_points_temp)
    global_points.extend(global_points_temp)
    geometry.extend(geometry_temp)
    result.extend(results_temp)

    #projected_points = np.array(projected_points)
    #p1 = projected_points[0, :]
    #p2 = projected_points[1, :]

    #min_edge_length = min(lengths)
    # hack to get rid of dot intersections
    #dist = np.linalg.norm(p2 - p1)
    #if dist < min_edge_length / 2.:
        ##print(projected_points)
        #print('removing dot...inid1=%s inid2=%s d=%s mel=%s' % (
            #inid1, inid2, dist, min_edge_length))
        #for unused_i in range(2):
            #global_points.pop()
            #local_points.pop()
            #geometry.pop()
            #result.pop()
            #return eid_new, nid_new

    #print('   cut_edge =', cut_edge)
    # if there are 3 nodes in the cut edge, it's fine
    # we'll take the first two
    conrod = ['CONROD', eid, nid_a_prime, nid_b_prime, mid, area, J]
    #print('  ', conrod)
    fbdf.write(print_card_8(conrod))

    eid_new += 1
    nid_new += 2
    return eid_new, nid_new

def _is_dot(ivalues, percent_values, plane_atol):
    """
    we don't want dots
    """
    percent_array = np.array(percent_values)
    if ivalues == [0, 1]:  # source
        is_dot = np.allclose(percent_array, 0., atol=plane_atol)
        print('source; percents=%s is_dot=%s' % (percent_array, is_dot))
    elif ivalues == [0, 2]:  # corner
        shifted_array = np.abs(percent_array - 0.5)
        is_dot = np.allclose(shifted_array, 0.5, atol=plane_atol)
        print('corner; percents=%s is_dot=%s' % (percent_array, is_dot))
    elif ivalues == [1, 2]:  # sink
        is_dot = np.allclose(percent_array, 1., atol=plane_atol)
        print('sink; percents=%s is_dot=%s' % (percent_array, is_dot))
    else:
        raise RuntimeError('incorrect ivalues=%s' % ivalues)
    return is_dot
