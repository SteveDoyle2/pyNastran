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
    nid_cp_cd, xyz_cid0, xyz_cp, unused_icd_transform, unused_icp_transform = out
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
    nid_cp_cd, xyz_cid0, xyz_cp, unused_icd_transform, unused_icp_transform = out
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
            face_eids.append(eid)
            faces.append(elem.node_ids)

    #out = model._get_maps(eids=None, map_names=None,
                          #consider_0d=False, consider_0d_rigid=False,
                          #consider_1d=False, consider_2d=True, consider_3d=False)
    #edge_to_eid_map = out['edge_to_eid_map']
    #faces = iterkeys(edge_to_eid_map)
    return nids, xyz_cid0, faces, face_eids

def cut_edge_model_by_coord(bdf_filename, coord, tol,
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
    nids, xyz_cid0, edges = _setup_edges(bdf_filename)
    local_points_array, global_points_array, result_array = _cut_edge_model_by_coord(
        nids, xyz_cid0, edges, coord, tol,
        nodal_result, plane_atol=plane_atol)
    return local_points_array, global_points_array, result_array

def cut_face_model_by_coord(bdf_filename, coord, tol,
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
    nids, xyz_cid0, faces, face_eids = _setup_faces(bdf_filename)
    local_points_array, global_points_array, result_array = _cut_face_model_by_coord(
        nids, xyz_cid0, faces, face_eids, coord, tol,
        nodal_result, plane_atol=plane_atol)
    return local_points_array, global_points_array, result_array

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
        raise NotImplementedError(method)
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
        the plane tolerance
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
    #aaa

    iclose_faces_array = np.searchsorted(nids, close_faces_array.ravel()).reshape(
        close_faces_array.shape)

    #print('iclose_edges_array:')
    #print(iclose_edges_array)

    local_points_array, global_points_array, result_array = slice_faces(
        xyz_cid0, xyz_cid, iclose_faces_array, close_face_eids_array,
        nodal_result, plane_atol=plane_atol)

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

def get_close_faces(faces, face_eids, nids_close):
    """this seems like it could be a lot faster"""
    #n1 = edges[:, 0]
    #n2 = edges[:, 1]
    close_faces = []
    close_face_eids = []
    for eid, face in zip(face_eids, faces):
        if not all((nid for nid in face)):
            continue
        close_faces.append(face)
        close_face_eids.append(eid)
    return close_faces, close_face_eids

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
            result1 = nodal_result[nid1]
            result2 = nodal_result[nid2]
            result.append([nid1, x1, y1, z1, x1g, y1g, z1g, result1])
            result.append([nid2, x2, y2, z2, x2g, y2g, z2g, result2])

            nid_new += 2
            eid_new += 2

        elif is_same_sign:  # Labs == Lpos
            # same sign, so no crossing
            #print('*edge =', edge)
            #print("  xyz1_global=%s xyz2_global=%s" % (xyz1_global, xyz2_global))
            #print("  xyz1_local =%s xyz2_local =%s" % (xyz1_local, xyz2_local))
            continue
        else:
            nid_new, eid_new = _interpolate_bar_to_node(
                nid_new, eid_new, mid, area, J, fbdf,
                nid1, nid2,
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

def _interpolate_bar_to_node(nid_new, eid_new, mid, area, J, fbdf,
                             nid1, nid2,
                             xyz1_local, xyz2_local,
                             xyz1_global, xyz2_global,
                             nodal_result,
                             local_points, global_points, result):
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
    y1_local = xyz1_local[1]
    y2_local = xyz2_local[1]
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
    return nid_new, eid_new

def slice_faces(xyz_cid0, xyz_cid, faces, face_eids, nodal_result, plane_atol=1e-5):
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
    cid = 0
    local_points = []
    global_points = []
    nid_new = 1
    eid_new = 1
    mid = 1
    area = 1.0
    J = 1.0

    result = []
    tri_faces = faces
    for eid, face in zip(face_eids, tri_faces):
        if len(face) == 4:
            continue
        (nid1, nid2, nid3) = face
        xyz1_local = xyz_cid[nid1]
        xyz2_local = xyz_cid[nid2]
        xyz3_local = xyz_cid[nid3]
        xyz1_global = xyz_cid0[nid1]
        xyz2_global = xyz_cid0[nid2]
        xyz3_global = xyz_cid0[nid3]
        y1_local = xyz1_local[1]
        y2_local = xyz2_local[1]
        y3_local = xyz3_local[1]
        y_local = xyz_cid0[face, 1]
        abs_y_local = np.abs(y_local)
        abs_y1_local = np.abs(y1_local)
        abs_y2_local = np.abs(y2_local)
        abs_y3_local = np.abs(y3_local)
        is_same_sign = np.sign(y1_local) == np.sign(y2_local) == np.sign(y3_local)
        is_far_from_plane = all(np.greater(abs_y_local, plane_atol))
        #print("  y_local = %s" % y_local)
        #print('  xyz1-local=%s xyz2-local=%s' % (xyz1_local, xyz2_local))
        #print('  xyz1-global=%s xyz2-global=%s' % (xyz1_global, xyz2_global))
        #print('  is_same_sign=%s is_far_from_plane=%s' % (is_same_sign, is_far_from_plane))
        if is_far_from_plane and is_same_sign:
            #print('  far-face=%s' % (face))
            #print('skip y1_local=%.3f y2_local=%.3f plane_atol=%.e' % (
                #y1_local, y2_local, plane_atol))
            continue
        elif np.allclose(y_local, y2_local, atol=plane_atol):
            print('  on_edge?-face=%s' % (face))
            raise NotImplementedError('on edge?')
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
            result.append([nid2, x2, y2, z2, x2g, y2g, z2g, result2])

            nid_new += 2
            eid_new += 2

        elif is_same_sign:  # Labs == Lpos
            print('  same sign-face=%s' % (face))
            # same sign, so no crossing
            #print('*edge =', edge)
            #print("  xyz1_global=%s xyz2_global=%s" % (xyz1_global, xyz2_global))
            #print("  xyz1_local =%s xyz2_local =%s" % (xyz1_local, xyz2_local))
            continue
        else:
            print('  intersection-eid=%s face=%s' % (eid, face))
            print('  is_same_sign=%s is_far_from_plane=%s' % (is_same_sign, is_far_from_plane))
            nid_new, eid_new = _interpolate_face_to_bar(
                eid, nid_new, eid_new, mid, area, J, fbdf,
                nid1, nid2, nid3,
                xyz1_local, xyz2_local, xyz3_local,
                xyz1_global, xyz2_global, xyz3_global,
                nodal_result,
                local_points, global_points, result, plane_atol)
        assert len(local_points) == len(result)
        fbdf.write('$------\n')
    fbdf.close()
    local_points_array = np.array(local_points)
    global_points_array = np.array(global_points)
    result_array = np.array(result)
    return local_points_array, global_points_array, result_array

def _interpolate_face_to_bar(eid, nid_new, eid_new, mid, area, J, fbdf,
                             nid1, nid2, nid3,
                             xyz1_local, xyz2_local, xyz3_local,
                             xyz1_global, xyz2_global, xyz3_global,
                             nodal_result,
                             local_points, global_points, result, plane_atol):
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
    edgesi = (
        ((nid1, xyz1_local, xyz1_global), (nid2, xyz2_local, xyz2_global)), # edge 1-2
        ((nid2, xyz2_local, xyz2_global), (nid3, xyz3_local, xyz3_global)), # edge 2-3
        ((nid1, xyz1_local, xyz1_global), (nid3, xyz3_local, xyz3_global)), # edge 1-3
    )
    cut_edge = []
    nid_a_prime = nid_new
    nid_b_prime = nid_new + 1
    for ((nid_a, p1_local, p1_global), (nid_b, p2_local, p2_global)) in edgesi:
        print('  nid_a=%s, p1_local=%s, p1_global=%s' % (nid_a, p1_local, p1_global))
        print('  nid_b=%s, p2_local=%s, p2_global=%s' % (nid_b, p2_local, p2_global))
        py1_local = p1_local[1]
        py2_local = p2_local[1]
        dy = py2_local - py1_local
        if np.allclose(dy, 0.0, atol=plane_atol):
            continue

        #dy21 = y2_local - y1_local
        #dy31 = y3_local - y1_local
        #dy23 = y2_local - y3_local

        # the second number is on the top
        percent = (0. - py1_local) / dy
        abs_percent_shifted = abs(percent - 0.5)
        print('  percent = %s' % percent)
        print('  abs_percent_shifted = %s' % abs_percent_shifted)

        # catching the case where all edges will intersect with the plane
        # if the edges are extended
        #
        # a "valid" percent is ranged from [0.-tol, 1.+tol], so:
        # b = [0.-tol, 1.+tol] - 0.5 = [-0.5-tol, 0.5+tol] # is the same thing
        # in_range = abs(b) < 0.5+tol
        #
        in_range = abs_percent_shifted < 0.5 + plane_atol
        if not in_range:
            #print('  **too big...\n')
            continue
        cut_edge.append((nid_a, nid_b))
        #percent31 = (0. - y1_local) / dy31
        #percent23 = (0. - y3_local) / dy23
        avg_local = p2_local  * percent + p1_local  * (1 - percent)
        avg_global = p2_global * percent + p1_global * (1 - percent)
        xl, yl, zl = avg_local
        xg, yg, zg = avg_global

        local_points.append(avg_local)
        global_points.append(avg_global)

        result1 = nodal_result[nid_a]
        result2 = nodal_result[nid_b]
        resulti = result2  * percent + result1  * (1 - percent)

        #print('  nid1=%s nid2=%s edge=%s' % (nid1, nid2, str(edge)))
        #print('    xyz1_local=%s xyz2_local=%s' % (xyz1_local, xyz2_local))
        #print('    avg_local=%s' % avg_local)
        #print('    avg_global=%s' % avg_global)
        sid = 1
        out_grid = ['GRID', nid_new, None, ] + list(avg_local)
        out_temp = ['TEMP', sid, nid_new, resulti] #+ resulti.tolist()
        fbdf.write(print_card_8(out_grid))
        fbdf.write(print_card_8(out_temp))
        print('  ', out_grid)
        print('  plane_atol=%s dy=%s\n' % (plane_atol, dy))

        result.append([eid, nid_new, xl, yl, zl, xg, yg, zg, resulti])
        nid_new += 1

    print('   cut_edge =', cut_edge)
    # if there are 3 nodes in the cut edge, it's fine
    # we'll take the first two
    conrod = ['CONROD', eid, nid_a_prime, nid_b_prime, mid, area, J]
    print('  ', conrod)
    #conm2 = ['CONM2', eid_new, nid_new, 0, 100.]

    fbdf.write(print_card_8(conrod))
    #fbdf.write(print_card_8(conm2))

    xl, yl, zl = avg_local
    xg, yg, zg = avg_global
    result1 = nodal_result[nid_a_prime]
    result2 = nodal_result[nid_b_prime]
    resulti = result2  * percent + result1  * (1 - percent)

    nid_new += 2
    eid_new += 1
    return nid_new, eid_new

