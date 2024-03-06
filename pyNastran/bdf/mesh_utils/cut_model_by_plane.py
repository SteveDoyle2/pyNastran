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
from __future__ import annotations
import os
from itertools import count
from typing import Union, Any, TYPE_CHECKING

import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model
from pyNastran.bdf.mesh_utils.cut_edge_model_by_plane import cut_edge_model_by_coord
from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.bdf.cards.coordinate_systems import Coord, xyz_to_rtz_array, rtz_to_xyz_array

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, CTRIA3, CQUAD4
    from pyNastran.bdf.cards.coordinate_systems import Coord
    from pyNastran.nptyping_interface import NDArrayNint, NDArray3float, NDArrayNfloat

def get_nid_cd_xyz_cid0(model: BDF) -> tuple[NDArrayNint, NDArrayNint,
                                             dict[int, NDArrayNint], NDArray3float]:
    out = model.get_displacement_index_xyz_cp_cd()
    icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
    nids = nid_cp_cd[:, 0]
    nid_cd = nid_cp_cd[:, [0, 2]]
    xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
        xyz_cp, nids, icp_transform,
        cid=0)
    return nids, nid_cd, icd_transform, xyz_cid0

def get_element_centroids(model: BDF,
                          idtype: str='int32',
                          fdtype: str='float32') -> tuple[NDArrayNint, NDArray3float]:
    """gets the element ids and their centroids"""
    eids_list = []
    element_centroids_cid0_list = []
    springs = (
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
    )
    for eid, elem in sorted(model.elements.items()):
        eids_list.append(eid)
        if elem.type in springs:
            centroid = np.zeros(3)
            nnodes = 0
            for node_ref in elem.nodes_ref:
                try:
                    centroid += node_ref.get_position()
                    nnodes += 1
                except AttributeError:
                    # spoint hack
                    pass
            if nnodes:
                centroid /= nnodes
        elif elem.type in {'CVISC'}:
            centroid = np.zeros(3)
        else:
            centroid = elem.Centroid()
        element_centroids_cid0_list.append(centroid)

    eids = np.array(eids_list, dtype=idtype)
    element_centroids_cid0 = np.array(element_centroids_cid0_list, dtype=fdtype)
    return eids, element_centroids_cid0

def get_stations(model: BDF,
                 p1: NDArray3float, p2: NDArray3float, p3: NDArray3float,
                 zaxis: NDArray3float,
                 method: str='Vector',
                 cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                 nplanes: int=20) -> tuple[NDArray3float, NDArray3float, NDArray3float,
                                           NDArray3float, NDArray3float,
                                           CORD2R,
                                           NDArray3float, NDArrayNfloat]:
    """
    Gets the axial stations

    Parameters
    ----------
    p1: (3,) float ndarray
        defines the starting point for the shear, moment, torque plot
    p3: (3,) float ndarray
        defines the end point for the shear, moment, torque plot
    p2: (3,) float ndarray
        defines the XZ plane for the shears/moments (depends on method)
        'Vectors':
           i = p2
        'Z-Axis Projection':
           i = p2 - p1
    zaxis: (3,) float ndarray
        the direction of the z-axis
    cid_p1 / cid_p2 / cid_p3 : int
        the coordinate systems for p1, p2, and p3
    method : str
       'CORD2R':
          zaxis: point on the z-axis
          p2:     point on the xz-plane
       'Vector':
          zaxis:  k vector
          p2:     xz-plane vector
        'Z-Axis Projection':
          zaxis:  point on the z-axis
          p2:     p2 is a point on the xz-plane

    Returns
    -------
    xyz1 / xyz2 / xyz3 : (3,) float ndarray
        the 1=starting 2=ending, 3=normal coordinates of the
        coordinate frames to create in the cid=0 frame
    i / k : (3,) float ndarray
        the i and k vectors of the coordinate system
    coord_out : Coord
        the output coordinate system
    iaxis_march : (3,) float ndarray
        the normalized x-axis that defines the direction to march
    stations : (n,) float ndarray
        the coordinates in the x-axis that will be marched down

    Example
    -------
    For the BWB example, we to calculate an SMT down the global x-axis

    1--------> y
    |
    |
    |
    2, 3
    |
    v x

    # axial
    p1 = np.array([0., 0., 0.]) # origin
    p2 = np.array([1600., 0., 0.]) # xaxis
    p3 = np.array([1600., 0., 0.]) # end
    zaxis = np.array([0., 0., 1.])
    method = 'Z-Axis Projection'

    xyz1, xyz2, xyz3, i, k, coord_out, stations = get_stations(
        model, p1, p2, p3, zaxis,
        method=method, cid_p1=0, cid_p2=0, cid_p3=0,
        cid_zaxis=0, nplanes=100)
    print(stations)

    """
    p1 = np.asarray(p1) # start
    p2 = np.asarray(p2) # xz-plane
    p3 = np.asarray(p3) # end point
    zaxis = np.asarray(zaxis)

    # define a local coordinate system
    xyz1, xyz2, unused_z_global, i, k, origin, zaxis2, xzplane = _p1_p2_zaxis_to_cord2r(
        model, p1, p2, zaxis,
        cid_p1=cid_p1, cid_p2=cid_p2, cid_zaxis=cid_zaxis,
        method=method)
    xyz3 = model.coords[cid_p3].transform_node_to_global(p3)

    try:
        coord_out = CORD2R(-1, origin=origin, zaxis=zaxis2, xzplane=xzplane)
    except Exception:
        msg = f'Cannot create ouput coordinate system.  origin={origin} zaxis={zaxis} xzplane={xzplane}\n'
        #msg += coord_out.get_stats()
        raise ValueError(msg)

    #coord_march = coord_out
    xaxis_march = xyz3 - xyz1
    xaxis_march_norm = np.linalg.norm(xaxis_march)
    if xaxis_march_norm == 0.:
        msg = f'Coincident starting and end points.  dx={xaxis_march_norm} xyz1={xyz1} xyz3={xyz3}\n'
        #msg += coord_out.get_stats()
        raise ValueError(msg)

    iaxis_march = xaxis_march / xaxis_march_norm
    # k has been rotated into the output coordinate frame, so we'll maintain that
    # k is length=1
    assert np.allclose(np.linalg.norm(k), 1.0)
    jaxis_march = np.cross(k, iaxis_march)
    jaxis_march_norm = np.linalg.norm(jaxis_march)
    if jaxis_march_norm == 0.:
        msg = f'Equal k axis and iaxis.  k={str(k)} iaxis_march={str(iaxis_march)}\n'
        #msg += coord_out.get_stats()
        raise ValueError(msg)

    kaxis_march = np.cross(iaxis_march, jaxis_march)
    kaxis_march_norm = np.linalg.norm(kaxis_march)
    if kaxis_march_norm == 0.:
        msg = f'Equal iaxis and jaxis.  k={str(k)} iaxis_march={str(iaxis_march)} jaxis_march={str(jaxis_march)}\n'
        #msg += coord_out.get_stats()
        raise ValueError(msg)


    coord_march = CORD2R(-1, origin=origin, zaxis=origin+kaxis_march, xzplane=origin+iaxis_march)
    #coord_march = CORD2R(None, origin=origin, zaxis=axis_march, xzplane=xzplane)
    #print(coord_out.get_stats())

    xyz1p = coord_march.transform_node_to_local(xyz1) # start
    xyz3p = coord_march.transform_node_to_local(xyz3) # end
    xaxis = xyz3p - xyz1p

    # we want to give the dx the right sign in the coord_out frame
    #i_abs = np.abs(coord_march.i)
    #i_abs_max = i_abs.max()
    #idir = np.where(i_abs == i_abs_max)[0][0]
    #isign = np.sign(coord_march.i[idir])
    dx = xaxis[0] # * isign

    if abs(dx) == 0.:
        msg = f'Coincident starting and end points.  dx={dx} xyz1={xyz1} xyz3={xyz3}\n'
        msg += coord_out.get_stats()
        raise ValueError(msg)
    x_stations_march = np.linspace(0., dx, num=nplanes, endpoint=True)
    assert x_stations_march.shape == (nplanes, ), x_stations_march.shape
    #stations.sort()

    return xyz1, xyz2, xyz3, i, k, coord_out, iaxis_march, x_stations_march

def _setup_faces(bdf_filename: Union[str, BDF]) -> tuple[Any, Any, Any, Any]:
    """helper method"""
    model = get_bdf_model(bdf_filename, xref=False, log=None, debug=False)
    out = model.get_xyz_in_coord_array(cid=0, fdtype='float64', idtype='int32')
    nid_cp_cd, xyz_cid0, unused_xyz_cp, unused_icd_transform, unused_icp_transform = out
    nids = nid_cp_cd[:, 0]
    #eid_to_edge_map, nid_to_edge_map, edge_to_eid_map = create_maps(model)
    #model = BDF()
    face_eids = []
    faces = []
    shells = {
        'CTRIA3', 'CTRIAX', 'CTRIA6', 'CTRIAX6',
        'CQUAD4', 'CQUAD', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUADX8',
        'CSHEAR'}
    for eid, elem in model.elements.items():
        if elem.type in shells:
            if elem.type == 'CQUAD4':
                # split to 2 faces
                n1, n2, n3, n4 = elem.node_ids
                face_eids.append(eid)
                face_eids.append(-eid)
                faces.append((n1, n2, n3))
                faces.append((n1, n3, n4))
            elif elem.type == 'CTRIA3':
                face_eids.append(eid)
                faces.append(elem.node_ids)
            else:
                model.log.debug('skipping %s' % elem.type)

    #out = model._get_maps(eids=None, map_names=None,
                          #consider_0d=False, consider_0d_rigid=False,
                          #consider_1d=False, consider_2d=True, consider_3d=False)
    #edge_to_eid_map = out['edge_to_eid_map']
    return nids, xyz_cid0, faces, face_eids

def cut_face_model_by_coord(bdf_filename: Union[str, BDF], coord: CORD2R, tol,
                            nodal_result, plane_atol=1e-5, skip_cleanup=True,
                            csv_filename=None,
                            plane_bdf_filename='plane_face.bdf',
                            plane_bdf_filename2='plane_face2.bdf',
                            plane_bdf_offset=0.):
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
    plane_bdf_filename : str; default='plane_face.bdf'
        the path to the simplified conrod model

    """
    assert isinstance(tol, float), tol
    nids, xyz_cid0, faces, face_eids = _setup_faces(bdf_filename)
    unique_geometry_array, unique_results_array, rods_array = _cut_face_model_by_coord(
        nids, xyz_cid0, faces, face_eids, coord, tol,
        nodal_result, plane_atol=plane_atol, skip_cleanup=skip_cleanup,
        plane_bdf_filename=plane_bdf_filename,
        plane_bdf_filename2=plane_bdf_filename2,
        plane_bdf_offset=plane_bdf_offset)
    if csv_filename and unique_geometry_array is not None:
        export_face_cut(csv_filename, unique_geometry_array, unique_results_array)
    #print('unique_geometry_array=%s unique_results_array=%s' % (
        #unique_geometry_array, unique_results_array))
    return unique_geometry_array, unique_results_array, rods_array


def export_face_cut(csv_filename: str, geometry_arrays, results_arrays, header: str='') -> None:
    """
    Writes a face cut file of the format:

        header
        Curve 1
        eid, nid1, nid2, x, y, z, Cp
        ...

        Curve 2
        eid, nid1, nid2, x, y, z, Cp
        ...
    """
    assert geometry_arrays is not None
    assert results_arrays is not None
    with open(csv_filename, 'w') as csv_file:
        if header:
            csv_file.write(header)
        for i, geometry_array, results_array in zip(count(), geometry_arrays, results_arrays):
            geometry_array2 = geometry_array[:, [0, 2, 3]]
            assert geometry_array.shape[1] == 4, geometry_array.shape
            nints = geometry_array2.shape[1]
            nfloats = results_array.shape[1]
            #max_int = geometry_array.max()
            #len_max_int = len(str(max_int))
            #fmt = ('%%%ii,' % (len_max_int)) * nints + '%19.18e,' * nfloa# ts
            fmt = '%i,' * nints + '%19.18e,' * nfloats
            fmt = fmt.rstrip(',')
            # eid, nid, nid1, nid2 -> eid, nid1, nid2
            X = np.concatenate((geometry_array2, results_array), axis=1)
            header2 = 'Curve %i\n' % (i+1)
            header2 += 'eid, nid1, nid2, x, y, z, Cp'
            np.savetxt(csv_file, X, fmt=fmt, newline='\n', header=header2,
                       footer='', comments='# ') # , encoding=None # numpy 1.14
            csv_file.write('\n')

def _determine_cord2r(origin, zaxis, xzplane):
    k = zaxis / np.linalg.norm(zaxis)
    iprime = xzplane / np.linalg.norm(xzplane)
    j = np.cross(k, iprime)
    j /= np.linalg.norm(j)
    i = np.cross(j, k)
    return i, k

def _project_z_axis(p1: NDArray3float,
                    p2: NDArray3float,
                    z_global: NDArray3float) -> tuple[NDArray3float, NDArray3float, NDArray3float,
                                                      NDArray3float, NDArray3float]:
    """
    p1-p2 defines the x-axis
    k is defined by the z-axis
    k = z / |z|
    xz = p2 - p1
    xz /= |xz|
    j = z × xz
    j /= |j|
    i = k × j
    """
    x = p2 - p1
    norm_x = np.linalg.norm(x)
    if norm_x == 0.:
        raise RuntimeError(f'p1={p1} and p2={p2} are coincident; distance={norm_x}')
    iprime = x / norm_x
    k = z_global / np.linalg.norm(z_global)
    j = np.cross(k, iprime)
    jhat = j / np.linalg.norm(j)
    i = np.cross(jhat, k)

    origin = p1
    zaxis = p1 + k
    xzplane = p1 + i
    return i, k, origin, zaxis, xzplane

def _project_vectors(p1: NDArray3float,
                     p2: NDArray3float,
                     z_global: NDArray3float) -> tuple[NDArray3float, NDArray3float, NDArray3float,
                                                       NDArray3float, NDArray3float, NDArray3float]:
    """
    p1 defines the origin
    p2 defines the x-axis
    k is defined by the z-axis
    k = z / |z|
    xz = p2
    xz /= |xz|
    j = k × xz
    j /= |j|
    i = k × j
    """
    x = p2
    origin = p1

    norm_z =  np.linalg.norm(z_global)
    if norm_z == 0.:
        raise RuntimeError(f'p1={p1} is has no length')

    norm_x = np.linalg.norm(x)
    if norm_x == 0.:
        raise RuntimeError(f'p2={p2} is has no length')

    iprime = x / norm_x
    k = z_global / norm_z
    j = np.cross(k, iprime)
    jhat = j / np.linalg.norm(j)
    i = np.cross(jhat, k)

    # we want x to be normal to the plane
    origin_zaxis = origin + k
    origin_xzplane = origin + j

    # add the origin
    xyz2 = origin + i
    return xyz2, i, k, origin, origin_zaxis, origin_xzplane

def _p1_p2_zaxis_to_cord2r(model: BDF,
                           p1: NDArray3float, p2: NDArray3float, zaxis: NDArray3float,
                           method: str='Z-Axis Projection',
                           cid_p1: int=0, cid_p2: int=0,
                           cid_zaxis: int=0) -> tuple[NDArray3float, NDArray3float, NDArray3float,
                                                      NDArray3float, NDArray3float,
                                                      NDArray3float, NDArray3float, NDArray3float]:
    """
    Creates the coordinate system that will define the cutting plane

    Parameters
    ----------
    model : BDF
        the bdf object
    method: str; default='Z-Axis Projection'
        method = 'CORD2R'
           p1:    origin
           zaxis: origin + zaxis (zaxis)
           p3:    origin + xzplane
        method = 'Z-Axis Projection'
           p1:    origin
           zaxis: zaxis
           p2:    origin + xzplane
    cid_p1/p2/zaxis: int; default=0
        the coordinate system for the points

    Returns
    -------
    p1 / p2 / p3 : (3,) float ndarray
        A, B, C in the CORD2R in the rid=0 frame
    i, k : the i/k vectors for the coord_out frame
    origin, zaxis, xzplane : (3,) float ndarray
        the CORD2R-ready coordinate system
    """
    p1 = np.asarray(p1) # origin
    p2 = np.asarray(p2) # xz-plane
    zaxis = np.asarray(zaxis)
    #print("coord:")
    #print('  p1 =', p1)
    #print('  p2 =', p2)
    #print('  zaxis =', zaxis)

    coords = model.coords
    xyz1 = coords[cid_p1].transform_node_to_global(p1)
    method_lower = method.lower().strip()
    if method_lower == 'vector':
        origin = xyz1
        xyz2 = origin + p2
        z_global = origin + zaxis
        method_lower = 'cord2r'
    elif method_lower == 'coord id':
        coord = coords[cid_p2]
        origin = xyz1
        xyz2 = origin + coord.i
        z_global = origin + coord.k
        method_lower = 'cord2r'
    else:
        xyz2 = coords[cid_p2].transform_node_to_global(p2)
        z_global = coords[cid_zaxis].transform_node_to_global(zaxis)

    if method_lower == 'cord2r':
        origin = xyz1
        xzplane_ = xyz2 - xyz1
        zaxis_ = z_global - xyz1
        i, k = _determine_cord2r(origin, zaxis_, xzplane_)
        origin_zaxis = z_global # origin + i
        origin_xzplane = xyz2 # origin + k
    #elif method == 'Vectors':
        #xyz2, i, k, origin, origin_zaxis, origin_xzplane = _project_vectors(xyz1, xyz2, z_global)
    elif method_lower == 'z-axis projection':
        i, k, origin, origin_zaxis, origin_xzplane = _project_z_axis(xyz1, xyz2, z_global)
    else:
        raise NotImplementedError(f"method={method!r}; valid_methods=['CORD2R', 'Vector', 'Z-Axis Projection']")
    #print(f'origin={origin}')
    #print(f'origin_zaxis={origin_zaxis}')
    #print(f'origin_xzplane={origin_xzplane}')
    return xyz1, xyz2, z_global, i, k, origin, origin_zaxis, origin_xzplane

#def _merge_bodies(local_points_array, global_points_array, result_array):
    #local_points_dict = {}
    #global_points_dict = {}
    #result_dict = {}
    #return local_points_dict, global_points_dict, result_dict
    #return NotImplementedError()

def _cut_face_model_by_coord(nids, xyz_cid0, faces, face_eids, coord: Coord, tol: float,
                             nodal_result, plane_atol: float=1e-5, skip_cleanup: bool=True,
                             plane_bdf_filename: str='plane_face.bdf',
                             plane_bdf_filename2: str='plane_face2.bdf',
                             plane_bdf_offset: float=0.) -> tuple[Any, Any, Any]:
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
    face_eids : list[int]
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
    plane_bdf_offset : float; default=0.
        ???

    Returns
    -------
    unique_geometry_array : ???
        ???
    unique_results_array : ???
        ???
    rods_array : ???
        ???

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
    unique_geometry_array, unique_results_array, rods_array = slice_faces(
        nids, xyz_cid0, xyz_cid, iclose_faces_array, close_face_eids_array,
        nodal_result, coord, plane_atol=plane_atol, skip_cleanup=skip_cleanup,
        plane_bdf_filename=plane_bdf_filename,
        plane_bdf_filename2=plane_bdf_filename2,
        plane_bdf_offset=plane_bdf_offset)

    #print(coord)
    return unique_geometry_array, unique_results_array, rods_array

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

#def faces_to_tri_faces(face_eids, faces):
    #"""splits quads into tris"""
    #tri_face_eids = []
    #tri_faces = []
    #for eid, face in zip(face_eids, faces):
        #if len(face) == 4:
            #n1, n2, n3, n4 = face
            #tri_face_eids.append(eid)
            #tri_face_eids.append(-eid)
            #tri_faces.append((n1, n2, n3))
            #tri_faces.append((n1, n3, n4))
        #else:
            #tri_faces.append(face)
    #return tri_face_eids, tri_faces

def split_to_trias(model: BDF) -> None:
    elements2 = {}
    neids = len(model.elements)
    for eid, elem in model.elements.items():
        elem_a, elem_b = elem.split_to_ctria3(eid, eid + neids)
        elements2[elem_a.eid] = elem_a
        elements2[elem_b.eid] = elem_b
    model.elements = elements2

def slice_faces(nodes, xyz_cid0, xyz_cid, faces, face_eids, nodal_result,
                coord: Coord, # CORD2R,
                plane_atol: float=1e-5,
                skip_cleanup: bool=True,
                plane_bdf_filename: str='plane_face.bdf',
                plane_bdf_filename2: str='plane_face2.bdf',
                plane_bdf_offset: float=0.):
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
    skip_cleanup : bool; default=True
        cleanup C and O loops
    plane_bdf_offset : float; default=0.
        shifts the plane ylocation...remove...

    Returns
    -------
    unique_geometry_array : ???
        ???
    unique_results_array : ???
        ???
    """
    assert len(face_eids) > 0, face_eids
    rod_elements: list[tuple[int, int, int]] = []
    rod_nids: list[Any] = []
    rod_xyzs: list[np.ndarray] = []

    #cid = 0
    local_points: list[np.ndarray] = []
    global_points: list[np.ndarray] = []
    nid_new = 1
    eid_new = 1
    mid = 1
    area = 1.0
    J = 1.0

    #tri_face_eids, tri_faces = faces_to_tri_faces(face_eids, faces)
    tri_face_eids = face_eids
    tri_faces = faces
    result: list[np.ndarray] = []
    geometry: list[np.ndarray] = []

    #base, ext = os.path.splitext(plane_bdf_filename)
    #plane_bdf_filename2 = base + '2' + ext
    #with open(plane_bdf_filename, 'w') as fbdf, open(plane_bdf_filename2, 'w') as fbdf2:
    fbdf = open(plane_bdf_filename, 'w')
    fbdf.write('$pyNastran: punch=True\n')
    fbdf.write('MAT1,1,3.0e7,,0.3\n')
    fbdf.write('MAT1,2,3.0e7,,0.3\n')

    fbdf2 = open(plane_bdf_filename2, 'w')
    fbdf2.write('$pyNastran: punch=True\n')
    fbdf2.write('MAT1,1,3.0e7,,0.3\n')
    fbdf2.write('MAT1,2,3.0e7,,0.3\n')

    for eid, face in zip(tri_face_eids, tri_faces):
        #if len(face) == 4:
            #print('skipping face=%s' % face)
            #continue
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
            #print('  far-face=%s' % str(face))
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
            #print('  same sign-face=%s' % (face))
            # same sign, so no crossing
            #print('*edge =', edge)
            #print("  xyz1_global=%s xyz2_global=%s" % (xyz1_global, xyz2_global))
            #print("  xyz1_local =%s xyz2_local =%s" % (xyz1_local, xyz2_local))
            continue
        else:
            # a crossing
            #print('  intersection-eid=%s face=%s' % (eid, face))
            #print('  is_same_sign=%s is_far_from_plane=%s' % (is_same_sign, is_far_from_plane))
            eid_new, nid_new = _interpolate_face_to_bar(
                nodes, eid, eid_new, nid_new, mid, area, J,
                fbdf, fbdf2,
                inid1, inid2, inid3,
                xyz1_local, xyz2_local, xyz3_local,
                xyz1_global, xyz2_global, xyz3_global,
                nodal_result,
                local_points, global_points,
                geometry, result,
                rod_elements, rod_nids, rod_xyzs,
                coord, plane_atol, plane_bdf_offset=plane_bdf_offset)

        if len(local_points) != len(result):
            msg = 'lengths are not equal; local_points=%s result=%s' % (
                len(local_points), len(result))
            raise RuntimeError(msg)
        fbdf.write('$------\n')
        fbdf2.write('$------\n')
        #print('----------------------')
    fbdf.close()
    fbdf2.close()

    if len(rod_elements) == 0:
        os.remove(plane_bdf_filename)

    if len(geometry) == 0:
        return None, None, (None, None, None)
    unused_local_points_array = np.array(local_points)
    unused_global_points_array = np.array(global_points)
    results_array = np.array(result)
    #print('*result', result)
    #print('*results_array', results_array, type(results_array))

    #for g in geometry:
        #print(g)
    geometry_array = np.array(geometry, dtype='int32')
    rods_elements_array = np.array(rod_elements, dtype='int32')
    rod_nids_array = np.array(rod_nids, dtype='int32')
    rod_xyzs_array = np.array(rod_xyzs, dtype='float64')

    unique_geometry_array, unique_results_array = _unique_face_rows(
        geometry_array, results_array, nodes, skip_cleanup=skip_cleanup)
    #print('*unique_results_array', unique_results_array, type(unique_results_array))
    return unique_geometry_array, unique_results_array, (rods_elements_array, rod_nids_array, rod_xyzs_array)

def _unique_face_rows(geometry_array: np.ndarray,
                      results_array: np.ndarray,
                      nodes: np.ndarray,
                      skip_cleanup: bool=True) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns
    -------
    unique_geometry_array : (nelements, 3?)
        eid, ???, ???
    unique_results_array : (nelements, 6+nresults)
        xl, yl, zl, xg, yg, zg, resulti
    skip_cleanup : bool; default=True
        cleanup C and O loops

    """
    #print(geometry_array)
    #iedges = geometry_array[:, 1:]
    #geometry_array[:, 1:] = nodes[iedges.flatten()].reshape(iedges.shape)
    #myrow = None

    # eid, nid, inid1, inid2
    # cut it as [eid, nid, | inid, inid2]
    # | = 2
    n = 2
    for irow, row in enumerate(geometry_array):
        #if row[0] == 11029:
            #myrow = irow
            #print('ids=%s' % row[1:])
            #print('nids=%s' % nodes[row[1:]])

        out = nodes[row[n:]]
        geometry_array[irow, n:] = out

    #print('geometry_array:')
    #print(geometry_array)
    # eid, nid, n1, n2
    edges = geometry_array[:, n:]
    #print(edges)
    #print(','.join([str(val) for val in np.unique(edges)]))
    #print('geom =', geometry_array[myrow, :])

    # axis added in numpy 1.13
    unused_unique_edges, unique_edge_index = np.unique(
        edges,
        return_index=True, return_inverse=False, return_counts=False,
        axis=0)
    unique_geometry_array = geometry_array[unique_edge_index, :]
    unique_results_array = results_array[unique_edge_index, :]
    #for i in unique_edge_index:
        #(eid, nid1, nid2) = geometry_array[i, :]
        #if nid1 == 192 or nid2 == 192:
            #print('i=%s eid=%s nid1=%s nid2=%s' % (i, eid, nid1, nid2))
        #print(geometry_array[i, :])

    unused_iedges, unique_geometry_array2, unique_results_array2 = connect_face_rows(
        unique_geometry_array, unique_results_array, skip_cleanup=skip_cleanup)
    #print('iedges =', unused_iedges)
    return unique_geometry_array2, unique_results_array2

def connect_face_rows(geometry_array: np.ndarray,
                      results_array: np.ndarray,
                      skip_cleanup: bool=True,
                      ) -> tuple[list[Any], list[np.ndarray], list[np.ndarray]]:
    """
    Connects the faces by taking the count of how many times each node
    is used.  If a node is not used twice, then it is a starting/ending point,
    so we can find the C-shaped loops and O-shaped loops.  This is not intended
    to handle 3+ connected points, only 1 or 2.

    Parameters
    ----------
    skip_cleanup: bool; default=True
        cleanup C and O loops

    """
    # temp
    nedges = geometry_array.shape[0]
    #print('skip_cleanup=%s' % skip_cleanup)
    if skip_cleanup:
        iedges = np.arange(0, nedges)
        return [iedges], [geometry_array], [results_array]

    # TODO: need to handle the dot case first

    # grab (nid1, nid2) columns
    eid_edges = geometry_array
    eids_backup = geometry_array[:, 0]
    geometry_array[:, 0] = np.arange(0, nedges)
    edges = geometry_array[:, 1:]  # nodes, not inodes

    nodes, counts = np.unique(edges, return_counts=True)
    ibad = np.where(counts != 2)[0]
    include_end = False
    if len(ibad) == 0:
        ibad = np.array([0])
        include_end = True
    #print('ibad = ', ibad)
    node_start_end = nodes[ibad].tolist()
    #print('node_start_end = ', node_start_end)
    if len(ibad) == 2 or 1:
        iedges: list[int] = []
        iedges_all = []
        #print("ibad = %s" % ibad)
        #print("edges:")
        #print(edges)
        #print('nodes =', nodes)
        #print('counts =', counts)
        nid_start = node_start_end[0]
        node_start_end.remove(nid_start)
        #print('nid_start = %s' % nid_start)

        nedges = edges.shape[0]
        all_irows = list(range(nedges)) #np.arange(0, nedges)
        #print('------------')
        while len(all_irows):
            #print('all_irows =',all_irows)
            #print(eid_edges)
            irows = np.where(nid_start == edges)[0]
            if len(irows) == 0:
                # finished an O-ring or a C-ring
                if nid_start in node_start_end:
                    node_start_end.remove(nid_start)
                    nid_start = node_start_end[0]
                    node_start_end.remove(nid_start)
                else:
                    #print("**missing 1")
                    # include_end=True
                    iedges.append(iedges[0])
                    nid_start = edges[0, 0]
                    node_start_end = []

                #print('nid_end=%s' % nid_start)
                #print('*iedges=%s' % iedges)
                iedges_all.append(iedges)
                iedges = []
                continue
            irow = irows[0]
            #print('irow =', irow)

            eid_nodes_row = eid_edges[irow]
            eid = eid_nodes_row[0]
            nodes_list = eid_nodes_row[1:].tolist()
            #print('eid = %s' % eid)
            #print('nodes = %s' % nodes_list)
            nodes_list.remove(nid_start)

            nid_start = nodes_list[0]
            iedges.append(eid)
            #print('nid = %s' % nid_start)
            all_irows.remove(irow)
            #print('iedges =', iedges)
            #print('all_irows =', all_irows)

            eid_edges = eid_edges[all_irows, :]
            edges = edges[all_irows, :]
            nedges = edges.shape[0]
            all_irows = list(range(nedges))
            #print('------------')
        #print('iedges =', iedges)
        if include_end:
            iedges.append(iedges[0])
        if iedges:
            iedges_all.append(iedges)
        geometry_array[:, 0] = eids_backup
        #print('end!!!')
    else:  # pragma: no cover
        print("ibad = %s" % ibad)
        print("ibad.shape", ibad.shape)
        print("len(ibad) = %s" % len(ibad))
        print('nodes =', nodes)
        print('counts =', counts)
        print('unique = ', nodes[ibad])
        print(','.join([str(val) for val in nodes[ibad]]))

    geometry, results = iedges_to_geometry_results(
        iedges_all, geometry_array, results_array)
    #print(results, type(results))
    return iedges_all, geometry, results

def iedges_to_geometry_results(iedges_all,
                               geometry_array: np.ndarray,
                               results_array: np.ndarray,
                               ) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """
    Takes the iedges and slices the geometry/results to isolate the rings
    or C shapes.
    """
    #print("iedges_all =", iedges_all)
    geometry = []
    results = []
    for iedges in iedges_all:
        geometryi = geometry_array[iedges, :]
        resultsi = results_array[iedges]
        geometry.append(geometryi)
        results.append(resultsi)
    return geometry, results

def _face_on_edge(eid: int, eid_new: int,
                  nid_new: int, mid: int,
                  area: float, J: float, fbdf,
                  inid1: int, inid2: int, unused_inid3: int,
                  xyz1_local: np.ndarray, xyz2_local: np.ndarray, unused_xyz3_local: np.ndarray,
                  xyz1_global: np.ndarray, xyz2_global: np.ndarray, unused_xyz3_global: np.ndarray,
                  nodal_result: np.ndarray,
                  local_points, global_points,
                  geometry, result,
                  rod_elements, rod_nids, rod_xyzs,
                  unused_plane_atol: float) -> tuple[int, int]:
    """is this function necessary?"""
    #raise NotImplementedError('on edge-y1')
    #print('  y-sym; nid1=%s nid2=%s edge=%s' % (nid1, nid2, str(edge)))
    #print('     xyz1=%s xyz2=%s' % (xyz1_global, xyz2_global))

    cid = 0

    # TODO: hardcoded and wrong
    nid_a = inid1
    nid_b = inid2

    out_grid1: list[Union[str, int, float]] = ['GRID', nid_new, cid, ] + xyz1_local.tolist()
    out_grid2: list[Union[str, int, float]] = ['GRID', nid_new + 1, cid, ] + xyz2_local.tolist()
    conrod: list[Union[str, int, float]] = ['CONROD', eid_new, nid_new, nid_new + 1, mid, area, J]
    conm2: list[Union[str, int, float]] = ['CONM2', eid_new+1, nid_new, 0, 100.]
    fbdf.write(print_card_8(out_grid1))
    fbdf.write(print_card_8(out_grid2))
    fbdf.write(print_card_8(conrod))
    fbdf.write(print_card_8(conm2))
    rod_elements.append([eid_new, nid_new, nid_new + 1])
    rod_nids.extend([nid_new, nid_new + 1])
    rod_xyzs.append([xyz1_local, xyz2_local])

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

def _interpolate_face_to_bar(nodes,
                             eid: int, eid_new: int, nid_new: int, mid: int,
                             area: float, J: float,
                             fbdf, fbdf2,
                             inid1: int, inid2: int, inid3: int,
                             xyz1_local: np.ndarray, xyz2_local: np.ndarray, xyz3_local: np.ndarray,
                             xyz1_global: np.ndarray, xyz2_global: np.ndarray, xyz3_global: np.ndarray,
                             nodal_result,
                             local_points, global_points,
                             geometry, result,
                             rod_elements, rod_nids, rod_xyzs: list[np.ndarray],
                             coord: CORD2R,
                             plane_atol: float,
                             plane_bdf_offset: float=0.):
    """
    Parameters
    ----------
    plane_bdf_offset : float; default=0.
        ???

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

    As mentioned previously, only two vectors are used (e.g., e12 and e13).
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
    #if eid == 11029:
        #print('eid=%s inid1=%s, inid2=%s, inid3=%s' % (eid, inid1, inid2, inid3))
        #print('nid1=%s, nid2=%s, nid3=%s' % (nodes[inid1], nodes[inid2], nodes[inid3]))

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
    msg2 = ''
    results_temp = []
    geometry_temp = []
    i_values = []
    percent_values = []
    local_points_temp = []
    global_points_temp = []
    is_result = nodal_result is not None
    for i, (edge1, edge2) in enumerate(edgesi):
        (inid_a, p1_local, p1_global) = edge1
        (inid_b, p2_local, p2_global) = edge2
        #print('  inid_a=%s, p1_local=%s, p1_global=%s' % (inid_a, p1_local, p1_global))
        #print('  inid_b=%s, p2_local=%s, p2_global=%s' % (inid_b, p2_local, p2_global))
        py1_local = p1_local[1]
        py2_local = p2_local[1]
        #length = np.linalg.norm(p2_global - p1_global)
        #lengths.append(length)
        dy = py2_local - py1_local
        if np.allclose(dy, 0.0, atol=plane_atol):
            # We choose to ignore the triangle edge on/close to the symmetry plane.
            # Instead, we use the neighboring projected edges as it's more correct.
            # Also, that way do things in a more consistent way.
            #
            continue

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

        #print('  inid1=%s inid2=%s edge1=%s' % (inid1, inid2, str(edge1)))
        #print('    xyz1_local=%s xyz2_local=%s' % (xyz1_local, xyz2_local))
        #print('    avg_local=%s' % avg_local)
        #print('    avg_global=%s' % avg_global)
        sid = 1
        # ['GRID', nid_new, cp=None, x, y, z]
        out_grid = ['GRID', nid_new, None, ] + list(avg_local)
        out_grid2 = ['GRID', nid_new, None, ] + list(avg_global)
        #rod_elements, rod_nids, rod_xyzs
        rod_nids.append(nid_new)
        rod_xyzs.append(avg_local)
        out_grid[4] += plane_bdf_offset
        msg += print_card_8(out_grid)
        msg2 += print_card_8(out_grid2)

        #print('  ', out_grid)
        #print('  plane_atol=%s dy=%s\n' % (plane_atol, dy))
        if is_result:
            result1 = nodal_result[inid_a]
            result2 = nodal_result[inid_b]
            resulti = result2  * percent + result1  * (1 - percent)
            out_temp = ['TEMP', sid, nid_new, resulti] #+ resulti.tolist()
            msg += print_card_8(out_temp)
            msg2 += print_card_8(out_temp)

            geometry_temp.append([eid, nid_new] + cut_edgei)
            # TODO: doesn't handle results of length 2+
            results_temp.append([xl, yl, zl, xg, yg, zg, resulti])
        else:
            geometry_temp.append([eid, nid_new] + cut_edgei)
            results_temp.append([xl, yl, zl, xg, yg, zg])
        i_values.append(i)
        percent_values.append(percent)
        nid_new += 1

    #p1 = global_points[-2]
    #p2 = global_points[-1]
    #dxyz = np.linalg.norm(p2 - p1)
    if _is_dot(i_values, percent_values, plane_atol):
        #print('dot!!!')
        mid = 2
        return eid_new, nid_new

    fbdf.write(msg)
    fbdf2.write(msg2)
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
    fbdf2.write(print_card_8(conrod))
    rod_elements.append([eid, nid_a_prime, nid_b_prime])

    eid_new += 1
    nid_new += 2
    return eid_new, nid_new

def _is_dot(ivalues, percent_values, plane_atol: float) -> np.ndarray:
    """we don't want dots"""
    percent_array = np.array(percent_values)
    if ivalues == [0, 1]:
        #dot_type = 'source'
        is_dot = np.allclose(percent_array, 0., atol=plane_atol)
    elif ivalues == [0, 2]:
        #dot_type = 'corner'
        shifted_array = np.abs(percent_array - 0.5)
        is_dot = np.allclose(shifted_array, 0.5, atol=plane_atol)
    elif ivalues == [1, 2]:
        #dot_type = 'sink'
        is_dot = np.allclose(percent_array, 1., atol=plane_atol)
    else:
        raise RuntimeError('incorrect ivalues=%s' % ivalues)
    #print('%s; percents=%s is_dot=%s' % (dot_type, percent_array, is_dot))
    return is_dot

def calculate_area_moi(model: BDF, rods, normal_plane, thetas,
                       moi_filename=None) -> tuple[Any, Any, Any, Any]:
    """
    Parameters
    ----------
    model : BDF
        the model object
    rods : (eids, nids, xyzs)
        eids : list[int]
            the element id that was split
        nids : (nelements, 2) int ndarray
            the n1, n2 in xyzs that define the cut shell element
        xyzs : (nnodes, 3) float ndarray
            the xyz of the nodes
    normal_plane : (3, ) float ndarray
        the direction of the cut plane
    thetas : dict???
        thetas[eid] = (thetad, Ex, Ey, Gxy)
    moi_filename : str; default=None
        writes a csv file

    Returns
    -------
    total_area
    Isum
    Jsum
    EIsum
    GJsum
    avg_centroid
    """
    rod_elements, rod_nids, rod_xyzs = rods
    eids = np.abs(rod_elements[:, 0])
    neids = len(eids)
    all_nids = rod_nids
    n1 = rod_elements[:, 1]
    n2 = rod_elements[:, 2]
    inid1 = np.searchsorted(all_nids, n1)
    inid2 = np.searchsorted(all_nids, n2)
    xyz1 = rod_xyzs[inid1, :]
    xyz2 = rod_xyzs[inid2, :]
    centroid = (xyz1 + xyz2) / 2.
    length = np.linalg.norm(xyz2 - xyz1, axis=1)
    assert len(length) == neids

    centroid, area, thickness, E = get_element_inertias(
        model, normal_plane, thetas,
        eids, length, centroid)

    # [Ixx, Iyy, Izz, Ixy, Iyz, Ixz]
    I = np.zeros((len(area), 6), dtype='float64') # type: np.ndarray

    # (Ex, Ey, Gxy)
    Ex = E[:, 0]

    total_area = area.sum()
    avg_centroid = (centroid * area[:, np.newaxis]) .sum(axis=0) / total_area
    assert len(avg_centroid) == 3, len(avg_centroid)
    # y corresponds to the station in the plane of the coordinate system
    # and is 0. because we're in the local plane
    x = centroid[:, 0] - avg_centroid[0]
    y = centroid[:, 1] - avg_centroid[1]
    z = centroid[:, 2] - avg_centroid[2]
    xmin = x.min()
    xmax = x.max()
    ixmin = np.where(x == xmin)[0][0]
    ixmax = np.where(x == xmax)[0][0]
    xyz_min = centroid[ixmin, :]
    xyz_max = centroid[ixmax, :]
    d = xyz_max - xyz_min
    dx = d[0]
    dz = d[2]
    theta = np.arctan2(dx, dz)
    thetad = np.degrees(theta)

    nnodes = len(x)
    delta = np.zeros((nnodes, 3))
    delta[:, 1] = thetad

    xyz = np.zeros((nnodes, 3))
    # we're swapping what axes we have to make the transform easier
    xyz[:, 0] = x
    xyz[:, 1] = z
    xyz[:, 2] = 0.
    rtz = xyz_to_rtz_array(xyz)
    rtz2 = rtz + delta
    xyz2 = rtz_to_xyz_array(rtz2)
    x2 = xyz2[:, 0]
    y2 = xyz2[:, 1]
    #z2 = xyz2[:, 2]

    #origin = d
    #zaxis = np.array([0., 1., 0.])
    #xzplane = d
    dxi = x2.max() - x2.min()
    dyi = y2.max() - y2.min()
    #dzi = z2.max() - z2.min()  # zero by definition

    I[:, 0] = area * (x * x)  # Ixx
    I[:, 1] = area * (y * y)  # Iyy
    I[:, 2] = area * (z * z)  # Izz
    I[:, 3] = area * (x * y)  # Ixy
    I[:, 4] = area * (y * z)  # Iyz
    I[:, 5] = area * (x * z)  # Ixz

    Isum = I.sum(axis=0)
    ExIsum = (Ex[:, np.newaxis] * I).sum(axis=0)
    assert len(Isum) == 6, len(Isum)

    if moi_filename is not None:
        dirname = os.path.dirname(moi_filename)
        eid_filename = os.path.join(dirname, 'eid_file.csv')
        _write_moi_file(
            moi_filename, eid_filename,
            eids, n1, n2, xyz1, xyz2, length, thickness, area,
            centroid, avg_centroid, I, E
        )
    return dxi, dyi, total_area, Isum, ExIsum, avg_centroid

def _write_moi_file(moi_filename: str, eid_filename: str,
                    eids, n1, n2, xyz1, xyz2, length, thickness, area,
                    centroid, avg_centroid, I, E) -> None:
    eidi = 1
    mid = 1
    nid0 = max(n1.max(), n2.max()) + 1
    with open(moi_filename, 'w') as bdf_file, open(eid_filename, 'w') as eid_file:
        bdf_file.write('$ pyNastran: punch=True\n')
        bdf_file.write('MAT1,1,3.0e7,,0.3\n')
        conm2 = ['CONM2', 1, nid0]
        grid = ['GRID', nid0, 0, avg_centroid[0], avg_centroid[2], 0.]
        bdf_file.write(print_card_8(grid))
        bdf_file.write(print_card_8(conm2))

        fmt = ('%s,' * 7)[:-1] + '\n'
        eid_file.write('# eid(%i),pid(%i),area,thickness,Ixx,Izz,Ixz\n')
        for eid, n1i, n2i, xyz1i, xyz2i, lengthi, thicknessi, areai, centroidi, Ii, Ei in zip(
                eids, n1, n2, xyz1, xyz2, length, thickness, area, centroid, I, E):
            actual_eid = abs(eid)

            assert nid0 not in [n1i, n2i], (n1i, n2i)
            pidi = actual_eid
            #pid = eidi
            grid1 = ['GRID', n1i, None] + xyz1i.tolist()
            grid2 = ['GRID', n2i, None] + xyz2i.tolist()
            #crod = ['CROD', eidi, pid, n1i, n2i]
            A, J, nsm = Ii
            #prod = ['PROD', pid, mid, A, J, 0., nsm]
            assert eidi > 0, eidi
            conrod = ['CONROD', eidi, n1i, n2i, mid, A, J, 0., nsm]
            bdf_file.write(print_card_8(grid1))
            bdf_file.write(print_card_8(grid2))
            #bdf_file.write(print_card_8(crod))
            #bdf_file.write(print_card_8(prod))
            bdf_file.write(print_card_8(conrod))
            eidi += 1
            #PID | MID |  A  |  J  |  C  | NSM
            eid_file.write(fmt % (eidi, pidi, areai, thicknessi, Ii[0], Ii[1], Ii[2]))

def get_element_inertias(model: BDF, normal_plane: np.ndarray,
                         thetas: dict[int, tuple[float, float, float, float]],
                         eids: list[int],
                         length: list[float],
                         centroid: list[np.ndarray],
                         ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    normal_plane_vector = normal_plane.copy().reshape((3, 1))
    cg_list: list[np.ndarray] = []
    area_list: list[float] = []
    thickness_list: list[float] = []
    E_list: list[tuple[float, float, float]] = []

    for eid, lengthi, centroidi in zip(eids, length, centroid):
        #print(eid, lengthi)
        element = model.elements[eid]
        if element.type in ['CTRIA3', 'CQUAD4']:
            thicknessi, areai, thetad, Ex, Ey, Gxy, nu_xy = _get_shell_inertia(
                element, normal_plane, normal_plane_vector, lengthi)
            thetas[eid] = (thetad, Ex, Ey, Gxy)
            thickness_list.append(thicknessi)
            area_list.append(areai)
            cg_list.append(centroidi)
            E_list.append((Ex, Ey, Gxy))
        else:
            print(element)

    centroid = np.array(cg_list, dtype='float64')
    area = np.array(area_list, dtype='float64')
    thickness = np.array(thickness_list, dtype='float64')
    E = np.array(E_list, dtype='float64')
    return centroid, area, thickness, E

def _get_shell_inertia(element: Union[CTRIA3, CQUAD4],
                       normal_plane: np.ndarray,
                       normal_plane_vector: np.ndarray,
                       lengthi: float,) -> tuple[float, float, float,
                                                 float, float, float, float]:
    """
    Parameters
    ----------
    element : CTRIA3 / CQUAD4
        the object to cut
    normal_plane : (3,) float ndarray
        the normal vector of the cutting plane (should be roughly normal to the element face)
    normal_plane_vector : (3,1) float ndarray
        the normal vector of the cutting plane (should be roughly normal to the element face)
    lengthi : float
        the length the cutting plane makes with the element

    Returns
    -------
    thicknessi : float
        the total thickness of the element
    areai : float
        the cut area of the element
    imat_rotation_angle_deg : float
        the angle between the cutting plane and the normal_plane / normal_plane_vector
        this is NOT the angle of the fiber
    Ex : float
        the moduli normal to the cut plane
    Ey : float
        the moduli parallel to the cut plane (normal to Ex)
    Gxy : float
        the inplane shear moduli
    nu_xy : float
        the correlary to in-plane nu12

    """
    pid_ref = element.pid_ref
    thicknessi = element.Thickness()
    dxyz, centroid, imat, unused_jmat, element_normal = element.material_coordinate_system()
    #print('imat = ', imat)
    #print('normal = ', normal)
    n1, n2, n3 = element_normal
    n12 = n1 * n2
    n13 = n1 * n3
    n23 = n2 * n3
    #  http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf
    # expanding eq 20 into
    #  R(n,theta) = R0 + R1*sin(theta) + R2*cos(theta)
    #R0 = np.array([
        #[n1 ** 2, n12, n13],
        #[n12, n2 ** 2, n23],
        #[n13, n23, n3 ** 2],
    #], dtype='float64')
    R1 = np.array([
        [0., -n3, n2],
        [n3, 0., -n1],
        [-n2, n1, 0.],
    ], dtype='float64')
    R2 = np.array([
        [1 - n1 ** 2, -n12, -n13],
        [-n12, 1 - n2 ** 2, -n23],
        [-n13, -n23, 1 - n3 ** 2],
    ])
    imat = imat.reshape(3, 1)
    #print(normal_plane.shape, R1.shape, imat.shape)
    #a = np.linalg.multi_dot([normal_plane.T, R0, imat])
    b = np.linalg.multi_dot([normal_plane_vector.T, R1, imat])
    c = np.linalg.multi_dot([normal_plane_vector.T, R2, imat])

    #  maximize m' dot p = p.T dot m
    # m' = R dot m
    #    = a + b*sin(theta) + c*cos(theta)
    #  d/d(theta) = b*cos(theta)*sin(theta) = 0
    #
    #  d/d(theta) = b*cos(theta) - c*sin(theta) = 0
    #  b*cos(theta) = c*sin(theta)
    #  tan(theta) = b/c
    #
    # the theta to rotate by in order to orient imat with the normal
    #print(b, c)
    imat_rotation_angle = np.arctan2(b, c).item()
    imat_rotation_angle_deg = np.degrees(imat_rotation_angle)
    if imat_rotation_angle_deg <= -90.:
        imat_rotation_angle_deg += 180.
    elif imat_rotation_angle_deg > 90.:
        imat_rotation_angle_deg -= 180.

    #element_normal = element.Normal()
    # cos(theta) = a o b / (|a| * |b|)
    # |a| = length of normal vector = 1.0
    # |b| = length of normal_plane vector = 1.0
    #
    # cos(theta) = a o b
    # then we take the absolute value because we don't care if the element is +/- theta off

    abs_cos_theta = abs(normal_plane @ element_normal)
    assert isinstance(imat_rotation_angle, float), imat_rotation_angle
    if abs_cos_theta > 0.9:  # <25.8 degrees
        # filter out elements that are in-plane
        thicknessi = 0.
        areai = 0.
        Ex = 0.
        Ey = 0.
        Gxy = 0.
        nu_xy = 0.
    else:
        Ex, Ey, Gxy, nu_xy = pid_ref.get_Ainv_equivalent_pshell(
            imat_rotation_angle_deg, thicknessi, degrees=True)

        #thicknessi = prop.Thickness()
        areai = thicknessi * lengthi

    #pid = pid_ref.pid
    #if pid == 10:
        #import copy
        #pid_ref45 = copy.deepcopy(pid_ref)
        #pid_ref45.mids_ref = [copy.deepcopy(pid_ref.mids_ref[0])]
        #pid_ref45.thetas = [copy.deepcopy(pid_ref.thetas[0])]
        #pid_ref45.thicknesses = [copy.deepcopy(pid_ref.thicknesses[0])]
        #pid_ref45.mids = [copy.deepcopy(pid_ref.mids[0])]
        #pid_ref45.get_thetas()
        #Ex45, Ey45, Gxy45, nu_xy45 = pid_ref45.get_Ainv_equivalent_pshell(
            #imat_rotation_angle_deg, thicknessi)

        #pid_ref0 = copy.deepcopy(pid_ref)
        #pid_ref0.mids_ref = [copy.deepcopy(pid_ref.mids_ref[1])]
        #pid_ref0.thetas = [copy.deepcopy(pid_ref.thetas[1])]
        #pid_ref0.thicknesses = [copy.deepcopy(pid_ref.thicknesses[1])]
        #pid_ref0.mids = [copy.deepcopy(pid_ref.mids[1])]
        #pid_ref0.get_thetas()
        #Ex0, Ey0, Gxy0, nu_xy0 = pid_ref0.get_Ainv_equivalent_pshell(
            #imat_rotation_angle_deg, thicknessi)
        #x = 1

    return thicknessi, areai, imat_rotation_angle_deg, Ex, Ey, Gxy, nu_xy
