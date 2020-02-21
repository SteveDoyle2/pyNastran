"""
defines:
 - local_points_array, global_points_array, result_array = cut_edge_model_by_axes(
        bdf_filename, view_up, p1, p2, tol,
        nodal_result, plane_atol=1e-5,
        plane_bdf_filename=None)
 - local_points_array, global_points_array, result_array = cut_edge_model_by_coord(
        bdf_filename, coord, tol,
        nodal_result, plane_atol=1e-5,
        plane_bdf_filename=None)
 - slice_edges(xyz_cid0, xyz_cid, edges, nodal_result, plane_atol=1e-5,
               plane_bdf_filename=None)

"""
from io import StringIO
from typing import Optional, Tuple, Any
import numpy as np
from pyNastran.nptyping import NDArrayNint, NDArrayN2int, NDArray3float, NDArrayN3float
from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model


def cut_edge_model_by_coord(bdf_filename: str,
                            coord,
                            tol: float,
                            nodal_result,
                            plane_atol: float=1e-5,
                            csv_filename: Optional[str]=None,
                            header: str='x, y, z, Cp') -> Tuple[Any, Any, Any]:
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
    header : str
        the header for the csv; x, y, z, and nodal result names

    """
    assert isinstance(tol, float), tol
    nids, xyz_cid0, edges = _setup_edges(bdf_filename)
    local_points_array, global_points_array, result_array = _cut_edge_model_by_coord(
        nids, xyz_cid0, edges, coord, tol,
        nodal_result, plane_atol=plane_atol)
    if csv_filename:
        export_edge_cut(csv_filename, result_array, header=header)
    return local_points_array, global_points_array, result_array

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

def _setup_edges(bdf_filename: str) -> Tuple[NDArrayNint, NDArrayN3float, NDArrayN2int]:
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
    edges = edge_to_eid_map.keys()
    return nids, xyz_cid0, edges

def _cut_model(nids: NDArrayNint,
               xyz_cp: NDArrayN3float,
               edges: NDArrayN2int,
               view_up: NDArray3float, p1: NDArray3float, p2: NDArray3float,
               tol: float,
               nodal_result, plane_atol=1e-5, plane_bdf_filename=None):
    """
    Helper method for cut_edge_model_by_axes

    Parameters
    ----------
    nids : (nnodes,) int ndarray
        the nodes
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
    plane_bdf_filename : str; default=None
        optionally write a BDF of the cut

   Returns
    -------
    local_points_array : (N, 3) float ndarray
        the xyz points in the cutting plane coordinate system
    global_points_array : (N, 3) float ndarray
        the xyz points in the global xyz coordinate system
    result_array : (N, 7) float ndarray
        inid, x, y, z, xg, yg, zg, result

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
        plane_atol=plane_atol, plane_bdf_filename=plane_bdf_filename)
    return local_points_array, global_points_array, result_array

def _cut_edge_model_by_coord(nids, xyz_cid0, edges, coord, tol,
                             nodal_result, plane_atol=1e-5,
                             plane_bdf_filename=None):
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
    plane_bdf_filename : str; default=None
        optionally write a BDF of the cut

   Returns
    -------
    local_points_array : (N, 3) float ndarray
        the xyz points in the cutting plane coordinate system
    global_points_array : (N, 3) float ndarray
        the xyz points in the global xyz coordinate system
    result_array : (N, 7) float ndarray
        inid, x, y, z, xg, yg, zg, result

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
        xyz_cid0, xyz_cid, iclose_edges_array, nodal_result, plane_atol=plane_atol,
        plane_bdf_filename=plane_bdf_filename)

    #print(coord)
    return local_points_array, global_points_array, result_array

def slice_edges(xyz_cid0: NDArrayN3float, xyz_cid: NDArrayN3float, edges, nodal_result, plane_atol=1e-5,
                plane_bdf_filename: Optional[str]=None) -> Tuple[Any, Any, Any]:
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
    plane_bdf_filename : str; default=None
        optionally write a BDF of the cut

    Returns
    -------
    local_points_array : (N, 3) float ndarray
        the xyz points in the cutting plane coordinate system
    global_points_array : (N, 3) float ndarray
        the xyz points in the global xyz coordinate system
    result_array : (N, 7) float ndarray
        inid, x, y, z, xg, yg, zg, result

    TODO: split result_array, so we don't have mixed int/float/complex
          results being all casted to the highest data type

    """
    #plane_bdf_filename = 'plane_edge.bdf'

    fbdf = StringIO()
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
            if plane_bdf_filename:
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
            if nodal_result is not None:
                result1 = nodal_result[inid1]
                result2 = nodal_result[inid2]
                result.append([inid1, x1, y1, z1, x1g, y1g, z1g, result1])
                result.append([inid2, x2, y2, z2, x2g, y2g, z2g, result2])
            else:
                result.append([inid1, x1, y1, z1, x1g, y1g, z1g])
                result.append([inid2, x2, y2, z2, x2g, y2g, z2g])

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
        if plane_bdf_filename:
            fbdf.write('$------\n')

    if plane_bdf_filename:
        fbdf.close()

    if 0:  # pragma: no cover
        plane_bdf_filename = 'plane_edge2.bdf'
        if plane_bdf_filename:
            cid = 0
            nid_new = 1
            eid_new = 1
            mid = 1
            area = 1.0
            J = 1.0

            ipoint = 0
            with open(plane_bdf_filename, 'w') as fbdf:
                fbdf.write('$pyNastran: punch=True\n')
                fbdf.write('MAT1,1,3.0e7,,0.3\n')

                while ipoint < len(local_points):
                    xyz1_local = local_points[ipoint]
                    #resulti = result[ipoint]
                    print(xyz1_local)
                    nid = xyz1_local[0]
                    if nid == 0:
                        ipoint += 1
                        continue

                    xyz2_local = local_points[ipoint+1]
                    out_grid1 = ['GRID', nid_new, cid, ] + list(xyz1_local)
                    out_grid2 = ['GRID', nid_new + 1, cid, ] + list(xyz2_local)
                    conrod = ['CONROD', eid_new, nid_new, nid_new + 1, mid, area, J]
                    conm2 = ['CONM2', eid_new+1, nid_new, 0, 100.]
                    fbdf.write(print_card_8(out_grid1))
                    fbdf.write(print_card_8(out_grid2))
                    fbdf.write(print_card_8(conrod))
                    fbdf.write(print_card_8(conm2))
                    fbdf.write('$------\n')
                    ipoint += 2
                    nid_new += 2
                    eid_new += 1

    # a hack to avoid making complicated code to sometimes write to the bdf
    if plane_bdf_filename:
        with open(plane_bdf_filename, 'w') as bdf_file:
            bdf_file.write(fbdf.getvalue())

    local_points_array = np.array(local_points)
    global_points_array = np.array(global_points)
    result_array = np.array(result)
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

    #if fbdf is not None:
    out_grid1 = ['GRID', nid_new, None, ] + list(xyz1_local)
    out_grid2 = ['GRID', nid_new+1, None, ] + list(xyz2_local)
    conrod = ['CONROD', eid_new, nid_new, nid_new + 1, mid, area, J]
    #conm2 = ['CONM2', eid_new, nid_new, 0, 100.]

    #fbdf.write('$ interpolated\n')
    fbdf.write(print_card_8(out_grid1))
    fbdf.write(print_card_8(out_grid2))
    fbdf.write(print_card_8(conrod))
    #fbdf.write(print_card_8(conm2))

    local_points.append(avg_local)
    global_points.append(avg_global)

    xl, yl, zl = avg_local
    xg, yg, zg = avg_global
    if nodal_result is not None:
        result1 = nodal_result[inid1]
        result2 = nodal_result[inid2]
        resulti = result2  * percent + result1  * (1 - percent)
        result.append([0, xl, yl, zl, xg, yg, zg, resulti])
    else:
        result.append([0, xl, yl, zl, xg, yg, zg])

    nid_new += 2
    eid_new += 1
    return eid_new, nid_new

def export_edge_cut(csv_filename, result_array, header='x, y, z, Cp'):
    """
    Writes a face cut file of the format:

        x, y, z, Cp
        ...
    """
    #nints = geometry_array.shape[1]
    #nfloats = results_array.shape[1]
    #max_int = geometry_array.max()
    #len_max_int = len(str(max_int))
    #fmt = ('%%%ii,' % (len_max_int)) * nints + '%19.18e,' * nfloats
    #fmt = fmt.rstrip(',')
    #print(fmt)
    #X = np.concatenate((geometry_array, results_array), axis=1)
    np.savetxt(csv_filename, result_array, delimiter=',', newline='\n', header=header,
               footer='', comments='# ') # , encoding=None # numpy 1.14
