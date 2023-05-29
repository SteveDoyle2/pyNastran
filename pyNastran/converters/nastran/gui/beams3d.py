"""creates 3d beams"""
from __future__ import annotations
from collections import defaultdict
from typing import Union, Optional, TYPE_CHECKING

import numpy as np
from numpy.linalg import norm
from vtk import vtkPoints, VTK_FLOAT, vtkIdList

from pyNastran.gui.vtk_interface import vtkUnstructuredGrid, VTK_POLYHEDRON
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk
from pyNastran.bdf.cards.elements.bars import rotate_v_wa_wb
from pyNastran.bdf.cards.elements.beam_connectivity import (
    bar_setup, rod_setup, tube_setup,
    chan_setup, chan1_setup, box_setup, i_setup, i1_setup,
    t_setup, t1_setup, t2_setup,
    h_setup, l_setup, z_setup,
    hexa_setup, hat_setup,
    _transform_points as transform_points,
    Faces,
)
BEAM_SETUP_MAP = {
    'BAR' : bar_setup,
    'ROD' : rod_setup,
    'TUBE' : tube_setup,
    #'TUBE2' : tube2_setup,
    'BOX' : box_setup,
    #'BOX1' : box1_setup,
    'L' : l_setup,
    'CHAN' : chan_setup,
    'CHAN1' : chan1_setup,
    #'CHAN2' : chan2_setup,
    #'CROSS' : cross_setup,
    'T' : t_setup,
    'T1' : t1_setup,
    'T2' : t2_setup,
    'I' : i_setup,
    'I1' : i1_setup,
    'H' : h_setup,
    'HAT' : hat_setup,
    #'HAT1' : hat1_setup,
    'HEXA' : hexa_setup,
    'Z' : z_setup,
    #'DBOX' : dbox_setup,
}

if TYPE_CHECKING:  # pragma: no cover
    #from pyNastran.nptyping_interface import NDArray3float
    from pyNastran.bdf.bdf import BDF, CBAR, CBEAM

def get_bar_nids(model: BDF, bar_beam_eids: list[int]) -> tuple[list[int],
                                                                dict[int, tuple[int, int]]]:
    """gets the bar nids"""
    nids_set = set([])
    nid_release_map_default = defaultdict(list)
    for eid in bar_beam_eids:
        elem = model.elements[eid]  # type: Union[CBAR, CBEAM]
        nid1, nid2 = elem.node_ids
        nids_set.update([nid1, nid2])

        if elem.pa != 0:
            nid_release_map_default[nid1].append((eid, elem.pa))
        if elem.pb != 0:
            nid_release_map_default[nid2].append((eid, elem.pb))

    nids = list(nids_set)
    nids.sort()
    nid_release_map = dict(nid_release_map_default)
    return nids, nid_release_map

def get_beam_sections_map(model: BDF,
                          bar_beam_eids: list[int]) -> dict[int, list[int]]:
    """gets the beams sorted by property_id that can be represented as 3d elements"""
    bar_pid_to_eids_default = defaultdict(list)
    int_offts = []
    for eid in bar_beam_eids:
        elem = model.elements[eid]  # type: Union[CBAR, CBEAM]
        if isinstance(elem.offt, int):
            int_offts.append(eid)
            continue

        pid_ref = elem.pid_ref
        if pid_ref is None:
            pid_ref = model.Property(elem.pid)
        assert not isinstance(pid_ref, int), elem
        ptype = pid_ref.type
        if ptype in {'PBEAML', 'PBARL'}:
            pass
        elif ptype in {'PBEAM', 'PBAR'}:
            continue
        else:
            raise NotImplementedError(pid_ref)
        beam_type = pid_ref.Type
        if beam_type not in BEAM_SETUP_MAP:
            model.log.warning(f'The following beam section is not supported:\n{pid_ref}')
            continue
        pid = pid_ref.pid
        bar_pid_to_eids_default[pid].append(eid)
    bar_pid_to_eids = dict(bar_pid_to_eids_default)
    if int_offts:
        int_offts.sort()
        model.log.warning(f'The following CBAR/CBEAMs have integer OFFTs: {int_offts}')

    return bar_pid_to_eids

def create_3d_beams(model: BDF,
                    bar_pid_to_eids: dict[int, list[int]]) -> Optional[vtkUnstructuredGrid]:
    if len(bar_pid_to_eids) == 0:
        return None
    ugrid = vtkUnstructuredGrid()
    node0 = 0
    points_list = []
    eids_bad = []
    for pid, eids in bar_pid_to_eids.items():
        pid_ref = model.properties[pid]
        ptype = pid_ref.type
        bar_type = pid_ref.beam_type
        if bar_type == {'BAR', 'TUBE', 'TUBE2', 'ROD', 'DBOX', 'HAT1', 'BOX1'}:# TODO: wut?
            continue

        if ptype == 'PBARL':
            dim1 = dim2 = pid_ref.dim
        elif ptype == 'PBEAML':
            dim1 = pid_ref.dim[0, :]
            dim2 = pid_ref.dim[-1, :]
        else:
            raise NotImplementedError(pid_ref)
            #dim1 = dim2 = None
            #return node0

        try:
            func = BEAM_SETUP_MAP[bar_type]
        except KeyError:
            raise NotImplementedError(pid_ref)
            #print('skipping 3d bar_type = %r' % bar_type)
            #return node0
        faces, points1, points2 = func(dim1, dim2)
        for eid in eids:
            elem = model.elements[eid]
            (nid1, nid2) = elem.node_ids
            #bar_nids.update([nid1, nid2])
            node1 = model.nodes[nid1]
            node2 = model.nodes[nid2]
            n1 = node1.get_position()
            n2 = node2.get_position()

            # wa/wb are not considered in i_offset
            # they are considered in ihat
            i = n2 - n1
            Li = norm(i)
            ihat = i / Li

            #if elem.pa != 0:
                #nid_release_map[nid1].append((eid, elem.pa))
            #if elem.pb != 0:
                #nid_release_map[nid2].append((eid, elem.pb))

            unused_v, wa, wb, xform = rotate_v_wa_wb(
                model, elem,
                n1, n2, node1, node2,
                ihat, i, eid, Li, model.log)
            if wb is None:
                # one or more of v, wa, wb are bad
                eids_bad.append(eid)
                continue

            #yhat = xform[1, :]
            #zhat = xform[2, :]
            pointsi = transform_points(n1+wa, n2+wb, points1, points2, xform)
            face_idlist = faces_to_element_facelist(faces, node0)
            ugrid.InsertNextCell(VTK_POLYHEDRON, face_idlist)

            dnode = points1.shape[0] * 2
            node0 += dnode
            points_list.append(pointsi)
            #--------------------------------------------
            #bar_typei = get_bar_type(ptype, pid_ref)
            #centroid = (n1 + n2) / 2.
            #bar_types[bar_typei][0].append(eid)
            #bar_types[bar_typei][1].append((centroid, centroid + yhat * Li * scale))
            #bar_types[bar_typei][2].append((centroid, centroid + zhat * Li * scale))

    if eids_bad:
        eids_bad.sort()
        model.log.warning(f'failed transform PBARL/PBEAML for eids={eids_bad}')
    points = _create_vtk_points_from_list(points_list)
    ugrid.SetPoints(points)
    ugrid.Modified()
    return ugrid

def _create_vtk_points_from_list(points_list: list[np.ndarray]) -> vtkPoints:
    #points_array = _make_points_array(points_list)
    points_array = np.vstack(points_list)
    points = vtkPoints()
    vtk_points = numpy_to_vtk(
        num_array=points_array,
        deep=1,
        array_type=VTK_FLOAT,
    )
    points.SetData(vtk_points)
    points.Modified()
    return points

def update_3d_beams(ugrid: vtkUnstructuredGrid,
                    model: BDF,
                    bar_pid_to_eids: dict[int, list[int]]) -> None:
    node0 = 0
    points_list = []
    for pid, eids in bar_pid_to_eids.items():
        pid_ref = model.properties[pid]
        ptype = pid_ref.type
        bar_type = pid_ref.bar_type
        #if bar_type == {'BAR', 'TUBE', 'TUBE2', 'ROD'}:
            #continue

        if ptype == 'PBARL':
            dim1 = dim2 = pid_ref.dim
        elif ptype == 'PBEAML':
            dim1 = pid_ref.dim[0, :]
            dim2 = pid_ref.dim[-1, :]
        else:
            raise NotImplementedError(pid_ref)
            #dim1 = dim2 = None
            #return node0

        try:
            func = BEAM_SETUP_MAP[bar_type]
        except KeyError:
            raise NotImplementedError(pid_ref)
            #print('skipping 3d bar_type = %r' % bar_type)
            #return node0
        faces, points1, points2 = func(dim1, dim2)
        del faces

        for eid in eids:
            elem = model.elements[eid]
            (nid1, nid2) = elem.node_ids
            #bar_nids.update([nid1, nid2])
            node1 = model.nodes[nid1]
            node2 = model.nodes[nid2]
            n1 = node1.get_position()
            n2 = node2.get_position()

            # wa/wb are not considered in i_offset
            # they are considered in ihat
            i = n2 - n1
            Li = norm(i)
            ihat = i / Li

            unused_v, wa, wb, xform = rotate_v_wa_wb(
                model, elem,
                n1, n2, node1, node2,
                ihat, i, eid, Li, model.log)
            if wb is None:
                # one or more of v, wa, wb are bad
                continue

            pointsi = transform_points(n1+wa, n2+wb, points1, points2, xform)
            dnode = points1.shape[0] * 2
            node0 += dnode
            points_list.append(pointsi)
    if node0:
        points = _create_vtk_points_from_list(points_list)
        ugrid.SetPoints(points)
        ugrid.Modified()
    return

def faces_to_element_facelist(faces: Faces, node0: int) -> vtkIdList:
    """creates a series of faces for the custom elements"""
    face_idlist = vtkIdList()

    nfaces = len(faces)
    face_idlist.InsertNextId(nfaces) # Number faces that make up the cell.
    for face in faces: # Loop over all the faces
        #print(face)
        face_idlist.InsertNextId(len(face)) # Number of points in face

        # Insert the pointIds for the face
        #for i in face:
            #face_idlist.InsertNextId(i + node0)
        [face_idlist.InsertNextId(i + node0) for i in face]
    return face_idlist

def get_bar_type(ptype: str, pid_ref):
    """helper method for _get_bar_yz_arrays"""
    if ptype in ['PBAR', 'PBEAM']:
        bar_type = 'bar'
    #if ptype == 'PBAR':
        #bar_type = 'pbar'
    #elif ptype == 'PBEAM':
        #bar_type = 'pbeam'
    elif ptype in ['PBARL', 'PBEAML']:
        bar_type = pid_ref.Type
    elif ptype == 'PBCOMP':
        bar_type = 'pbcomp'
    else:  # pragma: no cover
        raise NotImplementedError(pid_ref)
    return bar_type
