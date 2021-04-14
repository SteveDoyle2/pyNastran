from __future__ import annotations
import os
import numpy as np
import h5py
import vtk
import vtkmodules
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from pyNastran.dev.h5.read_h5 import pyNastranH5
from pyNastran.utils import object_methods, object_stats, object_attributes
from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk, numpy_to_vtkIdTypeArray

def numpy_to_vtk_points(nodes, points=None, dtype='<f', deep=1):
    """common method to account for vtk endian quirks and efficiently adding points"""
    assert isinstance(nodes, np.ndarray), type(nodes)
    if points is None:
        points = vtk.vtkPoints()
        nnodes, ndim = nodes.shape
        assert ndim == 3, nodes.shape
        points.SetNumberOfPoints(nnodes)

        # if we're in big endian, VTK won't work, so we byte swap
        nodes = np.asarray(nodes, dtype=np.dtype(dtype))

    points_array = numpy_to_vtk(
        num_array=nodes,
        deep=deep,
        array_type=vtk.VTK_FLOAT,
    )
    points.SetData(points_array)
    return points


def fill_vtk_unstructured_grid_aero(geom_model) -> vtk.vtkUnstructuredGrid:
    grid_aero = vtk.vtkUnstructuredGrid()
    #grid_aero = vtk.vtk.vtkStructuredGrid()

def fill_vtk_unstructured_grid(geom_model: BDF,
                               vtk_ugrid: vtk.vtkUnstructuredGrid,
                               add_property=True,
                               add_material=True):
    #cell_type_point = 1 # vtk.vtkVertex().GetCellType()
    #cell_type_line = 3 # vtk.vtkLine().GetCellType()
    #cell_type_tri3 = 5
    #cell_type_tri6 = 22
    #cell_type_quad4 = 9
    #cell_type_quad8 = 23
    #cell_type_tetra4 = 10
    #cell_type_tetra10 = 24
    #cell_type_pyram5 = 14 # vtk.vtkPyramid().GetCellType()
    #cell_type_pyram13 = 27 # vtk.vtkQuadraticPyramid().GetCellType()
    #cell_type_hexa8 = 12
    #cell_type_hexa20 = 25
    #cell_type_penta6 = 13
    #cell_type_penta15 = 26
    #nodes, node_ids, nid_map, idtype = _load_nodes(geom_model)
    nodes = geom_model._nodes
    node_ids = geom_model._node_ids
    #nid_map = geom_model._nid_map
    idtype = geom_model._idtype
    vtk_points = numpy_to_vtk_points(nodes, points=None, dtype='<f', deep=1)

    vtk_ugrid.SetPoints(vtk_points)

    etype_nids, cell_offsets_array, cell_types_array, eids, pids = _load_elements(
        geom_model, node_ids, idtype=idtype)
    nelements = len(eids)

    # build the grid
    _elements_to_vtk(vtk_ugrid, etype_nids, cell_offsets_array, cell_types_array)

    # fill the grid with results
    #point_data = vtk_ugrid.GetPointData()

    cell_data = vtk_ugrid.GetCellData()
    eid_array = numpy_to_vtk(eids, deep=0, array_type=None)
    eid_array.SetName('ElementID')
    cell_data.AddArray(eid_array)

    if add_property:
        pid_array = numpy_to_vtk(pids, deep=0, array_type=None)
        pid_array.SetName('PropertyID')
        cell_data.AddArray(pid_array)

    #if add_property:
        #pid_array = numpy_to_vtk(pids, deep=0, array_type=None)
        #pid_array.SetName('PropertyID')
        #cell_data.AddArray(pid_array)

    if add_material:
        psolid_mids = np.full(nelements, -1, dtype='int64')
        pshell_mids = np.full((nelements, 4), -1, dtype='int64')
        thickness = np.full(nelements, np.nan, dtype='float32')
        upids = np.unique(pids)
        is_solid = False
        is_shell = False
        for pid in upids:
            ipid = np.where(pid == pids)
            prop = geom_model.properties[pid]
            if prop.type == 'PSOLID':
                is_solid = True
                psolid_mids[ipid] = prop.mid
            elif prop.type == 'PSHELL':
                is_shell = True
                thickness[ipid] = prop.t
                for imid, mid in enumerate([prop.mid1, prop.mid2, prop.mid3, prop.mid4]):
                    if mid is None:
                        continue
                    pshell_mids[ipid, imid] = mid
            else:
                geom_model.log.warning(f'skipping:\n{prop}')

        if is_solid:
            mid_array = numpy_to_vtk(psolid_mids, deep=0, array_type=None)
            mid_array.SetName('Solid Material')
            cell_data.AddArray(mid_array)

        if is_shell:
            thickness_array = numpy_to_vtk(thickness, deep=0, array_type=None)
            thickness_array.SetName('Shell Thickness')
            cell_data.AddArray(thickness_array)
            for imid in range(4):
                mids = pshell_mids[:, imid]
                shell_mid_array = numpy_to_vtk(mids, deep=0, array_type=None)
                thickness_array.SetName(f'Shell Material {imid+1}')
                cell_data.AddArray(shell_mid_array)
    #print(point_data)
    return eids

def _load_nodes(geom_model: BDF) -> Tuple[np.ndarray,
                                          Dict[int, int],
                                          str]:
    idtype = 'int64'
    nid_map = {}
    if 0:
        nodes = geom_model.nodes
        nnodes = len(nodes)

        i = 0
        node_ids = np.full(nnodes, -1, dtype=idtype)
        points = []
        for nid, grid in sorted(nodes.items()):
            #vtk_points.SetPoint(i, *grid.xyz)
            points.append(grid.xyz)
            nid_map[nid] = i
            node_ids[i] = nid
            i += 1
        del nid, grid
        nodes = np.array(points, dtype='float32')
        del points
    else:
        # TODO: consider CP and CD
        GRID = geom_model.GRID
        node_ids = GRID['ID']
        nodes = GRID['X']
    return nodes, node_ids, nid_map, idtype

def _load_elements(geom_model: BDF, node_ids: np.ndarray, idtype: str='int32'):
    missing_types = set()
    ielem = 0
    ioffset = 0
    etypes = {
        'CTRIA3', 'CQUAD4', 'CTRIA6', 'CQUAD8',
        'CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM',
    }
    nelements = get_number_of_elements(etypes, geom_model)

    cell_types_all = []
    etype_nids_all = []
    cell_offsets_all = []
    eids = np.full(nelements, -1, dtype=idtype)
    pids = np.full(nelements, -1, dtype=idtype)

    card_count = geom_model.card_count
    for card_type in etypes:
        if card_type not in card_count:
            continue
        nelementsi = card_count[card_type]
        set_nnodes = True
        set_pid = True
        # 1D
        #if card_type == 'CONROD':
            #etype = 5 # line
            #elem = geom_model.CONROD
        # 2D
        if card_type == 'CTRIA3':
            etype = 5 # CTRIA3
            elem = geom_model.CTRIA3
        elif card_type == 'CQUAD4':
            etype = 9 # CQUAD4
            elem = geom_model.CQUAD4
        elif card_type == 'CTRIA6':
            etype = 22 # CTRIA6
            elem = geom_model.CTRIA6
        elif card_type == 'CQUAD8':
            etype = 23 # CQUAD8
            elem = geom_model.CQUAD8

        # 3D
        elif card_type == 'CTETRA':
            set_nnodes = False # we're just going to draw the first order geometry
            etype = 10 # CTETRA
            nnodes = 4
            elem = geom_model.CTETRA
            G = elem['G'][:, :nnodes]

            #EID = group['EID']
            #PID = group['PID']
            #NIDS = group['G']
            #DOMAIN_ID = group['DOMAIN_ID']
            #quadratic_nids = NIDS[:, nnodes:]
            #is_large = quadratic_nids.max(axis=1) == 1

        elif card_type == 'CPENTA':
            set_nnodes = False # we're just going to draw the first order geometry
            etype = 13 # CPENTA
            nnodes = 6
            elem = geom_model.CPENTA
            G = elem['G'][:, :nnodes]
        elif card_type == 'CHEXA':
            set_nnodes = False # we're just going to draw the first order geometry
            etype = 12 # CHEXA
            nnodes = 8
            elem = geom_model.CHEXA
            G = elem['G'][:, :nnodes]
        else:
            missing_types.add(elem.type)
            print(card_type, nelementsi)
            continue
        #print(card_type)
        EID = elem['EID']
        if set_pid:
            PID = elem['PID']

        if set_nnodes:
            G = elem['G']
        G2 = np.searchsorted(node_ids, G, side='left', sorter=None)
        if set_nnodes:
            nnodes = G.shape[1]
        etype_nids = np.zeros((nelementsi, 1+nnodes), dtype=idtype)
        etype_nids[:, 0] = nnodes # CQUAD4
        etype_nids[:, 1:] = G2

        eids[ielem:ielem+nelementsi] = EID.flatten()
        pids[ielem:ielem+nelementsi] = PID.flatten()
        cell_typesi = np.ones(nelementsi, dtype=idtype) * etype
        cell_offsetsi = np.arange(ioffset, ioffset+nelementsi*nnodes, nnodes, dtype=idtype)
        ielem += nelementsi
        ioffset += nelementsi*nnodes
        cell_types_all.append(cell_typesi)
        cell_offsets_all.append(cell_offsetsi)
        etype_nids_all.append(etype_nids.ravel())

    if len(etype_nids_all) == 1:
        etype_nids = etype_nids_all[0]
        cell_offsets_array = cell_offsets_all[0]
        cell_types_array = cell_types_all[0]
    else:
        etype_nids = np.hstack(etype_nids_all)
        cell_offsets_array = np.hstack(cell_offsets_all)
        cell_types_array = np.hstack(cell_types_all)
    print('missing_types = ', missing_types)
    return etype_nids, cell_offsets_array, cell_types_array, eids, pids

def _elements_to_vtk(vtk_ugrid: vtk.vtkUnstructuredGrid,
                     etype_nids: np.ndarray,
                     cell_offsets_array: np.ndarray,
                     cell_types_array: np.ndarray):
    # Create the array of cells
    nelements = len(cell_offsets_array)
    cells_id_type = numpy_to_vtkIdTypeArray(etype_nids, deep=1)
    vtk_cells = vtk.vtkCellArray()
    vtk_cells.SetCells(nelements, cells_id_type)

    # Cell types
    deep = False
    vtk_cell_types = numpy_to_vtk(
        cell_types_array, deep=deep,
        array_type=vtk.vtkUnsignedCharArray().GetDataType())

    vtk_cell_offsets = numpy_to_vtk(cell_offsets_array, deep=deep,
                                    array_type=vtk.VTK_ID_TYPE)

    #grid = vtk.vtkUnstructuredGrid()
    vtk_ugrid.SetCells(vtk_cell_types, vtk_cell_offsets, vtk_cells)

def get_number_of_elements(etypes, geom_model: BDF) -> int:
    nelements = 0
    for card_type, nelementsi in geom_model.card_count.items():
        if card_type in etypes:
            nelements += nelementsi
    return nelements
