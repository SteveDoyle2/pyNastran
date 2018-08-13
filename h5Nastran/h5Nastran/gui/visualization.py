from __future__ import print_function, absolute_import

"""
Quick and easy way to view model.  Not a full blown gui.
"""

from typing import List, Dict
from six import iteritems, itervalues, iterkeys
import numpy as np
from itertools import chain

from pyNastran.bdf.bdf import BDF

import vtk


class VTKData(object):
    def __init__(self, ugrid, nid_list, nid_dict, eid_list, eid_dict):
        self.ugrid = ugrid  # type: vtk.vtkUnstructuredGrid
        self.nid_list = nid_list  # type: List[int]
        self.nid_dict = nid_dict  # type: Dict[int, int]
        self.eid_list = eid_list  # type: List[int]
        self.eid_dict = eid_dict  # type: Dict[int, int]

    def visualize(self):
        ugrid = self.ugrid

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(ugrid)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        renderer = vtk.vtkRenderer()
        rw = vtk.vtkRenderWindow()
        rw.AddRenderer(renderer)

        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(rw)
        renderer.AddActor(actor)

        interactor_style = vtk.vtkInteractorStyleTrackballCamera()

        iren.SetInteractorStyle(interactor_style)

        lut = vtk.vtkLookupTable()
        lut.SetHueRange(2 / 3, 0)
        lut.Build()

        mapper.SetLookupTable(lut)
        # actor.GetProperty().SetRepresentationToWireframe()

        iren.Initialize()
        rw.Render()
        iren.Start()


def _basic_elem(elem, cells, nid_dict):
    nids = elem.nodes
    _cells = [nid_dict[nid] for nid in nids]
    cells.append(len(nids))
    cells.extend(_cells)
    return elem.eid


def _grid(grid, cells, nid_dict):
    nid = grid.nid
    _cells = [nid_dict[nid]]
    cells.append(1)
    cells.extend(_cells)
    return grid.nid


def _rbe2(elem, cells, nid_dict):
    independent = elem.gn

    dependent = elem.Gmi

    nids = []

    for i in range(len(dependent)):
        nids.append(nid_dict[independent])
        nids.append(nid_dict[dependent[i]])

    cells.append(len(nids))
    cells.extend(nids)
    return elem.eid


def _rbe3(elem, cells, nid_dict):
    independent = elem.refgrid
    _dependent = elem.Gijs

    dependent = []

    for dep in _dependent:
        if isinstance(dep, list):
            dependent.extend(dep)
        else:
            dependent.append(dep)

    nids = []

    for i in range(len(dependent)):
        nids.append(nid_dict[independent])
        nids.append(nid_dict[dependent[i]])

    cells.append(len(nids))
    cells.extend(nids)
    return elem.eid


def to_vtk(bdf):
    # type: (BDF) -> VTKData
    import vtk
    # TODO: use pynastran numpy_to_vtk
    from vtk.util.numpy_support import numpy_to_vtk

    nid_list = []
    nid_dict = {}

    node_pos = np.empty((len(bdf.nodes), 3), dtype=np.float64)

    i = 0
    for node in itervalues(bdf.nodes):
        node_pos[i] = node.get_position()
        nid = node.nid
        nid_list.append(nid)
        nid_dict[nid] = i
        i += 1

    _points = vtk.vtkPoints()
    _points.SetData(numpy_to_vtk(node_pos))

    points = vtk.vtkPoints()
    points.DeepCopy(_points)

    cells = []
    cell_types = []
    cell_count = 0
    elem_types = []

    eid_list = []
    eid_dict = {}

    _nastran_to_vtk = nastran_to_vtk

    bdf_data_to_plot = chain(
        itervalues(bdf.nodes),
        itervalues(bdf.elements),
        itervalues(bdf.rigid_elements)
    )

    category_list = []

    for elem in bdf_data_to_plot:
        elem_type = elem.type
        cell_type, add_method, category = _nastran_to_vtk.get(elem_type, (None, None, None))

        if cell_type is None:
            continue

        cell_types.append(cell_type)
        elem_types.append(elem_type)

        eid = add_method(elem, cells, nid_dict)  # returns element/grid id

        eid_list.append(eid)
        eid_dict[eid] = cell_count

        category_list.append(categories[category])

        cell_count += 1

    cells = np.array(cells, dtype=np.int64)

    id_array = vtk.vtkIdTypeArray()
    id_array.SetVoidArray(cells, len(cells), 1)

    vtk_cells = vtk.vtkCellArray()
    vtk_cells.SetCells(cell_count, id_array)

    cell_types = np.array(cell_types, 'B')
    vtk_cell_types = numpy_to_vtk(cell_types)

    cell_locations = np.array([i for i in range(cell_count)])

    vtk_cell_locations = numpy_to_vtk(cell_locations, deep=1, array_type=vtk.VTK_ID_TYPE)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    ugrid.SetCells(vtk_cell_types, vtk_cell_locations, vtk_cells)
    
    vtk_cell_locations.SetName('index')
    ugrid.GetCellData().AddArray(vtk_cell_locations)
    
    _elem_types = vtk.vtkStringArray()
    _elem_types.SetName('element_type')
    _elem_types.SetNumberOfValues(len(elem_types))
    
    for i in range(len(elem_types)):
        _elem_types.SetValue(i, elem_types[i])

    ugrid.GetCellData().AddArray(_elem_types)

    _cat_arr = np.array(category_list, dtype=np.int64)
    _cat = vtk.vtkIdTypeArray()
    _cat.SetNumberOfValues(len(category_list))
    _cat.SetVoidArray(_cat_arr, len(_cat_arr), 1)
    _cat.SetName('category')

    ugrid.GetCellData().AddArray(_cat)

    id_array = numpy_to_vtk(np.array(eid_list), deep=1, array_type=vtk.VTK_ID_TYPE)
    
    # id_array = vtk.vtkIdTypeArray()
    # id_array.SetVoidArray(eid_list, len(eid_list), 1)
    id_array.SetName('element_id')
    
    ugrid.GetCellData().AddArray(id_array)

    copy = vtk.vtkUnstructuredGrid()
    copy.DeepCopy(ugrid)

    vtk_data = VTKData(copy, nid_list, nid_dict, eid_list, eid_dict)

    return vtk_data


categories = {
    'grid': 0,
    'element': 1,
    'mpc': 2,
    'force': 3,
    'disp': 4,
    'cord': 5
}


nastran_to_vtk = {
    'GRID': (vtk.VTK_VERTEX, _grid, 'grid'),
    'CQUAD4': (vtk.VTK_QUAD, _basic_elem, 'element'),
    'CTRIA3': (vtk.VTK_TRIANGLE, _basic_elem, 'element'),
    'CBEAM': (vtk.VTK_LINE, _basic_elem, 'element'),
    'CBUSH': (vtk.VTK_LINE, _basic_elem, 'element'),
    'CBAR': (vtk.VTK_LINE, _basic_elem, 'element'),
    'RBE2': (vtk.VTK_POLY_LINE, _rbe2, 'mpc'),
    'RBE3': (vtk.VTK_POLY_LINE, _rbe3, 'mpc')
}