from __future__ import annotations
from typing import Union, TYPE_CHECKING
import numpy as np

#from vtk import (
    #vtkLODActor,
    #vtkLabeledDataMapper,
    #vtkCellCenters, vtkIdFilter,
    #vtkUnstructuredGridGeometryFilter,
    #vtkVertexGlyphFilter,
#)

from vtkmodules.vtkCommonCore import vtkPoints, vtkDataArray
from vtkmodules.vtkCommonDataModel import vtkCellData, vtkPointData, vtkPolyData
from vtkmodules.vtkRenderingCore import (
    vtkActor, vtkDataSetMapper, vtkActor2D,
    vtkPolyDataMapper, vtkRenderer)
from vtkmodules.vtkRenderingLOD import vtkLODActor
from vtkmodules.vtkRenderingLabel import vtkLabeledDataMapper
from vtkmodules.vtkFiltersCore import vtkCellCenters, vtkIdFilter
from vtkmodules.vtkFiltersGeometry import vtkUnstructuredGridGeometryFilter
from vtkmodules.vtkFiltersGeneral import vtkVertexGlyphFilter

#from pyNastran.gui.vtk_rendering_core import vtkActor, vtkActor2D, vtkRenderer, vtkPolyDataMapper
from pyNastran.gui.vtk_util import vtk_to_numpy
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.utils.vtk.vtk_utils import set_vtk_id_filter_name
#from pyNastran.gui.menus.menu_utils import eval_float_from_string

from pyNastran.gui.gui_objects.settings import Settings
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_selection_node_by_cell_ids,
    extract_selection_node_from_grid_to_ugrid,
    create_unstructured_point_grid, numpy_to_vtk_points, numpy_to_vtk)
if TYPE_CHECKING:
    from pyNastran.gui.main_window import MainWindow


def get_ids_filter(grid: Union[vtkUnstructuredGrid, vtkPolyData],
                   idsname: str='Ids',
                   is_nids: bool=True,
                   is_eids: bool=True) -> vtkIdFilter:
    """
    get the vtkIdFilter associated with a grid and either
    nodes/elements or both

    """
    ids = vtkIdFilter()
    if isinstance(grid, vtkUnstructuredGrid):
        # this is typically what's called in the gui
        ids.SetInputData(grid)
    elif isinstance(grid, vtkPolyData):  # pragma: no cover
        # this doesn't work...
        cell_data: vtkCellData = grid.GetCellData()
        point_data: vtkPointData = grid.GetPointData()
        ids.SetCellIds(cell_data)
        ids.SetPointIds(point_data)
    else:  # pragma: no cover
        raise NotImplementedError(ids)

    #self.is_eids = False
    #print('is_eids=%s is_nids=%s' % (is_eids, is_nids))
    #ids.FieldDataOn()

    if is_nids:
        set_vtk_id_filter_name(ids, idsname, is_points=True)
        ids.PointIdsOn()
    else:
        ids.PointIdsOff()

    if is_eids:
        set_vtk_id_filter_name(ids, idsname, is_cells=True)
        ids.CellIdsOn()
    else:
        ids.CellIdsOff()
    return ids


def create_cell_label_actor(actor: vtkLODActor, rend: vtkRenderer,
                            label_size: float) -> vtkActor2D:
    mapper = actor.GetMapper()
    mygrid = mapper.GetInput()

    element_id_filter = get_ids_filter(
        mygrid, idsname='Ids_cells', is_nids=False, is_eids=True)
    element_id_filter.SetFieldData(1)
    element_id_filter.SetCellIds(0)
    element_id_filter.FieldDataOn()

    # Create labels for cells
    cell_centers = vtkCellCenters()
    cell_centers.SetInputConnection(element_id_filter.GetOutputPort())

    cell_mapper = vtkLabeledDataMapper()
    cell_mapper.SetInputConnection(cell_centers.GetOutputPort())
    cell_mapper.SetLabelModeToLabelScalars()

    label_actor = vtkActor2D()
    label_actor.SetMapper(cell_mapper)
    return label_actor


def create_node_label_actor(actor: vtkLODActor, rend: vtkRenderer,
                            label_size: float) -> vtkActor2D:
    mapper = actor.GetMapper()
    mygrid = mapper.GetInput()

    point_id_filter = get_ids_filter(
        mygrid, idsname='Ids_points', is_nids=True, is_eids=False)
    point_id_filter.SetFieldData(1)
    point_id_filter.SetPointIds(0)
    point_id_filter.FieldDataOn()

    label_actor = create_node_labels(
        point_id_filter, mygrid, rend, label_size=label_size)
    return label_actor

def create_node_labels(point_id_filter: vtkIdFilter,
                       grid: vtkUnstructuredGrid,
                       rend: vtkRenderer,
                       label_size: float=10.0) -> vtkActor2D:
    """creates the node labels"""
    # filter inner points, so only surface points will be available
    geo = vtkUnstructuredGridGeometryFilter()
    geo.SetInputData(grid)

    # points
    vertex_filter = vtkVertexGlyphFilter()
    vertex_filter.SetInputConnection(geo.GetOutputPort())
    vertex_filter.Update()

    points_mapper = vtkPolyDataMapper()
    points_mapper.SetInputConnection(vertex_filter.GetOutputPort())
    points_mapper.ScalarVisibilityOn()
    points_actor = vtkActor()
    points_actor.SetMapper(points_mapper)
    prop = points_actor.GetProperty()
    prop.SetPointSize(label_size)

    # point labels
    label_mapper = vtkLabeledDataMapper()
    label_mapper.SetInputConnection(point_id_filter.GetOutputPort())
    label_mapper.SetLabelModeToLabelFieldData()
    label_actor = vtkActor2D()
    label_actor.SetMapper(label_mapper)
    return label_actor

def create_highlighted_ugrids(gui, grid: vtkUnstructuredGrid,
                              all_nodes=None, nodes=None, set_node_scalars: bool=True,
                              all_elements=None, elements=None, set_element_scalars: bool=True,
                              ) -> tuple[Optional[vtkUnstructuredGrid],
                                         Optional[vtkUnstructuredGrid]]:
    """
    Creates nodes & element highlighted objects

    Parameters
    ----------
    all_nodes / all_elements : None / (n,) np.ndarray
       all the nodes/elements
    nodes / elements : None / (n,) np.ndarray
       the subset of nodes/elements
    set_node_scalars / set_element_scalars : bool; default=True
        sets the node/element result because ???



    """
    ugrids = []
    nnodes = 0
    nelements = 0
    if nodes is not None:
        nnodes = len(nodes)
        assert len(all_nodes) >= nnodes

    if elements is not None:
        nelements = len(elements)
        assert len(all_elements) >= nelements
    assert nnodes + nelements > 0

    ugrid_node = None
    ugrid_element = None
    if nnodes:
        point_ids = np.searchsorted(all_nodes, nodes)
        point_data: vtkPoints = grid.GetPoints()
        output_data: vtkDataArray = point_data.GetData()
        points_array = vtk_to_numpy(output_data)  # yeah!

        point_array2 = points_array[point_ids, :]
        points2: vtkPoints = numpy_to_vtk_points(point_array2)

        ugrid = create_unstructured_point_grid(points2, nnodes)
        if set_node_scalars:
            point_ids_array: vtkDataArray = numpy_to_vtk(nodes)
            point_data2: vtkPointData = ugrid.GetPointData()
            point_data2.SetScalars(point_ids_array)
        ugrid_node = ugrid

    if nelements:
        cell_ids = np.searchsorted(all_elements, elements)

        selection_node = create_vtk_selection_node_by_cell_ids(cell_ids)
        ugrid = extract_selection_node_from_grid_to_ugrid(grid, selection_node)
        if set_element_scalars:
            element_ids_array: vtkDataArray = numpy_to_vtk(elements)
            point_data: vtkPointData = ugrid.GetPointData()
            cell_data: vtkCellData = ugrid.GetCellData()
            point_data.SetScalars(None)
            cell_data.SetScalars(element_ids_array)
        ugrid_element = ugrid
    return ugrid_node, ugrid_element

def create_highlighted_actors(gui, grid: vtkUnstructuredGrid,
                              all_nodes=None, nodes=None, set_node_scalars: bool=True,
                              all_elements=None, elements=None, set_element_scalars: bool=True,
                              add_actors: bool=False) -> list[vtkLODActor]:
    actors = []
    ugrid_node, ugrid_element = create_highlighted_ugrids(
        gui, grid,
        all_nodes=all_nodes, nodes=nodes, set_node_scalars=set_node_scalars,
        all_elements=all_elements, elements=elements, set_element_scalars=set_element_scalars)

    if nodes is not None:
        actor = create_highlighted_actor(
            gui, ugrid_node, representation='points',
            add_actor=add_actors)
        actors.append(actor)

    if elements is not None:
        actor = create_highlighted_actor(
            gui, ugrid_element, representation='wire',
            add_actor=add_actors)
        actors.append(actor)
    return actors

def create_highlighted_actors_old(gui, grid: vtkUnstructuredGrid,
                                  all_nodes=None, nodes=None, set_node_scalars: bool=True,
                                  all_elements=None, elements=None, set_element_scalars: bool=True,
                                  add_actors: bool=False) -> list[vtkLODActor]:
    """
    Creates nodes & element highlighted objects

    Parameters
    ----------
    all_nodes / all_elements : None / (n,) np.ndarray
       all the nodes/elements
    nodes / elements : None / (n,) np.ndarray
       the subset of nodes/elements
    set_node_scalars / set_element_scalars : bool; default=True
        sets the node/element result because ???
    add_actors : bool; default=True
        adds the actors to the renderer

    """
    actors = []
    nnodes = 0
    nelements = 0
    if nodes is not None:
        nnodes = len(nodes)
        assert len(all_nodes) >= nnodes

    if elements is not None:
        nelements = len(elements)
        assert len(all_elements) >= nelements
    assert nnodes + nelements > 0

    if nnodes:
        point_ids = np.searchsorted(all_nodes, nodes)
        point_data: vtkPoints = grid.GetPoints()
        output_data = point_data.GetData()
        points_array = vtk_to_numpy(output_data)  # yeah!

        point_array2 = points_array[point_ids, :]
        points2: vtkPoints = numpy_to_vtk_points(point_array2)

        ugrid = create_unstructured_point_grid(points2, nnodes)
        if set_node_scalars:
            point_ids_array: vtkDataArray = numpy_to_vtk(nodes)
            point_data: vtkPointData = ugrid.GetPointData()
            point_data.SetScalars(point_ids_array)
        actor = create_highlighted_actor(
            gui, ugrid, representation='points',
            add_actor=add_actors)
        actors.append(actor)

    if nelements:
        cell_ids = np.searchsorted(all_elements, elements)

        selection_node = create_vtk_selection_node_by_cell_ids(cell_ids)
        ugrid = extract_selection_node_from_grid_to_ugrid(grid, selection_node)
        if set_element_scalars:
            element_ids_array: vtkDataArray = numpy_to_vtk(elements)
            point_data: vtkPointData = ugrid.GetPointData()
            cell_data: vtkCellData = ugrid.GetCellData()
            point_data.SetScalars(None)
            cell_data.SetScalars(element_ids_array)
        actor = create_highlighted_actor(gui, ugrid, representation='wire',
                                         add_actor=add_actors)
        actors.append(actor)
    return actors

def create_highlighted_actor(gui: MainWindow, ugrid: vtkUnstructuredGrid,
                             representation: str='wire',
                             add_actor: bool=True) -> vtkLODActor:
    """creates a highlighted actor given a vtkUnstructuredGrid"""
    actor = vtkLODActor()
    mapper = vtkDataSetMapper()
    mapper.SetInputData(ugrid)
    # don't use a single color; makes setting prop values work
    mapper.ScalarVisibilityOff()
    actor.SetMapper(mapper)

    settings: Settings = gui.settings
    prop = actor.GetProperty()
    prop.SetColor(settings.highlight_color)
    prop.SetOpacity(settings.highlight_opacity)
    if representation == 'surface':
        pass
    elif representation == 'points':
        prop.SetRepresentationToPoints()
        prop.RenderPointsAsSpheresOn()
        prop.SetLighting(False)
        #prop.SetInterpolationToFlat()
        prop.SetPointSize(settings.highlight_point_size)
    elif representation == 'wire':
        prop.SetRepresentationToWireframe()
        prop.SetLineWidth(settings.highlight_line_width)
    else:  # pragma: no cover
        raise NotImplementedError(f'representation={representation!r} and must be [points, wire, surface]')

    if add_actor:
        gui.rend.AddActor(actor)
    return actor
