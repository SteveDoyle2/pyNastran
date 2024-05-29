from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

#from pyNastran.gui.vtk_interface import
#from vtk import (
    #vtkInteractorStyleRubberBandZoom,
    #vtkSelection, vtkSelectionNode, vtkPlanes,
    #vtkIdFilter,
    #vtkExtractPoints,
    #vtkExtractSelectedFrustum,
    #vtkExtractSelection,
    #vtkRenderedAreaPicker,
#)
from vtkmodules.vtkCommonDataModel import (
    vtkCellData, vtkPointData, vtkSelection, vtkSelectionNode, vtkPlanes)
from vtkmodules.vtkFiltersPoints import vtkExtractPoints
from vtkmodules.vtkFiltersGeneral import vtkExtractSelectedFrustum
from vtkmodules.vtkFiltersExtraction import vtkExtractSelection
from vtkmodules.vtkRenderingCore import vtkActor

from pyNastran.gui.vtk_rendering_core import vtkActor
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid

from pyNastran.gui.vtk_util import vtk_to_numpy
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_unstructured_point_grid, numpy_to_vtk_points,
)
from pyNastran.gui.menus.highlight.vtk_utils import (
    create_highlighted_ugrids,
    create_highlighted_actors, get_ids_filter)

if TYPE_CHECKING:
    from vtkmodules.vtkCommonCore import vtkDataArray, vtkPoints, vtkIdTypeArray
    from vtkmodules.vtkRenderingCore import vtkAreaPicker
    from pyNastran.gui.main_window import MainWindow

def get_actors_by_area_picker(gui: MainWindow, area_picker: vtkAreaPicker,
                              model_name: str,
                              is_nids: bool=True, is_eids: bool=True,
                              representation: str='points',
                              add_actors: bool=False) -> tuple[list[vtkActor], list[int], list[int]]:
    """doesn't handle multiple actors yet..."""
    frustum: vtkPlanes = area_picker.GetFrustum()

    ugrid, eids, nids = get_depth_ids(
        gui, frustum, model_name=model_name,
        is_nids=is_nids, is_eids=is_eids,
        representation=representation)

    actors = []
    is_points = 'points' in representation
    is_wire = 'wire' in representation
    is_surface = 'surface' in representation

    if is_nids and is_points:
        actor = gui.create_highlighted_actor(
            ugrid, representation='points', add_actor=add_actors)
        actors.append(actor)

    if is_eids and is_wire:
        actor = gui.create_highlighted_actor(  # or surface
            ugrid, representation='wire', add_actor=add_actors)
        actors.append(actor)
    elif is_eids and is_surface:
        actor = gui.create_highlighted_actor(  # or surface
            ugrid, representation='surface', add_actor=add_actors)
        actors.append(actor)

    return actors, eids, nids


def get_depth_ids(gui: MainWindow,
                  frustum: vtkPlanes,
                  model_name: str='main',
                  is_nids: bool=True,
                  is_eids: bool=True,
                  representation: str='points') -> tuple[list[vtkActor], list[int], list[int]]:
    """
    Picks the nodes and/or elements.  Only one grid (e.g., the elements)
    is currently returned.
    """
    #len(gui.node_ids)
    #2675
    #len(gui.element_ids)
    #2657

    #Number Of Points: 2675
    #Number Of Cells: 2657
    grid = gui.get_grid_selected(model_name)

    #Number Of Points: 2675
    #Number Of Cells: 2657
    #grid = gui.get_grid_selected(model_name)

    #extract_ids = vtkExtractSelectedIds()
    #extract_ids.AddInputData(grid)

    ids = get_ids_filter(
        grid, idsname='Ids',
        is_nids=is_nids, is_eids=is_eids)

    # TODO: this ugrid includes all elements in the box, including the ones
    #       not visible...it's an issue.  The eids/nids  are right tho
    ugrid, ugrid_flipped = grid_ids_frustum_to_ugrid_ugrid_flipped(
        grid, ids, frustum)

    eids = None
    if is_eids:
        #Number Of Points: 46
        #Number Of Cells: 31
        cells: vtkCellData = ugrid.GetCellData()
        if cells is not None:
            ids = cells.GetArray('Ids')  # local ids, not global...
            if ids is not None:
                cell_ids = vtk_to_numpy(ids)
                assert len(cell_ids) == len(np.unique(cell_ids))
                eids = gui.get_element_ids(model_name, cell_ids)

    nids = None
    if is_nids:
        # TODO: the highlighted node ugrid is a problem for groups
        ugrid_points, nids = get_inside_point_ids(
            gui, ugrid, ugrid_flipped, model_name,
            representation=representation)
        ugrid = ugrid_points

    if is_eids and eids is not None and len(eids):
        # the highlighted element ugrid works for groups
        ugrid_node, ugrid_element = create_highlighted_ugrids(
            gui, grid,
            #all_nodes=None, nodes=None, set_node_scalars=True,
            all_elements=gui.element_ids, elements=eids, set_element_scalars=True)
        del ugrid_node
        ugrid = ugrid_element

    return ugrid, eids, nids


def get_inside_point_ids(gui, ugrid: vtkUnstructuredGrid,
                         ugrid_flipped: vtkUnstructuredGrid,
                         model_name: str,
                         representation: str='points') -> tuple[vtkUnstructuredGrid, list[int]]:
    """
    The points that are returned from the frustum, despite being
    defined as inside are not all inside.  The cells are correct
    though.  If you determine the cells outside the volume and the
    points associated with that, and boolean the two, you can find
    the points that are actually inside.

    In other words, ``points`` corresponds to the points inside the
    volume and barely outside.  ``point_ids_flipped`` corresponds to
    the points entirely outside the volume.

    Parameters
    ==========
    ugrid : vtkUnstructuredGrid()
        the "inside" grid
    ugrid_flipped : vtkUnstructuredGrid()
        the outside grid

    Returns
    =======
    ugrid : vtkUnstructuredGrid()
        an updated grid that has the correct points
    nids : (n, ) int ndarray
        the node_ids

    """
    nids = None
    points: vtkPointData = ugrid.GetPointData()
    if points is None:
        return ugrid, nids

    ids: vtkIdTypeArray = points.GetArray('Ids')
    if ids is None:
        return ugrid, nids

    # all points associated with the correctly selected cells are returned
    # but we get extra points for the cells that are inside and out
    point_ids = vtk_to_numpy(ids)

    # these are the points outside the box/frustum (and also include the bad point)
    points_flipped: vtkPointData = ugrid_flipped.GetPointData()
    ids_flipped: vtkIdTypeArray = points_flipped.GetArray('Ids')
    point_ids_flipped = vtk_to_numpy(ids_flipped)

    # setA - setB
    point_ids2 = np.setdiff1d(point_ids, point_ids_flipped)

    # convert form point_id to node_id
    nids2 = gui.get_node_ids(model_name, point_ids2)

    #------------------
    if representation == 'points':
        # we need to filter the nodes that were filtered by the
        # numpy setdiff1d, so we don't show extra points
        nids = gui.get_node_ids(model_name, point_ids)
        ugrid = create_filtered_point_ugrid_from_nids(ugrid, nids, nids2)

    nids = nids2
    return ugrid, nids


def grid_ids_frustum_to_ugrid_ugrid_flipped(grid: vtkUnstructuredGrid,
                                            ids,
                                            frustum: vtkPlanes):
    if 1:
        selected_frustum = vtkExtractSelectedFrustum()
        #selected_frustum.ShowBoundsOn()
        #selected_frustum.SetInsideOut(1)
        selected_frustum.SetFrustum(frustum)
        # PreserveTopologyOn: return an insidedness array
        # PreserveTopologyOff: return a ugrid
        selected_frustum.PreserveTopologyOff()
        #selected_frustum.PreserveTopologyOn()
        selected_frustum.SetInputConnection(ids.GetOutputPort())  # was grid?
        selected_frustum.Update()
        ugrid = selected_frustum.GetOutput()

        # we make a second frustum to remove extra points
        selected_frustum_flipped = vtkExtractSelectedFrustum()
        selected_frustum_flipped.SetInsideOut(1)
        selected_frustum_flipped.SetFrustum(frustum)
        selected_frustum_flipped.PreserveTopologyOff()
        selected_frustum_flipped.SetInputConnection(ids.GetOutputPort())  # was grid?
        selected_frustum_flipped.Update()
        ugrid_flipped = selected_frustum_flipped.GetOutput()
    else:  # pragma: no cover
        unused_extract_points = vtkExtractPoints()
        selection_node = vtkSelectionNode()
        selection = vtkSelection()
        #selection_node.SetContainingCellsOn()
        selection_node.Initialize()
        selection_node.SetFieldType(vtkSelectionNode.POINT)
        selection_node.SetContentType(vtkSelectionNode.INDICES)

        selection.AddNode(selection_node)

        extract_selection = vtkExtractSelection()
        extract_selection.SetInputData(0, grid)
        extract_selection.SetInputData(1, selection) # vtk 6+
        extract_selection.Update()

        ugrid = extract_selection.GetOutput()
        ugrid_flipped = None
    return ugrid, ugrid_flipped

def create_filtered_point_ugrid_from_nids(
    ugrid: vtkUnstructuredGrid,
    nids: np.ndarray,
    nids2: np.ndarray) -> vtkUnstructuredGrid:
    """
    We need to filter the nodes that were filtered by the
    numpy setdiff1d, so we don't show extra points

    """
    #unused_pointsu = ugrid.GetPoints()
    point_data: vtkPoints = ugrid.GetPoints()
    output_data: vtkDataArray = point_data.GetData()
    points_array = vtk_to_numpy(output_data)  # yeah!

    isort_nids = np.argsort(nids)
    nids = nids[isort_nids]
    inids = np.searchsorted(nids, nids2)

    points_array_sorted = points_array[isort_nids, :]
    point_array2 = points_array_sorted[inids, :]
    points2: vtkPoints = numpy_to_vtk_points(point_array2)

    npoints = len(nids2)
    ugrid = create_unstructured_point_grid(points2, npoints)
    return ugrid

#def create_filtered_point_ugrid_from_inids(ugrid: vtkUnstructuredGrid,
                                           #inids: np.ndarray) -> vtkUnstructuredGrid:
