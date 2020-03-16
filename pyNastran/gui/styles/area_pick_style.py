"""
defines the AreaPick class

http://www.paraview.org/Wiki/Selection_Implementation_in_VTK_and_ParaView_III
http://ruby-vtk.rubyforge.org/svn/trunk/VTK/Rendering/Testing/Cxx/TestAreaSelections.cxx
http://vtk.1045678.n5.nabble.com/picking-objects-from-a-subset-of-a-grid-td5143877.html
http://www.vtk.org/Wiki/VTK/Examples/Cxx/Picking/HighlightSelectedPoints
http://www.vtk.org/doc/nightly/html/classvtkExtractSelectedFrustum.html
http://www.vtk.org/doc/nightly/html/classvtkUnstructuredGridAlgorithm.html
http://www.vtk.org/doc/nightly/html/classvtkExtractCells.html
http://www.vtk.org/Wiki/VTK/Examples/Cxx/Picking/HighlightSelection
http://public.kitware.com/pipermail/vtkusers/2012-January/072046.html
http://vtk.1045678.n5.nabble.com/Getting-the-original-cell-id-s-from-vtkExtractUnstructuredGrid-td1239667.html
"""
from typing import Tuple, List, Union
import numpy as np
import vtk
#from vtk.util import numpy_support
from vtk.util.numpy_support import vtk_to_numpy
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_unstructured_point_grid, numpy_to_vtk_points)
from pyNastran.gui.utils.vtk.gui_utils import add_actors_to_gui

#class AreaPickStyle(vtk.vtkInteractorStyleRubberBandPick):
    #"""Custom Rubber Band Picker"""
    #def __init__(self, parent=None):
        #"""creates the AreaPickStyle instance"""
        #pass
        ##super(AreaPickStyle, self).__init__()

        ## for vtk.vtkInteractorStyleRubberBandZoom
        #self.AddObserver("LeftButtonPressEvent", self._left_button_press_event)
        #self.AddObserver("LeftButtonReleaseEvent", self._left_button_release_event)
        ##self.AddObserver("RightButtonPressEvent", self.right_button_press_event)
        #self.parent = parent
        #self.area_pick_button = self.parent.actions['area_pick']
        #self.picker_points = []
        #self.parent.area_picker.SetRenderer(self.parent.rend)

#vtkInteractorStyleRubberBandPick # sucks?
#class AreaPickStyle(vtk.vtkInteractorStyleDrawPolygon):  # not sure how to use this one...
class AreaPickStyle(vtk.vtkInteractorStyleRubberBandZoom):  # works
    """Picks nodes & elements with a visible box widget"""
    def __init__(self, parent=None, is_eids=True, is_nids=True, representation='wire',
                 name=None, callback=None, cleanup=True):
        """creates the AreaPickStyle instance"""
        # for vtk.vtkInteractorStyleRubberBandZoom
        self.AddObserver("LeftButtonPressEvent", self._left_button_press_event)
        self.AddObserver("LeftButtonReleaseEvent", self._left_button_release_event)
        self.AddObserver("RightButtonPressEvent", self.right_button_press_event)
        self.parent = parent
        self.area_pick_button = self.parent.actions['area_pick']
        self.picker_points = []
        self.parent.area_picker.SetRenderer(self.parent.rend)
        self.is_eids = is_eids
        self.is_nids = is_nids
        self.representation = representation
        assert is_eids or is_nids, 'is_eids=%r is_nids=%r, must not both be False' % (is_eids, is_nids)
        self.callback = callback
        self.cleanup = cleanup
        self._pick_visible = False
        self.name = name
        assert name is not None
        self.actors = []

    def _left_button_press_event(self, obj, event):
        """gets the first point"""
        #print('area_picker - left_button_press_event')
        self.OnLeftButtonDown()
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        self.picker_points.append((pixel_x, pixel_y))

    def _left_button_release_event(self, obj, event):
        """
        gets the second point and zooms

        TODO: doesn't handle panning of the camera to center the image
              with respect to the selected limits

        """
        #self.OnLeftButtonUp()
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        #selector = vtk.vtkVisibleCellSelector()

        self.picker_points.append((pixel_x, pixel_y))

        #print(self.picker_points)
        if len(self.picker_points) == 2:
            p1x, p1y = self.picker_points[0]
            p2x, p2y = self.picker_points[1]
            self.picker_points = []
            xmin = min(p1x, p2x)
            ymin = min(p1y, p2y)
            xmax = max(p1x, p2x)
            ymax = max(p1y, p2y)
            #print(self.picker_points)
            #print('_area_pick_left_button_release', cell_id)

            dx = abs(p1x - p2x)
            dy = abs(p1y - p2y)
            self.picker_points = []
            if dx > 0 and dy > 0:
                #self.remove_actors()
                if self._pick_visible:
                    self._pick_visible_ids(xmin, ymin, xmax, ymax)
                else:
                    self._pick_depth_ids(xmin, ymin, xmax, ymax)
            self.parent.vtk_interactor.Render()
        self.picker_points = []

    def _pick_visible_ids(self, xmin, ymin, xmax, ymax):
        """Does an area pick of all the visible ids inside the box"""
        #vtk.vtkSelectVisiblePoints()
        #vcs = vtk.vtkVisibleCellSelector()
        area_picker = vtk.vtkRenderedAreaPicker()
        area_picker.AreaPick(xmin, ymin, xmax, ymax, self.parent.rend)
        #area_picker.Pick()

    def _pick_depth_ids(self, xmin, ymin, xmax, ymax):
        """
        Does an area pick of all the ids inside the box, even the ones
        behind the front elements
        """
        area_picker = self.parent.area_picker
        #area_picker.Pick()  # double pick?

        area_picker.AreaPick(xmin, ymin, xmax, ymax, self.parent.rend)

        gui = self.parent
        actors, eids, nids = get_actors_by_area_picker(
            gui, area_picker, self.name,
            is_nids=self.is_nids, is_eids=self.is_eids,
            representation=self.representation, add_actors=False)

        self.actors = actors
        add_actors_to_gui(gui, actors, render=False)

        if self.callback is not None:
            self.callback(eids, nids, self.name)

        self.area_pick_button.setChecked(False)

        # TODO: it would be nice if you could do a rotation without
        #       destroying the highlighted actor
        if self.cleanup:
            self.cleanup_observer = self.parent.setup_mouse_buttons(
                mode='default', left_button_down_cleanup=self.cleanup_callback)

    def remove_actors(self):
        if len(self.actors) == 0:
            return
        for actor in self.actors:
            self.parent.rend.RemoveActor(actor)
        self.actors = []

    def cleanup_callback(self, obj, event):
        """this is the cleanup step to remove the highlighted actor"""
        self.remove_actors()
        #self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
        self.parent.vtk_interactor.RemoveObserver(self.cleanup_observer)
        #cleanup_observer = None

    def right_button_press_event(self, obj, event):
        """cancels the button"""
        self.area_pick_button.setChecked(False)
        self.parent.setup_mouse_buttons(mode='default')
        self.parent.vtk_interactor.Render()


def get_actors_by_area_picker(gui, area_picker,
                              model_name: str,
                              is_nids: bool=True, is_eids: bool=True,
                              representation: str='points',
                              add_actors: bool=False) -> Tuple[List[vtk.vtkActor], List[int], List[int]]:
    """doesn't handle multiple actors yet..."""
    frustum = area_picker.GetFrustum() # vtkPlanes

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


def get_depth_ids(gui, frustum: vtk.vtkPlanes, model_name: str='main',
                  is_nids: bool=True, is_eids: bool=True,
                  representation: str='points') -> Tuple[List[vtk.vtkActor], List[int], List[int]]:
    """
    Picks the nodes and/or elements.  Only one grid (e.g., the elements)
    is currently returned.
    """
    grid = gui.get_grid(model_name)

    #extract_ids = vtk.vtkExtractSelectedIds()
    #extract_ids.AddInputData(grid)

    ids = get_ids_filter(
        grid, idsname='Ids',
        is_nids=is_nids, is_eids=is_eids)
    ugrid, ugrid_flipped = grid_ids_frustum_to_ugrid_ugrid_flipped(
        grid, ids, frustum)

    eids = None
    if is_eids:
        cells = ugrid.GetCellData()
        if cells is not None:
            ids = cells.GetArray('Ids')
            if ids is not None:
                cell_ids = vtk_to_numpy(ids)
                assert len(cell_ids) == len(np.unique(cell_ids))
                eids = gui.get_element_ids(model_name, cell_ids)

    nids = None
    if is_nids:
        ugrid_points, nids = get_inside_point_ids(
            gui, ugrid, ugrid_flipped, model_name,
            representation=representation)
        ugrid = ugrid_points

    return ugrid, eids, nids


def get_inside_point_ids(gui, ugrid: vtk.vtkUnstructuredGrid,
                         ugrid_flipped: vtk.vtkUnstructuredGrid,
                         model_name: str,
                         representation: str='points') -> Tuple[vtk.vtkUnstructuredGrid, List[int]]:
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
    ugrid : vtk.vtkUnstructuredGrid()
        the "inside" grid
    ugrid_flipped : vtk.vtkUnstructuredGrid()
        the outside grid

    Returns
    =======
    ugrid : vtk.vtkUnstructuredGrid()
        an updated grid that has the correct points
    nids : (n, ) int ndarray
        the node_ids

    """
    nids = None
    points = ugrid.GetPointData()
    if points is None:
        return ugrid, nids

    ids = points.GetArray('Ids')
    if ids is None:
        return  ugrid, nids

    # all points associated with the correctly selected cells are returned
    # but we get extra points for the cells that are inside and out
    point_ids = vtk_to_numpy(ids)
    nids = gui.get_node_ids(model_name, point_ids)

    # these are the points outside the box/frustum (and also include the bad point)
    points_flipped = ugrid_flipped.GetPointData()
    ids_flipped = points_flipped.GetArray('Ids')
    point_ids_flipped = vtk_to_numpy(ids_flipped)
    nids_flipped = gui.get_node_ids(model_name, point_ids_flipped)
    #nids = gui.gui.get_reverse_node_ids(model_name, point_ids_flipped)

    # setA - setB
    nids2 = np.setdiff1d(nids, nids_flipped, assume_unique=True)

    #narrays = points.GetNumberOfArrays()
    #for iarray in range(narrays):
        #name = points.GetArrayName(iarray)
        #print('iarray=%s name=%r' % (iarray, name))

    #------------------
    if representation == 'points':
        # we need to filter the nodes that were filtered by the
        # numpy setdiff1d, so we don't show extra points
        ugrid = create_filtered_point_ugrid(ugrid, nids, nids2)

    nids = nids2
    return ugrid, nids


def get_ids_filter(grid: Union[vtk.vtkUnstructuredGrid, vtk.vtkPolyData],
                   idsname: str='Ids',
                   is_nids: bool=True, is_eids: bool=True) -> vtk.vtkIdFilter:
    """
    get the vtkIdFilter associated with a grid and either
    nodes/elements or both

    """
    ids = vtk.vtkIdFilter()
    if isinstance(grid, vtk.vtkUnstructuredGrid):
        # this is typically what's called in the gui
        ids.SetInputData(grid)
    elif isinstance(grid, vtk.vtkPolyData):  # pragma: no cover
        # this doesn't work...
        ids.SetCellIds(grid.GetCellData())
        ids.SetPointIds(grid.GetPointData())
    else:
        raise NotImplementedError(ids)

    #self.is_eids = False
    ids.CellIdsOn()
    ids.PointIdsOn()

    #print('is_eids=%s is_nids=%s' % (is_eids, is_nids))
    if not is_eids:
        ids.CellIdsOff()
    if not is_nids:
        ids.PointIdsOff()
    #ids.FieldDataOn()
    ids.SetIdsArrayName(idsname)
    return ids

def grid_ids_frustum_to_ugrid_ugrid_flipped(grid, ids, frustum):
    if 1:
        selected_frustum = vtk.vtkExtractSelectedFrustum()
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
        selected_frustum_flipped = vtk.vtkExtractSelectedFrustum()
        selected_frustum_flipped.SetInsideOut(1)
        selected_frustum_flipped.SetFrustum(frustum)
        selected_frustum_flipped.PreserveTopologyOff()
        selected_frustum_flipped.SetInputConnection(ids.GetOutputPort())  # was grid?
        selected_frustum_flipped.Update()
        ugrid_flipped = selected_frustum_flipped.GetOutput()
    else:  # pragma: no cover
        unused_extract_points = vtk.vtkExtractPoints()
        selection_node = vtk.vtkSelectionNode()
        selection = vtk.vtkSelection()
        #selection_node.SetContainingCellsOn()
        selection_node.Initialize()
        selection_node.SetFieldType(vtk.vtkSelectionNode.POINT)
        selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)

        selection.AddNode(selection_node)

        extract_selection = vtk.vtkExtractSelection()
        extract_selection.SetInputData(0, grid)
        extract_selection.SetInputData(1, selection) # vtk 6+
        extract_selection.Update()

        ugrid = extract_selection.GetOutput()
        ugrid_flipped = None
    return ugrid, ugrid_flipped

def create_filtered_point_ugrid(ugrid: vtk.vtkUnstructuredGrid,
                                nids, nids2) -> vtk.vtkUnstructuredGrid:
    """
    We need to filter the nodes that were filtered by the
    numpy setdiff1d, so we don't show extra points

    """
    #unused_pointsu = ugrid.GetPoints()
    output_data = ugrid.GetPoints().GetData()
    points_array = vtk_to_numpy(output_data)  # yeah!

    isort_nids = np.argsort(nids)
    nids = nids[isort_nids]
    inids = np.searchsorted(nids, nids2)

    points_array_sorted = points_array[isort_nids, :]
    point_array2 = points_array_sorted[inids, :]
    points2 = numpy_to_vtk_points(point_array2)

    npoints = len(nids2)
    ugrid = create_unstructured_point_grid(points2, npoints)
    return ugrid
