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
from __future__ import print_function, division
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


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
    def __init__(self, parent=None, is_eids=True, is_nids=True, name=None, callback=None):
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
        assert is_eids or is_nids, 'is_eids=%r is_nids=%r, must not both be False' % (is_eids, is_nids)
        self.callback = callback
        self._pick_visible = False
        self.name = name
        assert name is not None

    def _left_button_press_event(self, obj, event):
        """
        gets the first point
        """
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
                if self._pick_visible:
                    self._pick_visible_ids(xmin, ymin, xmax, ymax)
                else:
                    self._pick_depth_ids(xmin, ymin, xmax, ymax)
            self.parent.vtk_interactor.Render()
        self.picker_points = []

    def _pick_visible_ids(self, xmin, ymin, xmax, ymax):
        """
        Does an area pick of all the visible ids inside the box
        """
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
        frustum = area_picker.GetFrustum() # vtkPlanes
        #frustum = create_box_frustum(xmin, ymin, xmax, ymax, self.parent.rend)

        grid = self.parent.get_grid(self.name)

        #extract_ids = vtk.vtkExtractSelectedIds()
        #extract_ids.AddInputData(grid)

        idsname = "Ids"
        ids = vtk.vtkIdFilter()
        if isinstance(grid, vtk.vtkUnstructuredGrid):
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

        #if not self.is_eids:
            #ids.CellIdsOff()
        #if not self.is_nids:
            #ids.PointIdsOff()
        #ids.FieldDataOn()
        ids.SetIdsArrayName(idsname)

        if 1:
            selected_frustum = vtk.vtkExtractSelectedFrustum()
            selected_frustum.SetFrustum(frustum)
            selected_frustum.PreserveTopologyOff() #  don't make an unstructured grid
            selected_frustum.SetInputConnection(ids.GetOutputPort())  # was grid?
            selected_frustum.Update()

            ugrid = selected_frustum.GetOutput()
        else:  # pragma: no cover
            extract_points = vtk.vtkExtractPoints()
            selection_node =vtk.vtkSelectionNode()
            selection = vtk.vtkSelection()
            selection_node.SetContainingCellsOn()
            selection_node.Initialize()
            selection_node.SetFieldType(vtkSelectionNode.POINT)
            selection_node.SetContentType(vtkSelectionNode.INDICES)

            selection.AddNode(selection_node)

            extract_selection = vtk.vtkExtractSelection()
            extract_selection.SetInputData(0, grid)
            extract_selection.SetInputData(1, selection) # vtk 6+
            extract_selection.Update()

            ugrid = extract_selection.GetOutput()

        eids = None
        nids = None

        msg = ''
        if self.is_eids:
            cells = ugrid.GetCellData()
            if cells is not None:
                ids = cells.GetArray('Ids')
                if ids is not None:
                    cell_ids = vtk_to_numpy(ids)
                    assert len(cell_ids) == len(np.unique(cell_ids))
                    eids = self.parent.get_element_ids(self.name, cell_ids)
        if self.is_nids:
            points = ugrid.GetPointData()
            if points is not None:
                ids = points.GetArray('Ids')
                if ids is not None:
                    point_ids = vtk_to_numpy(ids)
                    nids = self.parent.get_node_ids(self.name, point_ids)

        if self.callback is not None:
            self.callback(eids, nids, self.name)

        self.area_pick_button.setChecked(False)
        self.parent.setup_mouse_buttons(mode='default')

    def right_button_press_event(self, obj, event):
        """cancels the button"""
        self.area_pick_button.setChecked(False)
        self.parent.setup_mouse_buttons(mode='default')
        self.parent.vtk_interactor.Render()

def create_box_frustum(x0_, y0_, x1_, y1_, renderer):  # pragma: no cover

    if x0_ < x1_:
        x0 = x0_
        x1 = x1_
    else:
        x0 = x1_
        x1 = x0_

    if y0_ < y1_:
        y0 = y0_
        y1 = y1_
    else:
        y0 = y1_
        y1 = y0_

    if x0 == x1:
        x1 += 1.

    if y0 == y1:
        y1 += 1.

    verts = []

    renderer.SetDisplayPoint(x0, y0, 0)
    renderer.DisplayToWorld()
    verts.extend(renderer.GetWorldPoint()[:4])

    renderer.SetDisplayPoint(x0, y0, 1)
    renderer.DisplayToWorld()
    verts.extend(renderer.GetWorldPoint()[:4])

    renderer.SetDisplayPoint(x0, y1, 0)
    renderer.DisplayToWorld()
    verts.extend(renderer.GetWorldPoint()[:4])

    renderer.SetDisplayPoint(x0, y1, 1)
    renderer.DisplayToWorld()
    verts.extend(renderer.GetWorldPoint()[:4])

    renderer.SetDisplayPoint(x1, y0, 0)
    renderer.DisplayToWorld()
    verts.extend(renderer.GetWorldPoint()[:4])

    renderer.SetDisplayPoint(x1, y0, 1)
    renderer.DisplayToWorld()
    verts.extend(renderer.GetWorldPoint()[:4])

    renderer.SetDisplayPoint(x1, y1, 0)
    renderer.DisplayToWorld()
    verts.extend(renderer.GetWorldPoint()[:4])

    renderer.SetDisplayPoint(x1, y1, 1)
    renderer.DisplayToWorld()
    verts.extend(renderer.GetWorldPoint()[:4])

    extract_selected_frustum = vtk.vtkExtractSelectedFrustum()
    extract_selected_frustum.CreateFrustum(verts)

    return extract_selected_frustum.GetFrustum()
