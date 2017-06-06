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
from pyNastran.bdf.utils import write_patran_syntax_dict


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
    def __init__(self, parent=None, is_eids=True, is_nids=True, callback=None):
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
        self.callback = callback
        self._pick_visible = False

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
        pass


    def _pick_depth_ids(self, xmin, ymin, xmax, ymax):
        """
        Does an area pick of all the ids inside the box, even the ones
        behind the front elements
        """
        area_picker = self.parent.area_picker
        area_picker.Pick()

        area_picker.AreaPick(xmin, ymin, xmax, ymax, self.parent.rend)
        frustrum = area_picker.GetFrustum() # vtkPlanes
        grid = self.parent.grid

        #extract_ids = vtk.vtkExtractSelectedIds()
        #extract_ids.AddInputData(grid)

        idsname = "Ids"
        ids = vtk.vtkIdFilter()
        ids.SetInputData(grid)
        if self.is_eids:
            ids.CellIdsOn()
        if self.is_nids:
            ids.PointIdsOn()
        #ids.FieldDataOn()
        ids.SetIdsArrayName(idsname)

        selected_frustrum = vtk.vtkExtractSelectedFrustum()
        selected_frustrum.SetFrustum(frustrum)
        selected_frustrum.PreserveTopologyOff()
        selected_frustrum.SetInputConnection(ids.GetOutputPort())  # was grid?
        selected_frustrum.Update()

        ugrid = selected_frustrum.GetOutput()
        eids = None
        nids = None

        msg = ''
        if self.is_eids:
            cells = ugrid.GetCellData()
            if cells is not None:
                cell_ids = vtk_to_numpy(cells.GetArray('Ids'))
                assert len(cell_ids) == len(np.unique(cell_ids))
                eids = self.parent.element_ids[cell_ids]
        if self.is_nids:
            points = ugrid.GetPointData()
            if points is not None:
                point_ids = vtk_to_numpy(points.GetArray('Ids'))
                nids = self.parent.node_ids[point_ids]

        if self.callback is not None:
            self.callback(eids, nids)

        self.area_pick_button.setChecked(False)
        self.parent.setup_mouse_buttons(mode='default')

    def right_button_press_event(self, obj, event):
        """cancels the button"""
        self.area_pick_button.setChecked(False)
        self.parent.setup_mouse_buttons(mode='default')
        self.parent.vtk_interactor.Render()
