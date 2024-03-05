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
from __future__ import annotations
from typing import Callable, Optional

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
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleRubberBandZoom
from vtkmodules.vtkRenderingCore import vtkActor, vtkRenderedAreaPicker

from pyNastran.gui.utils.vtk.gui_utils import add_actors_to_gui
from pyNastran.gui.styles.area_pick_utils import (
    get_actors_by_area_picker)

#if TYPE_CHECKING:
    #from pyNastran.gui.main_window import MainWindow

#class AreaPickStyle(vtkInteractorStyleRubberBandPick):
    #"""Custom Rubber Band Picker"""
    #def __init__(self, parent=None):
        #"""creates the AreaPickStyle instance"""
        #pass
        ##super(AreaPickStyle, self).__init__()

        ## for vtkInteractorStyleRubberBandZoom
        #self.AddObserver("LeftButtonPressEvent", self._left_button_press_event)
        #self.AddObserver("LeftButtonReleaseEvent", self._left_button_release_event)
        ##self.AddObserver("RightButtonPressEvent", self.right_button_press_event)
        #self.parent = parent
        #self.area_pick_button = self.parent.actions['area_pick']
        #self.picker_points = []
        #self.parent.area_picker.SetRenderer(self.parent.rend)

#vtkInteractorStyleRubberBandPick # sucks?
#class AreaPickStyle(vtkInteractorStyleDrawPolygon):  # not sure how to use this one...
class AreaPickStyle(vtkInteractorStyleRubberBandZoom):  # works
    """Picks nodes & elements with a visible box widget"""
    def __init__(self, parent=None,
                 is_eids: bool=True,
                 is_nids: bool=True,
                 representation: str='wire',
                 name=None,
                 callback: Optional[Callable]=None,
                 cleanup: bool=True):
        """creates the AreaPickStyle instance"""
        # for vtkInteractorStyleRubberBandZoom
        self.AddObserver("LeftButtonPressEvent", self._left_button_press_event)
        self.AddObserver("LeftButtonReleaseEvent", self._left_button_release_event)
        self.AddObserver("RightButtonPressEvent", self.right_button_press_event)
        self.parent = parent
        self.area_pick_button = self.parent.actions['area_pick']
        self.picker_points: list[int] = []
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
        self.actors: list[vtkActor] = []

    def _left_button_press_event(self, obj, event) -> None:
        """gets the first point"""
        #print('area_picker - left_button_press_event')
        self.OnLeftButtonDown()
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        self.picker_points.append((pixel_x, pixel_y))

    def _left_button_release_event(self, obj, event) -> None:
        """
        gets the second point and zooms

        TODO: doesn't handle panning of the camera to center the image
              with respect to the selected limits

        """
        #self.OnLeftButtonUp()
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        #selector = vtkVisibleCellSelector()

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

    def _pick_visible_ids(self,
                          xmin: float, ymin: float,
                          xmax: float, ymax: float) -> None:
        """Does an area pick of all the visible ids inside the box"""
        #vtkSelectVisiblePoints()
        #vcs = vtkVisibleCellSelector()
        area_picker = vtkRenderedAreaPicker()
        area_picker.AreaPick(xmin, ymin, xmax, ymax, self.parent.rend)
        #area_picker.Pick()

    def _pick_depth_ids(self,
                        xmin: float, ymin: float,
                        xmax: float, ymax: float) -> None:
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

    def remove_actors(self) -> None:
        if len(self.actors) == 0:
            return
        for actor in self.actors:
            self.parent.rend.RemoveActor(actor)
        self.actors = []

    def cleanup_callback(self, obj, event) -> None:
        """this is the cleanup step to remove the highlighted actor"""
        self.remove_actors()
        #self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
        self.parent.vtk_interactor.RemoveObserver(self.cleanup_observer)
        #cleanup_observer = None

    def right_button_press_event(self, obj, event) -> None:
        """cancels the button"""
        self.area_pick_button.setChecked(False)
        self.parent.setup_mouse_buttons(mode='default')
        self.parent.vtk_interactor.Render()
