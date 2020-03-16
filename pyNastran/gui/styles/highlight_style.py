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
import numpy as np
import vtk
#from vtk.util.numpy_support import vtk_to_numpy
#from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.gui.utils.vtk.vtk_utils import (
    extract_selection_node_from_grid_to_ugrid,
    find_point_id_closest_to_xyz,
    create_vtk_selection_node_by_point_ids,
    create_vtk_selection_node_by_cell_ids)

#vtkInteractorStyleRubberBandPick # sucks?
#class AreaPickStyle(vtk.vtkInteractorStyleDrawPolygon):  # not sure how to use this one...
class HighlightStyle(vtk.vtkInteractorStyleTrackballCamera):  # works
    """Highlights nodes & elements"""
    def __init__(self, parent=None, is_eids=True, is_nids=True, representation='wire',
                 name=None, callback=None, cleanup=True):
        """
        Creates the HighlightStyle instance

        Parameters
        ----------
        is_eids/is_nids : bool; default=True
            should elements/nodes be highlighted
        representation : str; default='wire'
            allowed = {'wire', 'points', 'surface'}
        name : str; default=None
            the name of the actor
        callback : function
            fill up a QLineEdit or some other custom action
        cleanup : bool; default=True
            should the actor be removed when the camera is moved

        """
        self.AddObserver("LeftButtonPressEvent", self._left_button_press_event)
        #self.AddObserver("LeftButtonReleaseEvent", self._left_button_release_event)
        #self.AddObserver("RightButtonPressEvent", self.right_button_press_event)
        self.parent = parent
        self.highlight_button = self.parent.actions['highlight']
        #self.area_pick_button = self.parent.actions['area_pick']
        #self.picker_points = []
        #self.parent.area_picker.SetRenderer(self.parent.rend)
        self.is_eids = is_eids
        self.is_nids = is_nids
        self.representation = representation
        assert is_eids or is_nids, 'is_eids=%r is_nids=%r, must not both be False' % (is_eids, is_nids)
        self.callback = callback
        self.cleanup = cleanup
        #self._pick_visible = False
        self.model_name = name
        assert self.model_name is not None
        self.actors = []

    def _left_button_press_event(self, obj, event):
        """
        gets the first point
        """
        self.OnLeftButtonDown()

        picker = self.parent.cell_picker
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        picker.Pick(pixel_x, pixel_y, 0, self.parent.rend)

        cell_id = picker.GetCellId()

        if cell_id < 0:
            return

        #icase = self.gui.icase_fringe
        #if icase is None:
            #return

        world_position = picker.GetPickPosition()

        grid = self.parent.get_grid_selected(self.model_name)
        cell_ids = [cell_id]
        point_ids = []
        if self.is_eids: # highlight_style = 'centroid
            actor = self._highlight_picker_cell(cell_id, grid)
        elif self.is_nids: # highlight_style = 'node'
            actor, point_ids = self._highlight_picker_node(cell_id, grid, world_position)
        else:
            raise RuntimeError('invalid highlight_style=%r' % self.highlight_style)


        #print('highlight_style  point_id=', point_id)
        #self.remove_actors()
        self.actors = [actor]
        self.parent.vtk_interactor.Render()

        if self.callback is not None:
            eids, nids = map_cell_point_to_model(self.parent.gui, cell_ids, point_ids,
                                                 model_name=self.model_name)
            self.callback(eids, nids, self.model_name)

        self.highlight_button.setChecked(False)

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

    def _highlight_picker_node(self, cell_id, grid, node_xyz, add_actor=True):
        """won't handle multiple cell_ids/node_xyz"""
        point_id = find_point_id_closest_to_xyz(grid, cell_id, node_xyz)
        selection_node = create_vtk_selection_node_by_point_ids([point_id])
        ugrid = extract_selection_node_from_grid_to_ugrid(grid, selection_node)
        actor = self.parent.create_highlighted_actor(
            ugrid, representation='points', add_actor=add_actor)
        return actor, [point_id]

    def _highlight_picker_cell(self, cell_ids, grid, add_actor=True):
        """won't handle multiple cell_ids/node_xyz"""
        selection_node = create_vtk_selection_node_by_cell_ids(cell_ids)
        ugrid = extract_selection_node_from_grid_to_ugrid(grid, selection_node)
        actor = self.parent.create_highlighted_actor(
            ugrid, representation='surface', add_actor=add_actor)
        return actor

def map_cell_point_to_model(gui, cell_ids, point_ids, model_name=None):
    eids = []
    nids = []
    if cell_ids:
        cell_array = np.asarray(cell_ids, dtype='int32')
        eids = gui.get_element_ids(model_name, cell_array)

    if point_ids:
        point_array = np.asarray(point_ids, dtype='int32')
        nids = gui.get_node_ids(model_name, point_array)
    return eids, nids
