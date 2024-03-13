from __future__ import annotations
from typing import Optional, Callable, TYPE_CHECKING
import numpy as np

#from vtk import (
    #vtkInteractorStyleDrawPolygon,
    #vtkLODActor,
    #vtkSelection,
    #vtkExtractSelection,
#)
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleDrawPolygon
from vtkmodules.vtkRenderingLOD import vtkLODActor
from vtkmodules.vtkCommonDataModel import vtkSelection
from vtkmodules.vtkFiltersExtraction import vtkExtractSelection

from pyNastran.bdf.utils import write_patran_syntax_dict

from pyNastran.gui.vtk_common_core import vtkIdTypeArray
from pyNastran.gui.vtk_rendering_core import vtkRenderer
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid, vtkSelectionNode

from pyNastran.gui.gui_objects.settings import Settings

from pyNastran.gui.menus.highlight.vtk_utils import create_highlighted_actor

from pyNastran.gui.styles.highlight_style import HighlightStyle
from pyNastran.gui.styles.area_pick_style import AreaPickStyle
from pyNastran.gui.styles.zoom_style import ZoomStyle
#from pyNastran.gui.styles.probe_style import ProbeResultStyle
from pyNastran.gui.styles.rotation_center_style import RotationCenterStyle
from pyNastran.gui.styles.trackball_style_camera import TrackballStyleCamera
from pyNastran.gui.utils.vtk.vtk_utils import (
        find_point_id_closest_to_xyz, create_vtk_selection_node_by_cell_ids)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui import MainWindow


class MouseActions:
    def __init__(self, gui: MainWindow):
        self.gui = gui
        #self._camera_mode = None
        self._camera_mode = 'default'
        self.pick_state = 'node/centroid'
        self._zoom = []
        self._picker_points = []
        self._measure_distance_pick_points = []
        self.revert = None
        self.style = None
        self.cleanup_observer = None
        self.left_button_down_cleanup = None
        self._depress_when_done = []

    def setup_mouse_buttons(self,
                            mode: str='',
                            revert: bool=False,
                            left_button_down: Optional[Callable]=None,
                            left_button_up: Optional[Callable]=None,
                            right_button_down: Optional[Callable]=None,
                            end_pick=None,
                            style=None,
                            force: bool=False,
                            left_button_down_cleanup: Optional[Callable]=None):
        """
        Remaps the mouse buttons temporarily

        Parameters
        ----------
        mode : str
            lets you know what kind of mapping this is
        revert : bool; default=False
            does the button revert when it's finished
        left_button_down : function; default=None
            the callback function
            None: depends on the mode
        left_button_up : function; default=None
            the callback function
            None: depends on the mode
        right_button_down : function; default=None
            the callback function
            None: depends on the mode
        left_button_down_cleanup : function; default=None
            the callback function
            None: depends on the mode
        style : vtkInteractorStyle (default=None)
            a custom vtkInteractorStyle
            None -> keep the same style, but overwrite the left mouse button
        force : bool; default=False
            override the mode=camera_mode check

        """
        assert isinstance(mode, str), mode
        assert revert in [True, False], revert

        #print('setup_mouse_buttons mode=%r _camera_mode=%r' % (mode, self._camera_mode))
        if mode == self._camera_mode and not force:
            #print('auto return from set mouse mode')
            return
        self._camera_mode = mode

        if mode is None:
            # same as default
            #print('auto return 2 from set mouse mode')
            return
        elif mode == 'default':
            #print('set mouse mode as default')

            # standard rotation
            # Disable default left mouse click function (Rotate)
            self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
            self.vtk_interactor.RemoveObservers('EndPickEvent')
            self.vtk_interactor.AddObserver('EndPickEvent', self._probe_picker)
            # there should be a cleaner way to revert the trackball Rotate command
            # it apparently requires an (obj, event) argument instead of a void...
            self.set_style_as_trackball()

            # the more correct-ish way to reset the 'LeftButtonPressEvent' to Rotate
            # that doesn't work...
            #
            # Re-assign left mouse click event to custom function (Point Picker)
            #self.vtk_interactor.AddObserver('LeftButtonPressEvent', self.style.Rotate)

        elif mode == 'measure_distance': # 'rotation_center',
            # hackish b/c the default setting is so bad
            self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
            self.vtk_interactor.AddObserver('LeftButtonPressEvent', left_button_down)

            self.vtk_interactor.RemoveObservers('EndPickEvent')
            self.vtk_interactor.AddObserver('EndPickEvent', left_button_down)
        elif mode in {'probe_result', 'probe_result_all', 'highlight_cell', 'highlight_node'}:
            # hackish b/c the default setting is so bad
            self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
            self.vtk_interactor.AddObserver('LeftButtonPressEvent', left_button_down)

            self.vtk_interactor.RemoveObservers('EndPickEvent')
            self.vtk_interactor.AddObserver('EndPickEvent', left_button_down)
            #self.vtk_interactor.AddObserver('LeftButtonPressEvent', func, 1) # on press down
            #self.vtk_interactor.AddObserver('LeftButtonPressEvent', func, -1) # on button up

        elif mode == 'zoom':
            assert style is not None, style
            self.vtk_interactor.SetInteractorStyle(style)

            # on press down
            self.vtk_interactor.AddObserver('LeftButtonPressEvent', left_button_down)

            # on button up
            self.vtk_interactor.AddObserver('LeftButtonReleaseEvent', left_button_up, -1)
            if right_button_down:
                self.vtk_interactor.AddObserver('RightButtonPressEvent', right_button_down)


        #elif mode == 'node_pick':
            #self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
            #self.vtk_interactor.AddObserver('LeftButtonPressEvent', self.on_node_pick_event)
        #elif mode == 'cell_pick':
            #self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
            #self.vtk_interactor.AddObserver('LeftButtonPressEvent', self.on_cell_pick_event)
        elif mode == 'cell_pick':
            #print('set mouse mode as cell_pick')
            self.vtk_interactor.SetPicker(self.cell_picker)
        elif mode == 'node_pick':
            #print('set mouse mode as node_pick')
            self.vtk_interactor.SetPicker(self.node_picker)
        elif mode == 'style':
            self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
            self.vtk_interactor.RemoveObservers('RightButtonPressEvent')
            self.vtk_interactor.SetInteractorStyle(style)

        #elif mode == 'area_cell_pick':
            #self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
            #self.vtk_interactor.AddObserver('LeftButtonPressEvent',
                                            #self.on_area_cell_pick_event)
        #elif mode == 'area_node_pick':
            #self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
            #self.vtk_interactor.AddObserver('LeftButtonPressEvent',
                                            #self.on_area_cell_pick_event)
        #elif mode == 'polygon_cell_pick':
            #self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
            #self.vtk_interactor.AddObserver('LeftButtonPressEvent',
                                            #self.on_polygon_cell_pick_event)
        #elif mode == 'polygon_node_pick':
            #self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
            #self.vtk_interactor.AddObserver('LeftButtonPressEvent',
                                            #self.on_polygon_cell_pick_event)

        #elif mode == 'pan':
            #pass
        else:
            raise NotImplementedError('camera_mode = %r' % self._camera_mode)

        if left_button_down_cleanup:
            self.left_button_down_cleanup = left_button_down_cleanup
            self.cleanup_observer = self.vtk_interactor.AddObserver(
                'LeftButtonPressEvent', left_button_down_cleanup)
        #if left_button_up_cleanup:
            #self.left_button_down_cleanup = left_button_down_cleanup
            #self.cleanup_observer = self.vtk_interactor.AddObserver(
                #'LeftButtonPressEvent', left_button_down_cleanup)

        self.revert = revert
        return self.cleanup_observer

    def cleanup(self) -> None:
        """
        this cleanup gets called when the user immediately goes from an
        area pick into another area pick without left clicking.
        """
        self.left_button_down_cleanup(None, None)
        self.vtk_interactor.RemoveObserver(self.cleanup_observer)
        self.cleanup_observer = None

    def revert_pressed(self, active_name: str) -> None:
        if self.cleanup_observer is not None:
            self.cleanup()

        if active_name != 'probe_result':
            probe_button = self.actions['probe_result']
            is_checked = probe_button.isChecked()
            if is_checked:  # revert probe_result
                probe_button.setChecked(False)
                self.setup_mouse_buttons(mode='default')
                return
        if active_name != 'probe_result_all':
            probe_button = self.actions['probe_result_all']
            is_checked = probe_button.isChecked()
            if is_checked:  # revert probe_result_all
                probe_button.setChecked(False)
                self.setup_mouse_buttons(mode='default')
                return

        if active_name != 'highlight_cell':
            highlight_cell = self.actions['highlight_cell']
            is_checked = highlight_cell.isChecked()
            if is_checked:  # revert highlight_cell
                highlight_cell.setChecked(False)
                self.setup_mouse_buttons(mode='default')
                return

        if active_name != 'rotation_center':
            rotation_button = self.actions['rotation_center']
            is_checked = rotation_button.isChecked()
            if is_checked:  # revert rotation_center
                rotation_button.setChecked(False)
                self.setup_mouse_buttons(mode='default')
                return

        if active_name != 'measure_distance':
            measure_distance_button = self.actions['measure_distance']
            is_checked = measure_distance_button.isChecked()
            if is_checked:
                # revert on_measure_distance
                measure_distance_button.setChecked(False)
                self._measure_distance_pick_points = []
                self.setup_mouse_buttons(mode='default')
                return

        if active_name != 'zoom':
            zoom_button = self.actions['zoom']
            is_checked = zoom_button.isChecked()
            if is_checked:
                # revert on_measure_distance
                zoom_button.setChecked(False)
                self._zoom = []
                self.setup_mouse_buttons(mode='default')
                return

    def set_style_as_trackball(self) -> None:
        """sets the default rotation style"""
        #self._simulate_key_press('t') # change mouse style to trackball
        self.style = TrackballStyleCamera(self.vtk_interactor, self)
        self.vtk_interactor.SetInteractorStyle(self.style)

    def on_highlight_node(self) -> None:
        """highlights a single node (when it's added)"""
        #self.on_highlight(is_eids=True, is_nids=True, representation='surface',
                          #name=None, callback=None, force=True)

        self.highlight_style = 'node'
        self.revert_pressed('highlight_node')
        is_checked = self.actions['highlight_node'].isChecked()
        if not is_checked:
            # revert probe_result
            self.setup_mouse_buttons(mode='default')
            return
        self._depress_when_done = ['highlight_node']
        self.setup_mouse_buttons('highlight_node', left_button_down=self._highlight_picker,
                                 revert=True)

    def on_highlight_cell(self) -> None:
        """highlights a single cell"""
        #self.on_highlight(is_eids=True, is_nids=True, representation='surface',
                          #name=None, callback=None, force=True)
        self.highlight_style = 'centroid'
        self.revert_pressed('highlight_cell')
        is_checked = self.actions['highlight_cell'].isChecked()
        if not is_checked:
            # revert probe_result
            self.setup_mouse_buttons(mode='default')
            return
        self._depress_when_done = ['highlight_cell']
        self.setup_mouse_buttons('highlight_cell', left_button_down=self._highlight_picker,
                                 revert=True)

    def on_probe_result(self) -> None:
        #print('on_probe_result')
        self.revert_pressed('probe_result')
        is_checked = self.actions['probe_result'].isChecked()
        if not is_checked:
            # revert probe_result
            self.setup_mouse_buttons(mode='default')
            return
        self.setup_mouse_buttons('probe_result', left_button_down=self._probe_picker)

        #style = ProbeResultStyle(parent=self)
        #self.vtk_interactor.SetInteractorStyle(style)

    def on_quick_probe_result(self) -> None:
        """is this used?"""
        #print('on_quick_probe_result')
        self.revert_pressed('probe_result')
        unused_is_checked = self.actions['probe_result'].isChecked()
        self.setup_mouse_buttons('probe_result',
                                 left_button_down=self._probe_picker, revert=True)

    def on_probe_result_all(self) -> None:
        """called by button press, then click till done"""
        #print('on_probe_result_all')
        self.revert_pressed('probe_result_all')
        is_checked = self.actions['probe_result_all'].isChecked()
        if not is_checked:
            # revert probe_result
            self.setup_mouse_buttons(mode='default')
            return
        self.setup_mouse_buttons('probe_result_all', left_button_down=self._probe_picker_all)

    def on_quick_probe_result_all(self) -> None:
        """click and THEN press a"""
        #print('on_quick_probe_result_all')
        self.revert_pressed('probe_result_all')
        unused_is_checked = self.actions['probe_result_all'].isChecked()
        self.setup_mouse_buttons('probe_result_all',
                                 left_button_down=self._probe_picker_all, revert=True)

    def on_area_pick_callback(self, eids: list[int],
                              nids: list[int],
                              name: str) -> None:
        """prints the message when area_pick succeeds"""
        msg = ''
        if eids is not None and len(eids):
            msg += write_patran_syntax_dict({'Elem' : eids})
        if nids is not None and len(nids):
            msg += '\n' + write_patran_syntax_dict({'Node' : nids})
        if msg:
            self.gui.log_info('\n%s' % msg.lstrip())

    def on_area_pick(self, is_eids: bool=True,
                     is_nids: bool=True,
                     representation: str='wire+points',
                     name=None,
                     callback: Optional[Callable]=None,
                     cleanup: bool=True,
                     force: bool=False) -> Optional[AreaPickStyle]:
        """Creates a box picker"""
        style = None
        if name is None:
            name = self.gui.name
        self.revert_pressed('area_pick')

        is_checked = self.actions['area_pick'].isChecked()
        if not is_checked:
            # revert area_pick
            self.setup_mouse_buttons(mode='default')
            if not force:
                return style

        #self.gui.log_info('on_area_pick')
        self._picker_points = []

        if callback is None:
            callback = self.on_area_pick_callback
        style = AreaPickStyle(parent=self, is_eids=is_eids, is_nids=is_nids,
                              representation=representation,
                              name=name, callback=callback, cleanup=cleanup)
        self.setup_mouse_buttons(mode='style', revert=True,
                                 style=style) #, style_name='area_pick'
        return style


    def on_highlight(self, is_eids: bool=True,
                     is_nids: bool=True,
                     representation: str='wire',
                     name: Optional[str]=None,
                     callback: Optional[Callable]=None,
                     cleanup: bool=True,
                     force: bool=False) -> HighlightStyle:
        """
        Selects a single point/cell

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
        force : bool; default=False
            handles gui interaction; don't mess with this

        """
        if name is None:
            name = self.gui.name
        self.revert_pressed('highlight')

        is_checked = self.actions['highlight'].isChecked()
        if not is_checked:
            # revert probe_result
            self.setup_mouse_buttons(mode='default')
            if not force:
                return
        style = HighlightStyle(parent=self, is_eids=is_eids, is_nids=is_nids,
                               representation=representation,
                               name=name, callback=callback, cleanup=cleanup)
        self.setup_mouse_buttons(mode='style', revert=True,
                                 style=style)
        return style

    def on_area_pick_not_square(self) -> None:
        self.revert_pressed('area_pick')
        is_checked = self.actions['area_pick'].isChecked()
        if not is_checked:
            # revert area_pick
            self.setup_mouse_buttons(mode='default')
            return

        self.gui.log_info('on_area_pick')
        self.vtk_interactor.SetPicker(self.area_picker)

        def _area_picker_up(*args):
            pass
        style = vtkInteractorStyleDrawPolygon()
        self.setup_mouse_buttons('area_pick',
                                 #left_button_down=self._area_picker,
                                 left_button_up=_area_picker_up,
                                 #end_pick=self._area_picker_up,
                                 style=style)
        #self.area_picker = vtkAreaPicker()  # vtkRenderedAreaPicker?
        #self.rubber_band_style = vtkInteractorStyleRubberBandPick()
        #vtkInteractorStyleRubberBand2D
        #vtkInteractorStyleRubberBand3D
        #vtkInteractorStyleRubberBandZoom
        #vtkInteractorStyleAreaSelectHover
        #vtkInteractorStyleDrawPolygon

    def on_zoom(self) -> None:
        """creates a Rubber Band Zoom"""
        #self.revert_pressed('zoom')
        is_checked = self.actions['zoom'].isChecked()
        if not is_checked:
            # revert zoom
            self.setup_mouse_buttons(mode='default')
            return
        style = ZoomStyle(parent=self)
        self.setup_mouse_buttons(mode='style', revert=True, style=style)
        #self.vtk_interactor.SetInteractorStyle(style)

    def on_rotation_center(self) -> None:
        """
        http://osdir.com/ml/lib.vtk.user/2002-09/msg00079.html
        """
        self.revert_pressed('rotation_center')
        is_checked = self.actions['rotation_center'].isChecked()
        if not is_checked:
            # revert on_rotation_center
            self.setup_mouse_buttons(mode='default')
            return

        style = RotationCenterStyle(parent=self)
        self.setup_mouse_buttons('style', revert=True, style=style)

    def _get_closest_node_xyz(self,
                              cell_id: int,
                              world_position: np.ndarray) -> np.ndarray:
        unused_duplicate_key = None
        out = self.gui.get_result_by_xyz_cell_id(world_position, cell_id)
        if out is None:
            ## TODO: I think this happens when you do a focal point
            print(f'MouseActions._get_closest_node_xyz bug')
            return
        (result_name, unused_result_value, unused_node_id, xyz) = out
        assert self.gui.icase in self.gui.label_actors, result_name
        assert not isinstance(xyz, int), xyz
        return xyz

    def on_measure_distance(self) -> None:
        self.revert_pressed('measure_distance')
        measure_distance_button = self.actions['measure_distance']
        is_checked = measure_distance_button.isChecked()
        if not is_checked:
            # revert on_measure_distance
            self._measure_distance_pick_points = []
            self.setup_mouse_buttons(mode='default')
            return
        self._measure_distance_pick_points = []
        self.setup_mouse_buttons('measure_distance',
                                 left_button_down=self._measure_distance_picker)

    def _measure_distance_picker(self, unused_obj, unused_event) -> None:
        """picking is based on the two closest nodes"""
        picker = self.cell_picker
        pixel_x, pixel_y = self.vtk_interactor.GetEventPosition()
        picker.Pick(pixel_x, pixel_y, 0, self.rend)

        cell_id = picker.GetCellId()
        #print('_measure_distance_picker', cell_id)

        if cell_id < 0:
            #self.picker_textActor.VisibilityOff()
            pass
        else:
            world_position = picker.GetPickPosition()
            closest_point = self._get_closest_node_xyz(cell_id, world_position)

            if len(self._measure_distance_pick_points) == 0:
                self._measure_distance_pick_points.append(closest_point)
                self.gui.log_info('point1 = %s' % str(closest_point))
            else:
                self.gui.log_info('point2 = %s' % str(closest_point))
                p1 = self._measure_distance_pick_points[0]
                dxyz = closest_point - p1
                mag = np.linalg.norm(dxyz)

                self._measure_distance_pick_points = []
                self.gui.log_info('Node-Node: dxyz=%s mag=%s' % (str(dxyz), str(mag)))

                measure_distance_button = self.actions['measure_distance']
                measure_distance_button.setChecked(False)
                self.setup_mouse_buttons(mode='default')

    def _highlight_picker_node(self, cell_id: int, grid, node_xyz) -> vtkLODActor:
        """won't handle multiple cell_ids/node_xyz"""
        point_id = find_point_id_closest_to_xyz(grid, cell_id, node_xyz)

        ids = vtkIdTypeArray()
        ids.SetNumberOfComponents(1)
        ids.InsertNextValue(point_id)

        selection_node = vtkSelectionNode()
        #selection_node.SetContainingCellsOn()
        #selection_node.Initialize()
        selection_node.SetFieldType(vtkSelectionNode.POINT)
        selection_node.SetContentType(vtkSelectionNode.INDICES)
        selection_node.SetSelectionList(ids)
        actor = self._highlight_picker_by_selection_node(
            grid, selection_node, representation='points')
        return actor

    def _highlight_picker_cell(self,
                               cell_ids: list[int],
                               grid: vtkUnstructuredGrid) -> vtkLODActor:
        """won't handle multiple cell_ids/node_xyz"""
        selection_node = create_vtk_selection_node_by_cell_ids(cell_ids)
        actor = self._highlight_picker_by_selection_node(
            grid, selection_node, representation='surface')
        return actor

    def _highlight_picker_by_selection_node(self, grid: vtkUnstructuredGrid,
                                            selection_node: vtkSelectionNode,
                                            representation: str='surface',
                                            add_actor: bool=True) -> vtkLODActor:
        selection = vtkSelection()
        selection.AddNode(selection_node)

        extract_selection = vtkExtractSelection()
        extract_selection.SetInputData(0, grid)
        extract_selection.SetInputData(1, selection)
        extract_selection.Update()

        ugrid = extract_selection.GetOutput()
        actor = self.create_highlighted_actor(
            ugrid, representation=representation, add_actor=add_actor)
        return actor

        #-----------------------------------------------

    def create_highlighted_actor(self, ugrid: vtkUnstructuredGrid,
                                 representation: str='wire',
                                 add_actor: bool=True) -> list[vtkLODActor]:
        """creates a highlighted actor given a vtkUnstructuredGrid"""
        actor = create_highlighted_actor(
            self.gui, ugrid, representation=representation,
            add_actor=add_actor)
        return actor

    def _highlight_picker(self, unused_obj, unused_event) -> None:
        """
        pick a point/cell and highlight it
        """
        picker = self.cell_picker
        pixel_x, pixel_y = self.vtk_interactor.GetEventPosition()
        picker.Pick(pixel_x, pixel_y, 0, self.rend)

        cell_id = picker.GetCellId()
        #print('_probe_picker', cell_id)

        if cell_id < 0:
            pass
        else:
            #icase = self.gui.icase_fringe
            #if icase is None:
                #return

            world_position = picker.GetPickPosition()

            grid = self.gui.grid_selected
            if self.highlight_style == 'centroid':
                actor = self._highlight_picker_cell(cell_id, grid)
            elif self.highlight_style == 'node':
                actor = self._highlight_picker_node(cell_id, grid, world_position)
            else:
                raise RuntimeError('invalid highlight_style=%r' % self.highlight_style)
            self.actor = actor
            self.vtk_interactor.Render()

        if self.revert:
            self.cleanup_observer = self.setup_mouse_buttons(
                mode='default', left_button_down_cleanup=self._highlight_cleanup_callback)
            self.depress_buttons()

    def _highlight_cleanup_callback(self, obj, event) -> None:
        """this is the cleanup step to remove the highlighted actor"""
        if hasattr(self, 'actor'):
            self.rend.RemoveActor(self.actor)
        #self.vtk_interactor.RemoveObservers('LeftButtonPressEvent')
        self.vtk_interactor.RemoveObserver(self.cleanup_observer)
        #cleanup_observer = None

    def depress_buttons(self) -> None:
        """buttons may still be clicked after reverting, unpress them"""
        for button_name in self._depress_when_done:
            self.actions[button_name].setChecked(False)

    def _probe_picker(self, unused_obj, unused_event) -> None:
        """
        pick a point and apply the label based on the
        current displayed result
        """
        picker = self.cell_picker
        pixel_x, pixel_y = self.vtk_interactor.GetEventPosition()
        picker.Pick(pixel_x, pixel_y, 0, self.rend)

        cell_id = picker.GetCellId()
        #print('_probe_picker', cell_id)

        if cell_id < 0:
            pass
        else:
            gui = self.gui
            icase = gui.icase_fringe
            if icase is None:
                return

            world_position = picker.GetPickPosition()
            mark_actions = gui.mark_actions
            if 0:  # pragma : no cover
                camera = self.rend.GetActiveCamera()
                #focal_point = world_position
                out = mark_actions.get_result_by_xyz_cell_id(
                    world_position, cell_id)
                #if out is None:
                    #return
                _result_name, result_value, unused_node_id, node_xyz = out
                focal_point = node_xyz
                gui.log_info('focal_point = %s' % str(focal_point))
                self.setup_mouse_buttons(mode='default')

                # now we can actually modify the camera
                camera.SetFocalPoint(focal_point[0], focal_point[1], focal_point[2])
                camera.OrthogonalizeViewUp()
                probe_result_button = self.actions['probe_result']
                probe_result_button.setChecked(False)


                world_position = picker.GetPickPosition()
                cell_id = picker.GetCellId()
                #ds = picker.GetDataSet()
                #select_point = picker.GetSelectionPoint()
                gui.log_command("self.annotate_cell_picker()")
                gui.log_info("XYZ Global = %s" % str(world_position))
                #self.log_info("cell_id = %s" % cell_id)
                #self.log_info("data_set = %s" % ds)
                #self.log_info("selPt = %s" % str(select_point))

                #method = 'get_result_by_cell_id()' # self.model_type
                #print('pick_state =', self.pick_state)

            key = gui.case_keys[icase]
            location = gui.get_case_location(key)

            if location is None:
                return  # Normals
            if location == 'centroid':
                out = self._cell_centroid_pick(cell_id, world_position)
            elif location == 'node':
                out = self._cell_node_pick(cell_id, world_position)
            else:  # pragma: no cover
                raise RuntimeError(f'probe_picker: invalid pick location={location!r}')

            return_flag, duplicate_key, result_value, unused_result_name, xyz = out
            if return_flag is True:
                return

            # prevent duplicate labels with the same value on the same cell
            if duplicate_key is not None and duplicate_key in gui.label_ids[icase]:
                return
            gui.label_ids[icase].add(duplicate_key)

            #if 0:
                #result_value2, xyz2 = self.convert_units(case_key, result_value, xyz)
                #result_value = result_value2
                #xyz2 = xyz
            #x, y, z = world_position
            x, y, z = xyz
            text = '(%.3g, %.3g, %.3g); %s' % (x, y, z, result_value)
            text = str(result_value)
            assert icase in gui.label_actors, icase
            annotation = mark_actions.create_annotation(text, x, y, z)
            gui.label_actors[icase].append(annotation)
            self.vtk_interactor.Render()
        if self.revert:
            self.setup_mouse_buttons(mode='default')

    def _probe_picker_all(self, unused_obj, unused_event) -> None:
        """pick a point and apply the appropriate label for all cases"""
        picker = self.cell_picker
        pixel_x, pixel_y = self.vtk_interactor.GetEventPosition()
        picker.Pick(pixel_x, pixel_y, 0, self.rend)

        cell_id = picker.GetCellId()
        #print('_probe_picker', cell_id)

        #print(f'_probe_picker_all; cell_id={cell_id}')
        if cell_id < 0:
            pass
        else:
            gui = self.gui
            icase_temp = gui.icase_fringe
            if icase_temp is None:
                #print(f'return icase_temp={icase_temp}')
                return
            mark_actions = gui.mark_actions
            world_position = picker.GetPickPosition()

            ncases = len(gui.result_cases)
            for icase in range(ncases):
                gui.cycle_results(icase, update=False)
                key = gui.case_keys[icase]
                location = gui.get_case_location(key)

                # RuntimeError: invalid pick location
                if location is None:
                    gui.log_warning(f'skipping probing on case {icase}/{ncases} key={key} location={location}')
                    continue

                try:
                    if location == 'centroid':
                        out = self._cell_centroid_pick(cell_id, world_position, icase=icase)
                    elif location == 'node':
                        out = self._cell_node_pick(cell_id, world_position, icase=icase)
                    else:
                        raise RuntimeError(f'probe_picker_all: invalid pick location={location!r}')
                except IndexError:
                    # IndexError: out of range?
                    gui.log_warning(f'failed probing on case {icase}/{ncases} key={key} location={location}')
                    continue

                return_flag, duplicate_key, result_value, result_name, xyz = out
                if return_flag is True:
                    continue

                # prevent duplicate labels with the same value on the same cell
                if duplicate_key is not None and duplicate_key in gui.label_ids[icase]:
                    continue
                gui.label_ids[icase].add(duplicate_key)

                #if 0:
                    #result_value2, xyz2 = self.convert_units(case_key, result_value, xyz)
                    #result_value = result_value2
                    #xyz2 = xyz
                #x, y, z = world_position
                x, y, z = xyz
                text = '(%.3g, %.3g, %.3g); %s' % (x, y, z, result_value)
                text = str(result_value)
                assert icase in gui.label_actors, icase
                #print(f'icase={icase}/{ncases}: text={text!r}')
                annotation = mark_actions.create_annotation(text, x, y, z)
                gui.label_actors[icase].append(annotation)
            gui.cycle_results(icase_temp, show_msg=False)
            self.vtk_interactor.Render()
        if self.revert:
            self.setup_mouse_buttons(mode='default')

    def _cell_node_pick(self, cell_id: int,
                        world_position: np.ndarray,
                        icase: Optional[int]=None):
        gui = self.gui
        if icase is None:
            icase = gui.icase
        assert icase in gui.label_actors, f'icase={icase!r} result_name={result_name!r}'

        duplicate_key = None
        pick_state = self.pick_state
        if pick_state == 'node/centroid':
            return_flag = False
            out = gui.mark_actions.get_result_by_xyz_cell_id(
                world_position, cell_id, icase=icase)
            if out is None:
                print(f'MouseActions._cell_node_pick bug for icase={icase}')
                return_flag = True
                return return_flag, None, None, None, None
            result_name, result_value, node_id, xyz = out
            assert not isinstance(xyz, int), xyz
            duplicate_key = node_id
        else:
            method = f'get_nodal_{gui.format}_result_pick_state_{pick_state}_by_xyz_cell_id'
            if hasattr(self, method):
                methodi = getattr(self, method)
                return_flag, unused_value = methodi(world_position, cell_id)
                if return_flag is True:
                    return return_flag, None, None, None, None
            else:
                msg = "pick_state is set to 'centroidal', but the result is 'nodal'\n"
                msg += f'  cannot find: self.{method}(xyz, cell_id)'
                gui.log_error(msg)
            return return_flag, None, None, None
        msg = f'{result_name} = {result_value}'
        #if gui.result_name in {'Node_ID', 'Node ID', 'NodeID'}:
        x1, y1, z1 = xyz
        x2, y2, z2 = world_position
        msg += f'; xyz=({x1}, {y1}, {z1}); pierce_xyz=({x2}, {y2}, {z2})'
        gui.log_info(msg)
        return return_flag, duplicate_key, result_value, result_name, xyz

    def _cell_centroid_pick(self, cell_id: int,
                            world_position: np.ndarray,
                            icase: Optional[int]=None):
        gui = self.gui
        if icase is None:
            icase = gui.icase
        assert icase in gui.label_actors, f'icase={icase!r} result_name={result_name!r}'

        duplicate_key = None
        pick_state = self.pick_state
        if pick_state == 'node/centroid':
            return_flag = False
            duplicate_key = cell_id
            out = gui.get_result_by_cell_id(cell_id, world_position, icase=icase)
            result_name, result_value, xyz = out
        else:
            #cell = self.grid.GetCell(cell_id)
            # get_nastran_centroidal_pick_state_nodal_by_xyz_cell_id()
            method = f'get_centroidal_{gui.format}_result_pick_state_{pick_state}_by_xyz_cell_id'
            if hasattr(self, method):
                methodi = getattr(self, method)
                return_flag, unused_value = methodi(world_position, cell_id)
                if return_flag is True:
                    return return_flag, None, None, None, None
            else:
                msg = "pick_state is set to 'nodal', but the result is 'centroidal'\n"
                msg += f'  cannot find: self.{method}(xyz, cell_id)'
                gui.log_error(msg)
            return return_flag, None, None, None
        gui.log_info(f'{result_name} = {result_value}')
        return return_flag, duplicate_key, result_value, result_name, xyz

    #---------------------------------------------------------------------------
    @property
    def actions(self):
        return self.gui.actions

    @property
    def rend(self) -> vtkRenderer:
        return self.gui.rend

    @property
    def vtk_interactor(self):
        return self.gui.vtk_interactor

    @property
    def iren(self):
        return self.vtk_interactor

    @property
    def area_picker(self):
        return self.gui.area_picker

    @property
    def cell_picker(self):
        return self.gui.cell_picker

    @property
    def node_picker(self):
        return self.gui.node_picker

    def get_grid(self, name: str):
        return self.grid

    def get_grid_selected(self, name: str):
        try:
            return self.grid_selected
        except Exception:
            return self.grid

    @property
    def grid(self):
        return self.gui.grid

    @property
    def node_ids(self) -> np.ndarray:
        return self.gui.node_ids

    @property
    def element_ids(self) -> np.ndarray:
        return self.gui.element_ids

    def set_focal_point(self, focal_point):
        return self.gui.set_focal_point(focal_point)

    def window(self):
        return self.gui.window()

    def zoom(self, zoom_factor: float):
        return self.gui.zoom(zoom_factor)

    @property
    def view_actions(self):
        return self.gui.view_actions
    def on_pan_left(self, event):
        """helper method for trackball camera"""
        self.view_actions.on_pan_left(event)

    def on_pan_right(self, event) -> None:
        """helper method for trackball camera"""
        self.view_actions.on_pan_right(event)

    def on_pan_up(self, event) -> None:
        """helper method for trackball camera"""
        self.view_actions.on_pan_up(event)

    def on_pan_down(self, event) -> None:
        """helper method for trackball camera"""
        self.view_actions.on_pan_down(event)

    def get_node_ids(self, model_name, ids=None):
        """wrapper around node_ids"""
        return self.gui.get_node_ids(model_name, ids)

    def get_element_ids(self, model_name, ids=None):
        """wrapper around node_ids"""
        return self.gui.get_element_ids(model_name, ids)

    @property
    def settings(self) -> Settings:
        return self.gui.settings
