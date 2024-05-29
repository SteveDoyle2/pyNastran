"""
defines the ProbeResultStyle class
"""
from __future__ import annotations
from typing import TYPE_CHECKING
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui import MainWindow


class ProbeResultStyle(vtkInteractorStyleTrackballCamera):
    """Custom TrackballCamera"""

    def __init__(self, parent=None):
        """creates the ProbeResultStyle instance"""
        self.AddObserver("LeftButtonPressEvent", self.left_button_press_event)
        self.parent = parent
        self.probe_result_button = self.parent.actions['probe_result']

    def left_button_press_event(self, obj, event) -> None:
        """pick a point and apply the label based on the current displayed result"""
        gui: MainWindow = self.parent
        picker = gui.cell_picker
        pixel_x, pixel_y = gui.vtk_interactor.GetEventPosition()
        picker.Pick(pixel_x, pixel_y, 0, gui.rend)

        cell_id = picker.GetCellId()
        if cell_id < 0:
            return

        #print('_rotation_center_cell_picker', cell_id)
        world_position = picker.GetPickPosition()
        if 0:
            camera = gui.rend.GetActiveCamera()
            #focal_point = world_position
            out = gui.get_result_by_xyz_cell_id(world_position, cell_id)
            result_name, result_value, node_id, focal_point = out
            gui.log_info('focal_point = %s' % str(focal_point))
            gui.mouse_actions.setup_mouse_buttons(mode='default')

            # now we can actually modify the camera
            camera.SetFocalPoint(focal_point[0], focal_point[1], focal_point[2])
            camera.OrthogonalizeViewUp()
            self.probe_result_button.setChecked(False)


            world_position = picker.GetPickPosition()
            cell_id = picker.GetCellId()
            #ds = picker.GetDataSet()
            #select_point = picker.GetSelectionPoint()
            gui.log_command("self.annotate_cell_picker()")
            gui.log_info("XYZ Global = %s" % str(world_position))
            #gui.log_info("cell_id = %s" % cell_id)
            #gui.log_info("data_set = %s" % ds)
            #gui.log_info("selPt = %s" % str(select_point))

            #method = 'get_result_by_cell_id()' # gui.model_type
            #print('pick_state =', gui.pick_state)

        icase = gui.icase
        key = gui.case_keys[icase]
        location = gui.get_case_location(key)

        if location == 'centroid':
            out = gui._cell_centroid_pick(cell_id, world_position)
        elif location == 'node':
            out = gui._cell_node_pick(cell_id, world_position)
        else:  # pragma: no cover
            raise RuntimeError(f'invalid pick location={location!r}')

        return_flag, duplicate_key, result_value, result_name, xyz = out
        if return_flag is True:
            return
        if xyz is None:
            gui.log.warning(f'ProbeResult has xyz=None from gui._cell_{location}_pick({cell_id}, {world_position})')
            return


        # prevent duplicate labels with the same value on the same cell
        if duplicate_key is not None and duplicate_key in gui.label_ids[result_name]:
            return
        gui.label_ids[result_name].add(duplicate_key)

        #if 0:
            #result_value2, xyz2 = gui.convert_units(result_name, result_value, xyz)
            #result_value = result_value2
            #xyz2 = xyz
        #x, y, z = world_position
        x, y, z = xyz
        text = '(%.3g, %.3g, %.3g); %s' % (x, y, z, result_value)
        text = str(result_value)
        assert result_name in gui.label_actors, result_name
        self.label_actors[result_name].append(gui.create_annotation(text, x, y, z))
        gui.vtk_interactor.Render()
        gui.vtk_interactor.Update()
        gui.Update()
        gui.log_command('update...')

    #def right_button_press_event(self, obj, event):
        #"""cancels the probe button"""
        #self.probe_result_button.setChecked(False)
        #self.parent.mouse_actions.setup_mouse_buttons(mode='default')
        #self.parent.vtk_interactor.Render()
