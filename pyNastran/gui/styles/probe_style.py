"""
defines the ProbeResultStyle class
"""
import vtk

class ProbeResultStyle(vtk.vtkInteractorStyleTrackballCamera):
    """Custom TrackballCamera"""

    def __init__(self, parent=None):
        """creates the ProbeResultStyle instance"""
        self.AddObserver("LeftButtonPressEvent", self.left_button_press_event)
        self.parent = parent
        self.probe_result_button = self.parent.actions['probe_result']

    def left_button_press_event(self, obj, event):
        """pick a point and apply the label based on the current displayed result"""
        picker = self.parent.cell_picker
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        picker.Pick(pixel_x, pixel_y, 0, self.parent.rend)

        cell_id = picker.GetCellId()
        #print('_rotation_center_cell_picker', cell_id)

        if cell_id < 0:
            pass
        else:
            world_position = picker.GetPickPosition()
            if 0:
                camera = self.parent.rend.GetActiveCamera()
                #focal_point = world_position
                out = self.parent.get_result_by_xyz_cell_id(world_position, cell_id)
                result_name, result_value, node_id, node_xyz = out
                focal_point = node_xyz
                self.parent.log_info('focal_point = %s' % str(focal_point))
                self.parent.mouse_actions.setup_mouse_buttons(mode='default')

                # now we can actually modify the camera
                camera.SetFocalPoint(focal_point[0], focal_point[1], focal_point[2])
                camera.OrthogonalizeViewUp()
                self.probe_result_button.setChecked(False)


                world_position = picker.GetPickPosition()
                cell_id = picker.GetCellId()
                #ds = picker.GetDataSet()
                #select_point = picker.GetSelectionPoint()
                self.parent.log_command("annotate_cell_picker()")
                self.parent.log_info("XYZ Global = %s" % str(world_position))
                #self.parent.log_info("cell_id = %s" % cell_id)
                #self.parent.log_info("data_set = %s" % ds)
                #self.parent.log_info("selPt = %s" % str(select_point))

                #method = 'get_result_by_cell_id()' # self.parent.model_type
                #print('pick_state =', self.parent.pick_state)

            icase = self.parent.icase
            key = self.parent.case_keys[icase]
            location = self.parent.get_case_location(key)

            if location == 'centroid':
                out = self.parent._cell_centroid_pick(cell_id, world_position)
            elif location == 'node':
                out = self.parent._cell_node_pick(cell_id, world_position)
            else:
                raise RuntimeError('invalid pick location=%r' % location)

            return_flag, duplicate_key, result_value, result_name, xyz = out
            if return_flag is True:
                return

            # prevent duplicate labels with the same value on the same cell
            if duplicate_key is not None and duplicate_key in self.parent.label_ids[result_name]:
                return
            self.parent.label_ids[result_name].add(duplicate_key)

            #if 0:
                #result_value2, xyz2 = self.parent.convert_units(result_name, result_value, xyz)
                #result_value = result_value2
                #xyz2 = xyz
            #x, y, z = world_position
            x, y, z = xyz
            text = '(%.3g, %.3g, %.3g); %s' % (x, y, z, result_value)
            text = str(result_value)
            assert result_name in self.parent.label_actors, result_name
            self.label_actors[result_name].append(self.parent.create_annotation(text, x, y, z))
            self.parent.vtk_interactor.Render()
            self.parent.vtk_interactor.Update()
            self.parent.Update()
            self.parent.log_command('update...')

    #def right_button_press_event(self, obj, event):
        #"""cancels the probe button"""
        #self.probe_result_button.setChecked(False)
        #self.parent.mouse_actions.setup_mouse_buttons(mode='default')
        #self.parent.vtk_interactor.Render()
