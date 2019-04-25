"""
defines the RotationCenterStyle class
"""
import vtk

class RotationCenterStyle(vtk.vtkInteractorStyleTrackballCamera):
    """Custom TrackballCamera"""

    def __init__(self, parent=None):
        """creates the RotationCenterStyle instance"""
        self.AddObserver("LeftButtonPressEvent", self.left_button_press_event)
        self.parent = parent
        self.rotation_center_button = self.parent.actions['rotation_center']

    def left_button_press_event(self, obj, event):
        """pick a point and apply the label based on the current displayed result"""
        picker = self.parent.cell_picker
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        picker.Pick(pixel_x, pixel_y, 0, self.parent.rend)

        cell_id = picker.GetCellId()
        #print('_rotation_center_cell_picker', cell_id)

        if cell_id < 0:
            return
        camera = self.parent.rend.GetActiveCamera()
        world_position = picker.GetPickPosition()
        focal_point = self.parent._get_closest_node_xyz(cell_id, world_position)

        self.parent.set_focal_point(focal_point)
        self.parent.setup_mouse_buttons(mode='default')
        self.rotation_center_button.setChecked(False)

    def right_button_press_event(self, obj, event):
        """cancels the probe button"""
        self.rotation_center_button.setChecked(False)
        self.parent.setup_mouse_buttons(mode='default')
        self.parent.vtk_interactor.Render()
