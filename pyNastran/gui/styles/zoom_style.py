"""defines the ZoomStyle class"""
import vtk

#left_button_down=self._zoom_picker,
#left_button_up=self._zoom_picker,
#right_button_down=self._zoom_reset,

class ZoomStyle(vtk.vtkInteractorStyleRubberBandZoom):
    """Custom Rubber Band Zoom"""
    def __init__(self, parent=None):
        """creates the ZoomStyle instance"""
        self.AddObserver("LeftButtonPressEvent", self.left_button_press_event)
        self.AddObserver("LeftButtonReleaseEvent", self.left_button_release_event)
        self.AddObserver("RightButtonPressEvent", self.right_button_press_event)
        self.parent = parent
        self.zoom_button = self.parent.actions['zoom']
        self.picker_points = []

    #def leftButtonPressEvent(self, obj, event):
        #pass

    def left_button_press_event(self, obj, event):
        """
        gets the first point
        """
        self.OnLeftButtonDown()
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        self.picker_points.append((pixel_x, pixel_y))

    def left_button_release_event(self, obj, event):
        """
        gets the second point and zooms

        TODO: doesn't handle panning of the camera to center the image
              with respect to the selected limits
        """
        self.OnLeftButtonUp()
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        self.picker_points.append((pixel_x, pixel_y))

        camera = self.parent.rend.GetActiveCamera()
        x, y, z = camera.GetPosition()
        p1x, p1y = self.picker_points[0]
        p2x, p2y = self.picker_points[1]

        dx = abs(p1x - p2x)
        dy = abs(p1y - p2y)
        #x_avg = (p1x + p2x) / 2.
        #y_avg = (p1y + p2y) / 2.

        main_window = self.parent.window()
        width = main_window.frameGeometry().width()
        height = main_window.frameGeometry().height()
        #print('dx=%s dy=%s' % (dx, dy))

        # otherwise it's a failed zoom (they didn't hold the button down)
        self.picker_points = []
        if dx > 0 and dy > 0:
            #xmin = min(p1x, p2x)
            #ymin = min(p1y, p2y)
            #xmax = max(p1x, p2x)
            #ymax = max(p1y, p2y)

            aspect_ratio_x = width / dx
            aspect_ratio_y = height / dy
            zoom_factor = min([aspect_ratio_x, aspect_ratio_y])

            #distance = camera.GetDistance()
            #a = vtk.vtkCamera()


              # +---------+ --- ymax
              # |         |
              # |         |
              # |         |
              # +---------+ --- ymin
              #
            #camera.SetScreenBottomLeft(xmin, ymin)
            #camera.SetScreenBottomRight(float, float)
            #camera.SetScreenTopRight(float, float)

            #print('  p1 =', p1x, p1y)
            #print('  p2 =', p2x, p2y)
            #print('  z=%s distance=%s' % (z, distance))
            #print('  zoom_factor = %s\n' % aspect_ratio_x, aspect_ratio_y)
            #camera.SetPosition(x, y, z)
            self.parent.zoom(zoom_factor)

            self.zoom_button.setChecked(False)
            self.parent.setup_mouse_buttons(mode='default')
            #self.parent.actions['zoom'].SetChecked(False)


    def right_button_press_event(self, obj, event):
        """cancels the zoom button"""
        self.zoom_button.setChecked(False)
        self.parent.setup_mouse_buttons(mode='default')
        self.parent.vtk_interactor.Render()
