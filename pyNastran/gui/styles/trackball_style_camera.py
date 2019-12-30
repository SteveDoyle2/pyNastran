from __future__ import annotations
from typing import TYPE_CHECKING
import vtk
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.qt_files.mouse_actions import MouseActions

class TrackballStyleCamera(vtk.vtkInteractorStyleTrackballCamera):
    #https://stackoverflow.com/questions/33108670/arrow-key-events-in-vtk-on-windows
    def __init__(self, iren, parent: MouseActions):
        self.parent = parent

        # works
        vtk.vtkInteractorStyleTrackballCamera.__init__(self)

        # TypeError: object.__init__() takes no parameters
        #super(TrackballStyleCamera, self).__init__(self, iren)

        # TypeError: object.__init__() takes no parameters
        #super(TrackballStyleCamera, self).__init__(self)

        #self.AddObserver("CharEvent", self.onKeyPressEvent)

        self.AddObserver("KeyPressEvent", self.keyPressEvent)

    def keyPressEvent(self, unused_obj, event):
        key = self.parent.iren.GetKeySym()
        is_control = self.parent.iren.GetControlKey()
        if is_control:
            view = self.parent.gui.view_actions
            if key == 'Left':
                view.yaw(-5)
            elif key == 'Right':
                view.yaw(5)
            elif key == 'Up':
                # azimuth failed...
                view.pitch(5.)
            elif key == 'Down':
                view.pitch(-5.)
        else:
            if key == 'Left':
                self.parent.on_pan_left(event)
            elif key == 'Right':
                self.parent.on_pan_right(event)
            elif key == 'Up':
                self.parent.on_pan_up(event)
            elif key == 'Down':
                self.parent.on_pan_down(event)
