from __future__ import print_function
import vtk

class TrackballStyleCamera(vtk.vtkInteractorStyleTrackballCamera):
    #https://stackoverflow.com/questions/33108670/arrow-key-events-in-vtk-on-windows
    def __init__(self, iren, parent):
        self.parent = parent
        vtk.vtkInteractorStyleTrackballCamera.__init__(self, iren)
        #self.AddObserver("CharEvent", self.onKeyPressEvent)

        self.AddObserver("KeyPressEvent", self.keyPressEvent)

    def keyPressEvent(self, unused_obj, event):
        key = self.parent.iren.GetKeySym()
        if key == 'Left':
            self.parent.on_pan_left(event)
        elif key == 'Right':
            self.parent.on_pan_right(event)
        elif key == 'Up':
            self.parent.on_pan_up(event)
        elif key == 'Down':
            self.parent.on_pan_down(event)

