from __future__ import print_function
from copy import deepcopy
from pyNastran.gui.gui_interface.camera.camera import CameraWindow


def set_camera_menu(self):
    camera = self.rend.GetActiveCamera()
    #position = camera.GetPosition()
    #clip_range = camera.GetClippingRange()
    #focal_point = camera.GetFocalPoint()

    data = {'cameras' : self.cameras}
    window = CameraWindow(data, win_parent=self)
    window.show()
    window.exec_()

    if data['clicked_ok']:
        self.cameras = deepcopy(data['cameras'])
        #self._apply_camera(data)
    #self.log_info('position = %s' % str(position))
    #self.log_info('clip_range = %s' % str(clip_range))
    #self.log_info('focal_point = %s' % str(focal_point))
