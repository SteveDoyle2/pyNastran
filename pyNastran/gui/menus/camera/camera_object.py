"""
defines:
 - CameraObject

"""
from copy import deepcopy
from pyNastran.gui.menus.camera.camera import CameraWindow

class CameraObject:
    """defines CameraObject"""
    def __init__(self, gui):
        """creates CameraObject"""
        self.gui = gui
        self.cameras = {}
        self._camera_window_shown = False
        self._camera_window = None

    def set_font_size(self, font_size):
        """sets the font size for the camera window"""
        if self._camera_window_shown:
            self._camera_window.set_font_size(font_size)

    def set_camera_menu(self):
        """loads the camera window"""
        #camera = self.gui.rend.GetActiveCamera()
        #position = camera.GetPosition()
        #clip_range = camera.GetClippingRange()
        #focal_point = camera.GetFocalPoint()

        data = {'cameras' : self.cameras}
        window = CameraWindow(data, win_parent=self.gui)
        window.show()
        window.exec_()

        if data['clicked_ok']:
            self.cameras = deepcopy(data['cameras'])
            #self._apply_camera(data)
        #self.log_info('position = %s' % str(position))
        #self.log_info('clip_range = %s' % str(clip_range))
        #self.log_info('focal_point = %s' % str(focal_point))

    #def _apply_camera(self, data):
        #name = data['name']
        #self.cameras = deepcopy(data['cameras'])
        #self.on_set_camera(name)

    def get_camera_data(self):
        """see ``set_camera_data`` for arguments"""
        camera = self.gui.rend.GetActiveCamera()
        position = camera.GetPosition()
        focal_point = camera.GetFocalPoint()
        view_angle = camera.GetViewAngle()
        view_up = camera.GetViewUp()

        # TODO: do I need clip_range and parallel_scale?
        clip_range = camera.GetClippingRange()
        parallel_scale = camera.GetParallelScale()
        parallel_proj = None
        if hasattr(camera, 'GetParralelProjection'):
            parallel_proj = camera.GetParralelProjection()
        distance = camera.GetDistance()

        # clip_range, view_up, distance
        camera_data = {
            'position' : position,
            'focal_point' : focal_point,
            'view_angle' : view_angle,
            'view_up' : view_up,
            'clip_range' : clip_range,
            'parallel_scale' : parallel_scale,
            'prallel_proj' : parallel_proj,
            'distance' : distance,
        }
        return camera_data

    def on_set_camera_data(self, camera_data, show_log=True):
        """
        Sets the current camera

        Parameters
        ----------
        camera_data : Dict[key] : value
            defines the camera
            position : (float, float, float)
                where am I is xyz space
            focal_point : (float, float, float)
                where am I looking
            view_angle : float
                field of view (angle); perspective only?
            view_up : (float, float, float)
                up on the screen vector
            clip_range : (float, float)
                start/end distance from camera where clipping starts
            parallel_scale : float
                ???
            parallel_projection : bool (0/1)
                flag?
                TODO: not used
            distance : float
                distance to the camera

        i_vector = focal_point - position
        j'_vector = view_up

        use:
           i x j' -> k
           k x i -> j
           or it's like k'
        """
        position = camera_data['position']
        focal_point = camera_data['focal_point']
        view_angle = camera_data['view_angle']
        view_up = camera_data['view_up']
        clip_range = camera_data['clip_range']
        parallel_scale = camera_data['parallel_scale']
        unused_parallel_proj = camera_data['prallel_proj']
        distance = camera_data['distance']

        camera = self.gui.rend.GetActiveCamera()
        camera.SetPosition(position)
        camera.SetFocalPoint(focal_point)
        camera.SetViewAngle(view_angle)
        camera.SetViewUp(view_up)
        camera.SetClippingRange(clip_range)

        camera.SetParallelScale(parallel_scale)
        #parallel_proj

        camera.SetDistance(distance)

        camera.Modified()
        self.gui.vtk_interactor.Render()
        if show_log:
            self.gui.log_command('on_set_camera_data(%s)' % str(camera_data))

    def on_set_camera(self, name, show_log=True):
        """see ``set_camera_data`` for arguments"""
        camera_data = self.cameras[name]
        self.on_set_camera_data(camera_data, show_log=show_log)
