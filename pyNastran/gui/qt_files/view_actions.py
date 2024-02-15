# coding: utf-8
from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
#import vtkmodules
from pyNastran.gui.vtk_common_core import vtkMath
from pyNastran.gui.vtk_rendering_core import vtkCamera, vtkRenderer
from pyNastran.gui.gui_objects.settings import Settings

from pyNastran.gui.utils.vtk.gui_utils import numpy_array_to_vtk_array, flip_actor_visibility

#from pyNastran.gui.gui_objects.coord_properties import CoordProperties
if TYPE_CHECKING:  # pragma: no cover
    import numpy as np


class ViewActions:
    def __init__(self, gui):
        self.gui = gui
        self.is_wireframe = False

    def log_command(self, msg: str) -> None:
        return self.gui.log_command

    def on_reset_camera(self) -> None:
        self.log_command('self.on_reset_camera()')
        self._simulate_key_press('r')
        self.vtk_interactor.Render()

    def on_pan_left(self, event):
        """https://semisortedblog.wordpress.com/2014/09/04/building-vtk-user-interfaces-part-3c-vtk-interaction"""
        camera, cam, focal = self._setup_pan()

        # Create a vector that points upward, i.e. (0, 1, 0)
        up = camera.GetViewUp() # We don't want roll
        vec = [0, 0, 0]
        new_cam = [0, 0, 0]
        new_focal = [0, 0, 0]

        # Calculate the forward pointing unit-vector 'vec' again in the same way,
        # i.e. the normalized vector of focal point – camera position
        vtkMath.Subtract(focal, cam, vec)

        vec[1] = 0 #We don't want roll
        vtkMath.Normalize(vec)

        # Calculate the cross product of the forward vector by the up vector,
        # which will give us an orthogonal vector pointing right relative to
        #the camera
        vtkMath.Cross(vec, up, vec)

        # Add this to the camera position and focal point to move it right
        # new_cam = cam + vec
        vtkMath.Add(cam, vec, new_cam)

        # new_focal = focal + vec
        vtkMath.Add(focal, vec, new_focal)
        self._set_camera_position_focal_point(camera, new_cam, new_focal)

    def on_pan_right(self, event):
        """https://semisortedblog.wordpress.com/2014/09/04/building-vtk-user-interfaces-part-3c-vtk-interaction"""
        camera, cam, focal = self._setup_pan()

        # Create a vector that points upward, i.e. (0, 1, 0)
        up = camera.GetViewUp() # We don't want roll
        vec = [0, 0, 0]
        new_cam = [0, 0, 0]
        new_focal = [0, 0, 0]

        # Calculate the forward pointing unit-vector 'vec' again in the same way,
        # i.e. the normalized vector of focal point – camera position
        vtkMath.Subtract(focal, cam, vec)

        vec[1] = 0 #We don't want roll
        vtkMath.Normalize(vec)

        # Calculate the cross product of the forward vector by the up vector,
        # which will give us an orthogonal vector pointing right relative to
        #the camera
        #vec = up x vec
        vtkMath.Cross(vec, up, vec)

        # Subtract vec from the camera position and focal point to move it right
        # new_cam = cam - vec
        vtkMath.Subtract(cam, vec, new_cam)

        # new_focal = focal - vec
        vtkMath.Subtract(focal, vec, new_focal)
        self._set_camera_position_focal_point(camera, new_cam, new_focal)

    def on_pan_up(self, event):
        """not 100% on this"""
        camera, cam, focal = self._setup_pan()

        # Create a 'vec' vector that will be the direction of movement
        # (numpad 8 and 5 generate movement along the z-axis; numpad 4
        # and 6 along the x-axis; numpad 7 and 9 along the y-axis)
        vec = camera.GetViewUp() # We don't want roll
        new_cam = [0, 0, 0]
        new_focal = [0, 0, 0]

        # Add the movement to the current camera position and focal point,
        # and save these in 'new_cam' and 'new_focal' respectively
        vtkMath.Subtract(cam, vec, new_cam)

        # new_focal = focal - vec
        vtkMath.Subtract(focal, vec, new_focal)
        self._set_camera_position_focal_point(camera, new_cam, new_focal)

    def on_pan_down(self, event):
        """not 100% on this"""
        camera, cam, focal = self._setup_pan()

        # Create a 'vec' vector that will be the direction of movement
        # (numpad 8 and 5 generate movement along the z-axis; numpad 4
        # and 6 along the x-axis; numpad 7 and 9 along the y-axis)
        vec = camera.GetViewUp() # We don't want roll
        new_cam = [0, 0, 0]
        new_focal = [0, 0, 0]

        # Add the movement to the current camera position and focal point,
        # and save these in 'new_cam' and 'new_focal' respectively
        vtkMath.Add(cam, vec, new_cam)

        # new_focal = focal + vec
        vtkMath.Add(focal, vec, new_focal)
        self._set_camera_position_focal_point(camera, new_cam, new_focal)

    def _setup_pan(self):
        camera = self.rend.GetActiveCamera()
        cam = camera.GetPosition()
        focal = camera.GetFocalPoint()
        return camera, cam, focal

    def _set_camera_position_focal_point(self, camera, new_cam, new_focal):
        """Set the camera position and focal point to the new vectors"""
        camera.SetPosition(new_cam)
        camera.SetFocalPoint(new_focal)

        # Update the clipping range of the camera
        self.rend.ResetCameraClippingRange()
        self.Render()

    #---------------------------------------------------------------------------

    def on_show_hide_min_actor(self) -> None:
        """flips the status of the min label actor"""
        self.settings.is_min_visible = not self.settings.is_min_visible
        actor = self.gui.min_max_actors[0]
        flip_actor_visibility(actor)
        self.Render()

    def on_show_hide_max_actor(self) -> None:
        """flips the status of the max label actor"""
        self.settings.is_max_visible = not self.settings.is_max_visible
        actor = self.gui.min_max_actors[1]
        flip_actor_visibility(actor)
        self.Render()

    def on_increase_magnification(self) -> None:
        """zoom in"""
        self.zoom(1.1)

    def on_decrease_magnification(self) -> None:
        """zoom out"""
        self.zoom(1.0 / 1.1)

    def on_rotate_clockwise(self) -> None:
        """rotate clockwise"""
        self.rotate(15.)

    def on_rotate_cclockwise(self) -> None:
        """rotate counter clockwise"""
        self.rotate(-15.)

    def azimuth(self, azimuth_deg: float, render: bool=True) -> None:
        """see the gui"""
        camera = self.GetCamera()
        camera.Azimuth(-azimuth_deg)
        camera.Modified()
        if render:
            self.vtk_interactor.Render()
        self.log_command('self.azimuth(%s)' % azimuth_deg)

    def pitch(self, pitch_deg: float, render: bool=True) -> None:
        """see the gui"""
        camera = self.GetCamera()
        camera.Pitch(-pitch_deg)
        camera.Modified()
        if render:
            self.vtk_interactor.Render()
        self.log_command('self.pitch(%s)' % pitch_deg)

    def yaw(self, yaw_deg: float, render: bool=True) -> None:
        """see the gui"""
        camera = self.GetCamera()
        camera.Yaw(-yaw_deg)
        camera.Modified()
        if render:
            self.vtk_interactor.Render()
        self.log_command('self.yaw(%s)' % yaw_deg)

    def rotate(self, rotate_deg: float, render: bool=True) -> None:
        """see the gui"""
        camera = self.GetCamera()
        camera.Roll(-rotate_deg)
        camera.Modified()
        if render:
            self.vtk_interactor.Render()
        self.log_command('self.rotate(%s)' % rotate_deg)

    def zoom(self, value: float, render: bool=True) -> None:
        camera = self.GetCamera()
        camera.Zoom(value)
        camera.Modified()
        if render:
            self.vtk_interactor.Render()
        self.log_command('self.zoom(%s)' % value)

    def set_focal_point(self, focal_point: np.ndarray, render: bool=True) -> None:
        """
        Parameters
        ----------
        focal_point : (3, ) float ndarray
            The focal point
            [ 188.25109863 -7. -32.07858658]
        """
        camera = self.rend.GetActiveCamera()
        self.log_command("self.set_focal_point(focal_point=%s)" % str(focal_point))

        # now we can actually modify the camera
        camera.SetFocalPoint(focal_point[0], focal_point[1], focal_point[2])
        camera.OrthogonalizeViewUp()
        if render:
            self.vtk_interactor.Render()

    def on_surface(self, render: bool=True) -> None:
        """sets the main/toggle actors to surface"""
        if self.is_wireframe:
            self.log_command('self.on_surface()')
            for name, actor in self.gui.geometry_actors.items():
                #if name != 'main':
                    #print('name: %s\nrep: %s' % (
                        #name, self.geometry_properties[name].representation))
                representation = self.gui.geometry_properties[name].representation
                if name == 'main' or representation in ['main', 'toggle']:
                    prop = actor.GetProperty()

                    prop.SetRepresentationToSurface()
            self.is_wireframe = False
            if render:
                self.vtk_interactor.Render()

    def on_wireframe(self, render: bool=True) -> None:
        """sets the main/toggle actors to wirefreme"""
        if not self.is_wireframe:
            self.log_command('self.on_wireframe()')
            for name, actor in self.gui.geometry_actors.items():
                #if name != 'main':
                    #print('name: %s\nrep: %s' % (
                        #name, self.geometry_properties[name].representation))
                representation = self.gui.geometry_properties[name].representation
                if name == 'main' or representation in ['main', 'toggle']:
                    prop = actor.GetProperty()
                    prop.SetRepresentationToWireframe()
                #prop.SetRepresentationToPoints()
                #prop.RenderPointsAsSpheresOn()
                #prop.SetLighting(False)
                #prop.SetInterpolationToFlat()
                #prop.GetPointSize()
                #prop.SetPointSize(5.0)
                #prop.ShadingOff()
            if render:
                self.vtk_interactor.Render()
            self.is_wireframe = True

    #---------------------------------------------------------------------------
    # camera
    def update_camera(self, code: str) -> None:
        camera = self.GetCamera()
        #print("code =", code)
        if code == '+x':  # set x-axis
            # +z up
            # +y right
            # looking forward
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(1., 0., 0.)
        elif code == '-x':  # set x-axis
            # +z up
            # +y to the left (right wing)
            # looking aft
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(-1., 0., 0.)

        elif code == '+y':  # set y-axis
            # +z up
            # +x aft to left
            # view from right wing
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(0., 1., 0.)
        elif code == '-y':  # set y-axis
            # +z up
            # +x aft to right
            # view from left wing
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(0., -1., 0.)

        elif code == '+z':  # set z-axis
            # +x aft
            # +y up (right wing up)
            # top view
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 1., 0.)
            camera.SetPosition(0., 0., 1.)
        elif code == '-z':  # set z-axis
            # +x aft
            # -y down (left wing up)
            # bottom view
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., -1., 0.)
            camera.SetPosition(0., 0., -1.)
        else:
            self.log_error('invalid camera code...%r' % code)
            return
        self._update_camera(camera)
        self.rend.ResetCamera()
        self.log_command('self.update_camera(%r)' % code)

    def _update_camera(self, camera: Optional[vtkCamera]=None) -> None:
        if camera is None:
            camera = self.GetCamera()
        camera.Modified()
        self.vtk_interactor.Render()

    #---------------------------------------------------------------------------
    def Render(self) -> None:
        self.vtk_interactor.GetRenderWindow().Render()

    def GetCamera(self) -> vtkCamera:
        return self.rend.GetActiveCamera()

    #@property
    #def settings(self):
        #return self.gui.settings

    @property
    def rend(self) -> vtkRenderer:
        return self.gui.rend

    @property
    def vtk_interactor(self) -> Any:
        return self.gui.vtk_interactor

    @property
    def settings(self) -> Settings:
        return self.gui.settings
