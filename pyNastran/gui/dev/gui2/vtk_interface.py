from __future__ import annotations
from typing import Optional, TYPE_CHECKING

from vtkmodules.vtkRenderingCore import vtkRenderer
from pyNastran.gui.qt_files.colors import BLACK_FLOAT
from pyNastran.gui.qt_files.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from pyNastran.gui.styles.trackball_style_camera import TrackballStyleCamera

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    from pyNastran.gui.dev.gui2.gui2 import MainWindow2


class ScalarBar:
    def __init__(self):
        pass
    def Modified(self):
        pass
    def VisibilityOn(self):
        pass
    def VisibilityOff(self):
        pass


class VtkInterface:
    def __init__(self, gui: MainWindow2, parent):
        self.gui = gui
        self.scalar_bar_actor = ScalarBar()


        self.vtk_interactor = QVTKRenderWindowInteractor(parent=parent)
        self.set_style_as_trackball()

        self.rend = vtkRenderer()
        #self.vtk_interactor.GetRenderWindow().AddRenderer(self.rend)
        fill_render_window(self.vtk_interactor, self.rend, nframes=1)

        camera = self.rend.GetActiveCamera()
        if self.settings.use_parallel_projection:
            camera.ParallelProjectionOn()
        #else:
            #camera.ParallelProjectionOff()

    @property
    def settings(self) -> Settings:
        return self.gui.settings

    @property
    def log(self):
        return self.gui.log

    def set_style_as_trackball(self):
        """sets the default rotation style"""
        #self._simulate_key_press('t') # change mouse style to trackball
        self.style = TrackballStyleCamera(self.vtk_interactor, self)
        self.vtk_interactor.SetInteractorStyle(self.style)

    def set_quad_grid(self, box_name: str,
                      nodes: np.ndarray, elements: np.ndarray,
                      color: Optional[list[float]]=None,
                      line_width: float=1, opacity: float=1.) -> None:
        if color is None:
            color = BLACK_FLOAT
        self.log.warning('set_quad_grid')

    def create_global_axes(self, dim_max: float) -> None:
        self.log.warning('create_global_axes')


def fill_render_window(vtk_interactor,
                       rend: vtkRenderer,
                       nframes: int=1) -> list[vtkRenderer]:
    assert nframes in [1, 2, 4], nframes

    render_window = vtk_interactor.GetRenderWindow()
    if nframes == 1:
        render_window.AddRenderer(rend)
        return [rend]

    if nframes == 2:
        # +-----+-----+
        # |     |     |
        # |  A  |  B  |
        # |     |     |
        # +-----+-----+
        # xmin, xmax, ymin, ymax
        #
        # xmin, ymin, xmax, ymax
        frame1 = [0.0, 0.0, 0.5, 1.0]
        frame2 = [0.5, 0.0, 1.0, 1.0]
    elif nframes == 4:
        # +-----+-----+
        # |     |     |
        # |  C  |  D  |
        # |     |     |
        # +-----+-----+
        # |     |     |
        # |  A  |  B  |
        # |     |     |
        # +-----+-----+
        frame1 = [0.0, 0.0, 0.5, 0.5]
        frame2 = [0.5, 0.0, 1.0, 0.5]
        frame3 = [0.5, 0.5, 1.0, 1.0]
        frame4 = [0.5, 0.5, 1.0, 1.0]
    else:
        raise ValueError(nframes)

    if nframes > 1:
        rend.SetViewport(*frame1)
    render_window.AddRenderer(rend)

    if nframes == 2:
        rend2 = vtkRenderer()
        rend.SetViewport(*frame2)
        render_window.AddRenderer(rend2)
        return [rend, rend2]
    elif nframes == 4:
        rend2 = vtkRenderer()
        rend3 = vtkRenderer()
        rend4 = vtkRenderer()
        rend2.SetViewport(*frame2)
        rend3.SetViewport(*frame3)
        rend4.SetViewport(*frame4)
        render_window.AddRenderer(rend2)
        render_window.AddRenderer(rend3)
        render_window.AddRenderer(rend4)
        return [rend, rend2, rend3, rend4]
    raise ValueError(nframes)

