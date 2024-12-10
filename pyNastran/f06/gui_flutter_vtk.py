from __future__ import annotations
import os
from pathlib import Path
from typing import Callable, Optional, Any, TYPE_CHECKING

from qtpy.QtWidgets import (
    QMessageBox, QHBoxLayout,
    QMainWindow, QDockWidget, QFrame, QToolBar,
)

from pyNastran.dev.bdf_vectorized3.nastran_io3 import Nastran3
from pyNastran.gui.vtk_rendering_core import (
    vtkActor, vtkDataSetMapper, vtkRenderer,
    vtkRenderWindow, vtkRenderWindowInteractor)
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.qt_files.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

import pyNastran

from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.gui_objects.settings import Settings
from pyNastran.gui.qt_files.view_actions import ViewActions
from pyNastran.gui.qt_files.tool_actions import ToolActions

from pyNastran.gui.dev.gui2.vtk_interface import VtkInterface, ScalarBar
from pyNastran.gui.utils.qt.qsettings import QSettingsLike2
if TYPE_CHECKING:  # pragma: no cover
    #from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor
    #from pyNastran.gui.menus.results_sidebar import ResultsSidebar
    from pyNastran.gui.typing import ColorInt


from pyNastran.gui.styles.trackball_style_camera import TrackballStyleCamera
#RED_FLOAT = (1.0, 0., 0.)

PKG_PATH = Path(pyNastran.__path__[0])
#AERO_PATH = PKG_PATH / '..' / 'models' / 'aero'
#assert os.path.exists(AERO_PATH), AERO_PATH
#BDF_FILENAME = AERO_PATH / 'flutter_bug' / 'msc' / 'wing_b1.bdf'

class NewWindow(QMainWindow): # was QFame
    def __init__(self, bdf_filename: str):
        super().__init__()
        self.run_vtk = True
        self.eid_maps = {}
        self.nid_maps = {}
        self.alt_grids = {}
        self.geometry_actors = {}
        self.geometry_properties = {}
        self.settings = Settings(self)
        qsettings = QSettingsLike2()
        if hasattr(qsettings, 'load_json'):
            qsettings.load_json()

        self.scalar_bar_actor = ScalarBar()
        #self.load_actions = LoadActions(self)
        self.view_actions = ViewActions(self)
        self.tool_actions = ToolActions(self)

        self.setWindowTitle("New Window")

        if 1:
            self.vtk_frame = QFrame(self)

            self.vtk_interface = VtkInterface(self, self.vtk_frame)

            # put the vtk_interactor inside the vtk_frame
            self.set_vtk_frame_style()
            renderer = self.vtk_interface.rend
        else:
            # Setup renderer
            renderer = vtkRenderer()

            # Setup render window
            render_window = vtkRenderWindow()
            render_window.AddRenderer(renderer)

            # Setup render window interactor
            #render_window_interactor = vtkRenderWindowInteractor()
            render_window_interactor = QVTKRenderWindowInteractor(parent=self)

            # Render and start interaction
            render_window_interactor.SetRenderWindow(render_window)
            #render_window_interactor.Initialize()

            render_window_interactor.Start()
            #self.iren = render_window_interactor
            self.vtk_interactor = render_window_interactor
            self.setCentralWidget(render_window_interactor)
            self.set_style_as_trackball(render_window_interactor)
            self.rend = renderer

        # put the corner axis into the renderer
        self.tool_actions.create_corner_axis()

        grid_mapper = vtkDataSetMapper()
        ugrid = vtkUnstructuredGrid()
        grid_mapper.SetInputData(ugrid)

        geom_actor = vtkActor()
        geom_actor.SetMapper(grid_mapper)

        renderer.SetBackground(*self.settings.background_color)
        renderer.AddActor(geom_actor)
        # if make_glyphs:
        # renderer.AddActor(arrow_actor)
        renderer.ResetCamera()

        #self.setCentralWidget(self.vtk_frame)
        #self.setMinimumSize(400, 400)
        self.grid = ugrid

        #assert BDF_FILENAME.exists(), BDF_FILENAME
        bdf_filename = str(bdf_filename)
        analysis = Nastran3(self)
        analysis.load_nastran3_geometry(bdf_filename)

        # Render again to set the correct view
        self.render()

    def set_style_as_trackball(self, vtk_interactor) -> None:
        """sets the default rotation style"""
        #self._simulate_key_press('t') # change mouse style to trackball
        self.style = TrackballStyleCamera(vtk_interactor, self)
        vtk_interactor.SetInteractorStyle(self.style)

    def set_vtk_frame_style(self):
        """uses the vtk objects to set up the window (frame)"""
        vtk_hbox = QHBoxLayout()
        vtk_hbox.setContentsMargins(2, 2, 2, 2)

        vtk_hbox.addWidget(self.vtk_interactor)
        vtk_frame = self.vtk_frame
        vtk_frame.setLayout(vtk_hbox)
        vtk_frame.setFrameStyle(QFrame.NoFrame | QFrame.Plain)
        # this is our main, 'central' widget
        self.setCentralWidget(self.vtk_frame)
        #print('build_vtk_frame')

    # @property
    # def vtk_interactor(self) -> QVTKRenderWindowInteractor:
    #     return self.vtk_interface.vtk_interactor

    @property
    def vtk_interactor(self) -> QVTKRenderWindowInteractor:
        return self.vtk_interface.vtk_interactor
    @property
    def rend(self) -> vtkRenderer:
        return self.vtk_interface.rend
    @property
    def iren(self) -> QVTKRenderWindowInteractor:
        return self.vtk_interface.vtk_interactor
    @property
    def render_window(self) -> vtkRenderWindow:
        return self.vtk_interactor.GetRenderWindow()
    def render(self) -> None:
        self.vtk_interactor.GetRenderWindow().Render()

    def log_info(self, msg: str) -> None:
        print(msg)
    def log_debug(self, msg: str) -> None:
        print(msg)
    def create_global_axes(self, dim_max: float):
        return
    def _add_alt_actors(self,
                        grids_dict: dict[str, vtkUnstructuredGrid],
                        names_to_ignore=None):
        """
        Parameters
        ----------
        ignore_names : list[str]; default=None -> [main]
            add the actors to
        """
        if names_to_ignore is None:
            names_to_ignore = ['main']

        names = set(list(grids_dict.keys()))
        names_old = set(list(self.geometry_actors.keys()))
        names_old = names_old - set(names_to_ignore)
        #print('names_old1 =', names_old)

        #names_to_clear = names_old - names
        #self._remove_alt_actors(names_to_clear)
        #print('names_old2 =', names_old)
        #print('names =', names)
        for name in names:
            grid = grids_dict[name]
            self.tool_actions.add_alt_geometry(grid, name)

    def get_new_icase(self):
        return 0
    def _finish_results_io2(self, name, form, cases):
        pass
    def create_alternate_vtk_grid(self, name: str,
                                  color: ColorInt=None,
                                  line_width: int=5,
                                  opacity: float=1.0,
                                  point_size: int=1,
                                  bar_scale: float=0.0,
                                  representation: Optional[str]=None,
                                  display: Optional[str]=None,
                                  is_visible: bool=True,
                                  follower_nodes=None,
                                  follower_function: Optional[Callable]=None,
                                  is_pickable: bool=False,
                                  ugrid: vtkUnstructuredGrid=None,
                                  visible_in_geometry_properties: bool=True,
                                  ) -> None:
        """
        Creates an AltGeometry object

        Parameters
        ----------
        line_width : int
            the width of the line for 'surface' and 'main'
        color : ColorInt
            the RGB colors as ints
        opacity : float
            0.0 -> solid
            1.0 -> transparent
        point_size : int
            the point size for 'point'
        bar_scale : float
            the scale for the CBAR / CBEAM elements
        representation : str
            main - change with main mesh
            wire - always wireframe
            point - always points
            surface - always surface
            bar - can use bar scale
        is_visible : bool; default=True
            is this actor currently visible
        is_pickable : bool; default=False
            can you pick a node/cell on this actor
        follower_nodes : list[int]
            the nodes that are brought along with a deflection
        follower_function : function
            a custom follower_node update function
        ugrid : vtkUnstructuredGrid(); default=None
            the grid object; one will be created that you can fill
            if None is passed in
        visible_in_geometry_properties : bool; default=True
            True: show up in ``Edit Geometry Properties`` menu
            False: don't show up

        """
        if ugrid is None:
            ugrid = vtkUnstructuredGrid()
        self.alt_grids[name] = ugrid

        if name not in self.geometry_properties:
            self.geometry_properties[name] = AltGeometry(
                self, name, color=color,
                line_width=line_width, opacity=opacity,
                point_size=point_size, bar_scale=bar_scale,
                representation=representation, display=display,
                is_visible=is_visible, is_pickable=is_pickable,
                visible_in_geometry_properties=visible_in_geometry_properties,
            )
        if follower_nodes is not None:
            self.follower_nodes[name] = follower_nodes
        if follower_function is not None:
            self.follower_functions[name] = follower_function


class ScalarBar:
    def __init__(self):
        pass
    def VisibilityOn(self):
        pass
    def VisibilityOff(self):
        pass
    def Modified(self):
        pass
