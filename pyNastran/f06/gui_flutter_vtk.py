"""
TODO: add control over the update time
TODO: add label in lower left cornre

TODO: add displacement fringe
TODO: disable rotational modes to have fewer results
TODO: add way to select another case (j/k keys or GUI case form?)
TODO: show/hide set of elements (copy from femap, paste as list?)
TODO: fix weirdness with 0012 model (green lines)
TODO: add strain energy fringe
"""
from __future__ import annotations
import os
from pathlib import Path
from typing import Callable, Optional, Any, cast, TYPE_CHECKING
import numpy as np

from qtpy.QtWidgets import (
    QHBoxLayout,
    QMainWindow, QDockWidget, QFrame, QToolBar,
)
from qtpy.QtCore import QTimer
#from pyNastran.utils import object_attributes
from pyNastran.utils import PathLike
from pyNastran.dev.bdf_vectorized3.nastran_io3 import (
    Nastran3, DisplacementResults)
from pyNastran.gui.vtk_rendering_core import (
    vtkActor, vtkDataSetMapper, vtkRenderer,
    vtkRenderWindow, vtkRenderWindowInteractor)
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.qt_files.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from pyNastran.gui.utils.vtk.vtk_utils import (
    numpy_to_vtk_points)

import pyNastran

from pyNastran.gui.gui_objects.settings import NastranSettings
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


class VtkWindow(QMainWindow):
    def __init__(self, parent: QMainWindow,
                 bdf_filename: str, op2_filename: str=''):
        super().__init__()
        self.gui = parent
        self.log = parent.log
        self.run_vtk = True
        self.eid_maps = {}
        self.nid_maps = {}
        self.alt_grids = {}
        self.geometry_actors = {}
        self.geometry_properties = {}
        self.result_cases = {}
        self.icase = 0
        self.is_deflected = False

        self.settings = Settings(self)
        qsettings = QSettingsLike2()
        if hasattr(qsettings, 'load_json'):
            qsettings.load_json()

        self.scalar_bar_actor = ScalarBar()
        #self.load_actions = LoadActions(self)
        self.view_actions = ViewActions(self)
        self.tool_actions = ToolActions(self)

        title = f'FlutterGui - {bdf_filename}'
        if op2_filename:
            title += f' - {op2_filename}'
        self.setWindowTitle(title)

        #---------------------------------------------------
        self.vtk_frame = QFrame(self)

        self.vtk_interface = VtkInterface(self, self.vtk_frame)

        # put the vtk_interactor inside the vtk_frame
        self.set_vtk_frame_style()
        renderer = self.vtk_interface.rend

        # put the corner axis into the renderer
        self.tool_actions.create_corner_axis()

        grid_mapper = vtkDataSetMapper()
        ugrid = vtkUnstructuredGrid()
        grid_mapper.SetInputData(ugrid)

        geom_actor = vtkActor()
        geom_actor.SetMapper(grid_mapper)

        camera = renderer.GetActiveCamera()
        self._update_settings()
        if self.settings.use_parallel_projection:
            camera.ParallelProjectionOn()

        renderer.SetBackground(*self.settings.background_color)
        renderer.AddActor(geom_actor)
        # if make_glyphs:
        # renderer.AddActor(arrow_actor)

        #self.setCentralWidget(self.vtk_frame)
        self.grid = ugrid
        self._load_model(bdf_filename, op2_filename)
        renderer.ResetCamera()

        # Render again to set the correct view
        self.dphase = 0.0
        self.plot_deformation(icase=self.icase, dphase=self.dphase)

        self.start_animation_timer()
        self.render()

    def stop_animation_timer(self):
        self.timer.stop()

    def setup_animation_timer(self) -> None:
        self.iphase = 0
        self.nphase = 10
        # Update every N milliseconds
        self.dt_ms = 100  # ms
        self.timer = QTimer()

    def start_animation_timer(self) -> None:
        if not hasattr(self, 'timer'):
            self.setup_animation_timer()
        dphases = np.linspace(0.0, 360.0, num=self.nphase + 1)[:-1]
        nphase = len(dphases)
        def update_vtk():
            if self.iphase == nphase:
                self.iphase = 0
            dphase = dphases[self.iphase]
            self.plot_deformation(icase=self.icase, dphase=dphase)
            self.render()
            self.iphase += 1
        self.timer.timeout.connect(update_vtk)
        self.timer.start(self.dt_ms)

    def plot_deformation(self, icase: int, dphase: int) -> None:
        case, case_tuple = self.result_cases[icase]
        case = cast(DisplacementResults, case)
        i, name = case_tuple
        scale = case.get_scale(i, name)
        phase = case.get_phase(i, name) + dphase
        #print(f'phase={phase}; scale={scale:g}')
        xyz, deflected_xyz = case.get_vector_result_by_scale_phase(
            i, name, scale, phase=phase)
        z = deflected_xyz[:, 2]
        #print(f'{z.min():g}, {z.max():g}')
        #print(case)
        self.is_deflected = True
        self._update_grid(self.grid, deflected_xyz)

    def _update_grid(self, grid: vtkUnstructuredGrid,
                     nodes: np.ndarray) -> None:
        vtk_points = grid.GetPoints()
        #inan = np.where(nodes.ravel() == np.nan)[0]
        #if len(inan) > 0:
            #raise RuntimeError('nan in nodes...')
        vtk_points = numpy_to_vtk_points(nodes, points=vtk_points, dtype='<f', deep=1)
        grid.SetPoints(vtk_points)
        grid.Modified()
        #self.grid_selected.Modified()
        #self._update_follower_grids(nodes)
        #self._update_follower_grids_complex(nodes)

    def _update_settings(self):
        """we took the pyNastranGUI settings and are hacking on them"""
        nastran_settings: NastranSettings = self.settings.nastran_settings
        nastran_settings.is_bar_axes = False
        nastran_settings.is_shell_mcids = False
        nastran_settings.is_element_quality = False
        nastran_settings.is_rbe = False
        nastran_settings.is_constraints = False
        nastran_settings.is_3d_bars = False
        nastran_settings.is_3d_bars_update = False
        nastran_settings.is_mass_update = False
        nastran_settings.is_aero = False
        return

    def _load_model(self, bdf_filename: PathLike,
                   op2_filename: PathLike='') -> None:
        bdf_filename = str(bdf_filename)
        analysis = Nastran3(self)
        analysis.save_results_model = True
        self.is_geom = True

        #['bar_eids', 'bar_lines', 'card_index', 'data_map', 'element_cards', 'element_id', 'gui', 'gui_elements',
        # 'include_mass_in_geometry', 'mean_edge_length', 'model', 'node_id', 'save_results_model', 'xyz_cid0']
        analysis.load_nastran3_geometry(bdf_filename)
        #print(object_attributes(analysis))

        # self.inormal = -1
        # for key, (case, case_tag) in self.result_cases.items():
        #     print(case_tag)

        self.is_geom = False
        if os.path.exists(op2_filename):
            analysis.load_nastran3_results(op2_filename)
            for key, (case, case_tag) in self.result_cases.items():
                # if cae.titles[0] == ['Eigenvector']
                i, name = case_tag
                print(f'{key}: is_complex={case.is_complex} headers={case.headers[i]!r}')
                # print(case_tag)
            self.icase = 0  # plunge
            #self.icase = 1  # pitch
            #self.icase = 2  # flutter? pitch-plunge
            #self.icase = 3 # pitch
            return

        if len(analysis.model.eigenvectors) == 0:
            self.gui.mode2_pulldown.clear()
            return
        eigs = analysis.model.eigenvectors
        self.log.info(f'analysis.model.eigenvectors = {eigs}')

        #(1, 2, 1, 0, 0, '', '')
        #(1, 9, 1, 0, 0, '', '')
        #(subcase, flag, 1, 0, 0, '', '')
        subcase_dict = {}
        for key, case in eigs.items():
            subcase = key[0]
            flag = key[1]
            if subcase not in subcase_dict:
                subcase_dict[subcase] = {}
            if flag not in subcase_dict[subcase]:
                subcase_dict[subcase][flag] = (key, case)
            print(key)

        subcase_dict2 = {}
        for subcase, flag_dict in subcase_dict.items():
            if 2 in flag_dict and 9 in flag_dict:
                subcase_dict2[subcase] = flag_dict
            else:
                keys = list(flag_dict.keys())
                self.log.info(f'skipping flag_dict because 2/9 not in keys; keys={keys}')
        print(subcase_dict2)
        subcase = 1
        subcases = list(self.gui.responses)
        self.log.info(f'subcases = {subcases}')
        response = self.gui.responses[subcase]
        eigenvector_mpfs = response.eigenvectors
        self.log.info(f'eigenvector_mpfs = {eigenvector_mpfs}')
        asdf

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
        self.gui.log_info(msg)
    def log_debug(self, msg: str) -> None:
        self.gui.log_debug(msg)
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
    def get_form(self) -> list:
        return []
    def _finish_results_io2(self, name: str, form: list, cases: dict[int, Any]):
        if self.is_geom:
            self.result_cases = {}
            # self.result_cases = cases
            return
        for case_id, case in cases.items():
            print(case)
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
