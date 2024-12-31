"""
done
----
easy:
 - preferences window can now hide
 - add control over the update time (not tested)
 - added icase to Preferences Menu
 - add label in lower left cornre
 - add fast way to:
    - select another case (k/l keys)
    - enable/disable fringe (a key for apply_fringe)
    - update nphase (n key)
    - update scale factor (s key)

medium:
 - add displacement fringe

hard:
minor:

not done
--------
easy:
TODO: nphase minus crashes
TODO: nphase +/- doesn't scale right (seems to drop phases, but doesn't)

medium:
TODO: fix vtk window show/hide
TODO: add easier way to select another case (GUI case form?)
TODO: delink icase_disp and icase_fringe
TODO: -> add strain energy fringe

hard:
TODO: show/hide set of elements by:
TODO:   - element ids (copy from femap, paste as list?)
TODO:   - property ids
TODO:   - material ids
TODO:   - combined

minor:
TODO: fix up k/l key swapping (not quite right)
TODO: make vtk font support unicode
TODO: fix weirdness with 0012 model (green lines)
TODO: disable rotational modes to have fewer results (are these on?)
"""
from __future__ import annotations
import os
from typing import Callable, Optional, Any, cast, TYPE_CHECKING
import numpy as np

from qtpy.QtWidgets import (
    QHBoxLayout, QMenu,
    QMainWindow, QDockWidget, QFrame, QToolBar,
)
from qtpy.QtCore import QTimer
from pyNastran.utils import PathLike
from pyNastran.dev.bdf_vectorized3.nastran_io3 import (
    Nastran3, DisplacementResults)
from pyNastran.f06.dev.flutter.actions_builder import (
    Actions, Action, build_menus)

from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.vtk_rendering_core import (
    vtkActor, vtkDataSetMapper, vtkRenderer,
    vtkRenderWindow, vtkTextActor)
from pyNastran.gui.qt_files.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from pyNastran.gui.utils.vtk.vtk_utils import (
    numpy_to_vtk_points)
from pyNastran.gui.utils.vtk.gui_utils import numpy_array_to_vtk_array

import pyNastran

from pyNastran.gui.gui_objects.settings import NastranSettings
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.gui_objects.settings import Settings
from pyNastran.gui.qt_files.view_actions import ViewActions
from pyNastran.gui.qt_files.tool_actions import ToolActions

from pyNastran.gui.dev.gui2.vtk_interface import VtkInterface, ScalarBar
from pyNastran.gui.utils.qt.qsettings import QSettingsLike2
from pyNastran.gui.styles.trackball_style_camera import TrackballStyleCamera


if TYPE_CHECKING:  # pragma: no cover
    from vtkmodules.vtkCommonDataModel import vtkCellData, vtkPointData
    #from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor
    #from pyNastran.gui.menus.results_sidebar import ResultsSidebar
    from pyNastran.gui.typing import ColorInt
    from pyNastran.f06.dev.flutter.vtk_window_object import VtkWindowObject


class VtkWindow(QMainWindow):
    def __init__(self,
                 gui_obj: VtkWindowObject,
                 gui: QMainWindow,
                 data: dict[str, Any],
                 bdf_filename: str, op2_filename: str=''):
        """
        Sets up the vtk window

        Parameters
        ----------
        gui: QMainWindow
            the real gui
        gui_obj: VtkWindowObject
            the interface object

        """
        super().__init__()
        self.gui = gui
        self.gui_obj = gui_obj
        self.log = gui.log
        self.run_vtk = True
        self.global_apply_fringe_plot = True
        self.global_scale_factor = 1.0

        icon_path = self.gui_obj.icon_path
        self._create_menus(icon_path)

        # avoid redoing the fringe
        self.apply_fringe_plot = True

        self.set_data(data)
        self.eid_maps = {}
        self.nid_maps = {}
        self.alt_grids = {}

        self.vtk_mappers = {}
        self.vtk_actors = {}
        self.geometry_actors = self.vtk_actors
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
        self.vtk_mappers['main'] = grid_mapper

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

        self.corner_text_actors = setup_text_actors(
            renderer, 10, num=3)

        # if make_glyphs:
        # renderer.AddActor(arrow_actor)

        #self.setCentralWidget(self.vtk_frame)
        self.grid = ugrid
        self._load_model(bdf_filename, op2_filename)
        renderer.ResetCamera()

        # Render again to set the correct view
        self.plot_deformation(icase=self.icase, dphase=0.0)

        self.start_animation_timer()
        self.render()
        self.show()

    #---------------------------------------------------------------------
    ## teardown
    def closeEvent(self, event):
        # print('closeEvent')
        # event.accept()
        self.gui_obj.on_close()
        # print('back...')
    #---------------------------------------------------------------------
    # basic setup
    def _create_menus(self, icon_path: str) -> None:
        actions_dict = {
            # 'cycle_results': 'L',
            # 'rcycle_results': 'k',
            'apply_fringe_plot': Action(name='apply_fringe_plot', text='apply_fringe_plot...', icon='',
                                        shortcut='a', func=self.on_apply_fringe_plot),
            'cycle_icase': Action(name='cycle_icase', text='cycle_icase...', icon='',
                                  shortcut='l', func=self.on_cycle_icase),
            'rcycle_icase': Action(name='rcycle_icase', text='rcycle_icase...', icon='',
                                   shortcut='k', func=self.on_rcycle_icase),
            'scale_up': Action(name='scale_up', text='scale_up...', icon='',
                               shortcut='S', func=self.on_scale_up),
            'scale_down': Action(name='scale_down', text='scale_down...', icon='',
                                 shortcut='Ctrl+S', func=self.on_scale_down),
            'nphase_up': Action(name='nphase_up', text='nphase_up...', icon='',
                                shortcut='N', func=self.on_nphase_up),
            'nphase_down': Action(name='nphase_down', text='scale_down...', icon='',
                                  shortcut='Ctrl+N', func=self.on_nphase_down),
        }
        actions_input = Actions(icon_path, actions_dict)  # , load_icon=False
        self.qactions = actions_input.build_qactions(self)

        self.menubar = self.menuBar()
        self.menu_hidden = self.menubar.addMenu('File')
        self.menu_hidden.menuAction().setVisible(False)
        hidden_tools = ('apply_fringe_plot',
                        'cycle_icase', 'rcycle_icase',
                        'scale_up', 'scale_down',
                        'nphase_up', 'nphase_down',
                        )
        menus_dict = {
            'hidden': (self.menu_hidden, hidden_tools)
        }
        build_menus(menus_dict, self.qactions)

    def on_nphase_up(self):
        """changes the number of frames"""
        print('on_nphase_up')
        nphase = self.nphase + 1
        self.gui_obj.set_preferences(nphase=nphase)
        #self.update_dphases()

    def on_nphase_down(self):
        """changes the number of frames"""
        print('on_nphase_down')
        # self.nphase -= 1
        # if self.nphase == 0:
        #     self.nphase -= 1
        # self.update_dphases()
        nphase = self.nphase - 1
        if nphase == 0:
            return
            #nphase -= 1
        self.gui_obj.set_preferences(nphase=nphase)

    def on_scale_up(self) -> None:
        """changes the global scale factor"""
        self.global_scale_factor *= 1.1
    def on_scale_down(self) -> None:
        """changes the global scale factor"""
        self.global_scale_factor /= 1.1
    def on_rcycle_icase(self) -> None:
        icase = self._check_icase(self.icase - 1)
        self.gui_obj.set_preferences(icase=icase)

    def on_apply_fringe_plot(self):
        self.global_apply_fringe_plot = not self.global_apply_fringe_plot
        ugrid: vtkUnstructuredGrid = self.grid

        case, case_tuple = self.result_cases[self.icase]
        case = cast(DisplacementResults, case)
        i, name = case_tuple
        self.iphase = 0
        self.apply_fringe_plot = True
        self._update_fringe(i, name, case, ugrid)

    def on_cycle_icase(self) -> None:
        icase = self._check_icase(self.icase+1)
        self.gui_obj.set_preferences(icase=icase)

    def set_data(self, data: dict[str, int]) -> None:
        #print('setting data', data)
        self.dt_ms = data['dt_ms']
        self.iphase = 0
        self.nphase = data['nphase']

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
    #---------------------------------------------------------------------
    def stop_animation_timer(self):
        self.timer.stop()

    def update_dphases(self):
        num = self.nphase + 1
        self.dphases = np.linspace(0.0, 360.0, num=num)[:-1]
        print(num, self.dphases)
        self.scales = np.cos(np.radians(self.dphases))

    def start_animation_timer(self) -> None:
        if not hasattr(self, 'timer'):
            self.timer = QTimer()
        self.update_dphases()
        nphase = len(self.dphases)
        def update_vtk():
            if self.iphase == nphase:
                self.iphase = 0
            dphase = self.dphases[self.iphase]
            self.plot_deformation(icase=self.icase, dphase=dphase)
            self.render()
            self.iphase += 1
        self.timer.timeout.connect(update_vtk)
        self.timer.start(self.dt_ms)

    def _check_icase(self, icase: int) -> int:
        ncase = len(self.result_cases)
        if icase >= ncase:
            icase = 0
            self.icase = 0
            self.gui_obj.reset_icase_ncase(icase, ncase)
        elif icase < 0:
            print(f'  icase1={icase}')
            icase = ncase + icase
            print(f'  icase2={icase}')
            self.icase = icase
            self.gui_obj.reset_icase_ncase(icase, ncase)
        return icase

    def plot_deformation(self, icase: int, dphase: float) -> None:
        #ncase = max(self.result_cases)
        #print(f'plot_deformation; icase={icase}/{ncase}; dphase={dphase}')
        case, case_tuple = self.result_cases[icase]
        #print('found case')
        case = cast(DisplacementResults, case)
        # print(case)
        i, name = case_tuple
        scale = case.get_scale(i, name) * self.global_scale_factor
        phase = case.get_phase(i, name)
        #print(f'  scale={scale}; phase={phase}')
        if phase is None:
            # simple sinusoid
            scale = scale * self.scales[self.iphase]
        else:
            phase += dphase
        #print(f'phase={phase}; scale={scale:g}')
        xyz, deflected_xyz = case.get_vector_result_by_scale_phase(
            i, name, scale, phase=phase)
        ugrid = self.grid
        self._update_fringe(i, name, case, ugrid)

        #z = deflected_xyz[:, 2]
        #print(f'{z.min():g}, {z.max():g}')
        #print(case)
        self.is_deflected = True
        self._update_grid(ugrid, deflected_xyz)

    def _update_fringe(self, i: int, name: str,
                       case: DisplacementResults,
                       ugrid: vtkUnstructuredGrid) -> None:
        """
        Apply the fringe to the model only if:
         - global_apply_fringe_plot is True
        otherwise clear it.

        To reduce work, updates only happen when:
         - iphase = 0
         - apply_fringe_plot is True

        """
        if self.iphase != 0:
            return
        if not self.apply_fringe_plot:
            return

        # update corner text
        subcase = case.subcase_id[0]
        is_complex = case.is_complex
        complex_real_word = 'Complex' if is_complex else 'Real'
        legend_title = case.get_legend_title(i, name)
        header = case.headers[i]
        texts = [
            f'Subcase {subcase}',
            f'{complex_real_word} {legend_title}',
            header]
        set_corner_text(self.corner_text_actors, texts)

        vtk_point_data: vtkPointData = ugrid.GetPointData()
        if not self.global_apply_fringe_plot:
            vtk_point_data.SetScalars(None)
            vtk_point_data.Modified()
            self.apply_fringe_plot = False
            return

        fringe = case.get_fringe_result(i, name)
        dfringe = fringe.max() - fringe.min()
        dfringe = 1.0 if dfringe == 0.0 else dfringe
        fringe = (fringe - fringe.min()) / dfringe

        grid_mapper = self.vtk_mappers['main']
        vtk_fringe = numpy_array_to_vtk_array(
            grid_mapper, fringe, vector_size=1, phase=0)

        #vtk_point_data: vtkPointData = ugrid.GetPointData()
        vtk_point_data.SetScalars(vtk_fringe)
        vtk_point_data.Modified()
        self.apply_fringe_plot = False

    def _update_grid(self, grid: vtkUnstructuredGrid,
                     nodes: np.ndarray) -> None:
        vtk_points = grid.GetPoints()
        #inan = np.where(nodes.ravel() == np.nan)[0]
        #if len(inan) > 0:
            #raise RuntimeError('nan in nodes...')
        vtk_points = numpy_to_vtk_points(nodes, points=vtk_points,
                                         dtype='<f', deep=1)
        grid.SetPoints(vtk_points)
        grid.Modified()
        #self.grid_selected.Modified()
        #self._update_follower_grids(nodes)
        #self._update_follower_grids_complex(nodes)

    def _update_settings(self):
        """
        We took the pyNastranGUI settings and are hacking on them
        to speed up loading, minimize RAM usage, and focus the software.
        """
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
        """creates result_cases"""
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
            # creates result_cases
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
    def _finish_results_io2(self, name: str, form: list,
                            cases: dict[int, Any]):
        if self.is_geom:
            self.result_cases = {}
            # self.result_cases = cases
            return
        for case_id, case in cases.items():
            print(case)
        self.gui_obj.ncase = len(cases)

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

def setup_text_actors(renderer: vtkRenderer,
                      size: int, num: int=4):
    """
    Parameters
    ----------
    size: int
       the size of the text
    num: int; default=4
       the number of rows

    Returns
    -------
    vtk_text_actors : list[vtkTextActor]
       the actors

    Title:    ???
    Label:    ???
    Subtitle: ???
    Result:   Eigenvector Txyz
    Header:   mode=1 eigr=-0.0982041 eigi=17.2172 f=2.74021 Hz ζ=-0.0114076
    """
    vtk_text_actors = []
    text_size = 25
    #color = (1., 1., 1.)  # WHITE
    color = (0., 0., 0.)  # BLACK
    dtext_size = text_size + 5
    for i in range(num):
        text_actor = vtkTextActor()
        text_prop = text_actor.GetTextProperty()
        text_prop.SetFontSize(text_size)
        text_prop.SetBold(True)
        text_prop.SetColor(color)

        position = [5, 5 + (num-i-1) * dtext_size]
        text_actor.SetDisplayPosition(*position)

        vtk_text_actors.append(text_actor)
        renderer.AddActor(text_actor)

    texts = ['Subcase', 'Result: ', 'Note; ']
    assert len(texts) == num, (texts, num)
    set_corner_text(vtk_text_actors, texts)
    return vtk_text_actors

def set_corner_text(vtk_text_actors: list[vtkTextActor],
                    texts: list[str]) -> None:
    #print('set_corner_text', texts)
    for text, text_actor in zip(texts, vtk_text_actors):
        text_actor.SetInput(text)
