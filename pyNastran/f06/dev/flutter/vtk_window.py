"""
easy:
 - preferences window can now hide
 - control over the update time (TODO: untested)
 - icase to Preferences Menu
 - edit_geometry_properties menu
   - TODO: some actors are busted and crash the gui (e.g., no main, show/hide masses)
 - label in lower left corner
 - fast way to:
    - select another case (k/l keys)
    - enable/disable fringe (a key for apply_fringe)
    - update nphase (n key)
    - update scale factor (u key for up)
    - set the model to points (d key)
      - change point size (i key for increase)

medium:
 - nodal displacement fringe
   - TODO: enable individual components
 - logic for creating a pickle file when the bdf
   is loaded (model.obj)
TODO: add export gif for set of views and modes
 - {model}_top_subcase-1_real-mode_f=-100 Hz_g=0.05.png
TODO: fix vtk window show/hide
TODO: add easier way to select another case (GUI case form?)
TODO: delink icase_disp and icase_fringe
TODO: -> add strain energy centroidal fringe

hard:
 - groups menu (looks like FEMAP):
   - picking:
     - one at a time
     - TODO: support multi-picking and turn off the checks after flipped
   - checkbox list for each material/property type
     - comments (TODO: support FEMAP syntax -> split by FEMAP property syntax)
   - checkbox list for each INCLUDE file
     - ifile logic done in shells/bars; warnings for others
     - intention is only elements use ifile logic

 - groups (still preliminary)
     - logic for user setting ifile_filter -> on_update_groups with g key
   - TODO: Use a vtkFilter for minimum code solution
     - Much harder code (more buggy) though, but potentially faster...
   - "Reload" the model for simplicity?
     - It's actually semi-decent...
     - it doesn't work cause it's a bad test case, but it's close...
     - I think CELAS2 from the 0012 are skipped
     - TODO: The code doesn't handle empty elements -> no update...

    - show/hide set of elements by:
        - TODO: element type (CTRIA3, CQUAD4, CBAR, ..., no RBEs/CONM2 since no results)
        - TODO: element ids (copy from femap, paste as list?)
        - TODO: ifile (untested)
        - property ids
        - TODO: determine properties based on materials
        - combined???

 - TODO: good picking
   - can I do a highlight an element when picked?
   - what about a window that only shows a the active element?
   - what about the info menu from femap?
   - currently have p key that will give element/property/nodes assoicated with an element
      - assumes it's only GRID
   - TODO: be able to pick bars
   - TODO: be able to pick nodes

minor:
 - vtk corner text font supports unicode
 - TODO: disable timer if it's not doing anything

minor:
TODO: fix up k/l key swapping (not quite right)
TODO: fix RBE2 bug ib 0012 model (green lines)
TODO: disable rotational modes to have fewer results (are these on?)
"""
from __future__ import annotations
import os
import pickle
import traceback
from typing import Callable, Optional, Any, cast, TYPE_CHECKING
#from pyNastran.dev.op2_vectorized3.op2_hdf5 import OP2, OP2Geom
import numpy as np

from qtpy.QtWidgets import (
    QHBoxLayout, QVBoxLayout, QMenu,
    QMainWindow, QDockWidget, QFrame, QToolBar,
    QTableWidget, QTableWidgetItem,
    #QListWidget, QListWidgetItem,
    #QTreeWidget, QTreeWidgetItem,
    QTreeView
)
from qtpy.QtGui import QStandardItemModel, QStandardItem
from qtpy.QtCore import QTimer, Qt

from pyNastran.utils import PathLike
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.properties.bars import _bar_areaL
from pyNastran.dev.bdf_vectorized3.bdf import BDF
from pyNastran.dev.bdf_vectorized3.nastran_io3 import (
    Nastran3, DisplacementResults)
from pyNastran.f06.dev.flutter.actions_builder import (
    Actions, Action, build_menus)
from pyNastran.f06.dev.flutter.scalar_bar import ScalarBar

from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.vtk_rendering_core import (
    vtkActor, vtkDataSetMapper, vtkRenderer,
    vtkRenderWindow, vtkTextActor)
from pyNastran.gui.utils.qt.qsettings import QSettingsLike2
from pyNastran.gui.utils.qt.pydialog import QFloatEdit, make_font
from pyNastran.gui.utils.vtk.vtk_utils import (
    numpy_to_vtk_points)
from pyNastran.gui.utils.vtk.gui_utils import numpy_array_to_vtk_array

import pyNastran
from pyNastran.gui.menus.edit_geometry_properties.edit_geometry_properties_object import (
    EditGeometryPropertiesObject)

from pyNastran.gui.gui_objects.settings import NastranSettings
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.gui_objects.settings import Settings
from pyNastran.gui.qt_files.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from pyNastran.gui.qt_files.view_actions import ViewActions
from pyNastran.gui.qt_files.tool_actions import ToolActions
from pyNastran.gui.qt_files.tool_actions import set_vtk_property_to_unicode
from pyNastran.gui.styles.trackball_style_camera import TrackballStyleCamera

from pyNastran.gui.dev.gui2.vtk_interface import VtkInterface, ScalarBar

from pyNastran.gui import font_file

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
                 bdf_filename: str,
                 op2_filename: str=''):
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
        self.bdf_filename = bdf_filename
        self.op2_filename = op2_filename
        self.load_results = True
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.pid_to_ptype: dict[int, str] = {}

        self.vtk_tools = VtkTools(gui_obj, self)
        self.log = gui.log
        self.edit_geometry_properties_obj = EditGeometryPropertiesObject(self)

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
        self.follower_nodes = {}
        self.follower_functions = {}

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
        self.set_vtk_frame_style(self.vtk_interactor, self.vtk_frame)
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
        self.vtk_actors['main'] = geom_actor

        self.corner_text_actors = setup_text_actors(
            renderer, 10, num=3)

        # if make_glyphs:
        # renderer.AddActor(arrow_actor)

        #self.setCentralWidget(self.vtk_frame)
        self.grid = ugrid
        self._load_model(bdf_filename, op2_filename)
        self.case_keys = list(self.result_cases)
        renderer.ResetCamera()

        if self.load_results:
            # Render again to set the correct view
            self.update_dphases()
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
        vtk_tools = self.vtk_tools
        actions_dict = {
            # 'cycle_results': 'L',
            # 'rcycle_results': 'k',
            'setup_pick': Action(name='setup_pick', text='Pick...', icon='',
                                 shortcut='p', func=self.setup_pick, tip='pick'),
            'update_groups': Action(name='update_groups', text='Update Groups...', icon='',
                                   shortcut='g', func=self.on_update_groups, tip='Update the groups'),
            'apply_fringe_plot': Action(name='apply_fringe_plot', text='apply_fringe_plot...', icon='',
                                        shortcut='a', func=self.on_apply_fringe_plot, tip='Turn the fringe on/off'),
            'main_point': Action(name='main_point', text='Set to Points', icon='',
                                        shortcut='d', func=self.vtk_tools.on_main_point, tip='Sets the main actor to points'),
            'point_size_up': Action(name='point_size_up', text='Increase Point Size', icon='',
                                        shortcut='i', func=self.vtk_tools.on_point_size_up, tip='Increases the point size'),
            'point_size_down': Action(name='point_size_down', text='Decrease Point Size', icon='',
                                        shortcut='Ctrl+i', func=self.vtk_tools.on_point_size_down, tip='Decreases the point size'),
            'cycle_icase': Action(name='cycle_icase', text='cycle_icase', icon='',
                                  shortcut='l', func=vtk_tools.on_cycle_icase, tip='Increases icase'),
            'rcycle_icase': Action(name='rcycle_icase', text='rcycle_icase', icon='',
                                   shortcut='k', func=vtk_tools.on_rcycle_icase, tip='Decreases icase'),
            'scale_up': Action(name='scale_up', text='Increase Global Deformation Scale', icon='',
                               shortcut='u', func=vtk_tools.on_scale_up, tip='Increases the deformation scale factor by 1.1x'),
            'scale_down': Action(name='scale_down', text='Decrease Global Deformation Scale', icon='',
                                 shortcut='Ctrl+U', func=vtk_tools.on_scale_down, tip='Decreases the deformation scale factor by 1.1x'),
            'nphase_up': Action(name='nphase_up', text='nphase_up...', icon='',
                                shortcut='N', func=vtk_tools.on_nphase_up, tip='Increase the number of frame'),
            'nphase_down': Action(name='nphase_down', text='scale_down...', icon='',
                                  shortcut='Ctrl+N', func=vtk_tools.on_nphase_down, tip='Decreases the number of frame'),
            'geo_properties': Action(name='geo_properties', text='Edit Geometry Properties...',
                                     icon='', tip='Change Model Color/Opacity/Line Width',
                                     shortcut='Ctrl+E',
                                     func=self.edit_geometry_properties_obj.edit_geometry_properties,
                                     show=True),
        }
        actions_input = Actions(icon_path, actions_dict)  # , load_icon=False
        self.qactions = actions_input.build_qactions(self)

        self.menubar = self.menuBar()
        self.menu_hidden = self.menubar.addMenu('Hidden')
        self.menu_view = self.menubar.addMenu('View')
        # self.menu_hidden.menuAction().setVisible(False)
        view_tools = (
            'geo_properties', '',
            'main_point', 'point_size_up', 'point_size_down',
        )
        hidden_tools = (
            'apply_fringe_plot', 'update_groups',
            'cycle_icase', 'rcycle_icase',
            'scale_up', 'scale_down',
            'nphase_up', 'nphase_down',
            'setup_pick',
        )
        menus_dict = {
            'view': (self.menu_view, view_tools),
            'hidden': (self.menu_hidden, hidden_tools),
        }
        build_menus(menus_dict, self.qactions)

    def on_apply_fringe_plot(self):
        self.global_apply_fringe_plot = not self.global_apply_fringe_plot
        ugrid: vtkUnstructuredGrid = self.grid
        vtk_mapper = self.vtk_mappers['main']
        if len(self.result_cases) == 0:
            return
        case, case_tuple = self.result_cases[self.icase]
        case = cast(DisplacementResults, case)
        i, name = case_tuple
        self.iphase = 0
        self.apply_fringe_plot = True
        self.apply_fringe_plot = _update_fringe(
            i, name, case, ugrid, vtk_mapper,
            self.corner_text_actors,
            iphase=self.iphase,
            global_apply_fringe_plot=self.global_apply_fringe_plot,
            apply_fringe_plot=self.apply_fringe_plot)
        print('done')

    def set_data(self, data: dict[str, int]) -> None:
        #print('setting data', data)
        self.dt_ms = data['dt_ms']
        self.iphase = 0
        self.nphase = data['nphase']
        self.animate = True #data['animate']
        self.point_size = 2

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

    def log_debug(self, msg: str) -> None:
        self.gui.log_debug(msg)
    def log_info(self, msg: str) -> None:
        self.gui.log_info(msg)
    def log_command(self, msg: str) -> None:
        self.gui.log_command(msg)
    def log_warning(self, msg: str) -> None:
        self.gui.log_warning(msg)
    def log_error(self, msg: str) -> None:
        self.gui.log_error(msg)
    #---------------------------------------------------------------------
    def stop_animation_timer(self):
        self.timer.stop()

    def update_dphases(self):
        num = self.nphase + 1
        self.dphases = np.linspace(0.0, 360.0, num=num)[:-1]
        #print(num, self.dphases)
        self.scales = np.cos(np.radians(self.dphases))

    def start_animation_timer(self) -> None:
        if not hasattr(self, 'timer'):
            self.timer = QTimer()
        self.update_dphases()
        def update_vtk():
            nphase = len(self.dphases)
            if self.iphase == nphase:
                self.iphase = 0
            dphase = self.dphases[self.iphase]
            self.plot_deformation(icase=self.icase, dphase=dphase)
            self.render()
            self.iphase += 1
        self.timer.timeout.connect(update_vtk)
        self.timer.start(self.dt_ms)

    def plot_deformation(self, icase: int, dphase: float) -> None:
        #ncase = max(self.result_cases)
        #print(f'plot_deformation; icase={icase}/{ncase}; dphase={dphase}')
        if len(self.result_cases) == 0:
            return
        case, case_tuple = self.result_cases[icase]
        #print('found case')
        case = cast(DisplacementResults, case)
        i, name = case_tuple
        scale = case.get_scale(i, name) * self.global_scale_factor
        phase = case.get_phase(i, name)
        #print(f'  scale={scale}; phase={phase}')
        if phase is None:
            # simple sinusoid
            scale = scale * self.scales[self.iphase]
        else:
            phase = dphase
            #print('phase =', phase)
        #print(f'phase={phase}; scale={scale:g}')
        xyz, deflected_xyz = case.get_vector_result_by_scale_phase(
            i, name, scale, phase=phase)
        ugrid = self.grid
        vtk_mapper = self.vtk_mappers['main']
        _update_fringe(i, name, case, ugrid, vtk_mapper,
                       self.corner_text_actors,
                       iphase=self.iphase,
                       global_apply_fringe_plot=self.global_apply_fringe_plot,
                       apply_fringe_plot=self.apply_fringe_plot)

        #print(case)
        self.is_deflected = True
        _update_grid(ugrid, deflected_xyz)

    def _update_settings(self):
        """
        We took the pyNastranGUI settings and are hacking on them
        to speed up loading, minimize RAM usage, and focus the software.
        """
        nastran_settings: NastranSettings = self.settings.nastran_settings
        _update_nastran_settings(nastran_settings)

    def _load_model(self, bdf_filename: PathLike,
                   op2_filename: PathLike='') -> None:
        """creates result_cases"""
        bdf_filename = str(bdf_filename)
        log = self.gui.log

        self.analysis = Nastran3(self)
        self.analysis.save_results_model = True
        self.is_geom = True

        print('reload_model...')
        model = self._reload_model(
            bdf_filename, use_obj_file=False)
        self.fill_table_tree(model)

        # self.inormal = -1
        # for key, (case, case_tag) in self.result_cases.items():
        #     print(case_tag)

        self.op2_filename = None
        if not self.load_results:
            return
        if os.path.exists(op2_filename):
            self.op2_filename = op2_filename
            self._reload_results(op2_filename, log)
            return

        if len(self.analysis.model.eigenvectors) == 0:
            self.gui.mode2_pulldown.clear()
            return
        eigs = self.analysis.model.eigenvectors
        log.info(f'analysis.model.eigenvectors = {eigs}')

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

    def fill_table_tree(self, model: BDF) -> None:
        for card in model.element_cards:
            if card.n == 0:
                continue
            if not len(card.ifile) == card.n:
                self.gui.log_warning(f'{card.type} ifile not created')
        self.ifile_name_dict = get_ifile_name_dict(
            self.analysis.model)
        #fill_table_tree(self.table_tree, self.ifile_name_dict)

        tree = self.tree
        pids, properties, pid_to_type_name = get_property_table(model)
        mids, materials, mid_to_type_name = get_material_table(model)
        self.pid_to_type_name = pid_to_type_name
        self.mid_to_type_name = mid_to_type_name

        self.ifile_name_dict = get_ifile_name_dict(
            self.analysis.model)
        groups = []
        for ifile, fname in self.ifile_name_dict.items():
            groups.append((f'{ifile}: File: {fname}', True, ifile, []))

        words = [
            ('Properties', False, 0, properties),
            ('Materials', False, 1, materials),
            ('Groups', False, 2, groups),
        ]

        self.materials_filter_set = set(mids.tolist())
        self.properties_filter_set = set(pids.tolist())
        self.ifile_filter_set = set(list(self.ifile_name_dict))
        tree = self.tree
        self.tree_model = QStandardItemModel()
        self.tree_model.setHorizontalHeaderLabels(['Model'])
        tree_model_add_items(self.tree_model, words)
        tree.setModel(self.tree_model)
        self.tree_model.itemChanged.connect(self.on_tree_item_changed)

    def on_tree_item_changed(self, item: QStandardItem) -> None:
        if self._updating:
            return
        #texti, idi, typei, commenti = _qitem_text_to_sline(item.text())
        #print(f'id={idi}; type={typei!r} other={commenti!r}')

        is_checked = (item.checkState() == Qt.Checked)
        text0 = item.text()
        selected_items = self.tree.selectedIndexes()

        # disables autoupdating when the check marks are flipped
        self._updating = True
        for itemii in selected_items:
            text = self.tree_model.data(itemii)
            textii, idii, typeii, commentii = _qitem_text_to_sline(text)
            print(f'id={idii}; type={typeii!r} other={commentii!r}')

            if typeii.startswith('MAT'):  # MAT1, MAT8, ...
                myset = self.materials_filter_set
            elif typeii.startswith('P'):  # PSHELL, PCOMP, ...
                myset = self.properties_filter_set
            elif typeii == 'File':
                myset = self.ifile_filter_set
            else:  # pragma: no cover
                raise RuntimeError(f'type = {typeii!r}')

            skip_checkstate = (text0 == textii)
            if is_checked:
                print(f"{textii} is checked")
                #if skip_checkstate:
                    #item.setCheckState(Qt.Checked)
                    #item.setCheckState(Qt.CheckState.Checked)
                myset.add(idii)
            else:
                print(f"{textii} is unchecked")
                #if skip_checkstate:
                    #item.setCheckState(Qt.Unchecked)
                    #item.setCheckState(Qt.CheckState.Unchecked)
                try:
                    myset.remove(idii)
                except KeyError:  # can't remove a thing that doesn't exist
                    pass
        #print(f'{type}: {myset}')
        self._updating = False

    def on_update_groups(self) -> None:
        print('on_update_groups')
        print(self.table_tree)
        # ifile_filter = get_ifile_groups_from_table(
        #     self.table_tree, self.model)
        print('self.ifile_filter_set =', self.ifile_filter_set)
        ifile_filter = np.array(list(self.ifile_filter_set), dtype='int32')
        iprop_filter = np.array(list(self.properties_filter_set), dtype='int32')
        imat_filter = np.array(list(self.materials_filter_set), dtype='int32')
        ifile_filter.sort()
        iprop_filter.sort()
        imat_filter.sort()
        print(f'ifile_filter = {ifile_filter}')
        print(f'iprop_filter = {iprop_filter}')
        print(f'imat_filter = {imat_filter}')

        self.analysis.materials_filter = imat_filter
        self.analysis.properties_filter = iprop_filter
        self.analysis.ifile_filter = ifile_filter
        self._reload_model(self.bdf_filename, use_obj_file=True)
        if self.load_results:
            self._reload_results(self.op2_filename, self.gui.log)
        self.render()

    def setup_pick(self):
        # Create a picker
        from vtkmodules.vtkRenderingCore import vtkCellPicker #, vtkPointPicker, vtkAreaPicker, vtkDataSetMapper
        picker = vtkCellPicker()
        # Add the observer for the picking event
        # def on_pick(self):
        #     cell_id = picker.GetCellId()
        #     if cell_id != -1:
        #         print("Hovering over the sphere!; cell_id={cell_id}")
        # self.iren.AddObserver("LeftButtonPressEvent", on_pick)
        #self.iren.AddObserver("KeyPressEvent", on_pick)
        pixel_x, pixel_y = self.vtk_interactor.GetEventPosition()
        picker.Pick(pixel_x, pixel_y, 0, self.rend)
        cell_id = picker.GetCellId()
        print('cell_id = {cell_id}')

        if cell_id != -1:
            i = cell_id
            etype = self.analysis.element_types[i]
            eid = self.element_id[i]
            pid = self.property_id[i]
            #ptype = self.pid_to_ptype[pid]
            ptype = self.pid_to_type_name[pid][0]
            ptype_lower = ptype.lower()
            etype_lower = etype.lower()
            elem_str = ''
            prop_str = ''
            nid_str = ''
            log = self.log
            if not hasattr(self, 'analysis_model'):
                log.error('missing analysis model')
                return
            if hasattr(self.analysis_model, etype_lower):
                elem = getattr(self.analysis_model, etype_lower)
                try:
                    elemi = elem.slice_card_by_element_id([eid])
                    elem_str = elemi.write()
                except Exception:
                    self.log_error(str(traceback.format_exc()))
                    log.warning(f'eid={eid} is missing from {etype_lower}.element_id = {elem.element_id}')
            else:
                log.warning(f"no such element type {etype_lower}")
            if hasattr(self.analysis_model, ptype_lower):
                prop = getattr(self.analysis_model, ptype_lower)
                try:
                    assert len(prop.ifile) == len(prop.property_id)
                    propi = prop.slice_card_by_property_id([pid])
                    prop_str = propi.write()
                except Exception:
                    self.log_error(str(traceback.format_exc()))
                    # self.log_error('\n' + ''.join(traceback.format_stack()))
                    log.warning(f'pid={pid} is missing from {ptype_lower}.property_id = {prop.property_id}')
            else:
                log.warning(f"no such property type {ptype_lower}")

            nids = elemi.nodes.flatten()
            nids.sort()
            try:
                #self.analysis_model.grid.slice_card_by_id(nids)
                nidsi = self.analysis_model.grid.slice_card_by_node_id(nids)
                nid_str = nidsi.write()
            except Exception as e:
                print('failed on nodes...')
                pass
            log.info(f'i={i}: eid={etype} {eid}; pid={pid} {ptype}\n'
                     f'{nid_str}{elem_str}{prop_str}')

    def _reload_model(self, bdf_filename: PathLike,
                      use_obj_file: bool=False) -> None:
        log = self.gui.log
        #log.info('start of _reload_model')
        self.is_geom = True
        self.case_keys = []
        self.result_cases = {}

        bdf_filename = str(bdf_filename)
        obj_filename = bdf_filename + '.obj'
        bdf_filename_lower = bdf_filename.lower()
        #if use_obj_file and os.path.exists(obj_filename):
            #log.info(f'loading obj {obj_filename}')
            #model = read_obj(obj_filename)
            #model.log = log
        if use_obj_file:
            model = self.analysis_model
        # elif os.path.exists(obj_filename):
        #     log.info(f'  loading obj {obj_filename}')
        #     model = read_obj(obj_filename)
        #     model.log = log
        elif bdf_filename_lower.endswith('.op2'):
            op2_read
            #geo = self.load_op2_geometry(bdf_filename)
        elif bdf_filename_lower.endswith('.h5'):
            h5_read
            #geo = self.load_h5_geometry(bdf_filename)
        else:
            self.bdf_filename = bdf_filename
            model = self.analysis.get_bdf_geometry(
                bdf_filename, log=log)
            #model.load_mode = 'bdf'
            #if self.use_obj_file:
            # print(type(model))
            log.info(f'saving obj {obj_filename}')
            write_obj(model, obj_filename)
        self.analysis_model = model
        self.set_pid_to_ptype(model)

        if isinstance(model, BDF):
            self.analysis.load_nastran3_geometry(bdf_filename, name='main')
        # elif isinstance(model, OP2):
        #     raise RuntimeError(type(model))
        else:
            raise RuntimeError(type(model))
        self.is_geom = False
        #log.info('end of _reload_model')
        return model

    def set_pid_to_ptype(self, model: BDF) -> None:
        for prop in model.property_cards:
            if prop.n == 0:
                continue
            ptype = prop.type
            for pid in prop.property_id:
                self.pid_to_ptype[int(pid)] = str(ptype)

    def _reload_results(self, op2_filename, log):
        # creates result_cases
        log.level = 'info'
        self.analysis.load_nastran3_results(op2_filename)
        log.level = 'debug'
        for key, (case, case_tag) in self.result_cases.items():
            # if cae.titles[0] == ['Eigenvector']
            i, name = case_tag
            print(f'{key}: is_complex={case.is_complex} headers={case.headers[i]!r}')
            # print(case_tag)
        self.icase = 0
        return

    def set_style_as_trackball(self, vtk_interactor) -> None:
        """sets the default rotation style"""
        #self._simulate_key_press('t') # change mouse style to trackball
        self.style = TrackballStyleCamera(vtk_interactor, self)
        vtk_interactor.SetInteractorStyle(self.style)

    def set_vtk_frame_style(self, vtk_interactor: QVTKRenderWindowInteractor,
                            qt_frame: QFrame) -> None:
        """uses the vtk objects to set up the window (frame)"""
        vtk_hbox = QHBoxLayout()
        vtk_hbox.setContentsMargins(2, 2, 2, 2)

        vbox = QVBoxLayout()
        font_size = self.gui.font_size
        font = make_font(font_size)
        self.ifile_name_dict = {}

        if 0:
            table_tree = QTableWidget(self)
            table_tree.setRowCount(3)
            table_tree.setColumnCount(1)
            table_tree.setHorizontalHeaderLabels([self.tr('Groups')])
            table_tree.setFont(font)
        else:
            self.tree = QTreeView()
            self._updating = False
            table_tree = self.tree
            #table_tree.setSelectionMode(QTreeView.MultiSelection)

        vbox.addWidget(table_tree)
        self.table_tree = table_tree

        vtk_hbox.addWidget(vtk_interactor)
        vtk_hbox.addLayout(vbox)
        qt_frame.setLayout(vtk_hbox)
        qt_frame.setFrameStyle(QFrame.NoFrame | QFrame.Plain)
        # this is our main, 'central' widget
        self.setCentralWidget(qt_frame)
        self.set_font_size(self.gui.font_size)
    def set_font_size(self, font_size: int) -> None:
        font = make_font(font_size)
        # self.tree.setFont(font)
        self.setFont(font)

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

    def get_new_icase(self) -> int:
        return 0
    def get_form(self) -> list:
        return []
    def _finish_results_io2(self, name: str, form: list,
                            cases: dict[int, Any]):
        for i, (case, i_name) in cases.items():
            if len(i_name) == 2:
                i, name = i_name
                if name == 'ElementID':
                    self.element_id = case.get_fringe_result(*i_name)
                elif name == 'PropertyID':
                    self.property_id = case.get_fringe_result(*i_name)
            else:
                raise RuntimeError(i_name)
            #print(i_name)

        if self.is_geom:
            self.result_cases = {}
            # self.result_cases = cases
            return
        #for case_id, case in cases.items():
            #print(case)
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
        set_vtk_property_to_unicode(text_prop, font_file)
        text_prop.SetFontSize(text_size)
        text_prop.SetBold(True)
        text_prop.SetColor(color)

        position = [5, 5 + (num-i-1) * dtext_size]
        text_actor.SetDisplayPosition(*position)

        vtk_text_actors.append(text_actor)
        renderer.AddActor(text_actor)

    texts = ['Subcase', 'Result: ', 'Note: ']
    assert len(texts) == num, (texts, num)
    set_corner_text(vtk_text_actors, texts)
    return vtk_text_actors

def set_corner_text(vtk_text_actors: list[vtkTextActor],
                    texts: list[str]) -> None:
    #print('set_corner_text', texts)
    for text, text_actor in zip(texts, vtk_text_actors):
        text_actor.SetInput(text)

class VtkTools:
    def __init__(self,
                 gui_obj: VtkWindowObject,
                 vtk_window: VtkWindow):
        self.gui_obj = gui_obj
        self.vtk_window = vtk_window

    def on_nphase_up(self):
        """changes the number of frames"""
        nphase = self.vtk_window.nphase + 1
        self.gui_obj.set_preferences(nphase=nphase)

    def on_nphase_down(self):
        """changes the number of frames"""
        nphase = self.vtk_window.nphase - 1
        if nphase == 0:
            return
            #nphase -= 1
        self.gui_obj.set_preferences(nphase=nphase)

    def on_scale_up(self) -> None:
        """changes the global scale factor"""
        self.vtk_window.global_scale_factor *= 1.1
    def on_scale_down(self) -> None:
        """changes the global scale factor"""
        self.vtk_window.global_scale_factor /= 1.1

    def on_rcycle_icase(self) -> None:
        icase = self._check_icase(self.vtk_window.icase - 1)
        self.gui_obj.set_preferences(icase=icase)

    def on_point_size_up(self) -> None:
        point_size = self.vtk_window.point_size + 1
        self.set_point_size(point_size)

    def on_point_size_down(self) -> None:
        point_size = self.vtk_window.point_size - 1
        self.set_point_size(point_size)

    def set_point_size(self, point_size: int) -> None:
        point_size = max(point_size, 1)
        actor: vtkActor = self.vtk_window.geometry_actors['main']
        prop = actor.GetProperty()
        prop.SetPointSize(point_size)
        self.vtk_window.point_size = point_size

    def on_main_point(self):
        actor: vtkActor = self.vtk_window.geometry_actors['main']
        prop = actor.GetProperty()
        prop.SetRepresentationToPoints()
        self.vtk_window.render()

    def on_cycle_icase(self) -> None:
        icase = self._check_icase(self.vtk_window.icase+1)
        self.gui_obj.set_preferences(icase=icase)

    def _check_icase(self, icase: int) -> int:
        ncase = len(self.vtk_window.result_cases)
        if icase >= ncase:
            icase = 0
            self.vtk_window.icase = 0
            self.gui_obj.reset_icase_ncase(icase, ncase)
        elif icase < 0:
            icase = ncase + icase
            self.vtk_window.icase = icase
            self.gui_obj.reset_icase_ncase(icase, ncase)
        return icase

def _update_nastran_settings(nastran_settings: NastranSettings) -> None:
    """
    We took the pyNastranGUI settings and are hacking on them
    to speed up loading, minimize RAM usage, and focus this tool.
    """
    nastran_settings.is_mass = False
    nastran_settings.is_mass_update = False
    nastran_settings.is_shell_mcids = False
    nastran_settings.is_element_quality = False
    nastran_settings.is_rbe = False
    nastran_settings.is_constraints = False
    nastran_settings.is_3d_bars = False
    nastran_settings.is_3d_bars_update = False
    nastran_settings.is_bar_axes = False
    nastran_settings.is_aero = False

def _update_fringe(i: int, name: str,
                   case: DisplacementResults,
                   ugrid: vtkUnstructuredGrid,
                   vtk_mapper: vtkDataSetMapper,
                   corner_text_actors: list[vtkTextActor],
                   iphase: int,
                   global_apply_fringe_plot: bool,
                   apply_fringe_plot: bool) -> bool:
    """
    Apply the fringe to the model only if:
     - global_apply_fringe_plot is True
    otherwise clear it.

    To reduce work, updates only happen when:
     - iphase = 0
     - apply_fringe_plot is True

    """
    assert isinstance(vtk_mapper, vtkDataSetMapper), type(vtk_mapper)
    if iphase != 0:
        return apply_fringe_plot
    if not apply_fringe_plot:
        return apply_fringe_plot

    # update corner text
    if isinstance(case.subcase_id, integer_types):
        subcase = case.subcase_id
    else:
        subcase = case.subcase_id[0]
    is_complex = case.is_complex
    complex_real_word = 'Complex' if is_complex else 'Real'
    legend_title = case.get_legend_title(i, name)
    header = case.headers[i]
    texts = [
        f'Subcase {subcase}',
        f'{complex_real_word} {legend_title}',
        header]
    set_corner_text(corner_text_actors, texts)

    vtk_point_data: vtkPointData = ugrid.GetPointData()
    if not global_apply_fringe_plot:
        vtk_point_data.SetScalars(None)
        vtk_point_data.Modified()
        apply_fringe_plot = False
        return apply_fringe_plot

    fringe = case.get_fringe_result(i, name)
    dfringe = fringe.max() - fringe.min()
    dfringe = 1.0 if dfringe == 0.0 else dfringe
    fringe = (fringe - fringe.min()) / dfringe

    vtk_fringe = numpy_array_to_vtk_array(
        vtk_mapper, fringe, vector_size=1, phase=0)

    #vtk_point_data: vtkPointData = ugrid.GetPointData()
    vtk_point_data.SetScalars(vtk_fringe)
    vtk_point_data.Modified()
    apply_fringe_plot = False
    return apply_fringe_plot

def _update_grid(grid: vtkUnstructuredGrid,
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

def fill_table_tree(table: QTableWidget,
                    ifile_name_dict: dict[int, str],
                    right_clicked_func=None) -> None:
    data = []
    for i, name in ifile_name_dict.items():
        data.append((name, i, []))
    nrows = len(data)
    table.setMaximumWidth(250)
    table.setRowCount(nrows)
    table.setColumnCount(1)
    for irow, element in enumerate(data):
        text, i, children = element
        assert len(children) == 0, children
        item = QTableWidgetItem(text)
        table.setItem(irow, 0, item)
        item.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
        item.setCheckState(Qt.CheckState.Checked)
        if right_clicked_func is not None:
            func = partial(right_clicked_func, i)
            item.clicked.connect(func)

    # Resize all columns to fit their contents
    table.resizeColumnsToContents()
    return table

def tree_model_add_items(tree_model: QStandardItemModel,
                         elements: dict,
                         level: int=0,
                         count_check: bool=False) -> None:
    nelements = len(elements)
    redo = False
    for element in elements:
        if not len(element) == 4:
            print('element = %r' % str(element))

        text, is_checkable, i, children = element
        nchildren = len(children)
        # print('text=%r' % text)
        item = QStandardItem(text)
        item.setEditable(False)
        if is_checkable:
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked)
            #item.setCheckable(True)
        tree_model.appendRow(item)
        if nelements == 1 and nchildren == 0 and level == 0:
            # self.result_data_window.setEnabled(False)
            item.setEnabled(False)
            # print(dir(self.tree_view))
            # self.tree_view.setCurrentItem(self, 0)
            # item.mousePressEvent(None)
            redo = True
        # else:
        # pass
        # print('item=%s count_check=%s nelements=%s nchildren=%s' % (
        # text, count_check, nelements, nchildren))
        if children:
            assert isinstance(children, list), children
            tree_model_add_items(item, children, level + 1,
                                 count_check=count_check)
            # print('*children = %s' % children)
    is_single = redo
    return is_single


def get_ifile_name_dict(model: BDF) -> dict[int, str]:
    active_filenames = model.active_filenames
    ifile_name_dict = {}
    for ifile, fname in enumerate(active_filenames):
        ifile_name_dict[ifile] = os.path.basename(fname)
    return ifile_name_dict

def read_obj(obj_filename: PathLike):
    with open(obj_filename, 'rb') as obj_file:
        model = pickle.load(obj_file)
    return model

def write_obj(obj, obj_filename: PathLike):
    with open(obj_filename, 'wb') as obj_file:
        pickle.dump(obj, obj_file)


def get_ifile_groups_from_table(table_tree: QTreeView,
                                model: QStandardItemModel,
                                level: int=0) -> np.ndarray:
    print(f'get_ifile_groups_from_table; level={level}')
    nrows = model.rowCount()
    print('nrows =', nrows)
    icol = 0
    ifile_groups_root = {}
    ifile_groups = []
    text_flags = ['Properties', 'Materials', 'Groups']
    is_root = (level == 0)
    for irow in range(nrows):
        print(f'irow = {irow}')
        if level == 0:
            print('a')
            item = model.item(irow, icol)
        else:
            print('b')
            item = model.item(irow)
        print(f'item = {item}')
        texti = item.text()
        if is_root and texti in text_flags:
            ifile_groups_root[texti] = []
            out = get_ifile_groups_from_table(table_tree, item, level=level+1)
            print(f'root! text={texti!r}; out={out}')
        else:
            print(f'    text={texti!r}')
        if item is None:
            continue
        is_checked = (item.checkState() == Qt.Checked)
        if is_checked:
            ifile_groups.append(irow)
    print(f'ifile_groups = {ifile_groups}')
    return ifile_groups
    return np.array(ifile_groups, dtype='int32')

def get_property_table(model: BDF) -> tuple[np.ndarray, list, dict[int, tuple[str, str]]]:
    properties = []
    pids = set([])
    pid_to_type_name: dict[int, tuple[str, str]] = {}
    log = model.log
    for card in model.property_cards:
        if card.n == 0:
            continue
        card_type = card.type
        for i, pid in enumerate(card.property_id):
            pids.add(pid)
            comment = card.comment.get(pid, '')
            if comment == '':
                if card_type == 'PBARL':
                    idim0, idim1 = card.idim[i]
                    beam_type =card.Type[i]
                    dim = card.dims[idim0:idim1]
                    area, i1, i2, i12 = _bar_areaL('PBARL', beam_type, dim, card)
                    if i1 is None:
                        i1 = 0.
                        i2 = 0.
                        i12 = 0.
                    comment = (f'Type={beam_type} '
                               f'area={engieering_format_str(area)} '
                               f'I1={engieering_format_str(i1)} '
                               f'I2={engieering_format_str(i2)} '
                               f'I12={engieering_format_str(i12)}')
                elif card_type =='PBEAML':
                    comment = f'Type={card.Type[i]}'
                elif card_type == 'PSHELL':
                    comment = f't={card.t[i]} mid={card.material_id[i, :]}'
                elif card.type == 'PCOMP':
                    ilayer0, ilayer1 = card.ilayer[i]
                    thickness = sum(card.thickness[ilayer0:ilayer1])
                    theta_str = ply_format(card.theta[ilayer0:ilayer1])
                    comment = f't={thickness:g} theta={theta_str}'
                else:
                    log.warning(f'card.type={card.type} is not supported for comments')
            pid_to_type_name[pid] = (card_type, comment)

    pids_array = np.array(list(pids), dtype='int32')
    pids_array.sort()
    pid_to_type_name = {key: value for key, value in sorted(pid_to_type_name.items())}
    for pid, (prop_type, name) in pid_to_type_name.items():
        properties.append((f'{pid}: {prop_type}: {name}', True, pid, []))
    return pids_array, properties, pid_to_type_name

def get_material_table(model: BDF) -> tuple[np.ndarray, list, dict[int, tuple[str, str]]]:
    materials = []
    mids = set([])
    mid_to_type_name: dict[int, tuple[str, str]] = {}
    log = model.log
    for card in model.material_cards:
        if card.n == 0:
            continue
        card_type = card.type
        for i, mid in enumerate(card.material_id):
            mids.add(mid)
            comment = card.comment.get(mid, '')
            if card_type == 'MAT1':
                comment = (f'E={engieering_format_str(card.E[i])} '
                           f'G={engieering_format_str(card.G[i])} '
                           f'nu={card.nu[i]}')
            elif card.type == 'MAT8':
                comment = (f'E11={engieering_format_str(card.E11[i])} '
                           f'E22={engieering_format_str(card.E22[i])} '
                           f'nu12={card.nu12[i]}')
            else:
                log.warning(f'card.type={card.type} is not supported for comments')
            mid_to_type_name[mid] = (card_type, comment)

    mids_array = np.array(list(mids), dtype='int32')
    mids_array.sort()
    mid_to_type_name = {key: value for key, value in sorted(mid_to_type_name.items())}
    for mid, (prop_type, name) in mid_to_type_name.items():
        materials.append((f'{mid}: {card_type}: {name}', True, mid, []))
    return mids_array, materials, mid_to_type_name


def ply_format(theta: np.ndarray) -> str:
    thetai = theta.astype('int32')
    if np.allclose(theta, thetai):
        theta = thetai

    ntheta = len(theta)
    if ntheta > 4:
        is_even = (ntheta % 2 == 0)
        ntheta1 = ntheta // 2
        flag = ''
        if is_even:
            theta1 = theta[:ntheta1]
            theta2 = theta[ntheta1:]
            if np.allclose(theta1, theta2[::-1]):
                theta = theta1
                flag = ' Sym'
    theta_str = '/'.join(str(ti) for ti in theta)
    return theta_str + flag

def engieering_format_str(value: float) -> str:
    if value > 1.0:
        if value < 1000.0:
            return f'{value:g}'
        else:
            base, exp_str = f'{value:e}'.split('e')
            exp = int(exp_str)
            exp_remainder3 = exp % 3
            exp3 = exp - exp_remainder3
            base2 = float(base) * exp_remainder3
            return f'{base2}e+{exp3}'
    else:
        if value == 0:
            return '0'
        base, exp_str = f'{value:e}'.split('e')
        exp = int(exp_str)
        exp_remainder3 = -exp % 3
        exp3 = exp - exp_remainder3
        base2 = float(base) * exp_remainder3
        return f'{base2}e-{exp3}'

def _qitem_text_to_sline(text: str) -> tuple[str, int, str, str]:
    sline = text.split(':', 2)
    #print(f'sline = {sline}')
    idi_str, type, comment = sline
    type = type.strip()
    idi = int(idi_str)
    return text, idi, type, comment
