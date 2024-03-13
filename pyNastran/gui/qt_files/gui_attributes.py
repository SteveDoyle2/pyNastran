"""
defines GuiAttributes, which defines Gui getter/setter methods
and is inherited from many GUI classes
"""
from __future__ import annotations
import os
import sys
import traceback
from typing import Callable, Optional, Any, TYPE_CHECKING

import numpy as np
from vtkmodules.vtkRenderingCore import (
    vtkColorTransferFunction, vtkDataSetMapper, vtkTextActor, vtkRenderer)
from vtkmodules.vtkRenderingLOD import vtkLODActor

from qtpy import QtGui
from qtpy.QtWidgets import QMainWindow, QAction

try:
    import matplotlib
    IS_MATPLOTLIB = True
except ModuleNotFoundError:
    IS_MATPLOTLIB = False

import pyNastran
from pyNastran import DEV
from pyNastran.gui.typing import ColorFloat, Format
from pyNastran.gui.vtk_rendering_core import vtkPolyDataMapper
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.gui_objects.settings import Settings, FONT_SIZE_MIN, FONT_SIZE_MAX, force_ranged

from pyNastran.gui.qt_files.vtk_actor_actions import VtkActorActions
from pyNastran.gui.qt_files.tool_actions import ToolActions
from pyNastran.gui.qt_files.view_actions import ViewActions
from pyNastran.gui.qt_files.group_actions import GroupActions
from pyNastran.gui.qt_files.mouse_actions import MouseActions
from pyNastran.gui.qt_files.load_actions import LoadActions
from pyNastran.gui.qt_files.mark_actions import MarkActions
from pyNastran.gui.utils.qt.pydialog import make_font

from pyNastran.gui.menus.legend.legend_object import LegendObject
from pyNastran.gui.menus.highlight.highlight_object import HighlightObject, MarkObject
from pyNastran.gui.menus.preferences.preferences_object import PreferencesObject
IS_SHEAR_MOMENT_TORQUE = False
IS_CUTTING_PLANE = False
if IS_MATPLOTLIB:
    from pyNastran.gui.menus.cutting_plane.shear_moment_torque_object import ShearMomentTorqueObject
    IS_SHEAR_MOMENT_TORQUE = True
if IS_MATPLOTLIB and DEV:
    from pyNastran.gui.menus.cutting_plane.cutting_plane_object import CuttingPlaneObject
    IS_CUTTING_PLANE = True

from pyNastran.gui.menus.clipping.clipping_object import ClippingObject
from pyNastran.gui.menus.camera.camera_object import CameraObject
from pyNastran.gui.menus.edit_geometry_properties.edit_geometry_properties_object import (
    EditGeometryPropertiesObject)

from pyNastran.gui.qt_files.base_gui import BaseGui
from pyNastran.gui.utils.vtk.gui_utils import remove_actors_from_gui
from pyNastran.gui.utils.vtk.vtk_utils import (
    numpy_to_vtk_points, create_vtk_cells_of_constant_element_type)

from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.utils import print_bad_path
IS_TESTING = 'test' in sys.argv[0]
IS_OFFICIAL_RELEASE = 'dev' not in pyNastran.__version__
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow
    from pyNastran.gui.menus.results_sidebar import ResultsSidebar
    from pyNastran.gui.menus.groups_modify.groups_modify import Group
    from pyNastran.gui.qt_files.scalar_bar import ScalarBar
    #from vtkmodules.vtkFiltersGeneral import vtkAxes
    from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid, vtkPointData
    FollowerFunction = Callable[[dict[int, int], vtkUnstructuredGrid,
                                 vtkPointData, np.ndarray], None]


class GeometryObject(BaseGui):
    """
    """
    def __init__(self, gui: MainWindow):
        super().__init__(gui)
        #self.gui = parent

    def show(self):
        pass
    #def create(self):
        #"""
        #Create
         #- Point
         #- Line
         #- Surface
         #- Solid
         #- Coord
        #Modify
        #Delete
        #"""
        #pass
    #def create_point(self):
        #pass
    #def create_surface(self):
        #pass
    #def create_coord(self):
        #pass
    #def modify(self):
        #pass

class GuiAttributes:
    """All methods in this class must not require VTK"""
    def __init__(self, **kwds):
        """
        These variables are common between the GUI and
        the batch mode testing that fakes the GUI
        """
        inputs = kwds['inputs']
        res_widget = kwds['res_widget']
        self.dev = False
        self.log = None # it hasn't been initialized yet
        self.log_widget = None
        self._log_messages = []
        self._performance_mode = False
        #self.performance_mode = True

        # totally broken for solids
        self.make_contour_filter = False

        self.settings = Settings(self)
        self.vtk_actor_actions = VtkActorActions(self)
        self.tool_actions = ToolActions(self)
        self.view_actions = ViewActions(self)
        self.group_actions = GroupActions(self)
        self.mouse_actions = MouseActions(self)
        self.load_actions = LoadActions(self)
        self.mark_actions = MarkActions(self)

        self.legend_obj = LegendObject(self)
        self.camera_obj = CameraObject(self)
        self.clipping_obj = ClippingObject(self)
        self.highlight_obj = HighlightObject(self)
        self.mark_obj = MarkObject(self)
        self.preferences_obj = PreferencesObject(self)
        if IS_CUTTING_PLANE:
            self.cutting_plane_obj = CuttingPlaneObject(self)
        if IS_SHEAR_MOMENT_TORQUE:
            self.shear_moment_torque_obj = ShearMomentTorqueObject(self)

        self.edit_geometry_properties_obj = EditGeometryPropertiesObject(self)
        self.geometry_obj = GeometryObject(self)

        self.min_max_actors = []

        self.glyph_scale_factor = 1.0
        self.html_logging = False

        # the result type being currently shown
        # for a Nastran NodeID/displacement, this is 'node'
        # for a Nastran ElementID/PropertyID, this is 'element'
        self.result_location = None

        self.obj_names = []
        self.case_keys = []
        self.res_widget: ResultsSidebar = res_widget
        self._show_flag = True
        self.observers = {}

        # the gui is actually running
        # we set this to False when testing
        self.is_gui = True

        # testing enables additional checks
        # it's different than if we're just running tests
        if 'test' in inputs:
            self.is_testing_flag = inputs['test']
        else:
            self.is_testing_flag = False

        # just initializing the variable
        self.is_groups = False
        self._logo = None
        self._script_path = None
        self._icon_path = ''

        self.title = None
        self.min_value = None
        self.max_value = None
        self.blue_to_red = False
        self._is_axes_shown = True
        self.nvalues = 9
        #-------------

        # window variables
        self._modify_groups_window_shown = False
        #self._label_window = None
        #-------------
        # inputs dict
        #self.format = ''
        debug = inputs['debug']
        self.debug = debug
        assert debug in [True, False], 'debug=%s' % debug

        #-------------
        # format
        self.format = None
        self.format_class_map = {}
        self.supported_formats = []
        self.fmt_order = []

        self.infile_name = None
        self.out_filename = None

        # file
        self.menu_bar_format = None
        self.dirname = ''
        self.last_dir = '' # last visited directory while opening file
        self._default_python_file = None

        #-------------
        # internal params
        self.ncases = 0
        self.icase = 0
        self.icase_disp = None
        self.icase_vector = None
        self.icase_fringe = None

        self.nnodes = 0
        self.nelements = 0

        self.model_type = None

        self.tools = []
        self.checkables = []
        self.actions: dict[str, QAction] = {}
        self.modules = {}

        # actor_slots
        self.corner_text_actors: dict[int, vtkTextActor] = {}
        self.geometry_actors = {}
        self.alt_grids = {} #additional grids

        # coords
        self.transform = {}
        self.axes = {}

        #geom = Geom(color, line_thickness, etc.)
        #self.geometry_properties = {
        #    'name' : Geom(),
        #}

        self.model_data = ModelData(self)
        self.num_user_points = 0

        self._is_displaced = False
        self._is_forces = False
        self._is_fringe = False

        self._xyz_nominal = None

        self.nvalues = 9
        self.nid_maps: dict[str, dict[int, int]] = {}
        self.eid_maps: dict[str, dict[int, int]] = {}
        self.name = 'main'
        self.models = {}

        #if not isinstance(res_widget, MockResWidget):
            #if qt_version == 5:
                #super(QMainWindow, self).__init__()

        self.main_grids = {}
        self.main_grid_mappers = {}
        self.main_geometry_actors = {}

        self.main_edge_mappers = {}
        self.main_edge_actors = {}

        self.color_order = [
            (1.0, 0.145098039216, 1.0),
            (0.0823529411765, 0.0823529411765, 1.0),
            (0.0901960784314, 1.0, 0.941176470588),
            (0.501960784314, 1.0, 0.0941176470588),
            (1.0, 1.0, 0.117647058824),
            (1.0, 0.662745098039, 0.113725490196)
        ]

        self.color_function_black = vtkColorTransferFunction()
        self.color_function_black.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
        self.color_function_black.AddRGBPoint(1.0, 0.0, 0.0, 0.0)

    @property
    def geometry_properties(self):
        return self.model_data.geometry_properties
    @geometry_properties.setter
    def geometry_properties(self, geometry_properties):
        self.model_data.geometry_properties = geometry_properties

    @property
    def groups(self) -> dict[str, Any]:
        return self.model_data.groups
    @groups.setter
    def groups(self, groups: dict[str, Any]):
        self.model_data.groups = groups

    @property
    def group_active(self):
        return self.model_data.group_active
    @group_active.setter
    def group_active(self, group_active: str):
        self.model_data.group_active = group_active

    @property
    def follower_nodes(self) -> dict[str, list[int]]:
        return self.model_data.follower_nodes
    @follower_nodes.setter
    def follower_nodes(self, follower_nodes: dict[str, list[int]]) -> None:
        self.model_data.follower_nodes = follower_nodes

    @property
    def follower_functions(self) -> dict[str, FollowerFunction]:
        return self.model_data.follower_functions
    @follower_functions.setter
    def follower_functions(self, follower_functions: dict[str, FollowerFunction]):
        self.model_data.follower_functions = follower_functions

    @property
    def label_actors(self):
        return self.model_data.label_actors
    @label_actors.setter
    def label_actors(self, label_actors: list[Any]):
        self.model_data.label_actors = label_actors

    @property
    def label_ids(self):
        return self.model_data.label_ids
    @label_ids.setter
    def label_ids(self, label_ids: list[int]):
        self.model_data.label_ids = label_ids

    @property
    def label_scale(self) -> float:
        return self.model_data.label_scale
    @property
    def label_scale(self, label_scale: float) -> None:
        self.model_data.label_scale = label_scale

    @property
    def result_cases(self):
        return self.model_data.result_cases
    @result_cases.setter
    def result_cases(self, result_cases: dict[int, Any]):
        self.model_data.result_cases = result_cases

    @property
    def performance_mode(self) -> bool:
        """get the performance mode"""
        return self._performance_mode

    @performance_mode.setter
    def performance_mode(self, performance_mode: bool) -> None:
        """
        Set the performance mode.  If performance mode flips
        to False, we dump the log buffer.
        """
        if not performance_mode and self._log_messages:
            msg = ''.join(self._log_messages)
            #setUpdatesEnabled(False)
            #TxtBrows.append(SomeBigHTMLString)
            self._log_msg(msg)
            #setUpdatesEnabled(True)
            self._log_messages = []
        self._performance_mode = performance_mode

    def start_stop_performance_mode(func):
        """
        Suppresses logging.  If we started with logging suppressed,
        we won't unsuppress logging at the end of the function.
        """
        def new_func(self, *args, **kwargs):
            """The actual function exec'd by the decorated function."""
            performance_mode_initial = self.performance_mode
            if not performance_mode_initial:
                self.performance_mode = True
            try:
                n = func(self, *args, **kwargs)
            except Exception:
                if not performance_mode_initial:
                    self.performance_mode = False
                raise
            if not performance_mode_initial:
                self.performance_mode = False
            return n
        return new_func

    #-------------------------------------------------------------------
    # deprecated attributes
    def deprecated(self, old_name: str, new_name: str,
                   deprecated_version: Optional[list[int]]) -> None:
        """
        Throws a deprecation message and crashes if past a specific version.

        Parameters
        ----------
        old_name : str
            the old function name
        new_name : str
            the new function name
        deprecated_version : float
            the version the method was first deprecated in

        """
        deprecated(old_name, new_name, deprecated_version, levels=[0])

    #-------------------------------------------------------------------
    # geom
    def clear_actor(self, actor_name: str) -> None:
        if actor_name in self.gui.alt_grids:
            del self.alt_grids[actor_name]

        rend: vtkRenderer = self.rend
        if actor_name in self.geometry_actors:
            actor = self.geometry_actors[actor_name]
            rend.RemoveActor(actor)
            del self.geometry_actors[actor_name]

    @property
    def grid(self):
        """gets the active grid"""
        #print('get grid; %r' % self.name)
        return self.main_grids[self.name]

    @grid.setter
    def grid(self, grid):
        """sets the active grid"""
        #print('set grid; %r' % self.name)
        self.main_grids[self.name] = grid

    @property
    def grid_mapper(self):
        """gets the active grid_mapper"""
        return self.main_grid_mappers[self.name]

    @grid_mapper.setter
    def grid_mapper(self, grid_mapper):
        """sets the active grid_mapper"""
        self.main_grid_mappers[self.name] = grid_mapper

    @property
    def geom_actor(self):
        """gets the active geom_actor"""
        return self.main_geometry_actors[self.name]

    @geom_actor.setter
    def geom_actor(self, geom_actor):
        """sets the active geom_actor"""
        self.main_geometry_actors[self.name] = geom_actor

    #-------------------------------------------------------------------
    # edges
    @property
    def edge_mapper(self):
        return self.main_edge_mappers[self.name]

    @edge_mapper.setter
    def edge_mapper(self, edge_mapper):
        self.main_edge_mappers[self.name] = edge_mapper

    @property
    def edge_actor(self):
        """gets the active edge_actor"""
        return self.main_edge_actors[self.name]

    @edge_actor.setter
    def edge_actor(self, edge_actor):
        """sets the active edge_actor"""
        self.main_edge_actors[self.name] = edge_actor

    def set_glyph_scale_factor(self, scale):
        """sets the glyph scale factor"""
        if scale == np.nan:
            self.log.error('cannot set loads scale factor because no 1D, 2D, or 3D elements exist')
            return
        self.glyph_scale_factor = scale
        self.glyphs.SetScaleFactor(scale)

    @property
    def nid_map(self):
        """gets the node_id map"""
        return self.nid_maps[self.name]

    @nid_map.setter
    def nid_map(self, nid_map):
        """sets the node_id map"""
        self.nid_maps[self.name] = nid_map

    @property
    def eid_map(self):
        """gets the element_id map"""
        try:
            return self.eid_maps[self.name]
        except Exception:
            msg = 'KeyError: key=%r; keys=%s' % (self.name, list(self.eid_maps.keys()))
            raise KeyError(msg)

    @eid_map.setter
    def eid_map(self, eid_map):
        """sets the element_id map"""
        self.eid_maps[self.name] = eid_map

    #-------------------------------------------------------------------
    def set_point_grid(self, name: str,
                       nodes: np.ndarray, elements: np.ndarray,
                       color: ColorFloat,
                       point_size: int=5, opacity: int=1.,
                       add: bool=True) -> vtkUnstructuredGrid:
        """Makes a POINT grid"""
        self.create_alternate_vtk_grid(name, color=color, point_size=point_size,
                                       opacity=opacity, representation='point')

        nnodes = nodes.shape[0]
        if nnodes == 0:
            return

        assert isinstance(nodes, np.ndarray), type(nodes)

        points = numpy_to_vtk_points(nodes)
        grid = self.alt_grids[name]
        grid.SetPoints(points)

        etype = 9  # vtkQuad().GetCellType()
        create_vtk_cells_of_constant_element_type(grid, elements, etype)

        if add:
            self._add_alt_actors({name : self.alt_grids[name]})

            #if name in self.geometry_actors:
            self.geometry_actors[name].Modified()

    def set_quad_grid(self, name: str,
                      nodes: np.ndarray,
                      elements: np.ndarray,
                      color: ColorFloat,
                      point_size: int=5, line_width: int=5,
                      opacity: float=1.,
                      representation: str='wire',
                      add: bool=True,
                      visible_in_geometry_properties: bool=True) -> Optional[vtkUnstructuredGrid]:
        """Makes a quad grid"""
        etype = 9
        grid = self.vtk_actor_actions.create_grid_from_nodes_elements_etype(
            name, nodes, elements, etype, color,
            point_size=point_size, line_width=line_width,
            opacity=opacity,
            representation=representation,
            add=add,
            visible_in_geometry_properties=visible_in_geometry_properties,
        )
        return grid

    def _add_alt_actors(self, grids_dict: dict[str, vtkUnstructuredGrid],
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

    def _remove_alt_actors(self, names=None) -> None:
        if names is None:
            names = list(self.geometry_actors.keys())
            names.remove('main')
        for name in names:
            actor = self.geometry_actors[name]
            self.rend.RemoveActor(actor)
            del actor

    def _get_geometry_property_items(self, name: str,
                                     *property_name_defaults) -> list[Any]:
        """
        Matlab-esque way of accessing properties

        line_width = gui.get_geometry_property_items(
            LINE_NAME,
            'line_width', 5)
        line_width, opacity = gui.get_geometry_property_items(
            LINE_NAME,
            'line_width', 5,
            'opacity', 1.0)
        """
        length = len(property_name_defaults)
        assert length % 2 == 0, property_name_defaults
        out = []
        for i in range(0, length, 2):
            prop_name = property_name_defaults[i]
            assert prop_name in ['line_width', 'point_size', 'color', 'opacity'], prop_name
            if name in self.geometry_properties:
                prop = self.geometry_properties[name]
                outi = getattr(prop, prop_name)
            else:
                outi = property_name_defaults[i+1]
            out.append(outi)
        return out

    @property
    def displacement_scale_factor(self) -> float:
        """
        # dim_max = max_val * scale
        # scale = dim_max / max_val
        # 0.25 added just cause

        scale = self.displacement_scale_factor / tnorm_abs_max
        """
        #scale = dim_max / tnorm_abs_max * 0.25
        scale = self.settings.dim_max * 0.25
        return scale

    def set_script_path(self, script_path: str) -> None:
        """Sets the path to the custom script directory"""
        self._script_path = script_path

    def set_icon_path(self, icon_path: str) -> None:
        """Sets the path to the icon directory where custom icons are found"""
        self._icon_path = icon_path

    def form(self):
        formi = self.res_widget.get_form()
        return formi

    def get_form(self):
        return self._form

    def set_form(self, formi):
        self._form = formi
        data = []
        for key in self.case_keys:
            assert isinstance(key, int), key
            unused_obj, (i, resname) = self.result_cases[key]
            assert resname != 'main', resname
            form_tuple = (i, [])
            data.append(form_tuple)
        self.res_widget.update_results(formi, self.name)

        key = list(self.case_keys)[0]
        location = self.get_case_location(key)
        method = 'centroid' if location else 'nodal'

        data2 = [(method, None, [])]
        self.res_widget.update_methods(data2)

    def _remove_old_geometry(self, geom_filename: str) -> bool:
        skip_reading = False
        if self.dev:
            return skip_reading

        self.eid_map = {}
        self.nid_map = {}
        params_to_delete = (
            'case_keys', 'icase', 'isubcase_name_map',
            'result_cases', 'eid_map', 'nid_map',
        )
        if geom_filename is None or geom_filename == '':
            skip_reading = True
            return skip_reading
        else:
            self.turn_corner_text_off()
            self.grid.Reset()

            self.model_data.result_cases = {}
            self.ncases = 0
            for param in params_to_delete:
                if hasattr(self, param):  # TODO: is this correct???
                    try:
                        delattr(self, param)
                    except AttributeError:
                        msg = 'cannot delete %r; hasattr=%r' % (param, hasattr(self, param))
                        self.log.warning(msg)

            skip_reading = False
        #self.scalar_bar_actor.VisibilityOff()
        self.scalar_bar_actor.Modified()
        return skip_reading

    #---------------------------------------------------------------------------
    @start_stop_performance_mode
    def on_run_script(self, python_file=False) -> bool:
        """pulldown for running a python script"""
        return self.load_actions.on_run_script(python_file)

    def _execute_python_code(self, txt: str, show_msg: bool=True) -> bool:
        """executes python code"""
        is_passed = False
        if len(txt) == 0:
            return is_passed
        if show_msg:
            self.log_command(txt)
        try:
            exec(txt)
        except TypeError as error:
            self.log_error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)
            self.log_error(str(error))
            self.log_error(str(txt))
            self.log_error(str(type(txt)))
            return is_passed
        except Exception as error:
            #self.log_error(traceback.print_stack(f))
            self.log_error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)
            self.log_error(str(error))
            self.log_error(str(txt))
            return is_passed
        is_passed = True
        return is_passed

    #---------------------------------------------------------------------------
    def reset_labels(self, reset_minus1: bool=True) -> None:
        """
        Wipe all labels and regenerate the key slots based on the case keys.
        This is used when changing the model.
        """
        self._remove_labels()

        reset_minus1 = True
        # new geometry
        if reset_minus1:
            self.model_data.label_actors = {-1 : []}
        else:
            for idi in self.label_actors:
                if idi == -1:
                    continue
                self.label_actors[idi] = []
        self.label_ids = {}

        #self.case_keys = [
            #(1, 'ElementID', 1, 'centroid', '%.0f'),
            #(1, 'Region', 1, 'centroid', '%.0f')
        #]
        for icase in self.case_keys:
            #result_name = self.get_result_name(icase)
            self.label_actors[icase] = []
            self.label_ids[icase] = set()
        #print(self.label_actors)
        #print(self.label_ids)

    def _remove_labels(self) -> None:
        """
        Remove all labels from the current result case.
        This happens when the user explicitly selects the clear label button.
        """
        if len(self.label_actors) == 0:
            self.log.warning('No actors to remove')
            return

        # existing geometry
        for icase, actors in self.label_actors.items():
            if icase == -1:
                continue
            for actor in actors:
                self.rend.RemoveActor(actor)
                del actor
            self.label_actors[icase] = []
            self.label_ids[icase] = set()

    def clear_labels(self) -> None:
        """This clears out all labels from all result cases."""
        if len(self.label_actors) == 0:
            self.log.warning('No actors to clear')
            return

        # existing geometry
        icase = self.icase
        if icase not in self.label_actors:
            self.log.warning('No actors to clear')
            return

        actors = self.label_actors[icase]
        remove_actors_from_gui(self, actors, render=True)
        self.label_actors[icase] = []
        self.label_ids[icase] = set()

    def resize_labels(self, case_keys=None, show_msg: bool=True) -> None:
        """
        This resizes labels for all result cases.
        TODO: not done...
        """
        if case_keys is None:
            names = 'None)  # None -> all'
            case_keys = sorted(self.label_actors.keys())
        else:
            mid = '%s,' * len(case_keys)
            names = '[' + mid[:-1] + '])'

        count = 0
        for icase in case_keys:
            actors = self.label_actors[icase]
            for actor in actors:
                actor.VisibilityOff()
                count += 1
        if count and show_msg:
            self.log_command('self.resize_labels(%s)' % names)

    #---------------------------------------------------------------------------
    def on_update_clipping(self, min_clip=None, max_clip=None) -> None:
        self.clipping_obj.on_update_clipping(min_clip=min_clip, max_clip=max_clip)

        #---------------------------------------------------------------------------
    def hide_legend(self) -> None:
        """hides the legend"""
        self.scalar_bar.VisibilityOff()
        #self.scalar_bar.is_shown = False
        self.legend_obj.hide_legend()

    def show_legend(self) -> None:
        """shows the legend"""
        self.scalar_bar.VisibilityOn()
        #self.scalar_bar.is_shown = True
        self.legend_obj.show_legend()

    def clear_legend(self) -> None:
        """clears the legend"""
        self._is_fringe = False
        self.legend_obj.clear_legend()

    def _set_legend_fringe(self, is_fringe) -> None:
        self._is_fringe = is_fringe
        self.legend_obj._set_legend_fringe(is_fringe)

    def on_update_legend(self,
                         title: str='Title',
                         min_value: float=0., max_value: float=1.,
                         scale: float=0.0, phase: float=0.0,
                         arrow_scale: float=1.,
                         data_format: str='%.0f',
                         is_low_to_high: bool=True,
                         is_discrete: bool=True,
                         is_horizontal: bool=True,
                         nlabels: Optional[int]=None,
                         labelsize: Optional[int]=None,
                         ncolors: Optional[int]=None,
                         colormap: Optional[str]=None,
                         is_shown: bool=True, render: bool=True) -> None:
        """
        Updates the legend/model

        Parameters
        ----------
        scale : float
            displacemnt scale factor; true scale

        TODO: speed up by using existing values to skip update steps

        """
        self.legend_obj.on_update_legend(
            title=title, min_value=min_value, max_value=max_value,
            scale=scale, phase=phase,
            arrow_scale=arrow_scale,
            data_format=data_format,
            is_low_to_high=is_low_to_high, is_discrete=is_discrete, is_horizontal=is_horizontal,
            nlabels=nlabels, labelsize=labelsize, ncolors=ncolors, colormap=colormap,
            is_shown=is_shown, render=render)

    def update_scalar_bar(self, title: str, min_value: float, max_value: float,
                          data_format: str,
                          nlabels: Optional[int]=None,
                          labelsize: Optional[int]=None,
                          ncolors: Optional[int]=None,
                          colormap: Optional[str]=None,
                          is_horizontal: Optional[bool]=None,
                          is_shown: bool=True) -> None:
        """
        Updates the Scalar Bar

        Parameters
        ----------
        title : str
            the scalar bar title
        min_value : float
            the blue value
        max_value :
            the red value
        data_format : str
            '%g','%f','%i', etc.
        nlabels : int (default=None -> auto)
            the number of labels
        labelsize : int (default=None -> auto)
            the label size
        ncolors : int (default=None -> auto)
            the number of colors
        colormap : varies
            str : the name
            ndarray : (N, 3) float ndarry
                red-green-blue array
        is_shown : bool
            show the scalar bar

        """
        if colormap is None:
            colormap = self.settings.colormap
        #print("update_scalar_bar min=%s max=%s" % (min_value, max_value))
        scalar_bar: ScalarBar = self.scalar_bar

        if is_horizontal is None:
            is_horizontal = self.settings.is_horizontal_scalar_bar
        assert isinstance(is_horizontal, bool), is_horizontal

        scalar_bar.update(title, min_value, max_value, data_format,
                          nlabels=nlabels, labelsize=labelsize,
                          ncolors=ncolors, colormap=colormap,
                          is_low_to_high=self.legend_obj.is_low_to_high,
                          is_horizontal=is_horizontal,
                          is_shown=is_shown)

    #---------------------------------------------------------------------------
    def create_global_axes(self, dim_max: float) -> None:
        """creates the global axis"""
        cid = 0
        self.tool_actions.create_coordinate_system(
            cid, dim_max, label='', origin=None, matrix_3x3=None, coord_type='xyz')

    #def create_corner_axis(self) -> None:
        #"""creates the axes that sits in the corner"""
        #self.tool_actions.create_corner_axis()

    def get_corner_axis_visiblity(self) -> bool:
        """gets the visibility of the corner axis"""
        corner_axis = self.corner_axis
        axes_actor = corner_axis.GetOrientationMarker()
        is_visible = axes_actor.GetVisibility()
        return is_visible

    def set_corner_axis_visiblity(self, is_visible, render: bool=True) -> None:
        """sets the visibility of the corner axis"""
        corner_axis = self.corner_axis
        axes_actor = corner_axis.GetOrientationMarker()
        axes_actor.SetVisibility(is_visible)
        if render:
            self.Render()

    def update_axes_length(self, dim_max: float) -> None:
        """
        sets the driving dimension for:
          - picking?
          - coordinate systems
          - label size
        """
        self.settings.dim_max = dim_max
        dim = self.settings.dim_max * self.settings.coord_scale
        self.on_set_axes_length(dim)

    def on_set_axes_length(self, dim=None) -> None:
        """
        scale coordinate system based on model length
        """
        if dim is None:
            dim = self.settings.dim_max * self.settings.coord_scale
        for axes in self.axes.values():
            axes.SetTotalLength(dim, dim, dim)

    #---------------------------------------------------------------------------
    @property
    def window_title(self) -> str:
        return self.getWindowTitle()

    @window_title.setter
    def window_title(self, msg) -> None:
        #msg2 = "%s - "  % self.base_window_title
        #msg2 += msg
        self.setWindowTitle(msg)

    def build_fmts(self, fmt_order: list[str], stop_on_failure: bool=False) -> None:
        """populates the formats that will be supported"""
        fmts: list[Format] = []
        self.supported_formats = []
        #assert 'h5nastran' in fmt_order
        for fmt in fmt_order:
            geom_results_funcs = 'get_%s_wildcard_geometry_results_functions' % fmt

            if fmt in self.format_class_map:
                cls = self.format_class_map[fmt](self)
                data = getattr(cls, geom_results_funcs)()
            elif hasattr(self, geom_results_funcs):
                data = getattr(self, geom_results_funcs)()
            else:
                msg = 'get_%s_wildcard_geometry_results_functions does not exist' % fmt
                if stop_on_failure:
                    raise RuntimeError(msg)
                if not IS_OFFICIAL_RELEASE:
                    if self.log is None:
                        print('***', msg)
                    else:
                        self.log_error(msg)
            _add_fmt(self.supported_formats, fmts, fmt, geom_results_funcs, data)

        if len(fmts) == 0:
            RuntimeError('No formats...expected=%s' % fmt_order)
        self.fmts: list[Format] = fmts
        #print("fmts =", fmts)

        if not IS_TESTING:  # pragma: no cover
            print('supported_formats = %s' % self.supported_formats)
        #assert 'h5nastran' in self.supported_formats, self.supported_formats
        if len(fmts) == 0:
            print('supported_formats = %s' % self.supported_formats)
            raise RuntimeError('no modules were loaded...')

    @property
    def model(self):
        return self.models[self.name]
    @model.setter
    def model(self, model) -> None:
        self.models[self.name] = model

    def _reset_model(self, name: str) -> None:
        """resets the grids; sets up alt_grids"""
        if hasattr(self, 'main_grids') and name not in self.main_grids:
            grid = vtkUnstructuredGrid()
            grid_mapper = vtkDataSetMapper()
            grid_mapper.SetInputData(grid)

            geom_actor = vtkLODActor()
            geom_actor.DragableOff()
            geom_actor.SetMapper(grid_mapper)
            self.rend.AddActor(geom_actor)

            self.name = 'main'
            self.models = {}
            self.grid = grid
            self.grid_mapper = grid_mapper
            self.geom_actor = geom_actor
            self.grid.Modified()

            # link the current "main" to the scalar bar
            scalar_range = self.grid_selected.GetScalarRange()
            self.grid_mapper.ScalarVisibilityOn()
            self.grid_mapper.SetScalarRange(scalar_range)
            self.grid_mapper.SetLookupTable(self.color_function)

            self.edge_actor = vtkLODActor()
            self.edge_actor.DragableOff()
            self.edge_mapper = vtkPolyDataMapper()

            # create the edges
            self.get_edges()
        else:
            self.grid.Reset()
            self.grid.Modified()

        # reset alt grids
        alt_names = self.alt_grids.keys()
        for alt_name in alt_names:
            self.alt_grids[alt_name].Reset()
            self.alt_grids[alt_name].Modified()

    #---------------------------------------------------------------------------
    @start_stop_performance_mode
    def on_load_geometry(self, infile_name=None, geometry_format=None,
                         name: str='main',
                         plot: bool=True,
                         stop_on_failure: bool=False) -> None:
        """
        Loads a baseline geometry

        Parameters
        ----------
        infile_name : str; default=None -> popup
            path to the filename
        geometry_format : str; default=None
            the geometry format for programmatic loading
        name : str; default='main'
            the name of the actor; don't use this
        plot : bool; default=True
            Should the baseline geometry have results created and plotted/rendered?
            If you're calling the on_load_results method immediately after, set it to False
        stop_on_failure : bool; default=True
            stop the code if there's an error

        """
        if name is False:
            # fixing weird pyqt5, python 3.12 issue
            name = 'main'
        self.load_actions.on_load_geometry(
            infile_name=infile_name, geometry_format=geometry_format,
            name=name, plot=plot, stop_on_failure=stop_on_failure)

    @start_stop_performance_mode
    def on_load_results(self, out_filename=None) -> None:
        """
        Loads a results file.  Must have called on_load_geometry first.

        Parameters
        ----------
        out_filename : str / None
            the path to the results file
        """
        self.load_actions.on_load_results(out_filename=out_filename)

    @start_stop_performance_mode
    def on_load_custom_results(self,
                               out_filename=None, restype=None,
                               stop_on_failure: bool=False) -> None:
        """will be a more generalized results reader"""
        self.load_actions.on_load_custom_results(
            out_filename=out_filename, restype=restype, stop_on_failure=stop_on_failure)

    @start_stop_performance_mode
    def load_patran_nod(self, nod_filename) -> None:
        """reads a Patran formatted *.nod file"""
        self.load_actions.load_patran_nod(nod_filename)

    @start_stop_performance_mode
    def load_batch_inputs(self, inputs) -> None:
        print('load_batch_inputs', inputs)
        geom_script = inputs['geomscript']
        if geom_script is not None:
            self.on_run_script(geom_script)

        formats = inputs['format']
        if isinstance(formats, str):
            formats = [formats]
        if not formats:
            return
        input_filenames = inputs['input']
        results_filename = inputs['output']
        plot = True
        if results_filename:
            plot = False

        #print('input_filename =', input_filename)
        #print(formats)
        #print(input_filenames)
        assert len(formats) == len(input_filenames)
        if input_filenames is not None:
            for form, input_filename in zip(formats, input_filenames):
                form = form.lower()
                if not os.path.exists(input_filename):
                    msg = 'input filename: %s does not exist\n%s' % (
                        input_filename, print_bad_path(input_filename))
                    self.log.error(msg)
                    if self.html_logging:
                        print(msg)
                    return
            for results_filenamei in results_filename:
                #print('results_filenamei =', results_filenamei)
                if results_filenamei is not None:
                    if not os.path.exists(results_filenamei):
                        msg = '%s does not exist\n%s' % (
                            results_filenamei, print_bad_path(results_filenamei))
                        self.log.error(msg)
                        if self.html_logging:
                            print(msg)
                        return

        #unused_is_geom_results = input_filename == results_filename and len(input_filenames) == 1
        unused_is_geom_results = False
        is_failed = False
        print('input_filenames =', input_filenames)
        for i, input_filename in enumerate(input_filenames):
            if i == 0:
                name = 'main'
            else:
                name = input_filename
            self.name = name
            #form = inputs['format'].lower()
            #if is_geom_results:
            #    is_failed = self.on_load_geometry_and_results(
            #        infile_name=input_filename, name=name, geometry_format=form,
            #        plot=plot, stop_on_failure=True)
            #else:
            is_failed = self.on_load_geometry(
                infile_name=input_filename, name=name, geometry_format=form,
                plot=plot, stop_on_failure=True)
        self.name = 'main'
        #print('keys =', self.nid_maps.keys())

        if is_failed:
            return
        if results_filename:  #  and not is_geom_results
            self.on_load_results(results_filename)

        post_script = inputs['postscript']
        if post_script is not None:
            self.on_run_script(post_script)
        self.on_reset_camera()
        #self.log.debug('debug')
        #self.log.info('info')
        #self.log.warning('warning')
        #self.log.error('error')

        #self.log_debug('debug2')
        #self.log_info('info2')
        #self.log_warning('warning2')
        #self.log_command('command2')
        #self.log_error('error2')
        self.vtk_interactor.Modified()

    #---------------------------------------------------------------------------
    def on_increase_font_size(self) -> None:
        """used by the hidden_tools for Ctrl +"""
        self.on_set_font_size(self.settings.font_size + 1)

    def on_decrease_font_size(self) -> None:
        """used by the hidden_tools for Ctrl -"""
        self.on_set_font_size(self.settings.font_size - 1)

    def on_set_font_size(self, font_size: int, show_command: bool=True) -> bool:
        """changes the font size"""
        is_failed = True
        if not isinstance(font_size, int):
            self.log_error('font_size=%r must be an integer; type=%s' % (
                font_size, type(font_size)))
            return is_failed
        font_size = force_ranged(font_size, min_value=FONT_SIZE_MIN, max_value=FONT_SIZE_MAX)
        if self.settings.font_size == font_size:
            return False
        self.settings.font_size = font_size
        font = make_font(self.settings.font_size, is_bold=False)
        self.setFont(font)

        if isinstance(self, QMainWindow):
            #self.toolbar.setFont(font)
            self.menu_file.setFont(font)
            self.menu_view.setFont(font)
            self.menu_window.setFont(font)
            self.menu_help.setFont(font)

        self.legend_obj.set_font_size(font_size)
        self.camera_obj.set_font_size(font_size)
        self.highlight_obj.set_font_size(font_size)
        self.mark_obj.set_font_size(font_size)
        self.clipping_obj.set_font_size(font_size)
        if self._modify_groups_window_shown:
            self._modify_groups_window.set_font_size(font_size)
        self.preferences_obj.set_font_size(font_size)
        if hasattr(self, 'cutting_plane_obj'):
            self.cutting_plane_obj.set_font_size(font_size)
        if hasattr(self, 'shear_moment_torque_obj'):
            self.shear_moment_torque_obj.set_font_size(font_size)
        self.edit_geometry_properties_obj.set_font_size(font_size)

        #self.menu_scripts.setFont(font)
        self.log_command(f'self.settings.on_set_font_size({font_size})')
        return False

    def make_cutting_plane(self, data) -> None:
        model_name = data['model_name']
        unused_model = self.models[model_name]
        cid_p1, p1 = data['p1']
        cid_p2, p2 = data['p2']
        unused_method, cid_zaxis, zaxis = data['zaxis']
        unused_xyz1 = self.model.coords[cid_p1].transform_node_to_global(p1)
        unused_xyz2 = self.model.coords[cid_p2].transform_node_to_global(p2)
        unused_zaxis_xyz = self.model.coords[cid_zaxis].transform_node_to_global(zaxis)

    #---------------------------------------------------------------------------
    def get_result_by_xyz_cell_id(self, node_xyz: np.ndarray,
                                  cell_id: int,
                                  icase: Optional[int]=None) -> tuple[str, Any, int, np.ndarray]:
        """won't handle multiple cell_ids/node_xyz"""
        out = self.mark_actions.get_result_by_xyz_cell_id(node_xyz, cell_id, icase=icase)
        if out is None:
            print('attrs.get_result_by_xyz_cell_id bug')
            return None
        result_name, result_values, node_id, xyz = out
        return result_name, result_values, node_id, xyz

    def get_result_by_cell_id(self, cell_id: int,
                              world_position: np.ndarray,
                              icase: Optional[int]=None) -> tuple[str, Any, np.ndarray]:
        """should handle multiple cell_ids"""
        res_name, result_values, xyz = self.mark_actions.get_result_by_cell_id(
            cell_id, world_position, icase=icase)
        return res_name, result_values, xyz

    def mark_elements(self, eids: list[int],
                      stop_on_failure: bool=False,
                      show_command: bool=True) -> None:
        """mark the elements by the ElementID"""
        icase_result = 1 # ElementID
        icase_to_apply = self.icase
        self.mark_elements_by_different_case(
            eids, icase_result, icase_to_apply,
            stop_on_failure=stop_on_failure, show_command=False)
        self.log_command(f'self.mark_elements(eids={eids})')

    def mark_elements_by_case(self, eids: list[int],
                              stop_on_failure: bool=False,
                              show_command: bool=True):
        """mark the elements by the current case"""
        icase_result = self.icase
        icase_to_apply = self.icase
        self.mark_elements_by_different_case(
            eids, icase_result, icase_to_apply,
            stop_on_failure=stop_on_failure, show_command=False)
        self.log_command(f'self.mark_elements_by_case(eids={eids})')

    def mark_elements_by_different_case(self, eids: list[int],
                                        icase_result: int,
                                        icase_to_apply: int,
                                        stop_on_failure: bool=False,
                                        show_command: bool=False) -> None:
        """
        Marks a series of elements with custom text labels

        Parameters
        ----------
        eids : int, list[int]
            the elements to apply a message to
        icase_result : int
            the case to draw the result from
        icase_to_apply : int
            the key in label_actors to slot the result into

        TODO: fix the following
        correct   : applies to the icase_to_apply
        incorrect : applies to the icase_result

        Examples
        --------
        .. code-block::

          eids = [16563, 16564, 8916703, 16499, 16500, 8916699,
                  16565, 16566, 8916706, 16502, 16503, 8916701]
          icase_result = 22
          icase_to_apply = 25
          self.mark_elements_by_different_case(eids, icase_result, icase_to_apply)

        """
        self.mark_actions.mark_elements_by_different_case(
            eids, icase_result, icase_to_apply,
            stop_on_failure=stop_on_failure, show_command=show_command)

    #def mark_max_elements(self, neids, show_command: bool=True):
        #"""mark the elements by the top/btm x elements"""
    #def mark_min_elements(self, neids, show_command: bool=True):
        #"""mark the elements by the top/btm x elements"""

    def mark_nodes(self, nids, icase, text) -> None:
        """
        Marks a series of nodes with custom text labels

        Parameters
        ----------
        nids : int, list[int]
            the nodes to apply a message to
        icase : int
            the key in label_actors to slot the result into
        text : str, list[str]
            the text to display

        0 corresponds to the NodeID result
        self.mark_nodes(1, 0, 'max')
        self.mark_nodes(6, 0, 'min')
        self.mark_nodes([1, 6], 0, 'max')
        self.mark_nodes([1, 6], 0, ['max', 'min'])

        """
        self.mark_actions.mark_nodes(nids, icase, text)

    def create_annotation(self, text: str,
                          x: float, y: float, z: float) -> None:
        """
        Creates the actual annotation

        Parameters
        ----------
        text : str
            the text to display
                the annotation object
        x, y, z : float
            the position of the label

        Returns
        -------
        annotation : vtkBillboardTextActor3D
            the annotation object

        """
        annotation = self.mark_actions.create_annotation(text, x, y, z)
        return annotation

    #---------------------------------------------------------------------------
    def on_update_geometry_properties_window(self, geometry_properties) -> None:
        """updates the EditGeometryProperties window"""
        self.edit_geometry_properties_obj.on_update_geometry_properties_window(
            geometry_properties)

    @start_stop_performance_mode
    def on_update_geometry_properties(self, out_data, name=None,
                                      write_log: bool=True) -> None:
        """
        Applies the changed properties to the different actors if
        something changed.

        Note that some of the values are limited.  This prevents
        points/lines from being shrunk to 0 and also the actor being
        actually "hidden" at the same time.  This prevents confusion
        when you try to show the actor and it's not visible.
        """
        self.edit_geometry_properties_obj.on_update_geometry_properties(
            out_data, name=name, write_log=write_log)

    @start_stop_performance_mode
    def on_update_geometry_properties_override_dialog(self, geometry_properties) -> None:
        """
        Update the goemetry properties and overwrite the options in the
        edit geometry properties dialog if it is open.

        Parameters
        -----------
        geometry_properties : dict {str : CoordProperties or AltGeometry}
            Dictionary from name to properties object. Only the names included in
            ``geometry_properties`` are modified.
        """
        self.edit_geometry_properties_obj.on_update_geometry_properties_override_dialog(
            geometry_properties)

    #---------------------------------------------------------------------------

    def turn_corner_text_off(self) -> None:
        """turns all the text actors off"""
        self.tool_actions.turn_corner_text_off()

    def turn_corner_text_on(self) -> None:
        """turns all the text actors on"""
        self.tool_actions.turn_corner_text_on()

    @start_stop_performance_mode
    def export_case_data(self, icases=None) -> None:
        """exports CSVs of the requested cases"""
        self.tool_actions.export_case_data(icases=icases)

    @start_stop_performance_mode
    def on_load_user_geom(self, csv_filename=None, name=None, color=None) -> None:
        """
        Loads a User Geometry CSV File of the form:

        #    id  x    y    z
        GRID, 1, 0.2, 0.3, 0.3
        GRID, 2, 1.2, 0.3, 0.3
        GRID, 3, 2.2, 0.3, 0.3
        GRID, 4, 5.2, 0.3, 0.3
        grid, 5, 5.2, 1.3, 2.3  # case insensitive

        #    ID, nodes
        BAR,  1, 1, 2
        TRI,  2, 1, 2, 3
        # this is a comment

        QUAD, 3, 1, 5, 3, 4
        QUAD, 4, 1, 2, 3, 4  # this is after a blank line

        #RESULT,4,CENTROID,AREA(%f),PROPERTY_ID(%i)
        # in element id sorted order: value1, value2
        #1.0, 2.0 # bar
        #1.0, 2.0 # tri
        #1.0, 2.0 # quad
        #1.0, 2.0 # quad

        #RESULT,NODE,NODEX(%f),NODEY(%f),NODEZ(%f)
        # same difference

        #RESULT,VECTOR3,GEOM,DXYZ
        # 3xN

        Parameters
        ----------
        csv_filename : str (default=None -> load a dialog)
            the path to the user geometry CSV file
        name : str (default=None -> extract from fname)
            the name for the user points
        color : (float, float, float)
            RGB values as 0.0 <= rgb <= 1.0

        """
        self.load_actions.on_load_user_geom(csv_filename=csv_filename, name=name, color=color)

    @start_stop_performance_mode
    def on_save_vtk(self, vtk_filename=None) -> bool:
        is_failed = self.tool_actions.on_save_vtk(
            vtk_filename=vtk_filename)
        return is_failed

    @start_stop_performance_mode
    def on_load_csv_points(self, csv_filename=None, name=None, color=None) -> bool:
        """
        Loads a User Points CSV File of the form:

        1.0, 2.0, 3.0
        1.5, 2.5, 3.5

        Parameters
        -----------
        csv_filename : str (default=None -> load a dialog)
            the path to the user points CSV file
        name : str (default=None -> extract from fname)
            the name for the user points
        color : (float, float, float)
            RGB values as 0.0 <= rgb <= 1.0

        """
        is_failed = self.load_actions.on_load_csv_points(
            csv_filename=csv_filename, name=name, color=color)
        return is_failed

    #---------------------------------------------------------------------------
    def create_groups_by_visible_result(self, nlimit: int=50):
        """
        Creates group by the active result

        This should really only be called for integer results < 50-ish.
        """
        return self.group_actions.create_groups_by_visible_result(nlimit=nlimit)

    def create_groups_by_property_id(self) -> int:
        """
        Creates a group for each Property ID.

        As this is somewhat Nastran specific, create_groups_by_visible_result exists as well.
        """
        return self.group_actions.create_groups_by_property_id()

    def create_groups_by_model_group(self) -> int:
        if hasattr(self, 'model') and hasattr(self.model, 'model_groups'):
            return self.group_actions.create_groups_by_model_group(self.model.model_groups)
        return 0

    #---------------------------------------------------------------------------
    def update_camera(self, code) -> None:
        self.view_actions.update_camera(code)

    def _update_camera(self, camera=None) -> None:
        self.view_actions._update_camera(camera)

    def on_pan_left(self, event) -> None:
        self.view_actions.on_pan_left(event)

    def on_pan_right(self, event) -> None:
        self.view_actions.on_pan_right(event)

    def on_pan_up(self, event) -> None:
        self.view_actions.on_pan_up(event)

    def on_pan_down(self, event) -> None:
        self.view_actions.on_pan_down(event)

    #------------------------------
    def rotate(self, rotate_deg, render: bool=True):
        """rotates the camera by a specified amount"""
        self.view_actions.rotate(rotate_deg, render=render)

    def on_rotate_clockwise(self):
        """rotate clockwise"""
        self.view_actions.rotate(15.0)

    def on_rotate_cclockwise(self):
        """rotate counter clockwise"""
        self.view_actions.rotate(-15.0)

    #------------------------------
    def zoom(self, value):
        return self.view_actions.zoom(value)

    def on_increase_magnification(self):
        """zoom in"""
        self.view_actions.on_increase_magnification()

    def on_decrease_magnification(self):
        """zoom out"""
        self.view_actions.on_decrease_magnification()

    def set_focal_point(self, focal_point):
        """
        Parameters
        ----------
        focal_point : (3, ) float ndarray
            The focal point
            [ 188.25109863 -7. -32.07858658]

        """
        self.view_actions.set_focal_point(focal_point)

    def on_surface(self):
        """sets the main/toggle actors to surface"""
        self.view_actions.on_surface()

    def on_wireframe(self):
        """sets the main/toggle actors to wirefreme"""
        self.view_actions.on_wireframe()

    def on_take_screenshot(self, fname=None, magnify=None, show_msg=True):
        """
        Take a screenshot of a current view and save as a file

        Parameters
        ----------
        fname : str; default=None
            None : pop open a window
            str : bypass the popup window
        magnify : int; default=None
            None : use self.settings.magnify
            int : resolution increase factor
        show_msg : bool; default=True
            log the command

        """
        self.tool_actions.on_take_screenshot(fname=fname, magnify=magnify, show_msg=show_msg)

    def get_camera_data(self):
        """see ``set_camera_data`` for arguments"""
        return self.camera_obj.get_camera_data()

    def on_set_camera(self, name: str, show_log=True):
        """see ``set_camera_data`` for arguments"""
        self.camera_obj.on_set_camera(name, show_log=show_log)

    def on_set_camera_data(self, camera_data, show_log: bool=True) -> None:
        """
        Sets the current camera

        Parameters
        ----------
        camera_data : dict[key] : value
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
        self.camera_obj.on_set_camera_data(camera_data, show_log=show_log)

    @property
    def IS_GUI_TESTING(self) -> bool:
        return 'test_' in sys.argv[0]
    @property
    def iren(self):
        return self.vtk_interactor
    @property
    def render_window(self):
        return self.vtk_interactor.GetRenderWindow()

    #------------------------------
    def get_xyz_cid0(self, model_name=None) -> np.ndarray:
        xyz = self.xyz_cid0
        return xyz

    def get_element_ids(self, model_name: Optional[str]=None,
                        ids: Optional[np.ndarray]=None) -> np.ndarray:
        """wrapper around element_ids"""
        #if self.group_active == 'main':
        eids_all = self.element_ids
        if ids is None:
            eids = eids_all
        else:
            eids = eids_all[ids]
        # skin_eids=1:1632
        if self.group_active == 'main':
            return eids

        group: Group = self.groups[self.group_active]
        group_eids = group.element_ids
        eids2 = np.intersect1d(eids, group_eids)
        return eids2

    def get_node_ids(self, model_name: Optional[str]=None,
                     ids: Optional[np.ndarray]=None) -> np.ndarray:
        """wrapper around node_ids"""
        nids_all = self.node_ids
        if ids is None:
            nids = nids_all
        else:
            nids = nids_all[ids]

        if self.group_active == 'main':
            return nids

        group: Group = self.groups[self.group_active]
        group_nids = group.node_ids
        nids2 = np.intersect1d(nids, group_nids)
        return nids2

    def get_reverse_node_ids(self, model_name=None, ids=None):
        """wrapper around node_ids"""
        if ids is None:
            return np.array([])
        # include all but the indicies specified
        include_ids = np.ones(self.node_ids.shape, dtype=self.node_ids.dtype)
        include_ids[ids] = 0
        return self.node_ids[include_ids]

    #------------------------------
    # these are overwritten
    def log_debug(self, msg: str) -> None:
        """turns logs into prints to aide debugging"""
        if self.debug:
            print('DEBUG:  ', msg)

    def log_info(self, msg: str) -> None:
        """turns logs into prints to aide debugging"""
        if self.debug:
            print('INFO:  ', msg)

    def log_error(self, msg: str) -> None:
        """turns logs into prints to aide debugging"""
        #if self.debug:
        print('ERROR:  ', msg)

    def log_warning(self, msg: str) -> None:
        """turns logs into prints to aide debugging"""
        if self.debug:
            print('WARNING:  ', msg)

    def log_command(self, msg: str) -> None:
        """turns logs into prints to aide debugging"""
        if self.debug:
            print('COMMAND:  ', msg)

    def Render(self) -> None:  # pragma: no cover
        pass


class ModelData:
    def __init__(self, parent: GuiAttributes):
        self.geometry_properties: dict[str, Any] = {}

        self.groups: dict[str, Any] = {}
        self.group_active = 'main'

        self.follower_nodes: dict[str, list[int]] = {}
        self.follower_functions: dict[str, Callable] = {}

        self.label_actors: dict[int, list[int]] = {-1 : []}
        self.label_ids = {}
        self.label_scale = 1.0 # in percent

        self.result_cases = {}

    def __repr__(self):
        msg = ('ModelData:\n'
               f'result_cases.keys() = {self.result_cases.keys()}')
        return msg

def _add_fmt(supported_fmts: list[str],
             fmts: list[Format], fmt: str,
             geom_results_funcs: str, data: Callable) -> None:
    """
    Adds a format

    Parameters
    ----------
    supported_fmts : list[str]
        the names in fmts (without duplicates) in the same order
    fmts : list[formats]
        format : list[fmt, macro_name, geo_fmt, geo_func, res_fmt, res_func]
        macro_name : ???
            ???
        geo_fmt : ???
            ???
        geo_func : ???
            ???
        res_fmt : ???
            ???
        res_func : ???
            ???
    fmt : str
        nastran, cart3d, etc.
    geom_results_funcs : str
        'get_nastran_wildcard_geometry_results_functions'
        'get_cart3d_wildcard_geometry_results_functions'
    data : function
        the outputs from ``get_nastran_wildcard_geometry_results_functions()``
        so 1 or more formats (macro_name, geo_fmt, geo_func, res_fmt, res_func)

    """
    msg = 'macro_name, geo_fmt, geo_func, res_fmt, res_func = data\n'
    msg += 'data = %s'
    if isinstance(data, tuple):
        assert len(data) == 5, msg % str(data)
        macro_name, geo_fmt, geo_func, res_fmt, res_func = data
        fmts.append((fmt, macro_name, geo_fmt, geo_func, res_fmt, res_func))
        supported_fmts.append(fmt)
    elif isinstance(data, list):
        for datai in data:
            assert len(datai) == 5, msg % str(datai)
            macro_name, geo_fmt, geo_func, res_fmt, res_func = datai
            fmts.append((fmt, macro_name, geo_fmt, geo_func, res_fmt, res_func))
            if fmt not in supported_fmts:
                supported_fmts.append(fmt)
    else:  # pragma: no cover
        raise TypeError(data)
