"""
defines GuiAttributes, which defines Gui getter/setter methods
and is inherited from many GUI classes
"""
import os
import sys
import traceback
from collections import OrderedDict
from typing import List

import numpy as np
import vtk
from qtpy import QtGui
from qtpy.QtWidgets import QMainWindow

try:
    import matplotlib
    IS_MATPLOTLIB = True
except ImportError:
    IS_MATPLOTLIB = False

import pyNastran
from pyNastran.gui.gui_objects.settings import Settings

from pyNastran.gui.qt_files.tool_actions import ToolActions
from pyNastran.gui.qt_files.view_actions import ViewActions
from pyNastran.gui.qt_files.group_actions import GroupActions
from pyNastran.gui.qt_files.mouse_actions import MouseActions
from pyNastran.gui.qt_files.load_actions import LoadActions
from pyNastran.gui.qt_files.mark_actions import MarkActions

from pyNastran.gui.menus.legend.legend_object import LegendObject
from pyNastran.gui.menus.highlight.highlight_object import HighlightObject, MarkObject
from pyNastran.gui.menus.preferences.preferences_object import PreferencesObject
IS_CUTTING_PLANE = False
if IS_MATPLOTLIB:
    from pyNastran.gui.menus.cutting_plane.cutting_plane_object import CuttingPlaneObject
    from pyNastran.gui.menus.cutting_plane.shear_moment_torque_object import ShearMomentTorqueObject
    IS_CUTTING_PLANE = True
from pyNastran.gui.menus.clipping.clipping_object import ClippingObject
from pyNastran.gui.menus.camera.camera_object import CameraObject
from pyNastran.gui.menus.edit_geometry_properties.edit_geometry_properties_object import (
    EditGeometryPropertiesObject)

from pyNastran.gui.utils.vtk.gui_utils import remove_actors_from_gui
from pyNastran.gui.utils.vtk.vtk_utils import (
    numpy_to_vtk_points, create_vtk_cells_of_constant_element_type)

from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.utils import print_bad_path
IS_TESTING = 'test' in sys.argv[0]
IS_OFFICIAL_RELEASE = 'dev' not in pyNastran.__version__

class GeometryObject:
    """
    """
    def __init__(self, parent):
        self.gui = parent

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
    #def crate_surface(self):
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
        if IS_MATPLOTLIB:
            self.cutting_plane_obj = CuttingPlaneObject(self)
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
        self.res_widget = res_widget
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
        self.is_edges = False
        self.is_edges_black = self.is_edges

        #self.format = ''
        debug = inputs['debug']
        self.debug = debug
        assert debug in [True, False], 'debug=%s' % debug

        #-------------
        # format
        self.format = None
        self.format_class_map = {}
        self.supported_formats = []
        self.fmts = []

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
        self.actions = {}
        self.modules = OrderedDict()

        # actor_slots
        self.text_actors = {}
        self.geometry_actors = OrderedDict()
        self.alt_grids = {} #additional grids

        # coords
        self.transform = {}
        self.axes = {}

        #geom = Geom(color, line_thickness, etc.)
        #self.geometry_properties = {
        #    'name' : Geom(),
        #}
        self.geometry_properties = OrderedDict()
        self.follower_nodes = {}
        self.follower_functions = {}

        self.label_actors = {-1 : []}
        self.label_ids = {}
        self.label_scale = 1.0 # in percent

        self.result_cases = {}
        self.num_user_points = 0

        self._is_displaced = False
        self._is_forces = False
        self._is_fringe = False

        self._xyz_nominal = None

        self.nvalues = 9
        self.nid_maps = {}
        self.eid_maps = {}
        self.name = 'main'
        self.models = {}

        self.groups = {}
        self.group_active = 'main'

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

        self.color_function_black = vtk.vtkColorTransferFunction()
        self.color_function_black.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
        self.color_function_black.AddRGBPoint(1.0, 0.0, 0.0, 0.0)

    @property
    def performance_mode(self):
        """get the performance mode"""
        return self._performance_mode

    @performance_mode.setter
    def performance_mode(self, performance_mode):
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
        Supresses logging.  If we started with logging suppressed,
        we won't unsupress logging at the end of the function.
        """
        def new_func(self, *args, **kwargs):
            """The actual function exec'd by the decorated function."""
            performance_mode_initial = self.performance_mode
            if not performance_mode_initial:
                self.performance_mode = True
            try:
                n = func(self, *args, **kwargs)
            except:
                if not performance_mode_initial:
                    self.performance_mode = False
                raise
            if not performance_mode_initial:
                self.performance_mode = False
            return n
        return new_func

    #-------------------------------------------------------------------
    # deprecated attributes
    def deprecated(self, old_name, new_name, deprecated_version):
        # type: (str, str, str, Optional[List[int]]) -> None
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
    def clear_actor(self, actor_name):
        if actor_name in self.gui.alt_grids:
            del self.alt_grids[actor_name]
        if actor_name in self.geometry_actors:
            actor = self.geometry_actors[actor_name]
            self.rend.RemoveActor(actor)
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
        except:
            msg = 'KeyError: key=%r; keys=%s' % (self.name, list(self.eid_maps.keys()))
            raise KeyError(msg)

    @eid_map.setter
    def eid_map(self, eid_map):
        """sets the element_id map"""
        self.eid_maps[self.name] = eid_map

    #-------------------------------------------------------------------
    def set_point_grid(self, name, nodes, elements, color,
                       point_size=5, opacity=1., add=True):
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

        etype = 9  # vtk.vtkQuad().GetCellType()
        create_vtk_cells_of_constant_element_type(grid, elements, etype)

        if add:
            self._add_alt_actors({name : self.alt_grids[name]})

            #if name in self.geometry_actors:
            self.geometry_actors[name].Modified()

    def set_quad_grid(self, name, nodes, elements, color,
                      line_width=5, opacity=1., representation='wire', add=True):
        """Makes a CQUAD4 grid"""
        self.create_alternate_vtk_grid(name, color=color, line_width=line_width,
                                       opacity=opacity, representation=representation)

        nnodes = nodes.shape[0]
        nquads = elements.shape[0]
        if nnodes == 0:
            return
        if nquads == 0:
            return

        #print('adding quad_grid %s; nnodes=%s nquads=%s' % (name, nnodes, nquads))
        assert isinstance(nodes, np.ndarray), type(nodes)

        points = numpy_to_vtk_points(nodes)
        grid = self.alt_grids[name]
        grid.SetPoints(points)

        etype = 9  # vtk.vtkQuad().GetCellType()
        create_vtk_cells_of_constant_element_type(grid, elements, etype)

        if add:
            self._add_alt_actors({name : self.alt_grids[name]})

            #if name in self.geometry_actors:
        self.geometry_actors[name].Modified()

    def _add_alt_actors(self, grids_dict, names_to_ignore=None):
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
            self.tool_actions._add_alt_geometry(grid, name)

    def _remove_alt_actors(self, names=None):
        if names is None:
            names = list(self.geometry_actors.keys())
            names.remove('main')
        for name in names:
            actor = self.geometry_actors[name]
            self.rend.RemoveActor(actor)
            del actor

    @property
    def displacement_scale_factor(self):
        """
        # dim_max = max_val * scale
        # scale = dim_max / max_val
        # 0.25 added just cause

        scale = self.displacement_scale_factor / tnorm_abs_max
        """
        #scale = dim_max / tnorm_abs_max * 0.25
        scale = self.settings.dim_max * 0.25
        return scale

    def set_script_path(self, script_path):
        """Sets the path to the custom script directory"""
        self._script_path = script_path

    def set_icon_path(self, icon_path):
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
            unused_obj, (i, unused_name) = self.result_cases[key]
            form_tuple = (i, [])
            data.append(form_tuple)

        self.res_widget.update_results(formi, self.name)

        key = list(self.case_keys)[0]
        location = self.get_case_location(key)
        method = 'centroid' if location else 'nodal'

        data2 = [(method, None, [])]
        self.res_widget.update_methods(data2)

    def _remove_old_geometry(self, geom_filename):
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
            self.turn_text_off()
            self.grid.Reset()

            self.result_cases = OrderedDict()
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
    def _create_load_file_dialog(self, qt_wildcard, title, default_filename=None):
        wildcard_level, fname = self.load_actions.create_load_file_dialog(
            qt_wildcard, title, default_filename=default_filename)
        return wildcard_level, fname

    @start_stop_performance_mode
    def on_run_script(self, python_file=False):
        """pulldown for running a python script"""
        is_passed = False
        if python_file in [None, False]:
            title = 'Choose a Python Script to Run'
            wildcard = "Python (*.py)"
            infile_name = self._create_load_file_dialog(
                wildcard, title, self._default_python_file)[1]
            if not infile_name:
                return is_passed # user clicked cancel

            #python_file = os.path.join(script_path, infile_name)
            python_file = os.path.join(infile_name)

        if not os.path.exists(python_file):
            msg = 'python_file = %r does not exist' % python_file
            self.log_error(msg)
            return is_passed

        with open(python_file, 'r') as python_file_obj:
            txt = python_file_obj.read()
        is_passed = self._execute_python_code(txt, show_msg=False)
        if not is_passed:
            return is_passed
        self._default_python_file = python_file
        self.log_command('self.on_run_script(%r)' % python_file)
        print('self.on_run_script(%r)' % python_file)
        return is_passed

    def _execute_python_code(self, txt, show_msg=True):
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
    def reset_labels(self, reset_minus1=True):
        """
        Wipe all labels and regenerate the key slots based on the case keys.
        This is used when changing the model.
        """
        self._remove_labels()

        reset_minus1 = True
        # new geometry
        if reset_minus1:
            self.label_actors = {-1 : []}
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

    def _remove_labels(self):
        """
        Remove all labels from the current result case.
        This happens when the user explictly selects the clear label button.
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

    def clear_labels(self):
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

    def resize_labels(self, case_keys=None, show_msg=True):
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
            self.log_command('resize_labels(%s)' % names)

    #---------------------------------------------------------------------------
    def on_update_clipping(self, min_clip=None, max_clip=None):
        self.clipping_obj.on_update_clipping(min_clip=min_clip, max_clip=max_clip)

        #---------------------------------------------------------------------------
    def hide_legend(self):
        """hides the legend"""
        self.scalar_bar.VisibilityOff()
        #self.scalar_bar.is_shown = False
        self.legend_obj.hide_legend()

    def show_legend(self):
        """shows the legend"""
        self.scalar_bar.VisibilityOn()
        #self.scalar_bar.is_shown = True
        self.legend_obj.show_legend()

    def clear_legend(self):
        """clears the legend"""
        self._is_fringe = False
        self.legend_obj.clear_legend()

    def _set_legend_fringe(self, is_fringe):
        self._is_fringe = is_fringe
        self.legend_obj._set_legend_fringe(is_fringe)

    def on_update_legend(self,
                         title='Title', min_value=0., max_value=1.,
                         scale=0.0, phase=0.0,
                         arrow_scale=1.,
                         data_format='%.0f',
                         is_low_to_high=True, is_discrete=True, is_horizontal=True,
                         nlabels=None, labelsize=None, ncolors=None, colormap=None,
                         is_shown=True, render=True):
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

    def update_scalar_bar(self, title, min_value, max_value,
                          data_format,
                          nlabels=None, labelsize=None,
                          ncolors=None, colormap=None,
                          is_shown=True):
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
        self.scalar_bar.update(title, min_value, max_value, data_format,
                               nlabels=nlabels, labelsize=labelsize,
                               ncolors=ncolors, colormap=colormap,
                               is_low_to_high=self.legend_obj.is_low_to_high,
                               is_horizontal=self.legend_obj.is_horizontal_scalar_bar,
                               is_shown=is_shown)

    def on_update_scalar_bar(self, title, min_value, max_value, data_format):
        self.title = str(title)
        self.min_value = float(min_value)
        self.max_value = float(max_value)

        try:
            data_format % 1
        except:
            msg = ("failed applying the data formatter format=%r and "
                   "should be of the form: '%i', '%8f', '%.2f', '%e', etc.")
            self.log_error(msg)
            return
        #self.data_format = data_format
        self.log_command('on_update_scalar_bar(%r, %r, %r, %r)' % (
            title, min_value, max_value, data_format))

    #---------------------------------------------------------------------------
    def create_coordinate_system(self, coord_id: int, dim_max: float, label: str='',
                                 origin=None, matrix_3x3=None,
                                 coord_type: str='xyz'):
        """
        Creates a coordinate system

        Parameters
        ----------
        coord_id : int
            the coordinate system id
        dim_max : float
            the max model dimension; 10% of the max will be used for the coord length
        label : str
            the coord id or other unique label (default is empty to indicate the global frame)
        origin : (3, ) ndarray/list/tuple
            the origin
        matrix_3x3 : (3, 3) ndarray
            a standard Nastran-style coordinate system
        coord_type : str
            a string of 'xyz', 'Rtz', 'Rtp' (xyz, cylindrical, spherical)
            that changes the axis names

        .. todo::  coord_type is not supported ('xyz' ONLY)
        .. todo::  Can only set one coordinate system

        """
        self.tool_actions.create_coordinate_system(
            coord_id, dim_max, label=label,
            origin=origin, matrix_3x3=matrix_3x3,
            coord_type=coord_type)

    def create_global_axes(self, dim_max: float):
        """creates the global axis"""
        cid = 0
        self.tool_actions.create_coordinate_system(
            cid, dim_max, label='', origin=None, matrix_3x3=None, coord_type='xyz')

    def create_corner_axis(self):
        """creates the axes that sits in the corner"""
        self.tool_actions.create_corner_axis()

    def get_corner_axis_visiblity(self):
        """gets the visibility of the corner axis"""
        corner_axis = self.corner_axis
        axes_actor = corner_axis.GetOrientationMarker()
        is_visible = axes_actor.GetVisibility()
        return is_visible

    def set_corner_axis_visiblity(self, is_visible, render=True):
        """sets the visibility of the corner axis"""
        corner_axis = self.corner_axis
        axes_actor = corner_axis.GetOrientationMarker()
        axes_actor.SetVisibility(is_visible)
        if render:
            self.Render()

    def update_axes_length(self, dim_max):
        """
        sets the driving dimension for:
          - picking?
          - coordinate systems
          - label size
        """
        self.settings.dim_max = dim_max
        dim = self.settings.dim_max * self.settings.coord_scale
        self.on_set_axes_length(dim)

    def on_set_axes_length(self, dim=None):
        """
        scale coordinate system based on model length
        """
        if dim is None:
            dim = self.settings.dim_max * self.settings.coord_scale
        for axes in self.axes.values():
            axes.SetTotalLength(dim, dim, dim)

    #---------------------------------------------------------------------------
    @property
    def window_title(self):
        return self.getWindowTitle()

    @window_title.setter
    def window_title(self, msg):
        #msg2 = "%s - "  % self.base_window_title
        #msg2 += msg
        self.setWindowTitle(msg)

    def build_fmts(self, fmt_order: List[str], stop_on_failure: bool=False):
        """populates the formats that will be supported"""
        fmts = []
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
            _add_fmt(fmts, fmt, geom_results_funcs, data)

        if len(fmts) == 0:
            RuntimeError('No formats...expected=%s' % fmt_order)
        self.fmts = fmts
        #print("fmts =", fmts)

        self.supported_formats = [fmt[0] for fmt in fmts]
        if not IS_TESTING:  # pragma: no cover
            print('supported_formats = %s' % self.supported_formats)
        #assert 'cart3d' in self.supported_formats, self.supported_formats
        if len(fmts) == 0:
            print('supported_formats = %s' % self.supported_formats)
            raise RuntimeError('no modules were loaded...')

    @property
    def model(self):
        return self.models[self.name]
    @model.setter
    def model(self, model):
        self.models[self.name] = model

    def _reset_model(self, name: str):
        """resets the grids; sets up alt_grids"""
        if hasattr(self, 'main_grids') and name not in self.main_grids:
            grid = vtk.vtkUnstructuredGrid()
            grid_mapper = vtk.vtkDataSetMapper()
            grid_mapper.SetInputData(grid)

            geom_actor = vtk.vtkLODActor()
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

            self.edge_actor = vtk.vtkLODActor()
            self.edge_actor.DragableOff()
            self.edge_mapper = vtk.vtkPolyDataMapper()

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
    def on_load_geometry(self, infile_name=None, geometry_format=None, name='main',
                         plot=True, raise_error=False):
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
        raise_error : bool; default=True
            stop the code if True
        """
        self.load_actions.on_load_geometry(
            infile_name=infile_name, geometry_format=geometry_format,
            name=name, plot=plot, raise_error=raise_error)

    @start_stop_performance_mode
    def on_load_results(self, out_filename=None):
        """
        Loads a results file.  Must have called on_load_geometry first.

        Parameters
        ----------
        out_filename : str / None
            the path to the results file
        """
        self.load_actions.on_load_results(out_filename=out_filename)

    @start_stop_performance_mode
    def on_load_custom_results(self, out_filename=None, restype=None, stop_on_failure: bool=False):
        """will be a more generalized results reader"""
        self.load_actions.on_load_custom_results(
            out_filename=out_filename, restype=restype, stop_on_failure=stop_on_failure)

    @start_stop_performance_mode
    def load_patran_nod(self, nod_filename):
        """reads a Patran formatted *.nod file"""
        self.load_actions.load_patran_nod(nod_filename)

    @start_stop_performance_mode
    def load_batch_inputs(self, inputs):
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
            #        plot=plot, raise_error=True)
            #else:
            is_failed = self.on_load_geometry(
                infile_name=input_filename, name=name, geometry_format=form,
                plot=plot, raise_error=True)
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
    def on_increase_font_size(self):
        """used by the hidden_tools for Ctrl +"""
        self.on_set_font_size(self.settings.font_size + 1)

    def on_decrease_font_size(self):
        """used by the hidden_tools for Ctrl -"""
        self.on_set_font_size(self.settings.font_size - 1)

    def on_set_font_size(self, font_size: int, show_command: bool=True):
        """changes the font size"""
        is_failed = True
        if not isinstance(font_size, int):
            self.log_error('font_size=%r must be an integer; type=%s' % (
                font_size, type(font_size)))
            return is_failed
        if font_size < 6:
            font_size = 6
        if self.settings.font_size == font_size:
            return False
        self.settings.font_size = font_size
        font = QtGui.QFont()
        font.setPointSize(self.settings.font_size)
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
        self.log_command('settings.on_set_font_size(%s)' % font_size)
        return False

    def make_cutting_plane(self, data):
        model_name = data['model_name']
        unused_model = self.models[model_name]
        cid_p1, p1 = data['p1']
        cid_p2, p2 = data['p2']
        unused_method, cid_zaxis, zaxis = data['zaxis']
        unused_xyz1 = self.model.coords[cid_p1].transform_node_to_global(p1)
        unused_xyz2 = self.model.coords[cid_p2].transform_node_to_global(p2)
        unused_zaxis_xyz = self.model.coords[cid_zaxis].transform_node_to_global(zaxis)

    #---------------------------------------------------------------------------
    def get_result_by_xyz_cell_id(self, node_xyz, cell_id):
        """won't handle multiple cell_ids/node_xyz"""
        result_name, result_values, node_id, xyz = self.mark_actions.get_result_by_xyz_cell_id(
            node_xyz, cell_id)
        return result_name, result_values, node_id, xyz

    def get_result_by_cell_id(self, cell_id, world_position, icase=None):
        """should handle multiple cell_ids"""
        res_name, result_values, xyz = self.mark_actions.get_result_by_cell_id(
            cell_id, world_position, icase=icase)
        return res_name, result_values, xyz

    def mark_elements(self, eids,
                      stop_on_failure: bool=False, show_command: bool=True):
        """mark the elements by the ElementID"""
        icase_result = 1 # ElementID
        icase_to_apply = self.icase
        self.mark_elements_by_different_case(
            eids, icase_result, icase_to_apply,
            stop_on_failure=stop_on_failure, show_command=False)
        self.log_command(f'mark_elements(eids={eids})')

    def mark_elements_by_case(self, eids,
                              stop_on_failure: bool=False, show_command: bool=True):
        """mark the elements by the current case"""
        icase_result = self.icase
        icase_to_apply = self.icase
        self.mark_elements_by_different_case(
            eids, icase_result, icase_to_apply,
            stop_on_failure=stop_on_failure, show_command=False)
        self.log_command(f'mark_elements_by_case(eids={eids})')

    def mark_elements_by_different_case(self, eids, icase_result: int, icase_to_apply: int,
                                        stop_on_failure: bool=False, show_command: bool=False):
        """
        Marks a series of elements with custom text labels

        Parameters
        ----------
        eids : int, List[int]
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

    def mark_nodes(self, nids, icase, text):
        """
        Marks a series of nodes with custom text labels

        Parameters
        ----------
        nids : int, List[int]
            the nodes to apply a message to
        icase : int
            the key in label_actors to slot the result into
        text : str, List[str]
            the text to display

        0 corresponds to the NodeID result
        self.mark_nodes(1, 0, 'max')
        self.mark_nodes(6, 0, 'min')
        self.mark_nodes([1, 6], 0, 'max')
        self.mark_nodes([1, 6], 0, ['max', 'min'])

        """
        self.mark_actions.mark_nodes(nids, icase, text)

    def create_annotation(self, text, x, y, z):
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
    def on_update_geometry_properties_window(self, geometry_properties):
        """updates the EditGeometryProperties window"""
        self.edit_geometry_properties_obj.on_update_geometry_properties_window(
            geometry_properties)

    @start_stop_performance_mode
    def on_update_geometry_properties(self, out_data, name=None, write_log=True):
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
    def on_update_geometry_properties_override_dialog(self, geometry_properties):
        """
        Update the goemetry properties and overwite the options in the
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
    def update_text_actors(self, subcase_id, subtitle,
                           imin, min_value,
                           imax, max_value, label, location):
        """
        Updates the text actors in the lower left

        Max:  1242.3
        Min:  0.
        Subcase: 1 Subtitle:
        Label: SUBCASE 1; Static
        """
        self.tool_actions.update_text_actors(subcase_id, subtitle,
                                             imin, min_value,
                                             imax, max_value, label, location)

    def create_text(self, position, label, text_size=18):
        """creates the lower left text actors"""
        self.tool_actions.create_text(position, label, text_size=text_size)

    def turn_text_off(self):
        """turns all the text actors off"""
        self.tool_actions.turn_text_off()

    def turn_text_on(self):
        """turns all the text actors on"""
        self.tool_actions.turn_text_on()

    @start_stop_performance_mode
    def export_case_data(self, icases=None):
        """exports CSVs of the requested cases"""
        self.tool_actions.export_case_data(icases=icases)

    @start_stop_performance_mode
    def on_load_user_geom(self, csv_filename=None, name=None, color=None):
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
        self.tool_actions.on_load_user_geom(csv_filename=csv_filename, name=name, color=color)

    @start_stop_performance_mode
    def on_load_csv_points(self, csv_filename=None, name=None, color=None):
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
        is_failed = self.tool_actions.on_load_csv_points(
            csv_filename=csv_filename, name=name, color=color)
        return is_failed

    #---------------------------------------------------------------------------
    def create_groups_by_visible_result(self, nlimit=50):
        """
        Creates group by the active result

        This should really only be called for integer results < 50-ish.
        """
        return self.group_actions.create_groups_by_visible_result(nlimit=nlimit)

    def create_groups_by_property_id(self):
        """
        Creates a group for each Property ID.

        As this is somewhat Nastran specific, create_groups_by_visible_result exists as well.
        """
        return self.group_actions.create_groups_by_property_id()

    #---------------------------------------------------------------------------
    def update_camera(self, code):
        self.view_actions.update_camera(code)

    def _update_camera(self, camera=None):
        self.view_actions._update_camera(camera)

    def on_pan_left(self, event):
        self.view_actions.on_pan_left(event)

    def on_pan_right(self, event):
        self.view_actions.on_pan_right(event)

    def on_pan_up(self, event):
        self.view_actions.on_pan_up(event)

    def on_pan_down(self, event):
        self.view_actions.on_pan_down(event)

    #------------------------------
    def rotate(self, rotate_deg, render=True):
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

    def on_set_camera(self, name, show_log=True):
        """see ``set_camera_data`` for arguments"""
        self.camera_obj.on_set_camera(name, show_log=show_log)

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
        self.camera_obj.on_set_camera_data(camera_data, show_log=show_log)

    @property
    def IS_GUI_TESTING(self):
        return 'test_' in sys.argv[0]
    @property
    def iren(self):
        return self.vtk_interactor
    @property
    def render_window(self):
        return self.vtk_interactor.GetRenderWindow()

    #------------------------------
    def get_xyz_cid0(self, model_name=None):
        xyz = self.xyz_cid0
        return xyz

    def get_element_ids(self, model_name=None, ids=None):
        """wrapper around element_ids"""
        if ids is None:
            return self.element_ids
        return self.element_ids[ids]

    def get_node_ids(self, model_name=None, ids=None):
        """wrapper around node_ids"""
        if ids is None:
            return self.node_ids
        return self.node_ids[ids]

    def get_reverse_node_ids(self, model_name=None, ids=None):
        """wrapper around node_ids"""
        if ids is None:
            return np.array([])
        # include all but the indicies sepecified
        include_ids = np.ones(self.node_ids.shape, dtype=self.node_ids.dtype)
        include_ids[ids] = 0
        return self.node_ids[include_ids]

    #------------------------------
    # these are overwritten
    def log_debug(self, msg):
        """turns logs into prints to aide debugging"""
        if self.debug:
            print('DEBUG:  ', msg)

    def log_info(self, msg):
        """turns logs into prints to aide debugging"""
        if self.debug:
            print('INFO:  ', msg)

    def log_error(self, msg):
        """turns logs into prints to aide debugging"""
        #if self.debug:
        print('ERROR:  ', msg)

    def log_warning(self, msg):
        """turns logs into prints to aide debugging"""
        if self.debug:
            print('WARNING:  ', msg)

    def log_command(self, msg):
        """turns logs into prints to aide debugging"""
        if self.debug:
            print('COMMAND:  ', msg)

    def Render(self):  # pragma: no cover
        pass

def _add_fmt(fmts, fmt, geom_results_funcs, data):
    """
    Adds a format

    Parameters
    ----------
    fmts : List[formats]
        format : List[fmt, macro_name, geo_fmt, geo_func, res_fmt, res_func]
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
    elif isinstance(data, list):
        for datai in data:
            assert len(datai) == 5, msg % str(datai)
            macro_name, geo_fmt, geo_func, res_fmt, res_func = datai
            fmts.append((fmt, macro_name, geo_fmt, geo_func, res_fmt, res_func))
    else:
        raise TypeError(data)
