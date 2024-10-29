"""
defines:
 - Settings(parent)
   - reset_settings(self)
   - load(self, settings)
   - save(self, settings)
   - on_increase_text_size(self)
   - on_decrease_font_size(self)
   - on_set_font_size(self, font_size, show_command=True)
   - set_annotation_size_color(self, size=None, color=None)
   - set_annotation_size(self, size, render=True)
   - set_annotation_color(self, color, render=True)
   - set_background_color_to_white(self)
   - set_background_color(self, color)
   - set_corner_text_color(self, color)
   - update_text_size(self, magnify=1.0)

 - repr_settings(settings)

"""
from __future__ import annotations
import os
import copy
from typing import Optional, Any, TYPE_CHECKING
import numpy as np
from qtpy import QtGui

from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.gui_objects.coord_properties import CoordProperties
from pyNastran.gui.gui_objects.utils import get_setting
from pyNastran.gui.utils.colormaps import colormap_keys as COLORMAPS
from pyNastran.utils import object_attributes #, object_stats

from pyNastran.gui.qt_files.colors import (
    BLACK_FLOAT, WHITE_FLOAT, GREY_FLOAT, ORANGE_FLOAT, HOT_PINK_FLOAT,
    YELLOW_FLOAT, LIGHT_GREEN_FLOAT,
)
if TYPE_CHECKING:  # pragma: no cover
    from vtkmodules.vtkFiltersGeneral import vtkAxes
    from qtpy.QtCore import QSettings
    from pyNastran.gui.typing import ColorFloat
    from pyNastran.gui.main_window import MainWindow

from pyNastran.gui import (
    USE_OLD_SIDEBAR_OBJS_ as USE_OLD_SIDEBAR_OBJECTS,
    USE_NEW_SIDEBAR_OBJS_ as USE_NEW_SIDEBAR_OBJECTS,
    USE_NEW_TERMS_ as USE_NEW_TERMS)

USE_NEW_SIDEBAR = False
if USE_NEW_SIDEBAR_OBJECTS:
    USE_NEW_SIDEBAR = False


IS_WINDOWS = 'nt' in os.name

BACKGROUND_COLOR = GREY_FLOAT
BACKGROUND_COLOR2 = GREY_FLOAT

HIGHLIGHT_COLOR = ORANGE_FLOAT
HIGHLIGHT_OPACITY = 0.9
HIGHLIGHT_POINT_SIZE = 12.0
HIGHLIGHT_LINE_WIDTH = 5.0

SHEAR_MOMENT_TORQUE_COLOR = HOT_PINK_FLOAT
SHEAR_MOMENT_TORQUE_OPACITY = 0.9
SHEAR_MOMENT_TORQUE_POINT_SIZE = 12.0
SHEAR_MOMENT_TORQUE_LINE_WIDTH = 5.0

LINE_WIDTH_MIN = 0.1
LINE_WIDTH_MAX = 2000.

IS_TRACKBALL_CAMERA = True
USE_PARALLEL_PROJECTION = True
DEFAULT_COLORMAP = 'jet'
NFILES_TO_SAVE = 9
OPACITY_MIN = 0.1
OPACITY_MAX = 1.0

FONT_SIZE = 8
FONT_SIZE_MIN = 7
FONT_SIZE_MAX = 20

ANNOTATION_SIZE = 18
ANNOTATION_SIZE_MIN = 1
ANNOTATION_SIZE_MAX = 500
ANNOTATION_COLOR = BLACK_FLOAT

POINT_SIZE_MIN = 5.0
POINT_SIZE_MAX = 30.0

COORD_TEXT_SCALE = 0.5 # percent of nominal
COORD_TEXT_SCALE_MIN = 0.1
COORD_TEXT_SCALE_MAX = 2000.

CORNER_TEXT_SIZE = 14
CORNER_TEXT_SIZE_MIN = 7
CORNER_TEXT_SIZE_MAX = 30
CORNER_TEXT_COLOR = BLACK_FLOAT

COORD_SCALE = 0.05  # in percent of max dimension
COORD_SCALE_MIN = 0.1
COORD_SCALE_MAX = 1000.

MAGNIFY = 5
MAGNIFY_MIN = 1
MAGNIFY_MAX = 10

CAERO_COLOR = YELLOW_FLOAT
RBE_LINE_COLOR = LIGHT_GREEN_FLOAT
DISPLACEMENT_MODEL_SCALE = 0.1

NASTRAN_VERSIONS = ['Guess', 'MSC', 'NX', 'Optistruct']
NASTRAN_COLOR_KEYS = ['nastran_caero_color', 'nastran_rbe_line_color']
NASTRAN_STR_KEYS = ['nastran_version']
NASTRAN_BOOL_KEYS = [
    'nastran_create_coords',
    'nastran_is_properties',
    'nastran_is_element_quality',
    'nastran_is_bar_axes',
    'nastran_is_3d_bars', 'nastran_is_3d_bars_update',
    'nastran_is_mass_update',
    'nastran_is_constraints',
    'nastran_is_shell_mcids',
    'nastran_is_rbe',

    'nastran_displacement', 'nastran_velocity',
    'nastran_acceleration', 'nastran_eigenvector', 'nastran_temperature',
    'nastran_spc_force', 'nastran_mpc_force',
    'nastran_applied_load', 'nastran_heat_flux',

    'nastran_stress', 'nastran_plate_stress', 'nastran_composite_plate_stress',
    'nastran_strain', 'nastran_plate_strain', 'nastran_composite_plate_strain',
    'nastran_rod_stress', 'nastran_bar_stress', 'nastran_beam_stress',
    'nastran_rod_strain', 'nastran_bar_strain', 'nastran_beam_strain',
    'nastran_spring_stress', 'nastran_solid_stress',
    'nastran_spring_strain', 'nastran_solid_strain',

    'nastran_force',
    'nastran_bar_force', 'nastran_beam_force', 'nastran_plate_force',
    'nastran_spring_force', 'nastran_gap_force', 'nastran_cbush_force',

    'nastran_grid_point_force', 'nastran_strain_energy',

    # --------------------------------------------------------------
    'nastran_show_caero_sub_panels',
    'nastran_show_caero_actor',
    #'nastran_show_control_surfaces',
    #'nastran_show_conm',
]
NASTRAN_BOOL_STR_KEYS = NASTRAN_BOOL_KEYS + NASTRAN_STR_KEYS

OTHER_STRING_KEYS = [
    'units_length', 'units_force', 'units_moment',
    'units_stress', 'units_pressure',
    'units_displacement',
    'units_velocity', 'units_acceleration',
]
OTHER_LIST_INT_KEYS = [
    'cart3d_fluent_include', 'cart3d_fluent_remove',]
OTHER_LIST_STR_KEYS = [
    'units_model_in', ]


class OtherSettings:
    def __init__(self, parent):
        """
        Creates the OtherSettings object
        """
        self.parent = parent
        self.reset_settings()

    def reset_settings(self) -> None:
        self.cart3d_fluent_include = ()
        self.cart3d_fluent_remove = ()
        #('in', 'lbf', 's', 'psi')
        self.units_model_in = ('unitless','','','')
        self.units_length = 'in'
        #self.units_area = 'in^2'
        self.units_force = 'lbf'
        self.units_moment = 'in-lbf'
        self.units_stress = 'psi'
        self.units_pressure = 'psi'
        self.units_displacement = 'in'
        self.units_velocity = 'in/s'
        self.units_acceleration = 'in/s^2'

    def update(self, out_data: dict[str, Any]):
        """from preferences"""
        self.cart3d_fluent_include = out_data['cart3d_fluent_include']
        self.cart3d_fluent_remove = out_data['cart3d_fluent_remove']
        #self.units_model_in = out_data['units_model_in']
        #self.units_length = out_data['units_length']
        #self.units_area = out_data['units_area']
        #self.units_force  = out_data['units_force']
        #self.units_moment = out_data['units_moment']
        #self.units_stress = out_data['units_stress']
        #self.units_pressure = out_data['units_pressure']
        #self.units_displacement = out_data['units_displacement']
        #self.units_velocity = out_data['units_velocity']
        #self.units_acceleration = out_data['units_acceleration']

    def save(self, settings: QSettings) -> None:
        """save_json -> all keys"""
        keys = object_attributes(self, mode='public', keys_to_skip=['parent'])
        for key in keys:
            #base, key2 = key.split('_', 1)
            value = getattr(self, key)
            #print(f'save other: key={key!r} value={value!r}')
            settings.setValue(key, value)

    def __repr__(self) -> str:
        msg = '<OtherSettings>\n'
        keys = object_attributes(self, mode='public', keys_to_skip=['parent'])
        for key in keys:
            #if key.startswith('nastran'):
            #    raise RuntimeError(key)
            value = getattr(self, key)
            if isinstance(value, tuple):
                value = str(value)
            msg += '  %r = %r\n' % (key, value)
        return msg

class NastranSettings:
    def __init__(self, parent):
        """
        Creates the NastranSettings object
        """
        self.parent = parent
        self.reset_settings()

    def reset_settings(self):
        self.version = 'Guess'
        self.is_element_quality = True
        self.is_properties = True
        self.is_3d_bars = True
        self.is_3d_bars_update = True
        self.is_mass_update = True
        self.is_constraints = True
        self.create_coords = True
        self.is_bar_axes = True
        self.is_shell_mcids = True
        self.is_rbe = True

        self.stress = True
        self.spring_stress = True
        self.rod_stress = True
        self.bar_stress = True
        self.beam_stress = True
        self.plate_stress = True
        self.composite_plate_stress = True
        self.solid_stress = True

        self.strain = True
        self.spring_strain = True
        self.rod_strain = True
        self.bar_strain = True
        self.beam_strain = True
        self.plate_strain = True
        self.composite_plate_strain = True
        self.solid_strain = True

        self.force = True
        self.spring_force = True
        self.cbush_force = True
        self.gap_force = True
        self.bar_force = True
        self.beam_force = True
        self.plate_force = True

        self.eigenvector = True
        self.displacement = True
        self.velocity = True
        self.acceleration = True
        self.temperature = True

        self.spc_force = True
        self.mpc_force = True
        self.applied_load = True
        self.heat_flux = True

        #self.stress = True
        #self.stress = True
        #self.strain = True
        self.strain_energy = True
        self.grid_point_force = True

        #---------------------------------------0000000000000000
        # colors
        self.caero_color = CAERO_COLOR
        self.rbe_line_color = RBE_LINE_COLOR

        # ------------------------------------------------------

        #: flips the nastran CAERO subpaneling
        #:   False -> borders of CAEROs can be seen
        #:   True  -> individual subpanels can be seen
        self.show_caero_sub_panels = False

        self.show_caero_actor = True  # show the caero mesh
        #self.show_control_surfaces = True
        #self.show_conm = True

    def save(self, settings: QSettings) -> None:
        #print(nastran_settings)
        for key in NASTRAN_BOOL_STR_KEYS:
            base, key2 = key.split('_', 1)
            value = getattr(self, key2)
            settings.setValue(key, value)
            #print(f'*key={key!r} key2={key2!r} value={value!r}')

        for key in NASTRAN_COLOR_KEYS:
            base, key2 = key.split('_', 1)
            value = getattr(self, key2)
            settings.setValue(key, value)

    def set_caero_color(self, color: ColorFloat, render: bool=True) -> None:
        """
        Set the CAEROx color

        Parameters
        ----------
        color : ColorFloat
            RGB values as floats
        """
        self.caero_color = color
        parent = self.parent
        if render:
            parent.vtk_interactor.Render()
        parent.log_command('self.settings.nastran_settings.set_caero_color(%s, %s, %s)' % color)

    def set_rbe_line_color(self, color: ColorFloat, render: bool=True) -> None:
        """
        Set the RBE2/RBE3 line color

        Parameters
        ----------
        color : ColorFloat
            RGB values as floats
        """
        self.rbe_line_color = color
        parent = self.parent
        if render:
            parent.vtk_interactor.Render()
        parent.log_command('self.settings.nastran_settings.set_rbe_line_color(%s, %s, %s)' % color)

    def __repr__(self) -> str:
        msg = '<NastranSettings>\n'
        keys = object_attributes(self, mode='public', keys_to_skip=['parent'])
        for key in keys:
            if key.startswith('nastran'):
                raise RuntimeError(key)
            value = getattr(self, key)
            if isinstance(value, tuple):
                value = str(value)
            msg += '  %r = %r\n' % (key, value)
        return msg


class Settings:
    """storage class for various settings"""
    def __init__(self, parent: MainWindow):
        """
        Creates the Settings object

        Parameters
        ----------
        parent : MainWindow()
            used by the class to access the MainWindow
        """
        self.parent = parent

        self.recent_files = []

        #self.annotation_scale = 1.0

        self.reset_settings(resize=True, reset_dim_max=True)
        self.nastran_settings = NastranSettings(parent)
        self.other_settings = OtherSettings(parent)

    def reset_settings(self, resize: bool=True,
                       reset_dim_max: bool=True) -> None:
        """helper method for ``setup_gui``"""
        # rgb tuple
        self.use_gradient_background = True
        self.background_color = BACKGROUND_COLOR
        self.background_color2 = BACKGROUND_COLOR2

        # grab bag
        # int
        self.font_size = FONT_SIZE
        self.magnify = MAGNIFY

        # default directory: os.getcwd()
        # activates after the first directory selection
        self.startup_directory = ''
        # False  os.getcwd()
        # True:  startup directory
        self.use_startup_directory = True

        self.use_old_sidebar_objects = USE_OLD_SIDEBAR_OBJECTS
        self.use_new_sidebar_objects = USE_NEW_SIDEBAR_OBJECTS
        self.use_new_sidebar = USE_NEW_SIDEBAR
        self.use_new_terms = USE_NEW_TERMS

        # probe color
        # this includes:
        #   - min/max actors
        #   - probes
        self.annotation_size = ANNOTATION_SIZE    # int
        self.annotation_color = ANNOTATION_COLOR  # rgb floats

        # text in the lower left corner
        self.corner_text_size = CORNER_TEXT_SIZE    # int
        self.corner_text_color = CORNER_TEXT_COLOR  # rgb floats

        self.highlight_color = HIGHLIGHT_COLOR            # rgb floats
        self.highlight_opacity = HIGHLIGHT_OPACITY        # float
        self.highlight_point_size = HIGHLIGHT_POINT_SIZE  # int
        self.highlight_line_width = HIGHLIGHT_LINE_WIDTH  # float

        self.shear_moment_torque_color = SHEAR_MOMENT_TORQUE_COLOR           # float
        self.shear_moment_torque_opacity = SHEAR_MOMENT_TORQUE_OPACITY       # rgb floats
        self.shear_moment_torque_point_size = SHEAR_MOMENT_TORQUE_POINT_SIZE # float
        self.shear_moment_torque_line_width = SHEAR_MOMENT_TORQUE_LINE_WIDTH # float

        self.is_trackball_camera = IS_TRACKBALL_CAMERA
        self.use_parallel_projection = USE_PARALLEL_PROJECTION
        self.displacement_model_scale = DISPLACEMENT_MODEL_SCALE
        self.show_info = True
        self.show_debug = False
        self.show_command = True
        self.show_warning = True
        self.show_error = True

        self.is_edges_visible = True
        self.is_edges_black = True

        self.is_horizontal_scalar_bar = False
        self.is_min_visible = True
        self.is_max_visible = True

        # float
        self.coord_scale = COORD_SCALE            # float
        self.coord_text_scale = COORD_TEXT_SCALE  # float
        self.coord_linewidth = 2.0

        # string
        self.colormap = 'jet' # 'viridis'

        if resize:
            self.parent.resize(1100, 700)

        if reset_dim_max: # not stored
            self.dim_max = 1.0
        #self.annotation_scale = 1.0

        self.nastran_settings = NastranSettings(self.parent)

    def finish_startup(self):
        self.set_background_color(self.background_color, render=False, quiet=True)
        self.set_background_color2(self.background_color2, render=False, quiet=True)
        self.set_gradient_background(self.use_gradient_background, render=True, quiet=True)

    def add_model_settings_to_dict(self, data: dict[str, Any]):
        nastran_settings = self.nastran_settings
        for key in NASTRAN_BOOL_STR_KEYS:
            base, key2 = key.split('_', 1)
            data[key] = getattr(nastran_settings, key2)
        for key in NASTRAN_COLOR_KEYS:
            base, key2 = key.split('_', 1)
            data[key] = getattr(nastran_settings, key2)

    def load(self, settings: QSettings) -> bool:
        """helper method for ``setup_gui``"""
        #red = (1.0, 0.0, 0.0)
        screen_shape_default = (1100, 700)

        setting_keys = [str(key) for key in settings.childKeys()]

        # sets the window size/position
        main_window_geometry = get_setting(
            settings, setting_keys, ['main_window_geometry'], None)
        main_window_state = get_setting(
            settings, setting_keys, ['main_window_state'], None)
        if main_window_geometry is not None:
            self.parent.restoreGeometry(main_window_geometry)
        if main_window_state is not None:
            self.parent.restoreState(main_window_state)
        #settings.setValue('main_window_geometry', self.parent.saveGeometry())
        #settings.setValue('main_window_state', self.parent.saveState())

        # this is the gui font
        self._set_setting(settings, setting_keys, ['font_size'],
                          default=self.font_size,
                          save=True, auto_type=int)
        self.font_size = force_ranged(
            self.font_size, min_value=FONT_SIZE_MIN, max_value=None)

        # parallel/perspective
        self._set_setting(settings, setting_keys, ['is_trackball_camera'],
                          default=self.is_trackball_camera,
                          save=True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['use_parallel_projection'],
                          default=self.use_parallel_projection,
                          save=True, auto_type=bool)

        # launch in local or specified directory
        self._set_setting(settings, setting_keys, ['startup_directory'],
                          default=self.startup_directory,
                          save=True, auto_type=str)
        self._set_setting(settings, setting_keys, ['use_startup_directory'],
                          default=self.use_startup_directory,
                          save=True, auto_type=bool)
        if os.path.exists(self.startup_directory):
            self.parent.last_dir = self.startup_directory
        else:
            self.startup_directory = ''

        self._set_setting(settings, setting_keys, ['use_old_sidebar_objects'],
                          default=self.use_old_sidebar_objects,
                          save=True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['use_new_sidebar_objects'],
                          default=self.use_new_sidebar_objects,
                          save=True, auto_type=bool)
        self.use_new_sidebar = self.use_new_sidebar_objects
        #self._set_setting(settings, setting_keys, ['use_new_sidebar'], self.use_new_sidebar,
                          #USE_NEW_SIDEBAR, auto_type=bool)
        self._set_setting(settings, setting_keys, ['use_new_terms'],
                          default=self.use_new_terms,
                          save=True, auto_type=bool)

        # the info/debug/gui/command preferences
        self._set_setting(settings, setting_keys, ['show_debug'],
                          default=self.show_debug,
                          save=True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_info'],
                          default=self.show_info,
                          save=True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_command'],
                          default=self.show_command,
                          save=True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_warning'],
                          default=self.show_warning,
                          save=True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_error'],
                          default=self.show_error,
                          save=True, auto_type=bool)

        # edges
        self._set_setting(settings, setting_keys, ['is_edges_visible'],
                          default=self.is_edges_visible,
                          save=True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['is_edges_black'],
                          default=self.is_edges_black,
                          save=True, auto_type=bool)

        self._set_setting(settings, setting_keys, ['is_horizontal_scalar_bar'],
                          default=self.is_horizontal_scalar_bar,
                          save=True, auto_type=bool)

        # min/max
        self._set_setting(settings, setting_keys, ['is_min_visible'],
                          default=self.is_min_visible,
                          save=True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['is_max_visible'],
                          default=self.is_max_visible,
                          save=True, auto_type=bool)

        # the vtk panel background color
        self._set_setting(settings, setting_keys, ['use_gradient_background'],
                          default=False, save=True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['background_color'],
                          default=BACKGROUND_COLOR, save=True, auto_type=float)
        self._set_setting(settings, setting_keys, ['background_color2'],
                          default=BACKGROUND_COLOR2, save=True, auto_type=float)
        self.background_color = force_color_ranged(self.background_color, BACKGROUND_COLOR)
        self.background_color2 = force_color_ranged(self.background_color2, BACKGROUND_COLOR2)

        # scales the coordinate systems
        self._set_setting(settings, setting_keys, ['coord_scale'],
                          default=COORD_SCALE, save=True, auto_type=float)
        self._set_setting(settings, setting_keys, ['coord_text_scale'],
                          default=COORD_TEXT_SCALE, save=True, auto_type=float)
        self.coord_scale = force_ranged(
            self.coord_scale,
            min_value=COORD_SCALE_MIN, max_value=COORD_SCALE_MAX)
        self.coord_text_scale = force_ranged(
            self.coord_text_scale,
            min_value=COORD_TEXT_SCALE_MIN, max_value=COORD_TEXT_SCALE_MAX)

        # this is for the 3d annotation
        self._set_setting(settings, setting_keys, ['annotation_color'],
                          default=ANNOTATION_COLOR, save=True, auto_type=float)
        self._set_setting(settings, setting_keys, ['annotation_size'],
                          default=ANNOTATION_SIZE, save=True, auto_type=int) # int
        self.annotation_color = force_color_ranged(self.annotation_color, ANNOTATION_COLOR)
        self.annotation_size = force_ranged(
            self.annotation_size,
            min_value=ANNOTATION_SIZE_MIN, max_value=ANNOTATION_SIZE_MAX)

        #elif isinstance(self.annotation_size, int):
            #pass
        #else:
            #print('annotation_size = ', self.annotation_size)

        self._set_setting(settings, setting_keys, ['magnify'],
                          default=self.magnify, save=True, auto_type=int)
        self.magnify = force_ranged(
            self.magnify, min_value=MAGNIFY_MIN, max_value=MAGNIFY_MAX)

        # this is the text in the lower left corner
        self._set_setting(settings, setting_keys, ['corner_text_color'],
                          default=CORNER_TEXT_COLOR, save=True, auto_type=float)
        self._set_setting(settings, setting_keys, ['corner_text_size'],
                          default=CORNER_TEXT_SIZE, save=True, auto_type=int)
        self.corner_text_color = force_color_ranged(self.corner_text_color, CORNER_TEXT_COLOR)
        self.corner_text_size = force_ranged(
            self.corner_text_size, min_value=CORNER_TEXT_SIZE_MIN, max_value=CORNER_TEXT_SIZE_MAX)

        # highlight
        self._set_setting(settings, setting_keys, ['highlight_color'],
                          default=HIGHLIGHT_COLOR, save=True, auto_type=float)
        self._set_setting(settings, setting_keys, ['highlight_opacity'],
                          default=HIGHLIGHT_OPACITY, save=True, auto_type=float)
        self.highlight_color = force_color_ranged(self.highlight_color, HIGHLIGHT_COLOR)
        self.highlight_opacity = force_ranged(
            self.highlight_opacity, min_value=OPACITY_MIN, max_value=OPACITY_MAX)

        self._set_setting(settings, setting_keys, ['highlight_point_size'],
                          default=HIGHLIGHT_POINT_SIZE, save=True, auto_type=float)
        self._set_setting(settings, setting_keys, ['highlight_line_width', 'highlight_line_thickness'],
                          default=HIGHLIGHT_LINE_WIDTH, save=True, auto_type=float)
        self.highlight_point_size = force_ranged(
            self.highlight_point_size, min_value=POINT_SIZE_MIN, max_value=POINT_SIZE_MAX)
        self.highlight_line_width = force_ranged(
            self.highlight_line_width,
            min_value=LINE_WIDTH_MIN, max_value=LINE_WIDTH_MAX)
        #self._set_setting(settings, setting_keys, ['highlight_style'],
                          #HIGHLIGHT_OPACITY, auto_type=float)

        # shear moment torque
        self._set_setting(settings, setting_keys, ['shear_moment_torque_color'],
                          default=SHEAR_MOMENT_TORQUE_COLOR, save=True, auto_type=float)
        self._set_setting(settings, setting_keys, ['shear_moment_torque_opacity'],
                          default=SHEAR_MOMENT_TORQUE_OPACITY, save=True, auto_type=float)
        self.shear_moment_torque_color = force_color_ranged(
            self.shear_moment_torque_color, SHEAR_MOMENT_TORQUE_COLOR)
        self.shear_moment_torque_opacity = force_ranged(
            self.shear_moment_torque_opacity, min_value=OPACITY_MIN, max_value=OPACITY_MAX)

        self._set_setting(settings, setting_keys, ['shear_moment_torque_point_size'],
                          default=SHEAR_MOMENT_TORQUE_POINT_SIZE, save=True, auto_type=float)
        self._set_setting(settings, setting_keys, ['shear_moment_torque_line_width'],
                          default=SHEAR_MOMENT_TORQUE_LINE_WIDTH, save=True, auto_type=float)
        self.shear_moment_torque_point_size = force_ranged(
            self.shear_moment_torque_point_size, min_value=POINT_SIZE_MIN, max_value=POINT_SIZE_MAX)
        self.shear_moment_torque_line_width = force_ranged(
            self.shear_moment_torque_line_width,
            min_value=LINE_WIDTH_MIN, max_value=LINE_WIDTH_MAX)

        # displacement_model_scale - unused
        self._set_setting(settings, setting_keys, ['displacement_model_scale'],
                          default=DISPLACEMENT_MODEL_SCALE, save=True, auto_type=float)

        # default colormap for legend
        self._set_setting(settings, setting_keys, ['colormap'], default=DEFAULT_COLORMAP, save=True)
        if self.colormap not in COLORMAPS:
            self.colormap = DEFAULT_COLORMAP
        # general gui sizing
        screen_shape = self._set_setting(settings, setting_keys, ['screen_shape'],
                                         default=screen_shape_default, save=False, auto_type=int)

        #try:
            #screen_shape = settings.value("screen_shape", screen_shape_default)
        #except (TypeError, AttributeError):
            #screen_shape = screen_shape_default

        #if 'recent_files' in setting_keys:
        recent_files_og = copy.deepcopy(self.recent_files)
        try:
            recent_files = settings.value("recent_files", default=recent_files_og)
        except (TypeError, AttributeError):
            recent_files = recent_files_og

        if recent_files is None:
            recent_files = []
            self.recent_files = recent_files

        if not isinstance(recent_files, int):
            # yeah seriously...it just returns 0
            try:
                recent_files2 = [(fname, fmt) for (fname, fmt) in recent_files
                                 if os.path.exists(fname)]
            except TypeError:
                raise TypeError(f'recent_files = {recent_files}')

            #  only save 10 files
            self.recent_files = recent_files2[:NFILES_TO_SAVE]

        self.recent_files = filter_recent_files(self.recent_files)
        self._load_nastran_settings(settings, setting_keys)
        self._load_other_settings(settings, setting_keys)

        #w = screen_shape.width()
        #h = screen_shape.height()
        #try:
        if screen_shape:
            self.parent.resize(screen_shape[0], screen_shape[1])
            #width, height = screen_shape

        pos = self._set_setting(settings, setting_keys, ['screen_position'],
                          default=None, save=False)
        #if pos is not None:
            #x = 1
            #qpos = parent.pos()
            #pos = qpos.x(), qpos.y()

        self.python_dock_visible = self._set_setting(
            settings, setting_keys, ['python_dock_visible'],
            default=False, save=False)
        self.log_dock_visible = self._set_setting(
            settings, setting_keys, ['log_dock_visible'],
            default=True, save=False)

        font = QtGui.QFont()
        font.setPointSize(self.font_size)
        self.parent.setFont(font)

        #if 0:
            #pos_default = 0, 0
            #pos = settings.value("pos", pos_default)
            #x_pos, y_pos = pos
            #print(pos)
            #self.mapToGlobal(QtCore.QPoint(pos[0], pos[1]))
            #y_pos = pos_default[0]
            #self.parent.setGeometry(x_pos, y_pos, width, height)
        #except TypeError:
            #self.resize(1100, 700)
        is_loaded = True
        return is_loaded

    def _load_other_settings(self, settings: QSettings,
                             setting_keys: list[str]) -> None:
        other_settings: OtherSettings = self.other_settings
        # self.cart3d_fluent_include = ()
        # self.cart3d_fluent_remove = ()
        for key in OTHER_STRING_KEYS:
            default = getattr(other_settings, key)
            value = self._set_setting(
                settings, setting_keys, [key],
                default=default, save=False, auto_type=str)
            setattr(other_settings, key, value)

        for key in OTHER_LIST_INT_KEYS:
            default = getattr(other_settings, key)
            #print(key, default)
            value = self._set_setting(
                settings, setting_keys, [key],
                default=default, save=False, auto_type=int)
            setattr(other_settings, key, tuple(value))

        for key in OTHER_LIST_STR_KEYS:
            default = getattr(other_settings, key)
            #print(key, default)
            value = self._set_setting(
                settings, setting_keys, [key],
                default=default, save=False, auto_type=str)
            setattr(other_settings, key, tuple(value))


    def _load_nastran_settings(self, settings: QSettings,
                               setting_keys: list[str]) -> None:
        """
        loads the settings from 'nastran_displacement' (or similar)
        and save it to 'nastran_settings.displacement'
        """
        nastran_settings: NastranSettings = self.nastran_settings
        #print('-----default------')
        #print(nastran_settings)
        for key in NASTRAN_BOOL_KEYS:
            # nastran_is_properties -> nastran, is_properties
            base, key2 = key.split('_', 1)

            # we get default from the nastran_settings
            default = getattr(nastran_settings, key2)

            # pull it from the QSettings
            value = self._set_setting(settings, setting_keys, [key],
                                      default, save=True, auto_type=bool)
            #print(f'key={key!r} key2={key2!r} default={default!r} value={value!r}')
            setattr(nastran_settings, key2, value)

        for key in NASTRAN_STR_KEYS:
            base, key2 = key.split('_', 1)
            default = getattr(nastran_settings, key2)
            value = self._set_setting(settings, setting_keys, [key],
                                      default, save=True, auto_type=str)
            setattr(nastran_settings, key2, value)

        for key in NASTRAN_COLOR_KEYS:
            # nastran_is_properties -> nastran, is_properties
            base, key2 = key.split('_', 1)

            # we get default from the nastran_settings
            default = getattr(nastran_settings, key2)

            # pull it from the QSettings
            value = self._set_setting(settings, setting_keys, [key],
                                      default, save=True, auto_type=float)
            #print(f'key={key!r} key2={key2!r} default={default!r} value={value!r}')
            setattr(nastran_settings, key2, value)

    def _set_setting(self, settings: QSettings, setting_keys: list[str],
                     setting_names: list[str], default: Any,
                     save: bool=True, auto_type=None) -> Any:
        """
        helper method for ``reapply_settings``
        """
        assert isinstance(save, bool), save
        set_name = setting_names[0]
        value = get_setting(settings, setting_keys, setting_names, default,
                            auto_type=auto_type)
        if save:
            setattr(self, set_name, value)
        return value

    def save(self, settings: QSettings,
             is_testing: bool=False) -> None:
        """saves the settings"""
        #if not is_testing:
        parent = self.parent
        if hasattr(parent, 'saveGeometry'):
            settings.setValue('main_window_geometry', parent.saveGeometry())
        if hasattr(parent, 'saveState'):
            settings.setValue('main_window_state', parent.saveState())

        # booleans
        settings.setValue('is_trackball_camera', self.is_trackball_camera)
        settings.setValue('use_parallel_projection', self.use_parallel_projection)
        settings.setValue('use_gradient_background', self.use_gradient_background)

        # startup directory
        settings.setValue('startup_directory', self.startup_directory)
        settings.setValue('use_startup_directory', self.use_startup_directory)

        settings.setValue('recent_files', self.recent_files[:NFILES_TO_SAVE])

        settings.setValue('use_old_sidebar_objects', self.use_old_sidebar_objects)
        settings.setValue('use_new_sidebar_objects', self.use_new_sidebar_objects)
        settings.setValue('use_new_terms', self.use_new_terms)

        # rgb tuple
        settings.setValue('background_color', self.background_color)
        settings.setValue('background_color2', self.background_color2)
        settings.setValue('annotation_color', self.annotation_color)
        settings.setValue('corner_text_color', self.corner_text_color)

        settings.setValue('highlight_color', self.highlight_color)
        settings.setValue('highlight_opacity', self.highlight_opacity)
        settings.setValue('highlight_point_size', self.highlight_point_size)
        settings.setValue('highlight_line_width', self.highlight_line_width)

        settings.setValue('shear_moment_torque_color', self.shear_moment_torque_color)
        settings.setValue('shear_moment_torque_opacity', self.shear_moment_torque_opacity)
        settings.setValue('shear_moment_torque_point_size', self.shear_moment_torque_point_size)
        settings.setValue('shear_moment_torque_line_width', self.shear_moment_torque_line_width)

        # float
        settings.setValue('displacement_model_scale', self.displacement_model_scale)

        # logging
        settings.setValue('show_info', self.show_info)
        settings.setValue('show_debug', self.show_debug)
        settings.setValue('show_command', self.show_command)
        settings.setValue('show_warning', self.show_warning)
        settings.setValue('show_error', self.show_error)

        # edges
        settings.setValue('is_edges_visible', self.is_edges_visible)
        settings.setValue('is_edge_black', self.is_edges_black)

        settings.setValue('is_horizontal_scalar_bar', self.is_horizontal_scalar_bar)

        # min/max
        settings.setValue('is_min_visible', self.is_min_visible)
        settings.setValue('is_max_visible', self.is_max_visible)

        # int
        settings.setValue('font_size', self.font_size)
        settings.setValue('annotation_size', self.annotation_size)
        settings.setValue('magnify', self.magnify)

        # float
        settings.setValue('corner_text_size', self.corner_text_size)
        settings.setValue('coord_scale', self.coord_scale)
        settings.setValue('coord_text_scale', self.coord_text_scale)

        # str
        settings.setValue('colormap', self.colormap)

        # format-specific
        nastran_settings = self.nastran_settings
        nastran_settings.save(settings)

        other_settings = self.other_settings
        other_settings.save(settings)

        #screen_shape = QtGui.QDesktopWidget().screenGeometry()

        # checks because tests don't have these
        if hasattr(self.parent, 'python_dock_widget'):
            python_dock_visible = self.parent.python_dock_widget.isVisible()
            settings.setValue('python_dock_visible', python_dock_visible)
        if hasattr(self.parent, 'log_dock_widget'):
            log_dock_widget = self.parent.log_dock_widget.isVisible()
            settings.setValue('log_dock_visible', log_dock_widget)

        if not is_testing:
            main_window = self.parent.window()
            width = main_window.frameGeometry().width()
            height = main_window.frameGeometry().height()
            settings.setValue('screen_shape', (width, height))

            qpos = parent.pos()
            pos = qpos.x(), qpos.y()
            settings.setValue('screen_position', pos)

    #---------------------------------------------------------------------------
    # FONT SIZE
    def on_increase_font_size(self):
        """shrinks the overall GUI font size"""
        self.on_set_font_size(self.font_size + 1)

    def on_decrease_font_size(self) -> None:
        """shrinks the overall GUI font size"""
        self.on_set_font_size(self.font_size - 1)

    def on_set_font_size(self, font_size: int, show_command: bool=True) -> None:
        """updates the GUI font size"""
        return self.parent.on_set_font_size(font_size, show_command=show_command)

    #---------------------------------------------------------------------------
    # ANNOTATION SIZE/COLOR
    def set_annotation_size_color(self, size: Optional[float]=None,
                                  color: Optional[ColorFloat]=None) -> None:
        """
        Parameters
        ----------
        size : float
            annotation size
        color : ColorFloat
            RGB values as floats

        """
        if size is not None:
            assert isinstance(size, int), 'size=%r' % size
            self.set_annotation_size(size)
        if color is not None:
            assert len(color) == 3, color
            assert isinstance(color[0], float), 'color=%r' % color
            self.set_annotation_color(color)

    def set_annotation_size(self, size: int, render: bool=True) -> None:
        """Updates the size of all the annotations"""
        assert size >= 0, size
        assert isinstance(size, int), size
        if self.annotation_size == size:
            return
        self.annotation_size = size

        # min/max
        for actor in self.parent.min_max_actors:
            actor.GetTextProperty().SetFontSize(size)
            actor.Modified()

        # case attached annotations (typical)
        for follower_actors in self.parent.label_actors.values():
            for follower_actor in follower_actors:
                follower_actor.GetTextProperty().SetFontSize(size)
                follower_actor.Modified()

        # geometry property attached annotations (e.g., flaps)
        for obj in self.parent.geometry_properties.values():
            if isinstance(obj, CoordProperties):
                continue
            elif isinstance(obj, AltGeometry):
                pass
            else:
                raise NotImplementedError(obj)

            follower_actors = obj.label_actors
            for follower_actor in follower_actors:
                follower_actor.GetTextProperty().SetFontSize(size)
                follower_actor.Modified()

        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()
            self.parent.log_command(f'self.settings.set_annotation_size({size})')

    def set_coord_scale(self, coord_scale: float, render: bool=True) -> None:
        """sets the coordinate system size"""
        self.coord_scale = coord_scale
        self.update_coord_scale(coord_scale, render=render)

    def set_coord_text_scale(self, coord_text_scale: float,
                             render: bool=True) -> None:
        """sets the coordinate system text size"""
        self.coord_text_scale = coord_text_scale
        self.update_coord_text_scale(coord_text_scale, render=render)

    def update_coord_scale(self, coord_scale=None, coord_text_scale=None,
                           linewidth=None, render: bool=True) -> None:
        """internal method for updating the coordinate system size"""
        if coord_scale is None:
            coord_scale = self.coord_scale
        #if coord_text_scale:
            #self.update_coord_text_scale(coord_text_scale=coord_text_scale, render=False)

        dim_max = self.dim_max
        scale = coord_scale * dim_max

        for unused_coord_id, axes in self.parent.axes.items():
            axes.SetTotalLength(scale, scale, scale) # was coord_scale
            #axes.SetScale(magnify, magnify, magnify)
            #if linewidth:
                #xaxis = axes.GetXAxisShaftProperty()
                #yaxis = axes.GetXAxisShaftProperty()
                #zaxis = axes.GetXAxisShaftProperty()
                #lw = xaxis.GetLineWidth()  #  1.0
                #xaxis.SetLineWidth(linewidth)
                #yaxis.SetLineWidth(linewidth)
                #zaxis.SetLineWidth(linewidth)
            #print(f'coord_scale coord_id={unused_coord_id} scale={scale} lw={linewidth}')

        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()

    def scale_coord(self, magnify: float, render: bool=True) -> None:
        """internal method for scaling the coordinate system size"""
        for unused_coord_id, axes in self.parent.axes.items():
            axes.SetScale(magnify)
        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()

    def update_coord_text_scale(self, coord_text_scale: Optional[float]=None,
                                render: bool=True) -> None:
        """internal method for updating the coordinate system size"""
        if coord_text_scale is None:
            coord_text_scale = self.coord_text_scale

        update_axes_text_size(self.parent.axes, coord_text_scale,
                              width=1.0, height=0.25)
        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()

    def set_annotation_color(self, color: ColorFloat,
                             render: bool=True) -> None:
        """
        Set the annotation color

        Parameters
        ----------
        color : ColorFloat
            RGB values as floats
        """
        if np.allclose(self.annotation_color, color):
            return
        self.annotation_color = color

        # min/max
        for min_max_actor in self.parent.min_max_actors:
            #print(dir(min_max_actor))
            prop = min_max_actor.GetTextProperty()  # was GetProperty
            prop.SetColor(*color)

        # case attached annotations (typical)
        for follower_actors in self.parent.label_actors.values():
            for follower_actor in follower_actors:
                prop = follower_actor.GetTextProperty()
                prop.SetColor(*color)

        # geometry property attached annotations (e.g., flaps)
        for obj in self.parent.geometry_properties.values():
            if isinstance(obj, CoordProperties):
                continue
            elif isinstance(obj, AltGeometry):
                pass
            else:  # pragma: no cover
                raise NotImplementedError(obj)

            follower_actors = obj.label_actors
            for follower_actor in follower_actors:
                prop = follower_actor.GetTextProperty()
                prop.SetColor(*color)

        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()
        self.parent.log_command('self.settings.set_annotation_color(%s, %s, %s)' % color)

    #---------------------------------------------------------------------------
    def set_background_color_to_white(self, render: bool=True) -> None:
        """sets the background color to white; used by gif writing?"""
        self.set_gradient_background(use_gradient_background=False, render=False)
        self.set_background_color(WHITE_FLOAT, render=render)

    def set_gradient_background(self,
                                use_gradient_background: bool=False,
                                render: bool=True,
                                quiet: bool=False) -> None:
        """enables/disables the gradient background"""
        self.use_gradient_background = use_gradient_background
        self.parent.rend.SetGradientBackground(self.use_gradient_background)
        if render:
            self.parent.vtk_interactor.Render()

    def set_background_color(self, color: ColorFloat,
                             render: bool=True, quiet: bool=False) -> None:
        """
        Set the background color

        Parameters
        ----------
        color : ColorFloat
            RGB values as floats
        """
        self.background_color = color
        self.parent.rend.SetBackground(*color)
        if render:
            self.parent.vtk_interactor.Render()
        if not quiet:
            self.parent.log_command('self.settings.set_background_color(%s, %s, %s)' % color)

    def set_background_color2(self, color: ColorFloat,
                              render: bool=True, quiet: bool=False):
        """
        Set the background color

        Parameters
        ----------
        color : ColorFloat
            RGB values as floats
        """
        self.background_color2 = color
        self.parent.rend.SetBackground2(*color)
        if render:
            self.parent.vtk_interactor.Render()
        if not quiet:
            self.parent.log_command('self.settings.set_background_color2(%s, %s, %s)' % color)

    def set_highlight_color(self, color: ColorFloat, render: bool=True) -> None:
        """
        Set the highlight color

        Parameters
        ----------
        color : ColorFloat
            RGB values as floats
        """
        self.highlight_color = color
        if render:
            self.parent.vtk_interactor.Render()
        self.parent.log_command('self.settings.set_highlight_color(%s, %s, %s)' % color)

    def set_highlight_opacity(self, opacity: float) -> None:
        """
        Set the highlight opacity

        Parameters
        ----------
        opacity : float
            0.0 : invisible
            1.0 : solid
        """
        self.highlight_opacity = opacity
        self.parent.log_command(f'self.settings.set_highlight_opacity({opacity})')

    def set_highlight_point_size(self, point_size: int) -> None:
        """
        Set the highlight point size

        Parameters
        ----------
        point_size : float
            the point size
        """
        self.highlight_point_size = point_size
        self.parent.log_command(f'self.settings.set_highlight_point_size({point_size})')

    #---------------------------------------------------------------------------
    # TEXT ACTORS - used for lower left notes

    def set_corner_text_color(self, color: list[float],
                              render: bool=True) -> None:
        """
        Set the corner_text color

        Parameters
        ----------
        color : (float, float, float)
            RGB values as floats
        """
        self.text_color = color
        text_actors =  self.parent.corner_text_actors
        for text_actor in text_actors.values():
            text_actor.GetTextProperty().SetColor(color)
        if render:
            self.parent.vtk_interactor.Render()
        self.parent.log_command('self.settings.set_corner_text_color(%s, %s, %s)' % color)

    def set_corner_text_size(self, corner_text_size: int, render: bool=True) -> None:
        """
        Set the corner text size

        Parameters
        ----------
        corner_text_size : int
            the lower left text size (typical 14)

        """
        # we built these actors in reverse order,
        # so that's how we update their sizes
        text_actors =  self.parent.corner_text_actors
        i = len(text_actors) - 1
        dtext_size = corner_text_size + 1
        self.text_size = corner_text_size
        for text_actor in text_actors.values():
            text_prop = text_actor.GetTextProperty()
            text_prop.SetFontSize(corner_text_size)

            position = [5, 5 + i * dtext_size]
            text_actor.SetDisplayPosition(*position)
            i -= 1
        if render:
            self.parent.vtk_interactor.Render()
        self.parent.log_command(f'self.settings.set_corner_text_size({corner_text_size})')

    def update_corner_text_size(self, magnify: float=1.0) -> None:
        """Internal method for updating the bottom-left text when we go to take a picture"""
        text_size = int(14 * magnify)
        text_actors =  self.parent.corner_text_actors
        for text_actor in text_actors.values():
            text_prop = text_actor.GetTextProperty()
            text_prop.SetFontSize(text_size)

    def set_magnify(self, magnify: int=5) -> None:
        """sets the screenshot magnification factor"""
        self.magnify = magnify

    def set_trackball_camera(self, is_trackball_camera: bool,
                             render: bool=True) -> None:
        """sets the trackball_camera flag"""
        self.is_trackball_camera = is_trackball_camera
        self.parent.mouse_actions.set_style()
        if render:
            self.parent.vtk_interactor.Render()

    def set_parallel_projection(self, parallel_projection: bool,
                                render: bool=True) -> None:
        """sets the parallel_projection flag"""
        self.use_parallel_projection = parallel_projection
        camera = self.parent.rend.GetActiveCamera()
        if parallel_projection:
            camera.ParallelProjectionOn()
        else:
            camera.ParallelProjectionOff()
        if render:
            self.parent.vtk_interactor.Render()

    def __repr__(self) -> str:
        msg = '<Settings>\n'
        for key in object_attributes(self, mode='public', keys_to_skip=['parent']):
            value = getattr(self, key)
            if isinstance(value, tuple):
                value = str(value)
            msg += '  %r = %r\n' % (key, value)
        return msg

def update_axes_text_size(axes: dict[int, vtkAxes],
                          coord_text_scale: float,
                          width: float=1.0, height: float=0.25):
    """updates the coordinate system text size"""
    # width doesn't set the width
    # it being very large (old=0.1) makes the width constraint inactive

    for unused_coord_id, axis in axes.items():
        #print(f'coord_text_scale coord_id={unused_coord_id} coord_text_scale={coord_text_scale}')
        texts = [
            axis.GetXAxisCaptionActor2D(),
            axis.GetYAxisCaptionActor2D(),
            axis.GetZAxisCaptionActor2D(),
        ]
        # this doesn't set the width
        # this being very large (old=0.1) makes the width constraint inactive
        for text in texts:
            text.SetWidth(coord_text_scale * width)
            text.SetHeight(coord_text_scale * height)

def isfloat(value) -> bool:
    """is the value floatable"""
    try:
        float(value)
        return True
    except ValueError:
        return False

def repr_settings(settings: QSettings) -> str:
    """works on a QSettings, not a Settings"""
    msg = 'QSettings:\n'
    for key in sorted(settings.allKeys()):
        value = settings.value(key)
        msg += '    %r : %r\n' % (key, value)
    return msg

def filter_recent_files(recent_files: list[tuple[str, str]]):
    if IS_WINDOWS:
        #recent_files_lower = [fname.lower() for (fname, fmt) in recent_files]]
        recent_files_out = []
        stored_files_lower = set([])
        for i, (fname, geometry_format) in enumerate(recent_files):
            if geometry_format is None or not os.path.exists(fname):
                continue

            fname_abs = os.path.abspath(fname)
            fname_lower = fname.lower()
            if fname_lower in stored_files_lower:
                continue
            stored_files_lower.add(fname_lower)
            recent_files_out.append((fname_lower, geometry_format))
    else:
        recent_files_out = [(fname, fmt) for (fname, fmt) in recent_files
                            if os.path.exists(fname) and fmt is not None]
    return recent_files_out

def force_ranged(value, min_value=None, max_value=None):
    """make sure a value is in the proper range"""
    if min_value is not None and max_value is not None:
        out = max(min(value, max_value), min_value)
    elif min_value is not None:
        out = max(value, min_value)
    elif max_value is not None:
        out = max(value, min_value)
    else:  # pragma: no cover
        raise RuntimeError(value)
    #if out != value:
        #print(out, value)
    return out

def force_color_ranged(color: ColorFloat,
                       default_color: ColorFloat) -> ColorFloat:
    """
    make sure a color is in the proper range
    default if it's out of range
    """
    assert isinstance(color, tuple), color
    if min(color) < 0.0 or max(color) > 1.0:
        return default_color
    return color
