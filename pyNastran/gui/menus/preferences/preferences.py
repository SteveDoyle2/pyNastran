"""
The preferences menu handles:
 - Font Size
 - Background Color
 - Text Color
 - Annotation Color
 - Annotation Size
 - Clipping Min
 - Clipping Max

"""
from __future__ import annotations
from math import log10, ceil
from functools import partial
from typing import Optional, Any, TYPE_CHECKING

from qtpy import QtGui
from qtpy.QtWidgets import (
    QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
    QSpinBox, QDoubleSpinBox, QColorDialog, QLineEdit, QCheckBox,
    QTabWidget, QWidget, QComboBox,
)

from pyNastran.utils.locale import func_str, float_locale
from pyNastran.gui.utils.qt.pydialog import PyDialog, make_font, check_color
from pyNastran.gui.utils.qt.qcombobox import make_combo_box, get_combo_box_text
from pyNastran.gui.utils.qt.qpush_button_color import QPushButtonColor
from pyNastran.gui.utils.qt.checks.qlineedit import QLINEEDIT_GOOD, QLINEEDIT_ERROR

from pyNastran.gui.menus.menu_utils import eval_float_from_string
from pyNastran.gui.gui_objects.settings import (
    FONT_SIZE, FONT_SIZE_MIN, FONT_SIZE_MAX,
    MAGNIFY,
    COORD_SCALE, COORD_TEXT_SCALE,
    BACKGROUND_COLOR, BACKGROUND_COLOR2,
    ANNOTATION_COLOR, ANNOTATION_SIZE,
    CORNER_TEXT_COLOR, CORNER_TEXT_SIZE,
    HIGHLIGHT_COLOR, HIGHLIGHT_OPACITY, HIGHLIGHT_POINT_SIZE, HIGHLIGHT_LINE_WIDTH,
    SHEAR_MOMENT_TORQUE_COLOR, SHEAR_MOMENT_TORQUE_OPACITY, SHEAR_MOMENT_TORQUE_POINT_SIZE, SHEAR_MOMENT_TORQUE_LINE_WIDTH,
    OPACITY_MIN, OPACITY_MAX,
    USE_PARALLEL_PROJECTION, IS_TRACKBALL_CAMERA,
    NASTRAN_BOOL_KEYS,
    POINT_SIZE_MIN, POINT_SIZE_MAX,
    COORD_TEXT_SCALE_MIN, COORD_TEXT_SCALE_MAX,
    CORNER_TEXT_SIZE_MIN, CORNER_TEXT_SIZE_MAX,

    COORD_SCALE_MIN, COORD_SCALE_MAX,
    MAGNIFY_MIN, MAGNIFY_MAX,
    ANNOTATION_SIZE_MIN, ANNOTATION_SIZE_MAX,
    NASTRAN_VERSIONS,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui_objects.settings import Settings, NastranSettings
    from pyNastran.gui.typing import ColorInt


USE_TABS = True
IS_SMT = False
UNITS_MODEL_IN = [
    'in-lbf-s', 'm-kg-s', 'unitless',
]
LENGTH_UNITS = ['in', 'ft', 'm', 'cm', 'mm']
PRESSURE_UNITS = ['psi', 'ksi', 'Pa', 'MPa']
STRESS_UNITS = PRESSURE_UNITS
FORCE_UNITS = ['lbf', 'N', 'kN', 'MN', 'mN']
MOMENT_UNITS = ['in-lbf', 'ft-lbf', 'N-m']
AREA_UNITS = [f'{unit}^2' for unit in LENGTH_UNITS]

DISPLAMCENT_UNITS = LENGTH_UNITS
VELOCITY_UNITS = [f'{unit}/s' for unit in LENGTH_UNITS]
ACCELERATION_UNITS = [f'{unit}/s^2' for unit in LENGTH_UNITS] + ['g']


class PreferencesWindow(PyDialog):
    """
    +-------------+
    | Preferences |
    +---------------------------------+
    | Text Size        ______ Default |
    | Annotation Color ______         |
    | Annotation Size  ______         |
    | Picker Size      ______         |
    | Back Color       ______         |
    | Text Color       ______         |
    |                                 |
    |            Reset Defaults       |
    |        Apply OK Cancel          |
    +---------------------------------+
    """
    def __init__(self, data, win_parent=None):
        """
        Saves the data members from data and
        performs type checks
        """
        PyDialog.__init__(self, data, win_parent)

        self._updated_preference = False

        self.dim_max = data['dim_max']

        # font size for menu
        self._default_font_size = data['font_size']

        # corner text
        self._default_corner_text_size = CORNER_TEXT_SIZE

        self.use_startup_directory = data['use_startup_directory']

        # an annotation is the marked/probe label
        self._default_annotation_size = ANNOTATION_SIZE

        self._default_coord_scale = COORD_SCALE * 100. # * self.dim_max
        self._default_coord_text_scale = COORD_TEXT_SCALE * 100. # * self.dim_max
        self._default_clipping_min = data['min_clip']
        self._default_clipping_max = data['max_clip']
        #self._default_annotation_size = data['annotation_size'] # int
        #self.default_magnify = data['magnify']

        self._cart3d_fluent_include = data['cart3d_fluent_include']
        self._cart3d_fluent_remove = data['cart3d_fluent_remove']

        self._units_model_in = data['units_model_in']
        self._units_length = data['units_length']
        #self._units_area = data['units_area']
        self._units_force = data['units_force']
        self._units_moment = data['units_moment']
        self._units_pressure = data['units_pressure']
        self._units_stress = data['units_stress']
        self._units_displacement = data['units_displacement']
        self._units_velocity = data['units_velocity']
        self._units_acceleration = data['units_acceleration']

        self._is_trackball_camera = data['is_trackball_camera']
        self._parallel_projection = data['use_parallel_projection']
        self._use_gradient_background = data['use_gradient_background'] # bool
        self._show_corner_coord = data['show_corner_coord']
        self._annotation_size = data['annotation_size'] # int

        #self.out_data = data

        # doesn't include dim_max
        self._picker_size = data['picker_size'] * 100.
        self._coord_scale = data['coord_scale'] * 100.
        self._coord_text_scale = data['coord_text_scale'] * 100.

        self._magnify = data['magnify']
        self._corner_text_size = data['corner_text_size']
        self._highlight_opacity = data['highlight_opacity']
        self._highlight_point_size = data['highlight_point_size']

        self.annotation_color_float, self.annotation_color_int = check_color(
            data['annotation_color'])
        self.background_color_float, self.background_color_int = check_color(
            data['background_color'])
        self.background_color2_float, self.background_color2_int = check_color(
            data['background_color2'])
        self.corner_text_color_float, self.corner_text_color_int = check_color(
            data['corner_text_color'])
        self.highlight_color_float, self.highlight_color_int = check_color(
            data['highlight_color'])

        self.caero_color_float, self.caero_color_int = check_color(
            data['caero_color'])
        self.rbe_line_color_float, self.rbe_line_color_int = check_color(
            data['rbe_line_color'])
        self.plotel_color_float, self.plotel_color_int = check_color(
            data['plotel_color'])

        #self._shear_moment_torque_opacity = data['shear_moment_torque_opacity']
        #self._shear_moment_torque_point_size = data['shear_moment_torque_point_size']
        #self._shear_moment_torque_color_int = data['shear_moment_torque_color']
        #self._shear_moment_torque_line_thickness = data['shear_moment_torque_line_thickness']

        self._nastran_version = data['nastran_version']
        self._nastran_is_element_quality = data['nastran_is_element_quality']
        self._nastran_is_properties = data['nastran_is_properties']
        self._nastran_is_3d_bars = data['nastran_is_3d_bars']
        self._nastran_is_3d_bars_update = data['nastran_is_3d_bars_update']
        self._nastran_is_mass_update = data['nastran_is_mass_update']
        self._nastran_is_constraints = data['nastran_is_constraints']
        self._nastran_is_bar_axes = data['nastran_is_bar_axes']
        self._nastran_create_coords = data['nastran_create_coords']
        self._nastran_is_shell_mcids = data['nastran_is_shell_mcids']
        self._nastran_is_rbe = data['nastran_is_rbe']
        self._nastran_is_aero = data['nastran_is_aero']
        self._nastran_is_plotel = data['nastran_is_plotel']

        self._nastran_stress = data['nastran_stress']
        self._nastran_plate_stress = data['nastran_plate_stress']
        self._nastran_composite_plate_stress = data['nastran_composite_plate_stress']
        self._nastran_strain = data['nastran_strain']
        self._nastran_plate_strain = data['nastran_plate_strain']
        self._nastran_composite_plate_strain = data['nastran_composite_plate_strain']
        self._nastran_rod_stress = data['nastran_rod_stress']
        self._nastran_bar_stress = data['nastran_bar_stress']
        self._nastran_beam_stress = data['nastran_beam_stress']
        self._nastran_rod_strain = data['nastran_rod_strain']
        self._nastran_bar_strain = data['nastran_bar_strain']
        self._nastran_beam_strain = data['nastran_beam_strain']

        self._nastran_displacement = data['nastran_displacement']
        self._nastran_velocity = data['nastran_velocity']
        self._nastran_acceleration = data['nastran_acceleration']
        self._nastran_eigenvector = data['nastran_eigenvector']
        self._nastran_temperature = data['nastran_temperature']

        self._nastran_spc_force = data['nastran_spc_force']
        self._nastran_mpc_force = data['nastran_mpc_force']
        self._nastran_applied_load = data['nastran_applied_load']
        self._nastran_heat_flux = data['nastran_heat_flux']

        self._nastran_force = data['nastran_force']
        self._nastran_grid_point_force = data['nastran_grid_point_force']
        self._nastran_strain_energy = data['nastran_strain_energy']

        self.setWindowTitle('Preferences')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        self.on_font(self.font_size)
        self.on_gradient_scale()
        #self.show()

    def create_widgets(self):
        """creates the display window"""
        # window text size
        self.font_size_label = QLabel('Font Size:')
        self.font_size_edit = QSpinBox(self)
        self.font_size_edit.setValue(self._default_font_size)
        self.font_size_edit.setRange(FONT_SIZE_MIN, FONT_SIZE_MAX)

        #-----------------------------------------------------------------------
        # Corner Text Color
        self.corner_text_size_label = QLabel('Corner Text Size:')
        self.corner_text_size_edit = QSpinBox(self)
        self.corner_text_size_edit.setValue(self._default_corner_text_size)
        self.corner_text_size_edit.setRange(CORNER_TEXT_SIZE_MIN, CORNER_TEXT_SIZE_MAX)
        self.corner_text_size_edit.setToolTip('Sets the lower left corner text size')
        self.corner_text_size_button = QPushButton("Default")

        # Text Color
        self.corner_text_color_label = QLabel("Corner Text Color:")
        self.corner_text_color_edit = QPushButtonColor(self.corner_text_color_int)
        self.corner_text_color_edit.setToolTip('Sets the lower left corner text color')

        #-----------------------------------------------------------------------
        # Highlight Color
        self.highlight_opacity_label = QLabel("Highlight Opacity:")
        self.highlight_opacity_edit = create_opacity_edit(self, self._highlight_opacity)
        self.highlight_opacity_edit.setToolTip('Sets the highlight opacity (0=invisible, 1=solid)')
        self.highlight_opacity_button = QPushButton("Default")

        self.highlight_point_size_label = QLabel("Highlight Point Size:")
        self.highlight_point_size_edit = create_point_size_edit(self, self._highlight_point_size)
        self.highlight_point_size_edit.setToolTip('Sets the highlight node size')
        self.highlight_point_size_button = QPushButton("Default")

        # Text Color
        self.highlight_color_label = QLabel("Highlight Color:")
        self.highlight_color_edit = QPushButtonColor(self.highlight_color_int)
        self.highlight_color_edit.setToolTip('Sets the highlight color')

        #-----------------------------------------------------------------------
        # shear_moment_torque Color
        if IS_SMT:
            self.shear_moment_torque_label = QLabel("Shear-Moment-Torque:")

            opacity_edit, point_size_edit, line_width_edit, color_edit = create_shear_moment_torque_edits(
                self,
                self._shear_moment_torque_opacity,
                self._shear_moment_torque_point_size,
                self._shear_moment_torque_line_thickness,
                self._shear_moment_torque_color_int)
            self.shear_moment_torque_opacity_edit = opacity_edit
            self.shear_moment_torque_point_size_edit = point_size_edit
            self.shear_moment_torque_line_width_edit = line_width_edit
            self.shear_moment_torque_color_edit = color_edit

            self.shear_moment_torque_opacity_label = QLabel("Opacity:")
            self.shear_moment_torque_opacity_button = QPushButton("Default")

            self.shear_moment_torque_point_size_label = QLabel("Point Size:")
            self.shear_moment_torque_point_size_button = QPushButton("Default")

            self.shear_moment_torque_line_width_label = QLabel("Line Width:")
            self.shear_moment_torque_line_width_button = QPushButton("Default")

            # Text Color
            self.shear_moment_torque_color_label = QLabel("Color:")

        #-----------------------------------------------------------------------
        # Background Color
        self.background_color_label = QLabel("Btm Background Color:")
        self.background_color_edit = QPushButtonColor(self.background_color_int)
        self.background_color_edit.setToolTip('Sets the lower background color')

        # Background Color2
        self.gradient_scale_label = QLabel("Gradient Background:")
        self.gradient_scale_checkbox = QCheckBox()
        self.gradient_scale_checkbox.setChecked(self._use_gradient_background)

        self.background_color2_label = QLabel("Top Background Color:")
        self.background_color2_edit = QPushButtonColor(self.background_color2_int)
        self.background_color2_edit.setToolTip('Sets the upper background color')

        #-----------------------------------------------------------------------
        # Annotation Size
        self.annotation_size_label = QLabel("Annotation Size:")
        self.annotation_size_edit = QSpinBox(self)
        self.annotation_size_edit.setRange(ANNOTATION_SIZE_MIN, ANNOTATION_SIZE_MAX)
        self.annotation_size_edit.setValue(self._annotation_size)
        self.annotation_size_edit.setToolTip('Sets the "Probe" and Min/Max text size')
        #self.annotation_size_edit.setToolTip('Sets the hiannotation text size')
        self.annotation_size_button = QPushButton("Default")

        # Annotation Color - unused
        self.annotation_color_label = QLabel("Annotation Color:")
        self.annotation_color_edit = QPushButtonColor(self.annotation_color_int)
        self.annotation_color_edit.setToolTip('The "Probe" is an annotation')
        #self.annotation_color_label.hide()
        #self.annotation_color_edit.hide()

        #-----------------------------------------------------------------------
        # Picker Sizef
        self.picker_size_label = QLabel("Picker Size (% of Screen):")
        self.picker_size_edit = QDoubleSpinBox(self)
        self.picker_size_edit.setRange(0., 10.)

        log_dim = log10(self.dim_max)
        decimals = int(ceil(abs(log_dim)))

        decimals = max(6, decimals)
        self.picker_size_edit.setDecimals(decimals)
        self.picker_size_edit.setSingleStep(10. / 5000.)
        self.picker_size_edit.setValue(self._picker_size)
        self.picker_size_edit.setToolTip('Sets the picker tolerance')

        self.parallel_projection_label = QLabel('Parallel Projection:')
        self.parallel_projection_checkbox = QCheckBox()
        self.parallel_projection_checkbox.setChecked(self._parallel_projection)
        self.parallel_projection_checkbox.setToolTip('Checked: Typical engineering perspectivfe (default)\n'
                                                     'Unchecked: Distort the model like a real camera')

        self.trackball_camera_label = QLabel('Trackball Camera:')
        self.trackball_camera_checkbox = QCheckBox()
        self.trackball_camera_checkbox.setChecked(self._is_trackball_camera)
        self.trackball_camera_checkbox.setToolTip('Checked: Trackball Camera (default)\n'
                                                  'Unchecked: Joystick Camera (for 3d mice)')

        #-----------------------------------------------------------------------
        # Clipping Min
        self.clipping_min_label = QLabel('Clipping Min:')
        self.clipping_min_edit = QLineEdit(func_str(self._default_clipping_min))
        self.clipping_min_button = QPushButton('Default')

        # Clipping Max
        self.clipping_max_label = QLabel('Clipping Max:')
        self.clipping_max_edit = QLineEdit(func_str(self._default_clipping_max))
        self.clipping_max_button = QPushButton('Default')

        #-----------------------------------------------------------------------
        self.coord_scale_label = QLabel('Coordinate System Scale:')
        self.coord_scale_button = QPushButton("Default")

        self.coord_scale_edit = QDoubleSpinBox(self)
        self.coord_scale_edit.setRange(COORD_SCALE_MIN, COORD_SCALE_MAX)
        self.coord_scale_edit.setDecimals(3)
        self.coord_scale_edit.setSingleStep(1.0)
        self.coord_scale_edit.setValue(self._coord_scale)
        self.picker_size_edit.setToolTip('Increase/decrease the coordinate system size')

        self.coord_text_scale_label = QLabel('Coordinate System Text Scale:')
        self.coord_text_scale_button = QPushButton("Default")

        self.coord_text_scale_edit = QDoubleSpinBox(self)
        self.coord_text_scale_edit.setRange(COORD_TEXT_SCALE_MIN, COORD_TEXT_SCALE_MAX)
        self.coord_text_scale_edit.setDecimals(3)
        self.coord_text_scale_edit.setSingleStep(2.)
        self.coord_text_scale_edit.setValue(self._coord_text_scale)
        self.picker_size_edit.setToolTip('Increase/decrease the coordinate system text size')

        # Show corner coord
        self.corner_coord_label = QLabel('Show Corner Coordinate System:')
        self.corner_coord_checkbox = QCheckBox()
        self.corner_coord_checkbox.setChecked(self._show_corner_coord)

        #-----------------------------------------------------------------------
        self.magnify_label = QLabel('Screenshot Magnify:')
        self.magnify_edit = QSpinBox(self)
        self.magnify_edit.setMinimum(MAGNIFY_MIN)
        self.magnify_edit.setMaximum(MAGNIFY_MAX)
        self.magnify_edit.setValue(self._magnify)
        self.magnify_edit.setToolTip('1: Standard resolution; >1: high quality')
        self._set_widgets_nastran()
        self._set_widgets_other()

    def _set_widgets_nastran(self):
        self.nastran_version_label = QLabel('Version')
        self.nastran_version_pulldown = QComboBox(self)
        self.nastran_version_pulldown.addItems(NASTRAN_VERSIONS)
        nastran_versions_lower = [version.lower() for version in NASTRAN_VERSIONS]
        iversion = nastran_versions_lower.index(self._nastran_version.lower())
        self.nastran_version_pulldown.setItemText(iversion, self._nastran_version)

        self.nastran_is_element_quality_checkbox = QCheckBox('Element Quality')
        self.nastran_is_element_quality_checkbox.setToolTip('Cacluate Aspect Ratio, Skew Angle, Max/Min Interior Angle, etc.')
        self.nastran_is_element_quality_checkbox.setChecked(self._nastran_is_element_quality)

        self.nastran_is_properties_checkbox = QCheckBox('Properties')
        self.nastran_is_properties_checkbox.setToolTip('Breakdown each layer of a PCOMP/PSHELL')
        self.nastran_is_properties_checkbox.setChecked(self._nastran_is_properties)

        self.nastran_is_3d_bars_checkbox = QCheckBox('3D Bars')
        self.nastran_is_3d_bars_checkbox.setToolTip('Crete 3D Bar/Beam geometry')
        self.nastran_is_3d_bars_checkbox.setChecked(self._nastran_is_3d_bars)
        #self.nastran_is_3d_bars_checkbox.setDisabled(True)

        self.nastran_is_3d_bars_update_checkbox = QCheckBox('Update 3D Bars')
        self.nastran_is_3d_bars_update_checkbox.setToolTip('Update the 3D geometry (Bar/Beam cross-sections and CONMw) when deformations are applied')
        self.nastran_is_3d_bars_update_checkbox.setChecked(self._nastran_is_3d_bars_update)

        self.nastran_is_mass_update_checkbox = QCheckBox('Update Masses')
        self.nastran_is_mass_update_checkbox.setToolTip('Update the masses when nodes change')
        self.nastran_is_mass_update_checkbox.setChecked(self._nastran_is_mass_update)

        self.nastran_is_constraints_checkbox = QCheckBox('Constraints')
        self.nastran_is_constraints_checkbox.setToolTip('Create actors for the constraints (SPC, MPC, SUPORT)')
        self.nastran_is_constraints_checkbox.setChecked(self._nastran_is_constraints)


        self.nastran_is_shell_mcids_checkbox = QCheckBox('Shell MCIDs')
        self.nastran_is_shell_mcids_checkbox.setToolTip('Calculate the Material Coordinate Systems for Shells')
        self.nastran_is_shell_mcids_checkbox.setChecked(self._nastran_is_shell_mcids)

        self.nastran_is_rbe_checkbox = QCheckBox('RBEs')
        self.nastran_is_rbe_checkbox.setToolTip('Create MPC/RBE2/RBE3 dependent and indepdent nodes and lines')
        self.nastran_is_rbe_checkbox.setChecked(self._nastran_is_rbe)

        self.nastran_is_aero_checkbox = QCheckBox('Aero')
        self.nastran_is_aero_checkbox.setToolTip('Create aero panel (CAERO/SPLINE/SET) visualization')
        self.nastran_is_aero_checkbox.setChecked(self._nastran_is_aero)

        self.nastran_is_plotel_checkbox = QCheckBox('PLOTELs')
        self.nastran_is_plotel_checkbox.setToolTip('Create PLOTELs')
        self.nastran_is_plotel_checkbox.setChecked(self._nastran_is_plotel)

        self.nastran_create_coords_checkbox = QCheckBox('Coords')
        self.nastran_create_coords_checkbox.setChecked(self._nastran_create_coords)

        self.nastran_is_bar_axes_checkbox = QCheckBox('Bar Axes')
        self.nastran_is_bar_axes_checkbox.setChecked(self._nastran_is_bar_axes)
        #self.nastran_is_bar_axes_checkbox.setDisabled(True)

        if 1:
            self.nastran_displacement_checkbox = QCheckBox('Displacement')
            self.nastran_velocity_checkbox = QCheckBox('Velocity')
            self.nastran_acceleration_checkbox = QCheckBox('Acceleration')
            self.nastran_eigenvector_checkbox = QCheckBox('Eigenvector')
            self.nastran_temperature_checkbox = QCheckBox('Temperature')
            self.nastran_displacement_checkbox.setChecked(self._nastran_displacement)
            self.nastran_velocity_checkbox.setChecked(self._nastran_velocity)
            self.nastran_acceleration_checkbox.setChecked(self._nastran_acceleration)
            self.nastran_eigenvector_checkbox.setChecked(self._nastran_eigenvector)
            self.nastran_temperature_checkbox.setChecked(self._nastran_temperature)

            self.nastran_spc_force_checkbox = QCheckBox('SPC Force')
            self.nastran_mpc_force_checkbox = QCheckBox('MPC Force')
            self.nastran_applied_load_checkbox = QCheckBox('Applied Load')
            self.nastran_heat_flux_checkbox = QCheckBox('Heat Flux')
            self.nastran_spc_force_checkbox.setChecked(self._nastran_spc_force)
            self.nastran_mpc_force_checkbox.setChecked(self._nastran_mpc_force)
            self.nastran_applied_load_checkbox.setChecked(self._nastran_applied_load)
            self.nastran_heat_flux_checkbox.setChecked(self._nastran_heat_flux)

            self.nastran_force_checkbox = QCheckBox('Force')
            self.nastran_grid_point_force_checkbox = QCheckBox('Grid Point Force')
            self.nastran_strain_energy_checkbox = QCheckBox('Strain Energy')

            self.nastran_force_checkbox.setChecked(self._nastran_force)
            self.nastran_grid_point_force_checkbox.setChecked(self._nastran_grid_point_force)
            self.nastran_strain_energy_checkbox.setChecked(self._nastran_strain_energy)

            self.nastran_stress_checkbox = QCheckBox('Stress')
            #self.nastran_plate_stress_checkbox = QCheckBox('Plate Stress')
            #self.nastran_composite_plate_stress_checkbox = QCheckBox('Composite Plate Stress')
            #self.nastran_rod_stress_checkbox = QCheckBox('Rod Stress')
            #self.nastran_bar_stress_checkbox = QCheckBox('Bar Stress')
            #self.nastran_beam_stress_checkbox = QCheckBox('Beam Stress')

            self.nastran_stress_checkbox.setChecked(self._nastran_stress)
            #self.nastran_plate_stress_checkbox.setChecked(self._nastran_plate_stress)
            #self.nastran_composite_plate_stress_checkbox.setChecked(self._nastran_composite_plate_stress)
            #self.nastran_rod_stress_checkbox.setChecked(self._nastran_rod_stress)
            #self.nastran_bar_stress_checkbox.setChecked(self._nastran_bar_stress)
            #self.nastran_beam_stress_checkbox.setChecked(self._nastran_beam_stress)

            self.nastran_strain_checkbox = QCheckBox('Strain')
            #self.nastran_plate_strain_checkbox = QCheckBox('Plate Strain')
            #self.nastran_composite_plate_strain_checkbox = QCheckBox('Composite Plate Strain')
            #self.nastran_rod_strain_checkbox = QCheckBox('Rod Strain')
            #self.nastran_bar_strain_checkbox = QCheckBox('Bar Strain')
            #self.nastran_beam_strain_checkbox = QCheckBox('Beam Strain')

            self.nastran_strain_checkbox.setChecked(self._nastran_strain)
            #self.nastran_plate_strain_checkbox.setChecked(self._nastran_plate_strain)
            #self.nastran_composite_plate_strain_checkbox.setChecked(self._nastran_composite_plate_strain)
            #self.nastran_rod_strain_checkbox.setChecked(self._nastran_rod_strain)
            #self.nastran_bar_strain_checkbox.setChecked(self._nastran_bar_strain)
            #self.nastran_beam_strain_checkbox.setChecked(self._nastran_beam_strain)

        #-----------------------------------------------------------------------
        # colors
        self.caero_color_label = QLabel("Default CAERO panel color:")
        self.caero_color_edit = QPushButtonColor(self.caero_color_int)
        self.caero_color_edit.setToolTip('Sets the color for the caero/caero sub-panels actors')

        self.rbe_line_color_label = QLabel("Default RBE2/RBE3 line color:")
        self.rbe_line_color_edit = QPushButtonColor(self.rbe_line_color_int)
        self.rbe_line_color_edit.setToolTip('Sets the color for the RBE2/RBE3 lines')

        self.plotel_color_label = QLabel("Default PLOTEL color:")
        self.plotel_color_edit = QPushButtonColor(self.plotel_color_int)
        self.plotel_color_edit.setToolTip('Sets the color for the PLOTELs')

        #-----------------------------------------------------------------------
        # closing
        self.reset_defaults_button = QPushButton('Reset Defaults')
        self.apply_button = QPushButton('Apply')
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Cancel')

    def _set_widgets_other(self):
        self.cart3d_fluent_label = QLabel('Cart3d/Fluent')
        self.region_include_label = QLabel('Regions (Include):')
        self.region_remove_label = QLabel('Regions (Remove):')
        include_str = ' '.join(str(val) for val in self._cart3d_fluent_include)
        remove_str = ' '.join(str(val) for val in self._cart3d_fluent_remove)
        self.cart3d_fluent_regions_include = QLineEdit(include_str)
        self.cart3d_fluent_regions_remove = QLineEdit(remove_str)

        self.units_label = QLabel('Units:')
        self.model_in_label = QLabel('Input Units:')
        self.length_label = QLabel('Length:')

        self.stress_label = QLabel('Stress:')
        self.pressure_label = QLabel('Pressure:')
        self.force_label = QLabel('Force:')
        self.moment_label = QLabel('Moment:')
        self.displacement_label = QLabel('Displacement:')
        self.velocity_label = QLabel('Velocity:')
        self.acceleration_label = QLabel('Acceleration:')
        self.length_pulldown = make_combo_box(LENGTH_UNITS, self._units_length, partial(self.on_unit, 'length'))
        self.stress_pulldown = make_combo_box(STRESS_UNITS, self._units_stress, partial(self.on_unit, 'stress'))
        self.pressure_pulldown = make_combo_box(PRESSURE_UNITS, self._units_pressure, partial(self.on_unit, 'pressure'))
        #self.area_pulldown = make_combo_box(AREA_UNITS, self._units_area, partial(self.on_unit, 'area'))
        self.force_pulldown = make_combo_box(FORCE_UNITS, self._units_force, partial(self.on_unit, 'force'))
        self.moment_pulldown = make_combo_box(MOMENT_UNITS, self._units_moment, partial(self.on_unit, 'moment'))
        self.displacement_pulldown = make_combo_box(DISPLAMCENT_UNITS, self._units_displacement, partial(self.on_unit, 'displacement'))
        self.velocity_pulldown = make_combo_box(VELOCITY_UNITS, self._units_velocity, partial(self.on_unit, 'velocity'))
        self.acceleration_pulldown = make_combo_box(ACCELERATION_UNITS, self._units_acceleration, partial(self.on_unit, 'acceleration'))

        pulldowns = [
            #self.area_labe, self.area_pulldown,
            self.force_label, self.force_pulldown,
            self.moment_label, self.moment_pulldown,
            self.pressure_label, self.pressure_pulldown,
            self.stress_label, self.stress_pulldown,
            self.displacement_label, self.displacement_pulldown,
            self.velocity_label, self.velocity_pulldown,
            self.acceleration_label, self.acceleration_pulldown,
        ]
        for pulldown in pulldowns:
            pulldown.setVisible(False)

        #('in', 'lbf', 's', 'psi')
        #self.units_model_in = ('unitless','','','')
        self.units_model_label = QLabel('Model Units:')
        units_model_in_str = '-'.join(self._units_model_in[:3]).rstrip('-')
        self.units_model_pulldown = make_combo_box(UNITS_MODEL_IN, units_model_in_str,
                                                   self.on_units_model_in)

    #def create_legend_widgets(self):
        #"""
        #Creates the widgets for the legend control

        #Name    Itailic  Bold     Font
        #====    =======  =====  ========
        #Title    check   check  pulldown
        #Label    check   check  pulldown
        #"""
        #self.name_label = QLabel("Name:")
        #self.italic_label = QLabel("Italic:")
        #self.bold_label = QLabel("Bold:")
        #self.font_label = QLabel("Font:")
        #self.legend_label = QLabel("Legend:")

        #self.legend_title_name = QLabel("Title")
        #self.legend_title_italic_check = QCheckBox()
        #self.legend_title_bold_check = QCheckBox()
        #self.legend_title_font_edit = QComboBox()
        #self.legend_title_font_edit.addItems(['cat', 'dog', 'frog'])

        #self.legend_label_italic_name = QLabel("Label")
        #self.legend_label_italic_check = QCheckBox()
        #self.legend_label_bold_check = QCheckBox()
        #self.legend_label_font_edit = QComboBox()
        #self.legend_label_font_edit.addItems(['cat2', 'dog2', 'frog2'])

    #def create_legend_layout(self):
        #"""
        #Creates the layout for the legend control

        #Name    Italic  Bold     Font
        #====    ======  =====  ========
        #Title    check  check  pulldown
        #Label    check  check  pulldown
        #"""
        #grid2 = QGridLayout()
        #grid2.addWidget(self.legend_label, 0, 0)

        #grid2.addWidget(self.name_label, 1, 0)
        #grid2.addWidget(self.italic_label, 1, 1)
        #grid2.addWidget(self.bold_label, 1, 2)
        #grid2.addWidget(self.font_label, 1, 3)

        #grid2.addWidget(self.legend_title_name, 2, 0)
        #grid2.addWidget(self.legend_title_italic_check, 2, 1)
        #grid2.addWidget(self.legend_title_bold_check, 2, 2)
        #grid2.addWidget(self.legend_title_font_edit, 2, 3)

        #grid2.addWidget(self.legend_label_italic_name, 3, 0)
        #grid2.addWidget(self.legend_label_italic_check, 3, 1)
        #grid2.addWidget(self.legend_label_bold_check, 3, 2)
        #grid2.addWidget(self.legend_label_font_edit, 3, 3)
        #return grid2

    def create_layout(self):
        grid = QGridLayout()

        irow = 0
        grid.addWidget(self.font_size_label, irow, 0)
        grid.addWidget(self.font_size_edit, irow, 1)
        irow += 1

        grid.addWidget(self.gradient_scale_label, irow, 0)
        grid.addWidget(self.gradient_scale_checkbox, irow, 1)
        irow += 1

        grid.addWidget(self.background_color2_label, irow, 0)
        grid.addWidget(self.background_color2_edit, irow, 1)
        irow += 1

        grid.addWidget(self.background_color_label, irow, 0)
        grid.addWidget(self.background_color_edit, irow, 1)
        irow += 1

        grid.addWidget(self.highlight_color_label, irow, 0)
        grid.addWidget(self.highlight_color_edit, irow, 1)
        irow += 1

        grid.addWidget(self.highlight_opacity_label, irow, 0)
        grid.addWidget(self.highlight_opacity_edit, irow, 1)
        grid.addWidget(self.highlight_opacity_button, irow, 2)
        irow += 1

        grid.addWidget(self.highlight_point_size_label, irow, 0)
        grid.addWidget(self.highlight_point_size_edit, irow, 1)
        grid.addWidget(self.highlight_point_size_button, irow, 2)
        irow += 1

        grid.addWidget(self.corner_text_color_label, irow, 0)
        grid.addWidget(self.corner_text_color_edit, irow, 1)
        irow += 1

        grid.addWidget(self.corner_text_size_label, irow, 0)
        grid.addWidget(self.corner_text_size_edit, irow, 1)
        grid.addWidget(self.corner_text_size_button, irow, 2)
        irow += 1

        grid.addWidget(self.annotation_color_label, irow, 0)
        grid.addWidget(self.annotation_color_edit, irow, 1)
        irow += 1

        grid.addWidget(self.annotation_size_label, irow, 0)
        grid.addWidget(self.annotation_size_edit, irow, 1)
        grid.addWidget(self.annotation_size_button, irow, 2)
        irow += 1

        #grid.addWidget(self.clipping_min_label, irow, 0)
        #grid.addWidget(self.clipping_min_edit, irow, 1)
        #grid.addWidget(self.clipping_min_button, irow, 2)
        #irow += 1

        #grid.addWidget(self.clipping_max_label, irow, 0)
        #grid.addWidget(self.clipping_max_edit, irow, 1)
        #grid.addWidget(self.clipping_max_button, irow, 2)
        #irow += 1

        grid.addWidget(self.corner_coord_label, irow, 0)
        grid.addWidget(self.corner_coord_checkbox, irow, 1)
        irow += 1

        grid.addWidget(self.coord_scale_label, irow, 0)
        grid.addWidget(self.coord_scale_edit, irow, 1)
        grid.addWidget(self.coord_scale_button, irow, 2)
        irow += 1

        grid.addWidget(self.coord_text_scale_label, irow, 0)
        grid.addWidget(self.coord_text_scale_edit, irow, 1)
        grid.addWidget(self.coord_text_scale_button, irow, 2)
        irow += 1

        #-----------------------------------------------
        grid.addWidget(self.magnify_label, irow, 0)
        grid.addWidget(self.magnify_edit, irow, 1)
        irow += 1

        grid.addWidget(self.picker_size_label, irow, 0)
        grid.addWidget(self.picker_size_edit, irow, 1)
        irow += 1

        grid.addWidget(self.parallel_projection_label, irow, 0)
        grid.addWidget(self.parallel_projection_checkbox, irow, 1)
        irow += 1

        grid.addWidget(self.trackball_camera_label, irow, 0)
        grid.addWidget(self.trackball_camera_checkbox, irow, 1)
        irow += 1

        #--------------------------------------------------

        out = self._get_nastran_vboxs()
        vbox_nastran, vbox_nastran_results, vbox_nastran_actors = out

        #self.create_legend_widgets()
        #grid2 = self.create_legend_layout()
        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        if USE_TABS:
            vbox = QVBoxLayout()
            tabs = QTabWidget(self)
            general_tab_widget = QWidget(self)
            general_tab_widget.setLayout(grid)

            vbox_nastran_tab = QVBoxLayout()
            vbox_nastran_tab.addLayout(vbox_nastran)
            vbox_nastran_tab.addLayout(vbox_nastran_results)
            vbox_nastran_tab.addLayout(vbox_nastran_actors)
            vbox_nastran_tab.addStretch()

            vbox_other = self._get_grid_other()
            vbox_fluent_cart3d = QVBoxLayout()
            vbox_fluent_cart3d.addLayout(vbox_other)

            vbox_other_tab = QVBoxLayout()
            vbox_other_tab.addLayout(vbox_fluent_cart3d)
            #vbox_other_tab.addLayout(vbox_nastran_results)
            #vbox_other_tab.addLayout(vbox_nastran_actors)
            vbox_other_tab.addStretch()

            nastran_tab_widget = QWidget(self)
            nastran_tab_widget.setLayout(vbox_nastran_tab)

            other_tab_widget = QWidget(self)
            other_tab_widget.setLayout(vbox_other_tab)

            tabs.addTab(general_tab_widget, 'General')
            tabs.addTab(nastran_tab_widget, 'Nastran')
            tabs.addTab(other_tab_widget, 'Other')
            vbox.addWidget(tabs)

        else:
            vbox = QVBoxLayout()
            vbox.addLayout(grid)
            vbox.addLayout(vbox_nastran)
            vbox.addLayout(vbox_nastran_results)
            vbox.addLayout(vbox_nastran_actors)
        #vbox.addStretch()
        #vbox.addLayout(grid2)
        vbox.addStretch()

        vbox.addWidget(self.reset_defaults_button)
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def _get_grid_other(self):
        vbox_other = QGridLayout()
        irow = 1
        vbox_other.addWidget(self.cart3d_fluent_label, irow, 0)
        irow += 1
        vbox_other.addWidget(self.region_include_label, irow, 0)
        vbox_other.addWidget(self.cart3d_fluent_regions_include, irow, 1)
        irow += 1
        vbox_other.addWidget(self.region_remove_label, irow, 0)
        vbox_other.addWidget(self.cart3d_fluent_regions_remove, irow, 1)


        irow += 1
        vbox_other.addWidget(self.units_label, irow, 0)
        irow += 1
        vbox_other.addWidget(self.units_model_label, irow, 0)
        vbox_other.addWidget(self.units_model_pulldown, irow, 1)
        irow += 1
        vbox_other.addWidget(self.length_label, irow, 0)
        vbox_other.addWidget(self.length_pulldown, irow, 1)
        irow += 1
        vbox_other.addWidget(self.stress_label, irow, 0)
        vbox_other.addWidget(self.stress_pulldown, irow, 1)
        irow += 1
        vbox_other.addWidget(self.pressure_label, irow, 0)
        vbox_other.addWidget(self.pressure_pulldown, irow, 1)
        irow += 1
        vbox_other.addWidget(self.force_label, irow, 0)
        vbox_other.addWidget(self.force_pulldown, irow, 1)
        irow += 1
        vbox_other.addWidget(self.moment_label, irow, 0)
        vbox_other.addWidget(self.moment_pulldown, irow, 1)

        irow += 1
        vbox_other.addWidget(self.displacement_label, irow, 0)
        vbox_other.addWidget(self.displacement_pulldown, irow, 1)
        irow += 1
        vbox_other.addWidget(self.velocity_label, irow, 0)
        vbox_other.addWidget(self.velocity_pulldown, irow, 1)
        irow += 1
        vbox_other.addWidget(self.acceleration_label, irow, 0)
        vbox_other.addWidget(self.acceleration_pulldown, irow, 1)
        return vbox_other

    def _get_nastran_vboxs(self) -> tuple[QVBoxLayout, QVBoxLayout, QVBoxLayout]:
        grid_nastran = self._get_grid_nastran_layout()
        grid_nastran_results = self._get_grid_nastran_results_layout()

        #bold_font = make_font(self.font_size, is_bold=True)
        vbox_nastran = QVBoxLayout()
        self.nastran_label = QLabel('Nastran Geometry:')
        vbox_nastran.addWidget(self.nastran_label)
        vbox_nastran.addLayout(grid_nastran)

        vbox_nastran_results = QVBoxLayout()
        self.nastran_results_label = QLabel('Nastran Results:')
        vbox_nastran_results.addWidget(self.nastran_results_label)
        vbox_nastran_results.addLayout(grid_nastran_results)

        grid_nastran_actors = QGridLayout()
        irow = 1
        grid_nastran_actors.addWidget(self.caero_color_label, irow, 0)
        grid_nastran_actors.addWidget(self.caero_color_edit, irow, 1)
        irow += 1
        grid_nastran_actors.addWidget(self.rbe_line_color_label, irow, 0)
        grid_nastran_actors.addWidget(self.rbe_line_color_edit, irow, 1)
        irow += 1
        grid_nastran_actors.addWidget(self.plotel_color_label, irow, 0)
        grid_nastran_actors.addWidget(self.plotel_color_edit, irow, 1)

        vbox_nastran_actors = QVBoxLayout()
        self.nastran_actors_label = QLabel('Nastran Actors:')
        vbox_nastran_actors.addWidget(self.nastran_actors_label)
        vbox_nastran_actors.addLayout(grid_nastran_actors)
        return vbox_nastran, vbox_nastran_results, vbox_nastran_actors

    def _get_grid_nastran_layout(self) -> QGridLayout:
        grid_nastran = QGridLayout()
        irow = 0

        grid_nastran.addWidget(self.nastran_version_label, irow, 0)
        grid_nastran.addWidget(self.nastran_version_pulldown, irow, 1)
        irow += 1

        grid_nastran.addWidget(self.nastran_create_coords_checkbox, irow, 0)
        irow += 1

        grid_nastran.addWidget(self.nastran_is_element_quality_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_is_properties_checkbox, irow, 1)
        irow += 1

        grid_nastran.addWidget(self.nastran_is_constraints_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_is_rbe_checkbox, irow, 1)
        irow += 1

        grid_nastran.addWidget(self.nastran_is_shell_mcids_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_is_aero_checkbox, irow, 1)
        irow += 1

        grid_nastran.addWidget(self.nastran_is_bar_axes_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_is_3d_bars_checkbox, irow, 1)
        grid_nastran.addWidget(self.nastran_is_3d_bars_update_checkbox, irow, 2)
        irow += 1
        grid_nastran.addWidget(self.nastran_is_mass_update_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_is_plotel_checkbox, irow, 1)
        return grid_nastran

    def _get_grid_nastran_results_layout(self) -> QGridLayout:
        grid_nastran = QGridLayout()
        irow = 0
        # ------------------------
        grid_nastran.addWidget(self.nastran_displacement_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_velocity_checkbox, irow, 1)
        grid_nastran.addWidget(self.nastran_acceleration_checkbox, irow, 2)
        irow += 1

        grid_nastran.addWidget(self.nastran_eigenvector_checkbox, irow, 0)
        irow += 1

        # ------------------------
        grid_nastran.addWidget(self.nastran_spc_force_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_mpc_force_checkbox, irow, 1)
        grid_nastran.addWidget(self.nastran_applied_load_checkbox, irow, 2)
        irow += 1

        grid_nastran.addWidget(self.nastran_temperature_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_heat_flux_checkbox, irow, 1)
        irow += 1

        # ------------------------
        grid_nastran.addWidget(self.nastran_stress_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_strain_checkbox, irow, 1)
        grid_nastran.addWidget(self.nastran_force_checkbox, irow, 2)
        irow += 1

        grid_nastran.addWidget(self.nastran_strain_energy_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_grid_point_force_checkbox, irow, 1)
        irow += 1

        if IS_SMT:
            grid_nastran.addWidget(self.shear_moment_torque_label, irow, 0)
            irow += 1

            grid_nastran.addWidget(self.shear_moment_torque_point_size_label, irow, 0)
            grid_nastran.addWidget(self.shear_moment_torque_point_size_edit, irow, 1)
            grid_nastran.addWidget(self.shear_moment_torque_point_size_button, irow, 2)
            irow += 1

            grid_nastran.addWidget(self.shear_moment_torque_line_width_label, irow, 0)
            grid_nastran.addWidget(self.shear_moment_torque_line_width_edit, irow, 1)
            grid_nastran.addWidget(self.shear_moment_torque_line_width_button, irow, 2)
            irow += 1

            grid_nastran.addWidget(self.shear_moment_torque_opacity_label, irow, 0)
            grid_nastran.addWidget(self.shear_moment_torque_opacity_edit, irow, 1)
            grid_nastran.addWidget(self.shear_moment_torque_opacity_button, irow, 2)
            irow += 1

            grid_nastran.addWidget(self.shear_moment_torque_color_label, irow, 0)
            grid_nastran.addWidget(self.shear_moment_torque_color_edit, irow, 1)
            irow += 1



        #self.nastran_plate_stress_checkbox = QCheckBox('Plate Stress')
        #self.nastran_composite_plate_stress_checkbox = QCheckBox('Composite Plate Stress')
        #self.nastran_rod_stress_checkbox = QCheckBox('Rod Stress')
        #self.nastran_bar_stress_checkbox = QCheckBox('Bar Stress')
        #self.nastran_beam_stress_checkbox = QCheckBox('Beam Stress')
        #
        #self.nastran_plate_strain_checkbox = QCheckBox('Plate Strain')
        #self.nastran_composite_plate_strain_checkbox = QCheckBox('Composite Plate Strain')
        #self.nastran_rod_strain_checkbox = QCheckBox('Rod Strain')
        #self.nastran_bar_strain_checkbox = QCheckBox('Bar Strain')
        #self.nastran_beam_strain_checkbox = QCheckBox('Beam Strain')
        return grid_nastran

    def _set_other_connections(self):
        #self.cart3d_fluent_regions_include
        pass

    def _set_nastran_connections(self):
        # format-specific
        self.nastran_version_pulldown.currentIndexChanged.connect(self.on_nastran_version)
        self.nastran_is_element_quality_checkbox.clicked.connect(partial(on_nastran, self, 'is_element_quality'))
        self.nastran_is_properties_checkbox.clicked.connect(partial(on_nastran, self, 'is_properties'))
        self.nastran_is_3d_bars_checkbox.clicked.connect(partial(on_nastran, self, 'is_3d_bars'))
        self.nastran_is_3d_bars_update_checkbox.clicked.connect(partial(on_nastran, self, 'is_3d_bars_update'))
        self.nastran_is_mass_update_checkbox.clicked.connect(partial(on_nastran, self, 'is_mass_update'))

        self.nastran_is_bar_axes_checkbox.clicked.connect(partial(on_nastran, self, 'is_bar_axes'))
        self.nastran_create_coords_checkbox.clicked.connect(partial(on_nastran, self, 'create_coords'))
        self.nastran_is_shell_mcids_checkbox.clicked.connect(partial(on_nastran, self, 'is_shell_mcids'))
        self.nastran_is_aero_checkbox.clicked.connect(partial(on_nastran, self, 'is_aero'))
        self.nastran_is_plotel_checkbox.clicked.connect(partial(on_nastran, self, 'is_plotel'))
        self.nastran_is_rbe_checkbox.clicked.connect(partial(on_nastran, self, 'is_rbe'))
        self.nastran_is_constraints_checkbox.clicked.connect(partial(on_nastran, self, 'is_constraints'))

        #self.nastran_is_shell_mcid_checkbox.clicked.connect(self.on_nastran_is_shell_mcids)
        #self.nastran_is_shell_mcid_checkbox.clicked.connect(self.on_nastran_is_shell_mcids2)

        self.nastran_displacement_checkbox.clicked.connect(partial(on_nastran, self, 'displacement'))
        self.nastran_velocity_checkbox.clicked.connect(partial(on_nastran, self, 'acceleration'))
        self.nastran_acceleration_checkbox.clicked.connect(partial(on_nastran, self, 'acceleration'))
        self.nastran_eigenvector_checkbox.clicked.connect(partial(on_nastran, self, 'eigenvector'))
        self.nastran_temperature_checkbox.clicked.connect(partial(on_nastran, self, 'temperature'))

        self.nastran_spc_force_checkbox.clicked.connect(partial(on_nastran, self, 'spc_force'))
        self.nastran_mpc_force_checkbox.clicked.connect(partial(on_nastran, self, 'mpc_force'))
        self.nastran_applied_load_checkbox.clicked.connect(partial(on_nastran, self, 'applied_load'))
        self.nastran_grid_point_force_checkbox.clicked.connect(partial(on_nastran, self, 'grid_point_force'))

        self.nastran_force_checkbox.clicked.connect(partial(on_nastran, self, 'force'))
        self.nastran_strain_checkbox.clicked.connect(partial(on_nastran, self, 'strain'))
        self.nastran_stress_checkbox.clicked.connect(partial(on_nastran, self, 'stress'))
        self.nastran_strain_energy_checkbox.clicked.connect(partial(on_nastran, self, 'strain_energy'))

        # --------------------------
        # colors
        self.caero_color_edit.clicked.connect(self.on_caero_color)
        self.rbe_line_color_edit.clicked.connect(self.on_rbe_line_color)
        self.plotel_color_edit.clicked.connect(self.on_plotel_color)

    def set_connections(self):
        """creates the actions for the menu"""
        self.font_size_edit.valueChanged.connect(self.on_font)

        self.annotation_size_edit.editingFinished.connect(self.on_annotation_size)
        self.annotation_size_edit.valueChanged.connect(self.on_annotation_size)
        self.annotation_color_edit.clicked.connect(self.on_annotation_color)
        self.annotation_size_button.clicked.connect(self.on_default_annotation_size)

        self.background_color_edit.clicked.connect(self.on_background_color)
        self.background_color2_edit.clicked.connect(self.on_background_color2)
        self.gradient_scale_checkbox.clicked.connect(self.on_gradient_scale)

        self.highlight_color_edit.clicked.connect(self.on_highlight_color)
        self.highlight_opacity_edit.valueChanged.connect(self.on_highlight_opacity)
        self.highlight_point_size_edit.valueChanged.connect(self.on_highlight_point_size)

        self.corner_text_color_edit.clicked.connect(self.on_corner_text_color)
        self.corner_text_size_edit.valueChanged.connect(self.on_corner_text_size)
        self.corner_text_size_button.clicked.connect(self.on_default_corner_text_size)

        self.picker_size_edit.valueChanged.connect(self.on_picker_size)
        self.picker_size_edit.editingFinished.connect(self.on_picker_size)

        self.parallel_projection_checkbox.clicked.connect(self.on_parallel_projection)
        self.trackball_camera_checkbox.clicked.connect(self.on_trackball_camera)

        self.coord_scale_edit.valueChanged.connect(self.on_coord_scale)
        self.coord_scale_edit.editingFinished.connect(self.on_coord_scale)
        self.coord_scale_button.clicked.connect(self.on_default_coord_scale)
        self.corner_coord_checkbox.clicked.connect(self.on_corner_coord)

        self.coord_text_scale_edit.valueChanged.connect(self.on_coord_text_scale)
        self.coord_text_scale_edit.editingFinished.connect(self.on_coord_text_scale)
        self.coord_text_scale_button.clicked.connect(self.on_default_coord_text_scale)

        self.magnify_edit.valueChanged.connect(self.on_magnify)
        self.magnify_edit.editingFinished.connect(self.on_magnify)

        self.clipping_min_button.clicked.connect(self.on_default_clipping_min)
        self.clipping_max_button.clicked.connect(self.on_default_clipping_max)

        #------------------------------------
        self._set_nastran_connections()
        self._set_other_connections()
        #------------------------------------

        self.reset_defaults_button.clicked.connect(self.on_reset_defaults)
        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)
        # closeEvent

    def on_reset_defaults(self):
        """reset all of the preferences to their defaults"""
        self.on_font(FONT_SIZE)
        self.on_default_annotation_size()
        self.on_default_clipping_max()
        self.on_default_clipping_min()
        self.on_default_coord_scale()
        self.on_default_coord_text_scale()
        self.on_default_corner_text_size()

        self.magnify_edit.setValue(MAGNIFY)
        self.picker_size_edit.setValue(self._picker_size)

        self.highlight_opacity_edit.setValue(HIGHLIGHT_OPACITY)
        self.highlight_point_size_edit.setValue(HIGHLIGHT_POINT_SIZE)
        self.parallel_projection_checkbox.setChecked(USE_PARALLEL_PROJECTION)
        self.trackball_camera_checkbox.setChecked(IS_TRACKBALL_CAMERA)
        #self.highlight_line_width_edit.setValue(HIGHLIGHT_LINE_WIDTH)

        self.background_color1_float = BACKGROUND_COLOR
        self.background_color2_float = BACKGROUND_COLOR2
        self.highlight_color_float = HIGHLIGHT_COLOR
        self.corner_text_color_float = CORNER_TEXT_COLOR
        self.annotation_color_float = ANNOTATION_COLOR
        self.shear_moment_torque_color_float = SHEAR_MOMENT_TORQUE_COLOR

        self.gradient_scale_checkbox.setChecked(True)
        self.on_gradient_scale()
        self.corner_coord_checkbox.setChecked(True)

        self.background_color1_int = tuple([round(val * 255) for val in BACKGROUND_COLOR])
        self.background_color2_int = tuple([round(val * 255) for val in BACKGROUND_COLOR2])
        self.highlight_color_int = tuple([round(val * 255) for val in HIGHLIGHT_COLOR])
        self.corner_text_color_int = tuple([round(val * 255) for val in CORNER_TEXT_COLOR])
        self.annotation_color_int = tuple([round(val * 255) for val in ANNOTATION_COLOR])
        self.shear_moment_torque_color_int = tuple([round(val * 255) for val in SHEAR_MOMENT_TORQUE_COLOR])

        set_pushbutton_color(self.corner_text_color_edit, self.corner_text_color_int)
        set_pushbutton_color(self.highlight_color_edit, self.highlight_color_int)
        set_pushbutton_color(self.background_color_edit, self.background_color1_int)
        set_pushbutton_color(self.background_color2_edit, self.background_color2_int)
        set_pushbutton_color(self.annotation_color_edit, self.annotation_color_int)
        #set_pushbutton_color(self.shear_moment_torque_color_edit, self.shear_moment_torque_color_int)

        for key in NASTRAN_BOOL_KEYS:
            checkbox_name = f'{key}_checkbox'
            if hasattr(self, checkbox_name):
                checkbox = getattr(self, checkbox_name)
                checkbox.setChecked(True)

        if self.win_parent is not None:
            settings: Settings = self.settings
            is_min = settings.is_min_visible
            is_max = settings.is_max_visible
            is_edges_black = settings.is_edges_black
            is_edges_visible = settings.is_edges_visible
            settings.reset_settings(resize=False, reset_dim_max=False)

            self.win_parent.on_set_edge_visibility(is_edges_black, render=False)
            self.win_parent.on_set_edge_visibility(is_edges_black, render=False)

            settings.set_highlight_color(self.highlight_color_float, render=False)
            settings.set_corner_text_color(self.corner_text_color_float, render=False)
            settings.set_background_color(self.background_color1_float, render=False)
            settings.set_background_color2(self.background_color2_float, render=True)
        self.on_apply()

    def on_units_model_in(self):
        text = get_combo_box_text(self.units_model_pulldown)
        units_model_in_sline0 = text.split('-')
        units_model_in_sline = [''] * 4
        units_model_in_sline[:len(units_model_in_sline0)] = units_model_in_sline0
        if self.win_parent is not None:
            settings: Settings = self.settings
            other_settings = settings.other_settings
            assert len(units_model_in_sline) == len(other_settings.units_model_in)
            #print(f'on_units_model_in: units_model_in_sline={units_model_in_sline}')
            other_settings.units_model_in = tuple(units_model_in_sline)

    def on_unit(self, name: str) -> None:
        pulldown = getattr(self, f'{name}_pulldown')
        text = get_combo_box_text(pulldown)
        if self.win_parent is not None:
            settings: Settings = self.settings
            other_settings = settings.other_settings
            setattr(other_settings, f'units_{name}', text)
            #print(f'set units_{name} = {text}')

    @property
    def settings(self) -> Settings:
        return self.win_parent.settings

    @property
    def nastran_settings(self) -> NastranSettings:
        return self.settings.nastran_settings

    def on_nastran_version(self) -> None:
        version = get_combo_box_text(self.nastran_version_pulldown).lower()
        #iversion = self.nastran_version_pulldown.currentIndex()
        #version = NASTRAN_VERSIONS[iversion].lower()
        self.nastran_settings.version = version

    #def on_nastran_is_shell_mcids2(self):
        #"""set the nastran properties preferences"""
        #is_checked = self.nastran_is_shell_mcid_checkbox.isChecked()
        #if self.win_parent is not None:
            #self.nastran_settings.is_shell_mcids = is_checked

    def on_font(self, value=None):
        """update the font for the current window"""
        if value in (0, None):
            value = self.font_size_edit.value()
        font = make_font(value, is_bold=False)
        self.setFont(font)
        bold_font = make_font(value, is_bold=True)

        bold_labels = [
            self.nastran_label, self.nastran_results_label,
            self.nastran_actors_label,
            self.cart3d_fluent_label, self.units_label,
        ]
        if IS_SMT:
            bold_labels.append(self.shear_moment_torque_label)
        for label in bold_labels:
            label.setFont(bold_font)

    def on_annotation_size(self, value=None) -> None:
        """update the annotation size"""
        if value is None:
            value = int(self.annotation_size_edit.text())
        self._annotation_size = value
        #self.on_apply(force=True)
        #self.min_edit.setText(func_str(self._default_min))
        #self.min_edit.setStyleSheet(QLINEEDIT_GOOD)
        self.update_annotation_size_color()

    def update_annotation_size_color(self) -> None:
        if self.win_parent is not None:
            self.settings.set_annotation_size_color(
                size=self._annotation_size,
                color=self.annotation_color_float)

    def on_gradient_scale(self):
        is_checked = self.gradient_scale_checkbox.isChecked()
        self.background_color2_label.setEnabled(is_checked)
        self.background_color2_edit.setEnabled(is_checked)
        if self.win_parent is not None:
            self.settings.set_gradient_background(use_gradient_background=is_checked)

    def on_corner_coord(self):
        is_checked = self.corner_coord_checkbox.isChecked()
        if self.win_parent is not None:
            self.win_parent.set_corner_axis_visiblity(is_checked, render=True)

    def on_annotation_color(self):
        rgb_color_ints = self.annotation_color_int
        title = "Choose an annotation color"
        passed, rgb_color_ints, rgb_color_floats = self.on_color(
            self.annotation_color_edit, rgb_color_ints, title)
        if passed:
            self.annotation_color_int = rgb_color_ints
            self.annotation_color_float = rgb_color_floats
            self.update_annotation_size_color()

    def on_background_color(self):
        """ Choose a background color"""
        title = "Choose a primary background color"
        rgb_color_ints = self.background_color_int
        color_edit = self.background_color_edit
        func_name = 'set_background_color'
        passed, rgb_color_ints, rgb_color_floats = self._load_color(
            title, color_edit, rgb_color_ints, func_name)
        if passed:
            self.background_color_int = rgb_color_ints
            self.background_color_float = rgb_color_floats

    def on_background_color2(self):
        """ Choose a background color"""
        title = "Choose a secondary background color"
        rgb_color_ints = self.background_color2_int
        color_edit = self.background_color2_edit
        func_name = 'set_background_color2'
        passed, rgb_color_ints, rgb_color_floats = self._load_color(
            title, color_edit, rgb_color_ints, func_name)
        if passed:
            self.background_color2_int = rgb_color_ints
            self.background_color2_float = rgb_color_floats

    def on_caero_color(self) -> None:
        """ Choose a CAERO color"""
        title = "Choose a default CAERO color"
        rgb_color_ints = self.caero_color_int
        color_edit = self.caero_color_edit
        func_name = 'set_caero_color'
        passed, rgb_color_ints, rgb_color_floats = self._load_nastran_color(
            title, color_edit, rgb_color_ints, func_name)
        if passed:
            self.caero_color_int = rgb_color_ints
            self.caero_color_float = rgb_color_floats

    def on_rbe_line_color(self) -> None:
        """ Choose an RBE color"""
        title = "Choose a default RBE2/RBE3 line color"
        rgb_color_ints = self.rbe_line_color_int
        color_edit = self.rbe_line_color_edit
        func_name = 'set_rbe_line_color'
        passed, rgb_color_ints, rgb_color_floats = self._load_nastran_color(
            title, color_edit, rgb_color_ints, func_name)
        if passed:
            self.rbe_line_color_int = rgb_color_ints
            self.rbe_line_color_float = rgb_color_floats

    def on_plotel_color(self) -> None:
        """ Choose an RBE color"""
        title = "Choose a default PLOTEL line color"
        rgb_color_ints = self.plotel_color_int
        color_edit = self.plotel_color_edit
        func_name = 'set_plotel_color'
        passed, rgb_color_ints, rgb_color_floats = self._load_nastran_color(
            title, color_edit, rgb_color_ints, func_name)
        if passed:
            self.plotel_color_int = rgb_color_ints
            self.plotel_color_float = rgb_color_floats

    def on_highlight_color(self) -> None:
        """ Choose a highlight color"""
        title = "Choose a highlight color"
        rgb_color_ints = self.highlight_color_int
        color_edit = self.highlight_color_edit
        func_name = 'set_highlight_color'
        passed, rgb_color_ints, rgb_color_floats = self._load_color(
            title, color_edit, rgb_color_ints, func_name)
        if passed:
            self.highlight_color_int = rgb_color_ints
            self.highlight_color_float = rgb_color_floats

    def on_highlight_opacity(self, value=None):
        if value is None:
            value = self.highlight_opacity_edit.value()
        self._highlight_opacity = value
        if self.win_parent is not None:
            self.settings.set_highlight_opacity(value)

    def on_highlight_point_size(self, value=None):
        if value is None:
            value = self.highlight_point_size_edit.value()
        self._highlight_point_size = value
        if self.win_parent is not None:
            self.settings.set_highlight_point_size(value)

    def on_corner_text_color(self) -> None:
        """Choose a corner text color"""
        rgb_color_ints = self.corner_text_color_int
        title = "Choose a corner text color"
        passed, rgb_color_ints, rgb_color_floats = self.on_color(
            self.corner_text_color_edit, rgb_color_ints, title)
        if passed:
            self.corner_text_color_int = rgb_color_ints
            self.corner_text_color_float = rgb_color_floats
            if self.win_parent is not None:
                self.settings.set_corner_text_color(rgb_color_floats)

    def on_default_corner_text_size(self):
        self.corner_text_size_edit.setValue(self._default_corner_text_size)
        self.on_corner_text_size(self._default_corner_text_size)

    def on_corner_text_size(self, value=None):
        if value is None:
            value = self.corner_text_size_edit.value()
        self._corner_text_size = value
        if self.win_parent is not None:
            self.settings.set_corner_text_size(value)

    def _load_nastran_color(self, title: str, color_edit: QPushButtonColor,
                            rgb_color_ints: tuple[int, int, int],
                            func_name: str):
        """helper method for nastran colors"""
        passed, rgb_color_ints, rgb_color_floats = self.on_color(
            color_edit, rgb_color_ints, title)
        if passed:
            if self.win_parent is not None:
                settings = self.settings.nastran_settings
                func_color = getattr(settings, func_name)
                func_color(rgb_color_floats)
        return passed, rgb_color_ints, rgb_color_floats

    def _load_color(self, title: str, color_edit: QPushButtonColor,
                    rgb_color_ints: tuple[int, int, int],
                    func_name: str):
        """helper method for ``on_background_color`` and ``on_background_color2``"""
        passed, rgb_color_ints, rgb_color_floats = self.on_color(
            color_edit, rgb_color_ints, title)
        if passed:
            if self.win_parent is not None:
                settings = self.settings
                func_color = getattr(settings, func_name)
                func_color(rgb_color_floats)
        return passed, rgb_color_ints, rgb_color_floats

    def on_color(self, color_edit: QPushButtonColor,
                 rgb_color_ints: tuple[int, int, int],
                 title: str) -> tuple[bool, tuple[int, int, int], Any]:
        """pops a color dialog"""
        col = QColorDialog.getColor(QtGui.QColor(*rgb_color_ints), self,
                                    title)
        if not col.isValid():
            return False, rgb_color_ints, None

        color_float: list[float] = col.getRgbF()[:3]
        color_int = tuple([int(colori * 255) for colori in color_float])

        assert isinstance(color_float[0], float), color_float
        assert isinstance(color_int[0], int), color_int

        set_pushbutton_color(color_edit, color_int)
        return True, color_int, color_float

    def on_picker_size(self) -> None:
        self._picker_size = float_locale(self.picker_size_edit.text())
        if self.win_parent is not None:
            self.win_parent.element_picker_size = self._picker_size / 100.
        #self.on_apply(force=True)

    def on_parallel_projection(self) -> None:
        """set the nastran properties preferences"""
        is_checked = self.parallel_projection_checkbox.isChecked()
        if self.win_parent is not None:
            self.settings.set_parallel_projection(is_checked)

    def on_trackball_camera(self) -> None:
        """set the nastran properties preferences"""
        is_checked = self.trackball_camera_checkbox.isChecked()
        if self.win_parent is not None:
            self.settings.set_trackball_camera(is_checked)

    def on_magnify(self, value=None) -> None:
        if value is None:
            value = self.magnify_edit.value()
        self._magnify = value
        if self.win_parent is not None:
            self.settings.set_magnify(value)

    #---------------------------------------------------------------------------
    def on_coord_scale(self, value=None) -> None:
        if value is None:
            value = self.coord_scale_edit.value()
        self._coord_scale = value
        if self.win_parent is not None:
            self.settings.set_coord_scale(value / 100.)

    def on_default_coord_scale(self) -> None:
        self.coord_scale_edit.setValue(self._default_coord_scale)
        self.on_coord_scale(self._default_coord_scale)

    def on_coord_text_scale(self, value=None) -> None:
        if value is None:
            value = self.coord_text_scale_edit.value()
        self._coord_text_scale = value
        if self.win_parent is not None:
            self.settings.set_coord_text_scale(value / 100.)

    def on_default_coord_text_scale(self):
        self.coord_text_scale_edit.setValue(self._default_coord_text_scale)
        self.on_coord_text_scale(self._default_coord_text_scale)

    #---------------------------------------------------------------------------

    def on_default_annotation_size(self) -> None:
        self.annotation_size_edit.setValue(self._default_annotation_size)
        self.on_annotation_size(self._default_annotation_size)

    def on_default_clipping_min(self) -> None:
        self.clipping_min_edit.setText(func_str(self._default_clipping_min))
        self.clipping_min_edit.setStyleSheet(QLINEEDIT_GOOD)

    def on_default_clipping_max(self) -> None:
        self.clipping_max_edit.setText(func_str(self._default_clipping_max))
        self.clipping_max_edit.setStyleSheet(QLINEEDIT_GOOD)

    def on_validate(self) -> bool:
        font_size_value, flag0 = check_float(self.font_size_edit)
        annotation_size_value, flag1 = check_float(self.annotation_size_edit)

        assert isinstance(self.annotation_color_float[0], float), self.annotation_color_float
        assert isinstance(self.annotation_color_int[0], int), self.annotation_color_int
        picker_size_value, flag2 = check_float(self.picker_size_edit)

        clipping_min_value, flag3 = check_float(self.clipping_min_edit)
        clipping_max_value, flag4 = check_float(self.clipping_max_edit)

        cart3d_fluent_include, flag5 = check_tuple_ints(self.cart3d_fluent_regions_include)
        cart3d_fluent_remove, flag6 = check_tuple_ints(self.cart3d_fluent_regions_remove)
        is_regions = len(cart3d_fluent_include) and len(cart3d_fluent_remove)
        if is_regions and ([flag5, flag6]):
            # error
            self.cart3d_fluent_regions_include.setStyleSheet(QLINEEDIT_ERROR)
            self.cart3d_fluent_regions_remove.setStyleSheet(QLINEEDIT_ERROR)

        if all([flag0, flag1, flag2, flag3, flag4, flag5, flag6]):
            self._annotation_size = annotation_size_value
            self._picker_size = picker_size_value

            self.out_data['font_size'] = int(font_size_value)
            self.out_data['min_clip'] = min(clipping_min_value, clipping_max_value)
            self.out_data['max_clip'] = max(clipping_min_value, clipping_max_value)
            self.out_data['cart3d_fluent_include'] = cart3d_fluent_include
            self.out_data['cart3d_fluent_remove'] = cart3d_fluent_remove
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_apply(self, force=False):
        passed = self.on_validate()

        if (passed or force) and self.win_parent is not None:
            self.settings.on_set_font_size(self.out_data['font_size'])
            self.settings.other_settings.update(self.out_data)
            #self.settings.set_annotation_size_color(self._annotation_size)
            #self.win_parent.element_picker_size = self._picker_size / 100.
        if passed and self.win_parent is not None:
            self.win_parent.clipping_obj.apply_clipping(self.out_data)
        return passed

    def on_ok(self):
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.out_data['close'] = True
        self.close()


def create_shear_moment_torque_edits(
    parent,
    opacity: float, point_size, line_width,
    color: ColorInt) -> tuple[QDoubleSpinBox, QDoubleSpinBox, QDoubleSpinBox, QPushButtonColor]:
    opacity_edit = create_opacity_edit(parent, opacity)
    opacity_edit.setToolTip('Sets the shear-moment-torque opacity (0=invisible, 1=solid)')

    point_size_edit = create_point_size_edit(parent, point_size)
    point_size_edit.setToolTip('Sets the shear-moment-torque node size')

    line_width_edit = create_line_width_edit(parent, line_width)
    line_width_edit.setToolTip('Sets the shear-moment-torque line width')

    color_edit = QPushButtonColor(color)
    color_edit.setToolTip('Sets the shear-moment-torque color')

    return opacity_edit, point_size_edit, line_width_edit, color_edit

def create_point_size_edit(parent, value: float) -> QDoubleSpinBox:
    point_size_edit = QDoubleSpinBox(parent)
    point_size_edit.setValue(value)
    point_size_edit.setRange(POINT_SIZE_MIN, POINT_SIZE_MAX)
    point_size_edit.setDecimals(2)
    point_size_edit.setSingleStep(0.5)
    return point_size_edit

def create_line_width_edit(parent, value: float) -> QDoubleSpinBox:
    line_width_edit = QDoubleSpinBox(parent)
    line_width_edit.setValue(value)
    #line_width_edit.setRange(POINT_SIZE_MIN, POINT_SIZE_MAX)
    #line_width_edit.setDecimals(2)
    #line_width_edit.setSingleStep(0.25)
    return line_width_edit

def create_opacity_edit(parent, value: float) -> QDoubleSpinBox:
    opacity_edit = QDoubleSpinBox(parent)
    opacity_edit.setValue(value)
    opacity_edit.setRange(OPACITY_MIN, OPACITY_MAX)
    opacity_edit.setDecimals(2)
    opacity_edit.setSingleStep(0.05)
    return opacity_edit


def set_pushbutton_color(color_edit: QPushButtonColor,
                         color_tuple: tuple[int, int, int]) -> None:
    color_edit.setStyleSheet(
        "QPushButton {"
        "background-color: rgb(%s, %s, %s);" % tuple(color_tuple) +
        #"border:1px solid rgb(255, 170, 255); "
        "}")

def check_float(cell: QDoubleSpinBox) -> tuple[str, bool]:
    text = cell.text()
    value = float_locale(text)
    return value, True

def check_label_float(cell) -> tuple[Optional[float], bool]:
    text = cell.text()
    try:
        value = eval_float_from_string(text)
        cell.setStyleSheet(QLINEEDIT_GOOD)
        return value, True
    except ValueError:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        return None, False

def check_tuple_ints(cell) -> tuple[tuple[int, ...], bool]:
    text = cell.text()
    text2 = text.replace(',', ' ').strip()
    sline = text2.split()
    try:
        values = [int(value) for value in sline]
        cell.setStyleSheet(QLINEEDIT_GOOD)
        return tuple(values), True
    except ValueError:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        return tuple(), False

def on_nastran(self: PreferencesWindow,
               result_name: str) -> None:
    """
    Auto-checks to verify result name is correct.
    Used for all Nastran settings, it's a lot less code.

    self = PreferencesWindow
    checkbox = self.nastran_displacement_checkbox
    result_name = displacement
    is_checked = True
    """
    #self.nastran_displacement_checkbox,
    checkbox: QCheckBox = getattr(self, f'nastran_{result_name}_checkbox')
    is_checked = checkbox.isChecked()
    if self.win_parent is not None:
        assert hasattr(self.nastran_settings, result_name), result_name
        setattr(self.nastran_settings, result_name, is_checked)


def main():  # pragma: no cover
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)


    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    #The Main window
    data = {
        'font_size' : 9,
        'use_startup_directory': True,
        'use_gradient_background' : True,
        'background_color' : (0., 0., 0.), # black
        'background_color2' : (1., 0., 1.), # purple
        'coord_scale' : 0.05,
        'coord_text_scale' : 1.0,
        'show_corner_coord' : False,
        'magnify' : 5,

        'corner_text_size' : 12,
        'corner_text_color' : (0., 1., 0.), # green

        'highlight_color' : (1., 1., 0.), # yellow
        'highlight_opacity' : 0.8,
        'highlight_point_size' : 10.0,

        'is_trackball_camera': True,
        'use_parallel_projection': True,
        'annotation_color' : (1., 0., 0.), # red
        'annotation_size' : 11,
        'picker_size' : 10.,

        'caero_color': (0.2, 0.7, 0.4),
        'rbe_line_color': (0.5, 0.6, 0.7),
        'plotel_color': (0.5, 0.6, 0.7),
        'nastran_version' : 'MSC',

        'min_clip' : 0.,
        'max_clip' : 10,

        'dim_max' : 502.,
        #------------------------------------
        #other
        'cart3d_fluent_include': (),
        'cart3d_fluent_remove': (3,),
        'units_model_in': ('in', 'lbf', 's', 'psi'),
        'units_length': 'in',
        #'units_area': 'in^2',
        'units_force': 'lbf',
        'units_moment': 'in-lbf',
        'units_pressure': 'psi',
        'units_stress': 'psi',
        'units_displacement': 'in',
        'units_velocity': 'in/s',
        'units_acceleration': 'in/s^2',
    }
    for name in NASTRAN_BOOL_KEYS:
        data[name] = True
    main_window = PreferencesWindow(data)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":  # pragma: no cover
    main()
