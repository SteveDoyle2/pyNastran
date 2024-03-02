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
    QTabWidget, QWidget,
)

from pyNastran.utils.locale import func_str, float_locale
from pyNastran.gui.utils.qt.pydialog import PyDialog, make_font, check_color
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
    HIGHLIGHT_COLOR, HIGHLIGHT_OPACITY, HIGHLIGHT_POINT_SIZE, # HIGHLIGHT_LINE_THICKNESS,
    SHEAR_MOMENT_TORQUE_COLOR, SHEAR_MOMENT_TORQUE_OPACITY, SHEAR_MOMENT_TORQUE_POINT_SIZE, SHEAR_MOMENT_TORQUE_LINE_THICKNESS,
    OPACITY_MIN, OPACITY_MAX,
    USE_PARALLEL_PROJECTION,
    NASTRAN_BOOL_KEYS,
    POINT_SIZE_MIN, POINT_SIZE_MAX,
    COORD_TEXT_SCALE_MIN, COORD_TEXT_SCALE_MAX,
    CORNER_TEXT_SIZE_MIN, CORNER_TEXT_SIZE_MAX,

    COORD_SCALE_MIN, COORD_SCALE_MAX,
    MAGNIFY_MIN, MAGNIFY_MAX,
    ANNOTATION_SIZE_MIN, ANNOTATION_SIZE_MAX,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui_objects.settings import Settings, NastranSettings


USE_TABS = True
IS_SMT = True
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

        self._nastran_is_element_quality = data['nastran_is_element_quality']
        self._nastran_is_properties = data['nastran_is_properties']
        self._nastran_is_3d_bars = data['nastran_is_3d_bars']
        self._nastran_is_3d_bars_update = data['nastran_is_3d_bars_update']
        self._nastran_is_bar_axes = data['nastran_is_bar_axes']
        self._nastran_create_coords = data['nastran_create_coords']
        self._nastran_is_shell_mcids = data['nastran_is_shell_mcids']

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

        self._nastran_spc_force = data['nastran_spc_force']
        self._nastran_mpc_force = data['nastran_mpc_force']
        self._nastran_applied_load = data['nastran_applied_load']

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
        self.startup_directory_label = QLabel('Remember Last Directory:')
        self.startup_directory_checkbox = QCheckBox(self)
        self.startup_directory_checkbox.setChecked(self.use_startup_directory)
        self.startup_directory_checkbox.setToolTip('True: Remember the last directory when saving\n'
                                                   'False: Start from local directory')
        #self.startup_directory_button = QPushButton('...')

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
        self.highlight_opacity_edit = QDoubleSpinBox(self)
        self.highlight_opacity_edit.setValue(self._highlight_opacity)
        self.highlight_opacity_edit.setRange(OPACITY_MIN, OPACITY_MAX)
        self.highlight_opacity_edit.setDecimals(2)
        self.highlight_opacity_edit.setSingleStep(0.05)
        self.highlight_opacity_edit.setToolTip('Sets the highlight opacity (0=invisible, 1=solid)')
        self.highlight_opacity_button = QPushButton("Default")

        self.highlight_point_size_label = QLabel("Highlight Point Size:")
        self.highlight_point_size_edit = QDoubleSpinBox(self)
        self.highlight_point_size_edit.setValue(self._highlight_point_size)
        self.highlight_point_size_edit.setRange(POINT_SIZE_MIN, POINT_SIZE_MAX)
        self.highlight_point_size_edit.setDecimals(2)
        self.highlight_point_size_edit.setSingleStep(0.25)
        self.highlight_point_size_edit.setToolTip('Sets the highlight node size')
        self.highlight_point_size_button = QPushButton("Default")

        # Text Color
        self.highlight_color_label = QLabel("Highlight Color:")
        self.highlight_color_edit = QPushButtonColor(self.highlight_color_int)
        self.highlight_color_edit.setToolTip('Sets the highlight color')

        #-----------------------------------------------------------------------
        # shear_moment_torque Color
        if IS_SMT:
            self._shear_moment_torque_opacity = 0.8
            self._shear_moment_torque_point_size = 10.0
            self.shear_moment_torque_color_int = (0, 0, 0)
            self._shear_moment_torque_line_width = 5.0

            self.shear_moment_torque_label = QLabel("Shear-Moment-Torque:")

            self.shear_moment_torque_opacity_label = QLabel("Opacity:")
            self.shear_moment_torque_opacity_edit = QDoubleSpinBox(self)
            self.shear_moment_torque_opacity_edit.setValue(self._shear_moment_torque_opacity)
            self.shear_moment_torque_opacity_edit.setRange(OPACITY_MIN, OPACITY_MAX)
            self.shear_moment_torque_opacity_edit.setDecimals(2)
            self.shear_moment_torque_opacity_edit.setSingleStep(0.05)
            self.shear_moment_torque_opacity_edit.setToolTip('Sets the shear-moment-torque opacity (0=invisible, 1=solid)')
            self.shear_moment_torque_opacity_button = QPushButton("Default")

            self.shear_moment_torque_point_size_label = QLabel("Point Size:")
            self.shear_moment_torque_point_size_edit = QDoubleSpinBox(self)
            self.shear_moment_torque_point_size_edit.setValue(self._shear_moment_torque_point_size)
            self.shear_moment_torque_point_size_edit.setRange(POINT_SIZE_MIN, POINT_SIZE_MAX)
            self.shear_moment_torque_point_size_edit.setDecimals(2)
            self.shear_moment_torque_point_size_edit.setSingleStep(0.25)
            self.shear_moment_torque_point_size_edit.setToolTip('Sets the shear-moment-torque node size')
            self.shear_moment_torque_point_size_button = QPushButton("Default")

            self.shear_moment_torque_line_width_label = QLabel("Line Width:")
            self.shear_moment_torque_line_width_edit = QDoubleSpinBox(self)
            self.shear_moment_torque_line_width_edit.setValue(self._shear_moment_torque_line_width)
            #self.shear_moment_torque_line_width_edit.setRange(POINT_SIZE_MIN, POINT_SIZE_MAX)
            #self.shear_moment_torque_line_width_edit.setDecimals(2)
            #self.shear_moment_torque_line_width_edit.setSingleStep(0.25)
            self.shear_moment_torque_line_width_edit.setToolTip('Sets the shear-moment-torque line width')
            self.shear_moment_torque_line_width_button = QPushButton("Default")

            # Text Color
            self.shear_moment_torque_color_label = QLabel("Color:")
            self.shear_moment_torque_color_edit = QPushButtonColor(self.shear_moment_torque_color_int)
            self.shear_moment_torque_color_edit.setToolTip('Sets the shear-moment-torque color')


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
        self.parallel_projection_edit = QCheckBox()
        self.parallel_projection_edit.setChecked(self._parallel_projection)
        self.parallel_projection_edit.setToolTip('Checked: Typical engineering perspective\n'
                                                 'Unchecked: Distort the model like a real camera')

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

        #-----------------------------------------------------------------------
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
        self.nastran_is_3d_bars_update_checkbox.setToolTip('Update the 3D Bar/Beam cross-sections when deformations are applied')
        self.nastran_is_3d_bars_update_checkbox.setChecked(self._nastran_is_3d_bars_update)

        self.nastran_is_shell_mcids_checkbox = QCheckBox('Shell MCIDs')
        self.nastran_is_shell_mcids_checkbox.setToolTip('Calculate the Material Coordinate Systems for Shells')
        self.nastran_is_shell_mcids_checkbox.setChecked(self._nastran_is_shell_mcids)

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
            self.nastran_displacement_checkbox.setChecked(self._nastran_displacement)
            self.nastran_velocity_checkbox.setChecked(self._nastran_velocity)
            self.nastran_acceleration_checkbox.setChecked(self._nastran_acceleration)
            self.nastran_eigenvector_checkbox.setChecked(self._nastran_eigenvector)

            self.nastran_spc_force_checkbox = QCheckBox('SPC Force')
            self.nastran_mpc_force_checkbox = QCheckBox('MPC Force')
            self.nastran_applied_load_checkbox = QCheckBox('Applied Load')
            self.nastran_spc_force_checkbox.setChecked(self._nastran_spc_force)
            self.nastran_mpc_force_checkbox.setChecked(self._nastran_mpc_force)
            self.nastran_applied_load_checkbox.setChecked(self._nastran_applied_load)

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
        # closing
        self.reset_defaults_button = QPushButton('Reset Defaults')
        self.apply_button = QPushButton('Apply')
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Cancel')

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

        grid.addWidget(self.startup_directory_label, irow, 0)
        grid.addWidget(self.startup_directory_checkbox, irow, 1)
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
        grid.addWidget(self.parallel_projection_edit, irow, 1)
        irow += 1

        #--------------------------------------------------

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
            vbox_nastran_tab.addStretch()

            nastran_tab_widget = QWidget(self)
            nastran_tab_widget.setLayout(vbox_nastran_tab)

            tabs.addTab(general_tab_widget, 'General')
            tabs.addTab(nastran_tab_widget, 'Nastran')
            vbox.addWidget(tabs)

        else:
            vbox = QVBoxLayout()
            vbox.addLayout(grid)
            vbox.addLayout(vbox_nastran)
            vbox.addLayout(vbox_nastran_results)
        #vbox.addStretch()
        #vbox.addLayout(grid2)
        vbox.addStretch()

        vbox.addWidget(self.reset_defaults_button)
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def _get_grid_nastran_layout(self) -> QGridLayout:
        grid_nastran = QGridLayout()
        irow = 0

        grid_nastran.addWidget(self.nastran_create_coords_checkbox, irow, 0)
        irow += 1

        grid_nastran.addWidget(self.nastran_is_element_quality_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_is_properties_checkbox, irow, 1)
        irow += 1

        grid_nastran.addWidget(self.nastran_is_bar_axes_checkbox, irow, 0)
        irow += 1

        grid_nastran.addWidget(self.nastran_is_shell_mcids_checkbox, irow, 0)
        irow += 1

        grid_nastran.addWidget(self.nastran_is_3d_bars_checkbox, irow, 0)
        grid_nastran.addWidget(self.nastran_is_3d_bars_update_checkbox, irow, 1)
        irow += 1
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

    def _set_nastran_connections(self):
        # format-specific
        self.nastran_is_element_quality_checkbox.clicked.connect(partial(on_nastran, self, 'is_element_quality'))
        self.nastran_is_properties_checkbox.clicked.connect(partial(on_nastran, self, 'is_properties'))
        self.nastran_is_3d_bars_checkbox.clicked.connect(partial(on_nastran, self, 'is_3d_bars'))
        self.nastran_is_3d_bars_update_checkbox.clicked.connect(partial(on_nastran, self, 'is_3d_bars_update'))
        self.nastran_is_bar_axes_checkbox.clicked.connect(partial(on_nastran, self, 'is_bar_axes'))
        self.nastran_create_coords_checkbox.clicked.connect(partial(on_nastran, self, 'create_coords'))
        self.nastran_is_shell_mcids_checkbox.clicked.connect(partial(on_nastran, self, 'is_shell_mcids'))

        #self.nastran_is_shell_mcid_checkbox.clicked.connect(self.on_nastran_is_shell_mcids)
        #self.nastran_is_shell_mcid_checkbox.clicked.connect(self.on_nastran_is_shell_mcids2)

        self.nastran_displacement_checkbox.clicked.connect(partial(on_nastran, self, 'displacement'))
        self.nastran_velocity_checkbox.clicked.connect(partial(on_nastran, self, 'acceleration'))
        self.nastran_acceleration_checkbox.clicked.connect(partial(on_nastran, self, 'acceleration'))
        self.nastran_eigenvector_checkbox.clicked.connect(partial(on_nastran, self, 'eigenvector'))

        self.nastran_spc_force_checkbox.clicked.connect(partial(on_nastran, self, 'spc_force'))
        self.nastran_mpc_force_checkbox.clicked.connect(partial(on_nastran, self, 'mpc_force'))
        self.nastran_applied_load_checkbox.clicked.connect(partial(on_nastran, self, 'applied_load'))
        self.nastran_grid_point_force_checkbox.clicked.connect(partial(on_nastran, self, 'grid_point_force'))

        self.nastran_force_checkbox.clicked.connect(partial(on_nastran, self, 'force'))
        self.nastran_strain_checkbox.clicked.connect(partial(on_nastran, self, 'strain'))
        self.nastran_stress_checkbox.clicked.connect(partial(on_nastran, self, 'stress'))
        self.nastran_strain_energy_checkbox.clicked.connect(partial(on_nastran, self, 'strain_energy'))

    def set_connections(self):
        """creates the actions for the menu"""
        self.font_size_edit.valueChanged.connect(self.on_font)

        self.startup_directory_checkbox.clicked.connect(self.on_startup_directory_checked)

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

        self.parallel_projection_edit.clicked.connect(self.on_parallel_projection)

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
        #------------------------------------


        self.reset_defaults_button.clicked.connect(self.on_reset_defaults)
        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)
        # closeEvent

    def on_startup_directory_checked(self):
        """opens a folder dialog"""
        is_checked = self.startup_directory_checkbox.isCheckable()
        if self.win_parent is not None:
            self.nastran_settings.use_startup_directory = is_checked

    def on_reset_defaults(self):
        """reset all of the preferences to their defaults"""
        self.on_font(FONT_SIZE)
        self.on_default_annotation_size()
        self.on_default_clipping_max()
        self.on_default_clipping_min()
        self.on_default_coord_scale()
        self.on_default_coord_text_scale()
        self.on_default_corner_text_size()
        self.startup_directory_checkbox.setChecked(True)

        self.magnify_edit.setValue(MAGNIFY)
        self.picker_size_edit.setValue(self._picker_size)

        self.highlight_opacity_edit.setValue(HIGHLIGHT_OPACITY)
        self.highlight_point_size_edit.setValue(HIGHLIGHT_POINT_SIZE)
        self.parallel_projection_edit.setChecked(USE_PARALLEL_PROJECTION)
        #self.highlight_line_thickness_edit.setValue(HIGHLIGHT_LINE_THICKNESS)

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

        set_label_color(self.corner_text_color_edit, self.corner_text_color_int)
        set_label_color(self.highlight_color_edit, self.highlight_color_int)
        set_label_color(self.background_color_edit, self.background_color1_int)
        set_label_color(self.background_color2_edit, self.background_color2_int)
        set_label_color(self.annotation_color_edit, self.annotation_color_int)
        #set_label_color(self.shear_moment_torque_color_edit, self.shear_moment_torque_color_int)

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

    @property
    def settings(self) -> Settings:
        return self.win_parent.settings

    @property
    def nastran_settings(self) -> NastranSettings:
        return self.settings.nastran_settings

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
        self.nastran_label.setFont(bold_font)
        self.nastran_results_label.setFont(bold_font)
        if IS_SMT:
            self.shear_moment_torque_label.setFont(bold_font)


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
        passed, rgb_color_ints, rgb_color_floats = self._background_color(
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
        passed, rgb_color_ints, rgb_color_floats = self._background_color(
            title, color_edit, rgb_color_ints, func_name)
        if passed:
            self.background_color2_int = rgb_color_ints
            self.background_color2_float = rgb_color_floats

    def on_highlight_color(self) -> None:
        """ Choose a highlight color"""
        title = "Choose a highlight color"
        rgb_color_ints = self.highlight_color_int
        color_edit = self.highlight_color_edit
        func_name = 'set_highlight_color'
        passed, rgb_color_ints, rgb_color_floats = self._background_color(
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

    def _background_color(self, title: str, color_edit: QPushButtonColor,
                          rgb_color_ints: tuple[int, int, int],
                          func_name: str):
        """helper method for ``on_background_color`` and ``on_background_color2``"""
        passed, rgb_color_ints, rgb_color_floats = self.on_color(
            color_edit, rgb_color_ints, title)
        if passed:
            if self.win_parent is not None:
                settings = self.settings
                func_background_color = getattr(settings, func_name)
                func_background_color(rgb_color_floats)
        return passed, rgb_color_ints, rgb_color_floats

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

        set_label_color(color_edit, color_int)
        return True, color_int, color_float

    def on_picker_size(self) -> None:
        self._picker_size = float_locale(self.picker_size_edit.text())
        if self.win_parent is not None:
            self.win_parent.element_picker_size = self._picker_size / 100.
        #self.on_apply(force=True)

    def on_parallel_projection(self) -> None:
        """set the nastran properties preferences"""
        is_checked = self.parallel_projection_edit.isChecked()
        if self.win_parent is not None:
            self.settings.set_parallel_projection(is_checked)

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

        if all([flag0, flag1, flag2, flag3, flag4]):
            self._annotation_size = annotation_size_value
            self._picker_size = picker_size_value

            self.out_data['font_size'] = int(font_size_value)
            self.out_data['min_clip'] = min(clipping_min_value, clipping_max_value)
            self.out_data['max_clip'] = max(clipping_min_value, clipping_max_value)
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_apply(self, force=False):
        passed = self.on_validate()

        if (passed or force) and self.win_parent is not None:
            self.settings.on_set_font_size(self.out_data['font_size'])
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


def set_label_color(color_edit: QPushButtonColor,
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

def on_nastran(self: PreferencesWindow, result_name: str) -> None:
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
        'font_size' : 8,
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

        'use_parallel_projection': True,
        'annotation_color' : (1., 0., 0.), # red
        'annotation_size' : 11,
        'picker_size' : 10.,

        'min_clip' : 0.,
        'max_clip' : 10,

        'dim_max' : 502.,

    }
    for name in NASTRAN_BOOL_KEYS:
        data[name] = True
    main_window = PreferencesWindow(data)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":  # pragma: no cover
    main()
