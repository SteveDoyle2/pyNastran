from __future__ import annotations
import os
import sys
import copy
import warnings
import traceback
from pathlib import Path
from functools import wraps
from typing import Optional, Any, TYPE_CHECKING

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

try:
    import json5 as json
except ModuleNotFoundError:
    warnings.warn('couldnt find json5, using json')
    import json

from pyNastran.f06.dev.flutter.utils import get_raw_json
JSON_FILENAME, USE_VTK, USE_TABS = get_raw_json(allow_vtk=False)
from pyNastran.f06.dev.flutter.trade_layout import TradeLayout
from pyNastran.f06.dev.flutter.utils_qt import create_grid_from_list

from qtpy import QtCore
from qtpy.compat import getopenfilename  # getsavefilename
# from qtpy.QtGui import QIcon, QPixmap
from qtpy.QtWidgets import (
    QLabel, QWidget,
    QApplication, QVBoxLayout, QComboBox,
    QHBoxLayout, QPushButton, QGridLayout,
    QAction,
    QCheckBox, QLineEdit,
    QListWidgetItem, QAbstractItemView,
    QListWidget, QSpinBox, QTabWidget,)

from cpylog import SimpleLogger
import pyNastran
from pyNastran.gui.utils.qt.pydialog import QFloatEdit, make_font
from pyNastran.gui.qt_files.named_dock_widget import NamedDockWidget
from pyNastran.gui.qt_files.loggable_gui import LoggableGui

from pyNastran.f06.dev.flutter.gui_flutter import (
    export_flutter_results, get_list_float_or_none, get_float_or_none,
    get_selected_items_flat)

from pyNastran.f06.dev.flutter.actions_builder import Actions, Action, build_menus
from pyNastran.f06.dev.flutter.preferences_object import FlutterPreferencesObject
from pyNastran.f06.dev.flutter.preferences import (
    FLUTTER_BBOX_TO_ANCHOR_DEFAULT, LEGEND_LOC_DEFAULT,
    FONT_SIZE_DEFAULT, FLUTTER_NCOLUMNS_DEFAULT, FREQ_NDIGITS_DEFAULT, FREQ_DIVERGENCE_TOL)

from pyNastran.f06.flutter_response import Limit  # FlutterResponse
from pyNastran.f06.parse_flutter import get_flutter_units

from pyNastran.f06.dev.flutter.utils_qt import (
    load_lineedits, load_pulldowns, load_min_max_lineedits)
from pyNastran.f06.dev.flutter.utils import (
    validate_json,
    get_point_removal_str, get_noline_nopoints,
    point_removal_str_to_point_removal,
    _float_passed_to_default, get_plot_flags,
    update_ylog_style, get_png_filename,
    load_f06_op2, get_vlines, get_damping_crossings,
    X_PLOT_TYPES, PLOT_TYPES, UNITS_IN, UNITS_OUT,
    MODE_SWITCH_METHODS)

PKG_PATH = Path(pyNastran.__path__[0])

AERO_PATH = PKG_PATH / '..' / 'models' / 'aero'

from pyNastran.f06.dev.flutter.vtk_data import VtkData
from pyNastran.f06.parse_flutter import FlutterResponse

if TYPE_CHECKING:
    from pyNastran.f06.dev.flutter.gui_flutter_plot import FlutterGui

QLINEEDIT_WHITE = 'QLineEdit {background-color: white;}'
QLINEEDIT_RED = 'QLineEdit {background-color: red;}'
ICON_PATH = Path('')


def _set_modes_table(modes_widget: QListWidget,
                     modes: list[int],
                     freqs: list[float]) -> None:
    modes_widget.clear()
    for imode, freq in zip(modes, freqs):
        mode = QListWidgetItem(f'Mode {imode}; f={freq:.2f}')
        # mode.itemClicked.connect(self._on_update_mode)
        mode.setSelected(True)
        modes_widget.addItem(mode)


class PlotLayout(QHBoxLayout):
    def __init__(self, parent, ifile: int, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parent: FlutterGui = parent
        self.log = parent.log
        self.ifile = ifile
        self.responses = {}
        self.modes = []
        # self.log = parent.log
        #--------------------------------
        self.f06_filename_label = QLabel('F06 Filename:', parent)
        self.f06_filename_edit = QLineEdit(parent)
        self.f06_filename_browse = QPushButton('Browse...', parent)

        self.bdf_filename_checkbox = QCheckBox('BDF Filename:', parent)
        self.bdf_filename_edit = QLineEdit(parent)
        self.bdf_filename_browse = QPushButton('Browse...', parent)
        self.bdf_filename_checkbox.setChecked(False)
        self.bdf_filename_edit.setEnabled(False)
        self.bdf_filename_browse.setEnabled(False)
        self.bdf_filename_edit.setToolTip('Loads the Nastran Geometry')

        self.op2_filename_checkbox = QCheckBox('OP2 Filename:', parent)
        self.op2_filename_edit = QLineEdit(parent)
        self.op2_filename_browse = QPushButton('Browse...', parent)
        self.op2_filename_checkbox.setChecked(False)
        self.op2_filename_edit.setEnabled(False)
        self.op2_filename_browse.setEnabled(False)
        self.op2_filename_edit.setToolTip('Loads the Nastran Results (and geometry if BDF Filename is empty)')

        self.use_rhoref_checkbox = QCheckBox('Sea Level Rho Ref', parent)
        self.use_rhoref_checkbox.setChecked(False)

        self.log_xscale_checkbox = QCheckBox('Log Scale x', parent)
        self.log_yscale1_checkbox = QCheckBox('Log Scale y1', parent)
        self.log_yscale2_checkbox = QCheckBox('Log Scale y2', parent)
        self.log_xscale_checkbox.setChecked(False)
        self.log_yscale1_checkbox.setChecked(False)
        self.log_yscale2_checkbox.setChecked(False)

        self.show_points_checkbox = QCheckBox('Show Points', parent)
        self.show_mode_number_checkbox = QCheckBox('Show Mode Number', parent)
        self.show_detailed_mode_info_checkbox = QCheckBox('Show Detailed Mode Info', parent)
        self.point_spacing_label = QLabel('Point Spacing', parent)
        self.point_spacing_spinner = QSpinBox(parent)
        self.include_rigid_body_modes_checkbox = QCheckBox('Include Rigid Body Modes', parent)
        self.number_rigid_body_modes_label = QLabel('nRigid Body Modes', parent)
        self.number_rigid_body_modes_spinner = QSpinBox(parent)

        for obj in [
            self.number_rigid_body_modes_label,
            self.include_rigid_body_modes_checkbox,
            self.number_rigid_body_modes_spinner]:
            obj.setVisible(False)

        self.show_lines_checkbox = QCheckBox('Show Lines', parent)
        self.show_points_checkbox.setChecked(True)
        self.show_lines_checkbox.setChecked(True)
        self.show_points_checkbox.setToolTip('The points are symbols')
        self.show_mode_number_checkbox.setToolTip('The points are the mode number')
        self.show_detailed_mode_info_checkbox.setToolTip('Lists the 0% eas/freq range')
        self.point_spacing_spinner.setToolTip('Skip Every Nth Point; 0=Plot All')
        self.point_spacing_spinner.setValue(0)
        self.point_spacing_spinner.setMinimum(0)
        self.point_spacing_spinner.setMaximum(30)

        self.index_lim_label = QLabel('Index Limits:', parent)
        self.index_lim_edit_min = QFloatEdit('0', parent)
        self.index_lim_edit_max = QFloatEdit(parent)

        self.eas_lim_label = QLabel('EAS Limits:', parent)
        self.eas_lim_edit_min = QFloatEdit('0', parent)
        self.eas_lim_edit_max = QFloatEdit(parent)

        self.tas_lim_label = QLabel('TAS Limits:', parent)
        self.tas_lim_edit_min = QFloatEdit('0', parent)
        self.tas_lim_edit_max = QFloatEdit(parent)

        self.mach_lim_label = QLabel('Mach Limits:', parent)
        self.mach_lim_edit_min = QFloatEdit(parent)
        self.mach_lim_edit_max = QFloatEdit(parent)

        self.alt_lim_label = QLabel('Alt Limits:', parent)
        self.alt_lim_edit_min = QFloatEdit(parent)
        self.alt_lim_edit_max = QFloatEdit(parent)

        self.q_lim_label = QLabel('Q Limits:', parent)
        self.q_lim_edit_min = QFloatEdit(parent)
        self.q_lim_edit_max = QFloatEdit(parent)

        self.rho_lim_label = QLabel('Rho Limits:', parent)
        self.rho_lim_edit_min = QFloatEdit('0', parent)
        self.rho_lim_edit_max = QFloatEdit(parent)

        self.damp_lim_label = QLabel('Damping Limits (g):', parent)
        self.damp_lim_edit_min = QFloatEdit('-0.3', parent)
        self.damp_lim_edit_max = QFloatEdit('0.3', parent)

        self.freq_lim_label = QLabel('Freq Limits (Hz):', parent)
        self.freq_lim_edit_min = QFloatEdit('0', parent)
        self.freq_lim_edit_max = QFloatEdit(parent)

        self.kfreq_lim_label = QLabel('KFreq Limits:', parent)
        self.kfreq_lim_edit_min = QFloatEdit(parent)
        self.kfreq_lim_edit_max = QFloatEdit(parent)

        self.ikfreq_lim_label = QLabel('1/KFreq Limits:', parent)
        self.ikfreq_lim_edit_min = QFloatEdit(parent)
        self.ikfreq_lim_edit_max = QFloatEdit(parent)

        self.eigr_lim_label = QLabel('Real Eigenvalue:', parent)
        self.eigr_lim_edit_min = QFloatEdit(parent)
        self.eigr_lim_edit_max = QFloatEdit(parent)

        self.eigi_lim_label = QLabel('Imag Eigenvalue:', parent)
        self.eigi_lim_edit_min = QFloatEdit(parent)
        self.eigi_lim_edit_max = QFloatEdit(parent)

        # --------------------------------------------
        self.freq_tol_label = QLabel('dFreq Tol (Hz) Dash:', parent)
        self.freq_tol_edit = QFloatEdit('-1.0', parent)
        self.freq_tol_edit.setToolTip("Applies a dotted line for modes that don't change by more than some amount")

        self.freq_tol_remove_label = QLabel('dFreq Tol (Hz) Remove:', parent)
        self.freq_tol_remove_edit = QFloatEdit('-1.0', parent)
        self.freq_tol_remove_edit.setToolTip('Removes a mode if it meets dFreq Tol (Hz) Dash and Remove')

        self.mag_tol_label = QLabel('Magnitude Tol:', parent)
        self.mag_tol_edit = QFloatEdit('-1.0', parent)
        self.mag_tol_edit.setToolTip('Filters modal participation factors based on magnitude')

        self.subcase_label = QLabel('Subcase:', parent)
        self.subcase_edit = QComboBox(parent)

        units_msg = (
            "english_in: inch/s, slich/in^3\n"
            "english_ft: ft/s,   slug/ft^3\n"
            "english_kt: knots,  slug/ft^3\n"
            "si:         m/s,    kg/m^3\n"
            "si-mm:      mm/s,   Mg/mm^3\n"
        )
        self.x_plot_type_label = QLabel('X-Axis Plot Type:', parent)
        self.x_plot_type_pulldown = QComboBox(parent)
        self.x_plot_type_pulldown.addItems(X_PLOT_TYPES)
        self.x_plot_type_pulldown.setToolTip('sets the x-axis')

        self.plot_type_label = QLabel('Plot Type:', parent)
        self.plot_type_pulldown = QComboBox(parent)
        self.plot_type_pulldown.addItems(PLOT_TYPES)
        # self.plot_type_pulldown.setToolTip(units_msg)

        self.units_in_label = QLabel('Units In:', parent)
        self.units_in_pulldown = QComboBox(parent)
        self.units_in_pulldown.addItems(UNITS_IN)
        self.units_in_pulldown.setToolTip(units_msg)
        iunits_in = UNITS_IN.index('english_in')
        self.units_in_pulldown.setCurrentIndex(iunits_in)
        self.units_in_pulldown.setToolTip('Sets the units for the F06/OP2; set when loaded')

        self.units_out_label = QLabel('Units Out:', parent)
        self.units_out_pulldown = QComboBox(parent)
        self.units_out_pulldown.addItems(UNITS_OUT)
        self.units_out_pulldown.setToolTip(units_msg)
        iunits_out = UNITS_IN.index('english_kt')
        self.units_out_pulldown.setCurrentIndex(iunits_out)
        self.units_out_pulldown.setToolTip('Sets the units for the plot; may be updated')

        self.output_directory_label = QLabel('Output Directory:', parent)
        self.output_directory_edit = QLineEdit('', parent)
        self.output_directory_browse = QPushButton('Browse...', parent)
        self.output_directory_edit.setDisabled(True)
        self.output_directory_browse.setDisabled(True)

        self.vl_label = QLabel('VL, Limit:', parent)
        self.vl_edit = QFloatEdit('', parent)
        self.vl_edit.setToolTip('Makes a vertical line for VL')

        self.vf_label = QLabel('VF, Flutter:', parent)
        self.vf_edit = QFloatEdit('', parent)
        self.vf_edit.setToolTip('Makes a vertical line for VF')

        self.damping_required_label = QLabel('Damping Required, g:', parent)
        self.damping_required_edit = QFloatEdit('', parent)
        self.damping_required_edit.setToolTip('Enables the flutter crossing (e.g., 0.0 for 0%)')

        self.damping_required_tol_label = QLabel('Damping Required Tol, g:', parent)
        self.damping_required_tol_edit = QFloatEdit('', parent)
        self.damping_required_tol_edit.setToolTip('Tolerance for Damping Required. The crossing will be reported at the required value')

        self.damping_label = QLabel('Damping, g:', parent)
        self.damping_edit = QFloatEdit('', parent)
        self.damping_edit.setToolTip('Enables the flutter crossing (e.g., 0.03 for 3%)')

        self.eas_flutter_range_label = QLabel('EAS Flutter/Diverg Range:', parent)
        self.eas_flutter_range_edit_min = QFloatEdit('', parent)
        self.eas_flutter_range_edit_max = QFloatEdit('', parent)
        self.eas_flutter_range_edit_min.setToolTip('Defines the flutter/divergence crossing range')
        self.eas_flutter_range_edit_max.setToolTip('Defines the flutter/divergence crossing range')

        # self.eas_diverg_range_label = QLabel('EAS Diverg Range:', parent)
        # self.eas_diverg_range_edit_min = QFloatEdit('', parent)
        # self.eas_diverg_range_edit_max = QFloatEdit('', parent)
        # self.eas_diverg_range_edit_min.setToolTip('Defines the divergence crossing range')
        # self.eas_diverg_range_edit_max.setToolTip('Defines the divergence crossing range')

        self.point_removal_label = QLabel('Point Removal:', parent)
        self.point_removal_edit = QLineEdit('', parent)
        self.point_removal_edit.setToolTip('Remove bad points from a mode; "400:410,450:500"')

        self.mode_label = QLabel('Mode:', parent)
        self.mode_edit = QSpinBox(parent)
        self.mode_edit.setMinimum(1)
        # self.mode_edit.SetValue(3)
        self.mode_edit.setToolTip('Sets the mode')

        self.velocity_label = QLabel('Velocity Point:', parent)
        self.velocity_edit = QComboBox(parent)
        self.velocity_edit.setToolTip('Sets the velocity (input units)')

        self.f06_load_button = QPushButton('Load F06', parent)
        self.run_button = QPushButton('Run', parent)

        self.pop_vtk_gui_button = QPushButton('Open GUI', parent)
        self.solution_type_label = QLabel('Solution Type:', parent)
        self.solution_type_pulldown = QComboBox(parent)
        self.mode2_label = QLabel('Mode:', parent)
        self.mode2_pulldown = QComboBox(parent)

        self.mode_switch_method_label = QLabel('Mode Switch Method:', parent)
        self.mode_switch_method_pulldown = QComboBox(parent)
        self.mode_switch_method_pulldown.addItems(MODE_SWITCH_METHODS)

        self.modes_widget = QListWidget(self.parent)
        self.modes_widget.setMaximumWidth(200)  # was 100 when no freq
        self.modes_widget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        _set_modes_table(self.modes_widget, [0], [0.])

        self.on_plot_type()
        self.on_enable_bdf()
        self.on_enable_op2()
        self.on_hide_vtk()

        #--------------------------------
        self.setup_layout()
        self.setup_connections()

    def setup_layout(self) -> None:
        file_row = 0
        parent = self.parent
        hbox = QGridLayout(parent)
        hbox.addWidget(self.f06_filename_label, file_row, 0)
        hbox.addWidget(self.f06_filename_edit, file_row, 1)
        hbox.addWidget(self.f06_filename_browse, file_row, 2)
        file_row += 1

        grid = create_grid_from_list(parent, [
            (self.units_in_label, self.units_in_pulldown, self.use_rhoref_checkbox),
            (self.units_out_label, self.units_out_pulldown),
            (self.subcase_label, self.subcase_edit),
            (self.x_plot_type_label, self.x_plot_type_pulldown),
            (self.plot_type_label, self.plot_type_pulldown),
            # x-axis
            (self.index_lim_label, self.index_lim_edit_min, self.index_lim_edit_max),
            (self.eas_lim_label, self.eas_lim_edit_min, self.eas_lim_edit_max),
            (self.tas_lim_label, self.tas_lim_edit_min, self.tas_lim_edit_max),
            (self.mach_lim_label, self.mach_lim_edit_min, self.mach_lim_edit_max),
            (self.alt_lim_label, self.alt_lim_edit_min, self.alt_lim_edit_max),
            (self.q_lim_label, self.q_lim_edit_min, self.q_lim_edit_max),
            (self.rho_lim_label, self.rho_lim_edit_min, self.rho_lim_edit_max),
            (self.kfreq_lim_label, self.kfreq_lim_edit_min, self.kfreq_lim_edit_max),
            (self.ikfreq_lim_label, self.ikfreq_lim_edit_min, self.ikfreq_lim_edit_max),
            # y-axes
            (self.damp_lim_label, self.damp_lim_edit_min, self.damp_lim_edit_max),
            (self.freq_lim_label, self.freq_lim_edit_min, self.freq_lim_edit_max),
            (self.eigr_lim_label, self.eigr_lim_edit_min, self.eigr_lim_edit_max),
            (self.eigi_lim_label, self.eigi_lim_edit_min, self.eigi_lim_edit_max),
            (self.freq_tol_label, self.freq_tol_edit),
            (self.freq_tol_remove_label, self.freq_tol_remove_edit),
            (self.mode_label, self.mode_edit),
            (self.velocity_label, self.velocity_edit),
            (self.mag_tol_label, self.mag_tol_edit),
            (self.output_directory_label, self.output_directory_edit, self.output_directory_browse),

            (self.vl_label, self.vl_edit),
            (self.vf_label, self.vf_edit),
            (self.damping_required_label, self.damping_required_edit),
            (self.damping_required_tol_label, self.damping_required_tol_edit),
            (self.damping_label, self.damping_edit),
            (self.eas_flutter_range_label, self.eas_flutter_range_edit_min, self.eas_flutter_range_edit_max),
            # (self.eas_diverg_range_label, self.eas_diverg_range_edit_min, self.eas_diverg_range_edit_max),
            (self.point_removal_label, self.point_removal_edit),
            (self.mode_switch_method_label, self.mode_switch_method_pulldown),
        ])
        self.output_directory_label.setDisabled(True)
        self.output_directory_edit.setDisabled(True)
        self.output_directory_browse.setDisabled(True)

        self.output_directory_label.setVisible(False)
        self.output_directory_edit.setVisible(False)
        self.output_directory_browse.setVisible(False)

        grid_check = create_grid_from_list(parent, [
            (self.log_xscale_checkbox, self.log_yscale1_checkbox, self.log_yscale2_checkbox),
            (self.show_points_checkbox, self.show_mode_number_checkbox, self.show_detailed_mode_info_checkbox),
            (self.point_spacing_label, self.point_spacing_spinner),
            (self.include_rigid_body_modes_checkbox, self.number_rigid_body_modes_label, self.number_rigid_body_modes_spinner),
            (self.show_lines_checkbox,),
        ])

        ok_cancel_hbox = QHBoxLayout()
        ok_cancel_hbox.addWidget(self.f06_load_button)
        ok_cancel_hbox.addWidget(self.run_button)

        hbox_check = QHBoxLayout()
        hbox_check.addLayout(grid_check)
        hbox_check.addStretch(1)

        grid_modes = self._grid_modes()

        vbox = QVBoxLayout()
        vbox.addLayout(hbox)
        vbox.addLayout(grid)
        vbox.addLayout(hbox_check)
        vbox.addStretch(1)
        vbox.addLayout(ok_cancel_hbox)
        vbox.addWidget(self.pop_vtk_gui_button)
        vbox.addLayout(grid_modes)
        # log_widget = ApplicationLogWidget(parent)

        # vbox2 = QHBoxLayout()
        if parent.use_dock_widgets:
            # self.modes_dock_widget = NamedDockWidget('Modes', self.modes_widget, self.parent)
            # self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.modes_dock_widget)
            # vbox2 = vbox
            self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock_widget)
            self.addWidget(self.modes_widget)
        else:
            # self.log_dock_widget.hide()
            self.addWidget(self.modes_widget)
        self.addLayout(vbox)

    def _grid_modes(self) -> QGridLayout:
        irow = 0
        grid_modes = QGridLayout()
        grid_modes.addWidget(self.solution_type_label, irow, 0)
        grid_modes.addWidget(self.solution_type_pulldown, irow, 1)
        self.solution_type_pulldown.addItems(['Real Modes', 'Complex Modes'])
        irow += 1

        grid_modes.addWidget(self.mode2_label, irow, 0)
        grid_modes.addWidget(self.mode2_pulldown, irow, 1)
        # self.solution_type_pulldown.addItems(['Real Modes', 'Complex Modes'])
        irow += 1
        return grid_modes

    def setup_connections(self) -> None:
        self.f06_load_button.clicked.connect(self.on_load_f06)
        # self.bdf_load_button.clicked.connect(self.on_load_bdf)
        # self.op2_load_button.clicked.connect(self.on_load_op2)

        self.x_plot_type_pulldown.currentIndexChanged.connect(self.on_plot_type)
        self.plot_type_pulldown.currentIndexChanged.connect(self.on_plot_type)
        self.subcase_edit.currentIndexChanged.connect(self.on_subcase)
        self.bdf_filename_checkbox.stateChanged.connect(self.on_enable_bdf)
        self.op2_filename_checkbox.stateChanged.connect(self.on_enable_op2)
        self.f06_filename_browse.clicked.connect(self.on_browse_f06)
        # self.modes_widget.itemSelectionChanged.connect(self.on_modes)
        # self.modes_widget.itemClicked.connect(self.on_modes)
        # self.modes_widget.currentRowChanged.connect(self.on_modes)
        self.run_button.clicked.connect(self.on_run)
        self.units_out_pulldown.currentIndexChanged.connect(self.on_units_out)
        # self.include_rigid_body_modes_checkbox.clicked.connect(self.on_rigid_body_modes)

    def on_units_out(self):
        units_out = self.units_out_pulldown.currentText()
        units_out_dict = get_flutter_units(units_out)

        tas_units = units_out_dict['velocity']
        eas_units = units_out_dict['eas']
        alt_units = units_out_dict['altitude']
        q_units = units_out_dict['dynamic_pressure']
        rho_units = units_out_dict['density']

        self.tas_lim_label.setText(f'TAS Limits ({tas_units}):')
        self.eas_lim_label.setText(f'EAS Limits ({eas_units}):')
        self.alt_lim_label.setText(f'Alt Limits ({alt_units}):')
        self.rho_lim_label.setText(f'Rho Limits ({rho_units}):')
        self.q_lim_label.setText(f'Q Limits ({q_units}):')
        self.vl_label.setText(f'VL, Limit ({eas_units}):')
        self.vf_label.setText(f'VF, Flutter ({eas_units}):')

    def on_hide_vtk(self) -> None:
        if USE_VTK:
            return
        objs = [
            self.bdf_filename_checkbox, self.bdf_filename_edit, self.bdf_filename_browse,
            self.op2_filename_checkbox, self.op2_filename_edit, self.op2_filename_browse,
            self.pop_vtk_gui_button, self.solution_type_label, self.solution_type_pulldown,
            self.mode2_label, self.mode2_pulldown,
        ]
        for obj in objs:
            obj.setVisible(False)

    def on_plot_type(self) -> None:
        x_plot_type = self.x_plot_type_pulldown.currentText()
        plot_type = self.plot_type_pulldown.currentText()
        self.on_units_out()

        flags = get_plot_flags(plot_type, x_plot_type)
        show_index_lim = flags['show_index_lim']
        show_eas_lim = flags['show_eas_lim']
        show_tas_lim = flags['show_tas_lim']
        show_mach_lim = flags['show_mach_lim']
        show_alt_lim = flags['show_alt_lim']
        show_q_lim = flags['show_q_lim']
        show_rho_lim = flags['show_rho_lim']

        show_xlim = flags['show_xlim']
        show_freq = flags['show_freq']
        show_damp = flags['show_damp']
        show_kfreq = flags['show_kfreq']
        show_ikfreq = flags['show_ikfreq']
        show_root_locus = flags['show_root_locus']
        show_zimmerman = flags['show_zimmerman']
        show_modal_participation = flags['show_modal_participation']

        # print(f'x_plot_type={x_plot_type} show_damp={show_damp}; show_xlim={show_xlim}')
        # assert show_xlim is False, show_xlim

        show_eigenvalue = show_root_locus or show_modal_participation
        show_xaxis = not show_eigenvalue
        show_freq_tol = show_xaxis or show_root_locus
        show_crossing = show_xaxis
        show_mode_switch = show_xaxis or show_root_locus

        self.mode_switch_method_label.setVisible(show_mode_switch)
        self.mode_switch_method_pulldown.setVisible(show_mode_switch)

        self.x_plot_type_label.setVisible(show_xaxis)
        self.x_plot_type_pulldown.setVisible(show_xaxis)
        self.freq_tol_label.setVisible(show_freq_tol)
        self.freq_tol_edit.setVisible(show_freq_tol)
        self.freq_tol_remove_edit.setVisible(show_freq_tol)

        self.index_lim_label.setVisible(show_index_lim)
        self.index_lim_edit_min.setVisible(show_index_lim)
        self.index_lim_edit_max.setVisible(show_index_lim)

        self.eas_lim_label.setVisible(show_eas_lim)
        self.eas_lim_edit_min.setVisible(show_eas_lim)
        self.eas_lim_edit_max.setVisible(show_eas_lim)

        self.tas_lim_label.setVisible(show_tas_lim)
        self.tas_lim_edit_min.setVisible(show_tas_lim)
        self.tas_lim_edit_max.setVisible(show_tas_lim)

        self.mach_lim_label.setVisible(show_mach_lim)
        self.mach_lim_edit_min.setVisible(show_mach_lim)
        self.mach_lim_edit_max.setVisible(show_mach_lim)

        self.alt_lim_label.setVisible(show_alt_lim)
        self.alt_lim_edit_min.setVisible(show_alt_lim)
        self.alt_lim_edit_max.setVisible(show_alt_lim)

        self.q_lim_label.setVisible(show_q_lim)
        self.q_lim_edit_min.setVisible(show_q_lim)
        self.q_lim_edit_max.setVisible(show_q_lim)

        self.rho_lim_label.setVisible(show_rho_lim)
        self.rho_lim_edit_min.setVisible(show_rho_lim)
        self.rho_lim_edit_max.setVisible(show_rho_lim)

        self.damp_lim_label.setVisible(show_damp)
        self.damp_lim_edit_min.setVisible(show_damp)
        self.damp_lim_edit_max.setVisible(show_damp)

        self.damping_required_label.setVisible(show_crossing)
        self.damping_required_edit.setVisible(show_crossing)
        self.damping_required_tol_label.setVisible(show_crossing)
        self.damping_required_tol_edit.setVisible(show_crossing)

        self.damping_label.setVisible(show_crossing)
        self.damping_edit.setVisible(show_crossing)

        self.freq_lim_label.setVisible(show_freq)
        self.freq_lim_edit_min.setVisible(show_freq)
        self.freq_lim_edit_max.setVisible(show_freq)

        self.kfreq_lim_label.setVisible(show_kfreq)
        self.kfreq_lim_edit_min.setVisible(show_kfreq)
        self.kfreq_lim_edit_max.setVisible(show_kfreq)

        self.ikfreq_lim_label.setVisible(show_ikfreq)
        self.ikfreq_lim_edit_min.setVisible(show_ikfreq)
        self.ikfreq_lim_edit_max.setVisible(show_ikfreq)

        self.eigr_lim_label.setVisible(show_root_locus)
        self.eigr_lim_edit_min.setVisible(show_root_locus)
        self.eigr_lim_edit_max.setVisible(show_root_locus)

        self.eigi_lim_label.setVisible(show_root_locus)
        self.eigi_lim_edit_min.setVisible(show_root_locus)
        self.eigi_lim_edit_max.setVisible(show_root_locus)

        self.vl_label.setVisible(show_eas_lim)
        self.vl_edit.setVisible(show_eas_lim)
        self.vf_label.setVisible(show_eas_lim)
        self.vf_edit.setVisible(show_eas_lim)

        show_items = [
            (show_modal_participation, (
                self.mode_label, self.mode_edit,
                self.velocity_label, self.velocity_edit,
                self.mag_tol_label, self.mag_tol_edit)),
            (not show_modal_participation, (
                self.point_spacing_label, self.point_spacing_spinner,
                self.show_mode_number_checkbox,
                self.show_detailed_mode_info_checkbox,
                self.show_lines_checkbox,
                self.log_xscale_checkbox,
                self.log_yscale1_checkbox, self.log_yscale2_checkbox,
            ),),
        ]
        for show_hide, items in show_items:
            for item in items:
                item.setVisible(show_hide)

    def on_enable_bdf(self) -> None:
        state = self.bdf_filename_checkbox.isChecked()
        self.bdf_filename_edit.setEnabled(state)
        self.bdf_filename_browse.setEnabled(state)
        if state:
            self.bdf_filename_edit.setText(self.parent.bdf_filename)
        else:
            self.bdf_filename_edit.setText(self.parent._bdf_filename_default)

    def on_enable_op2(self) -> None:
        state = self.op2_filename_checkbox.isChecked()
        self.op2_filename_edit.setEnabled(state)
        self.op2_filename_browse.setEnabled(state)
        if state:
            self.op2_filename_edit.setText(self.parent.op2_filename)
        else:
            self.op2_filename_edit.setText(self.parent._op2_filename_default)

    def on_browse_f06(self) -> None:
        """pops a dialog to select the f06 file"""
        title = 'Load a Flutter (Nastran F06, Zona Out) File'
        qt_wildcard = 'F06 File (*.f06);; Zona File (*.out)'
        basedir = os.path.dirname(self.parent.f06_filename)
        if not os.path.exists(basedir):
            basedir = '.'
        fname, wildcard_level = getopenfilename(
            self, caption=title, basedir=basedir, filters=qt_wildcard,)
        if fname == '':
            return
        self.f06_filename_edit.setText(fname)
        self.run_button.setEnabled(False)

    # @dont_crash
    def on_load_f06(self, event) -> None:
        parent = self.parent
        log = parent.log
        f06_filename = os.path.abspath(self.f06_filename_edit.text())
        if not os.path.exists(f06_filename) or not os.path.isfile(f06_filename):
            self.f06_filename_edit.setStyleSheet(QLINEEDIT_RED)
            self.log.error(f"can't find {f06_filename}")
            return
        self.f06_filename_edit.setStyleSheet(QLINEEDIT_WHITE)
        f06_units = self.units_in_pulldown.currentText()
        out_units = self.units_out_pulldown.currentText()

        parent.use_rhoref = self.use_rhoref_checkbox.isChecked()
        model, self.responses = load_f06_op2(
            f06_filename, log,
            f06_units, out_units,
            parent.use_rhoref, stop_on_failure=False)

        subcases = list(self.responses.keys())
        if len(subcases) == 0:
            log.error('No subcases found')
            return
        # self.log.info(f'on_load_f06: subcases={subcases}')
        parent.f06_filename = f06_filename
        parent._units_in = f06_units
        parent._units_out = out_units
        parent.add_recent_file(f06_filename)
        self.update_subcases(subcases)
        self.run_button.setEnabled(True)




    def on_subcase(self) -> None:
        subcase, is_subcase_valid = self._get_subcase()
        if not is_subcase_valid:
            return
        response: FlutterResponse = self.responses[subcase]
        # self.log.info(f'on_subcase; response.results.shape={response.results.shape}')
        freqs = response.results[:, 0, response.ifreq].ravel()
        self._update_modal_participation_velocity(response)
        # self.log.info(f'on_subcase; freqs={freqs}')
        self.update_modes_table(response.modes, freqs)

    def _update_modal_participation_velocity(self, response: FlutterResponse) -> None:
        nvelocity = len(response.eigr_eigi_velocity)
        self.velocity_edit.clear()
        if nvelocity:
            velocity = response.eigr_eigi_velocity[:, -1]
            if np.all(np.isfinite(velocity)):
                items = [f'V{ipoint+1}={vel}' for ipoint, vel in enumerate(velocity)]
            else:
                items = [f'V{ipoint+1}' for ipoint in range(nvelocity)]
            self.velocity_edit.addItems(items)

    def _get_subcase(self) -> tuple[int, bool]:
        subcase_str = self.subcase_edit.currentText()
        # self.log.info(f'_get_subcase: subcase_str={subcase_str!r}')
        subcase_sline = subcase_str.split()
        # self.log.info(f'_get_subcase: subcase_sline={subcase_sline}')
        try:
            subcase = int(subcase_sline[1])
            is_valid = True
        except IndexError:
            self.log.error(f'failed parsing subcase={subcase_str!r}; subcase_sline={subcase_sline}')
            subcase = -1
            is_valid = False
        self.log.info(f'_get_subcase: subcase={subcase}; is_valid={is_valid}')
        return subcase, is_valid

    def update_subcases(self, subcases: list[int]) -> None:
        subcases_text = [f'Subcase {isubcase}' for isubcase in subcases]
        self.log.info(f'update_subcases={subcases_text}')
        self.subcase_edit.clear()
        # self.log.info(f'update_subcases setting...')
        self.subcase_edit.addItems(subcases_text)

    def update_modes_table(self, modes: list[int],
                           freqs: list[float]) -> None:
        self.modes = modes
        _set_modes_table(self.modes_widget, modes, freqs)
        self.run_button.setEnabled(True)
        self.log.info(f'modes = {self.modes}')

    def on_modes(self) -> None:
        self.on_run()
        # self.validate()
        # self.plot(self.modes)

    def _on_update_mode(self) -> None:
        if not self.is_valid:
            # self.log.warning('_on_update_mode')
            self.validate()
        self.plot()

    def get_xlim(self) -> tuple[Limit, Limit, Limit, Limit,
                                Limit, Limit, Limit, Limit, Limit,
                                Optional[float], Optional[float], Optional[float],
                                Optional[float], Optional[float], bool]:
        index_lim, is_passed0 = get_list_float_or_none(
            [self.index_lim_edit_min, self.index_lim_edit_max])
        eas_lim, is_passed1 = get_list_float_or_none(
            [self.eas_lim_edit_min, self.eas_lim_edit_max])
        tas_lim, is_passed2 = get_list_float_or_none(
            [self.tas_lim_edit_min, self.tas_lim_edit_max])

        mach_lim, is_passed3 = get_list_float_or_none(
            [self.mach_lim_edit_min, self.mach_lim_edit_max])
        alt_lim, is_passed4 = get_list_float_or_none(
            [self.alt_lim_edit_min, self.alt_lim_edit_max])
        q_lim, is_passed5 = get_list_float_or_none(
            [self.q_lim_edit_min, self.q_lim_edit_max])
        rho_lim, is_passed6 = get_list_float_or_none(
            [self.rho_lim_edit_min, self.rho_lim_edit_max])

        is_passed_x = all([
            is_passed0, is_passed1, is_passed2, is_passed3,
            is_passed4, is_passed5, is_passed6,
        ])
        damp_lim, is_passed_damp = get_list_float_or_none(
            [self.damp_lim_edit_min, self.damp_lim_edit_max])

        freq_lim, is_passed_freq = get_list_float_or_none(
            [self.freq_lim_edit_min, self.freq_lim_edit_max])
        kfreq_lim, is_passed_kfreq = get_list_float_or_none(
            [self.kfreq_lim_edit_min, self.kfreq_lim_edit_max])
        ikfreq_lim, is_passed_ikfreq = get_list_float_or_none(
            [self.ikfreq_lim_edit_min, self.ikfreq_lim_edit_max])
        eigr_lim, is_passed_eigr = get_list_float_or_none(
            [self.eigr_lim_edit_min, self.eigr_lim_edit_max])
        eigi_lim, is_passed_eigi = get_list_float_or_none(
            [self.eigr_lim_edit_min, self.eigr_lim_edit_max])
        is_passed_eig = all([is_passed_eigr, is_passed_eigi])

        freq_tol, is_passed_tol1 = get_float_or_none(self.freq_tol_edit)
        freq_tol_remove, is_passed_tol2 = get_float_or_none(self.freq_tol_remove_edit)
        mag_tol, is_passed_tol3 = get_float_or_none(self.mag_tol_edit)
        if is_passed_tol1 and freq_tol is None:
            freq_tol = -1.0
        if is_passed_tol2 and freq_tol_remove is None:
            freq_tol_remove = -1.0
        if is_passed_tol3 and mag_tol is None:
            mag_tol = -1.0

        vl, is_passed_vl = get_float_or_none(self.vl_edit)
        vf, is_passed_vf = get_float_or_none(self.vf_edit)
        damping_required, is_passed_damping_required = get_float_or_none(self.damping_required_edit)
        damping_required_tol, is_passed_damping_required_tol = get_float_or_none(self.damping_required_tol_edit)
        damping, is_passed_damping = get_float_or_none(self.damping_edit)
        eas_flutter_range, is_passed_flutter_range = get_list_float_or_none(
            [self.eas_flutter_range_edit_min, self.eas_flutter_range_edit_max])

        # 595:596,600:601
        point_removal_str = self.point_removal_edit.text().strip()
        point_removal = point_removal_str_to_point_removal(point_removal_str, self.log)

        vl = _float_passed_to_default(vl, is_passed_vl)
        vf = _float_passed_to_default(vf, is_passed_vf)
        log = self.log
        damping_required = _float_passed_to_default(damping_required, is_passed_damping_required)
        damping = _float_passed_to_default(damping, is_passed_damping)

        # if is_passed_eas_damping_lim1 and eas_damping_lim_min is None:
        #    eas_damping_lim_min = None
        # if is_passed_eas_damping_lim2 and eas_damping_lim_max is None:
        #    eas_damping_lim_max = None

        # self.log.info(f'XLim = {xlim}')
        is_passed_flags = [
            is_passed_x,
            is_passed_damp,
            is_passed_freq, is_passed_kfreq, is_passed_ikfreq,
            is_passed_eig,
            is_passed_tol1, is_passed_tol2, is_passed_tol3,
            is_passed_vl, is_passed_vf,
            is_passed_damping_required, is_passed_damping_required_tol,
            is_passed_damping,
            is_passed_flutter_range,
        ]
        is_passed = all(is_passed_flags)
        # if not is_passed:
        # self.log.warning(f'is_passed_flags = {is_passed_flags}')
        # print(f'freq_tol = {freq_tol}')
        out = (
            index_lim, eas_lim, tas_lim, mach_lim, alt_lim, q_lim, rho_lim,
            damp_lim, freq_lim, kfreq_lim, ikfreq_lim,
            eigr_lim, eigi_lim,
            freq_tol, freq_tol_remove, mag_tol,
            vl, vf, damping_required, damping_required_tol, damping,
            eas_flutter_range, point_removal, is_passed,
        )
        return out

    def get_selected_modes(self) -> list[int]:
        mode_strs = get_selected_items_flat(self.modes_widget)
        # self.log.info(f'mode_strs = {mode_strs}')
        modes = [int(mode_str.split(';')[0].split(' ')[1])
                 for mode_str in mode_strs]
        self.log.info(f'modes = {modes}')
        return modes

    # @dont_crash
    def get_lineedits(self) -> list[tuple[str, int, QLineEdit]]:
        out = [
            ('recent_files', 0, self.f06_filename_edit),
            ('freq_tol', -1, self.freq_tol_edit),
            ('freq_tol_remove', -1, self.freq_tol_remove_edit),
            ('mag_tol', -1, self.mag_tol_edit),
            ('vl', -1, self.vl_edit),
            ('vf', -1, self.vf_edit),
            ('damping', -1, self.damping_edit),
            ('damping_required', -1, self.damping_required_edit),
            ('damping_required_tol', -1, self.damping_required_tol_edit),
            ('output_directory', -1, self.output_directory_edit),
        ]
        return out

    def get_comboboxs(self) -> list[tuple[str, QComboBox, list[str]]]:
        out = [
            ('x_plot_type', self.x_plot_type_pulldown, X_PLOT_TYPES),
            ('plot_type', self.plot_type_pulldown, PLOT_TYPES),
            ('units_in', self.units_in_pulldown, UNITS_IN),
            ('units_out', self.units_out_pulldown, UNITS_OUT),
            ('mode_switch_method', self.mode_switch_method_pulldown, MODE_SWITCH_METHODS),
        ]
        return out

    def get_checkboxs(self) -> list[tuple[str, QCheckBox]]:
        out = [
            ('use_rhoref', self.use_rhoref_checkbox),
            ('show_points', self.show_points_checkbox),
            ('show_mode_number', self.show_mode_number_checkbox),
            ('show_detailed_mode_info', self.show_detailed_mode_info_checkbox),
            ('show_lines', self.show_lines_checkbox),
        ]
        return out

    def get_min_max_lineedits(self) -> list[tuple[str, QLineEdit, QLineEdit]]:
        out = [
            ('eas_lim', self.eas_lim_edit_min, self.eas_lim_edit_max),
            ('tas_lim', self.tas_lim_edit_min, self.tas_lim_edit_max),
            ('mach_lim', self.mach_lim_edit_min, self.mach_lim_edit_max),
            ('alt_lim', self.alt_lim_edit_min, self.alt_lim_edit_max),
            ('q_lim', self.q_lim_edit_min, self.q_lim_edit_max),
            ('rho_lim', self.rho_lim_edit_min, self.rho_lim_edit_max),
            ('ikfreq_lim', self.ikfreq_lim_edit_min, self.ikfreq_lim_edit_max),

            ('damp_lim', self.damp_lim_edit_min, self.damp_lim_edit_max),
            ('eas_flutter_range', self.eas_flutter_range_edit_min, self.eas_flutter_range_edit_max),
            # ('eas_diverg_range', self.eas_diverg_range_edit_min, self.eas_diverg_range_edit_max),
            ('freq_lim', self.freq_lim_edit_min, self.freq_lim_edit_max),
            ('kfreq_lim', self.kfreq_lim_edit_min, self.kfreq_lim_edit_max),
        ]
        return out

    def validate(self) -> bool:
        ifile = self.ifile
        parent = self.parent
        log = self.log
        # log.warning('validate')
        (index_lim, eas_lim, tas_lim, mach_lim, alt_lim, q_lim, rho_lim,
         ydamp_lim, freq_lim, kfreq_lim, ikfreq_lim,
         eigr_lim, eigi_lim,
         freq_tol, freq_tol_remove, mag_tol,
         vl, vf, damping_required, damping_required_tol, damping,
         eas_flutter_range, point_removal,
         is_valid_xlim) = self.get_xlim()

        selected_modes = []
        subcase, is_subcase_valid = self._get_subcase()
        log.warning(f'subcase={subcase}; is_subcase_valid={is_subcase_valid}')
        if is_subcase_valid:
            selected_modes = self.get_selected_modes()

        self.subcase = subcase
        parent.subcase = subcase
        self.selected_modes = selected_modes
        parent.selected_modes = selected_modes

        parent.index_lim = index_lim
        parent.eas_lim = eas_lim
        parent.tas_lim = tas_lim
        parent.mach_lim = mach_lim
        parent.alt_lim = alt_lim
        parent.q_lim = q_lim
        parent.rho_lim = rho_lim
        parent.ikfreq_lim = ikfreq_lim
        parent.ydamp_lim = ydamp_lim
        parent.kfreq_lim = kfreq_lim
        parent.freq_lim = freq_lim
        parent.eigi_lim = eigi_lim
        parent.eigr_lim = eigr_lim
        parent.freq_tol = freq_tol
        parent.freq_tol_remove = freq_tol_remove
        parent.mag_tol = mag_tol
        parent.vl = vl
        parent.vf = vf
        parent.damping_required = damping_required
        parent.damping_required_tol = damping_required_tol
        parent.damping = damping
        parent.eas_flutter_range = eas_flutter_range
        parent.point_removal = point_removal

        parent.x_plot_type = self.x_plot_type_pulldown.currentText()
        parent.plot_type = self.plot_type_pulldown.currentText()
        parent.mode_switch_method = self.mode_switch_method_pulldown.currentText()

        units_in = self.units_in_pulldown.currentText()
        units_out = self.units_out_pulldown.currentText()
        output_directory = self.output_directory_edit.text()

        parent.show_lines = self.show_lines_checkbox.isChecked()
        parent.show_points = self.show_points_checkbox.isChecked()
        parent.show_mode_number = self.show_mode_number_checkbox.isChecked()
        parent.show_detailed_mode_info = self.show_detailed_mode_info_checkbox.isChecked()
        parent.point_spacing = self.point_spacing_spinner.value()
        parent.use_rhoref = self.use_rhoref_checkbox.isChecked()

        is_passed_modal_partipation = False
        subcases = list(self.responses)
        if len(subcases):
            log.info(f'subcases={subcases}')
            subcase0 = subcases[0]
            response = self.responses[subcase0]

            failed_modal_partipation = (
                (parent.plot_type == 'modal-participation') and
                ((response.eigr_eigi_velocity is None) or
                 (response.eigenvector is None))
            )
            is_passed_modal_partipation = not failed_modal_partipation
        # (
        #     (self.plot_type == 'modal-participation') and
        #     (response.eigr_eigi_velocity is not None)
        # ) or (self.plot_type != 'modal-participation'))
        data = {
            # 'bdf_filename': self.bdf_filename,
            # 'op2_filename': self.op2_filename,
            'log_scale_x': self.log_xscale_checkbox.isChecked(),
            'log_scale_y1': self.log_yscale1_checkbox.isChecked(),
            'log_scale_y2': self.log_yscale2_checkbox.isChecked(),
            'use_rhoref': parent.use_rhoref,
            'show_points': parent.show_points,
            'show_mode_number': parent.show_mode_number,
            'show_detailed_mode_info': parent.show_detailed_mode_info,
            'point_spacing': parent.point_spacing,
            'show_lines': parent.show_lines,

            'recent_files': parent.recent_files,
            'subcase': subcase,
            # 'modes': modes,
            'selected_modes': selected_modes,
            'x_plot_type': parent.x_plot_type,
            'plot_type': parent.plot_type,
            'index_lim': index_lim,
            'eas_lim': eas_lim,
            'tas_lim': tas_lim,
            'mach_lim': mach_lim,
            'alt_lim': alt_lim,
            'q_lim': q_lim,
            'rho_lim': rho_lim,
            'ikfreq_lim': ikfreq_lim,

            'damp_lim': ydamp_lim,
            'freq_lim': freq_lim,
            'kfreq_lim': kfreq_lim,
            'eigr_lim': eigr_lim,
            'eigi_lim': eigi_lim,
            'output_directory': output_directory,
            'units_in': units_in,
            'units_out': units_out,
            'freq_tol': freq_tol,
            'freq_tol_remove': freq_tol_remove,
            'mag_tol': mag_tol,
            'vl': vl,
            'vf': vf,
            'damping': damping, # 0.03
            'damping_required': damping_required, # 0.0
            'damping_required_tol': damping_required_tol, # 0.0001
            'eas_flutter_range': eas_flutter_range,
            'point_removal': point_removal,
            'mode_switch_method': parent.mode_switch_method,
        }
        parent.units_in = units_in
        parent.units_out = units_out
        is_passed = all([is_valid_xlim, is_subcase_valid, is_passed_modal_partipation])
        if is_passed:
            parent.data = data
            # self.xlim = xlim
            # self.ylim = ydamp_lim
            # self.data = data
            # is_valid = validate_json(self.data, self.log)
            # if is_valid != is_passed:
            # self.log.info(f'passed data:\n{str(self.data)}')
        else:
            del data['recent_files']
            log.error(
                f'is_valid_xlim = {is_valid_xlim}\n'
                f'is_subcase_valid = {is_subcase_valid}\n'
                f'is_passed_modal_partipation = {is_passed_modal_partipation}\n'
                f'failed data:\n{str(data)}'
            )
            # self.log.error(f'failed data:\n{str(data)}')
        return is_passed

    def on_run(self) -> None:
        # self.log.warning('on_run')
        parent = self.parent
        is_valid = self.validate()
        if not is_valid:
            return

        modes = self.selected_modes
        if len(modes) == 0:
            self.log.warning(f'modes = {modes}; assuming all modes -> {self.modes}')
            modes = self.modes
            # return
        self.log.info(f'is_valid = {is_valid}\n')
        parent.is_valid = True
        self.plot(modes)
        # self.log.warning('on_run; _save')
        parent._save(parent.save_filename)

    def plot(self, modes: list[int]) -> None:
        parent = self.parent
        log = parent.log
        log.info(f'plot; modes = {modes}\n')
        if not parent.is_valid:
            log.warning('not valid\n')
            return
        if len(self.responses) == 0:
            log.warning('no subcases\n')
            return

        x_plot_type = parent.x_plot_type
        plot_type = parent.plot_type
        log.info(f'plot_type = {plot_type}\n')

        freq_tol = parent.freq_tol
        freq_tol_remove = parent.freq_tol_remove
        mag_tol = parent.mag_tol
        log.info(f'freq_tol = {freq_tol}\n')
        noline, nopoints = get_noline_nopoints(
            parent.show_lines, parent.show_points)

        if x_plot_type == 'index':
            xlim = parent.index_lim
        elif x_plot_type == 'eas':
            xlim = parent.eas_lim
        elif x_plot_type == 'tas':
            xlim = parent.tas_lim
        elif x_plot_type == 'mach':
            xlim = parent.mach_lim
        elif x_plot_type == 'alt':
            xlim = parent.alt_lim
        elif x_plot_type == 'q':
            xlim = parent.q_lim
        elif x_plot_type == 'kfreq':
            xlim = parent.kfreq_lim
        elif x_plot_type == 'ikfreq':
            xlim = parent.ikfreq_lim
        else:  # pragma: no cover
            log.error(f'x_plot_type={x_plot_type!r} is not supported')
            # raise RuntimeError(x_plot_type)
            xlim = (None, None)

        # log.info(f'xlim={xlim}\n')
        if plot_type == 'zimmerman':
            print('skipping xlim check')
        else:
            assert xlim[0] != '' and xlim[1] != '', (xlim, x_plot_type)

        v_lines = get_vlines(parent.vf, parent.vl)
        # log.info(f'v_lines={v_lines}\n')
        # log.info(f'kfreq_lim={self.kfreq_lim}\n')
        # log.info(f'ydamp_lim={self.ydamp_lim}\n')
        # log.info(f'freq_lim={self.freq_lim}\n')
        # log.info(f'damping={self.damping}\n')
        xlim_kfreq = parent.kfreq_lim
        ylim_damping = parent.ydamp_lim
        ylim_freq = parent.freq_lim

        damping_limit = parent.damping  # % damping
        eas_flutter_range = parent.eas_flutter_range
        # eas_diverg_range = self.eas_diverg_range
        if plot_type not in {'root-locus'}:
            damping_required = parent.damping_required
            damping_required_tol = parent.damping_required_tol
            assert isinstance(damping_required_tol, float), f'damping_required_tol={damping_required_tol}'
            if damping_required_tol is None:
                damping_required_tol = 0.01
            if damping_required_tol < 0.0:
                damping_required_tol = 0.0

        # changing directory so we don't make a long filename
        # in the plot header
        # log.info(f'damping_limit = {damping_limit}\n')
        dirname = os.path.abspath(os.path.dirname(parent.f06_filename))
        basename = os.path.basename(parent.f06_filename)

        base = os.path.splitext(basename)[0]
        current_directory = os.getcwd()
        sys.stdout.flush()
        os.chdir(dirname)

        fig = plt.figure(1)
        fig.clear()
        log.info(f'cleared plot\n')
        if plot_type not in {'root-locus', 'modal-participation', 'zimmerman'}:
            gridspeci = gridspec.GridSpec(2, 4)
            damp_axes = fig.add_subplot(gridspeci[0, :3])
            freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        response = self.responses[self.subcase]

        # you can change the output units without reloading
        if parent._units_out != parent.units_out:
            response.convert_units(parent.units_out)
            parent._units_out = parent.units_out

        response.noline = noline
        response.freq_ndigits = parent.freq_ndigits
        response.set_symbol_settings(
            nopoints, parent.show_mode_number, parent.point_spacing)
        # log.info(f'self.plot_font_size = {self.plot_font_size}')
        response.set_font_settings(parent.plot_font_size)
        response.log = log
        # print('trying plots...')

        # log.info(f'getting logs\n')
        data = parent.data
        log_scale_x = data['log_scale_x']
        log_scale_y1 = data['log_scale_y1']
        log_scale_y2 = data['log_scale_y2']
        # print(f'log_scale_x={log_scale_x}; log_scale_y1={log_scale_y1}; log_scale_y2={log_scale_y2}')
        # print(f'export_to_png={self.export_to_png}')

        # print(f'point_removal = {self.point_removal}')
        png_filename0, png_filename = get_png_filename(
            base, x_plot_type, plot_type,
            parent.export_to_png)
        print(f'png_filename={png_filename}')
        try:
            if plot_type == 'zimmerman':
                axes = fig.add_subplot(111)
                response.plot_zimmerman(fig=fig, axes=axes, modes=modes, show=True)
            elif plot_type == 'root-locus':
                axes = fig.add_subplot(111)
                # log.info(f'modes={modes}; eigr_lim={self.eigr_lim}; eigi_lim={self.eigi_lim}; freq_tol={freq_tol}')
                # log.info(f'png_filename={png_filename}')
                response.plot_root_locus(
                    fig=fig, axes=axes,
                    modes=modes, eigr_lim=parent.eigr_lim, eigi_lim=parent.eigi_lim,
                    freq_tol=freq_tol, freq_tol_remove=freq_tol_remove,
                    show=True, clear=False, close=False,
                    legend=True,
                    png_filename=png_filename,
                )
            elif plot_type == 'modal-participation':
                axes = fig.add_subplot(111)
                mode = self.mode_edit.value()
                ivel = self.velocity_edit.currentIndex()
                # print(f'ivel={ivel}; mode={mode}')
                response.plot_modal_participation(
                    ivel, mode,
                    fig=fig, axes=axes,
                    modes=modes,  # eigr_lim=self.eigr_lim, eigi_lim=self.eigi_lim,
                    freq_tol=freq_tol, freq_tol_remove=freq_tol_remove,
                    mag_tol=mag_tol,
                    show=True, clear=False, close=False,
                    legend=True,
                    png_filename=png_filename,
                )
            elif plot_type == 'x-damp-kfreq':
                # xlabel: eas
                # ylabel1 = r'Structural Damping; $g = 2 \gamma $'
                # ylabel2 = r'KFreq [rad]; $ \omega c / (2 V)$'
                # print('plot_kfreq_damping')
                response.plot_kfreq_damping(
                    fig=fig, damp_axes=damp_axes, freq_axes=freq_axes,
                    modes=modes, plot_type=x_plot_type,
                    xlim=xlim, ylim_damping=ylim_damping, ylim_kfreq=xlim_kfreq,
                    freq_tol=freq_tol, freq_tol_remove=freq_tol_remove,
                    show=True, clear=False, close=False,
                    legend=True,
                    png_filename=png_filename,
                )
                update_ylog_style(fig, log_scale_x, log_scale_y1, log_scale_y2)
            else:
                assert plot_type in 'x-damp-freq', plot_type
                # print('plot_vg_vf')
                # log.info(f'png_filename={png_filename!r}')
                # log.info(f'modes={modes!r}')
                # log.info(f'freq_tol={freq_tol!r}')
                # log.info(f'v_lines={v_lines!r}')
                damping_crossings = get_damping_crossings(
                    damping_required, damping_required_tol, damping_limit)

                response.plot_vg_vf(
                    fig=fig, damp_axes=damp_axes, freq_axes=freq_axes,
                    plot_type=x_plot_type,
                    modes=modes,
                    xlim=xlim, ylim_damping=ylim_damping, ylim_freq=ylim_freq,
                    eas_range=eas_flutter_range,
                    freq_tol=freq_tol, freq_tol_remove=freq_tol_remove,
                    show=True, clear=False, close=False,
                    legend=True,
                    v_lines=v_lines,
                    damping_limit=damping_limit,
                    damping_required=damping_required,
                    damping_crossings=damping_crossings,
                    png_filename=png_filename,
                    point_removal=parent.point_removal,
                    mode_switch_method=parent.mode_switch_method,
                    show_detailed_mode_info=parent.show_detailed_mode_info,
                    ncol=parent.flutter_ncolumns,
                    divergence_legend_loc=parent.divergence_legend_loc,
                    flutter_bbox_to_anchor=(parent.flutter_bbox_to_anchor_x, 1.),
                )
                update_ylog_style(fig, log_scale_x, log_scale_y1, log_scale_y2)
                fig.canvas.draw()
        except Exception as e:  # pragma: no cover
            log.error(f'plot_type={plot_type}')
            log.error(str(e))
            print(traceback.format_exc())
            # print(traceback.print_tb())
            print(traceback.print_exception(e))
            raise

        export_flutter_results(
            response, modes, png_filename0,
            parent.export_to_csv, parent.export_to_zaero, parent.export_to_f06,
            log)
        os.chdir(current_directory)
        if png_filename:
            log.info(f'saved {png_filename}')
        else:
            log.info(f'did not write file because export_to_png=False')
