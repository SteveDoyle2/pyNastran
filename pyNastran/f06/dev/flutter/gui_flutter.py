import os
import sys
import copy
import warnings
import traceback
from pathlib import Path
from typing import Optional, Any  #TYPE_CHECKING

ICON_PATH = Path('')
try:
    import json5 as json
except ModuleNotFoundError:
    warnings.warn('couldnt find json5, using json')
    import json
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from qtpy import QtCore
from qtpy.compat import getopenfilename #, getsavefilename
from qtpy.QtWidgets import (
    QLabel, QWidget,
    QApplication, QMenu, QVBoxLayout, QLineEdit, QComboBox,
    QHBoxLayout, QPushButton, QGridLayout,
    QAction,
    QCheckBox, #QRadioButton,
    QListWidgetItem, QAbstractItemView,
    QListWidget, QSpinBox,
)
# from qtpy.QtWidgets import (
#     QMessageBox,
#     QMainWindow, QDockWidget, QFrame, QToolBar,
#     QToolButton, QMenuBar,
# )
#from qtpy.QtGui import QIcon
QLINEEDIT_WHITE = 'QLineEdit {background-color: white;}'
QLINEEDIT_RED = 'QLineEdit {background-color: red;}'

from cpylog import SimpleLogger
from pyNastran.gui.utils.qt.pydialog import QFloatEdit, make_font
from pyNastran.gui.qt_files.named_dock_widget import NamedDockWidget
from pyNastran.gui.qt_files.loggable_gui import LoggableGui

from pyNastran.f06.dev.flutter.actions_builder import Actions, Action, build_menus
from pyNastran.f06.dev.flutter.preferences_object import PreferencesObject
from pyNastran.f06.dev.flutter.vtk_window_object import VtkWindowObject

from pyNastran.f06.flutter_response import FlutterResponse, Limit
from pyNastran.f06.parse_flutter import get_flutter_units

X_PLOT_TYPES = ['eas', 'tas', 'rho', 'q', 'mach', 'alt', 'kfreq', 'ikfreq']
PLOT_TYPES = ['x-damp-freq', 'x-damp-kfreq', 'root-locus', 'modal-participation']
UNITS_IN = ['english_in', 'english_kt', 'english_ft',
            'si', 'si_mm']
UNITS_OUT = UNITS_IN

#FONT_SIZE = 12
from pyNastran.f06.dev.flutter.utils import (
    get_plot_file, update_ylog_style, load_f06_op2,
    get_png_filename,)

import pyNastran
PKG_PATH = Path(pyNastran.__path__[0])
HOME_FILENAME = get_plot_file()

AERO_PATH = PKG_PATH / '..' / 'models' / 'aero'
if 0:  # pragma: no cover
    BASE_PATH = AERO_PATH / 'flutter_bug'
    BDF_FILENAME = BASE_PATH / 'nx' / 'wing_b1.bdf'
    OP2_FILENAME = BASE_PATH / 'wing_b1.op2'
elif 1:  # pragma: no cover
    BDF_FILENAME = PKG_PATH / '..' / 'models' / 'bwb' / 'bwb_saero.bdf'
    OP2_FILENAME = PKG_PATH / '..' / 'models' / 'bwb' / 'bwb_saero.op2'
else:
    BASE_PATH = AERO_PATH / '2_mode_flutter'
    BDF_FILENAME = BASE_PATH / '0012_flutter.bdf'
    OP2_FILENAME = BASE_PATH / '0012_flutter.op2'


class FlutterGui(LoggableGui):
    def __init__(self, f06_filename: str=''):
        super().__init__(html_logging=False)

        self._export_settings_obj = PreferencesObject(self)
        self._vtk_window_obj = VtkWindowObject(self, ICON_PATH)
        self.font_size = 10
        self.plot_font_size = 10
        self.show_lines = True
        self.show_points = True
        self.show_mode_number = False
        self.point_spacing = 0
        self._units_in = ''
        self._units_out = ''
        self.units_in = ''
        self.units_out = ''
        self.use_dock_widgets = self.html_logging
        self.qactions = {}
        self.nrecent_files_max = 20
        self.recent_files = []
        self.save_filename = HOME_FILENAME
        self.is_valid = False

        self.data = {}
        self.f06_filename = ''
        self.bdf_filename = ''
        self.op2_filename = ''
        self._bdf_filename_default = ''
        self._op2_filename_default = ''
        self.subcase = 0
        self.x_plot_type = 'eas'
        self.plot_type = 'x-damp-freq'
        self.eas_lim = []
        self.tas_lim = []
        self.mach_lim = []
        self.alt_lim = []
        self.q_lim = []
        self.rho_lim = []
        #self.x_lim = []
        self.freq_lim = [None, None]
        self.damping_lim = [None, None]
        self.kfreq_lim = [None, None]
        self.eigr_lim = [None, None]
        self.eigi_lim = [None, None]
        self.responses = {}
        self.modes = []
        self.selected_modes = []
        self.freq_tol = -1.0
        self.freq_tol_remove = -1.0
        self.mag_tol = -1.0
        self.damping = -1.0
        self.vf = -1.0
        self.vl = -1.0
        self.export_to_png = True
        self.export_to_csv = False
        self.export_to_f06 = False
        self.export_to_zona = False

        self.setup_widgets()
        self.setup_layout()
        self.on_load_settings()
        if f06_filename:
            self.f06_filename_edit.setText(f06_filename)
            self._set_f06_default_names(f06_filename)

        self.setup_toolbar()
        self._update_recent_files_actions()
        self.setup_connections()
        self._set_window_title()
        self.on_font_size()
        self.on_plot_type()
        self._set_f06_default_names(self.f06_filename_edit.text())
        #self.on_open_new_window()
        self.show()

    def setup_toolbar(self):
        #frame = QFrame(self)
        actions_dict = {
            #'file_load': Action(name='file_load', text='Load...', func=self.on_file_load, icon='folder.png'),
            #'file_save': Action(name='file_save', text='Save...', func=self.on_file_save, icon='save.png'),
            #'file_save_as': Action(name='file_save_as', text='Save As...', func=self.on_file_save_as),
            'file_exit':       Action(name='exit', text='Exit...', icon='exit2.jpg', func=self.on_file_exit),
            'export_settings': Action(name='Export Settings', text='Export Settings...', icon='preferences.jpg',
                                      shortcut='Ctrl+P',func=self.on_export_settings),
        }
        actions_input = Actions(ICON_PATH, actions_dict) # , load_icon=False
        recent_files = actions_input.build_recent_file_qactions(
            self, self.recent_files, self.set_f06)
        self.qactions = actions_input.build_qactions(self)

        file_actions = [
            #'file_load',
            # 'file_save', 'file_save_as',
            ] + recent_files + [
            'file_exit']
        view_actions = ['export_settings']

        self.menubar = self.menuBar()
        self.file_menu = self.menubar.addMenu('File')
        self.view_menu = self.menubar.addMenu('View')
        #self.help_menu = self.menubar.addMenu('Help')

        help_actions = []
        menus_dict = {
            'File': (self.file_menu, file_actions),
            'View': (self.view_menu, view_actions),
            #'Help': (self.help_menu, help_actions),
        }
        build_menus(menus_dict, self.qactions)
        #self.file_menu.addAction(actions['file_load'])
        #self.file_menu.addAction(actions['exit'])

        #self.toolbar = self.addToolBar('Show toolbar')
        #self.toolbar.setObjectName('main_toolbar')
        self.statusbar = self.statusBar()

    def set_f06(self, ifile: int) -> None:
        f06_filename = self.recent_files[ifile]
        self.f06_filename_edit.setText(f06_filename)
        self.on_load_f06()

    def setup_modes(self):
        self.modes_widget = QListWidget(self)
        self.modes_widget.setMaximumWidth(200)  # was 100 when no freq
        self.modes_widget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self._set_modes_table(self.modes_widget, [0], [0.])

    # def dontcrash(func):
    #     @wraps(func)
    #     def wrapper(self, *args, **kwargs):
    #         # do something before `sum`
    #         result = func(self) # , *args, **kwargs
    #         # do something after `sum`
    #         return result
    #     return wrapper

    # @dontcrash
    def on_export_settings(self):
        self._export_settings_obj.show()

    def on_file_exit(self):
        if hasattr(self, 'on_file_save') and hasattr(self, 'save_filename'):
            self.on_file_save()

    def on_file_save(self) -> None:
        if self.save_filename == '' or not os.path.exists(self.save_filename):
            self.on_file_save_as()
        else:
            #self.log.warning('on_file_save; _save')
            self._save(self.save_filename)

    def on_file_save_as(self) -> None:
        # if 0:
        #     #print('on_file_save')
        #     qt_wildcard = '*.json'
        #     basedir = str(DIRNAME)
        #     json_filename, wildcard = getsavefilename(
        #         self, caption=title, basedir=basedir,
        #         filters=qt_wildcard,
        #         #options=QFileDialog.setLabelText('data.json'),
        #     )
        #     self.log.info(f'json_filename={json_filename!r} wildcard={wildcard!r}')
        json_filename = self.save_filename
        #self.log.warning('on_file_save_as; _save')
        self._save(json_filename)

    def _save(self, json_filename: str):
        #self.log.warning('_save')
        is_valid = self.validate()
        #self.log.info(f'self.data = {self.data}')
        if json_filename == '' or len(self.data) == 0:
            return
        #print(f'json_filename={json_filename!r} wildcard={wildcard!r}')
        #print(f'self.data = {self.data}')
        out_data = copy.deepcopy(self.data)
        out_data['vtk'] = self._vtk_window_obj.data
        with open(json_filename, 'w') as json_file:
            json.dump(out_data, json_file, indent=4)
        #print(f'fname="{fname}"')
        self.log.info(f'finished saving {json_filename!r}\n')
        self.save_filename = json_filename
        self._set_window_title()

    def _apply_settings(self, data: dict[str, Any]) -> None:
        self._vtk_window_obj.apply_settings(data)
        font_size0 = self.font_size
        # radios = [
        #     ('show_points', self.show_points_radio),
        # ]
        # for (key, checkbox) in radios:
        #     if key not in data:
        #         continue
        #     val = data[key]
        #     assert isinstance(val, bool), (key, val)
        #     checkbox.setChecked(val)

        spinners = [
             #('plot_font_size', self.plot_font_size_edit),
        ]
        for (key, spinner) in spinners:
            if key not in data:
                continue
            val = data[key]
            assert isinstance(val, int), (key, val)
            spinner.setValue(val)

        type_names = [
            (int,  ('font_size', 'plot_font_size',)),
            (bool, ('export_to_png',
                    'export_to_f06', 'export_to_csv', 'export_to_zona')),
        ]
        for value_type, keys in type_names:
            for key in keys:
                if key not in data:
                    print(f'skipping {key!r}')
                    assert len(key) > 1, keys
                    continue
                value = data[key]
                assert isinstance(value, value_type), (key, value, value_type)
                assert hasattr(self, key), (key, value)
                setattr(self, key, value)

        checkboxs = [
            ('use_rhoref', self.use_rhoref_checkbox),
            ('show_points', self.show_points_checkbox),
            ('show_mode_number', self.show_mode_number_checkbox),
            ('show_lines', self.show_lines_checkbox),
        ]
        # attrs aren't stored
        for (key, checkbox) in checkboxs:
            if key not in data:
                continue
            val = data[key]
            assert isinstance(val, bool), (key, val)
            checkbox.setChecked(val)

        min_max_line_edits = [
            ('eas_lim', self.eas_lim_edit_min, self.eas_lim_edit_max),
            ('tas_lim', self.tas_lim_edit_min, self.tas_lim_edit_max),
            ('mach_lim', self.mach_lim_edit_min, self.mach_lim_edit_max),
            ('alt_lim', self.alt_lim_edit_min, self.alt_lim_edit_max),
            ('q_lim', self.q_lim_edit_min, self.q_lim_edit_max),
            ('rho_lim', self.rho_lim_edit_min, self.rho_lim_edit_max),
            ('xlim', self.xlim_edit_min, self.xlim_edit_max),

            ('damp_lim', self.damp_lim_edit_min, self.damp_lim_edit_max),
            ('freq_lim', self.freq_lim_edit_min, self.freq_lim_edit_max),
            ('kfreq_lim', self.kfreq_lim_edit_min, self.kfreq_lim_edit_max),
        ]
        for key, line_edit_min, line_edit_max in min_max_line_edits:
            if key not in data:
                #print(f'apply_settings: skipping key={key!r}')
                continue
            values = data[key]
            value0 = _to_str(values[0])
            value1 = _to_str(values[1])
            line_edit_min.setText(value0)
            line_edit_max.setText(value1)

        line_edits = [
            ('recent_files', 0, self.f06_filename_edit),
            ('freq_tol', -1, self.freq_tol_edit),
            ('freq_tol_remove', -1, self.freq_tol_remove_edit),
            ('mag_tol', -1, self.mag_tol_edit),
            ('vl', -1, self.VL_edit),
            ('vf', -1, self.VF_edit),
            ('damping', -1, self.damping_edit),
            ('output_directory', -1, self.output_directory_edit),
        ]
        for key, index, line_edit in line_edits:
            if key not in data:
                #print(f'apply_settings: skipping key={key!r}')
                continue
            values = data[key]
            if index != -1:
                value = values[index]
            else:
                value = values

            str_value = _to_str(value)

            #print('type(value) =', type(value))
            #print(f'{key+":":<10} values={values}[{index!r}]={value!r} -> {str_value!r}')
            line_edit.setText(str_value)

        pulldown_edits = [
            ('x_plot_type', self.x_plot_type_pulldown, X_PLOT_TYPES),
            ('plot_type', self.plot_type_pulldown, PLOT_TYPES),
            ('units_in', self.units_in_pulldown, UNITS_IN),
            ('units_out', self.units_out_pulldown, UNITS_OUT),
        ]
        for key, pulldown_edit, values in pulldown_edits:
            value = data[key]
            index = values.index(value)
            pulldown_edit.setCurrentIndex(index)

        self.recent_files = []
        for fname in data['recent_files']:
            abs_path = os.path.abspath(fname)
            if abs_path not in self.recent_files:
                self.recent_files.append(abs_path)
        self.f06_filename = self.recent_files[0]
        if self.font_size != font_size0:
            self.on_set_font_size(self.font_size)

    def on_enable_bdf(self) -> None:
        state = self.bdf_filename_checkbox.isChecked()
        self.bdf_filename_edit.setEnabled(state)
        self.bdf_filename_browse.setEnabled(state)
        if state:
            self.bdf_filename_edit.setText(self.bdf_filename)
        else:
            self.bdf_filename_edit.setText(self._bdf_filename_default)

    def on_enable_op2(self) -> None:
        state = self.op2_filename_checkbox.isChecked()
        self.op2_filename_edit.setEnabled(state)
        self.op2_filename_browse.setEnabled(state)
        if state:
            self.op2_filename_edit.setText(self.op2_filename)
        else:
            self.op2_filename_edit.setText(self._op2_filename_default)

    def on_browse_f06(self) -> None:
        """pops a dialog to select the f06 file"""
        title = 'Load Nastran Flutter F06 File'
        qt_wildcard = 'F06 File (*.f06)'
        basedir = os.path.dirname(self.f06_filename)
        fname, wildcard_level = getopenfilename(
            self, caption=title, basedir=basedir, filters=qt_wildcard,)
        if fname == '':
            return
        self.f06_filename_edit.setText(fname)
        self.ok_button.setEnabled(False)
        self._set_f06_default_names(fname)
    def _set_f06_default_names(self, f06_filename: str) -> None:
        base = os.path.splitext(f06_filename)[0]
        self._bdf_filename_default = base + '.bdf'
        self._op2_filename_default = base + '.op2'
        if self.bdf_filename == '':
            self.bdf_filename = self._bdf_filename_default
        if self.op2_filename == '':
            self.op2_filename = self._op2_filename_default
        self.bdf_filename_edit.setText(self._bdf_filename_default)
        self.op2_filename_edit.setText(self._op2_filename_default)

    def on_browse_bdf(self) -> None:
        """pops a dialog to select the bdf file"""
        title = 'Load Nastran Flutter BDF/DAT File'
        qt_wildcard = 'F06 File (*.bdf, *.dat)'
        basedir = os.path.dirname(self.bdf_filename)
        fname, wildcard_level = getopenfilename(
            self, caption=title, basedir=basedir, filters=qt_wildcard,)
        if fname == '':
            return
        self.bdf_filename_edit.setText(fname)
        self.bdf_filename = fname
        #self.ok_button.setEnabled(False)

    def on_browse_op2(self) -> None:
        """pops a dialog to select the op2 file"""
        title = 'Load Nastran Flutter OP2 File'
        qt_wildcard = 'OP2 File (*.o2p)'
        basedir = os.path.dirname(self.op2_filename)
        fname, wildcard_level = getopenfilename(
            self, caption=title, basedir=basedir, filters=qt_wildcard,)
        if fname == '':
            return
        self.op2_filename_edit.setText(fname)
        self.pop_vtk_gui_button.setEnabled(True)
        self.op2_filename = fname

    # @dontcrash
    def on_load_settings(self) -> None:
        json_filename = self.save_filename
        if not os.path.exists(json_filename):
            self.log.warning(f'unable to find {json_filename}')
            return
        # self.save_filename = json_filename
        try:
            with open(json_filename, 'r') as json_file:
                data = json.load(json_file)
            is_valid = validate_json(data, self.log)
            self._apply_settings(data)
        except Exception as e:
            self.log.error(f'failed to load {json_filename}\n{str(e)}')
            #raise
            return
        self.log.info(f'finished loading {json_filename!r}')
        #return wildcard_level, fname
        self._set_window_title()

    def _set_window_title(self) -> None:
        if self.save_filename == '':
            self.setWindowTitle('Flutter Plot')
        else:
            self.setWindowTitle(f'Flutter Plot: {self.save_filename}')

    def setup_widgets(self) -> None:
        self.f06_filename_label = QLabel('F06 Filename:')
        self.f06_filename_edit = QLineEdit()
        self.f06_filename_browse = QPushButton('Browse...')

        self.bdf_filename_checkbox = QCheckBox('BDF Filename:')
        self.bdf_filename_edit = QLineEdit()
        self.bdf_filename_browse = QPushButton('Browse...')
        self.bdf_filename_checkbox.setChecked(False)
        self.bdf_filename_edit.setEnabled(False)
        self.bdf_filename_browse.setEnabled(False)
        self.bdf_filename_edit.setToolTip('Loads the Nastran Geometry')

        self.op2_filename_checkbox = QCheckBox('OP2 Filename:')
        self.op2_filename_edit = QLineEdit()
        self.op2_filename_browse = QPushButton('Browse...')
        self.op2_filename_checkbox.setChecked(False)
        self.op2_filename_edit.setEnabled(False)
        self.op2_filename_browse.setEnabled(False)
        self.op2_filename_edit.setToolTip('Loads the Nastran Results (and geometry if BDF Filename is empty)')

        self.use_rhoref_checkbox = QCheckBox('Sea Level Rho Ref')
        self.use_rhoref_checkbox.setChecked(False)

        self.log_xscale_checkbox = QCheckBox('Log Scale x')
        self.log_yscale1_checkbox = QCheckBox('Log Scale y1')
        self.log_yscale2_checkbox = QCheckBox('Log Scale y2')
        self.log_xscale_checkbox.setChecked(False)
        self.log_yscale1_checkbox.setChecked(False)
        self.log_yscale2_checkbox.setChecked(False)

        self.show_points_checkbox = QCheckBox('Show Points')
        self.show_mode_number_checkbox = QCheckBox('Show Mode Number')
        self.point_spacing_label = QLabel('Point Spacing')
        self.point_spacing_spinner = QSpinBox()
        self.show_lines_checkbox = QCheckBox('Show Lines')
        self.show_points_checkbox.setChecked(True)
        self.show_lines_checkbox.setChecked(True)
        self.show_points_checkbox.setToolTip('The points are symbols')
        self.show_mode_number_checkbox.setToolTip('The points are the mode number')
        self.point_spacing_spinner.setToolTip('Skip Every Nth Point; 0=Plot All')
        self.point_spacing_spinner.setValue(0)
        self.point_spacing_spinner.setMinimum(0)
        self.point_spacing_spinner.setMaximum(10)

        self.eas_lim_label = QLabel('EAS Limits:')
        self.eas_lim_edit_min = QFloatEdit('0')
        self.eas_lim_edit_max = QFloatEdit()

        self.tas_lim_label = QLabel('TAS Limits:')
        self.tas_lim_edit_min = QFloatEdit('0')
        self.tas_lim_edit_max = QFloatEdit()

        self.mach_lim_label = QLabel('Mach Limits:')
        self.mach_lim_edit_min = QFloatEdit()
        self.mach_lim_edit_max = QFloatEdit()

        self.alt_lim_label = QLabel('Alt Limits:')
        self.alt_lim_edit_min = QFloatEdit()
        self.alt_lim_edit_max = QFloatEdit()

        self.q_lim_label = QLabel('Q Limits:')
        self.q_lim_edit_min = QFloatEdit()
        self.q_lim_edit_max = QFloatEdit()

        self.rho_lim_label = QLabel('Rho Limits:')
        self.rho_lim_edit_min = QFloatEdit('0')
        self.rho_lim_edit_max = QFloatEdit()

        self.xlim_label = QLabel('X Limits:')
        self.xlim_edit_min = QFloatEdit('0')
        self.xlim_edit_max = QFloatEdit()

        self.damp_lim_label = QLabel('Damping Limits (g):')
        self.damp_lim_edit_min = QFloatEdit('-0.3')
        self.damp_lim_edit_max = QFloatEdit('0.3')

        self.freq_lim_label = QLabel('Freq Limits (Hz):')
        self.freq_lim_edit_min = QFloatEdit('0')
        self.freq_lim_edit_max = QFloatEdit()

        self.kfreq_lim_label = QLabel('KFreq Limits:')
        self.kfreq_lim_edit_min = QFloatEdit()
        self.kfreq_lim_edit_max = QFloatEdit()


        self.eigr_lim_label = QLabel('Real Eigenvalue:')
        self.eigr_lim_edit_min = QFloatEdit()
        self.eigr_lim_edit_max = QFloatEdit()

        self.eigi_lim_label = QLabel('Imag Eigenvalue:')
        self.eigi_lim_edit_min = QFloatEdit()
        self.eigi_lim_edit_max = QFloatEdit()

        #--------------------------------------------
        self.freq_tol_label = QLabel('dFreq Tol (Hz) Dash:')
        self.freq_tol_edit = QFloatEdit('-1.0')
        self.freq_tol_edit.setToolTip("Applies a dotted line for modes that don't change by more than some amount")

        self.mag_tol_label = QLabel('Magnitude Tol:')
        self.mag_tol_edit = QFloatEdit('-1.0')
        self.mag_tol_edit.setToolTip("Filters modal participation factors based on magnitude")

        self.freq_tol_remove_label = QLabel('dFreq Tol (Hz) Hide:')
        self.freq_tol_remove_edit = QFloatEdit('-1.0')
        self.freq_tol_remove_edit.setToolTip("Completely remove modes that don't change by more than some amount")
        self.freq_tol_remove_label.setVisible(False)
        self.freq_tol_remove_edit.setVisible(False)

        self.subcase_label = QLabel('Subcase:')
        self.subcase_edit = QComboBox()

        units_msg = (
            "english_in: inch/s, slich/in^3\n"
            "english_ft: ft/s,   slug/ft^3\n"
            "english_kt: knots,  slug/ft^3\n"
            "si:         m/s,    kg/m^3\n"
            "si-mm:      mm/s,   Mg/mm^3\n"
        )
        self.x_plot_type_label = QLabel('X-Axis Plot Type:')
        self.x_plot_type_pulldown = QComboBox(self)
        self.x_plot_type_pulldown.addItems(X_PLOT_TYPES)
        self.x_plot_type_pulldown.setToolTip('sets the x-axis')

        self.plot_type_label = QLabel('Plot Type:')
        self.plot_type_pulldown = QComboBox(self)
        self.plot_type_pulldown.addItems(PLOT_TYPES)
        #self.plot_type_pulldown.setToolTip(units_msg)

        self.units_in_label = QLabel('Units In:')
        self.units_in_pulldown = QComboBox(self)
        self.units_in_pulldown.addItems(UNITS_IN)
        self.units_in_pulldown.setToolTip(units_msg)
        iunits_in = UNITS_IN.index('english_in')
        self.units_in_pulldown.setCurrentIndex(iunits_in)
        self.units_in_pulldown.setToolTip('Sets the units for the F06/OP2; set when loaded')

        self.units_out_label = QLabel('Units Out:')
        self.units_out_pulldown = QComboBox()
        self.units_out_pulldown.addItems(UNITS_OUT)
        self.units_out_pulldown.setToolTip(units_msg)
        iunits_out = UNITS_IN.index('english_kt')
        self.units_out_pulldown.setCurrentIndex(iunits_out)
        self.units_out_pulldown.setToolTip('Sets the units for the plot; may be updated')

        self.output_directory_label = QLabel('Output Directory:')
        self.output_directory_edit = QLineEdit('', self)
        self.output_directory_browse = QPushButton('Browse...', self)
        self.output_directory_edit.setDisabled(True)
        self.output_directory_browse.setDisabled(True)

        self.VL_label = QLabel('VL, Limit:')
        self.VL_edit = QFloatEdit('')
        self.VL_edit.setToolTip('Makes a vertical line for VL')

        self.VF_label = QLabel('VF, Flutter:')
        self.VF_edit = QFloatEdit('')
        self.VF_edit.setToolTip('Makes a vertical line for VF')

        self.damping_label = QLabel('Damping, g:')
        self.damping_edit = QFloatEdit('')
        self.damping_edit.setToolTip('Enables the flutter crossing (e.g., 0.03 for 3%)')

        self.mode_label = QLabel('Mode:')
        self.mode_edit = QSpinBox()
        self.mode_edit.setMinimum(1)
        #self.mode_edit.SetValue(3)
        self.mode_edit.setToolTip('Sets the mode')

        self.velocity_label = QLabel('Velocity Point:')
        self.velocity_edit = QComboBox()
        self.velocity_edit.setToolTip('Sets the velocity (input units)')

        self.f06_load_button = QPushButton('Load F06', self)
        self.ok_button = QPushButton('Run', self)


        self.pop_vtk_gui_button = QPushButton('Open GUI', self)
        self.solution_type_label = QLabel('Solution Type:')
        self.solution_type_pulldown = QComboBox(self)
        self.mode2_label = QLabel('Mode:')
        self.mode2_pulldown = QComboBox(self)
        self.setup_modes()
        self.on_plot_type()
        self.on_enable_bdf()
        self.on_enable_op2()

    def on_plot_type(self) -> None:
        x_plot_type = self.x_plot_type_pulldown.currentText()
        plot_type = self.plot_type_pulldown.currentText()

        show_eas_lim = False
        show_tas_lim = False
        show_mach_lim = False
        show_alt_lim = False
        show_q_lim = False
        show_rho_lim = False

        show_xlim = False
        show_freq = False
        show_damp = False
        show_root_locus = False
        show_modal_participation = False

        #PLOT_TYPES = ['x-damp-freq', 'x-damp-kfreq', 'root-locus']
        assert plot_type in PLOT_TYPES, plot_type
        self.on_units_out()

        if x_plot_type == 'kfreq':
            show_kfreq = True
        else:
            show_kfreq = False

        if plot_type == 'x-damp-freq':
            show_xlim = True
            show_damp = True
            show_freq = True
        elif plot_type == 'x-damp-kfreq':
            # kfreq-damp-kfreq not handled
            show_xlim = True
            show_damp = True
            show_kfreq = True
        elif plot_type == 'root-locus':
            show_root_locus = True
            #show_kfreq = False
        elif plot_type == 'modal-participation':
            show_modal_participation = True
            #show_kfreq = False
        else:  # pragma: no cover
            raise RuntimeError(f'plot_type={plot_type!r}')

        if show_xlim:
            if 'eas' == x_plot_type:
                show_eas_lim = True
                show_xlim = False
            elif 'tas' == x_plot_type:
                show_tas_lim = True
                show_xlim = False
            elif 'mach' == x_plot_type:
                show_mach_lim = True
                show_xlim = False
            elif 'alt' == x_plot_type:
                show_alt_lim = True
                show_xlim = False
            elif 'q' == x_plot_type:
                show_q_lim = True
                show_xlim = False
            elif 'rho' == x_plot_type:
                show_rho_lim = True
                show_xlim = False
        #print(f'x_plot_type={x_plot_type} show_damp={show_damp}; show_xlim={show_xlim}')
        #assert show_xlim is False, show_xlim

        show_eigenvalue = show_root_locus or show_modal_participation
        show_xaxis = not show_eigenvalue
        show_freq_tol = show_xaxis or show_root_locus
        show_crossing = show_xaxis

        self.x_plot_type_label.setVisible(show_xaxis)
        self.x_plot_type_pulldown.setVisible(show_xaxis)
        self.freq_tol_label.setVisible(show_freq_tol)
        self.freq_tol_edit.setVisible(show_freq_tol)

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

        self.xlim_label.setVisible(show_xlim)
        self.xlim_edit_min.setVisible(show_xlim)
        self.xlim_edit_max.setVisible(show_xlim)

        self.damp_lim_label.setVisible(show_damp)
        self.damp_lim_edit_min.setVisible(show_damp)
        self.damp_lim_edit_max.setVisible(show_damp)
        self.damping_label.setVisible(show_crossing)
        self.damping_edit.setVisible(show_crossing)

        self.freq_lim_label.setVisible(show_freq)
        self.freq_lim_edit_min.setVisible(show_freq)
        self.freq_lim_edit_max.setVisible(show_freq)

        self.kfreq_lim_label.setVisible(show_kfreq)
        self.kfreq_lim_edit_min.setVisible(show_kfreq)
        self.kfreq_lim_edit_max.setVisible(show_kfreq)

        self.eigr_lim_label.setVisible(show_root_locus)
        self.eigr_lim_edit_min.setVisible(show_root_locus)
        self.eigr_lim_edit_max.setVisible(show_root_locus)

        self.eigi_lim_label.setVisible(show_root_locus)
        self.eigi_lim_edit_min.setVisible(show_root_locus)
        self.eigi_lim_edit_max.setVisible(show_root_locus)

        self.VL_label.setVisible(show_xlim)
        self.VL_edit.setVisible(show_xlim)
        self.VF_label.setVisible(show_xlim)
        self.VF_edit.setVisible(show_xlim)

        show_items = [
            (show_modal_participation, (
                self.mode_label, self.mode_edit,
                self.velocity_label, self.velocity_edit,
                self.mag_tol_label, self.mag_tol_edit)),
            (not show_modal_participation, (
                self.point_spacing_label, self.point_spacing_spinner,
                self.show_mode_number_checkbox,
                self.show_lines_checkbox,
                self.log_xscale_checkbox,
                self.log_yscale1_checkbox, self.log_yscale2_checkbox,
            ),),
        ]
        for show_hide, items in show_items:
            for item in items:
                item.setVisible(show_hide)

    def setup_layout(self) -> None:
        if 0:
            hbox = QHBoxLayout()
            hbox.addWidget(self.f06_filename_label)
            hbox.addWidget(self.f06_filename_edit)
            hbox.addWidget(self.f06_filename_browse)
        else:
            file_row = 0
            hbox = QGridLayout()
            hbox.addWidget(self.f06_filename_label, file_row, 0)
            hbox.addWidget(self.f06_filename_edit, file_row, 1)
            hbox.addWidget(self.f06_filename_browse, file_row, 2)
            file_row += 1
            hbox.addWidget(self.bdf_filename_checkbox, file_row, 0)
            hbox.addWidget(self.bdf_filename_edit, file_row, 1)
            hbox.addWidget(self.bdf_filename_browse, file_row, 2)
            file_row += 1
            hbox.addWidget(self.op2_filename_checkbox, file_row, 0)
            hbox.addWidget(self.op2_filename_edit, file_row, 1)
            hbox.addWidget(self.op2_filename_browse, file_row, 2)
            file_row += 1

        grid = QGridLayout()
        irow = 0
        grid.addWidget(self.units_in_label, irow, 0)
        grid.addWidget(self.units_in_pulldown, irow, 1)
        grid.addWidget(self.use_rhoref_checkbox, irow, 2)
        irow += 1

        grid.addWidget(self.units_out_label, irow, 0)
        grid.addWidget(self.units_out_pulldown, irow, 1)
        irow += 1

        grid.addWidget(self.subcase_label, irow, 0)
        grid.addWidget(self.subcase_edit, irow, 1)
        irow += 1

        grid.addWidget(self.x_plot_type_label, irow, 0)
        grid.addWidget(self.x_plot_type_pulldown, irow, 1)
        irow += 1
        grid.addWidget(self.plot_type_label, irow, 0)
        grid.addWidget(self.plot_type_pulldown, irow, 1)
        irow += 1

        #--------------------------------------------------
        # x-axis
        grid.addWidget(self.eas_lim_label, irow, 0)
        grid.addWidget(self.eas_lim_edit_min, irow, 1)
        grid.addWidget(self.eas_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.tas_lim_label, irow, 0)
        grid.addWidget(self.tas_lim_edit_min, irow, 1)
        grid.addWidget(self.tas_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.mach_lim_label, irow, 0)
        grid.addWidget(self.mach_lim_edit_min, irow, 1)
        grid.addWidget(self.mach_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.alt_lim_label, irow, 0)
        grid.addWidget(self.alt_lim_edit_min, irow, 1)
        grid.addWidget(self.alt_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.q_lim_label, irow, 0)
        grid.addWidget(self.q_lim_edit_min, irow, 1)
        grid.addWidget(self.q_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.rho_lim_label, irow, 0)
        grid.addWidget(self.rho_lim_edit_min, irow, 1)
        grid.addWidget(self.rho_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.xlim_label, irow, 0)
        grid.addWidget(self.xlim_edit_min, irow, 1)
        grid.addWidget(self.xlim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.kfreq_lim_label, irow, 0)
        grid.addWidget(self.kfreq_lim_edit_min, irow, 1)
        grid.addWidget(self.kfreq_lim_edit_max, irow, 2)
        irow += 1
        #--------------------------------------------------
        # y-axes
        grid.addWidget(self.damp_lim_label, irow, 0)
        grid.addWidget(self.damp_lim_edit_min, irow, 1)
        grid.addWidget(self.damp_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.freq_lim_label, irow, 0)
        grid.addWidget(self.freq_lim_edit_min, irow, 1)
        grid.addWidget(self.freq_lim_edit_max, irow, 2)
        irow += 1
        #--------------------------------------------------

        grid.addWidget(self.eigr_lim_label, irow, 0)
        grid.addWidget(self.eigr_lim_edit_min, irow, 1)
        grid.addWidget(self.eigr_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.eigi_lim_label, irow, 0)
        grid.addWidget(self.eigi_lim_edit_min, irow, 1)
        grid.addWidget(self.eigi_lim_edit_max, irow, 2)
        irow += 1
        #------------------------------------------
        grid.addWidget(self.freq_tol_label, irow, 0)
        grid.addWidget(self.freq_tol_edit, irow, 1)
        irow += 1
        grid.addWidget(self.freq_tol_remove_label, irow, 0)
        grid.addWidget(self.freq_tol_remove_edit, irow, 1)
        irow += 1
        grid.addWidget(self.mode_label, irow, 0)
        grid.addWidget(self.mode_edit, irow, 1)
        irow += 1
        grid.addWidget(self.velocity_label, irow, 0)
        grid.addWidget(self.velocity_edit, irow, 1)
        irow += 1
        grid.addWidget(self.mag_tol_label, irow, 0)
        grid.addWidget(self.mag_tol_edit, irow, 1)
        irow += 1

        grid.addWidget(self.output_directory_label, irow, 0)
        grid.addWidget(self.output_directory_edit, irow, 1)
        grid.addWidget(self.output_directory_browse, irow, 2)
        self.output_directory_label.setDisabled(True)
        self.output_directory_edit.setDisabled(True)
        self.output_directory_browse.setDisabled(True)

        self.output_directory_label.setVisible(False)
        self.output_directory_edit.setVisible(False)
        self.output_directory_browse.setVisible(False)
        irow += 1

        grid.addWidget(self.VL_label, irow, 0)
        grid.addWidget(self.VL_edit, irow, 1)
        irow += 1

        grid.addWidget(self.VF_label, irow, 0)
        grid.addWidget(self.VF_edit, irow, 1)
        irow += 1

        grid.addWidget(self.damping_label, irow, 0)
        grid.addWidget(self.damping_edit, irow, 1)
        irow += 1

        jrow = 0
        grid_check = QGridLayout()
        grid_check.addWidget(self.log_xscale_checkbox, jrow, 0)
        grid_check.addWidget(self.log_yscale1_checkbox, jrow, 1)
        grid_check.addWidget(self.log_yscale2_checkbox, jrow, 2)
        jrow += 1

        grid_check.addWidget(self.show_points_checkbox, jrow, 0)
        grid_check.addWidget(self.show_mode_number_checkbox, jrow, 1)
        jrow += 1
        grid_check.addWidget(self.point_spacing_label, jrow, 0)
        grid_check.addWidget(self.point_spacing_spinner, jrow, 1)
        jrow += 1
        grid_check.addWidget(self.show_lines_checkbox, jrow, 0)
        jrow += 1

        ok_cancel_hbox = QHBoxLayout()
        ok_cancel_hbox.addWidget(self.f06_load_button)
        ok_cancel_hbox.addWidget(self.ok_button)

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
        #log_widget = ApplicationLogWidget(self)

        log_widget = self.setup_logging()
        if self.use_dock_widgets:
            self.modes_dock_widget = NamedDockWidget('Modes', self.modes_widget, self)
            self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.modes_dock_widget)
            self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock_widget)
            vbox2 = vbox
        else:
            #self.log_dock_widget.hide()
            vbox2 = QHBoxLayout()
            vbox2.addWidget(self.modes_widget)
            vbox2.addLayout(vbox)
        widget = QWidget()
        widget.setLayout(vbox2)
        self.setCentralWidget(widget)

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
        #self.bdf_load_button.clicked.connect(self.on_load_bdf)
        #self.op2_load_button.clicked.connect(self.on_load_op2)

        self.x_plot_type_pulldown.currentIndexChanged.connect(self.on_plot_type)
        self.plot_type_pulldown.currentIndexChanged.connect(self.on_plot_type)
        self.subcase_edit.currentIndexChanged.connect(self.on_subcase)
        self.bdf_filename_checkbox.stateChanged.connect(self.on_enable_bdf)
        self.op2_filename_checkbox.stateChanged.connect(self.on_enable_op2)
        self.f06_filename_browse.clicked.connect(self.on_browse_f06)
        self.bdf_filename_browse.clicked.connect(self.on_browse_bdf)
        self.op2_filename_browse.clicked.connect(self.on_browse_op2)
        #self.modes_widget.itemSelectionChanged.connect(self.on_modes)
        # self.modes_widget.itemClicked.connect(self.on_modes)
        # self.modes_widget.currentRowChanged.connect(self.on_modes)
        self.ok_button.clicked.connect(self.on_ok)
        self.units_out_pulldown.currentIndexChanged.connect(self.on_units_out)

        self.pop_vtk_gui_button.clicked.connect(self.on_open_new_window)

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
        self.VL_label.setText(f'VL, Limit ({eas_units}):')
        self.VF_label.setText(f'VF, Flutter ({eas_units}):')

    def on_font_size(self) -> None:
        #font_size = self.font_size_edit.value()
        self.on_set_font_size(self.font_size)

    def on_set_font_size(self, font_size: int) -> None:
        self.font_size = font_size
        font = make_font(font_size, is_bold=False)
        self.setFont(font)

    # @dontcrash
    def on_load_f06(self) -> None:
        f06_filename = os.path.abspath(self.f06_filename_edit.text())
        if not os.path.exists(f06_filename) or not os.path.isfile(f06_filename):
            self.f06_filename_edit.setStyleSheet(QLINEEDIT_RED)
            self.log.error(f"can't find {f06_filename}")
            return
        self.f06_filename_edit.setStyleSheet(QLINEEDIT_WHITE)
        f06_units = self.units_in_pulldown.currentText()
        out_units = self.units_out_pulldown.currentText()

        self.use_rhoref = self.use_rhoref_checkbox.isChecked()
        model, self.responses = load_f06_op2(
            f06_filename, self.log,
            f06_units, out_units,
            self.use_rhoref,
        )

        subcases = list(self.responses.keys())
        if len(subcases) == 0:
            self.log.error('No subcases found')
            return
        #self.log.info(f'on_load_f06: subcases={subcases}')
        self.f06_filename = f06_filename
        self._units_in = f06_units
        self._units_out = out_units
        self.add_recent_file(f06_filename)
        self.update_subcases(subcases)
        self.ok_button.setEnabled(True)

    def add_recent_file(self, f06_filename: str) -> None:
        path = os.path.abspath(f06_filename)
        if path in self.recent_files:
            self.recent_files.remove(path)
        self.recent_files = [path] + self.recent_files
        if len(self.recent_files) > self.nrecent_files_max:
            self.recent_files = self.recent_files[:self.nrecent_files_max]

        self._update_recent_files_actions()
        is_valid = self.validate()
        if is_valid:
            self.log.warning('add_recent_file; save')
            self._save(self.save_filename)

    def _update_recent_files_actions(self) -> None:
        current_directory = os.getcwd()
        for ifile in range(self.nrecent_files_max):
            name = f'file_{ifile}'
            try:
                fname = self.recent_files[ifile]
                show = True
            except IndexError:
                show = False

            action: QAction = self.qactions[name]
            active_text = action.text()

            if show != action.isVisible():
                action.setVisible(show)

            if show and active_text != fname:
                # only update text if it's visible
                relative_fname = os.path.relpath(fname, current_directory)
                text = fname
                if len(relative_fname) < len(fname):
                    text = relative_fname
                action.setText(text)

    def on_subcase(self) -> None:
        subcase, is_subcase_valid = self._get_subcase()
        if not is_subcase_valid:
            return
        response: FlutterResponse = self.responses[subcase]
        #self.log.info(f'on_subcase; response.results.shape={response.results.shape}')
        freqs = response.results[:, 0, response.ifreq].ravel()
        self._update_modal_participation_velocity(response)
        #self.log.info(f'on_subcase; freqs={freqs}')
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
        #self.log.info(f'_get_subcase: subcase_str={subcase_str!r}')
        subcase_sline = subcase_str.split()
        #self.log.info(f'_get_subcase: subcase_sline={subcase_sline}')
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
        #self.log.info(f'update_subcases setting...')
        self.subcase_edit.addItems(subcases_text)

    def update_modes_table(self, modes: list[int],
                           freqs: list[float]) -> None:
        self.modes = modes
        self._set_modes_table(self.modes_widget, modes, freqs)
        self.ok_button.setEnabled(True)
        self.log.info(f'modes = {self.modes}')

    def on_modes(self) -> None:
        self.on_ok()
        # self.validate()
        # self.plot(self.modes)

    # @dontcrash
    def _set_modes_table(self, modes_widget: QListWidget,
                         modes: list[int], freqs: list[float]):
        modes_widget.clear()
        for imode, freq in zip(modes, freqs):
            mode = QListWidgetItem(f'Mode {imode}; f={freq:.2f}')
            # mode.itemClicked.connect(self._on_update_mode)
            mode.setSelected(True)
            modes_widget.addItem(mode)

    def _on_update_mode(self):
        if not self.is_valid:
            #self.log.warning('_on_update_mode')
            self.validate()
        self.plot()

    def on_ok(self) -> None:
        #self.log.warning('on_ok')
        is_valid = self.validate()
        if not is_valid:
            return

        modes = self.selected_modes
        if len(modes) == 0:
            self.log.warning(f'modes = {modes}; assuming all modes -> {self.modes}')
            modes = self.modes
            #return
        self.log.info(f'is_valid = {is_valid}\n')
        self.is_valid = True
        self.plot(modes)
        #self.log.warning('on_ok; _save')
        self._save(self.save_filename)

    # @dontcrash
    def plot(self, modes: list[int]) -> None:
        self.log.info(f'plot; modes = {modes}\n')
        if not self.is_valid:
            self.log.warning(f'not valid\n')
            return
        if len(self.responses) == 0:
            self.log.warning(f'no subcases\n')
            return

        x_plot_type = self.x_plot_type
        plot_type = self.plot_type
        self.log.info(f'plot_type = {plot_type}\n')

        noline = not self.show_lines
        nopoints = not self.show_points

        freq_tol = self.freq_tol
        freq_tol_remove = self.freq_tol_remove
        mag_tol = self.mag_tol
        self.log.info(f'freq_tol = {freq_tol}\n')
        if noline and nopoints:
            noline = False
            nopoints = True

        if x_plot_type == 'eas':
            xlim = self.eas_lim
        elif x_plot_type == 'tas':
            xlim = self.tas_lim
        elif x_plot_type == 'mach':
            xlim = self.mach_lim
        elif x_plot_type == 'alt':
            xlim = self.alt_lim
        elif x_plot_type == 'q':
            xlim = self.q_lim
        else:  # pragma: no cover
            self.log.error(f'x_plot_type={x_plot_type!r} is not supported')
            #raise RuntimeError(x_plot_type)
            xlim = self.xlim

        #self.log.info(f'xlim={xlim}\n')
        assert xlim[0] != '' and xlim[1] != '', (xlim, x_plot_type)
        v_lines = []
        #self.log.info(f'vf={self.vf!r}; vl={self.vl!r}\n')
        if isinstance(self.vf, float) and self.vf > 0.:
            # name, velocity, color, linestyle
            v_lines.append(('VF', self.vf, 'r', '-'))
        # if self.vd:
        #     # name, velocity, color, linestyle
        #     x_limits.append(('VD', self.vd, 'k', '--'))
        #     x_limits.append(('1.15*VD', 1.15*self.vd, 'k', '-'))
        if isinstance(self.vl, float) and self.vl > 0.:
            # name, velocity, color, linestyle
            v_lines.append(('VL', self.vl, 'k', '--'))
            v_lines.append(('1.15*VL', 1.15*self.vl, 'k', '-'))

        #self.log.info(f'v_lines={v_lines}\n')
        #self.log.info(f'kfreq_lim={self.kfreq_lim}\n')
        #self.log.info(f'ydamp_lim={self.ydamp_lim}\n')
        #self.log.info(f'freq_lim={self.freq_lim}\n')
        #self.log.info(f'damping={self.damping}\n')
        xlim_kfreq = self.kfreq_lim
        ylim_damping = self.ydamp_lim
        ylim_freq = self.freq_lim
        damping_limit = self.damping  # % damping

        # changing directory so we don't make a long filename
        # in teh plot header
        #self.log.info(f'damping_limit = {damping_limit}\n')
        dirname = os.path.abspath(os.path.dirname(self.f06_filename))
        basename = os.path.basename(self.f06_filename)

        base = os.path.splitext(basename)[0]
        current_directory = os.getcwd()
        sys.stdout.flush()
        os.chdir(dirname)

        fig = plt.figure(1)
        fig.clear()
        self.log.info(f'cleared plot\n')
        if plot_type not in {'root-locus', 'modal-participation'}:
            gridspeci = gridspec.GridSpec(2, 4)
            damp_axes = fig.add_subplot(gridspeci[0, :3])
            freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        response = self.responses[self.subcase]

        # you can change the output units without reloading
        if self._units_out != self.units_out:
            response.convert_units(self.units_out)
            self._units_out  = self.units_out

        response.noline = noline
        response.set_symbol_settings(
            nopoints, self.show_mode_number, self.point_spacing)
        response.set_font_settings(self.plot_font_size)
        response.log = self.log
        #print('trying plots...')

        self.log.info(f'getting logs\n')
        log_scale_x = self.data['log_scale_x']
        log_scale_y1 = self.data['log_scale_y1']
        log_scale_y2 = self.data['log_scale_y2']
        print(f'log_scale_x={log_scale_x}; log_scale_y1={log_scale_y1}; log_scale_y2={log_scale_y2}')
        print(f'export_to_png={self.export_to_png}')

        self.export_to_png = False
        png_filename0, png_filename = get_png_filename(
            base, x_plot_type, plot_type,
            self.export_to_png)
        print(f'png_filename={png_filename}')
        try:
            if plot_type == 'root-locus':
                axes = fig.add_subplot(111)
                #self.log.info(f'modes={modes}; eigr_lim={self.eigr_lim}; eigi_lim={self.eigi_lim}; freq_tol={freq_tol}')
                #self.log.info(f'png_filename={png_filename}')
                response.plot_root_locus(
                    fig=fig, axes=axes,
                    modes=modes, eigr_lim=self.eigr_lim, eigi_lim=self.eigi_lim,
                    freq_tol=freq_tol,
                    show=True, clear=False, close=False,
                    legend=True,
                    png_filename=png_filename,
                )
            elif plot_type == 'modal-participation':
                axes = fig.add_subplot(111)
                mode = self.mode_edit.value()
                ivel = self.velocity_edit.currentIndex()
                print(f'ivel={ivel}; mode={mode}')
                response.plot_modal_participation(
                    ivel, mode,
                    fig=fig, axes=axes,
                    modes=modes, #eigr_lim=self.eigr_lim, eigi_lim=self.eigi_lim,
                    freq_tol=freq_tol,
                    mag_tol=mag_tol,
                    show=True, clear=False, close=False,
                    legend=True,
                    png_filename=png_filename,
                )
            elif plot_type == 'x-damp-kfreq':
                #xlabel: eas
                #ylabel1 = r'Structural Damping; $g = 2 \gamma $'
                #ylabel2 = r'KFreq [rad]; $ \omega c / (2 V)$'
                #print('plot_kfreq_damping')
                response.plot_kfreq_damping(
                    fig=fig, damp_axes=damp_axes, freq_axes=freq_axes,
                    modes=modes, plot_type=x_plot_type,
                    xlim=xlim, ylim_damping=ylim_damping, ylim_kfreq=xlim_kfreq,
                    freq_tol=freq_tol,
                    show=True, clear=False, close=False,
                    legend=True,
                    png_filename=png_filename,
                )
                update_ylog_style(fig, log_scale_x, log_scale_y1, log_scale_y2)
            else:
                assert plot_type in 'x-damp-freq', plot_type
                #print('plot_vg_vf')
                #self.log.info(f'png_filename={png_filename!r}')
                #self.log.info(f'modes={modes!r}')
                #self.log.info(f'freq_tol={freq_tol!r}')
                #self.log.info(f'v_lines={v_lines!r}')
                response.plot_vg_vf(
                    fig=fig, damp_axes=damp_axes, freq_axes=freq_axes,
                    plot_type=x_plot_type,
                    modes=modes,
                    xlim=xlim, ylim_damping=ylim_damping, ylim_freq=ylim_freq,
                    freq_tol=freq_tol,
                    show=True, clear=False, close=False,
                    legend=True,
                    v_lines=v_lines,
                    #vl_limit=vl,
                    #vd_limit=vd,
                    damping_limit=damping_limit,
                    png_filename=png_filename,
                )
                update_ylog_style(fig, log_scale_x, log_scale_y1, log_scale_y2)
        except Exception as e:  # pragma: no cover
            self.log.error(f'plot_type={plot_type}')
            self.log.error(str(e))
            print(traceback.format_exc())
            #print(traceback.print_tb())
            print(traceback.print_exception(e))
            raise

        base2 = os.path.splitext(png_filename0)[0]
        csv_filename = base2 + '.export.csv'
        veas_filename = base2 + '.export.veas'
        f06_filename = base2 + '.export.f06'
        if self.export_to_csv:
            self.log.debug(f'writing {csv_filename}')
            response.export_to_csv(csv_filename, modes=modes)
        if self.export_to_zona:
            self.log.debug(f'writing {veas_filename}')
            response.export_to_veas(veas_filename, modes=modes, xlim=None)
        if self.export_to_f06:
            self.log.debug(f'writing {f06_filename}')
            response.export_to_f06(f06_filename, modes=modes)
        os.chdir(current_directory)
        self.log.info(f'saved {png_filename}')

    def get_xlim(self) -> tuple[Limit, Limit, Limit, Limit,
                                Limit, Limit, Limit, Limit, Limit,
                                Optional[float], Optional[float],
                                Optional[float], Optional[float], bool]:
        eas_lim_min, is_passed1a = get_float_or_none(self.eas_lim_edit_min)
        eas_lim_max, is_passed1b = get_float_or_none(self.eas_lim_edit_max)
        tas_lim_min, is_passed2a = get_float_or_none(self.tas_lim_edit_min)
        tas_lim_max, is_passed2b = get_float_or_none(self.tas_lim_edit_max)
        mach_lim_min, is_passed3a = get_float_or_none(self.mach_lim_edit_min)
        mach_lim_max, is_passed3b = get_float_or_none(self.mach_lim_edit_max)
        alt_lim_min, is_passed4a = get_float_or_none(self.alt_lim_edit_min)
        alt_lim_max, is_passed4b = get_float_or_none(self.alt_lim_edit_max)
        q_lim_min, is_passed5a = get_float_or_none(self.q_lim_edit_min)
        q_lim_max, is_passed5b = get_float_or_none(self.q_lim_edit_max)
        rho_lim_min, is_passed6a = get_float_or_none(self.rho_lim_edit_min)
        rho_lim_max, is_passed6b = get_float_or_none(self.rho_lim_edit_max)
        xlim_min, is_passed7a = get_float_or_none(self.xlim_edit_min)
        xlim_max, is_passed7b = get_float_or_none(self.xlim_edit_max)

        is_passed_x = all([is_passed1a, is_passed1b,
                           is_passed2a, is_passed2b,
                           is_passed3a, is_passed3b,
                           is_passed4a, is_passed4b,
                           is_passed5a, is_passed5b,
                           is_passed6a, is_passed6b,
                           is_passed7a, is_passed7a,
                           ])

        damp_lim_min, is_passed_damp1 = get_float_or_none(self.damp_lim_edit_min)
        damp_lim_max, is_passed_damp2 = get_float_or_none(self.damp_lim_edit_max)
        freq_lim_min, is_passed_freq1 = get_float_or_none(self.freq_lim_edit_min)
        freq_lim_max, is_passed_freq2 = get_float_or_none(self.freq_lim_edit_max)
        kfreq_lim_min, is_passed_kfreq1 = get_float_or_none(self.kfreq_lim_edit_min)
        kfreq_lim_max, is_passed_kfreq2 = get_float_or_none(self.kfreq_lim_edit_max)

        eigr_lim_min, is_passed_eigr1 = get_float_or_none(self.eigr_lim_edit_min)
        eigr_lim_max, is_passed_eigr2 = get_float_or_none(self.eigr_lim_edit_max)
        eigi_lim_min, is_passed_eigi1 = get_float_or_none(self.eigi_lim_edit_min)
        eigi_lim_max, is_passed_eigi2 = get_float_or_none(self.eigi_lim_edit_max)
        is_passed_eig = all([is_passed_eigr1, is_passed_eigr2,
                             is_passed_eigi1, is_passed_eigi2])

        freq_tol, is_passed_tol1 = get_float_or_none(self.freq_tol_edit)
        freq_tol_remove, is_passed_tol2 = get_float_or_none(self.freq_tol_remove_edit)
        mag_tol, is_passed_tol3 = get_float_or_none(self.mag_tol_edit)
        if is_passed_tol1 and freq_tol is None:
            freq_tol = -1.0
        if is_passed_tol2 and freq_tol_remove is None:
            freq_tol_remove = -1.0
        if is_passed_tol3 and mag_tol is None:
            mag_tol = -1.0

        eas_lim = [eas_lim_min, eas_lim_max]
        tas_lim = [tas_lim_min, tas_lim_max]
        mach_lim = [mach_lim_min, mach_lim_max]
        alt_lim = [alt_lim_min, alt_lim_max]
        q_lim = [q_lim_min, q_lim_max]
        rho_lim = [rho_lim_min, rho_lim_max]
        xlim = [xlim_min, xlim_max]

        vl, is_passed_vl = get_float_or_none(self.VL_edit)
        vf, is_passed_vf = get_float_or_none(self.VF_edit)
        damping, is_passed_damping = get_float_or_none(self.damping_edit)
        if is_passed_vl and vl is None:
            vl = -1.0
        if is_passed_vf and vf is None:
            vf = -1.0
        if is_passed_damping and damping is None:
            damping = -1.0

        damp_lim = [damp_lim_min, damp_lim_max]
        freq_lim = [freq_lim_min, freq_lim_max]
        kfreq_lim = [kfreq_lim_min, kfreq_lim_max]
        eigr_lim = [eigr_lim_min, eigr_lim_max]
        eigi_lim = [eigi_lim_min, eigi_lim_max]
        #self.log.info(f'XLim = {xlim}')
        is_passed_flags = [
            is_passed_x,
            is_passed_damp1, is_passed_damp2,
            is_passed_freq1, is_passed_freq2,
            is_passed_kfreq1, is_passed_kfreq2,
            is_passed_eig,
            is_passed_tol1, is_passed_tol2, is_passed_tol3,
            is_passed_vl, is_passed_vf, is_passed_damping,
        ]
        is_passed = all(is_passed_flags)
        # if not is_passed:
        #self.log.warning(f'is_passed_flags = {is_passed_flags}')
        #print(f'freq_tol = {freq_tol}')
        out = (
            eas_lim, tas_lim, mach_lim, alt_lim, q_lim, rho_lim, xlim,
            damp_lim, freq_lim, kfreq_lim,
            eigr_lim, eigi_lim,
            freq_tol, freq_tol_remove, mag_tol,
            vl, vf, damping, is_passed,
        )
        return out

    def get_selected_modes(self) -> list[int]:
        mode_strs = get_selected_items_flat(self.modes_widget)
        #self.log.info(f'mode_strs = {mode_strs}')
        modes = [int(mode_str.split(';')[0].split(' ')[1])
                 for mode_str in mode_strs]
        self.log.info(f'modes = {modes}')
        return modes

    def validate(self) -> bool:
        #self.log.warning('validate')
        (eas_lim, tas_lim, mach_lim, alt_lim, q_lim, rho_lim, xlim,
         ydamp_lim, freq_lim, kfreq_lim,
         eigr_lim, eigi_lim,
         freq_tol, freq_tol_remove, mag_tol,
         vl, vf, damping, is_valid_xlim) = self.get_xlim()

        selected_modes = []
        subcase, is_subcase_valid = self._get_subcase()
        self.log.warning(f'subcase={subcase}; is_subcase_valid={is_subcase_valid}')
        if is_subcase_valid:
            selected_modes = self.get_selected_modes()

        self.subcase = subcase
        self.selected_modes = selected_modes
        self.eas_lim = eas_lim
        self.tas_lim = tas_lim
        self.mach_lim = mach_lim
        self.alt_lim = alt_lim
        self.q_lim = q_lim
        self.rho_lim = rho_lim
        self.xlim = xlim
        self.ydamp_lim = ydamp_lim
        self.kfreq_lim = kfreq_lim
        self.freq_lim = freq_lim
        self.eigi_lim = eigi_lim
        self.eigr_lim = eigr_lim
        self.freq_tol = freq_tol
        self.freq_tol_remove = freq_tol_remove
        self.mag_tol = mag_tol
        self.vl = vl
        self.vf = vf
        self.damping = damping

        self.x_plot_type = self.x_plot_type_pulldown.currentText()
        self.plot_type = self.plot_type_pulldown.currentText()
        units_in = self.units_in_pulldown.currentText()
        units_out = self.units_out_pulldown.currentText()
        output_directory = self.output_directory_edit.text()

        self.show_lines = self.show_lines_checkbox.isChecked()
        self.show_points = self.show_points_checkbox.isChecked()
        self.show_mode_number = self.show_mode_number_checkbox.isChecked()
        self.point_spacing = self.point_spacing_spinner.value()
        self.use_rhoref = self.use_rhoref_checkbox.isChecked()

        is_passed_modal_partipation = False
        subcases = list(self.responses)
        if len(subcases):
            self.log.info(f'subcases={subcases}')
            subcase0 = subcases[0]
            response = self.responses[subcase0]

            failed_modal_partipation = (
                (self.plot_type == 'modal-participation') and
                ((response.eigr_eigi_velocity is None) or
                 (response.eigenvector is None))
            )
            is_passed_modal_partipation = not failed_modal_partipation
        # (
        #     (self.plot_type == 'modal-participation') and
        #     (response.eigr_eigi_velocity is not None)
        # ) or (self.plot_type != 'modal-participation'))
        data = {
            'bdf_filename': self.bdf_filename,
            'op2_filename': self.op2_filename,
            'log_scale_x': self.log_xscale_checkbox.isChecked(),
            'log_scale_y1': self.log_yscale1_checkbox.isChecked(),
            'log_scale_y2': self.log_yscale2_checkbox.isChecked(),
            'use_rhoref': self.use_rhoref,
            'show_points': self.show_points,
            'show_mode_number': self.show_mode_number,
            'point_spacing': self.point_spacing,
            'show_lines': self.show_lines,
            'export_to_png': self.export_to_png,
            'export_to_csv': self.export_to_csv,
            'export_to_f06': self.export_to_f06,
            'export_to_zona': self.export_to_zona,

            'recent_files': self.recent_files,
            'font_size': self.font_size,
            'plot_font_size': self.plot_font_size,
            'subcase': subcase,
            #'modes': modes,
            'selected_modes': selected_modes,
            'x_plot_type': self.x_plot_type,
            'plot_type': self.plot_type,
            'eas_lim': eas_lim,
            'tas_lim': tas_lim,
            'mach_lim': mach_lim,
            'alt_lim': alt_lim,
            'q_lim': q_lim,
            'rho_lim': rho_lim,
            'xlim': xlim,

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
            'damping': damping,
        }
        self.units_in = units_in
        self.units_out = units_out
        is_passed = all([is_valid_xlim, is_subcase_valid, is_passed_modal_partipation])
        if is_passed:
            self.data = data
            #self.xlim = xlim
            #self.ylim = ydamp_lim
            #self.data = data
            # is_valid = validate_json(self.data, self.log)
            #if is_valid != is_passed:
            #self.log.info(f'passed data:\n{str(self.data)}')
        else:
            del data['recent_files']
            self.log.error(
                f'is_valid_xlim = {is_valid_xlim}\n'
                f'is_subcase_valid = {is_subcase_valid}\n'
                f'is_passed_modal_partipation = {is_passed_modal_partipation}\n'
                f'failed data:\n{str(data)}'
            )
            #self.log.error(f'failed data:\n{str(data)}')
        return is_passed

    def on_open_new_window(self):
        #bdf_filename = self.bdf_filename if not (self.bdf_filename and os.path.exists(self.bdf_filename)) else BDF_FILENAME
        #op2_filename = self.op2_filename if not (self.op2_filename and os.path.exists(self.op2_filename)) else OP2_FILENAME
        bdf_filename = self.bdf_filename if os.path.exists(self.bdf_filename) else ''
        op2_filename = self.op2_filename if os.path.exists(self.op2_filename) else ''
        self._vtk_window_obj.show(bdf_filename, op2_filename)
        return
        try:
            from pyNastran.f06.dev.flutter.gui_flutter_vtk import VtkWindow
        except ImportError as e:
            self.log.error(str(e))
            # print(traceback.print_tb(e))
            print(traceback.print_exception(e))
            self.log.error('cant open window')
            return
        self.new_window = VtkWindow(gui, BDF_FILENAME, OP2_FILENAME)
        self.new_window.show()

    def log_debug(self, msg: str) -> None:
        print(f'DEBUG: {msg}')
    def log_info(self, msg: str) -> None:
        print(f'INFO:  {msg}')
    def log_command(self, msg: str) -> None:
        print(f'COMMAND: {msg}')
    def log_warning(self, msg: str) -> None:
        print(f'WARNING: {msg}')
    def log_error(self, msg: str) -> None:
        print(f'ERROR:   {msg}')

def get_selected_items_flat(list_widget: QListWidget) -> list[str]:
    items = list_widget.selectedItems()
    texts = []
    for item in items:
        text = item.text()
        texts.append(text)
    return texts


def validate_json(data: dict[str, Any],
                  log: SimpleLogger) -> bool:
    is_valid = True
    #log.warning(f'keys = {list(data.keys())}')
    key_allowed_values = [
        ('units_in', UNITS_IN),
        ('units_out', UNITS_OUT),
        ('plot_type', PLOT_TYPES),
    ]
    for (key, allowed_values) in key_allowed_values:
        if key not in data:
            is_valid = False
            log.error(f'data[{key}] is missing; defaulting to {allowed_values[0]}')
            default_value = allowed_values[0]
            data[key] = default_value
            continue

        value = data[key]
        if value not in allowed_values:
            is_valid = False
            log.error(f'{key}={value!r} not in {allowed_values}; defaulting to {allowed_values[0]}')
            default_value = allowed_values[0]
            data[key] = default_value
    return is_valid

def get_float_or_none(line_edit: QLineEdit) -> tuple[Optional[float | str], bool]:
    # is_passed = False
    if not line_edit.isVisible():
        # just echo it back...who cares if it's wrong
        text = line_edit.text().strip()
        #try:
            #value = float(text)
        #except ValueError:
            #value = text
        #value = None
        value = text
        is_passed = True
        return value, is_passed

    text = line_edit.text().strip()
    if len(text) == 0:
        value = None
        is_passed = True
    else:
        try:
            value = float(text)
            is_passed = True
        except ValueError:
            value = None
            is_passed = False
    return value, is_passed


def _to_str(value: Optional[int | float]) -> str:
    if value is None:
        str_value = ''
    else:
        str_value = str(value)
    return str_value

def main(f06_filename: str='') -> None:  # pragma: no cover
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    # The Main window

    main_window = FlutterGui(f06_filename)
    # Enter the main loop
    app.exec_()


if __name__ == '__main__':
    main()
