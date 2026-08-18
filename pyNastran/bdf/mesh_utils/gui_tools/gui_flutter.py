# encoding: utf-8
import json
import os
import sys
import signal
from pathlib import Path

import numpy as np
from cpylog import SimpleLogger
from cpylog.html_utils import str_to_html

#from qtpy import QtGui
from qtpy import QtCore
from qtpy.QtWidgets import (
    QLabel, QPushButton, QGridLayout, QApplication,
    QSpinBox, QLineEdit, QCheckBox,
    QWidget, QComboBox,
    QRadioButton, QButtonGroup,
    QFileDialog,
    QHBoxLayout, QVBoxLayout,)

from pyNastran.bdf.mesh_utils.cmd_line.create_flutter import (
    create_flutter, VELOCITY_UNITS, ALT_UNITS)
from pyNastran.utils.atmosphere import atm_density

from pyNastran.gui.menus.cutting_plane.results_dialog import ResultsDialog
from pyNastran.gui.utils.qt.checks.qlineedit import check_float

#import pyNastran
from pyNastran.gui.menus.application_log import ApplicationLogWidget
from pyNastran.gui.utils.qt.pydialog import QFloatEdit, QIntEdit, PyDialog

from pyNastran.gui.utils.qt.qcombobox import get_combo_box_text
from pyNastran.utils.convert import convert_velocity

# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
signal.signal(signal.SIGINT, signal.SIG_DFL)

USE_WIN = False

UNIT_SYSTEMS_MAP = {
    'in-slinch-s (English-in)': 'english_in',
    'ft-slug-s (English-ft)': 'english_ft',
    'm-kg-s (SI)': 'si',
    'mm-Mg-s (SI-mm)': 'si_mm',
}
UNIT_SYSTEMS = list(UNIT_SYSTEMS_MAP.keys())

SETTINGS_PATH = Path.home() / 'bdf_flutter.json'

CONSTANT_TYPE_MAP = {
    'Mach': 'mach',
    'Altitude': 'alt',
    'Velocity': 'tas',
    'Equivalent Airspeed': 'eas',
}
# same as const formats
SWEEP_FORMATS = list(CONSTANT_TYPE_MAP.keys())


class FlutterGui(PyDialog):
    def __init__(self, data, win_parent=None):
        """
        Saves the data members from data and
        performs type checks
        """
        PyDialog.__init__(self, data, win_parent)
        self.log = SimpleLogger(level='debug', encoding='utf-8')
        self._updated_preference = False
        self.dim_max = 10

        self.setWindowTitle('Flutter Cards')
        self.create_widgets()

        self.log = None

        self.html_logging = True
        self._start_logging()
        #if self.html_logging is True:
        self.log_dock_widget = ApplicationLogWidget(self)
        self.log_widget = self.log_dock_widget.log_widget
        #self.win.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock_widget)

        self.create_layout()
        self.set_connections()
        self._load_settings()

        # self.on_font(self.font_size)
        # self.show()

    def on_sweep_pulldown(self) -> None:
        sweep = get_combo_box_text(self.sweep_pulldown)
        constant = get_combo_box_text(self.constant_pulldown)
        is_eas_enabled = (sweep != 'Equivalent Airspeed') and (constant != 'Equivalent Airspeed')
        self.eas_limit_edit.setEnabled(is_eas_enabled)
        self.eas_limit_unit_pulldown.setEnabled(is_eas_enabled)

        if sweep == 'Mach':
            self.sweep_unit_pulldown.clear()
            self.sweep_unit_pulldown.addItems(['-'])
            self.sweep_unit_pulldown.setItemText(0, '-')
            self.sweep_unit_pulldown.setEnabled(False)
        elif sweep in {'Equivalent Airspeed', 'Velocity'}:
            self.sweep_unit_pulldown.setEnabled(True)
            self.sweep_unit_pulldown.clear()
            self.sweep_unit_pulldown.addItems(VELOCITY_UNITS)
            self.sweep_unit_pulldown.setItemText(0, VELOCITY_UNITS[0])
        elif sweep == 'Altitude':
            self.sweep_unit_pulldown.setEnabled(True)
            self.sweep_unit_pulldown.clear()
            self.sweep_unit_pulldown.addItems(ALT_UNITS)
            self.sweep_unit_pulldown.setItemText(0, ALT_UNITS[0])
        else:  # pragma: no cover
            raise NotImplementedError(sweep)

    def on_constant_pulldown(self) -> None:
        sweep = get_combo_box_text(self.sweep_pulldown)
        constant = get_combo_box_text(self.constant_pulldown)
        is_eas_enabled = (sweep != 'Equivalent Airspeed') and (constant != 'Equivalent Airspeed')
        self.eas_limit_edit.setEnabled(is_eas_enabled)
        self.eas_limit_unit_pulldown.setEnabled(is_eas_enabled)

        # SWEEP_FORMATS = ['Mach', 'Equivalent Airspeed', 'Velocity', 'Altitude']
        if constant == 'Mach':
            self.constant_unit_pulldown.clear()
            self.constant_unit_pulldown.addItems(['-'])
            self.constant_unit_pulldown.setItemText(0, '-')
            self.constant_unit_pulldown.setEnabled(False)
        elif constant in {'Equivalent Airspeed', 'Velocity'}:
            self.constant_unit_pulldown.setEnabled(True)
            self.constant_unit_pulldown.clear()
            self.constant_unit_pulldown.addItems(VELOCITY_UNITS)
            self.constant_unit_pulldown.setItemText(0, VELOCITY_UNITS[0])
        elif constant == 'Altitude':
            self.constant_unit_pulldown.setEnabled(True)
            self.constant_unit_pulldown.clear()
            self.constant_unit_pulldown.addItems(ALT_UNITS)
            self.constant_unit_pulldown.setItemText(0, ALT_UNITS[0])
        else:  # pragma: no cover
            raise NotImplementedError(constant)

    def on_apply(self) -> None:
        sweep_method_unmapped = self.sweep_pulldown.currentText()
        constant_type_unmapped = self.constant_pulldown.currentText()
        unit_system_unmapped = self.unit_system_pulldown.currentText()
        constant_unit = self.constant_unit_pulldown.currentText()
        sweep_unit = self.sweep_unit_pulldown.currentText()
        eas_units = self.eas_limit_unit_pulldown.currentText()
        flutter_id = int(self.flutter_id_edit.text())
        npoints = self.n_value.value()

        # log = self.log
        log = SimpleLogger(level='debug')
        log.info(f'flutter_id = {flutter_id}')

        sweep_method = CONSTANT_TYPE_MAP[sweep_method_unmapped]
        constant_type = CONSTANT_TYPE_MAP[constant_type_unmapped]
        units_out = UNIT_SYSTEMS_MAP[unit_system_unmapped]

        # if sweep_method == 'mach':
        #     sweep_unit = ''
        value1, value1_flag = check_float(self.sweep1_value)
        value2, value2_flag = check_float(self.sweep2_value)
        const_value, const_value_flag = check_float(self.constant_value_edit)

        is_large = self.large_field_checkbox.isChecked()
        clean = self.clean_checkbox.isChecked()
        size = 16 if is_large else 8
        output_dir = self.output_dir_edit.text().strip()
        flutter_filename = self.aero_filename.text().strip()
        if output_dir:
            bdf_filename_out = os.path.join(output_dir, flutter_filename)
        else:
            bdf_filename_out = os.path.abspath(flutter_filename)

        # bdf flutter UNITS eas  EAS1  EAS2  SWEEP_UNIT N CONST_TYPE CONST_VAL CONST_UNIT
        # bdf flutter UNITS tas  TAS1  TAS2  SWEEP_UNIT N CONST_TYPE CONST_VAL CONST_UNIT [--eas_limit EAS EAS_UNITS]
        # bdf flutter UNITS alt  ALT1  ALT2  SWEEP_UNIT N CONST_TYPE CONST_VAL CONST_UNIT [--eas_limit EAS EAS_UNITS]
        # bdf flutter UNITS mach MACH1 MACH2            N CONST_TYPE CONST_VAL CONST_UNIT [--eas_limit EAS EAS_UNITS]
        # ...  [-o OUT_BDF_FILENAME] [--size SIZE | --clean]

        # the defaults for create_flutter
        eas_limit = 1_000_000
        eas_str = ''
        if (sweep_method != 'eas') and (constant_type != 'eas'):
            eas_limit, eas_limit_flag = check_float(self.eas_limit_edit)
            eas_str = f'--eas_limit {eas_limit} {eas_units}'
            is_passed = value1_flag and value2_flag and const_value_flag and eas_limit_flag
        else:
            is_passed = value1_flag and value2_flag and const_value_flag

        if sweep_method == constant_type:
            log.error('sweep_method=constant_type Error\n'
                     f'sweep_method={sweep_method} constant_type={constant_type} and must be different')
            # self.sweep_pulldown.setColor
            # self.constant_unit_pulldown.setColor
            if is_passed:
                return

        is_zaero = self._radio_zaero.isChecked()
        if not is_passed:
            log.error('Invalid parsing')
            return

        size_str = '--clean' if clean else f'--size {size}'
        sweep_unit2 = ''
        if sweep_unit != '-':
            sweep_unit2 = f' {sweep_unit}'

        cmd = (
            f'bdf flutter {units_out} '
            f'{sweep_method} {value1} {value2}{sweep_unit2} {npoints} '
            f'{constant_type} {const_value} {constant_unit} {eas_str} {size_str}'
        ).rstrip()
        if bdf_filename_out:
            cmd += f' --output {bdf_filename_out!r}'
        log.info(cmd)

        if constant_unit == '-':
            constant_unit = 'none'

        try:
            model, density_units, velocity_units = create_flutter(
                log,
                sweep_method, value1, value2, sweep_unit, npoints,
                constant_type, const_value, constant_unit,
                eas_limit=eas_limit, eas_units=eas_units,
                units_out=units_out, sid=flutter_id,
                size=size, clean=clean,
                bdf_filename_out=bdf_filename_out,
                comment=cmd, is_zaero=is_zaero)
        except Exception as error:
            log.error(str(error))
            # raise
            return

        self._save_settings()

        if is_zaero:
            atmos_id = flutter_id
            fix_id = flutter_id + 1
            if sweep_method == 'eas' and constant_type == 'mach':
                log.info(f'atmos = {list(model.zaero.atmos)}')
                atmos = model.zaero.atmos[atmos_id]
                alt = atmos.alt
                nalt = len(alt)
                density = atmos.density
                sos = atmos.sos
                log.info(f'flutter_table = {list(model.zaero.flutter_table)}')
                table = model.zaero.flutter_table[fix_id]
                if table.type == 'FIXMATM':
                    # alts = table.alts
                    mkaeroz_id = table.mkaeroz_id
                    log.info(f'mkaeroz = {list(model.zaero.mkaeroz)}')
                    mkaeroz = model.zaero.mkaeroz[mkaeroz_id]
                    machi = mkaeroz.mach
                    mach = machi * np.ones(nalt)
                    velocity = mach * sos
                else:
                    raise NotImplementedError((sweep_method, constant_type, table.type))

                density0 = atm_density(alt=0.0, density_units=density_units)
                eas_scale = convert_velocity(1., velocity_units, eas_units)
                eas = velocity * np.sqrt(density / density0) * eas_scale
                data = np.column_stack([alt, density, sos.round(1), mach, velocity.round(1), eas.round(1)])
                alt_units = 'ft'
                labels = [
                    f'Altitude ({alt_units})',
                    f'Density ({density_units})',
                    f'SOS ({velocity_units})',
                    'Mach',
                    f'Velocity ({velocity_units})',
                    f'EAS ({eas_units})',
                ]
            else:
                return
        else:
            # flutter_id = 1
            flfact_rho = flutter_id + 1
            flfact_mach = flutter_id + 2
            flfact_velocity = flutter_id + 3
            flfact_eas = flutter_id + 4
            rho = model.flfacts[flfact_rho].factors
            mach = model.flfacts[flfact_mach].factors
            velocity = model.flfacts[flfact_velocity].factors
            eas = model.flfacts[flfact_eas].factors
            data = np.column_stack([rho, mach, velocity, eas])
            labels = [
                f'Density ({density_units})',
                'Mach',
                f'Velocity ({velocity_units})',
                f'EAS ({eas_units})',
            ]

        dlg = ResultsDialog(self, data, labels, title='Atmosphere Table')
        dlg.show()

    def create_widgets(self):
        """creates the display window"""
        #self.win = QMainWindow(self)
        # window text size
        # self.font_size_label = QLabel('Font Size:')
        # self.font_size_edit = QSpinBox(self)

        # self.font_size_edit.setValue(self._default_font_size)
        # self.font_size_edit.setRange(FONT_SIZE_MIN, FONT_SIZE_MAX)
        self.sweep_label = QLabel('Sweep Method:')
        self.sweep_pulldown = QComboBox(self)
        self.sweep_pulldown.addItems(SWEEP_FORMATS)

        self.sweep1_label = QLabel('Value 1:', self)
        self.sweep1_value = QFloatEdit('0.1')
        self.sweep1_value.setToolTip('Starting Value')

        self.sweep2_label = QLabel('Value 2:', self)
        self.sweep2_value = QFloatEdit('0.99')
        self.sweep2_value.setToolTip('Ending Value')

        self.unit_system_label = QLabel('Output Unit System:', self)
        self.unit_system_pulldown = QComboBox(self)
        self.unit_system_pulldown.addItems(UNIT_SYSTEMS)
        self.unit_system_pulldown.setItemText(1, UNIT_SYSTEMS[1])

        self.sweep_unit_label = QLabel('Unit:', self)
        self.sweep_unit_pulldown = QComboBox(self)
        self.sweep_unit_pulldown.setToolTip('Units')
        self.sweep_unit_pulldown.addItems(['-'])
        self.sweep_unit_pulldown.setItemText(0, '-')

        self.n_label = QLabel('Number of Points:', self)
        self.n_value = QSpinBox(self)
        self.n_value.setValue(20)
        self.n_value.setMinimum(2)
        self.n_value.setMaximum(1001)
        self.n_value.setToolTip('Number of Points')

        self.constant_label = QLabel('Constant Type:')
        self.constant_pulldown = QComboBox(self)
        self.constant_pulldown.addItems(SWEEP_FORMATS)
        self.constant_pulldown.setItemText(1, SWEEP_FORMATS[1])

        self.constant_value_label = QLabel('Constant Value:', self)
        self.constant_value_edit = QFloatEdit('0')
        self.constant_value_edit.setToolTip('The constant value')

        # self.constant_unit_label = QLabel('Constant Value:', self)
        self.constant_unit_pulldown = QComboBox(self)

        self.output_dir_label = QLabel('Output Directory:')
        self.output_dir_edit = QLineEdit(self)
        self.output_dir_edit.setText(os.getcwd())
        self.output_dir_edit.setToolTip('Directory for the output flutter file')
        self.output_dir_browse = QPushButton('Browse...')

        self.aero_filename_label = QLabel('Flutter File:')
        self.aero_filename = QLineEdit(self)
        self.aero_filename.setText('flutter_cards.bdf')
        self.aero_filename.setToolTip('Path to the Flutter File')
        self.aero_filename_load = QPushButton('Load...')
        self.aero_filename_load.setEnabled(False)

        self.flutter_id_label = QLabel('Flutter ID:')
        self.flutter_id_edit = QIntEdit('10')
        self.flutter_id_edit.setToolTip('ID of the FLUTTER card')
        self.large_field_checkbox = QCheckBox('Large Field')
        self.clean_checkbox = QCheckBox('Clean')
        # ------------------------------------------------------------------
        self.eas_limit_label = QLabel('EAS Limit:', self)
        self.eas_limit_edit = QFloatEdit('1000')
        self.eas_limit_edit.setToolTip('Equivalent Airspeed Limit; V_EAS = V_TAS * sqrt(rho/rho0)')
        self.eas_limit_unit_pulldown = QComboBox(self)
        self.eas_limit_unit_pulldown.addItems(VELOCITY_UNITS)
        self.eas_limit_unit_pulldown.setItemText(0, VELOCITY_UNITS[0])

        self._radio_nastran = QRadioButton('Nastran')
        self._radio_zaero = QRadioButton('ZAero')
        plot_type_group = QButtonGroup(self)
        plot_type_group.addButton(self._radio_nastran)
        plot_type_group.addButton(self._radio_zaero)
        self._radio_nastran.setChecked(True)

        # ------------------------------------------------------------------
        # closing
        self.apply_button = QPushButton('Apply')
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Exit')
        self.on_sweep_pulldown()
        self.on_constant_pulldown()

    def create_layout(self):
        agrid = self._create_aero_grid()
        #
        awidget = QWidget(self)
        awidget.setLayout(agrid)

        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)
        ok_widget = QWidget(self)
        ok_widget.setLayout(ok_cancel_box)

        # ------------------------------
        vbox = QVBoxLayout(self)
        vbox.addWidget(awidget)

        vbox.addStretch()
        # vbox.addLayout(ok_cancel_box)
        vbox.addWidget(ok_widget)
        if USE_WIN:
            self.win.setLayout(vbox)
            main_vbox = QVBoxLayout(self)
            main_vbox.addWidget(self.win)
            self.setLayout(main_vbox)
        else:
            vbox.addWidget(self.log_dock_widget)
            self.setLayout(vbox)

    def _create_aero_grid(self):
        grid = QGridLayout(self)
        irow = 0
        if 0:
            grid.addWidget(self.sweep_label, irow, 0)
            grid.addWidget(self.sweep_pulldown, irow, 1)
            # grid.addWidget(self.font_size_edit, irow, 1)
            irow += 1

            grid.addWidget(self.sweep1_label, irow, 0)
            grid.addWidget(self.sweep1_value, irow, 1)
            irow += 1

            grid.addWidget(self.sweep2_label, irow, 0)
            grid.addWidget(self.sweep2_value, irow, 1)
            irow += 1

            grid.addWidget(self.n_label, irow, 0)
            grid.addWidget(self.n_value, irow, 1)
            irow += 1

            grid.addWidget(self.constant_label, irow, 0)
            grid.addWidget(self.constant_pulldown, irow, 1)
            grid.addWidget(self.constant_value_edit, irow, 2)
            grid.addWidget(self.constant_unit_pulldown, irow, 3)
            irow += 1
        else:
            grid.addWidget(self.sweep_label, irow, 0)
            grid.addWidget(self.sweep_unit_label, irow, 1)
            grid.addWidget(self.sweep1_label, irow, 2)
            grid.addWidget(self.sweep2_label, irow, 3)
            grid.addWidget(self.n_label, irow, 4)
            irow += 1

            grid.addWidget(self.sweep_pulldown, irow, 0)
            grid.addWidget(self.sweep_unit_pulldown, irow, 1)
            grid.addWidget(self.sweep1_value, irow, 2)
            grid.addWidget(self.sweep2_value, irow, 3)
            grid.addWidget(self.n_value, irow, 4)
            irow += 1

            unit_label = QLabel('Unit:', self)
            grid.addWidget(self.constant_label, irow, 0)
            grid.addWidget(self.constant_value_label, irow, 1)
            grid.addWidget(unit_label, irow, 2)
            irow += 1
            grid.addWidget(self.constant_pulldown, irow, 0)
            grid.addWidget(self.constant_value_edit, irow, 1)
            grid.addWidget(self.constant_unit_pulldown, irow, 2)
            irow += 1

            grid.addWidget(self.eas_limit_label, irow, 0)
            grid.addWidget(self.eas_limit_edit, irow, 1)
            grid.addWidget(self.eas_limit_unit_pulldown, irow, 2)
            irow += 1
            #-----------------------------------------------------
            grid.addWidget(self.flutter_id_label, irow, 0)
            grid.addWidget(self.flutter_id_edit, irow, 1)
            irow += 1

            grid.addWidget(self.large_field_checkbox, irow, 0)
            grid.addWidget(self.clean_checkbox, irow, 1)
            irow += 1

            grid.addWidget(self.unit_system_label, irow, 0)
            grid.addWidget(self.unit_system_pulldown, irow, 1)
            irow += 1

            grid.addWidget(self._radio_nastran, irow, 0)
            grid.addWidget(self._radio_zaero, irow, 1)
            irow += 1

            grid.addWidget(self.output_dir_label, irow, 0)
            grid.addWidget(self.output_dir_edit, irow, 1)
            grid.addWidget(self.output_dir_browse, irow, 2)
            irow += 1

            grid.addWidget(self.aero_filename_label, irow, 0)
            grid.addWidget(self.aero_filename, irow, 1)
            grid.addWidget(self.aero_filename_load, irow, 2)
            irow += 1

        return grid

    def set_connections(self):
        self.sweep_pulldown.currentIndexChanged.connect(self.on_sweep_pulldown)
        self.constant_pulldown.currentIndexChanged.connect(self.on_constant_pulldown)
        self.output_dir_browse.clicked.connect(self.on_browse_output_dir)
        self.apply_button.clicked.connect(self.on_apply)

    def on_browse_output_dir(self) -> None:
        current_dir = self.output_dir_edit.text().strip()
        if not os.path.isdir(current_dir):
            current_dir = os.getcwd()
        dirname = QFileDialog.getExistingDirectory(
            self, 'Select Output Directory', current_dir)
        if dirname:
            self.output_dir_edit.setText(dirname)

    def _save_settings(self) -> None:
        settings = {
            'sweep_method': self.sweep_pulldown.currentText(),
            'sweep_unit': self.sweep_unit_pulldown.currentText(),
            'sweep_value1': self.sweep1_value.text(),
            'sweep_value2': self.sweep2_value.text(),
            'npoints': self.n_value.value(),
            'constant_type': self.constant_pulldown.currentText(),
            'constant_value': self.constant_value_edit.text(),
            'constant_unit': self.constant_unit_pulldown.currentText(),
            'eas_limit': self.eas_limit_edit.text(),
            'eas_limit_unit': self.eas_limit_unit_pulldown.currentText(),
            'unit_system': self.unit_system_pulldown.currentText(),
            'large_field': self.large_field_checkbox.isChecked(),
            'clean': self.clean_checkbox.isChecked(),
            'nastran': self._radio_nastran.isChecked(),
            'output_dir': self.output_dir_edit.text(),
            'flutter_filename': self.aero_filename.text(),
            'flutter_id': self.flutter_id_edit.text(),
        }
        try:
            with open(SETTINGS_PATH, 'w') as json_file:
                json.dump(settings, json_file, indent=2)
            self.log.info(f'Settings saved to {SETTINGS_PATH}')
        except OSError as error:
            self.log.warning(f'Could not save settings: {error}')

    def _load_settings(self) -> None:
        if not SETTINGS_PATH.is_file():
            return
        try:
            with open(SETTINGS_PATH, 'r') as f:
                settings = json.load(f)
        except (OSError, json.JSONDecodeError):
            return

        def _set_combo(combo: QComboBox, value: str) -> None:
            idx = combo.findText(value)
            if idx >= 0:
                combo.setCurrentIndex(idx)

        combos = [
            (self.sweep_pulldown, 'sweep_method', self.on_sweep_pulldown),
            (self.sweep_unit_pulldown, 'sweep_unit', None),
        ]
        log = self.log
        for obj, name, action in combos:
            if name in settings:
                value = settings[name]
                _set_combo(obj, value)
                if action is not None:
                    action()
            else:
                log.warning(f'Could not find {name!r} in settings')

        texts = [
            (self.sweep1_value, 'sweep_value1'),
            (self.sweep2_value, 'sweep_value2'),
            (self.eas_limit_edit, 'eas_limit'),
            (self.output_dir_edit, 'output_dir'),
            (self.aero_filename, 'flutter_filename'),
            (self.constant_value_edit, 'constant_value'),
            (self.flutter_id_edit, 'flutter_id'),
        ]
        for obj, name in texts:
            if name in settings:
                value = settings[name]
                obj.setText(value)
            else:
                log.warning(f'Could not find {name!r} in settings')

        checkboxes = [
            (self.large_field_checkbox, 'large_field'),
            (self.clean_checkbox, 'clean'),
        ]
        for obj, name in checkboxes:
            if name in settings:
                value = settings[name]
                obj.setChecked(value)
            else:
                log.warning(f'Could not find {name!r} in settings')

        if 'npoints' in settings:
            self.n_value.setValue(settings['npoints'])
        if 'constant_type' in settings:
            _set_combo(self.constant_pulldown, settings['constant_type'])
            self.on_constant_pulldown()
        if 'constant_unit' in settings:
            _set_combo(self.constant_unit_pulldown, settings['constant_unit'])
        if 'eas_limit_unit' in settings:
            _set_combo(self.eas_limit_unit_pulldown, settings['eas_limit_unit'])
        if 'unit_system' in settings:
            _set_combo(self.unit_system_pulldown, settings['unit_system'])
        if 'nastran' in settings:
            if settings['nastran']:
                self._radio_nastran.setChecked(True)
            else:
                self._radio_zaero.setChecked(True)

    def _start_logging(self) -> None:
        if self.log is not None:
            return
        if self.html_logging is True:
            log = SimpleLogger(
                level='debug', encoding='utf-8',
                log_func=lambda w, x, y, z: self._logg_msg(w, x, y, z))
            # logging needs synchronizing, so the messages from different
            # threads would not be interleave
            self.log_mutex = QtCore.QReadWriteLock()
        else:
            log = SimpleLogger(
                level='debug', encoding='utf-8',
                #log_func=lambda x, y: print(x, y)  # no colorama
            )
        self.log = log

    def _logg_msg(self, log_type: str, filename: str, lineno: int, msg: str) -> None:
        """
        Add message to log widget trying to choose right color for it.

        Parameters
        ----------
        log_type : str
            {DEBUG, INFO, ERROR, COMMAND, WARNING} or prepend 'GUI '
        filename : str
            the active file
        lineno : int
            line number
        msg : str
            message to be displayed
        """
        if not self.html_logging:
            # standard logger
            name = '%-8s' % (log_type + ':')
            filename_n = '%s:%s' % (filename, lineno)
            msg2 = ' %-28s %s\n' % (filename_n, msg)
            print(name, msg2)
            return

        # if 'DEBUG' in log_type and not self.settings.show_debug:
        #     return
        # elif 'INFO' in log_type and not self.settings.show_info:
        #     return
        # elif 'COMMAND' in log_type and not self.settings.show_command:
        #     return
        # elif 'WARNING' in log_type and not self.settings.show_warning:
        #     return
        # elif 'ERROR' in log_type and not self.settings.show_error:
        #     return

        if log_type in ['GUI ERROR', 'GUI COMMAND', 'GUI DEBUG', 'GUI INFO', 'GUI WARNING']:
            log_type = log_type[4:] # drop the GUI

        html_msg = str_to_html(log_type, filename, lineno, msg)

        self._log_msg(html_msg)

    def _log_msg(self, msg: str) -> None:
        """prints an HTML log message"""
        self.log_mutex.lockForWrite()
        text_cursor = self.log_widget.textCursor()
        end = text_cursor.End
        text_cursor.movePosition(end)
        text_cursor.insertHtml(msg)
        self.log_widget.ensureCursorVisible() # new message will be visible
        self.log_mutex.unlock()



def cmd_line_gui():
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)

    data = {
        'font_size': 8,
    }
    # The Main window
    main_window = FlutterGui(data)
    main_window.show()
    # Enter the main loop
    app.exec_()
    return data

if __name__ == '__main__':  # pragma: no cover
    cmd_line_gui()
