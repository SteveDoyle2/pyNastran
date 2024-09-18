# encoding: utf-8
import os
from typing import Any

from PIL.ImageChops import constant
import numpy as np
from cpylog import SimpleLogger
#from qtpy import QtGui
from qtpy.QtWidgets import (
    QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
    QSpinBox, QDoubleSpinBox, QColorDialog, QLineEdit, QCheckBox,
    QTabWidget, QWidget, QComboBox, QHBoxLayout, QVBoxLayout,
)
from pyNastran.gui.menus.cutting_plane.results_dialog import ResultsDialog
from pyNastran.gui.utils.qt.checks.qlineedit import check_save_path, check_float, QLINEEDIT_GOOD

#import pyNastran
#from pyNastran.gui.utils.qt.pydialog import PyDialog, make_font, check_color
from pyNastran.gui.utils.qt.pydialog import QFloatEdit, QIntEdit
#from pyNastran.gui.utils.qt.checks.qlineedit import (
#    check_int, check_float, check_name_str, check_path, QLINEEDIT_GOOD, QLINEEDIT_ERROR)
#from pyNastran.utils import print_bad_path
#from pyNastran.converters.cart3d.cart3d import read_cart3d
#from pyNastran.converters.tecplot.tecplot import read_tecplot
#from pyNastran.dev.tools.pressure_map_aero_setup import get_aero_model

#import sys
#from typing import Any
#from cpylog import SimpleLogger
#from .utils import filter_no_args
#from pyNastran.utils.convert import convert_altitude, convert_velocity

from pyNastran.gui.utils.qt.pydialog import PyDialog, make_font, check_color
from pyNastran.gui.utils.qt.qcombobox import get_combo_box_text
from pyNastran.bdf.mesh_utils.cmd_line.create_flutter import create_flutter
from pyNastran.bdf.bdf import BDF


#from pyNastran.bdf.mesh_utils.cmd_line.create_flutter import VELOCITY_UNITS, ALTITUDE_UNITS
VELOCITY_UNITS = ['knots', 'ft/s', 'in/s', 'm/s', 'cm/s', 'mm/s']
ALTITUDE_UNITS = ['ft', 'm']
UNIT_SYSTEMS_MAP = {
    'in-slinch-s (English-in)': 'english_in',
    'ft-slug-s (English)' : 'english',
    'm-kg-s (SI)' : 'si',
    'mm-Mg-s (SI-mm)' : 'si_mm',
}
UNIT_SYSTEMS = list(UNIT_SYSTEMS_MAP.keys())[

CONSTANT_TYPE_MAP = {
    'Mach': 'mach',
    'Altitude': 'alt',
    'Velocity': 'vel',
    'Equivalent Airspeed' : 'eas',
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
        self._updated_preference = False
        self.dim_max = 8

        self.setWindowTitle('Flutter Gui')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        # self.on_font(self.font_size)
        # self.show()

    def on_sweep_pulldown(self) -> None:
        sweep = get_combo_box_text(self.sweep_pulldown)
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
            self.sweep_unit_pulldown.addItems(ALTITUDE_UNITS)
            self.sweep_unit_pulldown.setItemText(0, ALTITUDE_UNITS[0])
        else:  # pragma: no cover
            raise NotImplementedError(sweep)

    def on_constant_pulldown(self) -> None:
        constant = get_combo_box_text(self.constant_pulldown)

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
            self.constant_unit_pulldown.addItems(ALTITUDE_UNITS)
            self.constant_unit_pulldown.setItemText(0, ALTITUDE_UNITS[0])
        else:  # pragma: no cover
            raise NotImplementedError(constant)

    def on_apply(self) -> None:
        sweep_method_unmapped = self.sweep_pulldown.currentText()
        constant_type_unmapped = self.constant_pulldown.currentText()
        unit_system_unmapped = self.unit_system_pulldown.currentText()
        constant_unit = self.constant_unit_pulldown.currentText()
        sweep_unit = self.sweep_unit_pulldown.currentText()
        npoints = self.n_value.value()

        sweep_method = CONSTANT_TYPE_MAP[sweep_method_unmapped]
        constant_type = CONSTANT_TYPE_MAP[constant_type_unmapped]
        units_out = UNIT_SYSTEMS_MAP[unit_system_unmapped]

        value1, value1_flag = check_float(self.sweep1_value)
        value2, value2_flag = check_float(self.sweep2_value)
        const_value, const_value_flag = check_float(self.constant_value)
        eas_limit = 1000
        eas_units = 'm/s'
        log = SimpleLogger(level='debug', encoding='utf-8')
        size = 8
        clean = False
        bdf_filename_out = ''
        cmd = ''

        model, density_units, velocity_units = create_flutter(
            log,
            sweep_method, value1, value2, sweep_unit, npoints,
            constant_type, const_value, constant_unit,
            eas_limit=eas_limit, eas_units=eas_units,
            units_out=units_out,
            size=size, clean=clean,
            bdf_filename_out=bdf_filename_out,
            comment=cmd)
        sid = 1
        flfact_rho = sid + 1
        flfact_mach = sid + 2
        flfact_velocity = sid + 3
        flfact_eas = sid + 4
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
        dlg = ResultsDialog(self, data, labels,
                            title='Atmsospere Table')
        dlg.show()

    def create_widgets(self):
        """creates the display window"""
        # window text size
        # self.font_size_label = QLabel('Font Size:')
        # self.font_size_edit = QSpinBox(self)

        # self.font_size_edit.setValue(self._default_font_size)
        # self.font_size_edit.setRange(FONT_SIZE_MIN, FONT_SIZE_MAX)
        self.sweep_label = QLabel('Sweep:')
        self.sweep_pulldown = QComboBox(self)
        self.sweep_pulldown.addItems(SWEEP_FORMATS)

        self.sweep1_label = QLabel('Value 1:', self)
        self.sweep1_value = QFloatEdit('0.1')
        self.sweep1_value.setToolTip('Starting Value')

        self.sweep2_label = QLabel('Value 2:', self)
        self.sweep2_value = QFloatEdit('0.99')
        self.sweep2_value.setToolTip('Ending Value')

        self.unit_system_label = QLabel('Unit System:', self)
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
        self.n_value.setToolTip('Number of Points')

        self.constant_label = QLabel('Sweep:')
        self.constant_pulldown = QComboBox(self)
        self.constant_pulldown.addItems(SWEEP_FORMATS)
        self.constant_pulldown.setItemText(1, SWEEP_FORMATS[1])

        self.constant_value_label = QLabel('Constant Value:', self)
        self.constant_value = QFloatEdit('0')
        self.constant_value.setToolTip('The constant value')

        # self.constant_unit_label = QLabel('Constant Value:', self)
        self.constant_unit_pulldown = QComboBox(self)

        self.aero_filename_label = QLabel('Flutter File:')
        self.aero_filename = QLineEdit(self)
        self.aero_filename.setText('flutter_cards.bdf')
        self.aero_filename.setToolTip('Path to the Flutter File')
        self.aero_filename_load = QPushButton('Load...')

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
            grid.addWidget(self.constant_value, irow, 2)
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
            grid.addWidget(self.constant_value, irow, 1)
            grid.addWidget(self.constant_unit_pulldown, irow, 2)

            irow += 1

            grid.addWidget(self.aero_filename_label, irow, 0)
            grid.addWidget(self.aero_filename, irow, 1)
            grid.addWidget(self.aero_filename_load, irow, 2)
            irow += 1

            grid.addWidget(self.unit_system_label, irow, 0)
            grid.addWidget(self.unit_system_pulldown, irow, 1)
        return grid

    def set_connections(self):
        self.sweep_pulldown.currentIndexChanged.connect(self.on_sweep_pulldown)
        self.constant_pulldown.currentIndexChanged.connect(self.on_constant_pulldown)
        self.apply_button.clicked.connect(self.on_apply)



def cmd_line_gui():
    from qtpy.QtWidgets import (
        QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
        QSpinBox, QDoubleSpinBox, QColorDialog, QLineEdit, QCheckBox,
        QTabWidget, QWidget, QComboBox, QHBoxLayout, QVBoxLayout,
    )
    import pyNastran
    from pyNastran.gui.utils.qt.pydialog import PyDialog, make_font, check_color
    from pyNastran.gui.utils.qt.pydialog import QFloatEdit, QIntEdit
    from pyNastran.gui.utils.qt.checks.qlineedit import (
        check_int, check_float, check_name_str, check_path, QLINEEDIT_GOOD, QLINEEDIT_ERROR)
    from pyNastran.utils import print_bad_path
    from pyNastran.converters.cart3d.cart3d import read_cart3d
    from pyNastran.converters.tecplot.tecplot import read_tecplot
    from pyNastran.dev.tools.pressure_map_aero_setup import get_aero_model

    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
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

if __name__ == '__main__':
    cmd_line_gui()
