# encoding: utf-8
import os
from typing import Any

import numpy as np
from qtpy import QtGui
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


PKG_PATH = pyNastran.__path__[0]
CART3D_PATH = os.path.join(PKG_PATH, 'converters', 'cart3d', 'models', 'threePlugs.a.tri')
TECPLOT_PATH = os.path.join(PKG_PATH, 'converters', 'tecplot', 'models', 'ascii', 'point_febrick_3d_02.dat')
NASTRAN_PATH = os.path.join(PKG_PATH, '..', 'models', 'bwb', 'bwb_saero.bdf')

ELEMENT_ID_OPTIONS = ['CSV File', 'Load ID']
#AERO_FORMATS = ['Cart3D', 'Fund3D', 'Fluent Vrt', 'Fluent Press', 'Tecplot']
AERO_VARIABLE_STR = 'Define the Pressure or Cp (Coefficient of Pressure)'
AERO_FILE_EXTENSION = {
    'Cart3D' : ['triq'],
    'Fund3D': ['triq'],
    'Fluent Vrt': ['vrt'],
    'Fluent Press': ['pres'],
    'Tecplot': ['dat', 'plt'],
}
STRUCTURE_FILE_EXTENSION = {
    'Nastran' : ['bdf', 'pch'],
}
AERO_FORMATS = list(AERO_FILE_EXTENSION)


def get_aero_model(aero_filename: str, aero_format: str) -> tuple[Any, list[str]]:
    variables = []
    if aero_format == 'Cart3D':
        model = read_cart3d(aero_filename)
        variables = list(model.loads)
    # elif aero_format == 'Fund3D':
    #     return None, []
    # elif aero_format == 'Fluent Vrt':
    #     pass
    # elif aero_format == 'Fluent Press':
    #     pass
    elif aero_format == 'Tecplot':
        model = read_tecplot(aero_filename)
        #print(model.object_stats())
        variables = model.result_variables
    else:  # pragma: no cover
        return None, []
        #raise NotImplementedError(aero_format)
    return model, variables


class PressureMap(PyDialog):
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
        self.dim_max = 8

        self.setWindowTitle('Pressure Map')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.on_font(self.font_size)
        #self.show()

    def create_widgets(self):
        """creates the display window"""
        # window text size
        #self.font_size_label = QLabel('Font Size:')
        #self.font_size_edit = QSpinBox(self)

        #self.font_size_edit.setValue(self._default_font_size)
        #self.font_size_edit.setRange(FONT_SIZE_MIN, FONT_SIZE_MAX)

        self.aero_format_label = QLabel('Aero Format:')
        self.aero_format = QComboBox(self)
        self.aero_format.addItems(AERO_FORMATS)

        self.aero_filename_label = QLabel('Aero File:', self)
        self.aero_filename = QLineEdit(self)
        self.aero_filename.setToolTip('Path to the Aero File')
        self.aero_filename_load_button = QPushButton('ALoad...', self)
        self.aero_filename.setText(CART3D_PATH)

        self.aero_variable_label = QLabel('Aero Variable:')
        self.aero_variable = QComboBox(self)
        self.aero_variable.addItems(['Cp'])
        self.aero_variable.setToolTip(f'{AERO_VARIABLE_STR}\nLoad the Aero File')
        self.aero_variable.setEnabled(False)

        self.aero_xyz_scale_label = QLabel('XYZ Scale:', self)
        self.aero_xyz_scale = QFloatEdit('1.0')
        self.aero_xyz_scale_default = QPushButton('Default')
        self.aero_xyz_scale.setToolTip('Scales the Aero Model \n(xyz_new = scale * (offset + xyz))')

        self.aero_txyz_scale_label = QLabel('Offset:')
        self.aero_txyz_scale = QLineEdit('0.0,0.0,0.0')
        self.aero_txyz_scale_default = QPushButton('Default')
        self.aero_txyz_scale.setToolTip('Shifts the Aero Model')

        self.aero_dynamic_pressure_label = QLabel('Dynamic Pressure:')
        self.aero_dynamic_pressure = QFloatEdit('1.0')
        self.aero_dynamic_pressure_default = QPushButton('Default')
        self.aero_dynamic_pressure.setToolTip('Scales the Aero Pressure (q<sub>∞</sub>)\n'
                                              'q<sub>∞</sub> (p=Cp*q<sub>∞</sub>)')

        #-------------------------------------------------
        self.aero_model_label = QLabel('Aero Model:')
        self.structure_model_label = QLabel('Structural Model:')
        self.pressure_label = QLabel('Pressure Output:')

        self.structure_format_label = QLabel('Structure Format:')
        self.structure_format = QComboBox(self)
        self.structure_format.addItems(['Nastran'])
        self.structure_format.setEnabled(False)

        self.element_id_source_label = QLabel('Source of Structural Element IDs:')
        self.element_id_source = QComboBox(self)
        self.element_id_source.addItems(ELEMENT_ID_OPTIONS)

        self.load_aero_button = QPushButton('Load Aero Model')
        self.load_structure_button = QPushButton('Load Structure Model')
        # self.load_aero_button.setStyleSheet('QPushButton {background-color: #A3C1DA; border:  none}')
        # self.load_aero_button.setStyleSheet('QPushButton {background-color: red}')
        self.load_aero_button.setStyleSheet('QPushButton {color: red; font-weight: bold;}')
        self.load_structure_button.setStyleSheet('QPushButton {color: red; font-weight: bold;}')

        self.eid_csv_filename_label = QLabel('Element IDs: CSV Filename:')
        self.eid_csv_filename = QLineEdit(self)
        self.eid_csv_filename.setToolTip('Path to the CSV File')
        self.eid_csv_filename_button = QPushButton('CSV Load...')
        # self.eid_csv_filename.setEnabled(False)
        self.eid_csv_filename_button.setEnabled(False)

        self.eid_load_id_label = QLabel('Element IDs: Load ID:')
        self.eid_load_id_pulldown = QComboBox(self)
        self.eid_load_id_pulldown.addItems(['1'])
        self.eid_load_id_pulldown.setToolTip('The Load ID specifies the mapped elements based on '
                                    'the PLOAD, PLOAD2, or PLOAD4 entries.')
        self.eid_load_id_pulldown.setEnabled(False)

        self.map_type_label = QLabel('Structure Format:')
        self.map_type = QComboBox(self)
        self.map_type.addItems(['Pressure (PLOAD)', 'Force'])

        self.structure_filename_label = QLabel('Structure File:')
        self.structure_filename = QLineEdit(self)
        self.structure_filename.setToolTip('Path to the Structure File')
        self.structure_filename_load_button = QPushButton('SLoad...')
        self.structure_filename.setText(NASTRAN_PATH)

        self.pressure_filename_label = QLabel('Pressure File:')
        self.pressure_filename = QLineEdit(self)
        self.pressure_filename.setText('pressure.bdf')
        self.pressure_filename.setToolTip('Path to the Pressure File')
        self.pressure_filename_load_button = QPushButton('PLoad...')

        self.structure_load_id_label = QLabel('Load ID:')
        self.structure_load_id = QIntEdit('1')
        self.structure_load_id_button = QPushButton('Default')

        #------------------------------------------------------------------
        # closing
        self.apply_button = QPushButton('Apply')
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Cancel')

    def create_layout(self):
        agrid = self._create_aero_grid()
        sgrid, sgrid2 = self._create_structure_grid()
        #
        awidget = QWidget(self)
        swidget = QWidget(self)
        swidget2 = QWidget(self)
        awidget.setLayout(agrid)
        swidget.setLayout(sgrid)
        swidget2.setLayout(sgrid2)

        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)
        ok_widget = QWidget(self)
        ok_widget.setLayout(ok_cancel_box)

        #------------------------------
        vbox = QVBoxLayout(self)
        vbox.addWidget(self.aero_model_label)
        vbox.addWidget(awidget)
        vbox.addWidget(self.load_aero_button)

        vbox.addWidget(self.structure_model_label)
        vbox.addWidget(swidget)
        vbox.addWidget(self.load_structure_button)
        vbox.addWidget(self.pressure_label)
        vbox.addWidget(swidget2)

        vbox.addStretch()
        #vbox.addLayout(ok_cancel_box)
        vbox.addWidget(ok_widget)
        self.setLayout(vbox)

    def _create_aero_grid(self):
        grid = QGridLayout(self)
        irow = 0
        grid.addWidget(self.aero_format_label, irow, 0)
        grid.addWidget(self.aero_format, irow, 1)
        # grid.addWidget(self.font_size_edit, irow, 1)
        irow += 1

        grid.addWidget(self.aero_filename_label, irow, 0)
        grid.addWidget(self.aero_filename, irow, 1)
        grid.addWidget(self.aero_filename_load_button, irow, 2)
        irow += 1

        grid.addWidget(self.aero_xyz_scale_label, irow, 0)
        grid.addWidget(self.aero_xyz_scale, irow, 1)
        grid.addWidget(self.aero_xyz_scale_default, irow, 2)
        irow += 1

        grid.addWidget(self.aero_txyz_scale_label, irow, 0)
        grid.addWidget(self.aero_txyz_scale, irow, 1)
        grid.addWidget(self.aero_txyz_scale_default, irow, 2)
        irow += 1

        grid.addWidget(self.aero_variable_label, irow, 0)
        grid.addWidget(self.aero_variable, irow, 1)
        irow += 1

        grid.addWidget(self.aero_xyz_scale_label, irow, 0)
        grid.addWidget(self.aero_xyz_scale, irow, 1)
        grid.addWidget(self.aero_xyz_scale_default, irow, 2)
        irow += 1

        grid.addWidget(self.aero_txyz_scale_label, irow, 0)
        grid.addWidget(self.aero_txyz_scale, irow, 1)
        grid.addWidget(self.aero_txyz_scale_default, irow, 2)
        irow += 1

        grid.addWidget(self.aero_dynamic_pressure_label, irow, 0)
        grid.addWidget(self.aero_dynamic_pressure, irow, 1)
        grid.addWidget(self.aero_dynamic_pressure_default, irow, 2)
        irow += 1
        return grid

    def _create_structure_grid(self):
        grid = QGridLayout(self)
        grid2 = QGridLayout(self)
        irow = 0
        grid.addWidget(self.structure_format_label, irow, 0)
        grid.addWidget(self.structure_format, irow, 1)
        irow += 1

        grid.addWidget(self.structure_filename_label, irow, 0)
        grid.addWidget(self.structure_filename, irow, 1)
        grid.addWidget(self.structure_filename_load_button, irow, 2)
        irow += 1

        #-------------------------------------------------------------
        irow = 0

        grid2.addWidget(self.element_id_source_label, irow, 0)
        grid2.addWidget(self.element_id_source, irow, 1)
        irow += 1

        grid2.addWidget(self.eid_csv_filename_label, irow, 0)
        grid2.addWidget(self.eid_csv_filename, irow, 1)
        grid2.addWidget(self.eid_csv_filename_button, irow, 2)
        irow += 1

        grid2.addWidget(self.eid_load_id_label, irow, 0)
        grid2.addWidget(self.eid_load_id_pulldown, irow, 1)
        #grid.addWidget(self.eid_load_id_button, irow, 2)
        irow += 1

        grid2.addWidget(self.map_type_label, irow, 0)
        grid2.addWidget(self.map_type, irow, 1)
        irow += 1

        grid2.addWidget(self.pressure_filename_label, irow, 0)
        grid2.addWidget(self.pressure_filename, irow, 1)
        grid2.addWidget(self.pressure_filename_load_button, irow, 2)
        irow += 1

        grid2.addWidget(self.structure_load_id_label, irow, 0)
        grid2.addWidget(self.structure_load_id, irow, 1)
        grid2.addWidget(self.structure_load_id_button, irow, 2)
        irow += 1
        return grid, grid2

    def set_connections(self):
        self.load_aero_button.clicked.connect(self.on_aero)
        self.load_structure_button.clicked.connect(self.on_structure)
        self.aero_filename_load_button.clicked.connect(self.on_aero_load)
        self.structure_filename_load_button.clicked.connect(self.on_structure_load)
        self.pressure_filename_load_button.clicked.connect(self.on_pressure_load)

        return
        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)

    def on_aero_load(self):
        aero_format = self.aero_format.currentText()
        if aero_format == 'Cart3D':
            aero_filename = (CART3D_PATH)
        if aero_format == 'Tecplot':
            aero_filename = TECPLOT_PATH
        else:  # pragma: no cover
            return
            #raise NotImplementedError(aero_format)
        self.aero_filename.setText(aero_filename)
        #self._on_file_load(f'Load Existing {aero_format} Aero Model', AERO_FILE_EXTENSION)

    def on_structure_load(self):
        self._on_file_load('Load Existing Nastran Structure Model', STRUCTURE_FILE_EXTENSION)

    def on_pressure_load(self):
        self._on_file_load('Load Pressure Nastran File', STRUCTURE_FILE_EXTENSION)

    def _on_file_load(self, title: str, file_extensions: dict[str, list[str]]) -> str:
        """TODO: implement this..."""
        return ''

    def on_aero(self):
        aero_filename, aero_passed = check_path(self.aero_filename)
        if not aero_passed:
            self.load_aero_button.setStyleSheet('QPushButton {color: red; font-weight: bold;}')
            self.load_aero_button.setToolTip("Couldn't find Aero File")
            self.aero_variable.setToolTip(f'{AERO_VARIABLE_STR}\nLoad the Aero File')
            return

        aero_format = self.aero_format.currentText()
        aero_model, variables = get_aero_model(aero_filename, aero_format)

        if not list(variables):
            self.load_aero_button.setStyleSheet('QPushButton {color: red; font-weight: bold;}')
            self.load_aero_button.setToolTip('No variables found in Aero File')
            self.aero_variable.setToolTip(f'{AERO_VARIABLE_STR}\nLoad the Aero File')
            return
        self.aero_variable.clear()
        self.aero_variable.addItems(variables)
        self.aero_variable.setEnabled(True)
        self.load_aero_button.setStyleSheet('QPushButton {color: black; font-weight: bold;}')
        self.aero_variable.setToolTip(AERO_VARIABLE_STR)

    def on_structure(self):
        from pyNastran.bdf.bdf import read_bdf
        bdf_filename, structure_passed = check_path(self.structure_filename)
        if not structure_passed:
            self.load_structure_button.setToolTip("Couldn't find Structure File")
            #self.eid_csv_filename.setEnabled(False)
            self.eid_load_id_pulldown.setEnabled(False)
            return

        model = read_bdf(bdf_filename)
        load_ids = list(model.loads) + list(model.load_combinations)
        uload_ids = np.unique(load_ids)
        str_load_ids = [str(val) for val in uload_ids]
        self.eid_load_id_pulldown.clear()

        is_eid_load_id_pulldown = len(str_load_ids)
        if is_eid_load_id_pulldown:
            self.eid_load_id_pulldown.addItems(str_load_ids)
        else:
            self.element_id_source.clear()
            self.element_id_source.addItems([ELEMENT_ID_OPTIONS[0]])
            self.element_id_source.setEnabled(False)
        self.eid_load_id_pulldown.setEnabled(is_eid_load_id_pulldown)

        #self.eid_csv_filename.setEnabled(True)
        self.load_structure_button.setStyleSheet('QPushButton {color: black; font-weight: bold;}')

def main():  # pragma: no cover
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
    #The Main window
    main_window = PressureMap(data)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":  # pragma: no cover
    main()
