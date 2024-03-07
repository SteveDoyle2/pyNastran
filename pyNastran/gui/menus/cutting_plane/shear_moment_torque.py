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
import os
from typing import Callable, TYPE_CHECKING

import numpy as np

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
    QColorDialog, QLineEdit, QCheckBox, QComboBox, QSpinBox,
    QFrame, QTableWidget, QTableWidgetItem, QDialog, QHeaderView)

from qtpy.QtGui import QColor# , QHeaderView


from pyNastran.utils.locale import func_str
from pyNastran.gui.utils.qt.pydialog import PyDialog, QFloatEdit, make_font, check_color, check_patran_syntax
from pyNastran.gui.utils.qt.qelement_edit import (
    QElementLineEdit, QElementTextEdit, # QNodeLineEdit
)
from pyNastran.gui.utils.qt.resize_qtextedit import AutoResizingTextEdit
from pyNastran.gui.utils.qt.qcombobox import make_combo_box # get_combo_box_text # set_combo_box_text,
from pyNastran.gui.utils.qt.qpush_button_color import QPushButtonColor
from pyNastran.gui.utils.qt.dialogs import save_file_dialog
from pyNastran.gui.utils.qt.checks.qlineedit import check_save_path, check_float, QLINEEDIT_GOOD
from pyNastran.gui.utils.wildcards import wildcard_csv
from pyNastran.gui.menus.cutting_plane.cutting_plane import get_zaxis
from pyNastran.gui.menus.preferences.preferences import (
    create_shear_moment_torque_edits)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui_objects.settings import Settings
    from .shear_moment_torque_object import ShearMomentTorqueObject
    from pyNastran.gui.typing import ColorInt, ColorFloat
    from pyNastran.gui.main_window import MainWindow

IS_DEMO = False
#IS_DEMO = True  # just for testing
CID_GLOBAL_STR = '0/Global'

IS_TIME = False
MAX_LENGTH = 100_000
USE_LINE_EDIT = True

class ResultsDialog(QDialog):
    def __init__(self, win_parent,
                 data: np.ndarray,
                 labels: list[str],
                 title: str='Results'):
        super().__init__(win_parent)

        self.setWindowTitle(title)
        nrows, ncolumns = data.shape

        table_widget = QTableWidget(self)
        table_widget.setRowCount(nrows)
        table_widget.setColumnCount(ncolumns)
        table_widget.setHorizontalHeaderLabels(labels)
        self.table_widget = table_widget

        header = table_widget.horizontalHeader()
        for irow, row in enumerate(data):
            for jcol, value in enumerate(row):
                obj = QTableWidgetItem(str(value))
                table_widget.setItem(irow, jcol, obj)
            header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.Stretch)

        vbox = QVBoxLayout(self)
        vbox.addWidget(table_widget)
        self.setLayout(vbox)
        self.show()


class ShearMomentTorqueWindow(PyDialog):
    """
    +-------------------------+
    | ShearMomentTorqueWindow |
    +-------------------------+
    | Origin      cid  x y z  |
    | P2          cid  x y z  |
    | z-axis      cid  x y z  |
    | tol         cid  x y z  |
    |                         |
    |    Apply OK Cancel      |
    +-------------------------+
    """
    def __init__(self, data, win_parent=None):
        """
        Saves the data members from data and
        performs type checks
        """
        PyDialog.__init__(self, data, win_parent)

        self._updated_preference = False

        self._default_font_size = data['font_size']

        #self.dim_max = data['dim_max']
        self.model_name = data['model_name']
        self.icase = data['icase']
        self.cids = data['cids']
        self.elements_pound = data['elements_pound']
        self.gpforce = None
        #self.gpforce = data['gpforce']
        #self._origin = data['origin']
        #self._p1 = data['origin']
        #self._p2 = data['origin']

        #self.out_data = data

        self.plane_color_float, self.plane_color_int = check_color(
            data['plane_color'])
        self.plane_opacity = data['plane_opacity']
        self.vector_line_width = data['vector_line_width']
        self.vector_point_size = data['vector_point_size']

        self.methods = ['Vector', 'CORD2R', 'Coord ID']
        #self.zaxis_methods = ['Global Z', 'Camera Normal', 'Manual']
        self.zaxis_methods = ['Manual', 'Global Z']

        self._icord2r = self.methods.index('CORD2R')
        self._imanual = self.zaxis_methods.index('Manual')
        self._zaxis_method = 0  # Global Z - nope...

        self.setWindowTitle('Shear, Moment, Torque')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        self.on_font(self._default_font_size)
        #self.on_gradient_scale()
        #self.show()

    def on_font(self, value=None) -> None:
        """update the font for the current window"""
        if value in (0, None):
            value = self.font_size_edit.value()
        font = make_font(value, is_bold=False)
        self.setFont(font)

    def set_font_size(self, font_size: int) -> None:
        """
        Updates the font size of all objects in the PyDialog

        Parameters
        ----------
        font_size : int
            the font size
        """
        if self.font_size == font_size:
            return
        self.font_size = font_size
        font = make_font(font_size, is_bold=False)
        self.setFont(font)
        self.set_bold_font(font_size)

    def set_bold_font(self, font_size: int) -> None:
        """
        Updates the font size of all bolded objects in the dialog

        Parameters
        ----------
        font_size : int
            the font size
        """
        bold_font = make_font(font_size, is_bold=True)

        self.additional_params_label.setFont(bold_font)
        self.case_info_label.setFont(bold_font)
        self.plane_label.setFont(bold_font)

        self.location_label.setFont(bold_font)
        self.cid_label.setFont(bold_font)
        self.x_label.setFont(bold_font)
        self.y_label.setFont(bold_font)
        self.z_label.setFont(bold_font)

        self.plot_info.setFont(bold_font)
        self.unit_label.setFont(bold_font)
        self.scale_label.setFont(bold_font)

    def create_widgets(self) -> None:
        """creates the display window"""
        # CORD2R
        #self.origin_label = QLabel("Origin:")
        #self.zaxis_label = QLabel("Z Axis:")
        #self.xz_plane_label = QLabel("XZ Plane:")

        desc = AutoResizingTextEdit(
            'Creates a shear force/bending moment diagram by creating '
            'a series of section cuts.')
        desc.append('')
        desc.append(' 1. Create a vector to march down (from the origin/start '
                    'to end) by defining two points.')
        desc.append('')
        desc.append(
            ' 2. Define an output coordinate system about the origin/start.')
        desc.append(
            'The goal is to orient the cutting plane in the direction that '
            'want to expose forces/moments.')
        desc.append('In other words, point the x-axis roughly down the '
                    'march axis to define what torque is.')
        desc.append('')
        desc.setReadOnly(True)
        desc.viewport().setAutoFillBackground(False)
        desc.setFrameStyle(QFrame.NoFrame)

        self.description = desc

        # Z-Axis Projection
        self.p1_label = QLabel('Origin:')
        self.p3_label = QLabel('End:')
        self.p2_label = QLabel('XZ Plane:')
        self.p1_label.setToolTip('Defines the starting point for the shear, moment, torque plot')
        self.p3_label.setToolTip('Defines the end point for the shear, moment, torque plot')
        self.p2_label.setToolTip('Defines the XZ plane for the shears/moments')

        self.station_location_label = QLabel('Station Label:')
        self.station_location_pulldown = QComboBox()
        self.station_location_pulldown.addItems(['End-Origin', 'X', 'Y', 'Z'])
        self.station_location_pulldown.setToolTip('Sets the x-axis on the 2D plots')

        self.zaxis_label = QLabel('Z Axis:')

        self.method_pulldown = QComboBox()
        for method in self.methods:
            self.method_pulldown.addItem(method)
        self.method_pulldown.setToolTip(
            'Define the output coordinate system\n'
            ' - Vector:    Z-axis and xz-axis of plane  (x is normal to plane)\n'
            ' - CORD2R:    Points on z-axis and xz-axis (x is normal to plane)\n'
            ' - Coord ID : Direct output coordinate system (x is normal to plane)\n'
        )

        self.zaxis_method_pulldown = QComboBox()
        self.zaxis_method_pulldown.setToolTip(
            'Define the output coordinate system\n'
            ' - "Start to XZ-Plane" defines x-axis\n\n'
            'Options for Z-Axis:\n'
            ' - Global Z: z=<0, 0, 1>\n'
            #' - Camera Normal: depends on orientation of model (out of the page)\n'
            ' - Manual: Explicitly define the z-axis'
        )
        for method in self.zaxis_methods:
            self.zaxis_method_pulldown.addItem(method)

        self.cid_label = QLabel('Coordinate System:')
        self.p1_cid_pulldown = QComboBox()
        self.p2_cid_pulldown = QComboBox()
        self.p3_cid_pulldown = QComboBox()
        self.zaxis_cid_pulldown = QComboBox()

        for cid in sorted(self.cids):
            if cid == 0:
                cid_str = CID_GLOBAL_STR
            else:
                cid_str = str(cid)
            #print('cid_str = %r' % cid_str)
            self.p1_cid_pulldown.addItem(cid_str)
            self.p2_cid_pulldown.addItem(cid_str)
            self.p3_cid_pulldown.addItem(cid_str)
            self.zaxis_cid_pulldown.addItem(cid_str)

        self.p1_cid_pulldown.setCurrentIndex(0)
        self.p2_cid_pulldown.setCurrentIndex(0)
        self.p3_cid_pulldown.setCurrentIndex(0)
        self.zaxis_cid_pulldown.setCurrentIndex(0)
        if len(self.cids) == 1:
            self.p1_cid_pulldown.setEnabled(False)
            self.p2_cid_pulldown.setEnabled(False)
            self.p3_cid_pulldown.setEnabled(False)
            self.zaxis_cid_pulldown.setEnabled(False)

        #self.p1_cid_pulldown.setItemText(0, cid_str)
        #self.p2_cid_pulldown.setItemText(0, cid_str)
        #self.zaxis_cid_pulldown.setItemText(0, cid_str)

        self.p1_cid_pulldown.setToolTip('Defines the coordinate system for Point P1/starting point')
        self.p2_cid_pulldown.setToolTip('Defines the coordinate system for Point P2/xz-plane point')
        self.p3_cid_pulldown.setToolTip('Defines the coordinate system for Point P3/ending point')
        self.zaxis_cid_pulldown.setToolTip('Defines the coordinate system for the Z Axis')

        self.p1_x_edit = QFloatEdit('')
        self.p1_y_edit = QFloatEdit('')
        self.p1_z_edit = QFloatEdit('')

        self.p2_x_edit = QFloatEdit('')
        self.p2_y_edit = QFloatEdit('')
        self.p2_z_edit = QFloatEdit('')

        self.p3_x_edit = QFloatEdit('')
        self.p3_y_edit = QFloatEdit('')
        self.p3_z_edit = QFloatEdit('')

        self.zaxis_x_edit = QFloatEdit('')
        self.zaxis_y_edit = QFloatEdit('')
        self.zaxis_z_edit = QFloatEdit('')

        self.additional_params_label = QLabel('Plane Parameters:')
        self.case_info_label = QLabel('Case Info:')

        self.p2_label = QLabel('XZ Plane:')

        # Plane Color


        # ------------------------------------------------
        opacity_edit, point_size_edit, line_width_edit, color_edit = create_shear_moment_torque_edits(
            self,
            self.plane_opacity,
            self.vector_point_size,
            self.vector_line_width,
            self.plane_color_int)

        self.point_size_label = QLabel("Point Size:")
        self.point_size_edit = point_size_edit

        self.line_width_label = QLabel("Line Width:")
        self.line_width_edit = line_width_edit

        self.plane_color_label = QLabel('Plane Color:')
        self.plane_color_edit = color_edit

        self.plane_opacity_label = QLabel('Plane Opacity:')
        self.plane_opacity_edit = opacity_edit

        if 0:
            self.flip_coord_label = QLabel('Flip Coordinate System:')
            self.flip_coord_checkbox = QCheckBox()

        #-----------------------------------------------------------------------
        self.icase_label = QLabel('iCase:')
        self.icase_edit = QSpinBox()
        self.icase_edit.setMinimum(0)
        self.icase_edit.setMaximum(1000)
        self.icase_edit.setValue(self.icase)
        self.icase_edit.setToolTip('Defines the GridPointForces result to be analyzed.\n'
                                   'Check the log for the case id')

        if IS_TIME:
            self.time_label = QLabel('Time:')
            if self.gpforce is None:  # pragma: no cover
                # for debugging; not real
                times = ['0.', '0.5', '1.' , '1.5', '2.']
                time = '0.'
            else:
                times = [func_str(time) for time in self.gpforce._times]
                time = times[0]
            self.times_pulldown = make_combo_box(times, time)
            self.time_label.setEnabled(False)
            self.times_pulldown.setEnabled(False)


        #name = 'main'
        #win_parent = self
        #parent = self.gui
        #self.element_edit = QElementEdit(
            #win_parent, name, parent=parent, pick_style='area',
            #tab_to_next=True, cleanup=True, max_length=32767)
        #self.node_edit = QNodeLineEdit(
            #win_parent, name, parent=parent, pick_style='area',
            #tab_to_next=True, cleanup=True)

        #gui = self.gui
        gui = self.win_parent
        #gui = self
        #self.element_node_checkbox = QCheckBox('Limit Nodes/Elements:')
        self.element_node_checkbox = QCheckBox('Limit Elements:')
        #self.node_label = QLabel('Nodes:')
        #self.node_edit = QNodeLineEdit(
            #self.win_parent, self.model_name, parent=gui,
            #pick_style='area', tab_to_next=False, max_length=MAX_LENGTH)

        self.element_label = QLabel('Elements:')
        if USE_LINE_EDIT:
            self.element_edit = QElementLineEdit(
                self, self.model_name, parent=gui,
                pick_style='area', tab_to_next=False, max_length=MAX_LENGTH)
        else:
            self.element_edit = QElementTextEdit(
                self, self.model_name, parent=gui,
                pick_style='area', tab_to_next=False)
        self.on_element_node_checkbox()

        #self.node_element_label = QLabel('Nodes/Elements:')
        #self.node_element_edit = QLineEdit()
        #self.node_element_edit.setReadOnly(True)

        self.nplanes_label = QLabel('Num Planes:')
        self.nplanes_spinner = QSpinBox()
        self.nplanes_spinner.setMinimum(2)
        self.nplanes_spinner.setMaximum(500)
        self.nplanes_spinner.setValue(20)

        #-----------------------------------------------------------------------
        self.method_label = QLabel('Method:')
        self.plane_label = QLabel('Plane:')
        self.location_label = QLabel('Location:')
        self.zaxis_method_label = QLabel('Z-Axis Method:')
        self.cid_label = QLabel('Coordinate System:')
        self.x_label = QLabel('X')
        self.y_label = QLabel('Y')
        self.z_label = QLabel('Z')

        if 'Z-Axis Projection' not in self.methods:
            self.zaxis_method_label.setVisible(False)

        #self.location_label.setAlignment(Qt.AlignCenter)
        self.cid_label.setAlignment(Qt.AlignCenter)

        self.x_label.setAlignment(Qt.AlignCenter)
        self.y_label.setAlignment(Qt.AlignCenter)
        self.z_label.setAlignment(Qt.AlignCenter)

        self.export_checkbox = QCheckBox()
        self.csv_label = QLabel('CSV Filename:')
        self.csv_edit = QLineEdit()
        self.csv_button = QPushButton('Browse...')

        default_dirname = os.getcwd()
        if self.win_parent is not None:
            default_dirname = self.win_parent.last_dir
        default_filename = os.path.join(default_dirname, 'shear_moment_torque.csv')
        self.csv_edit.setText(default_filename)

        #self.csv_label.setEnabled(False)
        self.csv_edit.setEnabled(False)
        self.csv_button.setEnabled(False)
        #-----------------------------------------------------------------------
        # nodes
        self.add_button = QPushButton('Add')
        self.remove_button = QPushButton('Remove')

        # elements
        self.add2_button = QPushButton('Add')
        self.remove2_button = QPushButton('Remove')
        #-----------------------------------------------------------------------
        self.plot_info = QLabel('Plot Info:')
        self.force_label = QLabel('Force:')
        self.moment_label = QLabel('Moment:')
        self.length_label = QLabel('Length:')
        self.unit_label = QLabel('Unit')
        self.scale_label = QLabel('Scale')

        self.unit_label.setAlignment(Qt.AlignCenter)
        self.scale_label.setAlignment(Qt.AlignCenter)

        self.length_unit_edit = QLineEdit('')
        self.force_unit_edit = QLineEdit('')
        self.moment_unit_edit = QLineEdit('')
        self.length_scale_edit = QFloatEdit('1.0')
        self.force_scale_edit = QFloatEdit('1.0')
        self.moment_scale_edit = QFloatEdit('1.0')

        self.length_unit_edit.setToolTip('Define the length unit for the output')
        self.force_unit_edit.setToolTip('Define the force unit for the output')
        self.moment_unit_edit.setToolTip('Define the moment unit for the output')

        self.length_scale_edit.setToolTip('Scale the output length by this')
        self.force_scale_edit.setToolTip('Scale the output force by this')
        self.moment_scale_edit.setToolTip('Scale the output moment by this')
        #-----------------------------------------------------------------------
        # closing
        self.plot_plane_button = QPushButton('Plot Plane')
        self.clear_plane_button = QPushButton('Clear Plane')
        self.apply_button = QPushButton('Apply')
        self.cancel_button = QPushButton('Cancel')
        self.set_bold_font(self._default_font_size)

        if IS_DEMO:  # pragma: no cover
            if 0:  # bwb
                # start
                self.p1_x_edit.setText('1389')
                self.p1_y_edit.setText('1262')
                self.p1_z_edit.setText('87')

                # end
                self.p3_x_edit.setText('911')
                self.p3_y_edit.setText('0.1')
                self.p3_z_edit.setText('0.')

                # z-axis
                self.p2_x_edit.setText('0')
                self.p2_y_edit.setText('-1')
                self.p2_z_edit.setText('0')
                self.nplanes_spinner.setValue(50)

            elif 0:  # solid_shell_bar
                self.p1_x_edit.setText('0')
                self.p1_y_edit.setText('0')
                self.p1_z_edit.setText('-2')

                self.p3_x_edit.setText('0')
                self.p3_y_edit.setText('0')
                self.p3_z_edit.setText('3')

                self.p2_x_edit.setText('0')
                self.p2_y_edit.setText('0')
                self.p2_z_edit.setText('1')

                self.zaxis_x_edit.setText('1')
                self.zaxis_y_edit.setText('0')
                self.zaxis_z_edit.setText('0')

                self.nplanes_spinner.setValue(5)
            elif 1:  # wingbox
                self.p1_x_edit.setText('128')
                self.p1_y_edit.setText('0.1')
                self.p1_z_edit.setText('23')

                # end
                self.p3_x_edit.setText('128')
                self.p3_y_edit.setText('96')
                self.p3_z_edit.setText('23')

                # xz-plane
                self.p2_x_edit.setText('0')
                self.p2_y_edit.setText('1')
                self.p2_z_edit.setText('0')

                # z-axis
                self.zaxis_x_edit.setText('1')
                self.zaxis_y_edit.setText('0')
                self.zaxis_z_edit.setText('0')

                self.nplanes_spinner.setValue(20)
                self.element_edit.setText('1633:1856')

    @property
    def gui(self) -> MainWindow:
        if self.win_parent is None:
            return None
        return self.win_parent.parent.gui

    def create_layout(self) -> None:
        """sets up the window"""
        grid = self._make_grid_layout()

        #hbox_csv = QHBoxLayout()
        grid2 = QGridLayout()
        #irow = 0

        #grid2.addWidget(self.node_label, irow, 0)
        #grid2.addWidget(self.node_edit, irow, 1)
        #grid2.addWidget(self.add_button, irow, 2)
        #grid2.addWidget(self.remove_button, irow, 3)
        #irow += 1

        #grid2.addWidget(self.element_label, irow, 0)
        #grid2.addWidget(self.element_edit, irow, 1)
        #grid2.addWidget(self.add2_button, irow, 2)
        #grid2.addWidget(self.remove2_button, irow, 3)
        #irow += 1

        #grid2.addWidget(self.node_element_label, irow, 0)
        #grid2.addWidget(self.node_element_edit, irow, 1)
        #irow += 1

        hbox_csv = QHBoxLayout()
        hbox_csv.addWidget(self.export_checkbox)
        hbox_csv.addWidget(self.csv_label)
        hbox_csv.addWidget(self.csv_edit)
        hbox_csv.addWidget(self.csv_button)
        #----------------------------------------------

        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.plot_plane_button)
        ok_cancel_box.addWidget(self.clear_plane_button)
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.cancel_button)

        vbox = QVBoxLayout()

        vbox.addWidget(self.description)
        vbox.addLayout(grid)
        vbox.addLayout(grid2)
        #vbox.addStretch()


        grid_node_element = QGridLayout()
        #grid_node_element.addWidget(self.node_label, 0, 0)
        #grid_node_element.addWidget(self.node_edit, 0, 1)
        grid_node_element.addWidget(self.element_label, 1, 0)
        grid_node_element.addWidget(self.element_edit, 1, 1)

        element_node_box = QVBoxLayout()
        element_node_box.addWidget(self.element_node_checkbox)
        element_node_box.addLayout(grid_node_element)
        vbox.addLayout(element_node_box)

        vbox.addLayout(hbox_csv)

        #self.element_node_checkbox = QCheckBox('Limit Nodes/ELements:')
        #self.element_label = QLabel('Elements:')
        #self.node_label = QLabel('Nodes:')

        vbox.addStretch()

        #-----------------------
        #vbox.addLayout(add_remove_box)
        vbox.addLayout(ok_cancel_box)
        self.on_method(0)  # self._icord2r
        self.on_zaxis_method(self._imanual)
        self.setLayout(vbox)

    def on_export_checkbox(self) -> None:
        """this is called when the checkbox is clicked"""
        is_checked = self.export_checkbox.isChecked()
        self.csv_label.setEnabled(is_checked)
        self.csv_edit.setEnabled(is_checked)
        self.csv_button.setEnabled(is_checked)

    def on_browse_csv(self) -> None:
        """opens a file dialog"""
        default_filename = self.csv_edit.text()
        csv_filename, wildcard = save_file_dialog(
            self, 'Select the file name for export',
            default_filename, wildcard_csv)
        if not csv_filename:
            return

        if self.win_parent is not None:
            last_dir = os.path.dirname(csv_filename)
            self.win_parent.load_actions._set_last_dir(last_dir)
        self.csv_edit.setText(csv_filename)

    def _make_grid_layout(self) -> QGridLayout:
        """builds the QGridLayout"""
        grid = QGridLayout()
        irow = 0
        #-------------------------
        grid.addWidget(self.location_label, irow, 0)
        grid.addWidget(self.cid_label, irow, 1)
        grid.addWidget(self.x_label, irow, 2)
        grid.addWidget(self.y_label, irow, 3)
        grid.addWidget(self.z_label, irow, 4)
        irow += 1

        add_row(irow, grid,
                self.p1_label,
                self.p1_cid_pulldown,
                self.p1_x_edit, self.p1_y_edit, self.p1_z_edit)
        irow += 1

        add_row(irow, grid,
                self.p3_label,
                self.p3_cid_pulldown,
                self.p3_x_edit, self.p3_y_edit, self.p3_z_edit)
        irow += 1

        # -------------------------------
        grid.addWidget(self.plane_label, irow, 0)
        irow += 1

        grid.addWidget(self.method_label, irow, 0)
        grid.addWidget(self.method_pulldown, irow, 1)
        irow += 1

        grid.addWidget(self.zaxis_method_label, irow, 0)
        grid.addWidget(self.zaxis_method_pulldown, irow, 1)
        irow += 1

        add_row(irow, grid,
                self.zaxis_label,
                self.zaxis_cid_pulldown,
                self.zaxis_x_edit, self.zaxis_y_edit, self.zaxis_z_edit)
        irow += 1

        add_row(irow, grid,
                self.p2_label,
                self.p2_cid_pulldown,
                self.p2_x_edit, self.p2_y_edit, self.p2_z_edit)
        irow += 1

        #-----------------------------------------
        grid.addWidget(self.case_info_label, irow, 0)
        irow += 1

        grid.addWidget(self.icase_label, irow, 0)
        grid.addWidget(self.icase_edit, irow, 1)
        irow += 1

        if IS_TIME:
            grid.addWidget(self.time_label, irow, 0)
            grid.addWidget(self.times_pulldown, irow, 1)
            irow += 1

        grid.addWidget(self.nplanes_label, irow, 0)
        grid.addWidget(self.nplanes_spinner, irow, 1)
        irow += 1

        #-----------------------------------------
        grid.addWidget(self.additional_params_label, irow, 0)
        irow += 1

        grid.addWidget(self.plane_color_label, irow, 0)
        grid.addWidget(self.plane_color_edit, irow, 1)
        irow += 1

        grid.addWidget(self.plane_opacity_label, irow, 0)
        grid.addWidget(self.plane_opacity_edit, irow, 1)
        irow += 1

        grid.addWidget(self.point_size_label, irow, 0)
        grid.addWidget(self.point_size_edit, irow, 1)
        irow += 1

        grid.addWidget(self.line_width_label, irow, 0)
        grid.addWidget(self.line_width_edit, irow, 1)
        irow += 1

        # -----------------------------------------
        grid.addWidget(self.plot_info, irow, 0)
        grid.addWidget(self.unit_label, irow, 1)
        grid.addWidget(self.scale_label, irow, 2)
        irow += 1

        grid.addWidget(self.length_label, irow, 0)
        grid.addWidget(self.length_unit_edit, irow, 1)
        grid.addWidget(self.length_scale_edit, irow, 2)
        irow += 1

        grid.addWidget(self.force_label, irow, 0)
        grid.addWidget(self.force_unit_edit, irow, 1)
        grid.addWidget(self.force_scale_edit, irow, 2)
        irow += 1

        grid.addWidget(self.moment_label, irow, 0)
        grid.addWidget(self.moment_unit_edit, irow, 1)
        grid.addWidget(self.moment_scale_edit, irow, 2)
        irow += 1

        grid.addWidget(self.station_location_label, irow, 0)
        grid.addWidget(self.station_location_pulldown, irow, 1)
        irow += 1

        #----------------------------------------------
        return grid

    def set_connections(self) -> None:
        """creates the actions for the menu"""
        self.method_pulldown.currentIndexChanged.connect(self.on_method)
        self.zaxis_method_pulldown.currentIndexChanged.connect(self.on_zaxis_method)
        self.plane_color_edit.clicked.connect(self.on_plane_color)
        self.plane_opacity_edit.valueChanged.connect(self.on_plane_opacity)
        self.point_size_edit.valueChanged.connect(self.on_plane_point_size)
        self.line_width_edit.valueChanged.connect(self.on_plane_line_width)
        self.element_node_checkbox.clicked.connect(self.on_element_node_checkbox)


        self.export_checkbox.clicked.connect(self.on_export_checkbox)
        self.csv_button.clicked.connect(self.on_browse_csv)
        #self.csv_label.clicked.connect(self.on_export_checkbox)

        self.plot_plane_button.clicked.connect(self.on_plot_plane)
        self.clear_plane_button.clicked.connect(self.on_clear_plane)
        self.apply_button.clicked.connect(self.on_apply)
        self.cancel_button.clicked.connect(self.on_cancel)

    def on_method(self, method_int=None) -> None:
        method = get_pulldown_text(method_int, self.methods, self.method_pulldown)

        is_p2_cid_enabled = True
        is_zaxis_cid_enabled = True
        zaxis_method_visible = False
        show_zaxis_xyz = True
        show_p2_xyz = True
        if method == 'CORD2R':
            self._zaxis_method = self.zaxis_method_pulldown.currentIndex()
            # set to manual
            #self.on_zaxis_method(method_int=2)  # manual

            self.plane_label.setText('Points on Plane:')
            self.zaxis_label.setText('Origin + Z Axis:')
            self.p2_label.setText('Origin + XZ Plane:')

            #self.zaxis_method_label.setText('Origin + Z-Axis')
            self.zaxis_method_pulldown.setCurrentIndex(self._imanual)
            self.on_zaxis_method()  # update

        elif method == 'Vector':
            is_p2_cid_enabled = False
            is_zaxis_cid_enabled = False
            self.plane_label.setText('Vectors:')
            self.zaxis_label.setText('Z Axis:')
            self.p2_label.setText('XZ Plane Axis:')
            #self.zaxis_method_label.setText('Z-Axis')
            self.zaxis_method_pulldown.setCurrentIndex(self._imanual)
            self.on_zaxis_method()  # update

        elif method == 'Coord ID':
            is_p2_cid_enabled = True
            is_zaxis_cid_enabled = False
            self.zaxis_method_pulldown.setVisible(False)
            self.zaxis_method_label.setVisible(False)

            self.zaxis_method_pulldown.setVisible(False)

            zaxis_method_visible = False
            show_zaxis_xyz = False
            self.zaxis_method_label.setVisible(False)

            show_p2_xyz = False
            self.p2_label.setText('Output Coord:')

        elif method == 'Z-Axis Projection':
            #is_p2_cid_enabled = False
            is_zaxis_cid_enabled = False
            zaxis_method_visible = True
            self.plane_label.setText('Point on Plane/Vector:')
            self.zaxis_label.setText('Z Axis:')
            self.p2_label.setText('Origin + XZ Plane:')

            #self.zaxis_method_label.setText('Z-Axis')
            self.zaxis_method_pulldown.setCurrentIndex(self._zaxis_method)
            self.on_zaxis_method()  # update
        else:  # pragma: no cover
            raise NotImplementedError(method)

        self.p2_cid_pulldown.setEnabled(is_p2_cid_enabled)
        self.zaxis_cid_pulldown.setEnabled(is_zaxis_cid_enabled)

        self.zaxis_method_pulldown.setEnabled(zaxis_method_visible)
        self.zaxis_method_pulldown.setVisible(zaxis_method_visible)
        self.zaxis_method_label.setEnabled(zaxis_method_visible)

        self.p2_x_edit.setVisible(show_p2_xyz)
        self.p2_y_edit.setVisible(show_p2_xyz)
        self.p2_z_edit.setVisible(show_p2_xyz)

        self.zaxis_x_edit.setVisible(show_zaxis_xyz)
        self.zaxis_y_edit.setVisible(show_zaxis_xyz)
        self.zaxis_z_edit.setVisible(show_zaxis_xyz)

    def on_zaxis_method(self, method_int=None) -> None:
        method = get_pulldown_text(method_int, self.zaxis_methods,
                                   self.zaxis_method_pulldown)

        if method == 'Global Z':
            is_visible = False
        #elif method == 'Camera Normal':
            #is_visible = False
        elif method == 'Manual':
            is_visible = True
        else:  # pragma: no cover
            raise NotImplementedError(method)

        self.zaxis_cid_pulldown.setVisible(is_visible)
        self.zaxis_x_edit.setVisible(is_visible)
        self.zaxis_y_edit.setVisible(is_visible)
        self.zaxis_z_edit.setVisible(is_visible)

    def _update_plane_settings(self) -> None:
        obj: ShearMomentTorqueObject = self.win_parent.shear_moment_torque_obj
        obj.set_plane_properties()
        return

    def on_plane_opacity(self) -> None:
        """ Sets the plane opacity"""
        opacity = self.plane_opacity_edit.value()
        if self.win_parent is not None:
            settings: Settings = self.win_parent.settings
            settings.shear_moment_torque_opacity = opacity
            self._update_plane_settings()

    def on_plane_point_size(self) -> None:
        """ Sets the plane opacity"""
        point_size = self.point_size_edit.value()
        if self.win_parent is not None:
            settings: Settings = self.win_parent.settings
            settings.shear_moment_torque_point_size = point_size
            self._update_plane_settings()

    def on_plane_line_width(self) -> None:
        """ Sets the plane opacity"""
        line_width = self.line_width_edit.value()
        if self.win_parent is not None:
            settings: Settings = self.win_parent.settings
            settings.shear_moment_torque_line_width = line_width
            self._update_plane_settings()

    def on_element_node_checkbox(self) -> None:
        is_checked = self.element_node_checkbox.isChecked()
        #self.node_label.setEnabled(is_checked)
        #self.node_edit.setEnabled(is_checked)
        self.element_label.setEnabled(is_checked)
        self.element_edit.setEnabled(is_checked)

    def on_plane_color(self) -> None:
        """ Choose a plane color"""
        title = 'Choose a cutting plane color'
        rgb_color_ints = self.plane_color_int
        color_edit = self.plane_color_edit
        func_name = 'set_plane_color'
        passed, rgb_color_ints, rgb_color_floats = self._background_color(
            title, color_edit, rgb_color_ints, func_name)
        if passed:
            settings: Settings = self.win_parent.settings
            settings.shear_moment_torque_color = rgb_color_floats
            self.plane_color_int = rgb_color_ints
            self.plane_color_float = rgb_color_floats
            self._update_plane_settings()

    def _background_color(self, title: str,
                          color_edit: QPushButtonColor,
                          rgb_color_ints: ColorInt,
                          func_name: Callable) -> tuple[bool, ColorInt, ColorFloat]:
        """
        helper method for:
         - ``on_background_color``
         - ``on_background_color2``

        """
        passed, rgb_color_ints, rgb_color_floats = self.on_color(
            color_edit, rgb_color_ints, title)
        #if passed and 0:
            #if self.win_parent is not None:
                #settings = self.win_parent.settings
                #func_background_color = getattr(settings, func_name)
                #func_background_color(rgb_color_floats)
        return passed, rgb_color_ints, rgb_color_floats

    def on_color(self, color_edit: QPushButtonColor,
                 rgb_color_ints: ColorInt,
                 title: str) -> tuple[bool, ColorInt, ColorFloat]:
        """pops a color dialog"""
        qcolor = QColor(*rgb_color_ints)
        col = QColorDialog.getColor(qcolor, self, title)
        if not col.isValid():
            return False, rgb_color_ints, None

        color_float: ColorFloat = col.getRgbF()[:3]  # floats
        color_int = [int(colori * 255) for colori in color_float]

        assert isinstance(color_float[0], float), color_float
        assert isinstance(color_int[0], int), color_int

        color_edit.set_color(color_int)
        #color_edit.setStyleSheet(
            #'QPushButton {'
            #'background-color: rgb(%s, %s, %s);' % tuple(color_int) +
            ##"border:1px solid rgb(255, 170, 255); "
            #'}')
        return True, color_int, color_float


    #---------------------------------------------------------------------------

    def on_validate(self) -> bool:
        station_location = self.station_location_pulldown.currentText()

        method = self.method_pulldown.currentText()
        assert method in self.methods, f'method={method!r}'

        p1_cidi = self.p1_cid_pulldown.currentText()
        p2_cidi = self.p2_cid_pulldown.currentText()
        p3_cidi = self.p3_cid_pulldown.currentText()
        zaxis_cidi = self.zaxis_cid_pulldown.currentText()
        p1_cid = int(p1_cidi) if CID_GLOBAL_STR not in p1_cidi else 0
        p2_cid = int(p2_cidi) if CID_GLOBAL_STR not in p2_cidi else 0
        p3_cid = int(p3_cidi) if CID_GLOBAL_STR not in p3_cidi else 0
        zaxis_cid = int(zaxis_cidi) if CID_GLOBAL_STR not in zaxis_cidi else 0
        #print('p1_cidi=%r p2_cidi=%r p3_cidi=%r' % (p1_cidi, p2_cidi, zaxis_cidi))
        #print('p2_cid=%r p2_cid=%r p3_cidi=%r' % (p2_cid, p2_cid, zaxis_cid))

        p1_x, flag1 = check_float(self.p1_x_edit)
        p1_y, flag2 = check_float(self.p1_y_edit)
        p1_z, flag3 = check_float(self.p1_z_edit)
        p1_flag = all([flag1, flag2, flag3])
        p1 = [p1_x, p1_y, p1_z]

        p3_x, flag7 = check_float(self.p3_x_edit)
        p3_y, flag8 = check_float(self.p3_y_edit)
        p3_z, flag9 = check_float(self.p3_z_edit)
        p3_flag = all([flag7, flag8, flag9])
        p3 = [p3_x, p3_y, p3_z]

        if method == 'Coord ID':
            p2_flag = True
            p2 = [0., 0., 0.]

            zaxis_cid = -1
            zaxis = np.full(3, np.nan)
            zaxis_flag = True
        else:
            p2_x, flag4 = check_float(self.p2_x_edit)
            p2_y, flag5 = check_float(self.p2_y_edit)
            p2_z, flag6 = check_float(self.p2_z_edit)
            p2 = [p2_x, p2_y, p2_z]
            p2_flag = all([flag4, flag5, flag6])

            flag10, flag11, flag12, zaxis_cid, zaxis = get_zaxis(
                self.win_parent, # for camera
                self.zaxis_method_pulldown,
                self.zaxis_x_edit, self.zaxis_y_edit, self.zaxis_z_edit)
            zaxis_flag = all([flag10, flag11, flag12])

        nplanes = self.nplanes_spinner.value()

        csv_filename = None
        csv_flag = True
        if self.export_checkbox.isChecked():
            csv_filename, csv_flag = check_save_path(self.csv_edit)

        length_scale, length_flag = check_float(self.length_scale_edit)
        force_scale, force_flag = check_float(self.force_scale_edit)
        moment_scale, moment_flag = check_float(self.moment_scale_edit)

        element_ids = None
        eids_flag = True
        if self.element_node_checkbox.isChecked():
            element_ids, eids_flag = check_patran_syntax(
                self.element_edit, self.elements_pound)
        else:
            self.element_node_checkbox.setStyleSheet(QLINEEDIT_GOOD)


        flags = [
            p1_flag, p2_flag, p3_flag,
            zaxis_flag,
            csv_flag,
            eids_flag,
            length_flag, force_flag, moment_flag]

        length_unit = self.length_unit_edit.text()
        force_unit = self.force_unit_edit.text()
        moment_unit = self.moment_unit_edit.text()

        if all(flags):
            # Z-Axis Method
            # p1: origin
            # p2: xz_plane
            # p3: end
            self.out_data['method'] = method
            self.out_data['station_location'] = station_location
            self.out_data['p1'] = [p1_cid, p1]  # origin
            self.out_data['p2'] = [p2_cid, p2]  # xzplane
            self.out_data['p3'] = [p3_cid, p3]  # end
            self.out_data['zaxis'] = [zaxis_cid, zaxis]
            self.out_data['nplanes'] = nplanes
            self.out_data['csv_filename'] = csv_filename
            self.out_data['length'] = [length_scale, length_unit]
            self.out_data['force'] = [force_scale, force_unit]
            self.out_data['moment'] = [moment_scale, moment_unit]
            self.out_data['element_ids'] = element_ids
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_clear_plane(self) -> None:
        if self.win_parent is not None:
            obj: ShearMomentTorqueObject = self.win_parent.shear_moment_torque_obj
            obj.on_clear_plane_actors()
        return

    def on_plot_plane(self) -> bool:
        passed = self.on_validate()
        if passed and self.win_parent is not None:
            obj: ShearMomentTorqueObject = self.win_parent.shear_moment_torque_obj
            obj.make_plane_from_data(self.out_data)
            #self.win_parent.make_smt_from_data(self.out_data)
        return passed

    def on_apply(self) -> bool:
        passed = self.on_validate()
        if passed and self.win_parent is not None:
            obj: ShearMomentTorqueObject = self.win_parent.shear_moment_torque_obj
            obj.make_smt_from_data(self.out_data, show=True)
            #self.win_parent.make_smt_from_data(self.out_data)
        return passed

    def on_cancel(self) -> None:
        self.out_data['close'] = True
        self.close()

def add_row(irow: int,
            grid: QGridLayout,
            p1_label, p1_cid_pulldown,
            p1_x_edit, p1_y_edit, p1_z_edit) -> None:
    """adds the items to the grid"""
    grid.addWidget(p1_label, irow, 0)
    grid.addWidget(p1_cid_pulldown, irow, 1)
    grid.addWidget(p1_x_edit, irow, 2)
    grid.addWidget(p1_y_edit, irow, 3)
    grid.addWidget(p1_z_edit, irow, 4)


def get_pulldown_text(method_int: int,
                      methods: list[str],
                      pulldown: QComboBox):
    if method_int is None:
        #method = pulldown.getText()
        method = pulldown.currentText()
    else:
        method = methods[method_int]
    return method

def main() -> None:  # pragma: no cover
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)


    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    #The Main window

    #gpforce = None
    data = {
        'font_size' : 8,
        'icase': 42,
        #'cids' : [0, 1, 2, 3],
        'cids' : [0],
        'elements_pound': 5000,
        'plane_color' : (1., 0., 1.), # purple
        'plane_opacity' : 0.9,
        'vector_line_width' : 0.9,
        'vector_point_size' : 0.9,
        #'gpforce' : gpforce,
        #'itime' : 0,
        'word' : 'Static',
        'model_name' : 'main',
    }
    main_window = ShearMomentTorqueWindow(data)
    main_window.show()
    # Enter the main loop
    app.exec_()


if __name__ == '__main__':   # pragma: no cover
    main()
