"""
The highlight menu handles:
 - Nodes/Elements
 - Preferences
"""
from __future__ import print_function
from math import log10, ceil

#import PySide  # for local testing
from qtpy import QtGui
from qtpy.QtWidgets import (
    QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
    QSpinBox, QDoubleSpinBox, QColorDialog, QLineEdit, QCheckBox, QComboBox)
import vtk

from pyNastran.gui.utils.qt.pydialog import PyDialog, check_patran_syntax
from pyNastran.gui.utils.qt.qpush_button_color import QPushButtonColor
from pyNastran.gui.menus.menu_utils import eval_float_from_string
from pyNastran.gui.utils.qt.qelement_edit import QNodeEdit, QElementEdit, QNodeElementEdit

from pyNastran.gui.styles.highlight_style import (
    create_vtk_selection_node_by_cell_ids,
    create_vtk_selection_node_by_point_ids,
    extract_selection_node_from_grid_to_ugrid,
)
class HighlightWindow(PyDialog):
    """
    +-----------+
    | Highlight |
    +--------------------------+
    | Nodes     ______         |
    | Elements  ______         |
    |                          |
    |     Highlight Close      |
    +--------------------------+
    """
    def __init__(self, data, win_parent=None):
        """
        Saves the data members from data and
        performs type checks
        """
        PyDialog.__init__(self, data, win_parent)
        gui = win_parent

        self._updated_highlight = False

        self._default_font_size = data['font_size']
        self.model_name = data['model_name']
        #self._default_annotation_size = data['annotation_size'] # int
        #self.default_magnify = data['magnify']

        if 'nodes_pound' in data: # testing
            nodes_pound = data['nodes_pound']
            elements_pound = data['elements_pound']
            nodes = np.arange(1, nodes_pound+1)
            elements = np.arange(1, elements_pound+1)
        else:
            # gui
            nodes = gui.get_node_ids(model_name=self.model_name)
            elements = gui.get_element_ids(model_name=self.model_name)
            nodes_pound = nodes.max()
            elements_pound = elements.max()

        self.nodes = nodes
        self.elements = elements
        self._nodes_pound = nodes_pound
        self._elements_pound = elements_pound

        self.setWindowTitle('Highlight')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        self.on_font(self._default_font_size)
        #self.show()

    def create_widgets(self):
        """creates the display window"""
        # Text Size

        model_name = ''
        self.nodes_label = QLabel("Nodes:")
        self.nodes_edit = QNodeEdit(self, model_name, pick_style='area', tab_to_next=False)

        self.elements_label = QLabel("Elements:")
        self.elements_edit = QNodeEdit(self, model_name, pick_style='area', tab_to_next=False)

        #-----------------------------------------------------------------------
        # Highlight Color
        #self.highlight_opacity_label = QLabel("Highlight Opacity:")
        #self.highlight_opacity_edit = QDoubleSpinBox(self)
        #self.highlight_opacity_edit.setValue(self._highlight_opacity)
        #self.highlight_opacity_edit.setRange(0.1, 1.0)
        #self.highlight_opacity_edit.setDecimals(1)
        #self.highlight_opacity_edit.setSingleStep(0.1)
        #self.highlight_opacity_button = QPushButton("Default")

        # Text Color
        #self.highlight_color_label = QLabel("Highlight Color:")
        #self.highlight_color_edit = QPushButtonColor(self.highlight_color_int)

        #-----------------------------------------------------------------------
        # closing
        self.show_button = QPushButton("Show")
        self.close_button = QPushButton("Close")

    def create_layout(self):
        grid = QGridLayout()

        irow = 0
        grid.addWidget(self.nodes_label, irow, 0)
        grid.addWidget(self.nodes_edit, irow, 1)
        irow += 1

        grid.addWidget(self.elements_label, irow, 0)
        grid.addWidget(self.elements_edit, irow, 1)
        irow += 1

        #grid.addWidget(self.highlight_color_label, irow, 0)
        #grid.addWidget(self.highlight_color_edit, irow, 1)
        #irow += 1

        #grid.addWidget(self.highlight_opacity_label, irow, 0)
        #grid.addWidget(self.highlight_opacity_edit, irow, 1)
        #irow += 1

        #self.create_legend_widgets()
        #grid2 = self.create_legend_layout()
        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.show_button)
        ok_cancel_box.addWidget(self.close_button)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        #vbox.addStretch()
        #vbox.addLayout(grid2)
        vbox.addStretch()

        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        #self.highlight_color_edit.clicked.connect(self.on_highlight_color)
        #self.highlight_opacity_edit.valueChanged.connect(self.on_highlight_opacity)
        self.nodes_edit.textChanged.connect(self.on_validate)
        self.elements_edit.textChanged.connect(self.on_validate)
        self.show_button.clicked.connect(self.on_show)
        self.close_button.clicked.connect(self.on_close)
        # closeEvent

    def on_font(self, value=None):
        """update the font for the current window"""
        if value is None:
            value = self.font_size_edit.value()
        font = QtGui.QFont()
        font.setPointSize(value)
        self.setFont(font)

    def on_highlight_color(self):
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
            self.win_parent.settings.set_highlight_opacity(value)

    def _background_color(self, title, color_edit, rgb_color_ints, func_name):
        """helper method for ``on_background_color`` and ``on_background_color2``"""
        passed, rgb_color_ints, rgb_color_floats = self.on_color(
            color_edit, rgb_color_ints, title)
        if passed:
            if self.win_parent is not None:
                settings = self.win_parent.settings
                func_background_color = getattr(settings, func_name)
                func_background_color(rgb_color_floats)
        return passed, rgb_color_ints, rgb_color_floats

    def on_color(self, color_edit, rgb_color_ints, title):
        """pops a color dialog"""
        col = QColorDialog.getColor(QtGui.QColor(*rgb_color_ints), self, title)
        if not col.isValid():
            return False, rgb_color_ints, None

        color_float = col.getRgbF()[:3]  # floats
        color_int = [int(colori * 255) for colori in color_float]

        assert isinstance(color_float[0], float), color_float
        assert isinstance(color_int[0], int), color_int

        color_edit.setStyleSheet(
            "QPushButton {"
            "background-color: rgb(%s, %s, %s);" % tuple(color_int) +
            #"border:1px solid rgb(255, 170, 255); "
            "}")
        return True, color_int, color_float

    #---------------------------------------------------------------------------

    def on_validate(self):
        nodes, flag1 = check_patran_syntax(self.nodes_edit, pound=self._nodes_pound)
        elements, flag2 = check_patran_syntax(self.elements_edit, pound=self._elements_pound)
        print('%s nodes=%s' % (flag1, nodes))
        print('%s elements=%s' % (flag2, elements))
        if all([flag1, flag2]):
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_show(self):
        passed = self.on_validate()

        if passed and self.win_parent is not None:
            nodes, flag1 = check_patran_syntax(self.nodes_edit, pound=self._nodes_pound)
            elements, flag2 = check_patran_syntax(self.elements_edit, pound=self._elements_pound)

            import numpy as np
            nodes_filtered = np.union1d(self.nodes, nodes)
            elements_filtered = np.union1d(self.elements, elements)

            point_ids = np.searchsorted(self.nodes, nodes_filtered)
            cell_ids = np.searchsorted(self.elements, elements_filtered)
            selection_node_cells = create_vtk_selecion_node_by_cell_ids(cell_ids)
            selection_node_points = create_vtk_selecion_node_by_cell_ids(cell_ids)

            self.actor = actor
            self.parent.vtk_interactor.Render()


        return passed

    def on_close(self):
        self.out_data['close'] = True
        self.close()


def check_float(cell):
    text = cell.text()
    value = float(text)
    return value, True

def check_label_float(cell):
    text = cell.text()
    try:
        value = eval_float_from_string(text)
        cell.setStyleSheet("QLineEdit{background: white;}")
        return value, True
    except ValueError:
        cell.setStyleSheet("QLineEdit{background: red;}")
        return None, False

def _check_color(color_float):
    assert len(color_float) == 3, color_float
    assert isinstance(color_float[0], float), color_float
    color_int = [int(colori * 255) for colori in color_float]
    return color_float, color_int

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
        'model_name' : 'main',
        'nodes_pound' : 100,
        'elements_pound' : 200,
    }
    main_window = HighlightWindow(data)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":  # pragma: no cover
    main()
