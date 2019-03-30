"""
The highlight menu handles:
 - Nodes/Elements
 - Preferences
"""
from __future__ import print_function
#from math import log10, ceil
import numpy as np

#import PySide  # for local testing
from qtpy import QtGui
from qtpy.QtWidgets import (
    QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
    QDoubleSpinBox, QColorDialog)
#import vtk
from vtk.util.numpy_support import vtk_to_numpy

from pyNastran.gui.utils.qt.pydialog import PyDialog, check_patran_syntax
from pyNastran.gui.utils.qt.qpush_button_color import QPushButtonColor
#from pyNastran.gui.menus.menu_utils import eval_float_from_string
from pyNastran.gui.utils.qt.qelement_edit import QNodeEdit, QElementEdit#, QNodeElementEdit

from pyNastran.gui.qt_files.mouse_actions import create_highlighted_actor
from pyNastran.gui.styles.highlight_style import (
    create_vtk_selection_node_by_cell_ids,
    #create_vtk_selection_node_by_point_ids,
    #create_surface_actor_from_grid_and_cell_ids,
)
from pyNastran.gui.utils.vtk.vtk_utils import (
    extract_selection_node_from_grid_to_ugrid,
    create_unstructured_point_grid, numpy_to_vtk_points, )
from pyNastran.gui.utils.vtk.gui_utils import add_actors_to_gui, remove_actors_from_gui


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

        if gui is None:  # pragma: no cover
            self.highlight_color_float = [0., 0., 0.]
            self.highlight_color_int = [0, 0, 0]
            self._highlight_opacity = 0.9
        else:
            self.highlight_color_float = gui.settings.highlight_color
            self.highlight_color_int = [int(val) for val in self.highlight_color_float]
            #self.highlight_color_int = gui.settings.highlight_color_int
            self._highlight_opacity = gui.settings.highlight_opacity
        self._updated_highlight = False

        self.actors = []
        self._default_font_size = data['font_size']
        self.model_name = data['model_name']
        assert len(self.model_name) > 0, self.model_name

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

        model_name = self.model_name
        self.nodes_label = QLabel("Nodes:")
        self.nodes_edit = QNodeEdit(self, model_name, pick_style='area', tab_to_next=False)

        self.elements_label = QLabel("Elements:")
        self.elements_edit = QElementEdit(self, model_name, pick_style='area', tab_to_next=False)

        #-----------------------------------------------------------------------
        # Highlight Color
        self.highlight_opacity_label = QLabel("Highlight Opacity:")
        self.highlight_opacity_edit = QDoubleSpinBox(self)
        self.highlight_opacity_edit.setValue(self._highlight_opacity)
        self.highlight_opacity_edit.setRange(0.1, 1.0)
        self.highlight_opacity_edit.setDecimals(1)
        self.highlight_opacity_edit.setSingleStep(0.1)
        self.highlight_opacity_button = QPushButton("Default")

        # Text Color
        self.highlight_color_label = QLabel("Highlight Color:")
        self.highlight_color_edit = QPushButtonColor(self.highlight_color_int)
        #-----------------------------------------------------------------------
        # closing
        self.show_button = QPushButton("Show")
        self.clear_button = QPushButton("Clear")
        self.close_button = QPushButton("Close")

    def create_layout(self):
        """displays the menu objects"""
        grid = QGridLayout()

        irow = 0
        grid.addWidget(self.nodes_label, irow, 0)
        grid.addWidget(self.nodes_edit, irow, 1)
        irow += 1

        grid.addWidget(self.elements_label, irow, 0)
        grid.addWidget(self.elements_edit, irow, 1)
        irow += 1

        grid.addWidget(self.highlight_color_label, irow, 0)
        grid.addWidget(self.highlight_color_edit, irow, 1)
        self.highlight_color_label.setVisible(False)
        self.highlight_color_edit.setVisible(False)
        irow += 1

        grid.addWidget(self.highlight_opacity_label, irow, 0)
        grid.addWidget(self.highlight_opacity_edit, irow, 1)
        self.highlight_opacity_label.setVisible(False)
        self.highlight_opacity_edit.setVisible(False)
        irow += 1

        #self.create_legend_widgets()
        #grid2 = self.create_legend_layout()
        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.show_button)
        ok_cancel_box.addWidget(self.clear_button)
        ok_cancel_box.addWidget(self.close_button)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        #vbox.addStretch()
        #vbox.addLayout(grid2)
        vbox.addStretch()

        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        self.highlight_color_edit.clicked.connect(self.on_highlight_color)
        self.highlight_opacity_edit.valueChanged.connect(self.on_highlight_opacity)
        self.nodes_edit.textChanged.connect(self.on_validate)
        self.elements_edit.textChanged.connect(self.on_validate)
        self.show_button.clicked.connect(self.on_show)
        self.clear_button.clicked.connect(self.on_remove_actors)
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
        """
        Choose a highlight color

        TODO: not implemented
        """
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
        """
        update the highlight opacity

        TODO: not implemented
        """
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
        """makes sure that all attributes are valid before doing any actions"""
        unused_nodes, flag1 = check_patran_syntax(self.nodes_edit, pound=self._nodes_pound)
        unused_elements, flag2 = check_patran_syntax(self.elements_edit, pound=self._elements_pound)
        if all([flag1, flag2]):
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_show(self):
        """show the highlight"""
        passed = self.on_validate()
        self.parent().mouse_actions.get_grid_selected(self.model_name)

        if passed and self.win_parent is not None:
            nodes, unused_flag1 = check_patran_syntax(self.nodes_edit, pound=self._nodes_pound)
            elements, unused_flag2 = check_patran_syntax(
                self.elements_edit, pound=self._elements_pound)
            if len(nodes) == 0 and len(elements) == 0:
                return False
            nodes_filtered = np.intersect1d(self.nodes, nodes)
            elements_filtered = np.intersect1d(self.elements, elements)
            nnodes = len(nodes_filtered)
            nelements = len(elements_filtered)
            if nnodes == 0 and nelements == 0:
                return False
            self.on_remove_actors()


            gui = self.parent()
            mouse_actions = gui.mouse_actions
            grid = mouse_actions.get_grid_selected(self.model_name)

            actors = create_highlighted_actors(
                gui, grid,
                all_nodes=self.nodes, nodes=nodes_filtered,
                all_elements=self.elements, elements=elements_filtered)

            if actors:
                add_actors_to_gui(gui, actors, render=True)
                self.actors = actors
        return passed

    def on_remove_actors(self):
        """removes multiple vtk actors"""
        gui = self.parent()
        if gui is not None:
            remove_actors_from_gui(gui, self.actors, render=True)
        self.actors = []

    def on_close(self):
        """close the window"""
        self.on_remove_actors()
        self.out_data['close'] = True
        self.close()


def create_highlighted_actors(gui, grid,
                              all_nodes=None, nodes=None,
                              all_elements=None, elements=None, add_actors=False):
    """creates nodes & element highlighted objects"""
    actors = []
    nnodes = 0
    nelements = 0
    if nodes is not None:
        nnodes = len(nodes)
        assert len(all_nodes) >= nnodes

    if elements is not None:
        nelements = len(elements)
        assert len(all_elements) >= nelements
    assert nnodes + nelements > 0

    if nnodes:
        point_ids = np.searchsorted(all_nodes, nodes)
        output_data = grid.GetPoints().GetData()
        points_array = vtk_to_numpy(output_data)  # yeah!

        point_array2 = points_array[point_ids, :]
        points2 = numpy_to_vtk_points(point_array2)

        ugrid = create_unstructured_point_grid(points2, nnodes)
        actor = create_highlighted_actor(gui, ugrid, representation='points', add_actor=add_actors)
        actors.append(actor)

    if nelements:
        cell_ids = np.searchsorted(all_elements, elements)

        selection_node = create_vtk_selection_node_by_cell_ids(cell_ids)
        ugrid = extract_selection_node_from_grid_to_ugrid(grid, selection_node)
        actor = create_highlighted_actor(gui, ugrid, representation='wire', add_actor=False)
        actors.append(actor)
    return actors

def check_float(cell):
    """validate the value is floatable"""
    text = cell.text()
    value = float(text)
    return value, True

#def check_label_float(cell):
    #text = cell.text()
    #try:
        #value = eval_float_from_string(text)
        #cell.setStyleSheet("QLineEdit{background: white;}")
        #return value, True
    #except ValueError:
        #cell.setStyleSheet("QLineEdit{background: red;}")
        #return None, False

def _check_color(color_float):
    assert len(color_float) == 3, color_float
    assert isinstance(color_float[0], float), color_float
    color_int = [int(colori * 255) for colori in color_float]
    return color_float, color_int

def main():  # pragma: no cover
    """basic testing"""
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
