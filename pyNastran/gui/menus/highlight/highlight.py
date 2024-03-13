"""
The highlight menu handles:
 - Nodes/Elements
 - Preferences

"""
from __future__ import annotations
from typing import Any, TYPE_CHECKING
import numpy as np

from qtpy import QtGui
from qtpy.QtWidgets import (
    QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
    QSpinBox, QDoubleSpinBox, QColorDialog) # QCheckBox

#from vtk import (
    #vtkLODActor,
    #vtkLabeledDataMapper,
    #vtkCellCenters, vtkIdFilter,
    #vtkUnstructuredGridGeometryFilter,
    #vtkVertexGlyphFilter,
#)
from vtkmodules.vtkRenderingLOD import vtkLODActor

from pyNastran.gui.vtk_rendering_core import vtkActor2D, vtkRenderer
from pyNastran.gui.utils.qt.pydialog import PyDialog, make_font, check_patran_syntax, check_color
from pyNastran.gui.utils.qt.qpush_button_color import QPushButtonColor
#from pyNastran.gui.menus.menu_utils import eval_float_from_string
from pyNastran.gui.utils.qt.qelement_edit import (
    QNodeLineEdit, QElementLineEdit,
    QNodeTextEdit, QElementTextEdit)

from pyNastran.gui.gui_objects.settings import Settings, ANNOTATION_SIZE
from pyNastran.gui.utils.vtk.gui_utils import add_actors_to_gui, remove_actors_from_gui
from pyNastran.gui.menus.highlight.vtk_utils import (
    create_highlighted_actors, create_node_label_actor, create_cell_label_actor,
    get_ids_filter)
if TYPE_CHECKING:
    from pyNastran.gui.menus.groups_modify.groups import Group

USE_LINEEDIT = True

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
    def __init__(self, data: dict[str, Any], menu_type: str, win_parent=None):
        """
        Saves the data members from data and
        performs type checks

        Parameters
        ----------
        menu_type : str
            'highlight'
            'mark'

        """
        PyDialog.__init__(self, data, win_parent)
        gui = win_parent
        self.menu_type = menu_type

        self._default_annotation_size = ANNOTATION_SIZE
        if gui is None:  # pragma: no cover
            self.highlight_color_float = [0., 0., 0.]
            self.highlight_color_int = [0, 0, 0]
            self._highlight_opacity = 0.9
            self._point_size = 10
            self._label_size = 10.0
        else:
            settings: Settings = gui.settings
            self.highlight_color_float, self.highlight_color_int = check_color(
                settings.highlight_color)

            self._highlight_opacity = settings.highlight_opacity

            self._point_size = settings.highlight_point_size
            self._line_width = settings.highlight_line_width

            self._point_size = 10
            self._annotation_size = settings.annotation_size
        self._updated_window = False

        self.actors: list[vtkActor2D] = []
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

        if self.menu_type == 'highlight':
            self.setWindowTitle('Highlight')
        elif self.menu_type == 'mark':
            self.setWindowTitle('Mark')
        else:  # pragma: no cover
            raise NotImplementedError(self.menu_type)

        self.create_widgets()
        self.create_layout()
        self.set_connections()
        self.on_font(self._default_font_size)
        #self.show()

    def _create_node_element_edits(self, model_name: str) -> None:
        if USE_LINEEDIT:
            MAX_LENGTH = 100_000
            self.nodes_edit = QNodeLineEdit(
                self, model_name, pick_style='area',
                cleanup=True, tab_to_next=False, max_length=MAX_LENGTH)

            self.elements_edit = QElementLineEdit(
                self, model_name, pick_style='area',
                cleanup=True, tab_to_next=False, max_length=MAX_LENGTH)
        else:
            self.nodes_edit = QNodeTextEdit(
                self, model_name, pick_style='area',
                cleanup=True, tab_to_next=False)

            self.elements_edit = QElementTextEdit(
                self, model_name, pick_style='area',
                cleanup=True, tab_to_next=False)

    def create_widgets(self):
        """creates the display window"""
        # Text Size

        model_name = self.model_name
        self.nodes_label = QLabel('Nodes:')
        self.elements_label = QLabel('Elements:')
        self._create_node_element_edits(model_name)

        #-----------------------------------------------------------------------
        # Highlight Color
        if self.menu_type == 'highlight':
            self.highlight_opacity_label = QLabel('Highlight Opacity:')
            self.highlight_opacity_edit = QDoubleSpinBox(self)
            self.highlight_opacity_edit.setValue(self._highlight_opacity)
            self.highlight_opacity_edit.setRange(0.1, 1.0)
            self.highlight_opacity_edit.setDecimals(2)
            self.highlight_opacity_edit.setSingleStep(0.05)
            self.highlight_opacity_button = QPushButton('Default')

            # Text Color
            self.highlight_color_label = QLabel('Highlight Color:')
            self.highlight_color_edit = QPushButtonColor(self.highlight_color_int)

            self.point_size_label = QLabel('Point Size:')
            self.point_size_edit = QSpinBox(self)
            self.point_size_edit.setValue(self._point_size)
            self.point_size_edit.setRange(7, 30)
            #self.point_size_button = QPushButton("Default")
        elif self.menu_type == 'mark':
            self.label_size_label = QLabel('Label Size:')
            self.label_size_edit = QSpinBox(self)
            self.label_size_edit.setValue(self._default_annotation_size)
            self.label_size_edit.setRange(7, 30)
            #self.label_size_button = QPushButton("Default")
        else:
            raise NotImplementedError(self.menu_type)
        #-----------------------------------------------------------------------
        # closing
        #if self.menu_type == 'highlight':
        self.show_button = QPushButton('Show')
        #elif self.menu_type == 'mark':
            #self.mark_button = QCheckBox('Mark')
        #else:
            #raise NotImplementedError(self.menu_type)
        self.clear_button = QPushButton('Clear')
        self.close_button = QPushButton('Close')

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

        if self.menu_type == 'highlight':
            # TODO: enable me
            grid.addWidget(self.highlight_color_label, irow, 0)
            grid.addWidget(self.highlight_color_edit, irow, 1)
            self.highlight_color_label.setEnabled(False)
            self.highlight_color_edit.setEnabled(False)
            irow += 1

            # TODO: enable me
            grid.addWidget(self.highlight_opacity_label, irow, 0)
            grid.addWidget(self.highlight_opacity_edit, irow, 1)
            self.highlight_opacity_label.setEnabled(False)
            self.highlight_opacity_edit.setEnabled(False)
            irow += 1

            # TODO: enable me
            grid.addWidget(self.point_size_label, irow, 0)
            grid.addWidget(self.point_size_edit, irow, 1)
            self.point_size_label.setEnabled(False)
            self.point_size_edit.setEnabled(False)
            irow += 1
        elif self.menu_type == 'mark':
            # TODO: enable me
            grid.addWidget(self.label_size_label, irow, 0)
            grid.addWidget(self.label_size_edit, irow, 1)
            #self.label_size_label.setEnabled(False)
            #self.label_size_edit.setEnabled(False)
            irow += 1
            #self.mark_button.setEnabled(False)
        else:  # pragma: no cover
            raise NotImplementedError(self.menu_type)

        #self.create_legend_widgets()
        #grid2 = self.create_legend_layout()
        ok_cancel_box = QHBoxLayout()
        #if self.menu_type == 'highlight':
        ok_cancel_box.addWidget(self.show_button)
        #elif self.menu_type == 'mark':
            #ok_cancel_box.addWidget(self.mark_button)
        #else:
            #raise NotImplementedError(self.menu_type)
        ok_cancel_box.addWidget(self.clear_button)
        ok_cancel_box.addWidget(self.close_button)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        #vbox.addStretch()
        #vbox.addLayout(grid2)
        if USE_LINEEDIT:
            vbox.addStretch()

        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        if self.menu_type == 'highlight':
            self._set_connections_highlight()
        elif self.menu_type == 'mark':
            self._set_connections_mark()
        else:  # pragma: no cover
            raise NotImplementedError(self.menu_type)
        self._set_connections_end()

    def _set_connections_highlight(self):
        """creates the actions for the menu"""
        self.highlight_color_edit.clicked.connect(self.on_highlight_color)
        self.highlight_opacity_edit.valueChanged.connect(self.on_highlight_opacity)
        self.show_button.clicked.connect(self.on_show)

    def _set_connections_mark(self) -> None:
        """creates the actions for the menu"""
        self.show_button.clicked.connect(self.on_show)

    def _set_connections_end(self) -> None:
        """creates the actions for the menu"""
        self.nodes_edit.textChanged.connect(self.on_validate)
        self.elements_edit.textChanged.connect(self.on_validate)
        self.clear_button.clicked.connect(self.on_remove_actors)
        self.close_button.clicked.connect(self.on_close)

    def on_font(self, value=None) -> None:
        """update the font for the current window"""
        if value in (0, None):
            value = self.font_size_edit.value()
        font = make_font(value, is_bold=False)
        self.setFont(font)

    def on_highlight_color(self) -> None:
        """
        Choose a highlight color

        TODO: not implemented
        """
        title = "Choose a highlight color"
        rgb_color_ints = self.highlight_color_int
        color_edit = self.highlight_color_edit
        func_name = 'set_highlight_color'
        passed, rgb_color_ints, rgb_color_floats = create_color_menu(
            self, self.win_parent, title, color_edit, rgb_color_ints, func_name)
        if passed:
            self.highlight_color_int = rgb_color_ints
            self.highlight_color_float = rgb_color_floats

    def on_highlight_opacity(self, value=None) -> None:
        """
        update the highlight opacity

        TODO: not implemented
        """
        if value is None:
            value = self.highlight_opacity_edit.value()
        self._highlight_opacity = value
        if self.win_parent is not None:
            self.win_parent.settings.set_highlight_opacity(value)

    #---------------------------------------------------------------------------

    def on_validate(self) -> bool:
        """makes sure that all attributes are valid before doing any actions"""
        unused_nodes, flag1 = check_patran_syntax(self.nodes_edit, pound=self._nodes_pound)
        unused_elements, flag2 = check_patran_syntax(self.elements_edit, pound=self._elements_pound)
        if all([flag1, flag2]):
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_show(self) -> bool:
        """show the highlight"""
        passed = self.on_validate()
        if not passed or self.win_parent is None:
            return passed
        self.parent().mouse_actions.get_grid_selected(self.model_name)

        nodes, unused_flag1 = check_patran_syntax(self.nodes_edit, pound=self._nodes_pound)
        elements, unused_flag2 = check_patran_syntax(
            self.elements_edit, pound=self._elements_pound)
        if len(nodes) == 0 and len(elements) == 0:
            return False

        if self.win_parent is not None:
            gui = self.win_parent
            group: Group = gui.groups[gui.group_active]
            all_nodes = group.node_ids
            all_elements = group.element_ids
        nodes_filtered = np.intersect1d(all_nodes, nodes)
        elements_filtered = np.intersect1d(all_elements, elements)
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
            all_nodes=all_nodes, nodes=nodes_filtered, set_node_scalars=True,
            all_elements=all_elements, elements=elements_filtered, set_element_scalars=True,
            add_actors=False)

        #make_highlight = self.menu_type == 'highlight'
        make_labels = self.menu_type == 'mark'
        make_element_labels = True
        make_node_labels = True

        if make_labels:
            rend = gui.rend
            actors = create_mark_actors(
                rend, make_node_labels, make_element_labels,
                nnodes, nelements, actors, label_size=self._annotation_size)

        if actors:
            add_actors_to_gui(gui, actors, render=True)
            self.actors = actors
        gui.Render()
        return passed

    def on_remove_actors(self):
        """removes multiple vtk actors"""
        gui = self.parent()
        if len(self.actors) == 0:
            return
        if gui is not None:
            if self.nodes_edit.style is not None:
                self.nodes_edit.style.remove_actors()
            if self.elements_edit.style is not None:
                self.elements_edit.style.remove_actors()
            remove_actors_from_gui(gui, self.actors, render=True, force_render=True)
            gui.Render()
        self.actors = []

    def closeEvent(self, unused_event) -> None:
        """close the window"""
        self.on_close()

    def on_close(self) -> None:
        """close the window"""
        self.on_remove_actors()
        self.out_data['close'] = True
        self.close()


def create_mark_actors(rend: vtkRenderer,
                       make_node_labels: bool, make_element_labels: bool,
                       nnodes: int, nelements: int,
                       actors: list[vtkLODActor],
                       label_size: float=10.0) -> list[vtkActor2D]:
    """Replace the actors with labels

    For one actor (e.g., the main model), you may want to show many node/element ids.
    Thus, you'll make 1 label_actor per actor.
    """
    iactor = 0
    actors2 = []
    if make_node_labels and nnodes:
        actor = actors[iactor]
        label_actor = create_node_label_actor(actor, rend, label_size=label_size)
        #actors.append(label_actor)
        actors2.append(label_actor)
        iactor += 1

    if make_element_labels and nelements:
        actor = actors[iactor]
        label_actor = create_cell_label_actor(actor, rend, label_size=label_size)
        #actors.append(label_actor)
        actors2.append(label_actor)
        iactor += 1
    return actors2

def check_float(cell) -> tuple[float, bool]:
    """validate the value is floatable"""
    text = cell.text()
    value = float(text)
    return value, True

def create_color_menu(parent, win_parent, title: str,
                      color_edit: QPushButtonColor,
                      rgb_color_ints: list[int],
                      func_name: str):
    """helper method for ``on_background_color`` and ``on_background_color2``"""
    passed, rgb_color_ints, rgb_color_floats = _pop_color_dialog(
        parent, color_edit, rgb_color_ints, title)
    if passed:
        if win_parent is not None:
            settings = win_parent.settings
            func_background_color = getattr(settings, func_name)
            func_background_color(rgb_color_floats)
    return passed, rgb_color_ints, rgb_color_floats

def _pop_color_dialog(parent,
                      color_edit: QPushButtonColor,
                      rgb_color_ints: list[int], title: str) -> tuple[bool, list[int], list[float]]:
    """pops a color dialog"""
    col = QColorDialog.getColor(QtGui.QColor(*rgb_color_ints), parent, title)
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


#def check_label_float(cell):
    #text = cell.text()
    #try:
        #value = eval_float_from_string(text)
        #cell.setStyleSheet(QLINEEDIT_GOOD)
        #return value, True
    #except ValueError:
        #cell.setStyleSheet(QLINEEDIT_ERROR)
        #return None, False

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
    if 1:
        main_window = HighlightWindow(data, menu_type='highlight')
    else:
        main_window = HighlightWindow(data, menu_type='mark')
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":  # pragma: no cover
    main()
