import sys
from typing import List, Dict, Union, Any

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QPushButton, QApplication, QGridLayout,
    QComboBox, QLabel, QSpinBox, QLineEdit,
)
#from qtpy import QtCore
from pyNastran.gui.utils.qt.results_window import ResultsWindow
from pyNastran.gui.gui_objects.gui_result import GuiResult, NullResult
from pyNastran.gui.utils.qt.utils import (
    add_obj_to_vbox,
    create_hbox_with_widgets,
    add_line_widgets_to_grid)
from pyNastran.gui.utils.utils import find_next_value_in_sorted_list

SkippableSpinBox = QSpinBox
#class SkippableSpinBox(QSpinBox):
    #stepChanged = QtCore.Signal() #  tested in PySide2, pyqtSignal?

    #def stepBy(self, step):
        #keys = self.parent().case_keys
        #old = self.value()
        #new = old + step
        #find_next_value_in_sorted_list(keys, old, new)
        #super(SkippableSpinBox, self).stepBy(step)
        #if self.value() != value:
            #self.stepChanged.emit()

class Sidebar(QWidget):
    """
    +----------------------+
    |        Results       |
    +======================+
    |                      |
    |  Name = Main         |
    |                      |
    |  +----------------+  |
    |  | ResultsWindow  |  |
    |  +----------------+  |
    |  |                |  |
    |  |  +----------+  |  |
    |  |  | - a      |  |  |
    |  |  |  - b1    |  |  |
    |  |  |  - b2    |  |  |
    |  |  |  - b3    |  |  |
    |  |  +----------+  |  |
    |  |                |  |
    |  +----------------+  |
    |                      |
    |                      |
    |         Apply        |
    +----------------------+

    For Nastran:
      - a: Subcase 1
       - b1. Displacement
       - b2. Stress
       - b3. Strain

    For Cart3d:
      - a1. Geometry
       - b1. ElementID
       - b2. Region
      - a2. Results Case 1
       - b1. U
       - b2. V
       - b3. W

    +--------------+
    |  Sub-Result  | (pulldown)
    +==============+
    | - a1         |
    | - a2         |
    | - a3         |
    | - a4         |
    +--------------+

    For Nastran:
      - a1: Displacement X
      - a2. Displacement Y
      - a3. Displacmenet Z
      - a4. Displacmenet Mag

    For Cart3d:
      - NA (Greyed Out)

    +----------------+
    |     Plot       | (pulldown)
    +================+
    | - Fringe       |
    | - Marker       |
    | - Displacement |
    +----------------+
    - Cart3d -> Fringe (disabled)


    +---------------+
    | Scale Display | (text box)
    +===============+
    | 0 < x < 1000  | (not for fringe)
    +---------------+

    +--------------+
    |   Location   | (pulldown)
    +==============+
    | - nodal      |
    | - centroidal |
    +--------------+
    (disabled)

    +------------------+
    |  Complex Method  | (pulldown)
    +==================+
    | - real           | (usually set to real and disabled)
    | - imag           |
    | - mag            |
    | - phase          |
    | - max over phase |
    +------------------+

    +--------------+
    |    Derive    | (pulldown; only for nodal results)
    +==============+
    | - derive/avg |  (default?)
    | - avg/derive |
    +--------------+
    """
    def __init__(self, parent, debug=False, data=None, clear_data=True, name='main',
                 setup_dict=None,
                 # buttons
                 include_case_spinner=False,
                 include_deflection_scale=False,
                 include_vector_scale=False,
                 # actions
                 include_clear=True,
                 include_export_case=False,
                 include_delete=True,
                 include_results=True):
        """
        Creates the buttons in the Sidebar, not the actual layout

        Parameters
        ----------
        parent : MainWindow()
            the gui
        debug : bool; default=False
            flag for debug info
        data : List[tree]
            the tree
        clear_data : bool; default=True
            ???
        name : str; default='main'
            the active name
        setup_dict : Dict[irow] = List[QWidgets]
            a way to add additional widgets to the sidebar

        """
        #include_case_spinner = False
        QWidget.__init__(self)
        self.parent = parent
        self.debug = debug
        self.setup_dict = setup_dict
        self._update_case = True
        self.case_keys = []
        self.icase = 0  # default

        # buttons
        self.include_case_spinner = include_case_spinner
        self.include_deflection_scale = include_deflection_scale
        self.include_vector_scale = include_vector_scale


        choices = ['keys2', 'purse2', 'cellphone2', 'credit_card2', 'money2']
        if data is None:
            data = []

        self.result_case_windows = [
            ResultsWindow(self, 'Case/Results', data, choices,
                          include_clear=include_clear,
                          include_export_case=include_export_case,
                          include_delete=include_delete,
                          include_results=include_results)
        ]
        data = [
            ('A', 1, []),
            #('B', 2, []),
            #('C', 3, []),
        ]
        self.result_method_window = ResultsWindow(self, 'Method', data, choices)
        self.result_method_window.setVisible(False)
        #else:
            #self.result_method_window = None

        self.show_pulldown = False
        if self.show_pulldown:
            #combo_options = ['a1', 'a2', 'a3']
            self.pulldown = QComboBox()
            self.pulldown.addItems(choices)
            self.pulldown.activated[str].connect(self.on_pulldown)

        self.apply_button = QPushButton('Apply', self)
        self.apply_button.clicked.connect(self.on_apply)

        if name is None:
            self.name = None
            self.names = ['N/A']
            name = 'N/A'
        else:
            self.name = str(name)
            self.names = [name]

        self.name_label = QLabel("Name:")
        self.name_pulldown = QComboBox()
        self.name_pulldown.addItem(name)
        self.name_pulldown.setDisabled(True)

        if include_case_spinner:
            self.case_spinner_label = QLabel('Case:')
            self.case_spinner = SkippableSpinBox()
            self.case_spinner_label.setVisible(False)
            self.case_spinner.setVisible(False)
            self.case_spinner.lineEdit().setReadOnly(True)

            # -1 is actually invalid, but we'll correct it later
            self.case_spinner.setMinimum(-1)
            if self.has_cases:
                self.set_max_case(self.parent.result_cases)
        if include_deflection_scale:
            self.deflection_label = QLabel('Deflection Scale:')
            self.deflection_edit = QLineEdit()
        if include_vector_scale:
            self.vector_label = QLabel('Vector Scale:')
            self.vector_edit = QLineEdit()
        #if include_vector:

        self.setup_layout(data, choices, clear_data=clear_data)
        self.set_connections()

    def set_connections(self):
        """creates the actions for the menu"""
        self.name_pulldown.currentIndexChanged.connect(self.on_update_name)
        if self.include_case_spinner:
            self.case_spinner.valueChanged.connect(self._on_case)
        #if self.include_deflection_scale:
            #self.deflection_edit.valueChanged.connect(self.on_deflection_scale)
        #if self.include_vector_scale:
            #self.vector_scale.valueChanged.connect(self.on_vector_scale)

    def set_max_case(self, cases: Union[List[int], Dict[int, Any]]):
        """
        The max case id needs to be dynamic because additional results
        can be added

        """
        if self.has_cases:
            if len(cases) == 0:
                return
            self.case_spinner_label.setVisible(True)
            self.case_spinner.setVisible(True)
            # we add the +1, so we can wrap around
            self.case_spinner.setMaximum(max(cases) + 1)
            self._on_case()

    def set_case_keys(self, case_keys: List[int]):
        """set the availiable keys for the case spinner"""
        if not self.include_case_spinner:
            return
        self.case_keys = case_keys
        if self.icase == -1:
            self.icase = case_keys[0]
        self.set_max_case(case_keys)

    def update_icase(self, icase_frige):
        """callback for updating the case spinner"""
        if not self.include_case_spinner:
            return
        self._update_case = False
        self.case_spinner.setValue(icase_frige)
        self._update_case = True

    def _on_case(self):
        """callback for updating the GUI"""
        if not self._update_case:
            return
        icase = self.case_spinner.value()

        next_value = find_next_value_in_sorted_list(self.case_keys, self.icase, icase)
        if icase != next_value:
            self._update_case = False
            self.case_spinner.setValue(next_value)
        self._set_case(next_value)
        self._update_case = True

    def _set_case(self, icase):
        #print(f'changing from icase={self.icase} -> {icase} (next_value={next_value})')
        #icase = next_value
        self.icase = icase
        if self.has_cases:
            result_name = None
            self.parent._set_case(result_name, icase, explicit=True)
            (obj, (i, name)) = self.parent.result_cases[icase]

            deflection_is_visible = False
            vector_is_visible = False
            if isinstance(obj, (GuiResult, NullResult)):
                self._set_buttons(deflection_is_visible, vector_is_visible)
                return
            #DisplacementResults
            #elif isinstance(obj, ForceTableResults):
                # include_vector_scale
            #print(i, name)
            #print(obj)
            #if self.include_deflection_scale:
                #case =
            #if self.include_vector_scale:
                #case

    def _set_buttons(self, deflection_is_visible, vector_is_visible):
        """show/hide the additional buttons"""
        if self.include_deflection_scale:
            self.deflection_label.setVisible(deflection_is_visible)
            self.deflection_edit.setVisible(deflection_is_visible)
        if self.include_vector_scale:
            self.vector_label.setVisible(vector_is_visible)
            self.vector_edit.setVisible(vector_is_visible)

    #def on_deflection_scale(self):
        #pass
    #def on_vector_scale(self):
        #pass

    @property
    def result_case_window(self):
        if self.name is None:
            i = 0
        else:
            i = self.names.index(self.name)
        i = 0
        return self.result_case_windows[i]

    def setup_layout(self, data, choices, clear_data=True, init=True):
        """creates the sidebar visual layout"""
        #if not init:
            #self.frameGeometry().
            #width = self.frameGeometry().width()
            #height = self.frameGeometry().height()
            #print('width=%s height=%s' % (width, height))

        #print('init...')
        vbox = QVBoxLayout()

        irow = 0
        self._add_from_setup_dict(vbox, irow)

        hbox = create_hbox_with_widgets([self.name_label, self.name_pulldown])
        vbox.addLayout(hbox)

        irow += 1
        self._add_from_setup_dict(vbox, irow)

        nwindows = len(self.result_case_windows)
        #print('nwindows=%s self.names=%s' % (nwindows, self.names))
        for i in range(nwindows):
            #print('*using existing window')
            result_case_window = self.result_case_windows[i]
            vbox.addWidget(result_case_window)
            #result_case_window.setVisible(False)  # be very careful of this...

        nwindows = len(self.result_case_windows)
        for name in self.names[nwindows:]:
            #print('*creating a window')
            result_case_window = ResultsWindow(self, 'Case/Results', data, choices)
            result_case_window.setVisible(False)
            vbox.addWidget(result_case_window)
            self.result_case_windows.append(result_case_window)

        iname = 0
        #if self.name is None:
            #iname = 0
        #else:
            #iname = self.names.index(self.name)
        #for i in range(nwindows):
            #if i != iname:
                #self.result_case_windows[iname].setVisible(False)
        #self.result_case_windows[iname].setVisible(True)

        irow += 1
        self._add_from_setup_dict(vbox, irow)

        if self.result_method_window:
            vbox.addWidget(self.result_method_window)
        if self.show_pulldown:
            vbox.addWidget(self.pulldown)

        irow += 1
        self._add_from_setup_dict(vbox, irow)

        self._add_grid_to_vbox(vbox)

        vbox.addWidget(self.apply_button)

        irow += 1
        self._add_from_setup_dict(vbox, irow)

        self.setLayout(vbox)

        if clear_data:
            self.clear_data()

        #if not init:
            #self.frameGeometry().width()
            #self.frameGeometry().height()
            #self.resize(width, height)

    def _add_grid_to_vbox(self, vbox):
        """creates the grid for some optional buttons"""
        if self.include_case_spinner or self.include_deflection_scale or self.include_vector_scale:
            grid = QGridLayout()
            irow = 0

            if self.include_case_spinner:
                irow = add_line_widgets_to_grid(
                    grid, irow, [self.case_spinner_label, self.case_spinner])
            if self.include_deflection_scale:
                irow = add_line_widgets_to_grid(
                    grid, irow, [self.deflection_label, self.deflection_edit])
            if self.include_vector_scale:
                irow = add_line_widgets_to_grid(
                    grid, irow, [self.vector_label, self.vector_edit])
            vbox.addLayout(grid)

    def _add_from_setup_dict(self, vbox, irow):
        if self.setup_dict is None:
            return

        if irow in self.setup_dict:
            widgets = self.setup_dict[irow]
            assert widgets is not None, widgets
            if isinstance(widgets, list):
                for widget_layout in widgets:
                    add_obj_to_vbox(vbox, widget_layout)
            else:
                # scalar
                add_obj_to_vbox(vbox, widgets)

    def update_method(self, method):
        if isinstance(method, str):
            datai = self.result_method_window.data[0]
            self.result_method_window.data[0] = (method, datai[1], datai[2])
            #print('method=%s datai=%s' % (method, datai))
            self.result_method_window.update_data(self.result_method_window.data)
        else:
            return
             # pragma: no cover
            #datai = self.result_method_window.data[0]

    def get_form(self):
        """
        TODO: At this point, we should clear out the data block and refresh it
        """
        return self.result_case_window.data

    def update_results(self, data, name):
        """
        Updates the sidebar

        Parameters
        ----------
        data : List[tuple]
            the form data
        name : str
            the name that goes at the side
        """
        name = str(name)
        update_name = False
        setup_layout = False
        if self.name is None:
            self.names = [name]
            self.name_pulldown.clear()
            self.name_pulldown.addItems(self.names)
            #self.name_pulldown.setItemText(0, name)
            #self.name_pulldown.setCurrentIndex(1)
            update_name = True
            setup_layout = True

        else:
            if name in self.names:
                i = self.names.index(name)
                self.name_pulldown.setCurrentIndex(i)
            else:
                self.name_pulldown.addItem(name)
                self.names.append(name)
                setup_layout = True

            if len(self.names) >= 2:
                self.name_pulldown.setEnabled(True)
        self.name = name

        self.result_case_window.update_data(data)
        self.apply_button.setEnabled(True)
        if update_name:
            self.on_update_name(None)

        if setup_layout and 0:  # pragma: no cover
            #print('setup_layout******')
            ## TODO: screws up the width of the window
            choices = ['keys2', 'purse2', 'cellphone2', 'credit_card2', 'money2']
            data = [
                ('A', 1, []),
                #('B', 2, []),
                #('C', 3, []),
            ]
            self.setup_layout(data, choices, init=False, clear_data=True)

    def update_methods(self, data):
        """the methods is a hidden box"""
        if self.result_method_window is not None:
            self.result_method_window.update_data(data)
        self.apply_button.setEnabled(True)

    def clear_data(self):
        self.result_case_window.clear_data()
        if self.result_method_window is not None:
            self.result_method_window.clear_data()
        self.apply_button.setEnabled(False)

    def on_pulldown(self):  # pragma: no cover
        print('pulldown...')

    def on_update_name(self):  # pragma: no cover
        """user clicked the pulldown"""
        name = str(self.name_pulldown.currentText())
        data = self.parent._get_sidebar_data(name)
        #self.result_case_window.update_data(data)

    def on_apply(self):
        data = self.result_case_window.data
        valid_a, keys_a = self.result_case_window.treeView.get_row()

        data = self.result_method_window.data
        valid_b, keys_b = self.result_method_window.treeView.get_row()
        if valid_a and valid_b:
            if self.debug:  # pragma: no cover
                print('  rows1 = %s' % self.result_case_window.treeView.old_rows)
                print('        = %s' % str(keys_a))
                print('  rows2 = %s' % self.result_method_window.treeView.old_rows)
                print('        = %s' % str(keys_b))
            else:
                self.update_vtk_window(keys_a, keys_b)

    def update_vtk_window(self, keys_a, keys_b):
        if 0:  # pragma: no cover
            print('keys_a = %s' % str(keys_a))
            for i, key in enumerate(self.parent.case_keys):
                if key[1] == keys_a[0]:
                    break
            print('*i=%s key=%s' % (i, str(key)))
            #self.parent.update_vtk_window_by_key(i)
            #result_name = key[1]
            #self.parent.cycle_results_explicit(result_name=result_name, explicit=True)
            #j = self.parent._get_icase(result_name)
            #j = i
        i = keys_a

        #self.case_spinner.setValue(i)  # this might just work?

        # set the spinner, but don't take any actions
        if 0:
            if self._update_case:
                self._update_case = False
                self.case_spinner.setValue(i)
                self._update_case = True
            else:
                self.case_spinner.setValue(i)
        result_name = None
        #self._set_case(i)
        self.parent._set_case(result_name, i, explicit=True)

    @property
    def has_cases(self):
        return hasattr(self.parent, 'result_cases')

    def on_clear_results(self):
        if hasattr(self.parent, 'on_clear_results'):
            self.parent.on_clear_results()
    def on_fringe(self, icase):
        if hasattr(self.parent, 'on_fringe'):
            self.parent.on_fringe(icase)
    def on_disp(self, icase):
        if hasattr(self.parent, 'on_disp'):
            self.parent.on_disp(icase)
    def on_vector(self, icase):
        if hasattr(self.parent, 'on_vector'):
            self.parent.on_vector(icase)

    def get_clicked(self):
        self.result_method_window.treeView.get_clicked()


def main():  # pragma: no cover
    app = QApplication(sys.argv)

    form = [
        [u'Geometry', None, [
            (u'NodeID', 0, []),
            (u'ElementID', 1, []),
            (u'PropertyID', 2, []),
            (u'MaterialID', 3, []),
            (u'E', 4, []),
            (u'Element Checks', None, [
                (u'ElementDim', 5, []),
                (u'Min Edge Length', 6, []),
                (u'Min Interior Angle', 7, []),
                (u'Max Interior Angle', 8, [])],
            ),],
        ],
    ]
    #form = []
    res_widget = Sidebar(app, data=form, clear_data=False, debug=True)

    name = 'name'
    #res_widget.update_results(form, name)

    res_widget.show()
    sys.exit(app.exec_())

if __name__ == "__main__":  # pragma: no cover
    main()
