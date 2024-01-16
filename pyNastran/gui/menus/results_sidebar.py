import sys
from typing import Union, Any

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QPushButton, QApplication, QGridLayout,
    QComboBox, QLabel, QSpinBox, QLineEdit,
    QSplitter, QStyleFactory,
)
#from qtpy.QtCore import QSize
from qtpy import QtCore
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

SHOW_DEV = False
USE_NEW_SIDEBAR = False

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
    def __init__(self, parent,
                 debug: bool=False,
                 data=None, clear_data: bool=True,
                 name: str='main',
                 left_click_callback=None,
                 # actions
                 include_clear: bool=True,
                 include_export_case: bool=False,
                 include_delete: bool=True,
                 include_results: bool=True):
        """
        Creates the buttons in the Sidebar, not the actual layout

        Parameters
        ----------
        parent : MainWindow()
            the gui
        debug : bool; default=False
            flag for debug info
        data : list[tree]
            the tree
        clear_data : bool; default=True
            ???
        name : str; default='main'
            the active name
        setup_dict : dict[irow] = list[QWidgets]
            a way to add additional widgets to the sidebar

        """
        #include_case_spinner = False
        QWidget.__init__(self)
        self.parent = parent
        self.debug = debug
        self._update_case = True
        self.case_keys = []
        self.icase = 0  # default

        choices = ['keys2', 'purse2', 'cellphone2', 'credit_card2', 'money2']
        if data is None:
            data = []

        self.result_case_windows = [
            ResultsWindow(self, 'Case/Results', data, choices,
                          left_click_callback=left_click_callback,
                          right_click_actions=[],
                          is_single_select=True,
                          include_clear=include_clear,
                          include_export_case=include_export_case,
                          include_delete=include_delete,
                          include_results=include_results)
        ]
        data = [
            ('Top', 1, []),
            ('Bottom', 2, []),
            #('C', 3, []),
        ]

        if USE_NEW_SIDEBAR:
            self.result_method_window = ResultsWindow(
                self, 'Method', data, choices,
                left_click_callback=None,
                right_click_actions=None,
                is_single_select=False,
                include_clear=False,
                include_export_case=False,
                include_delete=False,
                include_results=False)
            self.result_method_window.setVisible(True)

        self.show_pulldown = False
        if self.show_pulldown:
            #combo_options = ['a1', 'a2', 'a3']
            self.pulldown = QComboBox()
            self.pulldown.addItems(choices)
            self.pulldown.activated[str].connect(self.on_pulldown)

        self.apply_button = QPushButton('Apply', self)
        self.apply_button.clicked.connect(self.on_apply)

        assert name is not None, name
        self.name = str(name)
        self.names = [name]

        self.name_label = QLabel("Name:")
        self.name_pulldown = QComboBox()
        self.name_pulldown.addItem(name)
        self.name_pulldown.setDisabled(True)

        self.sub_result_name_label = QLabel("Name:")
        self.sub_result_name_pulldown = QComboBox()
        #self.sub_result_name_pulldown.addItem(name)
        #self.sub_result_name_pulldown.setDisabled(True)

        self.case_spinner_label = QLabel('Case:')
        self.case_spinner = SkippableSpinBox()
        self.case_spinner.lineEdit().setReadOnly(True)

        # -1 is actually invalid, but we'll correct it later
        self.case_spinner.setMinimum(-1)
        if self.has_cases:
            self.set_max_case(self.parent.result_cases)

        if USE_NEW_SIDEBAR:
            self.transform_coords_label = QLabel('2: Transform:')
            self.transform_coords_pulldown = QComboBox()
            self.transform_coords_pulldown.addItem('Coord 0')

            self.min_max_average_label = QLabel('3: Derivation Method:')
            self.min_max_average_pulldown = QComboBox()
            self.min_max_average_pulldown.addItems(['Derive/Average', 'Max', 'Min'])

            self.combine_label = QLabel('4: Nodal combine:')
            self.combine_pulldown = QComboBox(self)
            self.combine_pulldown.addItems(['Across Neighbors', 'Centroid Absolute Max'])
            self.combine_pulldown.setEnabled(False)

        self.deflection_label = QLabel('Deflection Scale:')
        self.deflection_edit = QLineEdit()

        self.vector_label = QLabel('Vector Scale:')
        self.vector_edit = QLineEdit()

        self.setup_layout(data, choices, clear_data=clear_data)
        self.set_connections()
        if not SHOW_DEV:
            self.hide_dev()

    def set_coord_transform_visible(self, is_visible: bool) -> None:
        if USE_NEW_SIDEBAR:
            self.transform_coords_label.setVisible(is_visible)
            self.transform_coords_pulldown.setVisible(is_visible)

    def set_derivation_visible(self, is_visible: bool) -> None:
        if USE_NEW_SIDEBAR:
            self.min_max_average_label.setVisible(is_visible)
            self.min_max_average_pulldown.setVisible(is_visible)

    def set_nodal_combine_visible(self, is_visible: bool) -> None:
        if USE_NEW_SIDEBAR:
            self.combine_label.setVisible(is_visible)
            self.combine_pulldown.setVisible(is_visible)

    def hide_dev(self):
        objs = [
            self.case_spinner_label, self.case_spinner,
            self.deflection_label, self.deflection_edit,
            self.vector_label, self.vector_edit,
        ]
        for obj in objs:
            obj.hide()

    def set_connections(self):
        """creates the actions for the menu"""
        self.name_pulldown.currentIndexChanged.connect(self.on_update_name)
        self.case_spinner.valueChanged.connect(self._on_case)
        #self.deflection_edit.valueChanged.connect(self.on_deflection_scale)
        #self.vector_scale.valueChanged.connect(self.on_vector_scale)

    def set_max_case(self, cases: Union[list[int], dict[int, Any]]):
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

    def set_case_keys(self, case_keys: list[int]):
        """set the available keys for the case spinner"""
        self.case_keys = case_keys
        if self.icase == -1:
            self.icase = case_keys[0]
        self.set_max_case(case_keys)

    def update_icase(self, icase_frige: int):
        """callback for updating the case spinner"""
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

    def _set_case(self, icase: int) -> None:
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

    def _set_buttons(self, deflection_is_visible: bool,
                     vector_is_visible: bool) -> None:
        """show/hide the additional buttons"""
        self.deflection_label.setVisible(deflection_is_visible)
        self.deflection_edit.setVisible(deflection_is_visible)

        self.vector_label.setVisible(vector_is_visible)
        self.vector_edit.setVisible(vector_is_visible)

    #def on_deflection_scale(self):
        #pass
    #def on_vector_scale(self):
        #pass

    @property
    def result_case_window(self) -> str:
        i = self.names.index(self.name)
        return self.result_case_windows[i]

    def setup_layout(self, data, choices,
                     clear_data: bool=True,
                     init: bool=True) -> None:
        """creates the sidebar visual layout"""
        #if not init:
            #self.frameGeometry().
            #width = self.frameGeometry().width()
            #height = self.frameGeometry().height()
            #print('width=%s height=%s' % (width, height))

        vbox = QVBoxLayout()

        hbox_name = create_hbox_with_widgets([self.name_label, self.name_pulldown])
        if USE_NEW_SIDEBAR:
            hbox_avg = create_hbox_with_widgets([self.min_max_average_label, self.min_max_average_pulldown])
            hbox_coord = create_hbox_with_widgets([self.transform_coords_label, self.transform_coords_pulldown])
            hbox_combine = create_hbox_with_widgets([self.combine_label, self.combine_pulldown])

        #self.case_spinner_label, self.case_spinner,
        #self.deflection_label, self.deflection_edit,
        #self.vector_label, self.vector_edit,
        vbox.addLayout(hbox_name)

        if USE_NEW_SIDEBAR:
            vbox_top = QVBoxLayout(self)
            vbox_btm = QVBoxLayout(self)
        else:
            vbox_top = vbox

        nwindows = len(self.result_case_windows)
        #print('nwindows=%s self.names=%s' % (nwindows, self.names))
        for i in range(nwindows):
            #print('*using existing window')
            result_case_window = self.result_case_windows[i]
            vbox_top.addWidget(result_case_window)
            #result_case_window.setVisible(False)  # be very careful of this...

        nwindows = len(self.result_case_windows)
        for name in self.names[nwindows:]:
            #print('*creating a window')
            result_case_window = ResultsWindow(self, 'Case/Results', data, choices)
            result_case_window.setVisible(False)
            vbox_top.addWidget(result_case_window)
            self.result_case_windows.append(result_case_window)

        iname = 0
        if USE_NEW_SIDEBAR:
            vbox_top_widget = QWidget(self); vbox_top_widget.setLayout(vbox_top)
            vbox_btm_widget = QWidget(self); vbox_btm_widget.setLayout(vbox_btm)

            vbox_btm.addWidget(self.result_method_window)

            # TODO: make the border/split line thinner
            #vbox_top_widget.setMinimumHeight(40)
            #vbox_btm_widget.setMinimumHeight(40)

            splitter = QSplitter(QtCore.Qt.Vertical)
            splitter.addWidget(vbox_top_widget)
            splitter.addWidget(vbox_btm_widget)
            splitter.setCollapsible(0, False)
            splitter.setCollapsible(1, False)
            #splitter.setStyleSheet( "margin-btm: 1px" )
            splitter.setStretchFactor(2, 1)
            QApplication.setStyle(QStyleFactory.create('Cleanlooks'))
            vbox.addWidget(splitter)
            #if self.name is None:
                #iname = 0
            #else:
                #iname = self.names.index(self.name)
            #for i in range(nwindows):
                #if i != iname:
                    #self.result_case_windows[iname].setVisible(False)
            #self.result_case_windows[iname].setVisible(True)

        if self.show_pulldown:
            vbox.addWidget(self.pulldown)
        if USE_NEW_SIDEBAR:
            vbox.addLayout(hbox_coord)
            vbox.addLayout(hbox_avg)
            vbox.addLayout(hbox_combine)
        self._add_grid_to_vbox(vbox)

        vbox.addWidget(self.apply_button)
        self.setLayout(vbox)

        if clear_data:
            self.clear_data()

        #if not init:
            #self.frameGeometry().width()
            #self.frameGeometry().height()
            #self.resize(width, height)

    def _add_grid_to_vbox(self, vbox: QVBoxLayout) -> int:
        """creates the grid for some optional buttons"""
        grid = QGridLayout()
        irow = 0
        irow = add_line_widgets_to_grid(
            grid, irow, [self.case_spinner_label, self.case_spinner])
        irow = add_line_widgets_to_grid(
            grid, irow, [self.deflection_label, self.deflection_edit])
        irow = add_line_widgets_to_grid(
            grid, irow, [self.vector_label, self.vector_edit])
        vbox.addLayout(grid)
        return irow

    def update_method(self, method: list[str]) -> None:
        """called by cycle_results_explicit -> _set_case"""
        if isinstance(method, str):
            method = [method]
        assert isinstance(method, list), method
        assert len(method) > 0, method
        data = []
        for methodi in method:
            #datai = self.result_method_window.data[0]
            data.append((methodi, None, []))
            #self.result_method_window.data[0] = (methodi, datai[1], datai[2])
            #print('method=%s datai=%s' % (method, datai))
        if not USE_NEW_SIDEBAR:
            return
        self.result_method_window.data = data
        self.result_method_window.update_data(self.result_method_window.data)

        #else:
            # pragma: no cover
           #datai = self.result_method_window.data[0]

    def get_form(self):
        """
        TODO: At this point, we should clear out the data block and refresh it
        """
        return self.result_case_window.data

    def update_results(self, data: list[Any], name: str):
        """
        Updates the sidebar

        Parameters
        ----------
        data : list[tuple]
            the form data
        name : str
            the name that goes at the side

        """
        name = str(name)
        update_name = False
        setup_layout = False
        assert name is not None, name
        #if self.name is None:
            #self.names = [name]
            #self.name_pulldown.clear()
            #self.name_pulldown.addItems(self.names)
            ##self.name_pulldown.setItemText(0, name)
            ##self.name_pulldown.setCurrentIndex(1)
            #update_name = True
            #setup_layout = True

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
                ('Top', 1, []),
                ('Bottom', 2, []),
                #('C', 3, []),
            ]
            self.setup_layout(data, choices, init=False, clear_data=True)

    def update_methods(self, data):
        """the methods is a hidden box"""
        if USE_NEW_SIDEBAR:
            self.result_method_window.update_data(data)
            self.result_method_window.tree_view.selectAll()
        assert len(data), data
        self.apply_button.setEnabled(True)

    def clear_data(self):
        self.result_case_window.clear_data()
        if USE_NEW_SIDEBAR:
            self.result_method_window.clear_data()
        self.apply_button.setEnabled(False)

    def on_pulldown(self):  # pragma: no cover
        print('pulldown...')

    def on_update_name(self):  # pragma: no cover
        """user clicked the pulldown"""
        name = str(self.name_pulldown.currentText())
        data = self.parent._get_sidebar_data(name)
        #self.result_case_window.update_data(data)

    def on_click_result_name(self):
        asdf

    def on_apply(self):
        data = self.result_case_window.data
        valid_a, keys_a = self.result_case_window.tree_view.get_row()

        if USE_NEW_SIDEBAR:
            data = self.result_method_window.data
            valid_b, keys_b = self.result_method_window.tree_view.get_row()
            if valid_a and valid_b:
                if self.debug:  # pragma: no cover
                    print('  rows1 = %s' % self.result_case_window.tree_view.old_rows)
                    print('        = %s' % str(keys_a))
                    print('  rows2 = %s' % self.result_method_window.tree_view.old_rows)
                    print('        = %s' % str(keys_b))
                else:
                    self.update_vtk_window(keys_a, keys_b)
        else:
            keys_b = [0]
            self.update_vtk_window(keys_a, keys_b)

    def update_vtk_window(self,
                          keys_a: int,
                          keys_b: int) -> None:
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
        icase = keys_a
        if icase is None:
            #self.parent.log.error(f"icase={icase} and you're trying to set a result...")
            return

        #self.case_spinner.setValue(icase)  # this might just work?

        # set the spinner, but don't take any actions
        if 0:
            if self._update_case:
                self._update_case = False
                self.case_spinner.setValue(icase)
                self._update_case = True
            else:
                self.case_spinner.setValue(icase)
        result_name = None
        #self._set_case(i)
        self.parent._set_case(result_name, icase, explicit=True)

    @property
    def has_cases(self) -> bool:
        return hasattr(self.parent, 'result_cases')

    def on_clear_results(self) -> None:
        if hasattr(self.parent, 'on_clear_results'):
            self.parent.on_clear_results()
    def on_fringe(self, icase: int) -> None:
        if hasattr(self.parent, 'on_fringe'):
            self.parent.on_fringe(icase)
    def on_disp(self, icase: int) -> None:
        if hasattr(self.parent, 'on_disp'):
            self.parent.on_disp(icase)
    def on_vector(self, icase: int) -> None:
        if hasattr(self.parent, 'on_vector'):
            self.parent.on_vector(icase)

    def get_clicked(self) -> None:
        self.result_method_window.tree_view.get_clicked()


def main():  # pragma: no cover
    app = QApplication(sys.argv)

    form = [
        ['Geometry', None, [
            ('NodeID', 0, []),
            ('ElementID', 1, []),
            #('PropertyID', 2, []),
            #('MaterialID', 3, []),
            #('E', 4, []),
            ('Element Checks', None, [
                ('ElementDim', 2, []),
                ('Min Edge Length', 3, []),
                #('Min Interior Angle', 7, []),
                #('Max Interior Angle', 8, [])],
                ],
            ),],
        ],
        [
            'Results', None, [
                ('Centroidal Stress', 4, []),
                ('Nodal Stress', 5, []),
                ],
        ],
    ]
    #form = []
    def func(case: int):
        print(f'case={case}')
        x = 1

    res_widget = Sidebar(
        app, data=form,
        #left_click_callback=func,
        clear_data=False, debug=True)
    res_widget.update_results(data=form, name='main')
    res_widget.result_case_window.tree_view.expandAll()
    res_widget.case_keys = [1, 2, 3, 4, 5]

    name = 'name'
    #res_widget.update_results(form, name)

    res_widget.show()
    sys.exit(app.exec_())

if __name__ == "__main__":  # pragma: no cover
    main()
