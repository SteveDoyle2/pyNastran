from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QPushButton,
    QComboBox, QLabel, QHBoxLayout, QBoxLayout, )

from pyNastran.gui.utils.qt.results_window import ResultsWindow


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
    def __init__(self, parent, debug=False, data=None, actions=None,
                 results_window_title='Case/Results',
                 clear_data=True, name='main', setup_dict=None):
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
        results_window_title : str
            the title for "Case/Results"
        actions : Dict[int] = QWidget, List[QWidget]
            additional actions
        clear_data : bool; default=True
            ???
        name : str; default='main'
            the active name
        setup_dict : Dict[irow] = List[QWidgets]
            a way to add additional widgets to the sidebar
        """
        QWidget.__init__(self)
        self.parent = parent
        self.debug = debug
        self.setup_dict = setup_dict
        self.results_window_title = results_window_title

        choices = ['keys2', 'purse2', 'cellphone2', 'credit_card2', 'money2']
        if data is None:
            data = []

        if actions is None:
            actions = [
                ('Clear Results...', self.on_clear_results, False),
                ('Apply Results to Fringe...', self.on_fringe, True),
                ('Apply Results to Displacement...', self.on_disp, True),
                ('Apply Results to Vector...', self.on_vector, True),
            ]

        self.result_case_windows = [
            ResultsWindow(self, self.results_window_title, data, choices, actions=actions)
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

        self.name_label.setVisible(False)
        self.name_pulldown.setVisible(False)

        self.setup_layout(data, choices, clear_data=clear_data)
        self.name_pulldown.currentIndexChanged.connect(self.on_update_name)

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
        if not init:
            #self.frameGeometry().
            unused_width = self.frameGeometry().width()
            unused_height = self.frameGeometry().height()
            #print('width=%s height=%s' % (width, height))

        #print('init...')
        vbox = QVBoxLayout()
        hbox = QHBoxLayout()

        irow = 0
        self._add_from_setup_dict(vbox, irow)

        hbox.addWidget(self.name_label)
        hbox.addWidget(self.name_pulldown)
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
        for unused_name in self.names[nwindows:]:
            #print('*creating a window')
            result_case_window = ResultsWindow(
                self, self.results_window_title, data, choices)
            result_case_window.setVisible(False)
            vbox.addWidget(result_case_window)
            self.result_case_windows.append(result_case_window)

        unused_iname = 0
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

    def _add_from_setup_dict(self, vbox, irow):
        """
        0
        Name
        1
        Case/Results Window
        3
        Apply
        4
        """
        if self.setup_dict is None:
            return

        if irow in self.setup_dict:
            widgets = self.setup_dict[irow]
            assert widgets is not None, widgets
            if isinstance(widgets, list):
                for widget_layout in widgets:
                    if isinstance(widget_layout, QBoxLayout):
                        vbox.addLayout(widget_layout)
                    else:
                        vbox.addWidget(widget_layout)
            else:
                widget_layout = widgets
                if isinstance(widget_layout, QBoxLayout):
                    vbox.addLayout(widget_layout)
                else:
                    vbox.addWidget(widget_layout)

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

        if setup_layout and 0:
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
        unused_data = self.parent._get_sidebar_data(name)
        #self.result_case_window.update_data(data)

    def on_apply(self):
        unused_data = self.result_case_window.data
        valid_a, keys_a = self.result_case_window.treeView.get_row()

        unused_data = self.result_method_window.data
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
            result_name = key[1]
            #self.parent.cycle_results_explicit(result_name=result_name, explicit=True)
            #j = self.parent._get_icase(result_name)
            #j = i
        i = keys_a
        result_name = None
        self.parent._set_case(result_name, i, explicit=True)

    def on_clear_results(self):
        #print('*clear')
        if hasattr(self.parent, 'on_clear_results'):
            self.parent.on_clear_results()
    def on_fringe(self, icase):
        #print('*fringe', icase)
        if hasattr(self.parent, 'on_fringe'):
            self.parent.on_fringe(icase)
    def on_disp(self, icase):
        #print('*disp', icase)
        if hasattr(self.parent, 'on_disp'):
            self.parent.on_disp(icase)
    def on_vector(self, icase):
        #print('*vector', icase)
        if hasattr(self.parent, 'on_vector'):
            self.parent.on_vector(icase)

    def get_clicked(self):
        self.result_method_window.treeView.get_clicked()
