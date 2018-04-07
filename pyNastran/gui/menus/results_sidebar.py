from __future__ import print_function
import sys

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from six import string_types

#from qtpy import QtGui
#from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QPushButton, QApplication,
    QComboBox, QLabel, QHBoxLayout)
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
    def __init__(self, parent, debug=False, data=None, clear_data=True, name='main'):
        """creates the buttons in the Sidebar, not the actual layout"""
        QWidget.__init__(self)
        self.parent = parent
        self.debug = debug

        choices = ['keys2', 'purse2', 'cellphone2', 'credit_card2', 'money2']
        if data is None:
            data = []

        #self.result_case_windows = []
        self.result_case_windows = [ResultsWindow(self, 'Case/Results', data, choices)]
        data = [
            ('A', 1, []),
            #('B', 2, []),
            #('C', 3, []),
        ]
        self.result_method_window = ResultsWindow(self, 'Method', data, choices)
        self.result_method_window.setVisible(False)

        self.show_pulldown = False
        if self.show_pulldown:
            combo_options = ['a1', 'a2', 'a3']
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

        self.setup_layout(data, choices, clear_data=clear_data)
        self.name_pulldown.currentIndexChanged.connect(self.on_update_name)

    @property
    def result_case_window(self):
        #if self.name is None:
            #i = 0
        #else:
            #i = self.names.index(self.name)
        i = 0
        return self.result_case_windows[i]

    def setup_layout(self, data, choices, clear_data=True):
        """creates the sidebar visual layout"""
        vbox = QVBoxLayout()
        hbox = QHBoxLayout()

        hbox.addWidget(self.name_label)
        hbox.addWidget(self.name_pulldown)
        vbox.addLayout(hbox)

        nwindows = len(self.result_case_windows)
        print('nwindows=%s self.names=%s' % (nwindows, self.names))
        for i in range(nwindows):
            print('*using existing window')
            result_case_window = self.result_case_windows[i]
            vbox.addWidget(result_case_window)
            result_case_window.setVisible(False)

        for name in self.names[nwindows:]:
            print('*creating a window')
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
        self.result_case_windows[iname].setVisible(True)

        vbox.addWidget(self.result_method_window)
        if self.show_pulldown:
            vbox.addWidget(self.pulldown)
        vbox.addWidget(self.apply_button)
        self.setLayout(vbox)

        if clear_data:
            self.clear_data()

    def update_method(self, method):
        if isinstance(method, string_types):
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
            self.setup_layout(data, choices, clear_data=True)

    def update_methods(self, data):
        """the methods is a hidden box"""
        self.result_method_window.update_data(data)
        self.apply_button.setEnabled(True)

    def clear_data(self):
        self.result_case_window.clear_data()
        self.result_method_window.clear_data()
        self.apply_button.setEnabled(False)

    def on_pulldown(self, event):  # pragma: no cover
        print('pulldown...')

    def on_update_name(self, event):  # pragma: no cover
        """user clicked the pulldown"""
        name = str(self.name_pulldown.currentText())
        data = self.parent._get_sidebar_data(name)
        #self.result_case_window.update_data(data)

    def on_apply(self, event):
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
            result_name = key[1]
            #self.parent.cycle_results_explicit(result_name=result_name, explicit=True)
            #j = self.parent._get_icase(result_name)
            #j = i
        i = keys_a
        result_name = None
        self.parent._set_case(result_name, i, explicit=True)

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
