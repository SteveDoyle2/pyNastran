from __future__ import print_function
import sys
from copy import deepcopy

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from six import string_types

from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4 import QtGui
    from PyQt4.QtGui import (
        QTreeView, QWidget, QAbstractItemView, QVBoxLayout, QPushButton, QApplication,
        QComboBox, QLabel, QHBoxLayout)
elif qt_version == 5:
    from PyQt5 import QtGui
    from PyQt5.QtWidgets import (
        QTreeView, QWidget, QAbstractItemView, QVBoxLayout, QPushButton, QApplication,
        QComboBox, QLabel, QHBoxLayout)
elif qt_version == 'pyside':
    from PySide import QtGui
    from PySide.QtGui import (
        QTreeView, QWidget, QAbstractItemView, QVBoxLayout, QPushButton, QApplication,
        QComboBox, QLabel, QHBoxLayout)
else:
    raise NotImplementedError('qt_version = %r' % qt_version)


class QTreeView2(QTreeView):
    def __init__(self, data, choices):
        self.old_rows = []
        self.data = data
        self.choices = choices
        self.single = False
        QTreeView.__init__(self)

    def mousePressEvent(self, position):
        QTreeView.mousePressEvent(self, position)
        indexes = self.selectedIndexes()

        # trace the tree to find the selected item
        rows = []
        for index in indexes:
            row = index.row()
            rows.append(row)
            level = 0
            while index.parent().isValid():
                index = index.parent()
                row = index.row()
                rows.append(row)
                level += 1
        rows.reverse()

        # TODO: what is this for???
        if rows != self.old_rows:
            self.old_rows = rows
        valid, keys = self.get_row()
        if not valid:
            print('invalid=%s keys=%s' % (valid, keys))
        else:
            print('valid=%s keys=%s' % (valid, keys))
            #print('choice =', self.choices[keys])

    def get_row(self):
        """
        gets the row

        Returns
        -------
        is_valid : bool
            is this case valid
        row : None or tuple
            None : invalid case
            tuple : valid case
                ('centroid', None, [])
                0 - the location (e.g. node, centroid)
                1 - ???
                2 - ???
        """
        # if there's only 1 data member, we don't need to extract the data id
        if self.single:
            return True, self.data[0]

        # TODO: what is this for???
        irow = 0
        data = deepcopy(self.data)
        for row in self.old_rows:
            try:
                key = data[row][0]
            except IndexError:
                return False, irow
            irow = data[row][1]
            data = data[row][2]

        if data:
            return False, None
        return True, irow

    def set_single(self, single):
        self.single = single
        self.old_rows = [0]


#for subcase in subcases:
#    for time in times:
#        disp
#        stress
#        load

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
    def __init__(self, parent, debug=False):
        """creates the buttons in the Sidebar, not the actual layout"""
        QWidget.__init__(self)
        self.parent = parent
        self.debug = debug

        name = 'main'
        data = []
        data = [
            ("Alice", None, [
                ("Keys", 1, []),
                ("Purse", 2, [
                    ("Cellphone", 3, [])
                    ])
                ]),
            ("Bob", None, [
                ("Wallet", None, [
                    ("Credit card", 4, []),
                    ("Money", 5, [])
                    ])
                ]),
            ]

        choices = ['keys2', 'purse2', 'cellphone2', 'credit_card2', 'money2']
        self.result_case_window = ResultsWindow('Case/Results', data, choices)

        data = [
            ('A', 1, []),
            #('B', 2, []),
            #('C', 3, []),
        ]
        self.result_data_window = ResultsWindow('Method', data, choices)
        self.result_data_window.setVisible(False)

        self.show_pulldown = False
        if self.show_pulldown:
            combo_options = ['a1', 'a2', 'a3']
            self.pulldown = QComboBox()
            self.pulldown.addItems(choices)
            self.pulldown.activated[str].connect(self.on_pulldown)

        self.apply_button = QPushButton('Apply', self)
        self.apply_button.clicked.connect(self.on_apply)

        self.name = str(name)
        self.names = [name]
        self.name_label = QLabel("Name:")
        self.name_pulldown = QComboBox()
        self.name_pulldown.addItem(name)
        self.name_pulldown.setDisabled(True)
        self.name_pulldown.currentIndexChanged.connect(self.on_update_name)

        self.setup_layout()

    def setup_layout(self):
        """creates the sidebar visual layout"""
        vbox = QVBoxLayout()
        hbox = QHBoxLayout()

        hbox.addWidget(self.name_label)
        hbox.addWidget(self.name_pulldown)
        vbox.addLayout(hbox)
        vbox.addWidget(self.result_case_window)
        vbox.addWidget(self.result_data_window)
        if self.show_pulldown:
            vbox.addWidget(self.pulldown)
        vbox.addWidget(self.apply_button)
        self.setLayout(vbox)

        self.clear_data()

    def update_method(self, method):
        if isinstance(method, string_types):
            datai = self.result_data_window.data[0]
            self.result_data_window.data[0] = (method, datai[1], datai[2])
            print('method=%s datai=%s' % (method, datai))
            self.result_data_window.update_data(self.result_data_window.data)
        else:
            return
             # pragma: no cover
            datai = self.result_data_window.data[0]

    def get_form(self):
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
        if name in self.names:
            i = self.names.index(name)
            self.name_pulldown.setCurrentIndex(i)
        else:
            self.name_pulldown.addItem(name)
            self.names.append(name)
        if len(self.names) >= 2:
            self.name_pulldown.setEnabled(True)
        self.name = name

        self.result_case_window.update_data(data)
        self.apply_button.setEnabled(True)

    def update_methods(self, data):
        """the methods is a hidden box"""
        self.result_data_window.update_data(data)
        self.apply_button.setEnabled(True)

    def clear_data(self):
        self.result_case_window.clear_data()
        self.result_data_window.clear_data()
        self.apply_button.setEnabled(False)

    def on_pulldown(self, event):
        print('pulldown...')

    def on_update_name(self, event):
        """user clicked the pulldown"""
        name = str(self.name_pulldown.currentText())
        data = self.parent._get_sidebar_data(name)
        #self.result_case_window.update_data(data)

    def on_apply(self, event):
        data = self.result_case_window.data
        valid_a, keys_a = self.result_case_window.treeView.get_row()

        data = self.result_data_window.data
        valid_b, keys_b = self.result_data_window.treeView.get_row()
        if valid_a and valid_b:
            if self.debug:  # pragma: no cover
                print('  rows1 = %s' % self.result_case_window.treeView.old_rows)
                print('        = %s' % str(keys_a))
                print('  rows2 = %s' % self.result_data_window.treeView.old_rows)
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


class ResultsWindow(QWidget):
    """
    A ResultsWindow creates the box where we actually select our
    results case.  It does not have an apply button.
    """
    def __init__(self, name, data, choices):
        QWidget.__init__(self)
        self.name = name
        self.data = data
        self.choices = choices

        self.treeView = QTreeView2(self.data, choices)
        self.treeView.setEditTriggers(QAbstractItemView.NoEditTriggers)

        self.model = QtGui.QStandardItemModel()
        is_single = self.addItems(self.model, data)
        self.treeView.setModel(self.model)
        self.treeView.set_single(is_single)

        self.model.setHorizontalHeaderLabels([self.tr(self.name)])

        layout = QVBoxLayout()
        layout.addWidget(self.treeView)
        self.setLayout(layout)

    def update_data(self, data):
        self.clear_data()
        self.data = data
        try:
            self.addItems(self.model, data)
        except:
            adf
            if isinstance(data, string_types):
                self.addItems(self.model, data)
            else:
                self.addItems(self.model, *tuple(data))
        self.treeView.data = data
        #layout = QVBoxLayout()
        #layout.addWidget(self.treeView)
        #self.setLayout(layout)

    def clear_data(self):
        self.model.clear()
        self.treeView.data = []
        self.model.setHorizontalHeaderLabels([self.tr(self.name)])

    def addItems(self, parent, elements, level=0, count_check=False):
        nelements = len(elements)
        redo = False
        #print(elements[0])
        try:
            #if len(elements):
                #assert len(elements[0]) == 3, 'len=%s elements[0]=%s\nelements=\n%s\n' % (
                    #len(elements[0]), elements[0], elements)
            for element in elements:
                #if isinstance(element, str):
                    #print('elements = %r' % str(elements))

                #print('element = %r' % str(element))
                if not len(element) == 3:
                    print('element = %r' % str(element))
                text, i, children = element
                nchildren = len(children)
                #print('text=%r' % text)
                item = QtGui.QStandardItem(text)
                parent.appendRow(item)

                # TODO: count_check and ???
                if nelements == 1 and nchildren == 0 and level == 0:
                    #self.result_data_window.setEnabled(False)
                    item.setEnabled(False)
                    #print(dir(self.treeView))
                    #self.treeView.setCurrentItem(self, 0)
                    #item.mousePressEvent(None)
                    redo = True
                #else:
                    #pass
                    #print('item=%s count_check=%s nelements=%s nchildren=%s' % (
                        #text, count_check, nelements, nchildren))
                if children:
                    assert isinstance(children, list), children
                    self.addItems(item, children, level + 1, count_check=count_check)
            is_single = redo
            return is_single
        except ValueError:
            print()
            print('elements =', elements)
            print('element =', element)
            print('len(elements)=%s' % len(elements))
            for e in elements:
                print('  e = %s' % str(e))
            raise
        #if redo:
        #    data = [
        #        ('A', []),
        #        ('B', []),
        #    ]
        #    self.update_data(data)

def main():  # pragma: no cover
    app = QApplication(sys.argv)
    window = Sidebar(app, debug=True)
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":  # pragma: no cover
    main()

