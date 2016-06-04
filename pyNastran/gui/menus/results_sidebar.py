from __future__ import print_function
import sys
from copy import deepcopy
from PyQt4 import QtGui
#from PyQt4.QtCore import *
#from PyQt4.QtGui import *

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


class QTreeView2(QtGui.QTreeView):
    def __init__(self, data):
        self.old_rows = []
        self.data = data
        self.single = False
        QtGui.QTreeView.__init__(self)

    def mousePressEvent(self, position):
        QtGui.QTreeView.mousePressEvent(self, position)
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

    def get_row(self):
        # if there's only 1 data member, we don't need to extrac the data id
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

class Sidebar(QtGui.QWidget):
    """
    +--------------+
    | Case/Results |
    +==============+
    | - a          |
    |  - b1        |
    |  - b2        |
    |  - b3        |
    +--------------+

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
        QtGui.QWidget.__init__(self)
        self.parent = parent
        self.debug = debug

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

        self.result_case_window = ResultsWindow('Case/Results', data)

        data = [
            ('A', 1, []),
            #('B', 2, []),
            #('C', 3, []),
        ]
        self.result_data_window = ResultsWindow('Method', data)
        self.result_data_window.setVisible(False)

        if 0:
            combo_options = ['a1', 'a2', 'a3']
            self.pulldown = QtGui.QComboBox()
            self.pulldown.addItems(combo_options)
            self.pulldown.activated[str].connect(self.on_pulldown)

        self.apply_button = QtGui.QPushButton('Apply', self)
        self.apply_button.clicked.connect(self.on_apply)
        self.setup_layout()

    def setup_layout(self):
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.result_case_window)
        layout.addWidget(self.result_data_window)
        #layout.addWidget(self.pulldown)
        layout.addWidget(self.apply_button)
        self.setLayout(layout)

        self.clear_data()

    def update_method(self, method):
        datai = self.result_data_window.data[0]
        self.result_data_window.data[0] = (method, datai[1], datai[2])
        self.result_data_window.update_data(self.result_data_window.data)

    def get_form(self):
        return self.result_case_window.data

    def update_results(self, data):
        self.result_case_window.update_data(data)
        self.apply_button.setEnabled(True)

    def update_methods(self, data):
        self.result_data_window.update_data(data)
        self.apply_button.setEnabled(True)

    def clear_data(self):
        self.result_case_window.clear_data()
        self.result_data_window.clear_data()
        self.apply_button.setEnabled(False)

    def on_pulldown(self, event):
        print('pulldown...')

    def on_apply(self, event):
        data = self.result_case_window.data
        valid_a, keys_a = self.result_case_window.treeView.get_row()

        data = self.result_data_window.data
        valid_b, keys_b = self.result_data_window.treeView.get_row()
        if valid_a and valid_b:
            if self.debug:
                print('  rows1 = %s' % self.result_case_window.treeView.old_rows)
                print('        = %s' % str(keys_a))
                print('  rows2 = %s' % self.result_data_window.treeView.old_rows)
                print('        = %s' % str(keys_b))
            else:
                self.update_vtk_window(keys_a, keys_b)

    def update_vtk_window(self, keys_a, keys_b):
        if 0:
            #print('keys_a = %s' % str(keys_a))
            for i, key in enumerate(self.parent.case_keys):
                if key[1] == keys_a[0]:
                    break
            #print('*i=%s key=%s' % (i, str(key)))
            #self.parent.update_vtk_window_by_key(i)
            result_name = key[1]
            #self.parent.cycle_results_explicit(result_name=result_name, explicit=True)
            #j = self.parent._get_icase(result_name)
            #j = i
        i = keys_a
        result_name = None
        self.parent._set_case(result_name, i, explicit=True)


class ResultsWindow(QtGui.QWidget):
    def __init__(self, name, data):
        QtGui.QWidget.__init__(self)
        self.name = name
        self.data = data

        self.treeView = QTreeView2(self.data)
        self.treeView.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)

        self.model = QtGui.QStandardItemModel()
        is_single = self.addItems(self.model, data)
        self.treeView.setModel(self.model)
        self.treeView.set_single(is_single)

        self.model.setHorizontalHeaderLabels([self.tr(self.name)])

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.treeView)
        self.setLayout(layout)

    def update_data(self, data):
        self.clear_data()
        self.data = data
        self.addItems(self.model, data)
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
                else:
                    pass
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

def main():
    app = QtGui.QApplication(sys.argv)
    window = Sidebar(app, debug=True)
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()

