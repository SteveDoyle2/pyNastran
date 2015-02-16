import sys
from copy import deepcopy
from PyQt4 import QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


class QTreeView2(QTreeView):
    def __init__(self, data):
        self.old_rows = []
        self.data = data
        self.single = False
        QTreeView.__init__(self)

    def mousePressEvent(self, position):
    #def openMenu(self, position):
        QTreeView.mousePressEvent(self, position)
        #print('position = %s' % position)
        indexes = self.selectedIndexes()
        #print('indexes', indexes)
        rows = []
        for index in indexes:
            row = index.row()
            rows.append(row)
            #print('  %s' % row)
            level = 0
            while index.parent().isValid():
                index = index.parent()
                row = index.row()
                rows.append(row)
                level += 1
        rows.reverse()
        #print('rows = %s' % rows)

        if rows != self.old_rows:
            #print('rows = %s' % rows)
            self.old_rows = rows
        valid, keys = self.get_row()

    def get_row(self):
        #print('data=%s single=%s' % (self.data, self.single))
        if self.single:
            return True, self.data[0]

        keys = []
        data = deepcopy(self.data)
        for row in self.old_rows:
            key = data[row][0]
            #print('  %r' % key)
            data = data[row][1]
            keys.append(key)
        if data:
            return False, None
        return True, keys

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
    def __init__(self, parent):
        QWidget.__init__(self)
        self.parent = parent

        data = []
        data = [
            ("Alice", [
                ("Keys", []),
                ("Purse", [
                    ("Cellphone", [])
                    ])
                ]),
            ("Bob", [
                ("Wallet", [
                    ("Credit card", []),
                    ("Money", [])
                    ])
                ]),
            ]

        self.result_case_window = ResultsWindow('Case/Results', data)

        data = [
            ('A', []),
            #('B', []),
            #('C', []),
        ]
        self.result_data_window = ResultsWindow('Method', data)

        combo_options = ['a1', 'a2', 'a3']
        self.pulldown = QtGui.QComboBox()
        self.pulldown.addItems(combo_options)
        self.pulldown.activated[str].connect(self.on_pulldown)

        self.apply_button = QtGui.QPushButton('Apply', self)
        self.apply_button.clicked.connect(self.on_apply)

        layout = QVBoxLayout()
        layout.addWidget(self.result_case_window)
        layout.addWidget(self.result_data_window)
        #layout.addWidget(self.pulldown)
        layout.addWidget(self.apply_button)
        self.setLayout(layout)

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
        #print('Apply!')
        data = self.result_case_window.data
        validA, keysA = self.result_case_window.treeView.get_row()

        data = self.result_data_window.data
        validB, keysB = self.result_data_window.treeView.get_row()
        if validA and validB:
            #print('  rows1 = %s' % self.result_case_window.treeView.old_rows)
            #print('        = %s' % str(keysA))
            #print('  rows2 = %s' % self.result_data_window.treeView.old_rows)
            #print('        = %s' % str(keysB))
            self.update_vtk_window(keysA, keysB)

    def update_vtk_window(self, keysA, keysB):
        for i, key in enumerate(self.parent.caseKeys):
            if key[1] == keysA[0]:
                break
        #print('i=%s key=%s' % (i, key[1]))
        #self.parent.update_vtk_window_by_key(i)
        result_name = key[1]
        #self.parent.cycleResults_explicit(result_name=result_name, explicit=True)
        j = self.parent._get_icase(result_name)
        self.parent._set_case(result_name, j, explicit=True)

class ResultsWindow(QWidget):

    def __init__(self, name, data):
        QWidget.__init__(self)
        self.data = data
        self.name = name

        self.treeView = QTreeView2(self.data)
        self.treeView.setEditTriggers(QAbstractItemView.NoEditTriggers)

        self.model = QStandardItemModel()
        is_single = self.addItems(self.model, data)
        self.treeView.setModel(self.model)
        self.treeView.set_single(is_single)

        self.model.setHorizontalHeaderLabels([self.tr(self.name)])

        layout = QVBoxLayout()
        layout.addWidget(self.treeView)
        self.setLayout(layout)

    def update_data(self, data):
        self.clear_data()
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
        for text, children in elements:
            #print('text=%s' % text)
            item = QStandardItem(text)
            parent.appendRow(item)

            # count_check and
            if nelements == 1 and len(children)==0 and level==0:
                #self.result_data_window.setEnabled(False)
                item.setEnabled(False)
                #print(dir(self.treeView))
                #self.treeView.setCurrentItem(self, 0)
                #item.mousePressEvent(None)
                redo = True
            else:
                pass
                #print('item=%s count_check=%s nelements=%s len(children)=%s' % (
                    #text, count_check, nelements, len(children)))
            if children:
                #print('  children=%s' % children)
                self.addItems(item, children, level + 1, count_check=count_check)
        is_single = redo
        return is_single
        #if redo:
        #    data = [
        #        ('A', []),
        #        ('B', []),
        #    ]
        #    self.update_data(data)

def main():
    app = QApplication(sys.argv)
    window = Sidebar(app)
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()

