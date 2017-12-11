from __future__ import print_function
import sys
from copy import deepcopy

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from six import string_types

from qtpy import QtGui
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QTreeView, QWidget, QAbstractItemView, QVBoxLayout, QPushButton, QApplication,
    QComboBox, QLabel, QHBoxLayout, QMessageBox)


class QTreeView2(QTreeView):
    def __init__(self, parent, data, choices):
        self.parent = parent
        self.old_rows = []
        self.data = data
        self.cases_deleted = set([])
        self.choices = choices
        self.single = False
        QTreeView.__init__(self)
        #self.setAlternatingRowColors(True)

    def keyPressEvent(self, event): #Reimplement the event here, in your case, do nothing
        #if event.key() == QtCore.Qt.Key_Escape:
            #self.close()
        #return
        key = event.key()
        if key == Qt.Key_Delete:
            rows = self.find_list_index()
            self.remove_rows(rows)
            #self.model().removeRow(row)
            #self.parent().on_delete(index.row())

            #del self.data[row]
            #self.model().change_data(self.data)
        else:
            QTreeView.keyPressEvent(self, event)

    def remove_rows(self, rows):
        """
        We just hide the row to delete things to prevent refreshing
        the window and changing which items have been expanded

        Parameters
        ----------
        rows : List[int]
            the trace on the data/form block

        form = [
            ['Geometry', None, [
                ('NodeID', 0, []),
                ('ElementID', 1, []),
                ('PropertyID', 2, []),
                ('MaterialID', 3, []),
                ('E', 4, []),
                ('Element Checks', None, [
                    ('ElementDim', 5, []),
                    ('Min Edge Length', 6, []),
                    ('Min Interior Angle', 7, []),
                    ('Max Interior Angle', 8, [])],
                 ),],
             ],
        ]

        # delete Geometry
        data[0] = ('Geometry', None, [...])
        >>> remove_rows([0])

        # delete MaterialID
        data[0][3] = ('MaterialID', 3, [])
        >>> remove_rows([0, 3])

        # delete ElementChecks
        data[0][5] = ('Element Checks', None, [...])
        >>> remove_rows([0, 5])

        # delete Min Edge Length
        data[0][5][1] = ('Min Edge Length', 6, [])
        >>> remove_rows([0, 5, 1])
        """
        # find the row the user wants to delete
        data = self.data
        for row in rows[:-1]:
            data = data[row][2]

        # we got our data block
        # now we need to get 1+ results
        last_row = rows[-1]
        cases_to_delete = get_many_cases(data[last_row])

        cases_to_delete = list(set(cases_to_delete) - self.cases_deleted)
        cases_to_delete.sort()
        if len(cases_to_delete) == 0: # can this happen?
            # this happens when you cleared out a data block by
            # deleting to entries, but not the parent
            #
            # we'll just hide the row now
            msg = ''
            return
        elif len(cases_to_delete) == 1:
            msg = 'Are you sure you want to delete 1 result case load_case=%s' % cases_to_delete[0]
        else:
            msg = 'Are you sure you want to delete %s result cases %s' % (
                len(cases_to_delete), str(cases_to_delete))

        if msg:
            widget = QMessageBox()
            title = 'Delete Cases'
            result = QMessageBox.question(widget, title, msg,
                                          QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

            if result != QMessageBox.Yes:
                return

            self.cases_deleted.update(set(cases_to_delete))
            sidebar = self.parent.parent
            if hasattr(sidebar.parent, 'delete_cases'):
                gui = sidebar.parent
                gui.delete_cases(cases_to_delete, ask=False)

        # hide the line the user wants to delete
        row = rows[-1]
        indexes = self.selectedIndexes()
        self.setRowHidden(row, indexes[-1].parent(), True)
        self.update()

    def find_list_index(self):
        """
        trace the tree to find the selected item
        """
        indexes = self.selectedIndexes()

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
        return rows

    def mousePressEvent(self, position):
        """called when you click a result"""
        QTreeView.mousePressEvent(self, position)

        # trace the tree to find the selected item
        rows = self.find_list_index()

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
                1 - iCase
                2 - []
        """
        # if there's only 1 data member, we don't need to extract the data id
        if self.single:
            return True, self.data[0]

        irow = 0

        # TODO: what is this for???
        #     crashes some PyQt4 cases when clicking on the first
        #     non-results level of the sidebar
        #data = deepcopy(self.data)
        data = self.data
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
    def __init__(self, parent, debug=False, data=None, clear_data=True):
        """creates the buttons in the Sidebar, not the actual layout"""
        QWidget.__init__(self)
        self.parent = parent
        self.debug = debug

        name = 'main'

        choices = ['keys2', 'purse2', 'cellphone2', 'credit_card2', 'money2']
        if data is None:
            data = []

        self.result_case_window = ResultsWindow(self, 'Case/Results', data, choices)

        data = [
            ('A', 1, []),
            #('B', 2, []),
            #('C', 3, []),
        ]
        self.result_data_window = ResultsWindow(self, 'Method', data, choices)
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

        self.setup_layout(clear_data=clear_data)

    def setup_layout(self, clear_data=True):
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

        if clear_data:
            self.clear_data()

    def update_method(self, method):
        if isinstance(method, string_types):
            datai = self.result_data_window.data[0]
            self.result_data_window.data[0] = (method, datai[1], datai[2])
            #print('method=%s datai=%s' % (method, datai))
            self.result_data_window.update_data(self.result_data_window.data)
        else:
            return
             # pragma: no cover
            #datai = self.result_data_window.data[0]

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
    def __init__(self, parent, name, data, choices):
        QWidget.__init__(self)
        self.name = name
        self.data = data
        self.choices = choices
        self.parent = parent
        self.treeView = QTreeView2(self, self.data, choices)
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
            raise RuntimeError('cannot add data=\n%s' % data)
            #if isinstance(data, string_types):
                #self.addItems(self.model, data)
            #else:
                #self.addItems(self.model, *tuple(data))
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

def get_many_cases(data):
    """
    Get the result case ids that are a subset of the data/form list

    data = [
        (u'Element Checks', None, [
            (u'ElementDim', 5, []),
            (u'Min Edge Length', 6, []),
            (u'Min Interior Angle', 7, []),
            (u'Max Interior Angle', 8, [])],
        ),
    ]
    >>> get_many_cases(data)
    [5, 6, 7, 8]

    >>> data = [(u'Max Interior Angle', 8, [])]
    [8]
    """
    name, case, rows = data
    if case is None:
        # remove many results
        # (Geometry, None, [results...])
        cases = []
        for irow, row in enumerate(rows):
            name, row_id, data2 = row
            cases += get_many_cases(row)
    else:
        cases = [case]
    return cases

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

