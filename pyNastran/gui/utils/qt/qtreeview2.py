"""
creates:
 - QTreeView2 - allows for semi-easy access to the QTreeView
 - RightClickTreeView - adds some right click support for results
 - GenericRightClickTreeView - used for more general right click support

"""
from __future__ import annotations
from typing import Callable, TYPE_CHECKING
from functools import partial
from qtpy import QtGui
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QTreeView, QMessageBox, QMenu
ClickCallback = Callable[[int], bool]

if TYPE_CHECKING:
    from pyNastran.gui.menus.results_sidebar import Sidebar
    from pyNastran.gui.main_window import MainWindow


class QTreeView2(QTreeView):
    """
    creates a QTreeView with:
     - a nice-ish way to extract the location in the tree
     - key press delete support

    """
    def __init__(self, parent, data, choices):
        self.parent = parent
        self.old_rows = []
        self.data = data
        self.cases_deleted = set()
        self.choices = choices
        self.single = False
        QTreeView.__init__(self)
        #self.setAlternatingRowColors(True)

    def keyPressEvent(self, event):
        """
        Handles:
         - delete: delete result cases
         - enter/return: apply result
         - up/down/left/right: navigate the tree

        """
        #if event.key() == QtCore.Qt.Key_Escape:
            #self.close()
        #return
        key = event.key()
        if key == Qt.Key_Delete:
            self.on_delete()
        elif key in [Qt.Key_Enter, Qt.Key_Return]:
            return self.parent.parent.on_apply()
        elif key in [Qt.Key_Up, Qt.Key_Down]:
            QTreeView.keyPressEvent(self, event)
            self.set_rows()
            self.on_left_mouse_button()
        else:
            QTreeView.keyPressEvent(self, event)
        return None

    def remove_rows(self, rows):
        """
        We hide the row to delete things to prevent refreshing
        the window and changing which items have been expanded

        Parameters
        ----------
        rows : list[int]
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
            self.on_delete_parent_cases(cases_to_delete)

        # hide the line the user wants to delete
        row = rows[-1]
        indexes = self.selectedIndexes()
        self.setRowHidden(row, indexes[-1].parent(), True)
        self.update()

    def on_delete_parent_cases(self, cases_to_delete: list[int]) -> None:
        """delete cases from the parent"""
        sidebar = self.parent.parent
        if hasattr(sidebar.parent, 'delete_cases'):
            gui = sidebar.parent
            gui.delete_cases(cases_to_delete, ask=False)

    def find_list_index(self):
        """trace the tree to find the selected item"""
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

    def mousePressEvent(self, event):
        """called when you click a result"""
        QTreeView.mousePressEvent(self, event)
        if event.button() == Qt.LeftButton:
            self.on_left_mouse_button()
        elif event.button() == Qt.RightButton:
            self.on_right_mouse_button()

    def on_left_mouse_button(self):
        """overwrite this to make a right-click menu"""
        pass
        #self.set_rows()
        #unused_valid, unused_keys, imethods = self.get_row()
        #x = 1
        #if not valid:
            #print('invalid=%s keys=%s' % (valid, keys))
        #else:
            #print('valid=%s keys=%s' % (valid, keys))
            #print('choice =', self.choices[keys])

    def on_left_mouse_button(self):
        pass

    def on_right_mouse_button(self):
        """overwrite this to make a right-click menu"""
        pass

    def on_delete(self):
        """interface for deleting rows"""
        rows = self.find_list_index()
        self.remove_rows(rows)
        #self.model().removeRow(row)
        #self.parent().on_delete(index.row())

        #del self.data[row]
        #self.model().change_data(self.data)

    def set_rows(self):
        """trace the tree to find the selected item"""
        rows = self.find_list_index()

        # TODO: what is this for???
        #if rows != self.old_rows:
        self.old_rows = rows
        #print('rows =', rows)

    def get_row(self) -> tuple[bool, int]:
        """
        gets the row

        Returns
        -------
        is_valid : bool
            is this case valid
        irow : int
            the row index

        row : None or tuple
            None : invalid case
            tuple : valid case
                ('centroid', None, [])
                0 - the location (e.g. node, centroid)
                1 - icase
                2 - []
            index_list : list[int]
                the sorted??? list indices from the tree; 1d

        """
        # if there's only 1 data member, we don't need to extract the data id
        if self.single:
            return True, self.data[0]

        # TODO: what is this for???
        #     crashes some PyQt cases when clicking on the first
        #     non-results level of the sidebar
        #data = deepcopy(self.data)
        is_valid, irow, unused_res_name = self.get_trace()
        index_list = self.find_list_index()
        #index_list.sort()
        #if is_valid and 'NEW' not in unused_res_name:
            #indexs = self.selectedIndexes()
            #for index in indexs:
                #pass
                #x = 1
        return is_valid, irow, index_list

    def get_trace(self) -> tuple[bool, int, str]:
        """
        Returns
        -------
        is_valid : bool
            is this case valid
        row : None or tuple
            None : invalid case
            tuple : valid case
                ('centroid', None, [])
                0 - the location (e.g. node, centroid)
                1 - icase
                2 - []
        word : str
            the 0th index of the tuple; 'centroid' in this case

        """
        is_valid = True
        irow = 0
        data = self.data
        if len(self.old_rows) == 0:
            return False, irow, None

        for row in self.old_rows:
            data_old = data
            try:
                unused_key = data[row][0]
            except IndexError:
                return False, irow, None
            irow = data[row][1]
            data = data[row][2]
        return is_valid, irow, data_old[row][0]

    def set_single(self, single: bool) -> None:
        self.single = single
        self.old_rows = [0]


class ClickTreeView(QTreeView2):
    """
    creates a QTreeView with:
     - all the features of QTreeView2
     - a right click context menu with:
       - Clear Active Results
       - Apply Results to Fringe
       - Apply Results to Displacement
       - Apply Results to Vector
       - Delete Case

    """
    def __init__(self, parent, data, choices,
                 left_click_callback: ClickCallback,
                 right_click_actions: list[tuple[str, ClickCallback, bool]],
                 include_clear: bool=False,
                 include_export_case: bool=False,
                 include_delete: bool=False,
                 include_results: bool=False):
        if right_click_actions is None:
            right_click_actions = []
        QTreeView2.__init__(self, parent, data, choices)
        #
        # TODO: create a menu that only has clear/normals/fringe/delete
        #       if there is no transient result
        #
        right_click_menu = QMenu()

        if include_clear:
            self.clear = right_click_menu.addAction('Clear Results')
            self.clear.triggered.connect(self.on_clear_results)

        if include_export_case:
            self.export_case = right_click_menu.addAction('Export Case')
            self.export_case.triggered.connect(self.on_export_case)

        if include_results:
            self.fringe = right_click_menu.addAction("Apply Results to Fringe")
            self.fringe.triggered.connect(self.on_fringe)

            self.disp = right_click_menu.addAction('Apply Results to Displacement')
            self.vector = right_click_menu.addAction('Apply Results to Vector')

            self.disp.triggered.connect(self.on_disp)
            self.vector.triggered.connect(self.on_vector)

        if include_delete:
            self.delete = right_click_menu.addAction('Delete')
            self.delete.triggered.connect(self.on_delete)

        #self.fringe.setCheckable(True)
        #self.disp.setCheckable(True)
        #self.vector.setCheckable(True)

        self.left_click_callback = left_click_callback

        self._setup_right_click_menu_actions(
            right_click_menu, right_click_actions)
        self.right_click_menu = right_click_menu

    def _setup_right_click_menu_actions(
            self, right_click_menu: QMenu,
            right_click_actions: list[tuple[str, ClickCallback, bool]]) -> None:
        """populates the right click menu"""
        def false_callback(callback_func: ClickCallback):
            """dont validate the click"""
            #print('false')
            #print('  ', callback_func)
            #print('  ', menu)
            callback_func()

        def true_callback(callback_func: ClickCallback):
            """a validated callback returns the row"""
            #print('true')
            #print('  ', callback_func)
            #print('  ', menu)
            #print('  ', self)
            unused_is_valid, icase, imethods = self.get_row()
            #print('callback =', callback_func)
            callback_func(icase)

        for right_click_action in right_click_actions:
            (right_click_msg, callback_func, validate) = right_click_action
            action = right_click_menu.addAction(right_click_msg)

            true_false_callback = true_callback if validate else false_callback
            trigger_func = partial(true_false_callback, callback_func)
            action.triggered.connect(trigger_func)
            #self.clear = right_click_menu.addAction("Clear Results...")
            #self.clear.triggered.connect(self.on_clear_results)

    def on_left_mouse_button(self) -> None:
        """interfaces with other menus"""
        if self.left_click_callback is not None:
            self.set_rows()
            is_valid, icase, imethods = self.get_row()
            self.left_click_callback(icase)
            x = 1

    def on_right_mouse_button(self) -> None:
        """interfaces with the right click menu"""
        self.set_rows()
        is_valid, unused_icase, imethods = self.get_row()
        if not is_valid:
            return
        # TODO: check if we should show disp/vector
        self.right_click_menu.popup(QtGui.QCursor.pos())

    def get_clicked(self) -> dict[str, int]:
        """gets the state of the clickable buttons"""
        is_clicked = {
            'fringe' : self.fringe.isChecked(),
            'disp' : self.disp.isChecked(),
            'vector' : self.vector.isChecked(),
        }
        return is_clicked

    @property
    def sidebar(self) -> Sidebar:
        """gets the sidebar"""
        return self.parent.parent

    @property
    def gui(self) -> MainWindow:
        """get the MainWindow class"""
        return self.sidebar.parent

    def on_clear_results(self) -> None:
        """clears the active result"""
        self.sidebar.on_clear_results()

    def on_export_case(self) -> None:
        """exports the case to a file"""
        unused_is_valid, icase, imethods = self.get_row()
        self.gui.export_case_data(icase)

    def on_fringe(self) -> None:
        """applies a fringe result"""
        unused_is_valid, icase, imethods = self.get_row()
        # imethods are just the default
        self.sidebar.on_fringe(icase)

    def on_disp(self) -> None:
        """applies a displacement result"""
        unused_is_valid, icase, imethods = self.get_row()
        # imethods are just the default
        self.sidebar.on_disp(icase)

    def on_vector(self) -> None:
        """applies a vector result"""
        unused_is_valid, icase, imethods = self.get_row()
        # imethods are just the default
        self.sidebar.on_vector(icase)

    def on_right_mouse_button(self) -> None:
        """interfaces with the right click menu"""
        self.set_rows()
        is_valid, unused_icase, imethods = self.get_row()
        if not is_valid:
            return
        # TODO: check if we should show disp/vector
        self.right_click_menu.popup(QtGui.QCursor.pos())

def get_many_cases(data) -> list[int]:
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
    unused_name, case, rows = data
    if case is None:
        # remove many results
        # (Geometry, None, [results...])
        cases = []
        for unused_irow, row in enumerate(rows):
            unused_name, unused_row_id, unused_data2 = row
            cases += get_many_cases(row)
    else:
        cases = [case]
    return cases


#for subcase in subcases:
#    for time in times:
#        disp
#        stress
#        load
