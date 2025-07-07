from typing import Optional, Any
from qtpy import QtGui
from qtpy.QtWidgets import QWidget, QVBoxLayout, QAbstractItemView
#from qtpy.QtCore import QAbstractView
from pyNastran.gui.utils.qt.qtreeview2 import ClickTreeView, ClickCallback

class ResultsWindow(QWidget):
    """
    A ResultsWindow creates the box where we actually select our
    results case.  It does not have an apply button.
    """
    def __init__(self, parent, name: str,
                 data: list[Any],
                 choices,
                 is_single_select: bool,
                 left_click_callback: Optional[ClickCallback],
                 right_click_actions: list[tuple[str, ClickCallback, bool]]=None,
                 include_clear: bool=True,
                 include_export_case: bool=True,
                 include_delete: bool=True,
                 include_results: bool=True):
        """
        Parameters
        ----------
        right_click_actions : varies
            None:
                use the default right_click_actions
            list[tuple[str, Any, bool]]:
                action : (name, function, return_icase)
                name : str
                    the name of the action
                function : the callback function of the form:
                    def return_icase(icase):
                        pass
                    def dont_return_icase():
                        pass
                    the chosen function (return_icase/dont_return_icase) is determined by
                    return icase
                return_icase : bool
                    selects the corresponding function

        """
        QWidget.__init__(self)
        self.name = name
        self.data = data
        self.choices = choices
        self.parent = parent

        #def on_modify(icase: int):
            #print('modify...%i' % icase)
        #def on_case(icase: int):
            #print('case...%i' % icase)
        #def on_delete():
            #print('delete...')

        #right_click_actions = []
        #right_click_actions = [
            ## (right_click_msg, callback, validate?)
            ##('Clear Results...', self.on_clear_results, False),
            ##('Apply Results to Fringe...', 'fringe', self.on_fringe, True),
            ##('Apply Results to Displacement...', self.on_disp, True),
            ##('Apply Results to Vector...', self.on_vector, True),
            #('Delete...', on_delete, False),
            #('Modify...', on_modify, True),
        #]
        self.tree_view = ClickTreeView(
            self, self.data, choices,
            left_click_callback,
            right_click_actions,
            include_export_case=include_export_case,
            include_clear=include_clear,
            include_delete=include_delete,
            include_results=include_results)
        self.tree_view.setEditTriggers(QAbstractItemView.NoEditTriggers)

        self.model = QtGui.QStandardItemModel()
        is_single = self.addItems(self.model, data)
        self.tree_view.setModel(self.model)
        self.tree_view.set_single(is_single)

        assert isinstance(is_single_select, bool), is_single_select
        if is_single_select:
            self.tree_view.setSelectionMode(QAbstractItemView.SingleSelection)
        else:
            self.tree_view.setSelectionMode(QAbstractItemView.MultiSelection)
            #item = tree.invisibleRootItem()
            #def select_item(item)
                #item.setSelected(True)
                #for i in range(item.childCount()):
                    #child = item.child(i)
                    #select_item(child)
            #select_item(item)
            #widget.selectAll()
            self.tree_view.selectAll()
            #youQTreeWidget.setSelectionMode(QGui.QAbstractView.MultiSelection)

        self.model.setHorizontalHeaderLabels([self.tr(self.name)])

        layout = QVBoxLayout()
        layout.addWidget(self.tree_view)
        self.setLayout(layout)

    def update_data(self, data):
        self.clear_data()
        self.data = data
        try:
            self.addItems(self.model, data)
        except Exception:
            raise
            raise RuntimeError('cannot add data=\n%s' % data)
            #if isinstance(data, str):
                #self.addItems(self.model, data)
            #else:
                #self.addItems(self.model, *tuple(data))
        self.tree_view.data = data
        #layout = QVBoxLayout()
        #layout.addWidget(self.tree_view)
        #self.setLayout(layout)

    def clear_data(self):
        self.model.clear()
        self.tree_view.data = []
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
                try:
                    text, i, children = element
                except ValueError:
                    #  [
                    #     ('Point Data', None, [
                    #         ('NodeID', 0, []),
                    #         ('Displacement T_XYZ_subcase=1', 1, []),
                    #         ('Displacement R_XYZ_subcase=1', 2, []),
                    #     ])
                    #  ]
                    #
                    # should be:
                    #   ('Point Data', None, [
                    #       ('NodeID', 0, []),
                    #       ('Displacement T_XYZ_subcase=1', 1, []),
                    #       ('Displacement R_XYZ_subcase=1', 2, []),
                    #   ])
                    print('failed element = ', element)
                    raise
                nchildren = len(children)
                #print('text=%r' % text)
                item = QtGui.QStandardItem(text)
                parent.appendRow(item)

                # TODO: count_check and ???
                if nelements == 1 and nchildren == 0 and level == 0:
                    #self.result_data_window.setEnabled(False)
                    item.setEnabled(False)
                    #print(dir(self.tree_view))
                    #self.tree_view.setCurrentItem(self, 0)
                    #item.mousePressEvent(None)
                    redo = True
                #else:
                    #pass
                    #print('item=%s count_check=%s nelements=%s nchildren=%s' % (
                        #text, count_check, nelements, nchildren))
                if children:
                    assert isinstance(children, list), children
                    self.addItems(item, children, level + 1, count_check=count_check)
                    #print('*children = %s' % children)
            is_single = redo
            return is_single
        except ValueError:
            print()
            print(f'elements = {elements}')
            print(f'element = {element}')
            print(f'len(element) = {len(element)}')
            print(f'len(elements)={len(elements)}')
            for elem in elements:
                print('  e = %s' % str(elem))
            raise
        #if redo:
        #    data = [
        #        ('A', []),
        #        ('B', []),
        #    ]
        #    self.update_data(data)
