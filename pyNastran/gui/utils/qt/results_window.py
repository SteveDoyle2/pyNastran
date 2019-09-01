from qtpy import QtGui
from qtpy.QtWidgets import QWidget, QVBoxLayout, QAbstractItemView
from pyNastran.gui.utils.qt.qtreeview2 import RightClickTreeView, GenericRightClickTreeView


class ResultsWindow(QWidget):
    """
    A ResultsWindow creates the box where we actually select our
    results case.  It does not have an apply button.
    """
    def __init__(self, parent, name, data, choices, actions=None,
                 include_clear=True, include_export_case=True,
                 include_delete=True,
                 include_results=True):
        QWidget.__init__(self)
        self.name = name
        self.data = data
        self.choices = choices
        self.parent = parent

        def on_modify(icase):
            print('modify...%i' % icase)
        def on_case(icase):
            print('case...%i' % icase)
        def on_delete():
            print('delete...')

        if actions:
            #actions = [
                ## (right_click_msg, callback, validate?)
                ##('Clear Results...', self.on_clear_results, False),
                ##('Apply Results to Fringe...', 'fringe', self.on_fringe, True),
                ##('Apply Results to Displacement...', self.on_disp, True),
                ##('Apply Results to Vector...', self.on_vector, True),
                #('Delete...', on_delete, False),
                #('Modify...', on_modify, True),
            #]
            self.treeView = GenericRightClickTreeView(
                self, self.data, choices, actions,
                include_clear=include_clear, include_delete=include_delete,
                include_results=include_results)
        else:
            self.treeView = RightClickTreeView(
            self, self.data, choices,
            include_clear=include_clear, include_export_case=include_export_case,
            include_delete=include_delete,
            include_results=include_results, )
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
            raise
            raise RuntimeError('cannot add data=\n%s' % data)
            #if isinstance(data, str):
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
                try:
                    text, i, children = element
                except ValueError:
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
