"""
defines:
 - GroupsModify

"""
# -*- coding: utf-8 -*-
from numpy import setdiff1d, unique, hstack

from qtpy import QtGui
from qtpy.QtWidgets import (
    QLabel, QLineEdit, QPushButton, QApplication,
    QListWidget, QGridLayout, QHBoxLayout, QVBoxLayout,
)

from pyNastran.bdf.utils import parse_patran_syntax #, parse_patran_syntax_dict
#from pyNastran.gui.menus.manage_actors import Model
from pyNastran.gui.utils.qt.pydialog import PyDialog, check_patran_syntax
from pyNastran.gui.utils.qt.qelement_edit import QElementEdit

from pyNastran.gui.menus.groups_modify.groups import Group, _get_collapsed_text
#from .groups_modify.color_display import ColorDisplay
from pyNastran.gui.utils.vtk.gui_utils import add_actors_to_gui, remove_actors_from_gui
from pyNastran.gui.menus.highlight.highlight import create_highlighted_actors

#import pyNastran
from pyNastran.gui import ICON_PATH


class GroupsModify(PyDialog):
    """
    +-------------------------------+
    |     Groups : Modify           |
    +-------------------------------+
    |                               |
    |  Name        xxx Default      |
    |  Nodes       xxx Default Show |
    |  Color       xxx Default      |
    |  Add         xxx Add     Show |
    |  Remove      xxx Remove  Show |
    |                               |
    |      Set  OK Cancel           |
    +-------------------------------+
    """
    def __init__(self, data, win_parent=None, model_name=None, group_active='main'):
        super(GroupsModify, self).__init__(data, win_parent)
        self.set_font_size(data['font_size'])
        self._updated_groups = False
        self.model_name = model_name
        self.actors = []

        #self.out_data = data

        #print(data)
        self.keys = [group.name for key, group in sorted(data.items())
                     if isinstance(key, int)]
        self.active_key = self.keys.index(group_active)

        group_obj = data[self.active_key]
        unused_name = group_obj.name

        self.imain = 0
        self.nrows = len(self.keys)

        self._default_name = group_obj.name
        self._default_elements = group_obj.element_str
        self.elements_pound = group_obj.elements_pound

        self.table = QListWidget(parent=None)
        self.table.clear()
        self.table.addItems(self.keys)

        self.setWindowTitle('Groups: Modify')
        self.create_widgets()
        self.create_layout()
        self.set_connections()

        self.on_set_as_main()

    def create_widgets(self):
        """creates the menu objects"""
        #icon = QtGui.QPixmap(os.path.join(ICON_PATH, 'node.png'))
        #icon = QtGui.QPixmap(os.path.join(ICON_PATH, 'element.png'))
        # Name
        self.pick_style_label = QLabel('Pick Style:')
        #self.pick_node_button = QPushButton('Node')
        self.pick_element_button = QPushButton('Element')
        self.pick_area_button = QPushButton('Area')
        #self.pick_node_button.setIcon(icon)
        #self.pick_area_button.setIcon(icon)

        # Name
        self.name_label = QLabel('Name:')
        self.name_set = QPushButton('Set')
        self.name_edit = QLineEdit(str(self._default_name).strip())
        self.name_button = QPushButton('Default')

        # elements
        self.elements_label = QLabel('Element IDs:')
        self.elements_edit = QLineEdit(str(self._default_elements).strip())
        self.elements_button = QPushButton('Default')
        self.elements_highlight_button = QPushButton('Show')

        # add
        self.add_label = QLabel('Add Elements:')
        self.add_edit = QElementEdit(self, self.model_name, pick_style='area')
        self.add_button = QPushButton('Add')
        self.add_highlight_button = QPushButton('Show')

        # remove
        self.remove_label = QLabel('Remove Elements:')
        self.remove_edit = QElementEdit(self, self.model_name, pick_style='area')
        self.remove_button = QPushButton('Remove')
        self.remove_highlight_button = QPushButton('Show')

        # applies a unique implicitly
        self.eids = parse_patran_syntax(str(self._default_elements), pound=self.elements_pound)

        # closing
        #self.apply_button = QPushButton('Apply')
        self.ok_button = QPushButton('Close')
        #self.cancel_button = QPushButton('Cancel')

        self.set_as_main_button = QPushButton('Set As Main')
        self.create_group_button = QPushButton('Create New Group')
        self.delete_group_button = QPushButton('Delete Group')

        self.name_label.setEnabled(False)
        self.name_set.setEnabled(False)
        self.name_edit.setEnabled(False)
        self.name_button.setEnabled(False)

        self.elements_label.setEnabled(False)
        self.elements_button.setEnabled(False)
        self.elements_edit.setEnabled(False)
        self.elements_highlight_button.setEnabled(False)

        self.add_label.setEnabled(False)
        self.add_button.setEnabled(False)
        self.add_edit.setEnabled(False)
        self.add_highlight_button.setEnabled(False)

        self.remove_label.setEnabled(False)
        self.remove_button.setEnabled(False)
        self.remove_edit.setEnabled(False)
        self.remove_highlight_button.setEnabled(False)

        #self.apply_button.setEnabled(False)
        #self.ok_button.setEnabled(False)


    def create_layout(self):
        """displays the menu objects"""
        hbox = QHBoxLayout()
        hbox.addWidget(self.pick_style_label)
        hbox.addWidget(self.pick_element_button)
        hbox.addWidget(self.pick_area_button)
        hbox.addStretch()

        grid = QGridLayout()
        irow = 0
        grid.addWidget(self.name_label, irow, 0)
        grid.addWidget(self.name_edit, irow, 1)
        grid.addWidget(self.name_set, irow, 2)
        grid.addWidget(self.name_button, irow, 3)
        irow += 1

        grid.addWidget(self.elements_label, irow, 0)
        grid.addWidget(self.elements_edit, irow, 1)
        grid.addWidget(self.elements_button, irow, 2)
        grid.addWidget(self.elements_highlight_button, irow, 3)
        irow += 1

        grid.addWidget(self.add_label, irow, 0)
        grid.addWidget(self.add_edit, irow, 1)
        grid.addWidget(self.add_button, irow, 2)
        grid.addWidget(self.add_highlight_button, irow, 3)
        irow += 1

        grid.addWidget(self.remove_label, irow, 0)
        grid.addWidget(self.remove_edit, irow, 1)
        grid.addWidget(self.remove_button, irow, 2)
        grid.addWidget(self.remove_highlight_button, irow, 3)
        irow += 1

        ok_cancel_box = QHBoxLayout()
        #ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        #ok_cancel_box.addWidget(self.cancel_button)


        main_create_delete = QHBoxLayout()
        main_create_delete.addWidget(self.set_as_main_button)
        main_create_delete.addWidget(self.create_group_button)
        main_create_delete.addWidget(self.delete_group_button)

        vbox = QVBoxLayout()
        vbox.addWidget(self.table)
        vbox.addLayout(hbox)
        vbox.addLayout(grid)
        vbox.addLayout(main_create_delete)
        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def on_set_name(self):
        self.remove_highlight_actor()
        name = str(self.name_edit.text()).strip()
        if name not in self.keys:
            self.name_edit.setStyleSheet("QLineEdit{background: white;}")
            group = self.out_data[self.active_key]
            group.name = name
            self.keys[self.active_key] = name
            self.recreate_table()
        elif name != self.keys[self.active_key]:
            self.name_edit.setStyleSheet("QLineEdit{background: red;}")
        elif name == self.keys[self.active_key]:
            self.name_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_pick_element(self):
        self.add_edit.pick_style = 'single'
        self.remove_edit.pick_style = 'single'

    def on_pick_area(self):
        self.add_edit.pick_style = 'area'
        self.remove_edit.pick_style = 'area'


    def set_connections(self):
        """creates the actions for the menu"""
        self.pick_element_button.clicked.connect(self.on_pick_element)
        self.pick_area_button.clicked.connect(self.on_pick_area)

        self.name_set.clicked.connect(self.on_set_name)
        self.name_button.clicked.connect(self.on_default_name)
        self.elements_button.clicked.connect(self.on_default_elements)
        self.elements_highlight_button.clicked.connect(self.on_highlight_elements)

        self.add_button.clicked.connect(self.on_add)
        self.add_highlight_button.clicked.connect(self.on_highlight_add)

        self.remove_button.clicked.connect(self.on_remove)
        self.remove_highlight_button.clicked.connect(self.on_highlight_remove)

        self.table.itemClicked.connect(self.on_update_active_key)
        self.ok_button.clicked.connect(self.on_ok)

        self.set_as_main_button.clicked.connect(self.on_set_as_main)
        self.create_group_button.clicked.connect(self.on_create_group)
        self.delete_group_button.clicked.connect(self.on_delete_group)

    def on_create_group(self):
        self.remove_highlight_actor()

        irow = self.nrows
        new_key = 'Group %s' % irow
        while new_key in self.keys:
            irow += 1
            new_key = 'Group %s' % irow
        irow = self.nrows

        self.keys.append(new_key)
        group = Group(
            new_key,
            element_str='',
            elements_pound=self.elements_pound,
            editable=True)
        self.out_data[irow] = group

        self.table.reset()
        self.table.addItems(self.keys)
        self.nrows += 1

        #----------------------------------
        # update internal parameters
        #self.out_data = items
        if self.imain > self.active_key:
            self.imain += 1

        #make the new group the default
        self.active_key = self.nrows - 1

        self.keys = [group.name for key, group in sorted(self.out_data.items())
                     if isinstance(key, int)]
        self.recreate_table()

    def recreate_table(self):
        # update gui
        self.table.clear()
        self.table.addItems(self.keys)
        item = self.table.item(self.imain)

        bold = QtGui.QFont()
        bold.setBold(True)
        bold.setItalic(True)
        item.setFont(bold)
        self.table.update()

        # update key
        name = self.keys[self.active_key]
        self._update_active_key_by_name(name)

    def on_delete_group(self):
        self.remove_highlight_actor()
        if self.active_key == 0:
            return

        #self.deleted_groups.add(self.imain)
        items = {}
        j = 0
        for i, key in sorted(self.out_data.items()):
            if isinstance(i, int):
                continue
            if i != self.active_key:
                items[j] = key
                j += 1

        # update internal parameters
        self.out_data = items
        if self.imain >= self.active_key:
            self.imain = max(0, self.imain - 1)
        self.active_key = max(0, self.active_key - 1)
        self.nrows -= 1
        self.keys = [group.name for key, group in sorted(items.items())]

        self.recreate_table()

        # update key
        name = self.keys[self.active_key]
        self._update_active_key_by_name(name)

    def on_set_as_main(self):
        bold = QtGui.QFont()
        bold.setBold(True)
        bold.setItalic(True)

        normal = QtGui.QFont()
        normal.setBold(False)
        normal.setItalic(False)

        obj = self.table.item(self.imain)
        obj.setFont(normal)

        self.imain = self.active_key
        obj = self.table.item(self.imain)
        obj.setFont(bold)
        group = self.out_data[self.imain]
        self._default_elements = group.element_str
        self._default_name = group.name
        self.on_update_main()

    def on_update_main(self):
        """adds/removes the elements to the main actor when add/remove is pressed"""
        group = self.out_data[self.imain]
        if self._default_name == group.name and self.win_parent is not None:
            # we're not testing the menu
            self.win_parent.post_group(group, update_groups=True)

    def closeEvent(self, event):
        """closes the window"""
        self.remove_highlight_actor()
        self.out_data['close'] = True
        event.accept()

    def remove_highlight_actor(self):
        """removes the highlighted actor"""
        gui = self.win_parent
        remove_actors_from_gui(gui, self.actors, render=True)
        self.actors = []

    def on_highlight(self, nids=None, eids=None):
        """highlights the nodes"""
        gui = self.win_parent
        unused_name = self.keys[self.active_key]

        if gui is not None:
            self.remove_highlight_actor()
            ## TODO: super strange; doesn't work...
            mouse_actions = gui.mouse_actions
            grid = mouse_actions.get_grid_selected(self.model_name)
            #all_nodes = mouse_actions.node_ids
            all_elements = mouse_actions.element_ids

            actors = create_highlighted_actors(gui, grid,
                                               all_nodes=None, nodes=nids,
                                               all_elements=all_elements, elements=eids,
                                               add_actors=False)
            if actors:
                add_actors_to_gui(gui, actors, render=True)
                self.actors = actors

    def on_highlight_elements(self):
        """highlights the active elements"""
        self.on_highlight(eids=self.eids)

    def on_highlight_add(self):
        """highlights the elements to add"""
        eids, is_valid = check_patran_syntax(self.add_edit, pound=self.elements_pound)
        if not is_valid:
            #self.add_edit.setStyleSheet("QLineEdit{background: red;}")
            return
        self.on_highlight(eids=eids)

    def on_highlight_remove(self):
        """highlights the elements to remove"""
        eids, is_valid = check_patran_syntax(self.remove_edit)
        if not is_valid:
            #self.remove_edit.setStyleSheet("QLineEdit{background: red;}")
            return
        self.on_highlight(eids=eids)

    def on_add(self):
        self.remove_highlight_actor()
        eids, is_valid = check_patran_syntax(self.add_edit, pound=self.elements_pound)
        #adict, is_valid = check_patran_syntax_dict(self.add_edit)
        if not is_valid:
            #self.add_edit.setStyleSheet("QLineEdit{background: red;}")
            return

        self.eids = unique(hstack([self.eids, eids]))
        #self.eids = _add(adict, ['e', 'elem', 'element'], self.eids)
        #self.cids = _add(adict, ['c', 'cid', 'coord'], self.cids)
        self._apply_cids_eids()

        self.add_edit.clear()
        self.add_edit.setStyleSheet("QLineEdit{background: white;}")
        self.on_update_main()

    def _apply_cids_eids(self):
        #ctext = _get_collapsed_text(self.cids)
        etext = _get_collapsed_text(self.eids)

        #self.coords_edit.setText(str(ctext.lstrip()))
        self.elements_edit.setText(str(etext.lstrip()))
        self.out_data[self.active_key].element_ids = self.eids

    def on_remove(self):
        self.remove_highlight_actor()
        eids, is_valid = check_patran_syntax(self.remove_edit)
        #adict, is_valid = check_patran_syntax_dict(self.remove_edit)
        if not is_valid:
            #self.remove_edit.setStyleSheet("QLineEdit{background: red;}")
            return

        #self.eids = _remove(adict, ['e', 'elem', 'element'], self.eids)
        #self.cids = _remove(adict, ['c', 'cid', 'coord'], self.cids)
        self.eids = setdiff1d(self.eids, eids)
        self._apply_cids_eids()

        self.remove_edit.clear()
        self.remove_edit.setStyleSheet("QLineEdit{background: white;}")
        self.on_update_main()

    def on_default_name(self):
        self.remove_highlight_actor()
        name = str(self._default_name)
        self.name_edit.setText(name)
        self.name_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_elements(self):
        self.remove_highlight_actor()
        element_str = str(self._default_elements)
        self.elements_edit.setText(element_str)
        self.elements_edit.setStyleSheet("QLineEdit{background: white;}")
        group = self.out_data[self.active_key]
        group.element_str = element_str

    def check_name(self, cell):
        cell_value = cell.text()
        try:
            text = str(cell_value).strip()
        except UnicodeEncodeError:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

        if len(text):
            cell.setStyleSheet("QLineEdit{background: white;}")
            return text, True
        else:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

        if self._default_name != text:
            if self._default_name in self.out_data:
                cell.setStyleSheet("QLineEdit{background: white;}")
                return text, True
            else:
                cell.setStyleSheet("QLineEdit{background: red;}")
                return None, False

    def on_validate(self):
        name, flag0 = self.check_name(self.name_edit)
        unused_elements, flag1 = check_patran_syntax(
            self.elements_edit, pound=self.elements_pound)
        #coords_value, flag2 = check_patran_syntax(
            #self.coords_edit, pound=self.coords_pound)

        if all([flag0, flag1]):
            self._default_name = name
            self._default_elements = self.eids
            self.out_data['clicked_ok'] = True
            self.out_data['close'] = True
            return True
        return False

    def on_apply(self, force=False):
        passed = self.on_validate()
        if passed or force:
            self.win_parent._apply_modify_groups(self.out_data)
        return passed

    def on_ok(self):
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.remove_highlight_actor()
        self.out_data['close'] = True
        self.close()

    def on_update_active_key(self, index):
        self.update_active_key(index)
        #str(index.text())

    def update_active_key(self, index):
        #old_obj = self.out_data[self.imain]
        name = str(index.text())
        self._update_active_key_by_name(name)

    def _update_active_key_by_name(self, name):
        if name in self.keys:
            self.active_key = self.keys.index(name)
        else:
            # we (hopefully) just removed a row
            #self.active_key = self.keys[self.active_key]
            pass

        self.name_edit.setText(name)
        obj = self.out_data[self.active_key]

        self.eids = parse_patran_syntax(obj.element_str, pound=obj.elements_pound)
        self._default_elements = obj.element_str
        self._default_name = name
        self._apply_cids_eids()

        self.set_as_main_button.setEnabled(True)

        if name in ['main', 'anti-main']:
            enabled = False
            self.elements_edit.setEnabled(False)
            if name == 'anti-main':
                self.set_as_main_button.setEnabled(False)
            #self.apply_button.setEnabled(False)
            #self.ok_button.setEnabled(False)
        else:
            enabled = True
            self.elements_highlight_button.setEnabled(True)
            self.add_label.setEnabled(True)
            self.add_highlight_button.setEnabled(True)
            self.remove_highlight_button.setEnabled(True)
            #self.apply_button.setEnabled(True)
            #self.ok_button.setEnabled(True)

        self.name_label.setEnabled(enabled)
        self.name_set.setEnabled(enabled)
        self.name_edit.setEnabled(enabled)
        self.name_button.setEnabled(enabled)
        self.elements_label.setEnabled(enabled)
        self.elements_button.setEnabled(enabled)
        self.add_button.setEnabled(enabled)
        self.add_edit.setEnabled(enabled)

        self.remove_label.setEnabled(enabled)
        self.remove_button.setEnabled(enabled)
        self.remove_edit.setEnabled(enabled)
        self.delete_group_button.setEnabled(enabled)

        # TODO: call default
        #self.elements_edit # obj.eids



def _add(adict, keys, values_to_add):
    value_stack = []
    for key in keys:
        if key not in adict:
            continue
        values = adict[key]
        value_stack.append(values)
    if value_stack:
        value_stack.append(values_to_add)
        values_add = unique(hstack(value_stack))
        return values_add
    return values_to_add

#def _simple_add(values_to_add):
    #value_stack.append(values_to_add)
    #values_add = unique(hstack(value_stack))
    #return values_add

def _remove(adict, keys, values_to_remove):
    value_stack = []
    for key in keys:
        if key not in adict:
            continue
        value_stack.append(adict[key])
    if value_stack:
        values_remove = unique(hstack(value_stack))
        return setdiff1d(values_to_remove, values_remove)
    return values_to_remove


def main(): # pragma: no cover
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import os
    os.environ['QT_API'] = 'pyside'
    #os.environ['QT_API'] = 'pyqt5'

    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    #The Main window

    data = {
        'font_size' : 8,
        0 : Group(name='main', element_str='1:#', elements_pound='103', editable=False),
        1 : Group(name='wing', element_str='1:10', elements_pound='103', editable=True),
        2 : Group(name='fuselage1', element_str='50:60', elements_pound='103', editable=True),
        3 : Group(name='fuselage2', element_str='50:60', elements_pound='103', editable=True),
        4 : Group(name='fuselage3', element_str='50:60', elements_pound='103', editable=True),
        5 : Group(name='fuselage4', element_str='50:60', elements_pound='103', editable=True),
        6 : Group(name='fuselage5', element_str='50:60', elements_pound='103', editable=True),
        7 : Group(name='fuselage6', element_str='50:60', elements_pound='103', editable=True),
        8 : Group(name='fuselage7', element_str='50:60', elements_pound='103', editable=True),
        9 : Group(name='fuselage8', element_str='50:60', elements_pound='103', editable=True),
        10 : Group(name='fuselage9', element_str='50:60', elements_pound='103', editable=True),
        11 : Group(name='fuselage10', element_str='50:60', elements_pound='103', editable=True),
        12 : Group(name='fuselage11', element_str='50:60', elements_pound='103', editable=True),
        13 : Group(name='fuselage12', element_str='50:60', elements_pound='103', editable=True),
        14 : Group(name='fuselage13', element_str='50:60', elements_pound='103', editable=True),
        15 : Group(name='fuselage14', element_str='50:60', elements_pound='103', editable=True),
    }
    main_window = GroupsModify(data, win_parent=None, group_active='main')
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == '__main__':  # pragma: no cover
    main()
