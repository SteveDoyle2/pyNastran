# -*- coding: utf-8 -*-
from __future__ import print_function, unicode_literals
from six import iteritems

from PyQt4 import QtCore, QtGui
from numpy import setdiff1d, unique, array, hstack, ndarray

from pyNastran.bdf.utils import parse_patran_syntax, parse_patran_syntax_dict
from pyNastran.bdf.cards.collpase_card import collapse_colon_packs
from pyNastran.gui.menus.manage_actors import CustomQTableView, Model


class ColorDisplay(QtGui.QWidget):
    """
    http://stackoverflow.com/questions/4624985/how-simply-display-a-qcolor-using-pyqt
    """
    def __init__(self, parent, default_color=None):
        super(ColorDisplay, self).__init__(parent)
        self.color = default_color
        self.setColor(self.color)

    def setColor(self, color):
        if color is not None:
            color = [int(255 * i) for i in color]
        self.color = QtGui.QColor(*color)
        self.update()

    def paintEvent(self, event=None):
        painter = QtGui.QPainter(self)
        if self.color is not None:
            painter.setBrush(QtGui.QBrush(self.color))
            painter.drawRect(self.rect())

    def getColorName(self):
        return unicode(self.color.name())


class GroupsModify(QtGui.QDialog):
    """
    +--------------------------+
    |     Groups : Modify      |
    +--------------------------+
    |                          |
    |  Name        xxx Default |
    |  Coords      xxx Default |
    |  Elements    xxx Default |
    |  Color       xxx Default |
    |  Add         xxx Add     |
    |  Remove      xxx Remove  |
    |                          |
    |      Set  OK Cancel      |
    +--------------------------+
    """
    def __init__(self, data, win_parent=None, group_active='main'):
        self.win_parent = win_parent
        QtGui.QDialog.__init__(self, win_parent)

        self.win_parent = win_parent
        self.out_data = data

        #print(data)
        self.keys = [group.name for key, group in sorted(iteritems(data))]
        self.active_key = self.keys.index(group_active)

        group_obj = data[self.active_key]
        name = group_obj.name

        self.imain = 0
        self.nrows = len(self.keys)

        self._default_name = group_obj.name
        self._default_elements = group_obj.element_str
        self.elements_pound = group_obj.elements_pound

        self.table = QtGui.QListWidget(parent=None)
        self.table.clear()
        self.table.addItems(self.keys)

        # table
        self.setWindowTitle('Groups: Modify')
        self.create_widgets()
        self.create_layout()
        self.set_connections()

        self.on_set_as_main()
        #self.show()

    def create_widgets(self):
        # Name
        self.name = QtGui.QLabel("Name:")
        self.name_set = QtGui.QPushButton("Set")
        self.name_edit = QtGui.QLineEdit(str(self._default_name).strip())
        self.name_button = QtGui.QPushButton("Default")

        # elements
        self.elements = QtGui.QLabel("Element IDs:")
        self.elements_edit = QtGui.QLineEdit(str(self._default_elements).strip())
        self.elements_button = QtGui.QPushButton("Default")

        # add
        self.add = QtGui.QLabel("Add:")
        self.add_edit = QtGui.QLineEdit(str(''))
        self.add_button = QtGui.QPushButton("Add")

        # remove
        self.remove = QtGui.QLabel("Remove:")
        self.remove_edit = QtGui.QLineEdit(str(''))
        self.remove_button = QtGui.QPushButton("Remove")

        # applies a unique implicitly
        self.eids = parse_patran_syntax(str(self._default_elements), pound=self.elements_pound)

        # closing
        #self.apply_button = QtGui.QPushButton("Apply")
        self.ok_button = QtGui.QPushButton("Close")
        #self.cancel_button = QtGui.QPushButton("Cancel")

        self.set_as_main_button = QtGui.QPushButton("Set As Main")
        self.create_group_button = QtGui.QPushButton('Create New Group')
        self.delete_group_button = QtGui.QPushButton('Delete Group')

        self.name.setEnabled(False)
        self.name_set.setEnabled(False)
        self.name_edit.setEnabled(False)
        self.name_button.setEnabled(False)
        self.elements.setEnabled(False)
        self.elements_button.setEnabled(False)
        self.elements_edit.setEnabled(False)
        self.add.setEnabled(False)
        self.add_button.setEnabled(False)
        self.add_edit.setEnabled(False)
        self.remove.setEnabled(False)
        self.remove_button.setEnabled(False)
        self.remove_edit.setEnabled(False)
        self.delete_group_button.setEnabled(False)
        #self.apply_button.setEnabled(False)
        #self.ok_button.setEnabled(False)


    def create_layout(self):
        grid = QtGui.QGridLayout()
        grid.addWidget(self.name, 0, 0)
        grid.addWidget(self.name_edit, 0, 1)
        grid.addWidget(self.name_set, 0, 2)
        grid.addWidget(self.name_button, 0, 3)

        grid.addWidget(self.elements, 2, 0)
        grid.addWidget(self.elements_edit, 2, 1)
        grid.addWidget(self.elements_button, 2, 2)

        grid.addWidget(self.add, 4, 0)
        grid.addWidget(self.add_edit, 4, 1)
        grid.addWidget(self.add_button, 4, 2)

        grid.addWidget(self.remove, 5, 0)
        grid.addWidget(self.remove_edit, 5, 1)
        grid.addWidget(self.remove_button, 5, 2)


        ok_cancel_box = QtGui.QHBoxLayout()
        #ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        #ok_cancel_box.addWidget(self.cancel_button)


        main_create_delete = QtGui.QHBoxLayout()
        main_create_delete.addWidget(self.set_as_main_button)
        main_create_delete.addWidget(self.create_group_button)
        main_create_delete.addWidget(self.delete_group_button)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.table)
        vbox.addLayout(grid)
        vbox.addLayout(main_create_delete)
        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def on_set_name(self):
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

    def set_connections(self):
        self.connect(self.name_set, QtCore.SIGNAL('clicked()'), self.on_set_name)
        self.connect(self.name_button, QtCore.SIGNAL('clicked()'), self.on_default_name)
        self.connect(self.elements_button, QtCore.SIGNAL('clicked()'), self.on_default_elements)

        self.connect(self.add_button, QtCore.SIGNAL('clicked()'), self.on_add)
        self.connect(self.remove_button, QtCore.SIGNAL('clicked()'), self.on_remove)
        self.connect(self.table, QtCore.SIGNAL('itemClicked(QListWidgetItem *)'), self.on_update_active_key)

        #self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        #self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)

        self.connect(self.set_as_main_button, QtCore.SIGNAL('clicked()'), self.on_set_as_main)
        self.connect(self.create_group_button, QtCore.SIGNAL('clicked()'), self.on_create_group)
        self.connect(self.delete_group_button, QtCore.SIGNAL('clicked()'), self.on_delete_group)

    def on_create_group(self):
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

        self.keys = [group.name for key, group in sorted(iteritems(self.out_data))]
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
        if self.active_key == 0:
            return

        #self.deleted_groups.add(self.imain)
        items = {}
        j = 0
        for i, key in sorted(iteritems(self.out_data)):
            if i != self.active_key:
                items[j] = key
                j += 1

        # update internal parameters
        self.out_data = items
        if self.imain >= self.active_key:
            self.imain = max(0, self.imain - 1)
        self.active_key = max(0, self.active_key - 1)
        self.nrows -= 1
        self.keys = [group.name for key, group in sorted(iteritems(items))]

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
        if self.win_parent is not None:
            # we're not testing the menu
            self.win_parent.post_group(group)

    def closeEvent(self, event):
        event.accept()

    def on_add(self):
        eids, is_valid = self.check_patran_syntax(self.add_edit, pound=self.elements_pound)
        #adict, is_valid = self.check_patran_syntax_dict(self.add_edit)
        if not is_valid:
            #self.add_edit.setStyleSheet("QLineEdit{background: red;}")
            return

        self.eids = unique(hstack([self.eids, eids]))
        #self.eids = _add(adict, ['e', 'elem', 'element'], self.eids)
        #self.cids = _add(adict, ['c', 'cid', 'coord'], self.cids)
        self._apply_cids_eids()

        self.add_edit.clear()
        self.add_edit.setStyleSheet("QLineEdit{background: white;}")

    def _apply_cids_eids(self):
        #ctext = _get_collapsed_text(self.cids)
        etext = _get_collapsed_text(self.eids)

        #self.coords_edit.setText(str(ctext.lstrip()))
        self.elements_edit.setText(str(etext.lstrip()))
        self.out_data[self.active_key].element_ids = self.eids

    def on_remove(self):
        eids, is_valid = self.check_patran_syntax(self.remove_edit)
        #adict, is_valid = self.check_patran_syntax_dict(self.remove_edit)
        if not is_valid:
            #self.remove_edit.setStyleSheet("QLineEdit{background: red;}")
            return

        #self.eids = _remove(adict, ['e', 'elem', 'element'], self.eids)
        #self.cids = _remove(adict, ['c', 'cid', 'coord'], self.cids)
        self.eids = setdiff1d(self.eids, eids)
        self._apply_cids_eids()

        self.remove_edit.clear()
        self.remove_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_name(self):
        name = str(self._default_name)
        self.name_edit.setText(name)
        self.name_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_elements(self):
        element_str = str(self._default_elements)
        self.elements_edit.setText(element_str)
        self.elements_edit.setStyleSheet("QLineEdit{background: white;}")
        group = self.out_data[self.active_key]
        group.element_str = element_str

    def check_patran_syntax(self, cell, pound=None):
        text = str(cell.text())
        try:
            values = parse_patran_syntax(text, pound=pound)
            cell.setStyleSheet("QLineEdit{background: white;}")
            return values, True
        except ValueError as e:
            cell.setStyleSheet("QLineEdit{background: red;}")
            cell.setToolTip(str(e))
            return None, False

    def check_patran_syntax_dict(self, cell, pound=None):
        text = str(cell.text())
        try:
            value = parse_patran_syntax_dict(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
            cell.setToolTip('')
            return value, True
        except (ValueError, SyntaxError, KeyError) as e:
            cell.setStyleSheet("QLineEdit{background: red;}")
            cell.setToolTip(str(e))
            return None, False

    def check_float(self, cell):
        text = cell.text()
        try:
            value = float(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
            cell.setToolTip('')
            return value, True
        except ValueError:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    def check_name(self, cell):
        text = str(cell.text()).strip()
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

    def check_format(self, cell):
        text = str(cell.text())

        is_valid = True
        if len(text) < 2:
            is_valid = False
        elif 's' in text.lower():
            is_valid = False
        elif '%' not in text[0]:
            is_valid = False
        elif text[-1].lower() not in ['g', 'f', 'i', 'e']:
            is_valid = False

        try:
            text % 1
            text % .2
            text % 1e3
            text % -5.
            text % -5
        except ValueError:
            is_valid = False

        try:
            text % 's'
            is_valid = False
        except TypeError:
            pass

        if is_valid:
            cell.setStyleSheet("QLineEdit{background: white;}")
            return text, True
        else:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    def on_validate(self):
        name, flag0 = self.check_name(self.name_edit)
        elements, flag1 = self.check_patran_syntax(self.elements_edit,
                                                         pound=self.elements_pound)
        #coords_value, flag2 = self.check_patran_syntax(self.coords_edit,
        #pound=self.coords_pound)

        if flag0 and flag1:
            self._default_name = name
            self._default_elements = self.eids
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_apply(self, force=False):
        passed = self.on_validate()
        if passed or force:
            self.win_parent.on_modify_group(self.out_data)

    def on_ok(self):
        passed = self.on_validate()
        if passed:
            self.out_data['close'] = True
            self.out_data['clicked_ok'] = True
            self.out_data['clicked_cancel'] = False
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.out_data['close'] = True
        self.out_data['clicked_cancel'] = True
        self.close()

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.on_cancel()

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
            self.name.setEnabled(False)
            self.name_set.setEnabled(False)
            self.name_edit.setEnabled(False)
            self.name_button.setEnabled(False)
            self.elements.setEnabled(False)
            self.elements_button.setEnabled(False)
            self.elements_edit.setEnabled(False)
            self.add.setEnabled(False)
            self.add_button.setEnabled(False)
            self.add_edit.setEnabled(False)
            self.remove.setEnabled(False)
            self.remove_button.setEnabled(False)
            self.remove_edit.setEnabled(False)
            self.delete_group_button.setEnabled(False)
            if name == 'anti-main':
                self.set_as_main_button.setEnabled(False)
            #self.apply_button.setEnabled(False)
            #self.ok_button.setEnabled(False)
        else:
            self.name.setEnabled(True)
            self.name_set.setEnabled(True)
            self.name_edit.setEnabled(True)
            self.name_button.setEnabled(True)
            self.elements.setEnabled(True)
            self.elements_button.setEnabled(True)

            self.add.setEnabled(True)
            self.add_button.setEnabled(True)
            self.add_edit.setEnabled(True)
            self.remove.setEnabled(True)
            self.remove_button.setEnabled(True)
            self.remove_edit.setEnabled(True)
            self.delete_group_button.setEnabled(True)
            #self.apply_button.setEnabled(True)
            #self.ok_button.setEnabled(True)
        # TODO: call default
        #self.elements_edit # obj.eids


def _get_collapsed_text(values):
    singles, doubles = collapse_colon_packs(values)
    text = ' '.join([str(s) for s in singles]) + ' '
    text += ' '.join([''.join([str(doublei) for doublei in double]) for double in doubles])
    return text

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

def _simple_add(values_to_add):
    value_stack.append(values_to_add)
    values_add = unique(hstack(value_stack))
    return values_add

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


class Group(object):
    def __init__(self, name, element_str, elements_pound, editable=True):
        if len(name):
            assert len(name) > 0, name
            assert name[-1] != ' ', name
            assert '\n' not in name, name
            assert '\r' not in name, name
            assert '\t' not in name, name
        self.name = name
        #self.cids = [0]
        assert isinstance(element_str, (str, unicode)), element_str
        self.element_str = element_str
        self.elements_pound = elements_pound
        self.editable = editable

    @property
    def element_ids(self):
        return parse_patran_syntax(self.element_str, pound=self.elements_pound)

    @element_ids.setter
    def element_ids(self, eids):
        assert isinstance(eids, ndarray), eids
        self.element_str = _get_collapsed_text(eids).strip()

    def __repr__(self):
        msg = 'Group:\n'
        msg += '  name: %s\n' % self.name
        msg += '  editable: %s\n' % self.editable
        #msg += '  cids: [%s]\n' % _get_collapsed_text(self.cids).strip()
        msg += '  element_str: [%s]\n' % self.element_str
        #msg += '  elements: [%s]\n' % _get_collapsed_text(self.elements).strip()
        msg += '  elements_pound: %s\n' % self.elements_pound
        return msg

def main():
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QtGui.QApplication(sys.argv)

    d = {
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
    main_window = GroupsModify(d)
    main_window.show()
    app.exec_()

if __name__ == "__main__":
    main()
