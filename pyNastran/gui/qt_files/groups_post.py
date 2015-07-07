# -*- coding: utf-8 -*-
from __future__ import print_function, unicode_literals
from PyQt4 import QtCore, QtGui

from numpy import setdiff1d, unique, array, hstack, argsort

from pyNastran.bdf.utils import parse_patran_syntax, parse_patran_syntax_dict
from pyNastran.bdf.cards.baseCard import collapse_colon_packs
from pyNastran.gui.qt_files.groups_modify import ColorDisplay, _get_collapsed_text

class GroupsPostView(QtGui.QDialog):
    """
    +------------------------+
    |  Groups : Post/Delete  |
    +------------------------+
    |                        |
    |   check1      Name1    |
    |   check2      Name2    |
    |   check3      Name3    |
    |                        |
    |       SetAsMain        |
    |   Apply   OK   Close   |
    +------------------------+
    """
    def __init__(self, data, win_parent=None):
        self.win_parent = win_parent
        #Init the base class

        groups = data['groups']
        self.names = array([group.name for group in groups])
        self.inames = argsort(self.names)
        #print('inames =', inames)

        for iname, name in enumerate(self.names[self.inames]):
            print('name[%s] = %r' % (iname, name))

        # ignore these...
        self._default_name = data['name']
        self._default_coords = data['coords']
        self._default_elements = data['elements']
        self._default_color = data['color']

        self.coords_pound = data['coords_pound']
        self.elements_pound = data['elements_pound']

        #self._default_is_discrete = data['is_discrete']

        self.out_data = data

        QtGui.QDialog.__init__(self, win_parent)
        #self.setupUi(self)
        self.setWindowTitle('Groups: Post/View')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.show()

    def create_widgets(self):
        # Name
        self.name = QtGui.QLabel("Name:")
        self.name_edit = QtGui.QLineEdit(str(self._default_name).strip())
        self.name_button = QtGui.QPushButton("Default")

        # Name
        self.coords = QtGui.QLabel("Coord IDs:")
        self.coords_edit = QtGui.QLineEdit(str(self._default_coords).strip())
        self.coords_button = QtGui.QPushButton("Default")


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

        # color
        self.color = QtGui.QLabel("Color:")
        #self.color_edit = QtGui.QPushButton()
        #self.color_edit.setColor(self._default_color)
        #self.color_edit = ColorDisplay(self, self._default_color)
        #self.color_button = QtGui.QPushButton("Default")

        # applies a unique implicitly
        self.eids = parse_patran_syntax(str(self._default_elements), pound=self.elements_pound)
        self.cids = parse_patran_syntax(str(self._default_coords), pound=self.coords_pound)

        # continuous / discrete
        #self.checkbox_continuous = QtGui.QCheckBox("Continuous")
        #self.checkbox_discrete = QtGui.QCheckBox("Discrete")
        #self.checkbox_discrete.setChecked(self._default_is_discrete)

        # put these in a group
        #checkboxs2 = QtGui.QButtonGroup(self)
        #checkboxs2.addButton(self.checkbox_continuous)
        #checkboxs2.addButton(self.checkbox_discrete)

        # main/delete/supergroup
        self.set_as_main_button = QtGui.QPushButton("Set As Main")
        self.create_super_group_button = QtGui.QPushButton("Create Super Group")
        self.delete_groups_button = QtGui.QPushButton("Delete Groups")

        # closing
        self.apply_button = QtGui.QPushButton("Apply")
        self.ok_button = QtGui.QPushButton("OK")
        self.cancel_button = QtGui.QPushButton("Cancel")

        self.checks = []
        self.names_text = []
        for iname, name in enumerate(self.names[self.inames]):

            # TODO: create signal when clicked and SetAsMain is pressed
            check = QtGui.QTableWidgetItem()
            check.setCheckState(False)

            # TODO: create right click menu
            name_text = QtGui.QTableWidgetItem(str(name))

            self.checks.append(check)
            self.names_text.append(name_text)

    def create_layout(self):
        if 1:
            nrows = len(self.names)
            grid = QtGui.QTableWidget()
            grid.setRowCount(nrows)
            grid.setColumnCount(2)
            headers = [QtCore.QString('Is Active?'), QtCore.QString('Name')]
            grid.setHorizontalHeaderLabels(headers)
            #grid.resize(400, 250)
            for iname, name in enumerate(self.names[self.inames]):
                #print('name[%s] = %r' % (iname, name))
                check = self.checks[iname]
                name_text = self.names_text[iname]
                # row, col, value
                grid.setItem(iname, 0, check)
                grid.setItem(iname, 1, name_text)
            #grid.horizontalHeaderItem(1).setTextAlignment(QtCore.AlignHCenter)
        else:
            grid = QtGui.QGridLayout()
            for iname, name in enumerate(self.names[self.inames]):
                #print('name[%s] = %r' % (iname, name))
                grid.addWidget(check, iname, 0)
                grid.addWidget(name_text, iname, 1)


        #grid.addWidget(self.set_as_main_button, 6, 1)


        #= QtGui.QVBoxLayout()

        ok_cancel_box = QtGui.QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)


        vbox = QtGui.QVBoxLayout()
        if 1:
            vbox.addWidget(grid)
        else:
            vbox.addLayout(grid)
        vbox.addWidget(self.set_as_main_button)
        vbox.addWidget(self.create_super_group_button)
        vbox.addWidget(self.delete_groups_button)

        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)

        self.setLayout(vbox)

    def set_connections(self):
        #self.connect(self.name_button, QtCore.SIGNAL('clicked()'), self.on_default_name)
        #self.connect(self.coords_button, QtCore.SIGNAL('clicked()'), self.on_default_coords)
        #self.connect(self.elements_button, QtCore.SIGNAL('clicked()'), self.on_default_elements)

        #self.connect(self.add_button, QtCore.SIGNAL('clicked()'), self.on_add)
        #self.connect(self.remove_button, QtCore.SIGNAL('clicked()'), self.on_remove)

        #self.connect(self.color_edit, QtCore.SIGNAL('clicked()'), self.on_edit_color)
        #self.color_edit.clicked.connect(self.on_edit_color)


        #self.connect(self.color_button, QtCore.SIGNAL('clicked()'), self.on_default_color)
        self.connect(self.set_as_main_button, QtCore.SIGNAL('clicked()'), self.on_set_as_main)
        self.connect(self.delete_groups_button, QtCore.SIGNAL('clicked()'), self.on_delete_groups)
        self.connect(self.create_super_group_button, QtCore.SIGNAL('clicked()'), self.on_create_super_group)

        self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)

    def closeEvent(self, event):
        event.accept()

    def on_delete_groups(self):
        pass

    def on_create_super_group(self):
        pass

    def on_set_as_main(self):
        pass

    def on_add(self):
        adict, is_valid = self.check_patran_syntax_dict(self.add_edit)
        if not is_valid:
            #self.add_edit.setStyleSheet("QLineEdit{background: red;}")
            return

        self.eids = _add(adict, ['e', 'elem', 'element'], self.eids)
        self.cids = _add(adict, ['c', 'cid', 'coord'], self.cids)
        self._apply_cids_eids()

        self.add_edit.clear()
        self.add_edit.setStyleSheet("QLineEdit{background: white;}")

    def _apply_cids_eids(self):
        ctext = _get_collapsed_text(self.cids)
        etext = _get_collapsed_text(self.cids)

        self.coords_edit.setText(str(ctext.lstrip()))
        self.elements_edit.setText(str(etext.lstrip()))

    def on_remove(self):
        adict, is_valid = self.check_patran_syntax_dict(self.remove_edit)
        if not is_valid:
            #self.remove_edit.setStyleSheet("QLineEdit{background: red;}")
            return

        self.eids = _remove(adict, ['e', 'elem', 'element'], self.eids)
        self.cids = _remove(adict, ['c', 'cid', 'coord'], self.cids)

        self._apply_cids_eids()
        self.remove_edit.clear()
        self.remove_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_name(self):
        self.name_edit.setText(str(self._default_name))
        self.name_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_coords(self):
        self.coords_edit.setText(str(self._default_coords))
        self.coords_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_elements(self):
        self.elements_edit.setText(str(self._default_elements))
        self.elements_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_validate(self):
        flag0 = flag1 = flag2 = False
        if flag0 and flag1 and flag2:
            self.out_data['name'] = name_value
            self.out_data['elements'] = elements_value
            self.out_data['coords'] = coords_value
            self.out_data['clicked_ok'] = True

            #print("name = %r" % self.name_edit.text())
            #print("min = %r" % self.min_edit.text())
            #print("max = %r" % self.max_edit.text())
            #print("format = %r" % self.format_edit.text())
            return True
        return False

    def on_apply(self):
        passed = self.on_validate()
        if passed:
            self.win_parent.on_post_group(self.out_data)

    def on_ok(self):
        passed = self.on_validate()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.close()


class Group(object):
    def __init__(self, name):
        self.name = name
        self.cids = [0]
        self.eids = [1, 3, 8]

    def __repr__(self):
        msg = 'Group:\n'
        msg += '  cids: [%s]\n' % _get_collapsed_text(self.cids).strip()
        msg += '  eids: [%s]\n' % _get_collapsed_text(self.eids).strip()
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

    nice_blue = (0.1, 0.2, 0.4)

    group1 = Group('this is a really long name')
    group2 = Group('frog')
    group3 = Group('dog')
    print(group3)
    groups = [
        group1, group2, group3,
    ]
    #The Main window
    d = {
        'groups' : groups,
        'name' : 'asdf',
        'coords' : 0,
        'coords_pound' : 4,
        'elements_pound' : 103,
        'elements' : '1:#',
        'color' : nice_blue,
    }
    main_window = GroupsPostView(d)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
