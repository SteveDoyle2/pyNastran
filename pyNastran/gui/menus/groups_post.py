# -*- coding: utf-8 -*-
from numpy import array, arange

from pyNastran.gui.qt_version import qt_version
from qtpy import QtCore, QtGui
from qtpy.QtWidgets import (
    QDialog, QPushButton, QApplication,
    QHBoxLayout, QVBoxLayout, QTableWidget, QTableWidgetItem,
)

if qt_version == 5:
    QString =  str
else:
    raise NotImplementedError('qt_version = %r' % qt_version)

from pyNastran.gui.menus.groups_modify import Group


class GroupsPostView(QDialog):
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
        inames = data['inames']
        self.imain = data['imain']
        self.names = [group.name for group in groups]
        self.white = (255, 255, 255)
        self.light_grey = (211, 211, 211)
        self.inames = inames
        self.shown_set = data['shown']
        self.deleted_groups = set()
        #self.inames = argsort(self.names)
        #print('inames =', inames)

        anames = array(self.names)
        for iname, name in enumerate(anames[self.inames]):
            print('name[%s] = %r' % (iname, name))

        # ignore these...
        #self._default_name = data['name']
        #self._default_coords = data['coords']
        #self._default_elements = data['elements']
        #self._default_color = data['color']

        #self.coords_pound = data['coords_pound']
        #self.elements_pound = data['elements_pound']

        #self._default_is_discrete = data['is_discrete']

        self.out_data = data

        QDialog.__init__(self, win_parent)
        #self.setupUi(self)
        self.setWindowTitle('Groups: Post/View')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.show()

    def create_widgets(self):
        # main/delete/supergroup
        self.set_as_main_button = QPushButton("Set As Main")
        self.create_super_group_button = QPushButton("Create Super Group")
        self.delete_groups_button = QPushButton("Delete Groups")
        self.revert_groups_button = QPushButton("Revert Groups")

        self.show_groups_button = QPushButton("Show Groups")
        self.hide_groups_button = QPushButton("Hide Groups")

        # closing
        self.apply_button = QPushButton("Apply")
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")

        #table
        self.table = QTableWidget()
        self.checks = []
        self.names_text = []

        bold = QtGui.QFont()
        bold.setBold(True)
        bold.setItalic(True)
        bold.setWeight(75)
        anames = array(self.names)
        for iname, name in enumerate(anames[self.inames]):
            check = QTableWidgetItem()
            check.setCheckState(False)

            # TODO: create right click menu ???
            name_text = QTableWidgetItem(str(name))
            if iname == self.imain:
                name_text.setFont(bold)
                self.shown_set.add(iname)
                check.setCheckState(2)
                name_text.setBackground(QtGui.QColor(*self.light_grey))
            elif iname in self.shown_set:
                name_text.setBackground(QtGui.QColor(*self.light_grey))

            self.checks.append(check)
            self.names_text.append(name_text)

    def create_layout(self):
        nrows = len(self.names)
        table = self.table
        table.setRowCount(nrows)
        table.setColumnCount(2)
        headers = [QString('Operate On'), QString('Name')]
        table.setHorizontalHeaderLabels(headers)

        header = table.horizontalHeader()
        header.setStretchLastSection(True)

        #table.setAlternatingRowColors(True)

        #header = table.verticalHeader()
        #header.setStretchLastSection(True)
        #table.resize(400, 250)

        #heighti = table.rowHeight(0)
        #total_height = nrows * heighti
        #table.setMaximumHeight(total_height)
        #table.resize(total_height, None)

        #for iname, name in enumerate(self.names[self.inames]):
        #print('name[%s] = %r' % (iname, name))
        for iname in self.inames:
            check = self.checks[iname]
            name_text = self.names_text[iname]
            # row, col, value
            table.setItem(iname, 0, check)
            table.setItem(iname, 1, name_text)
        table.resizeRowsToContents()
        #table.horizontalHeaderItem(1).setTextAlignment(QtCore.AlignHCenter)

        #= QVBoxLayout()
        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        vbox = QVBoxLayout()
        vbox.addWidget(table)
        vbox.addWidget(self.set_as_main_button)
        #vbox.addWidget(self.create_super_group_button)

        vbox.addStretch()
        vbox.addWidget(self.show_groups_button)
        vbox.addWidget(self.hide_groups_button)
        vbox.addStretch()
        vbox.addWidget(self.delete_groups_button)
        vbox.addWidget(self.revert_groups_button)
        vbox.addStretch()

        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)

        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        self.set_as_main_button.clicked.connect(self.on_set_as_main)
        self.delete_groups_button.clicked.connect(self.on_delete_groups)
        self.revert_groups_button.clicked.connect(self.on_revert_groups)

        self.show_groups_button.clicked.connect(self.on_show_groups)
        self.hide_groups_button.clicked.connect(self.on_hide_groups)

        self.create_super_group_button.clicked.connect(self.on_create_super_group)

        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)

    def closeEvent(self, event):
        event.accept()

    @property
    def nrows(self):
        return self.table.rowCount()

    def on_hide_groups(self):
        self._set_highlight(self.white)

    def on_show_groups(self):
        self._set_highlight(self.light_grey)

    def _set_highlight(self, color):
        for irow in range(self.nrows):
            check = self.checks[irow]
            is_checked = check.checkState()

            # 0 - unchecked
            # 1 - partially checked (invalid)
            # 2 - checked
            if is_checked:
                name_text = self.names_text[irow]
                name_text.setBackground(QtGui.QColor(*color))

    def on_delete_groups(self):
        for irow in range(self.nrows):
            check = self.checks[irow]
            is_checked = check.checkState()

            # 0 - unchecked
            # 1 - partially checked (invalid)
            # 2 - checked
            if irow == 0 and is_checked:
                # TODO: change this to a log
                print('error deleting group ALL...change this to a log')
                #self.window_parent.log
                return
            if is_checked:
                self.table.hideRow(irow)
                self.deleted_groups.add(irow)
                check.setCheckState(0)

        if self.imain > 0 and self.shown_set == set([0]):
            bold = QtGui.QFont()
            bold.setBold(True)
            bold.setItalic(True)

            self.imain = 0
            irow = 0
            check = self.checks[irow]
            name_text = self.names_texts[irow]
            name_text.setFont(bold)
            name_text.setBackground(QtGui.QColor(*self.light_grey))

    def on_revert_groups(self):
        for irow in range(self.nrows):
            self.table.showRow(irow)
        self.deleted_groups = set()

    def on_create_super_group(self):
        inames = [iname for iname, check in enumerate(self.checks)
                  if bool(check.checkState())]

        if not len(inames):
            # TODO: add logging
            print('nothing is checked...')
            return
        if inames[0] == 0:
            # TODO: add logging
            print("cannot include 'ALL' in supergroup...")
            return

        name = 'SuperGroup'
        # popup gui and get a name

        irow = self.table.rowCount()
        self.table.insertRow(irow)


        check = QTableWidgetItem()
        check.setCheckState(False)
        name_text = QTableWidgetItem(str(name))

        self.names.extend(name)
        self.names_text.append(name_text)
        self.checks.append(check)

        self.table.setItem(irow, 0, check)
        self.table.setItem(irow, 1, name_text)


    def on_set_as_main(self):
        bold = QtGui.QFont()
        bold.setBold(True)
        bold.setItalic(True)

        normal = QtGui.QFont()
        normal.setBold(False)
        normal.setItalic(False)

        imain = None
        imain_set = False
        for irow in range(self.nrows):
            check = self.checks[irow]
            name_text = self.names_text[irow]
            is_checked = check.checkState()

            # 0 - unchecked
            # 1 - partially checked (invalid)
            # 2 - checked
            if is_checked and not imain_set:
                # TODO: change this to a log
                #self.window_parent.log
                imain_set = True
                imain = irow
                name_text.setFont(bold)
                name_text.setBackground(QtGui.QColor(*self.light_grey))
                self.shown_set.add(irow)
            elif irow == self.imain:
                name_text.setFont(normal)
                if irow == 0:
                    name_text.setBackground(QtGui.QColor(*self.white))
                    if irow in self.shown_set:
                        self.shown_set.remove(irow)
                elif imain == 0:
                    name_text.setBackground(QtGui.QColor(*self.white))
                    self.shown_set.remove(imain)
        self.imain = imain

    def get_main_group(self):
        return self.imain

    def get_shown_group(self):
        return self.shown_set

    def get_deleted_groups(self):
        return self.deleted_groups

    def on_validate(self):
        flag0 = flag1 = flag2 = True
        main_group_id = self.get_main_group()
        shown_groups_ids = self.get_shown_group()
        deleted_group_ids = self.get_deleted_groups()

        if flag0 and flag1 and flag2:
            self.out_data['imain'] = main_group_id
            self.out_data['shown'] = shown_groups_ids
            self.out_data['remove'] = deleted_group_ids
            self.out_data['clicked_ok'] = True
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

def on_post_group(data):
    print('hi')

def main():  # pragma: no cover
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    app.on_post_group = on_post_group

    group1 = Group('this is a really long name', [1, 2, 3], 4)
    group2 = Group('frog', [1, 3], 4)
    group3 = Group('dog', [1, 2, 3, 5], 4)
    all_group = Group('ALL', [1, 2, 3, 34], 4)
    print(group3)
    groups = [
        all_group, group1, group2, group3,
        all_group, group1, group2, group3,
        all_group, group1, group2, group3,
    ]
    #The Main window
    data = {
        'groups' : groups,
        'inames' : arange(len(groups)),
        'imain' : 0,
        'shown' : set([1, 2, 3]),
        'remove' : None,
    }
    win_parent = None
    main_window = GroupsPostView(data, win_parent=win_parent)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == '__main__':  # pragma: no cover
    main()
