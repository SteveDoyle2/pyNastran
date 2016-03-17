from copy import deepcopy
from PyQt4 import QtCore, QtGui


class CameraWindow(QtGui.QDialog):
    def __init__(self, data, win_parent=None):
        """
        +--------+
        | Camera |
        +--------+---------------+
        |  Camera Name           |
        |  +-------------------+ |
        |  |                   | |
        |  |                   | |
        |  |                   | |
        |  |                   | |
        |  |                   | |
        |  +-------------------+ |
        |                        |
        | Name xxx       Save    |
        | Delete   Set           |
        |                        |
        |    Apply   OK  Cancel  |
        +--------+---------------+
        """
        QtGui.QDialog.__init__(self, win_parent)
        self.win_parent = win_parent
        self.setWindowTitle('Camera Views')
        #self.setWindowIcon(view_icon)

        self._default_name = 'Camera'
        self.out_data = data
        self.out_data['clicked_ok'] = False

        self.cameras = deepcopy(data['cameras'])
        self.names = sorted(self.cameras.keys())

        self.name = QtGui.QLabel("Name:")
        self.name_edit = QtGui.QLineEdit(str(self._default_name))

        self.delete_button = QtGui.QPushButton("Delete")
        self.set_button = QtGui.QPushButton("Set")
        self.save_button = QtGui.QPushButton("Save")

        # closing
        self.apply_button = QtGui.QPushButton("Apply")
        #self.ok_button = QtGui.QPushButton("OK")
        self.close_button = QtGui.QPushButton("Close")
        self.cancel_button = QtGui.QPushButton("Cancel")

        self.table = QtGui.QTableWidget()
        names_text = []
        for iname, name in enumerate(self.names):
            name_text = QtGui.QTableWidgetItem(str(name))
            names_text.append(name_text)
        self.create_layout(names_text)
        self.set_connections()

    def create_layout(self, names_text):
        nrows = len(self.names)
        table = self.table
        table.setRowCount(nrows)
        table.setColumnCount(1)
        headers = [QtCore.QString('Camera Name')]
        table.setHorizontalHeaderLabels(headers)

        header = table.horizontalHeader()
        header.setStretchLastSection(True)

        for iname, name_text in enumerate(names_text):
            # row, col, value
            table.setItem(iname, 0, name_text)
        table.resizeRowsToContents()

        ok_cancel_box = QtGui.QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        #ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.close_button)
        ok_cancel_box.addWidget(self.cancel_button)

        grid = QtGui.QGridLayout()

        irow = 0
        grid.addWidget(self.name, irow, 0)
        grid.addWidget(self.name_edit, irow, 1)
        grid.addWidget(self.save_button, irow, 2)
        irow += 1

        grid.addWidget(self.delete_button, irow, 0)
        grid.addWidget(self.set_button, irow, 1)
        irow += 1


        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.table)

        vbox.addLayout(grid)
        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        self.connect(self.set_button, QtCore.SIGNAL('clicked()'), self.on_set)
        self.connect(self.save_button, QtCore.SIGNAL('clicked()'), self.on_save)
        self.connect(self.delete_button, QtCore.SIGNAL('clicked()'), self.on_delete)

        self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        #self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        self.connect(self.close_button, QtCore.SIGNAL('clicked()'), self.on_close)
        self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)

    def on_set(self):
        objs = self.table.selectedIndexes()
        if len(objs) == 1:
            obj = objs[0]
            irow = obj.row()
            name = self.names[irow]
            #print('name =', name)
            self.set_camera(name)
            return True
        return False

    def on_save(self):
        name = str(self.name_edit.text()).strip()
        if name in self.cameras:
            return
        irow = self.nrows
        if len(name):
            self.table.insertRow(irow)
            name_text = QtGui.QTableWidgetItem(str(name))
            self.table.setItem(irow, 0, name_text)
            self.name_edit.setText('')
            self.save_camera(name)

    def set_camera(self, name):
        camera_data = self.cameras[name]
        if self.win_parent is None:
            return
        self.win_parent.on_set_camera_data(camera_data)

    def save_camera(self, name):
        self.names.append(name)
        if self.win_parent is None:
            self.cameras[name] = None
            return

        self.cameras[name] = self.win_parent.get_camera_data()

    #@property
    #def camera(self):

    @property
    def nrows(self):
        return self.table.rowCount()

    def on_delete(self):
        irows = []
        for obj in self.table.selectedIndexes():
            irow = obj.row()
            irows.append(irow)
        irows.sort()

        for irow in reversed(irows):
            self.table.removeRow(irow)
            #print('delete', self.names)
            name = self.names.pop(irow)
            del self.cameras[name]
            #print('  removing irow=%s name=%r' % (irow, name))

    def closeEvent(self, event):
        event.accept()

    def check_name(self, cell):
        text = str(cell.text()).strip()
        if len(text):
            cell.setStyleSheet("QLineEdit{background: white;}")
            return text, True
        else:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    #def on_validate(self):
        #name_value, flag0 = self.check_name(self.name_edit)
        #if flag0:
            #self.out_data['cameras'] = self.cameras
            #self.out_data['clicked_ok'] = True
            #return True
        #return False

    def on_apply(self):
        passed = self.on_set()
        #if passed:
        #    self.win_parent.create_plane(self.out_data)
        return passed

    def on_close(self):
        self.out_data['clicked_ok'] = True
        self.out_data['cameras'] = self.cameras
        self.close()


    def on_ok(self):
        passed = self.on_apply()
        if passed:
            name = str(self.name_edit.text()).strip()
            self.out_data['name'] = name
            self.out_data['cameras'] = self.cameras
            self.out_data['clicked_ok'] = True
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.close()


def main():
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)


    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QtGui.QApplication(sys.argv)
    #The Main window
    a = [
        [1., 1., 1.],
        [1., 1., 1.],
        [1., 1., 1.],
    ]
    b = [
        [1., 1., 1.],
        [1., 1., 1.],
        [1., 1., 1.],
    ]
    c = [
        [1., 1., 1.],
        [1., 1., 1.],
        [1., 1., 1.],
    ]

    d = {
        'cameras' : {'a':a, 'b':b, 'c':c},
    }
    main_window = CameraWindow(d, win_parent=None)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
