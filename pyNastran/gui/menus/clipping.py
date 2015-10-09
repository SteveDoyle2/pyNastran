from six import string_types
from PyQt4 import QtCore, QtGui
from pyNastran.gui.qt_files.menu_utils import eval_float_from_string

class ClippingPropertiesWindow(QtGui.QDialog):

    def __init__(self, data, win_parent=None):
        self.win_parent = win_parent
        #Init the base class
        self._default_min = data['min']
        self._default_max = data['max']
        self.out_data = {}

        QtGui.QDialog.__init__(self, win_parent)
        #self.setupUi(self)
        self.setWindowTitle('Clipping Properties')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.show()

    def create_widgets(self):
        # Min
        self.min = QtGui.QLabel("Min:")
        self.min_edit = QtGui.QLineEdit(str(self._default_min))
        self.min_button = QtGui.QPushButton("Default")

        # Max
        self.max = QtGui.QLabel("Max:")
        self.max_edit = QtGui.QLineEdit(str(self._default_max))
        self.max_button = QtGui.QPushButton("Default")

        # closing
        self.apply_button = QtGui.QPushButton("Apply")
        self.ok_button = QtGui.QPushButton("OK")
        self.cancel_button = QtGui.QPushButton("Cancel")

    def create_layout(self):
        grid = QtGui.QGridLayout()

        grid.addWidget(self.min, 0, 0)
        grid.addWidget(self.min_edit, 0, 1)
        grid.addWidget(self.min_button, 0, 2)

        grid.addWidget(self.max, 1, 0)
        grid.addWidget(self.max_edit, 1, 1)
        grid.addWidget(self.max_button, 1, 2)

        ok_cancel_box = QtGui.QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        vbox = QtGui.QVBoxLayout()
        vbox.addLayout(grid)

        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        self.connect(self.min_button, QtCore.SIGNAL('clicked()'), self.on_default_min)
        self.connect(self.max_button, QtCore.SIGNAL('clicked()'), self.on_default_max)

        self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)
        self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.close()

    def closeEvent(self, event):
        self.out_data['close'] = True
        event.accept()

    def on_default_min(self):
        self.min_edit.setText(str(self._default_min))
        self.min_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_max(self):
        self.max_edit.setText(str(self._default_max))
        self.max_edit.setStyleSheet("QLineEdit{background: white;}")

    def check_float(self, cell):
        text = cell.text()
        try:
            value = eval_float_from_string(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
            return value, True
        except ValueError:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    def on_validate(self):
        min_value, flag0 = self.check_float(self.min_edit)
        max_value, flag1 = self.check_float(self.max_edit)

        if flag0 and flag1:
            self.out_data['min'] = min(min_value, max_value)
            self.out_data['max'] = max(min_value, max_value)
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_apply(self):
        passed = self.on_validate()
        if passed:
            self.win_parent._apply_clipping(self.out_data)
        return passed

    def on_ok(self):
        passed = self.on_apply()
        if passed:
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
    d = {
        'min' : 0.,
        'max' : 10,
    }
    main_window = ClippingPropertiesWindow(d)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
