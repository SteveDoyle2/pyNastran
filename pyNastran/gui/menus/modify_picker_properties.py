"""
based on manage_actors.py
"""
from __future__ import print_function
from math import log10, ceil, floor
from six import iteritems
from PyQt4 import QtCore, QtGui


class ModifyPickerPropertiesMenu(QtGui.QDialog):

    def __init__(self, data, win_parent=None):
        self.win_parent = win_parent

        self._size = data['size']
        self.out_data = data
        self.dim_max = data['dim_max']

        QtGui.QDialog.__init__(self, win_parent)
        self.setWindowTitle('Modify Picker Properties')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        width = 260
        height = 130
        self.resize(width, height)
        #self.show()

    def create_widgets(self):
        # Size
        self.size = QtGui.QLabel("Size:")
        self.size_edit = QtGui.QDoubleSpinBox(self)
        self.size_edit.setRange(0.0, self.dim_max)

        log_dim = log10(self.dim_max)
        decimals = int(ceil(abs(log_dim)))
        decimals = max(3, decimals)
        self.size_edit.setDecimals(decimals)
        self.size_edit.setSingleStep(self.dim_max / 100.)
        self.size_edit.setValue(self._size)

        # closing
        #self.apply_button = QtGui.QPushButton("Apply")
        #self.ok_button = QtGui.QPushButton("OK")
        self.cancel_button = QtGui.QPushButton("Close")

    def create_layout(self):
        grid = QtGui.QGridLayout()

        grid.addWidget(self.size, 1, 0)
        grid.addWidget(self.size_edit, 1, 1)

        ok_cancel_box = QtGui.QHBoxLayout()
        #ok_cancel_box.addWidget(self.apply_button)
        #ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        vbox = QtGui.QVBoxLayout()
        vbox.addLayout(grid)

        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)
        self.layout()

    #def on_color(self):
        #pass

    def set_connections(self):
        self.size_edit.valueChanged.connect(self.on_size)
        self.connect(self.size_edit, QtCore.SIGNAL('editingFinished()'), self.on_size)
        self.connect(self.size_edit, QtCore.SIGNAL('valueChanged()'), self.on_size)
        self.connect(self.size_edit, QtCore.SIGNAL('clicked()'), self.on_size)

        #self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        #self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)
        self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.close()

    def closeEvent(self, event):
        self.out_data['close'] = True
        event.accept()

    def on_size(self):
        self._size = float(self.size_edit.text())
        self.on_apply(force=True)

    def check_float(self, cell):
        text = cell.text()
        value = float(text)
        return value, True

    def on_validate(self):
        size_value, flag0 = self.check_float(self.size_edit)

        if flag0:
            self._size = size_value
            #self.out_data['min'] = min(min_value, max_value)
            #self.out_data['max'] = max(min_value, max_value)
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_apply(self, force=False):
        passed = self.on_validate()
        if (passed or Force) and self.win_parent is not None:
            self.win_parent.set_element_picker_size(self._size)
        return passed

    def on_ok(self):
        self.out_data['clicked_ok'] = True
        self.out_data['clicked_cancel'] = False
        self.out_data['close'] = True
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.out_data['clicked_cancel'] = True
        self.out_data['close'] = True
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
        'size' : 10.,
        'dim_max' : 502.
    }
    main_window = ModifyPickerPropertiesMenu(d)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
