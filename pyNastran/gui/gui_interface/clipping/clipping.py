from six import string_types

from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4 import QtCore
    from PyQt4.QtGui import (
        QLabel, QLineEdit, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout)
elif qt_version == 5:
    from PyQt5 import QtCore
    from PyQt5.QtWidgets import (
        QLabel, QLineEdit, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout)
elif qt_version == 'pyside':
    from PySide import QtCore
    from PySide.QtGui import (
        QLabel, QLineEdit, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout)
else:
    raise NotImplementedError('qt_version = %r' % qt_version)

from pyNastran.gui.gui_interface.common import PyDialog
from pyNastran.gui.qt_files.menu_utils import eval_float_from_string

class ClippingPropertiesWindow(PyDialog):
    """
    +---------------------+
    | Clipping Properties |
    +------------------------------+
    | Clipping Min  ______ Default |
    | Clipping Max  ______ Default |
    |                              |
    |       Apply OK Cancel        |
    +------------------------------+
    """

    def __init__(self, data, win_parent=None):
        #Init the base class
        PyDialog.__init__(self, data, win_parent)
        self.set_font_size(data['font_size'])

        self._updated_clipping = False

        self._default_min = data['clipping_min']
        self._default_max = data['clipping_max']

        #self.setupUi(self)
        self.setWindowTitle('Clipping Properties')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.show()

    def create_widgets(self):
        # Min
        self.min = QLabel("Min:")
        self.min_edit = QLineEdit(str(self._default_min))
        self.min_button = QPushButton("Default")

        # Max
        self.max = QLabel("Max:")
        self.max_edit = QLineEdit(str(self._default_max))
        self.max_button = QPushButton("Default")

        # closing
        self.apply_button = QPushButton("Apply")
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")

    def create_layout(self):
        grid = QGridLayout()

        grid.addWidget(self.min, 0, 0)
        grid.addWidget(self.min_edit, 0, 1)
        grid.addWidget(self.min_button, 0, 2)

        grid.addWidget(self.max, 1, 0)
        grid.addWidget(self.max_edit, 1, 1)
        grid.addWidget(self.max_button, 1, 2)

        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)

        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        if qt_version == 4:
            self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)
        #else:
            # closeEvent
        self.min_button.clicked.connect(self.on_default_min)
        self.max_button.clicked.connect(self.on_default_max)
        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)

    def on_default_min(self):
        self.min_edit.setText(str(self._default_min))
        self.min_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_max(self):
        self.max_edit.setText(str(self._default_max))
        self.max_edit.setStyleSheet("QLineEdit{background: white;}")

    @staticmethod
    def check_float(cell):
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
            self.out_data['clipping_min'] = min(min_value, max_value)
            self.out_data['clipping_max'] = max(min_value, max_value)
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
    app = QApplication(sys.argv)
    #The Main window
    d = {
        'clipping_min' : 0.,
        'clipping_max' : 10,
    }
    main_window = ClippingPropertiesWindow(d)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
