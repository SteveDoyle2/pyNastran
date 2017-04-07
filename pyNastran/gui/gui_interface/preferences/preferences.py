from __future__ import print_function
from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    #from PyQt4 import QtCore, QtGui
    from PyQt4 import QtGui
    from PyQt4.QtGui import (
        QDialog, QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
        QSpinBox)
elif qt_version == 5:
    #from PyQt5 import QtCore, QtGui
    from PyQt5 import QtGui
    from PyQt5.QtWidgets import (
        QDialog, QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
        QSpinBox)
elif qt_version == 'pyside':
    #from PySide import QtCore, QtGui
    from PySide import QtGui
    from PySide.QtGui import (
        QDialog, QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
        QSpinBox)
else:
    raise NotImplementedError('qt_version = %r' % qt_version)

from pyNastran.gui.gui_interface.common import PyDialog

class PreferencesWindow(PyDialog):
    """
    +-------------+
    | Preferences |
    +------------------------------+
    | Text Size     ______ Default |
    |                              |
    |       Apply OK Cancel        |
    +------------------------------+
    """

    def __init__(self, data, win_parent=None):
        #Init the base class
        PyDialog.__init__(self, data, win_parent)

        self._updated_preference = False

        self._default_textsize = data['font_size']

        #self.setupUi(self)
        self.setWindowTitle('Preferences')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        self.on_set_font(self._default_textsize)
        #self.show()

    def create_widgets(self):
        # Text Size
        self.textsize = QLabel("Text Size:")
        self.textsize_edit = QSpinBox(self)
        self.textsize_edit.setValue(self._default_textsize)
        self.textsize_edit.setRange(7, 20)
        self.textsize_button = QPushButton("Default")

        # closing
        self.apply_button = QPushButton("Apply")
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")

    def create_layout(self):
        grid = QGridLayout()

        grid.addWidget(self.textsize, 0, 0)
        grid.addWidget(self.textsize_edit, 0, 1)
        grid.addWidget(self.textsize_button, 0, 2)

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
        #if qt_version == 4:
            #self.connect(self.textsize_button, QtCore.SIGNAL('clicked()'), self.on_default_textsize)
            #self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
            #self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
            #self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)
            #self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)
        #else:
        self.textsize_button.clicked.connect(self.on_default_textsize)
        self.textsize_edit.valueChanged.connect(self.on_set_font)

        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)
        # closeEvent

    def on_set_font(self, value=None):
        if value is None:
            value = self.textsize_edit.value()
        font = QtGui.QFont()
        font.setPointSize(value)
        self.setFont(font)

    def on_default_textsize(self):
        self.textsize_edit.setValue(self._default_textsize)
        self.on_set_font(self._default_textsize)

    #@staticmethod
    #def check_float(cell):
        #text = cell.text()
        #try:
            #value = eval_float_from_string(text)
            #cell.setStyleSheet("QLineEdit{background: white;}")
            #return value, True
        #except ValueError:
            #cell.setStyleSheet("QLineEdit{background: red;}")
            #return None, False

    def on_validate(self):
        textsize_value, flag0 = self.check_float(self.textsize_edit)

        if flag0:
            self.out_data['font_size'] = int(textsize_value)
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_apply(self):
        passed = self.on_validate()
        if passed:
            self.win_parent._apply_preferences(self.out_data)
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
        'text_size' : 14,
    }
    main_window = PreferencesWindow(d)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
