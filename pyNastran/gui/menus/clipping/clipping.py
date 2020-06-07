"""
defines:
 - ClippingPropertiesWindow
"""
from pyNastran.gui.qt_version import qt_int as qt_version
from qtpy import QtCore
from qtpy.QtWidgets import (
    QLabel, QLineEdit, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout)

from pyNastran.gui.utils.qt.pydialog import PyDialog, check_float
from pyNastran.gui.menus.menu_utils import eval_float_from_string

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
        """creates ClippingPropertiesWindow"""
        #Init the base class
        PyDialog.__init__(self, data, win_parent)
        self.set_font_size(data['font_size'])

        self._updated_clipping = False

        self._default_min = data['min_clip']
        self._default_max = data['max_clip']

        #self.setupUi(self)
        self.setWindowTitle('Clipping Properties')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.show()

    def create_widgets(self):
        # Min
        self.min_label = QLabel("Min:")
        self.min_edit = QLineEdit(str(self._default_min))
        self.min_button = QPushButton("Default")

        # Max
        self.max_label = QLabel("Max:")
        self.max_edit = QLineEdit(str(self._default_max))
        self.max_button = QPushButton("Default")

        # closing
        self.ok_button = QPushButton("OK")

    def create_layout(self):
        grid = QGridLayout()

        grid.addWidget(self.min_label, 0, 0)
        grid.addWidget(self.min_edit, 0, 1)
        grid.addWidget(self.min_button, 0, 2)

        grid.addWidget(self.max_label, 1, 0)
        grid.addWidget(self.max_edit, 1, 1)
        grid.addWidget(self.max_button, 1, 2)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)

        vbox.addStretch()
        vbox.addLayout(self.ok_button)
        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        self.min_button.clicked.connect(self.on_default_min)
        self.max_button.clicked.connect(self.on_default_max)
        self.min_edit.textChanged.connect(self.on_apply)
        self.max_edit.textChanged.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)

    def on_default_min(self):
        self.min_edit.setText(str(self._default_min))
        self.min_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_max(self):
        self.max_edit.setText(str(self._default_max))
        self.max_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_validate(self):
        min_value, flag0 = check_float(self.min_edit)
        max_value, flag1 = check_float(self.max_edit)

        if flag0 and flag1:
            self.out_data['min_clip'] = min(min_value, max_value)
            self.out_data['max_clip'] = max(min_value, max_value)
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_apply(self):
        passed = self.on_validate()
        if passed:
            self.win_parent.clipping_obj.apply_clipping(self.out_data)
        return passed

    def on_ok(self):
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.out_data['close'] = True
        self.close()


def main():  # pragma: no cover
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
        'min_clip' : 0.,
        'max_clip' : 10,
    }
    main_window = ClippingPropertiesWindow(d)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == '__main__':   # pragma: no cover
    main()
