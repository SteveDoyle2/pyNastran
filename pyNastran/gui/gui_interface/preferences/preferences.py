"""
The preferences menu handles:
 - Font Size
 - Background Color
 - Text Color
 - Label Color
 - Label Size
 - Clipping Min
 - Clipping Max
"""
from __future__ import print_function
from math import log10, ceil
from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    #from PyQt4 import QtCore, QtGui
    from PyQt4 import QtGui
    from PyQt4.QtGui import (
        QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
        QSpinBox, QDoubleSpinBox, QColorDialog, QLineEdit)
elif qt_version == 5:
    #from PyQt5 import QtCore, QtGui
    from PyQt5 import QtGui
    from PyQt5.QtWidgets import (
        QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
        QSpinBox, QDoubleSpinBox, QColorDialog, QLineEdit)
elif qt_version == 'pyside':
    #from PySide import QtCore, QtGui
    from PySide import QtGui
    from PySide.QtGui import (
        QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
        QSpinBox, QDoubleSpinBox, QColorDialog, QLineEdit)
else:
    raise NotImplementedError('qt_version = %r' % qt_version)

from pyNastran.gui.gui_interface.common import PyDialog, QPushButtonColor
from pyNastran.gui.qt_files.menu_utils import eval_float_from_string


class PreferencesWindow(PyDialog):
    """
    +-------------+
    | Preferences |
    +------------------------------+
    | Text Size     ______ Default |
    | Label Color   ______         |
    | Label Size    ______         |
    | Picker Size   ______         |
    | Back Color    ______         |
    | Text Color    ______         |
    |                              |
    |       Apply OK Cancel        |
    +------------------------------+
    """
    def __init__(self, data, win_parent=None):
        """
        Saves the data members from data and
        performs type checks
        """
        PyDialog.__init__(self, data, win_parent)

        self._updated_preference = False

        self._default_font_size = data['font_size']
        self._default_clipping_min = data['clipping_min']
        self._default_clipping_max = data['clipping_max']

        #label_color_float = data['label_color']
        self._label_size = data['label_size']
        #self.out_data = data
        self.dim_max = data['dim_max']
        self._picker_size = data['picker_size'] * 100.

        self.label_color_float, self.label_color_int = _check_color(data['label_color'])
        self.background_color_float, self.background_color_int = _check_color(
            data['background_color'])
        self.text_color_float, self.text_color_int = _check_color(data['text_color'])

        #self.setupUi(self)
        self.setWindowTitle('Preferences')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        self.on_set_font(self._default_font_size)
        #self.show()

    def create_widgets(self):
        """creates the display window"""
        # Text Size
        self.font_size = QLabel("Text Size:")
        self.font_size_edit = QSpinBox(self)
        self.font_size_edit.setValue(self._default_font_size)
        self.font_size_edit.setRange(7, 20)
        self.font_size_button = QPushButton("Default")

        #-----------------------------------------------------------------------
        # Annotation Color
        self.label_color = QLabel("Label Color:")
        self.label_color_edit = QPushButtonColor(self.label_color_int)

        # Background Color
        self.background_color = QLabel("Background Color:")
        self.background_color_edit = QPushButtonColor(self.background_color_int)

        # Label Color
        self.text_color = QLabel("Text Color:")
        self.text_color_edit = QPushButtonColor(self.text_color_int)

        #-----------------------------------------------------------------------
        # Label Size
        self.label_size = QLabel("Label Size (3D Text):")
        self.label_size_edit = QDoubleSpinBox(self)
        self.label_size_edit.setRange(0.0, self.dim_max)

        log_dim = log10(self.dim_max)
        decimals = int(ceil(abs(log_dim)))
        decimals = max(6, decimals)
        self.label_size_edit.setDecimals(decimals)
        #self.label_size_edit.setSingleStep(self.dim_max / 100.)
        self.label_size_edit.setSingleStep(self.dim_max / 1000.)
        self.label_size_edit.setValue(self._label_size)

        #-----------------------------------------------------------------------
        # Picker Size
        self.picker_size = QLabel("Picker Size (% of Screen):")
        self.picker_size_edit = QDoubleSpinBox(self)
        self.picker_size_edit.setRange(0., 10.)

        log_dim = log10(self.dim_max)
        decimals = int(ceil(abs(log_dim)))

        decimals = max(6, decimals)
        self.picker_size_edit.setDecimals(decimals)
        self.picker_size_edit.setSingleStep(10. / 5000.)
        self.picker_size_edit.setValue(self._picker_size)

        #-----------------------------------------------------------------------
        # Clipping Min
        self.clipping_min = QLabel("Clipping Min:")
        self.clipping_min_edit = QLineEdit(str(self._default_clipping_min))
        self.clipping_min_button = QPushButton("Default")

        # Clipping Max
        self.clipping_max = QLabel("Clipping Max:")
        self.clipping_max_edit = QLineEdit(str(self._default_clipping_max))
        self.clipping_max_button = QPushButton("Default")

        #-----------------------------------------------------------------------
        # closing
        self.apply_button = QPushButton("Apply")
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")

    def create_layout(self):
        grid = QGridLayout()

        grid.addWidget(self.font_size, 0, 0)
        grid.addWidget(self.font_size_edit, 0, 1)
        grid.addWidget(self.font_size_button, 0, 2)

        grid.addWidget(self.background_color, 1, 0)
        grid.addWidget(self.background_color_edit, 1, 1)

        grid.addWidget(self.text_color, 2, 0)
        grid.addWidget(self.text_color_edit, 2, 1)

        grid.addWidget(self.label_color, 3, 0)
        grid.addWidget(self.label_color_edit, 3, 1)

        grid.addWidget(self.label_size, 4, 0)
        grid.addWidget(self.label_size_edit, 4, 1)

        grid.addWidget(self.picker_size, 5, 0)
        grid.addWidget(self.picker_size_edit, 5, 1)

        #grid.addWidget(self.clipping_min, 6, 0)
        #grid.addWidget(self.clipping_min_edit, 6, 1)
        #grid.addWidget(self.clipping_min_button, 6, 2)

        #grid.addWidget(self.clipping_max, 7, 0)
        #grid.addWidget(self.clipping_max_edit, 7, 1)
        #grid.addWidget(self.clipping_max_button, 7, 2)

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
        self.font_size_button.clicked.connect(self.on_default_font_size)
        self.font_size_edit.valueChanged.connect(self.on_set_font)

        self.label_size_edit.editingFinished.connect(self.on_label_size)
        self.label_size_edit.valueChanged.connect(self.on_label_size)
        self.label_color_edit.clicked.connect(self.on_label_color)

        self.background_color_edit.clicked.connect(self.on_background_color)
        self.text_color_edit.clicked.connect(self.on_text_color)

        self.picker_size_edit.valueChanged.connect(self.on_picker_size)
        self.picker_size_edit.editingFinished.connect(self.on_picker_size)
        self.picker_size_edit.valueChanged.connect(self.on_picker_size)

        self.clipping_min_button.clicked.connect(self.on_default_clipping_min)
        self.clipping_max_button.clicked.connect(self.on_default_clipping_max)

        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)
        # closeEvent

    def on_set_font(self, value=None):
        """update the font for the current window"""
        if value is None:
            value = self.font_size_edit.value()
        font = QtGui.QFont()
        font.setPointSize(value)
        self.setFont(font)

    def update_label_size_color(self):
        if self.win_parent is not None:
            self.win_parent.set_labelsize_color(self._label_size, self.label_color_float)

    def on_label_color(self):
        rgb_color_ints = self.label_color_int
        title = "Choose a label color"
        passed, rgb_color_ints, rgb_color_floats = self.on_color(
            self.label_color_edit, rgb_color_ints, title)
        if passed:
            self.label_color_int = rgb_color_ints
            self.label_color_float = rgb_color_floats
            self.update_label_size_color()

    def on_background_color(self):
        """ Choose a background color """
        rgb_color_ints = self.background_color_int
        title = "Choose a background color"
        passed, rgb_color_ints, rgb_color_floats = self.on_color(
            self.background_color_edit, rgb_color_ints, title)
        if passed:
            self.background_color_int = rgb_color_ints
            self.background_color_float = rgb_color_floats
            if self.win_parent is not None:
                self.win_parent.set_background_color(rgb_color_floats)

    def on_text_color(self):
        """ Choose a text color """
        rgb_color_ints = self.text_color_int
        title = "Choose a text color"
        passed, rgb_color_ints, rgb_color_floats = self.on_color(
            self.text_color_edit, rgb_color_ints, title)
        if passed:
            self.text_color_int = rgb_color_ints
            self.text_color_float = rgb_color_floats
            if self.win_parent is not None:
                self.win_parent.set_text_color(rgb_color_floats)

    def on_color(self, color_edit, rgb_color_ints, title):
        """pops a color dialog"""
        col = QColorDialog.getColor(QtGui.QColor(*rgb_color_ints), self,
                                    title)
        if not col.isValid():
            return False, rgb_color_ints, None

        color_float = col.getRgbF()[:3]  # floats
        color_int = [int(colori * 255) for colori in color_float]

        assert isinstance(color_float[0], float), color_float
        assert isinstance(color_int[0], int), color_int

        color_edit.setStyleSheet(
            "QPushButton {"
            "background-color: rgb(%s, %s, %s);" % tuple(color_int) +
            #"border:1px solid rgb(255, 170, 255); "
            "}")
        return True, color_int, color_float

    def on_label_size(self):
        self._label_size = float(self.label_size_edit.text())
        #self.on_apply(force=True)
        #self.min_edit.setText(str(self._default_min))
        #self.min_edit.setStyleSheet("QLineEdit{background: white;}")
        self.update_label_size_color()

    def on_picker_size(self):
        self._picker_size = float(self.picker_size_edit.text())
        if self.win_parent is not None:
            self.win_parent.element_picker_size = self._picker_size / 100.
        #self.on_apply(force=True)

    def on_default_font_size(self):
        self.font_size_edit.setValue(self._default_font_size)
        self.on_set_font(self._default_font_size)

    def on_default_clipping_min(self):
        self.clipping_min_edit.setText(str(self._default_clipping_min))
        self.clipping_min_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_clipping_max(self):
        self.clipping_max_edit.setText(str(self._default_clipping_max))
        self.clipping_max_edit.setStyleSheet("QLineEdit{background: white;}")

    @staticmethod
    def check_float(cell):
        text = cell.text()
        value = float(text)
        return value, True

    @staticmethod
    def check_label_float(cell):
        text = cell.text()
        try:
            value = eval_float_from_string(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
            return value, True
        except ValueError:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    def on_validate(self):
        font_size_value, flag0 = self.check_float(self.font_size_edit)
        label_size_value, flag1 = self.check_float(self.label_size_edit)
        assert isinstance(self.label_color_float[0], float), self.label_color_float
        assert isinstance(self.label_color_int[0], int), self.label_color_int
        picker_size_value, flag2 = self.check_float(self.picker_size_edit)

        clipping_min_value, flag3 = self.check_label_float(self.clipping_min_edit)
        clipping_max_value, flag4 = self.check_label_float(self.clipping_max_edit)

        if all([flag0, flag1, flag2, flag3, flag4]):
            self._label_size = label_size_value
            self._picker_size = picker_size_value

            self.out_data['font_size'] = int(font_size_value)
            self.out_data['clipping_min'] = min(clipping_min_value, clipping_max_value)
            self.out_data['clipping_max'] = max(clipping_min_value, clipping_max_value)
            self.out_data['clicked_ok'] = True
            return True
        return False

    def on_apply(self, force=False):
        passed = self.on_validate()
        if (passed or force) and self.win_parent is not None:
            self.win_parent.on_set_font_size(self.out_data['font_size'])

            #self.win_parent.set_labelsize_color(self._label_size, self.label_color_float)
            #self.win_parent.element_picker_size = self._picker_size / 100.
        if passed and self.win_parent is not None:
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


def _check_color(color_float):
    assert len(color_float) == 3, color_float
    assert isinstance(color_float[0], float), color_float
    color_int = [int(colori * 255) for colori in color_float]
    return color_float, color_int

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
    data = {
        'font_size' : 8,
        'background_color' : (0., 0., 0.), # black
        'text_color' : (0., 1., 0.), # green

        'label_color' : (1., 0., 0.), # red
        'label_size' : 10.,
        'picker_size' : 10.,

        'clipping_min' : 0.,
        'clipping_max' : 10,

        'dim_max' : 502.

    }
    main_window = PreferencesWindow(data)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
