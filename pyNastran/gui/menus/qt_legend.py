from __future__ import print_function
from six import string_types
from PyQt4 import QtCore, QtGui
from pyNastran.gui.qt_files.menu_utils import eval_float_from_string


class LegendPropertiesWindow(QtGui.QDialog):
    """
    +-------------------+
    | Legend Properties |
    +-----------------------+
    | Title  ______ Default |
    | Min    ______ Default |
    | Max    ______ Default |
    | Format ______ Default |
    | Scale  ______ Default |
    | Number of Colors ____ | (TODO)
    | Number of Labels ____ | (TODO)
    | Label Size       ____ | (TODO)
    |                       |
    | x Min/Max (Blue->Red) |
    | o Max/Min (Red->Blue) |
    |                       |
    | x Vertical/Horizontal |
    | x Show/Hide           |
    |                       |
    |    Apply OK Cancel    |
    +-----------------------+
    """

    def __init__(self, data, win_parent=None):
        self.win_parent = win_parent
        #Init the base class
        self._updated_legend = False
        self._icase = data['icase']
        self._default_name = data['name']
        self._default_min = data['min']
        self._default_max = data['max']
        self._format = data['format']
        self._scale = data['scale']
        self._default_is_blue_to_red = data['is_blue_to_red']
        self._default_is_discrete = data['is_discrete']
        self._default_is_horizontal = data['is_horizontal']
        self._default_is_shown = data['is_shown']

        self._default_format = data['default_format']
        self._default_scale = data['default_scale']
        self._default_icase = self._icase

        self.out_data = data

        QtGui.QDialog.__init__(self, win_parent)
        #self.setupUi(self)
        self.setWindowTitle('Legend Properties')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.show()

    def update_legend(self, icase, name,
                      min_value, max_value, data_format, scale,
                      default_title, default_min_value, default_max_value, default_data_format, default_scale,
                      is_blue_to_red, is_horizontal_scalar_bar):
        """
        We need to update the legend if there's been a result change request
        """
        if icase != self._default_icase:
            self._default_icase = icase
            self._default_name = default_title
            self._default_min = default_min_value
            self._default_max = default_max_value
            self._default_format = default_data_format
            self._default_is_blue_to_red = is_blue_to_red
            self._default_is_discrete = True
            self._default_is_horizontal = is_horizontal_scalar_bar
            self._default_scale = default_scale

            assert isinstance(scale, float), 'scale=%r' % scale
            assert isinstance(default_scale, float), 'default_scale=%r' % default_scale
            if self._default_scale == 0.0:
                self.scale_edit.setEnabled(False)
                self.scale_button.setEnabled(False)
            else:
                self.scale_edit.setEnabled(True)
                self.scale_button.setEnabled(True)

            #self.on_default_name()
            #self.on_default_min()
            #self.on_default_max()
            #self.on_default_format()
            #self.on_default_scale()
            # reset defaults
            self.name_edit.setText(name)
            self.name_edit.setStyleSheet("QLineEdit{background: white;}")

            self.min_edit.setText(str(min_value))
            self.min_edit.setStyleSheet("QLineEdit{background: white;}")

            self.max_edit.setText(str(max_value))
            self.max_edit.setStyleSheet("QLineEdit{background: white;}")

            self.format_edit.setText(str(data_format))
            self.format_edit.setStyleSheet("QLineEdit{background: white;}")

            self.scale_edit.setText(str(scale))
            self.scale_edit.setStyleSheet("QLineEdit{background: white;}")


    def create_widgets(self):
        # Name
        self.name = QtGui.QLabel("Title:")
        self.name_edit = QtGui.QLineEdit(str(self._default_name))
        self.name_button = QtGui.QPushButton("Default")

        # Min
        self.min = QtGui.QLabel("Min:")
        self.min_edit = QtGui.QLineEdit(str(self._default_min))
        self.min_button = QtGui.QPushButton("Default")

        # Max
        self.max = QtGui.QLabel("Max:")
        self.max_edit = QtGui.QLineEdit(str(self._default_max))
        self.max_button = QtGui.QPushButton("Default")

        # Format
        self.format = QtGui.QLabel("Format (e.g. %.3f, %g, %.6e):")
        self.format_edit = QtGui.QLineEdit(str(self._format))
        self.format_button = QtGui.QPushButton("Default")

        # Scale
        self.scale = QtGui.QLabel("Scale:")
        self.scale_edit = QtGui.QLineEdit(str(self._scale))
        self.scale_button = QtGui.QPushButton("Default")
        if self._default_scale == 0.0:
            self.scale_edit.setEnabled(False)
            self.scale_button.setEnabled(False)
        #tip = QtGui.QToolTip()
        #tip.setTe
        #self.format_edit.toolTip(tip)

        # red/blue or blue/red
        self.checkbox_blue_to_red = QtGui.QCheckBox("Min -> Blue; Max -> Red")
        self.checkbox_red_to_blue = QtGui.QCheckBox("Min -> Red; Max -> Blue")
        self.checkbox_blue_to_red.setChecked(self._default_is_blue_to_red)

        # continuous / discrete
        self.checkbox_continuous = QtGui.QCheckBox("Continuous")
        self.checkbox_discrete = QtGui.QCheckBox("Discrete")
        self.checkbox_discrete.setChecked(self._default_is_discrete)
        self.checkbox_continuous.setDisabled(True)
        self.checkbox_discrete.setDisabled(True)

        # horizontal / vertical
        self.checkbox_horizontal = QtGui.QCheckBox("Horizontal")
        self.checkbox_vertical = QtGui.QCheckBox("Vertical")
        self.checkbox_horizontal.setChecked(self._default_is_horizontal)
        self.checkbox_vertical.setChecked(not self._default_is_horizontal)

        # on / off
        self.checkbox_show = QtGui.QCheckBox("Show")
        self.checkbox_hide = QtGui.QCheckBox("Hide")
        print('_default_is_shown =', self._default_is_shown)
        self.checkbox_show.setChecked(self._default_is_shown)
        self.checkbox_hide.setChecked(not self._default_is_shown)

        #checkbox3.setChecked(False)

        # put these in a group
        checkboxs = QtGui.QButtonGroup(self)
        checkboxs.addButton(self.checkbox_blue_to_red)
        checkboxs.addButton(self.checkbox_red_to_blue)

        checkboxs2 = QtGui.QButtonGroup(self)
        checkboxs2.addButton(self.checkbox_continuous)
        checkboxs2.addButton(self.checkbox_discrete)

        checkboxs3 = QtGui.QButtonGroup(self)
        checkboxs3.addButton(self.checkbox_vertical)
        checkboxs3.addButton(self.checkbox_horizontal)

        checkboxs4 = QtGui.QButtonGroup(self)
        checkboxs4.addButton(self.checkbox_show)
        checkboxs4.addButton(self.checkbox_hide)

        # closing
        self.apply_button = QtGui.QPushButton("Apply")
        self.ok_button = QtGui.QPushButton("OK")
        self.cancel_button = QtGui.QPushButton("Cancel")

    def create_layout(self):
        grid = QtGui.QGridLayout()
        grid.addWidget(self.name, 0, 0)
        grid.addWidget(self.name_edit, 0, 1)
        grid.addWidget(self.name_button, 0, 2)

        grid.addWidget(self.min, 1, 0)
        grid.addWidget(self.min_edit, 1, 1)
        grid.addWidget(self.min_button, 1, 2)

        grid.addWidget(self.max, 2, 0)
        grid.addWidget(self.max_edit, 2, 1)
        grid.addWidget(self.max_button, 2, 2)

        grid.addWidget(self.format, 3, 0)
        grid.addWidget(self.format_edit, 3, 1)
        grid.addWidget(self.format_button, 3, 2)

        grid.addWidget(self.scale, 4, 0)
        grid.addWidget(self.scale_edit, 4, 1)
        grid.addWidget(self.scale_button, 4, 2)

        ok_cancel_box = QtGui.QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        if 0:
            vbox1 = QtGui.QVBoxLayout()
            vbox1.addWidget(self.checkbox_blue_to_red)
            vbox1.addWidget(self.checkbox_red_to_blue)

            vbox2 = QtGui.QVBoxLayout()
            vbox2.addWidget(self.checkbox_continuous)
            vbox2.addWidget(self.checkbox_discrete)

            vbox3 = QtGui.QVBoxLayout()
            vbox3.addWidget(self.checkbox_vertical)
            vbox3.addWidget(self.checkbox_horizontal)

            checkboxes = QtGui.QHBoxLayout()
            checkboxes.addLayout(vbox1)
            checkboxes.addLayout(vbox2)
            checkboxes.addLayout(vbox3)
        else:
            grid2 = QtGui.QGridLayout()

            title = QtGui.QLabel("Color Scale:")
            grid2.addWidget(title, 0, 0)
            grid2.addWidget(self.checkbox_blue_to_red, 1, 0)
            grid2.addWidget(self.checkbox_red_to_blue, 2, 0)

            grid2.addWidget(self.checkbox_continuous, 1, 1)
            grid2.addWidget(self.checkbox_discrete, 2, 1)

            grid2.addWidget(self.checkbox_vertical, 1, 2)
            grid2.addWidget(self.checkbox_horizontal, 2, 2)

            grid2.addWidget(self.checkbox_show, 1, 3)
            grid2.addWidget(self.checkbox_hide, 2, 3)
            #grid2.setSpacing(0)

        vbox = QtGui.QVBoxLayout()
        vbox.addLayout(grid)
        #vbox.addLayout(checkboxes)
        vbox.addLayout(grid2)
        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)

        #Create central widget, add layout and set
        #central_widget = QtGui.QWidget()
        #central_widget.setLayout(vbox)
        #self.setCentralWidget(central_widget)
        self.setLayout(vbox)

    def set_connections(self):
        self.connect(self.name_button, QtCore.SIGNAL('clicked()'), self.on_default_name)
        self.connect(self.min_button, QtCore.SIGNAL('clicked()'), self.on_default_min)
        self.connect(self.max_button, QtCore.SIGNAL('clicked()'), self.on_default_max)
        self.connect(self.format_button, QtCore.SIGNAL('clicked()'), self.on_default_format)
        self.connect(self.scale_button, QtCore.SIGNAL('clicked()'), self.on_default_scale)

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

    def on_default_name(self):
        self.name_edit.setText(str(self._default_name))
        self.name_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_min(self):
        self.min_edit.setText(str(self._default_min))
        self.min_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_max(self):
        self.max_edit.setText(str(self._default_max))
        self.max_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_format(self):
        self.format_edit.setText(str(self._default_format))
        self.format_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_scale(self):
        self.scale_edit.setText(str(self._default_scale))
        self.scale_edit.setStyleSheet("QLineEdit{background: white;}")

    def check_float(self, cell):
        text = cell.text()
        try:
            value = eval_float_from_string(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
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
        name_value, flag0 = self.check_name(self.name_edit)
        min_value, flag1 = self.check_float(self.min_edit)
        max_value, flag2 = self.check_float(self.max_edit)
        format_value, flag3 = self.check_format(self.format_edit)
        scale_value, flag4 = self.check_float(self.scale_edit)
        if 'i' in format_value:
            format_value = '%i'

        if flag0 and flag1 and flag2 and flag3 and flag4:
            self.out_data['name'] = name_value
            self.out_data['min'] = min_value
            self.out_data['max'] = max_value
            self.out_data['format'] = format_value
            self.out_data['scale'] = scale_value
            self.out_data['is_blue_to_red'] = self.checkbox_blue_to_red.isChecked()
            self.out_data['is_discrete'] = self.checkbox_discrete.isChecked()
            self.out_data['is_horizontal'] = self.checkbox_horizontal.isChecked()
            self.out_data['is_shown'] = self.checkbox_show.isChecked()
            self.out_data['clicked_ok'] = True
            self.out_data['close'] = True

            #print("name = %r" % self.name_edit.text())
            #print("min = %r" % self.min_edit.text())
            #print("max = %r" % self.max_edit.text())
            #print("format = %r" % self.format_edit.text())
            return True
        return False

    def on_apply(self):
        passed = self.on_validate()
        if passed:
            self.win_parent._apply_legend(self.out_data)
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
    app = QtGui.QApplication(sys.argv)
    #The Main window
    d = {
        'icase' : 1,
        'name' : 'asdf',
        'min' : 0.,
        'max' : 10,
        'scale' : 1e-12,
        'default_scale' : 1.0,

        'default_format' : '%s',
        'format' : '%g',
        'is_blue_to_red': True,
        'is_discrete' : False,
        'is_horizontal' : False,
        'is_shown' : True,
    }
    main_window = LegendPropertiesWindow(d)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
