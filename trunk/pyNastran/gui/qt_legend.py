from PyQt4 import QtCore, QtGui

class LegendPropertiesWindow(QtGui.QDialog):

    def __init__(self, data, win_parent=None):
        #Init the base class
        self._default_name = data['name']
        self._default_min = data['min']
        self._default_max = data['max']
        self._default_format = data['format']
        self._default_is_blue_to_red = data['is_blue_to_red']
        self._default_is_discrete = data['is_discrete']

        self.out_data = data

        QtGui.QDialog.__init__(self, win_parent)
        #self.setupUi(self)
        self.setWindowTitle('Legend Properties')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.show()

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
        self.format = QtGui.QLabel("Format (e.g. %s, %.3f, %i, %g, %.6e):")
        self.format_edit = QtGui.QLineEdit(str(self._default_format))
        self.format_button = QtGui.QPushButton("Default")
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

        #checkbox3.setChecked(False)

        # put these in a group
        checkboxs = QtGui.QButtonGroup(self)
        checkboxs.addButton(self.checkbox_blue_to_red)
        checkboxs.addButton(self.checkbox_red_to_blue)

        checkboxs2 = QtGui.QButtonGroup(self)
        checkboxs2.addButton(self.checkbox_continuous)
        checkboxs2.addButton(self.checkbox_discrete)

        # closing
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

        ok_cancel_box = QtGui.QHBoxLayout()
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        if 0:
            vbox1 = QtGui.QVBoxLayout()
            vbox1.addWidget(self.checkbox_blue_to_red)
            vbox1.addWidget(self.checkbox_red_to_blue)

            vbox2 = QtGui.QVBoxLayout()
            vbox2.addWidget(self.checkbox_continuous)
            vbox2.addWidget(self.checkbox_discrete)

            checkboxes = QtGui.QHBoxLayout()
            checkboxes.addLayout(vbox1)
            checkboxes.addLayout(vbox2)

        else:
            grid2 = QtGui.QGridLayout()

            title = QtGui.QLabel("Color Scale:")
            grid2.addWidget(title, 0, 0)
            grid2.addWidget(self.checkbox_blue_to_red, 1, 0)
            grid2.addWidget(self.checkbox_red_to_blue, 2, 0)

            grid2.addWidget(self.checkbox_continuous, 1, 1)
            grid2.addWidget(self.checkbox_discrete, 2, 1)
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

        self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)

    def closeEvent(self, event):
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

    def check_float(self, cell):
        text = cell.text()
        try:
            value = float(text)
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
            cell.setStyleSheet("QLineEdit{background: red;}");
            return None, False

    def on_validate(self):
        name_value, flag0 = self.check_name(self.name_edit)
        min_value, flag1 = self.check_float(self.min_edit)
        max_value, flag2 = self.check_float(self.max_edit)
        format_value, flag3 = self.check_format(self.format_edit)

        if flag0 and flag1 and flag2 and flag3:
            self.out_data['name'] = name_value
            self.out_data['min'] = min_value
            self.out_data['max'] = max_value
            self.out_data['format'] = format_value
            self.out_data['is_blue_to_red'] = self.checkbox_blue_to_red.isChecked()
            self.out_data['is_discrete'] = self.checkbox_discrete.isChecked()

            print("name = %r" % self.name_edit.text())
            print("min = %r" % self.min_edit.text())
            print("max = %r" % self.max_edit.text())
            print("format = %r" % self.format_edit.text())
            return True
        return False

    def on_ok(self):
        passed = self.on_validate()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.close()


def main():
    # Someone is launching this directly
    # Create the QApplication
    app = QtGui.QApplication(sys.argv)
    #The Main window
    main_window = LegendPropertiesWindow()
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
