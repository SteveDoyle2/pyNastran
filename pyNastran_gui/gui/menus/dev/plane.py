from qtpy import QtCore, QtGui
from pyNastran.gui.menus.menu_utils import eval_float_from_string
from pyNastran.gui.utils.qt.pydialog import check_float

class BCMap(QtGui.QDialog):
    def __init__(self, data, win_parent=None):
        """
        +--------+
        | BC Map |
        +--------+---------------+
        |  Origin                |
        |     X       Y      Z   |
        |                        |
        |  Normal                |
        |    NX      NY     NZ   |
        |                        |
        |  Apply Normal          |
        |    NX      NY     NZ   |
        |                        |
        |    Apply   OK  Cancel  |
        +--------+---------------+
        """
        QtGui.QDialog.__init__(self, win_parent)

        self._default_name = 'Plane'
        self.out_data = data
        self._axis = '?'
        self._plane = '?'

        ox_value, oy_value, oz_value = self.out_data['origin']
        nx_value, ny_value, nz_value = self.out_data['normal']
        ax_value, ay_value, az_value = self.out_data['a']
        bx_value, by_value, bz_value = self.out_data['b']

        self._origin_x_default = ox_value
        self._origin_y_default = oy_value
        self._origin_z_default = oz_value

        self._normal_x_default = nx_value
        self._normal_y_default = ny_value
        self._normal_z_default = nz_value

        self._ax_default = ax_value
        self._ay_default = ay_value
        self._az_default = az_value

        self._bx_default = bx_value
        self._by_default = by_value
        self._bz_default = bz_value

        self._default_is_apply = True

        self.name = QtGui.QLabel("Title:")
        self.name_edit = QtGui.QLineEdit(str(self._default_name))
        self.name_button = QtGui.QPushButton("Default")

        self.axis = QtGui.QLabel("Point on ? Axis:")
        self.combo_axis = QtGui.QComboBox(self)
        self.combo_axis.addItem("X")
        self.combo_axis.addItem("Y")
        self.combo_axis.addItem("Z")

        self.plane = QtGui.QLabel("Point on %s? Plane:" % self._axis)
        self.combo_plane = QtGui.QComboBox(self)
        self.combo_plane.addItem("X")
        self.combo_plane.addItem("Y")
        self.combo_plane.addItem("Z")

        self.origin = QtGui.QLabel("Origin:")
        self.origin_x_edit = QtGui.QLineEdit(str(self._origin_x_default))
        self.origin_y_edit = QtGui.QLineEdit(str(self._origin_y_default))
        self.origin_z_edit = QtGui.QLineEdit(str(self._origin_z_default))
        self.origin_default_button = QtGui.QPushButton("Default")
        #self.name_button = QtGui.QPushButton("Default")

        self.normal = QtGui.QLabel("Normal:")
        self.normal_x_edit = QtGui.QLineEdit(str(self._normal_x_default))
        self.normal_y_edit = QtGui.QLineEdit(str(self._normal_y_default))
        self.normal_z_edit = QtGui.QLineEdit(str(self._normal_z_default))
        self.normal_default_button = QtGui.QPushButton("Default")

        self.snap = QtGui.QLabel("Snap Normal:")
        self.snap_normal_xy_button = QtGui.QPushButton("XY Plane")
        self.snap_normal_yz_button = QtGui.QPushButton("YZ Plane")
        self.snap_normal_xz_button = QtGui.QPushButton("XZ Plane")


        self.point_a = QtGui.QLabel("Point on %s Axis:" % self._axis)
        self.ax_edit = QtGui.QLineEdit(str(self._ax_default))
        self.ay_edit = QtGui.QLineEdit(str(self._ay_default))
        self.az_edit = QtGui.QLineEdit(str(self._az_default))

        self.point_b = QtGui.QLabel("Point on %s%s Plane:" % (self._axis, self._plane))
        self.bx_edit = QtGui.QLineEdit(str(self._bx_default))
        self.by_edit = QtGui.QLineEdit(str(self._by_default))
        self.bz_edit = QtGui.QLineEdit(str(self._bz_default))


        self.check_apply = QtGui.QCheckBox("Automatically apply results")
        self.check_apply.setChecked(self._default_is_apply)

        # closing
        self.apply_button = QtGui.QPushButton("Apply")
        if self._default_is_apply:
            self.apply_button.setDisabled(True)

        self.ok_button = QtGui.QPushButton("OK")
        self.cancel_button = QtGui.QPushButton("Cancel")
        self.create_layout()
        self.set_connections()

    def create_layout(self):
        ok_cancel_box = QtGui.QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        # x axis; yplane
        # x axis; zplane
        mode = 2

        grid = QtGui.QGridLayout()

        irow = 0
        grid.addWidget(self.axis, irow, 0)
        grid.addWidget(self.combo_axis, irow, 1)
        irow += 1

        grid.addWidget(self.plane, irow, 0)
        grid.addWidget(self.combo_plane, irow, 1)
        irow += 1

        grid.addWidget(self.origin, irow, 0)
        grid.addWidget(self.origin_x_edit, irow, 1)
        grid.addWidget(self.origin_y_edit, irow, 2)
        grid.addWidget(self.origin_z_edit, irow, 3)
        grid.addWidget(self.origin_default_button, irow, 4)
        irow += 1

        grid.addWidget(self.normal, irow, 0)
        grid.addWidget(self.normal_x_edit, irow, 1)
        grid.addWidget(self.normal_y_edit, irow, 2)
        grid.addWidget(self.normal_z_edit, irow, 3)
        grid.addWidget(self.normal_default_button, irow, 4)
        irow += 1

        grid.addWidget(self.snap, irow, 0)
        grid.addWidget(self.snap_normal_xy_button, irow, 1)
        grid.addWidget(self.snap_normal_yz_button, irow, 2)
        grid.addWidget(self.snap_normal_xz_button, irow, 3)
        irow += 1

        grid.addWidget(self.point_a, irow, 0)
        grid.addWidget(self.ax_edit, irow, 1)
        grid.addWidget(self.ay_edit, irow, 2)
        grid.addWidget(self.az_edit, irow, 3)
        irow += 1

        grid.addWidget(self.point_b, irow, 0)
        grid.addWidget(self.bx_edit, irow, 1)
        grid.addWidget(self.by_edit, irow, 2)
        grid.addWidget(self.bz_edit, irow, 3)
        irow += 1




        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.name)

        vbox.addLayout(grid)
        vbox.addStretch()
        vbox.addWidget(self.check_apply)
        vbox.addLayout(ok_cancel_box)

        #Create central widget, add layout and set
        #central_widget = QtGui.QWidget()
        #central_widget.setLayout(vbox)
        #self.setCentralWidget(central_widget)
        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        self.connect(self.name_button, QtCore.SIGNAL('clicked()'), self.on_default_name)
        #combo.activated[str].connect(self.onActivated)
        self.combo_axis.activated[str].connect(self.on_axis)
        self.combo_plane.activated[str].connect(self.on_plane)
        #self.combo_axis.connect(self.on_axis)
        #self.combo_plane.connect(self.on_plane)

        self.connect(self.origin_default_button, QtCore.SIGNAL('clicked()'), self.on_default_origin)
        self.connect(self.normal_default_button, QtCore.SIGNAL('clicked()'), self.on_default_normal)

        self.connect(self.origin_x_edit, QtCore.SIGNAL('edited()'), self.on_origin_x)
        self.connect(self.origin_y_edit, QtCore.SIGNAL('edited()'), self.on_origin_y)
        self.connect(self.origin_z_edit, QtCore.SIGNAL('edited()'), self.on_origin_z)

        self.connect(self.snap_normal_xy_button, QtCore.SIGNAL('clicked()'), self.on_snap_xy)
        self.connect(self.snap_normal_yz_button, QtCore.SIGNAL('clicked()'), self.on_snap_yz)
        self.connect(self.snap_normal_xz_button, QtCore.SIGNAL('clicked()'), self.on_snap_xz)


        self.connect(self.check_apply, QtCore.SIGNAL('clicked()'), self.on_check_apply)

        self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)

    def on_axis(self, text):
        #print(self.combo_axis.itemText())
        self._axis = str(text)
        self.plane.setText('Point on %s? Plane:' % self._axis)
        self.point_a.setText('Point on %s Axis:' % self._axis)
        self.point_b.setText('Point on %s%s Plane:' % (self._axis, self._plane))

    def on_plane(self, text):
        self._plane = str(text)
        self.point_b.setText('Point on %s%s Plane:' % (self._axis, self._plane))

    def on_check_apply(self):
        is_checked = self.check_apply.isChecked()
        self.apply_button.setDisabled(is_checked)

    def on_snap_xy(self):
        self.normal_x_edit.setText(str(0.0))
        self.normal_y_edit.setText(str(0.0))
        self.normal_z_edit.setText(str(1.0))

    def on_snap_yz(self):
        self.normal_x_edit.setText(str(1.0))
        self.normal_y_edit.setText(str(0.0))
        self.normal_z_edit.setText(str(0.0))

    def on_snap_xz(self):
        self.normal_x_edit.setText(str(0.0))
        self.normal_y_edit.setText(str(1.0))
        self.normal_z_edit.setText(str(0.0))

    def on_origin_x(self):
        _on_float(self.origin_x_edit)

    def on_origin_y(self):
        _on_float(self.origin_y_edit)

    def on_origin_z(self):
        _on_float(self.origin_z_edit)

    def on_normal_x(self):
        _on_float(self.normal_x_edit)

    def on_normal_y(self):
        _on_float(self.normal_y_edit)

    def on_normal_z(self):
        _on_float(self.normal_z_edit)


    def closeEvent(self, event):
        event.accept()

    def on_default_name(self):
        self.name_edit.setText(str(self._default_name))
        self.name_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_origin(self):
        self.origin_x_edit.setText(str(self._origin_x_default))
        self.origin_y_edit.setText(str(self._origin_y_default))
        self.origin_z_edit.setText(str(self._origin_z_default))
        self.origin_x_edit.setStyleSheet("QLineEdit{background: white;}")
        self.origin_y_edit.setStyleSheet("QLineEdit{background: white;}")
        self.origin_z_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_normal(self):
        self.normal_x_edit.setText(str(self._normal_x_default))
        self.normal_y_edit.setText(str(self._normal_y_default))
        self.normal_z_edit.setText(str(self._normal_z_default))
        self.normal_x_edit.setStyleSheet("QLineEdit{background: white;}")
        self.normal_y_edit.setStyleSheet("QLineEdit{background: white;}")
        self.normal_z_edit.setStyleSheet("QLineEdit{background: white;}")

    def check_name(self, cell):
        text = str(cell.text()).strip()
        if text:
            cell.setStyleSheet("QLineEdit{background: white;}")
            return text, True
        else:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    def on_validate(self):
        name_value, flag0 = self.check_name(self.name_edit)
        ox_value, flag1 = check_float(self.origin_x_edit)
        oy_value, flag2 = check_float(self.origin_y_edit)
        oz_value, flag3 = check_float(self.origin_z_edit)

        nx_value, flag4 = check_float(self.normal_x_edit)
        ny_value, flag5 = check_float(self.normal_y_edit)
        nz_value, flag6 = check_float(self.normal_z_edit)

        if flag0 and flag1 and flag2 and flag3 and flag4 and flag5 and flag6:
            self.out_data['origin'] = [ox_value, oy_value, oz_value]
            self.out_data['normal'] = [nx_value, ny_value, nz_value]
            self.out_data['clicked_ok'] = True

            return True
        return False

    def on_apply(self):
        passed = self.on_validate()
        if passed:
            self.win_parent.create_plane(self.out_data)
        return passed

    def on_ok(self):
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.close()

def _on_float(field):
    try:
        eval_float_from_string(field.text())
        field.setStyleSheet("QLineEdit{background: white;}")
    except ValueError:
        field.setStyleSheet("QLineEdit{background: red;}")


def main():  # pragma: no cover
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
        'name' : 'plane_test',
        'origin': [1., 2., 3.],
        'normal' : [4., 5., 6.],
        'a' : [7., 8., 9.],
        'b' : [10., 11., 12.],
        'is_horizontal' : False,
    }
    main_window = BCMap(d, win_parent=None)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":  # pragma: no cover
    main()
