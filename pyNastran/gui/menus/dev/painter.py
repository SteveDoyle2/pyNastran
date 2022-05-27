from pyNastran.gui.menus.edit_geometry_properties.manage_actors import SingleChoiceQTableView, AltGeometry, Model
from qtpy import QtCore, QtGui
from qtpy.QtWidgets import (
    QLabel, QLineEdit, QPushButton, QCheckBox, QSpinBox,
    QDoubleSpinBox, QColorDialog, QApplication,
    QHBoxLayout, QGridLayout, QVBoxLayout, QButtonGroup)
#from pyNastran.utils.locale import func_str
from pyNastran.gui.utils.qt.pydialog import QDialog #, QFloatEdit
from pyNastran.gui.utils.qt.version import Background
#from pyNastran.gui.utils.qt.
#check_float


class EditBoundaryConditions(QDialog):
    def __init__(self, data, win_parent=None):
        """
        +---------+
        | Painter |
        +---------+---------+
        |  Name1            |
        |  Name2            |
        |  Name3            |
        |  Name4            |
        |                   |
        | x by cell   ___   |
        | x by angle  ___   |
        |                   |
        |   box elements    |
        |                   |
        |   Add    Remove   |
        |                   |
        | Apply  OK  Cancel |
        +-------------------+
        """
        QDialog.__init__(self, win_parent)
        self.setWindowTitle('Edit Boundary Conditions')

        #default
        self.win_parent = win_parent
        self.out_data = data

        self.keys = sorted(data.keys())
        keys = self.keys
        nrows = len(keys)
        self.active_key = keys[0]

        self._use_old_table = False
        items = keys
        header_labels = ['Groups']
        table_model = Model(items, header_labels, self)
        view = SingleChoiceQTableView(self) #Call your custom QTableView here
        view.setModel(table_model)
        #view.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)  #TODO: fixme
        self.table = view

        actor_obj = data[self.active_key]
        name = actor_obj.name
        line_width = actor_obj.line_width
        point_size = actor_obj.point_size
        opacity = actor_obj.opacity
        color = actor_obj.color
        show = actor_obj.is_visible
        self.representation = actor_obj.representation

        table = self.table
        #headers = [QtCore.QString('Groups')]

        header = table.horizontalHeader()
        header.setStretchLastSection(True)

        self._default_is_apply = False
        self.name = QLabel("Name:")
        self.name_edit = QLineEdit(str(name))
        self.name_edit.setDisabled(True)

        self.color = QLabel("Color:")
        self.color_edit = QPushButton()
        #self.color_edit.setFlat(True)

        color = self.out_data[self.active_key].color
        qcolor = QtGui.QColor()
        qcolor.setRgb(*color)
        #print('color =%s' % str(color))
        palette = QtGui.QPalette(self.color_edit.palette()) # make a copy of the palette
        #palette.setColor(QtGui.QPalette.Active, QtGui.QPalette.Base, \
                         #qcolor)
        palette.setColor(Background, QtGui.QColor('blue'))  # ButtonText
        self.color_edit.setPalette(palette)

        self.color_edit.setStyleSheet("QPushButton {"
                                      "background-color: rgb(%s, %s, %s);" % tuple(color) +
                                      #"border:1px solid rgb(255, 170, 255); "
                                      "}")


        self.opacity = QLabel('Opacity:')
        self.opacity_edit = QDoubleSpinBox(self)
        self.opacity_edit.setRange(0.1, 1.0)
        self.opacity_edit.setDecimals(1)
        self.opacity_edit.setSingleStep(0.1)
        self.opacity_edit.setValue(opacity)

        self.line_width = QLabel('Line Width:')
        self.line_width_edit = QSpinBox(self)
        self.line_width_edit.setRange(1, 10)
        self.line_width_edit.setSingleStep(1)
        self.line_width_edit.setValue(line_width)
        if self.representation in ['point', 'surface']:
            self.line_width.setEnabled(False)
            self.line_width_edit.setEnabled(False)

        self.point_size = QLabel('Point Size:')
        self.point_size_edit = QSpinBox(self)
        self.point_size_edit.setRange(1, 10)
        self.point_size_edit.setSingleStep(1)
        self.point_size_edit.setValue(point_size)
        if self.representation in ['wire', 'surface']:
            self.point_size.setEnabled(False)
            self.point_size_edit.setEnabled(False)

        # show/hide
        self.checkbox_show = QCheckBox("Show")
        self.checkbox_hide = QCheckBox("Hide")
        self.checkbox_show.setChecked(show)
        self.checkbox_hide.setChecked(not show)

        # closing
        self.apply_button = QPushButton("Apply")
        #if self._default_is_apply:
            #self.apply_button.setDisabled(True)

        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")

        self.create_layout()
        self.set_connections()

    def update_active_key(self, index):
        old_obj = self.out_data[self.active_key]
        old_obj.line_width = self.line_width_edit.value()
        old_obj.point_size = self.point_size_edit.value()
        old_obj.opacity = self.opacity_edit.value()
        old_obj.is_visible = self.checkbox_show.isChecked()

        name = index.data()
        #i = self.keys.index(self.active_key)

        self.active_key = name
        self.name_edit.setText(name)
        obj = self.out_data[name]
        line_width = obj.line_width
        point_size = obj.point_size
        opacity = obj.opacity
        representation = obj.representation
        is_visible = obj.is_visible

        self.color_edit.setStyleSheet("QPushButton {"
                                      "background-color: rgb(%s, %s, %s);" % tuple(obj.color) +
                                      #"border:1px solid rgb(255, 170, 255); "
                                      "}")
        self.line_width_edit.setValue(line_width)
        self.point_size_edit.setValue(point_size)
        if self.representation != representation:
            self.representation = representation

            if self.representation == 'point':
                self.point_size.setEnabled(True)
                self.point_size_edit.setEnabled(True)
            else:
                self.point_size.setEnabled(False)
                self.point_size_edit.setEnabled(False)

            if self.representation == ['wire']:
                self.line_width.setEnabled(True)
                self.line_width_edit.setEnabled(True)
            else:
                self.line_width.setEnabled(False)
                self.line_width_edit.setEnabled(False)
            #if self.representation in ['wire', 'surface']:

        self.opacity_edit.setValue(opacity)
        self.checkbox_show.setChecked(is_visible)
        self.checkbox_hide.setChecked(not is_visible)

    #def on_name_select(self):
        #print('on_name_select')
        #return

    def create_layout(self):
        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        grid = QGridLayout()

        irow = 0
        grid.addWidget(self.name, irow, 0)
        grid.addWidget(self.name_edit, irow, 1)
        irow += 1

        grid.addWidget(self.color, irow, 0)
        grid.addWidget(self.color_edit, irow, 1)
        irow += 1

        grid.addWidget(self.opacity, irow, 0)
        grid.addWidget(self.opacity_edit, irow, 1)
        irow += 1

        grid.addWidget(self.line_width, irow, 0)
        grid.addWidget(self.line_width_edit, irow, 1)
        irow += 1

        grid.addWidget(self.point_size, irow, 0)
        grid.addWidget(self.point_size_edit, irow, 1)
        irow += 1

        checkboxs = QButtonGroup(self)
        checkboxs.addButton(self.checkbox_show)
        checkboxs.addButton(self.checkbox_hide)

        vbox = QVBoxLayout()
        vbox.addWidget(self.table)
        vbox.addLayout(grid)

        if 0:
            vbox.addWidget(self.checkbox_show)
            vbox.addWidget(self.checkbox_hide)
        else:
            vbox1 = QVBoxLayout()
            vbox1.addWidget(self.checkbox_show)
            vbox1.addWidget(self.checkbox_hide)
            vbox.addLayout(vbox1)

        vbox.addStretch()
        #vbox.addWidget(self.check_apply)
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        self.opacity_edit.valueChanged.connect(self.on_opacity)        # clicked?
        self.line_width_edit.valueChanged.connect(self.on_line_width)  # clicked?
        self.point_size_edit.valueChanged.connect(self.on_point_size)  # clicked?
        self.color_edit.clicked.connect(self.on_color)
        self.checkbox_show.clicked.connect(self.on_show)
        self.checkbox_hide.clicked.connect(self.on_hide)
        #self.opacity_edit.clicked.connect(self.on_opacity)

        #self.connect(self.opacity_edit, QtCore.SIGNAL('clicked()'), self.on_opacity)
        #self.connect(self.line_width, QtCore.SIGNAL('clicked()'), self.on_line_width)
        #self.connect(self.point_size, QtCore.SIGNAL('clicked()'), self.on_point_size)
        #self.connect(self.color_edit, QtCore.SIGNAL('clicked()'), self.on_color)
        #self.connect(self.checkbox_show, QtCore.SIGNAL('clicked()'), self.on_show)
        #self.connect(self.checkbox_hide, QtCore.SIGNAL('clicked()'), self.on_hide)
        #self.connect(self.check_apply, QtCore.SIGNAL('clicked()'), self.on_check_apply)

        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)

        #self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        #self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        #self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)
        #self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.close()

    def closeEvent(self, event):
        self.out_data['close'] = True
        event.accept()

    def on_color(self):
        name = self.active_key
        obj = self.out_data[name]
        rgb_color_ints = obj.color

        msg = name
        col = QColorDialog.getColor(QtGui.QColor(*rgb_color_ints), self, "Choose a %s color" % msg)
        if col.isValid():
            color = col.getRgbF()[:3]
            obj.color = color
            #print('new_color =', color)
            self.color_edit.setStyleSheet("QPushButton {"
                                          "background-color: rgb(%s, %s, %s);" % tuple(obj.color) +
                                          #"border:1px solid rgb(255, 170, 255); "
                                          "}")

    def on_show(self):
        name = self.active_key
        is_checked = self.checkbox_show.isChecked()
        self.out_data[name].is_visible = is_checked

    def on_hide(self):
        name = self.active_key
        is_checked = self.checkbox_hide.isChecked()
        self.out_data[name].is_visible = not is_checked

    def on_line_width(self):
        name = self.active_key
        line_width = self.line_width_edit.value()
        self.out_data[name].line_width = line_width

    def on_point_size(self):
        name = self.active_key
        point_size = self.point_size_edit.value()
        self.out_data[name].point_size = point_size

    def on_opacity(self):
        name = self.active_key
        opacity = self.opacity_edit.value()
        self.out_data[name].opacity = opacity

    #def on_axis(self, text):
        ##print(self.combo_axis.itemText())
        #self._axis = str(text)
        #self.plane.setText('Point on %s? Plane:' % self._axis)
        #self.point_a.setText('Point on %s Axis:' % self._axis)
        #self.point_b.setText('Point on %s%s Plane:' % (self._axis, self._plane))

    #def on_plane(self, text):
        #self._plane = str(text)
        #self.point_b.setText('Point on %s%s Plane:' % (self._axis, self._plane))

    def on_check_apply(self):
        is_checked = self.check_apply.isChecked()
        self.apply_button.setDisabled(is_checked)

    def _on_float(self, field):
        try:
            eval_float_from_string(field.text())
            field.setStyleSheet("QLineEdit{background: white;}")
        except ValueError:
            field.setStyleSheet("QLineEdit{background: red;}")

    #def on_default_name(self):
        #self.name_edit.setText(str(self._default_name))
        #self.name_edit.setStyleSheet("QLineEdit{background: white;}")

    #def check_name(self, cell):
        #text = str(cell.text()).strip()
        #if len(text):
            #cell.setStyleSheet("QLineEdit{background: white;}")
            #return text, True
        #else:
            #cell.setStyleSheet("QLineEdit{background: red;}")
            #return None, False

    def on_validate(self):
        self.out_data['clicked_ok'] = True
        self.out_data['clicked_cancel'] = False

        old_obj = self.out_data[self.active_key]
        old_obj.line_width = self.line_width_edit.value()
        old_obj.point_size = self.point_size_edit.value()
        old_obj.opacity = self.opacity_edit.value()
        old_obj.is_visible = self.checkbox_show.isChecked()
        return True
        #name_value, flag0 = self.check_name(self.name_edit)
        #ox_value, flag1 = check_float(self.transparency_edit)
        #if flag0 and flag1:
            #self.out_data['clicked_ok'] = True
            #return True
        #return False

    def on_apply(self):
        passed = self.on_validate()
        if passed:
            self.win_parent.on_update_geometry_properties(self.out_data)
        return passed

    def on_ok(self):
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.out_data['clicked_ok'] = False
        self.out_data['clicked_cancel'] = True
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
    #g = GeometryHandle()
    #g.add('main', color=(0, 0, 0), line_thickness=0.0)
    #g.get_grid('name')
    #g.set_color('name')
    #g.set_grid('name')
    #g.set_grid('name')
    parent = app
    red = (255, 0, 0)
    blue = (0, 0, 255)
    green = (0, 255, 0)
    purple = (255, 0, 255)
    d = {
        'caero1' : AltGeometry(parent, 'caero', color=green, line_width=3, opacity=0.2),
        'caero2' : AltGeometry(parent, 'caero', color=purple, line_width=4, opacity=0.3),
        'caero' : AltGeometry(parent, 'caero', color=blue, line_width=2, opacity=0.1),
        'main' : AltGeometry(parent, 'main', color=red, line_width=1, opacity=0.0),
    }
    main_window = EditBoundaryConditions(d, win_parent=None)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":  # pragma: no cover
    main()
