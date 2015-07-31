from six import iteritems
from PyQt4 import QtCore, QtGui
from pyNastran.gui.qt_files.menu_utils import eval_float_from_string

#class Table1(QtGui.QTableWidget):
    #def __init__(self, parent, gui):
        #QtGui.QTableWidget.__init__(self, parent)
        #self.gui = gui

    #def mouseReleaseEvent(self, event):
        #if event.button() == QtCore.Qt.RightButton:
            #self.rightClickMenu(event)

    #def rightClickMenu(self,  event):
        #pos = event.pos
        #self.gui.ui.menuEdit.popup(QtGui.QCursor.pos())

#class TreeView(QtGui.QTreeView):
    #def edit(self, index, trigger, event):
        #if trigger == QtGui.QAbstractItemView.DoubleClicked:
            #print 'DoubleClick Killed!'
            #return False
        #return QtGui.QTreeView.edit(self, index, trigger, event)

class EditGroupProperties(QtGui.QDialog):
    def __init__(self, data, win_parent=None):
        """
        +------------------+
        | Edit Actor Props |
        +------------------+------+
        |  Name1                  |
        |  Name2                  |
        |  Name3                  |
        |  Name4                  |
        |                         |
        |  Color          box     |
        |  Line_Thickness 2.0     |
        |  Opacity        0.5     |
        |                         |
        |    Apply   OK   Cancel  |
        +-------------------------+
        """
        QtGui.QDialog.__init__(self, win_parent)
        self.setWindowTitle('Edit Group Properties')

        #default
        #self._default_name = 'Plane'
        self.out_data = data
        #self._axis = '?'
        #self._plane = '?'

        keys = sorted(data.keys())
        nrows = len(keys)
        self.active_key = keys[0]

        self.table = QtGui.QTableWidget()
        names_text = []
        for key, actor_obj in sorted(iteritems(data)):
            name_text = QtGui.QTableWidgetItem(str(key))
            #self.connect(lb, SIGNAL('doubleClicked()'), self.someMethod)
            #self.connect(name_text, QtCore.SIGNAL('doubleClicked()'), self.on_name_select)
            #self.connect(self.table, QtCore.SIGNAL("itemDoubleClicked (QListWidgetItem *)"), self.on_name_select)
            #self.connect(name_text, QtCore.SIGNAL('rightClicked()'), self.on_name_select)
            #QtCore.QObject.connect(
                #self.table,
                #QtCore.SIGNAL('itemChanged(QTableWidgetItem*)'),
                #self.on_name_select,
            #)
            #self.table.itemChanged.connect(someFunc)
            #item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled );
            names_text.append(name_text)

        self.table.cellDoubleClicked.connect(self.on_name_select)
        #self.table.doubleClicked.connect(self.on_name_select)

        actor_obj = data[self.active_key]
        name = actor_obj.name
        line_thickness = actor_obj.line_thickness
        transparency = actor_obj.transparency
        color = actor_obj.color

        table = self.table
        table.setRowCount(nrows)
        table.setColumnCount(1)
        headers = [QtCore.QString('Groups')]
        table.setHorizontalHeaderLabels(headers)

        header = table.horizontalHeader()
        header.setStretchLastSection(True)
        for iname, name_text in enumerate(names_text):
            # row, col, value
            table.setItem(iname, 0, name_text)
        table.resizeRowsToContents()


        self._default_is_apply = False

        self.name = QtGui.QLabel("Name:")
        self.name_edit = QtGui.QLineEdit(str(name))
        #self.name_button = QtGui.QPushButton("Default")

        self.color = QtGui.QLabel("Color:")
        self.color_edit = QtGui.QPushButton()
        #self.color_edit.setFlat(True)

        color = self.out_data[self.active_key].color
        qcolor = QtGui.QColor()
        qcolor.setRgb(*color)
        print('color =%s' % str(color))
        palette = QtGui.QPalette(self.color_edit.palette()) # make a copy of the palette
        #palette.setColor(QtGui.QPalette.Active, QtGui.QPalette.Base, \
                         #qcolor)
        palette.setColor(QtGui.QPalette.Background, QtGui.QColor('blue'))  # ButtonText
        #palette.setColor(QtGui.QPalette.ButtonText, qcolor)
        self.color_edit.setPalette(palette) # assign new palette

        self.color_edit.setStyleSheet("QPushButton {"
                                      "background-color: rgb(%s, %s, %s);" % tuple(color) +
                                      #"border:1px solid rgb(255, 170, 255); "
                                      "}")


        self.transparency = QtGui.QLabel("Transparency:")
        self.transparency_edit = QtGui.QDoubleSpinBox(self)
        self.transparency_edit.setRange(0.0, 1.0)
        self.transparency_edit.setDecimals(1)
        self.transparency_edit	.setSingleStep(0.1)

        self.line_thickness = QtGui.QLabel("Line Thickness:")
        self.line_thickness_edit = QtGui.QDoubleSpinBox(self)
        self.line_thickness_edit.setRange(0.0, 1.0)
        self.line_thickness_edit.setDecimals(1)
        self.line_thickness_edit.setSingleStep(0.1)

        # closing
        self.apply_button = QtGui.QPushButton("Apply")
        #if self._default_is_apply:
            #self.apply_button.setDisabled(True)

        self.ok_button = QtGui.QPushButton("OK")
        self.cancel_button = QtGui.QPushButton("Cancel")

        self.create_layout()
        self.set_connections()

    def on_name_select(self):
        print('on_name_select')
        return

    def create_layout(self):
        ok_cancel_box = QtGui.QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        grid = QtGui.QGridLayout()

        irow = 0
        grid.addWidget(self.name, irow, 0)
        grid.addWidget(self.name_edit, irow, 1)
        irow += 1

        grid.addWidget(self.color, irow, 0)
        grid.addWidget(self.color_edit, irow, 1)
        irow += 1

        grid.addWidget(self.transparency, irow, 0)
        grid.addWidget(self.transparency_edit, irow, 1)
        irow += 1

        grid.addWidget(self.line_thickness, irow, 0)
        grid.addWidget(self.line_thickness_edit, irow, 1)
        irow += 1

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.table)
        vbox.addLayout(grid)
        vbox.addStretch()
        #vbox.addWidget(self.check_apply)
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        ##self.connect(self.name_button, QtCore.SIGNAL('clicked()'), self.on_default_name)

        self.connect(self.transparency_edit, QtCore.SIGNAL('clicked()'), self.on_transparency)
        self.connect(self.line_thickness, QtCore.SIGNAL('clicked()'), self.on_line_thickness)
        self.connect(self.color_edit, QtCore.SIGNAL('clicked()'), self.on_color)


        #self.connect(self.check_apply, QtCore.SIGNAL('clicked()'), self.on_check_apply)

        self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)

    def on_color(self):
        name = self.active_key
        obj = self.out_data[name]
        rgb_color_floats = obj.color
        #c = [int(255 * i) for i in rgb_color_floats]
        c = rgb_color_floats
        msg = name
        col = QtGui.QColorDialog.getColor(QtGui.QColor(*c), self, "Choose a %s color" % msg)
        if col.isValid():
            color = col.getRgbF()[:3]
            obj.color = color
            print('new_color =', color)
            self.color_edit.setStyleSheet("QPushButton {"
                                          "background-color: rgb(%s, %s, %s);" % tuple(obj.color) +
                                          #"border:1px solid rgb(255, 170, 255); "
                                          "}")

    def on_line_thickness(self):
        name = self.active_key
        print('on_line_thickness')

    def on_transparency(self):
        name = self.active_key
        print('on_transparency')

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

    def _on_float(self, field):
        try:
            eval_float_from_string(field.text())
            field.setStyleSheet("QLineEdit{background: white;}")
        except ValueError:
            field.setStyleSheet("QLineEdit{background: red;}")

    def closeEvent(self, event):
        event.accept()

    def on_default_name(self):
        self.name_edit.setText(str(self._default_name))
        self.name_edit.setStyleSheet("QLineEdit{background: white;}")

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

    def on_validate(self):
        name_value, flag0 = self.check_name(self.name_edit)
        ox_value, flag1 = self.check_float(self.transparency_edit)
        if flag0 and flag1:
            #self.out_data['origin'] = [ox_value, oy_value, oz_value]
            #self.out_data['normal'] = [nx_value, ny_value, nz_value]
            self.out_data['clicked_ok'] = True

            return True
        return False

    def on_apply(self):
        passed = self.on_validate()
        if passed:
            self.win_parent.update_properties(self.out_data)
        return passed

    def on_ok(self):
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.close()


#class GeometryHandle(object):
    #def __init__(self, parent=None):
        #self.parent = parent
        #self.grids = []
        #self.mappers = []
        #self.geometry_actors = []
        #self.colors = []
        #self.line_thickness = []

class AltGeometry(object):
    def __init__(self, parent, name, grid=None, color=None, line_thickness=None, transparency=None):
        if color is None:
            color = (0., 0.2, 0.3)
        if line_thickness is None:
            line_thickness = 1.0
        if transparency is None:
            transparency = 0.0

        self.parent = parent
        self.name = name
        self.grid = grid
        self.mapper = None
        self.geometry_actor = None
        self.color = color
        self.line_thickness = line_thickness

        self._transparency = transparency

    @property
    def opacity(self):
        assert 0.0 <= self._transparency <= 1.0, self._transparency
        return 1.0 - self._transparency

    @opacity.setter
    def opacity(self, opacity):
        assert 0.0 <= opacity <= 1.0, opacity
        self._transparency = 1.0 - opacity

    @property
    def transparency(self):
        assert 0.0 <= self._transparency <= 1.0, self._transparency
        return self._transparency

    @transparency.setter
    def transparency(self, transparency):
        assert 0.0 <= transparency <= 1.0, transparency
        self._transparency = transparency

    @property
    def color(self):
        return self._color

    @color.setter
    def color(self, color):
        assert len(color) == 3, color
        if isinstance(color[0], int):
            assert isinstance(color[0], int), color[0]
            assert isinstance(color[1], int), color[1]
            assert isinstance(color[2], int), color[2]
            self._color = color
        else:
            assert isinstance(color[0], float), color[0]
            assert isinstance(color[1], float), color[1]
            assert isinstance(color[2], float), color[2]
            self._color = [color[0] * 255, color[1] * 255, color[2] * 255]


    def set_color(self, color, mode='rgb'):
        assert mode == 'rgb', mode
        self.color = color
        assert len(color)
        self.mode = 'rgb'

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
    #g = GeometryHandle()
    #g.add('main', color=(0, 0, 0), line_thickness=0.0)
    #g.get_grid('name')
    #g.set_color('name')
    #g.set_grid('name')
    #g.set_grid('name')
    parent = app
    red = (255, 0, 0)
    blue = (0, 0, 255)
    d = {
        'main' : AltGeometry(parent, 'main', color=red, line_thickness=1.0, transparency=0.0),
        'caero' : AltGeometry(parent, 'caero', color=blue, line_thickness=1.0, transparency=0.0),
        'caero1' : AltGeometry(parent, 'caero', color=blue, line_thickness=1.0, transparency=0.0),
        'caero2' : AltGeometry(parent, 'caero', color=blue, line_thickness=1.0, transparency=0.0),
    }
    main_window = EditGroupProperties(d, win_parent=None)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
