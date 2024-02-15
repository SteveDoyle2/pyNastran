"""
References:
-----------
https://wiki.python.org/moin/PyQt/Distinguishing%20between%20click%20and%20double%20click
http://www.saltycrane.com/blog/2007/12/pyqt-43-qtableview-qabstracttablemodel/
http://stackoverflow.com/questions/12152060/how-does-the-keypressevent-method-work-in-this-program
"""
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.menus.edit_geometry_properties.manage_actors import Model, SingleChoiceQTableView
from pyNastran.gui.qt_version import qt_int as qt_version
from pyNastran.gui.utils.qt.checks.qlineedit import QLINEEDIT_GOOD, QLINEEDIT_ERROR

from qtpy import QtCore
from qtpy.QtGui import QColor, QPalette
from qtpy.QtWidgets import (
    QDialog, QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
    QSpinBox, QCheckBox, QRadioButton, QButtonGroup, QDoubleSpinBox, QLineEdit,
    QColorDialog, QHeaderView)


class EditNodeProperties(QDialog):
    def __init__(self, data, win_parent=None):
        """
        +-----------------+
        | Edit Node Props |
        +-----------------+------+
        |  LEwingTip             |
        |  Node2                 |
        |  Node3                 |
        |  Node4                 |
        |                        |
        |  All Nodes:            |
        |    Color     red       |
        |    PointSize 3         |
        |    Opacity   0.3       |
        |    Show/Hide           |
        |                        |
        |  Name        LEwingTip |
        |  Location    X Y Z     |
        |  Coord       0         |
        |  CoordType   R, C, S   |
        |                        |
        |   Previous     Next    |
        |                        |
        |          Close         |
        +------------------------+
        """
        QDialog.__init__(self, win_parent)
        self.setWindowTitle('Edit Node Properties')

        #default
        self.win_parent = win_parent
        self.out_data = data

        point_properties = data['point_properties']
        print(point_properties)
        #name = point_properties.name
        point_size = point_properties.point_size
        opacity = point_properties.opacity
        color = point_properties.color
        show = point_properties.is_visible

        self.points = data['points']

        self.keys = sorted(self.points.keys())
        keys = self.keys
        #nrows = len(keys)

        active_point = data['active_point']
        #self.active_key = keys[0]
        self.active_key = active_point
        name = self.active_key
        description = self.points[self.active_key][0]

        self._use_old_table = False
        items = ['Node %i' % val for val in keys]
        header_labels = ['Nodes']
        table_model = Model(items, header_labels, self)
        view = SingleChoiceQTableView(self) #Call your custom QTableView here
        view.setModel(table_model)
        view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table = view

        #self.representation = actor_obj.representation
        #print('rep =', self.representation)

        table = self.table
        #headers = [QtCore.QString('Groups')]

        header = table.horizontalHeader()
        header.setStretchLastSection(True)

        #----------------------------------------------
        #self._default_is_apply = False

        self.color = QLabel("Color:")
        self.color_edit = QPushButton()
        #self.color_edit.setFlat(True)

        color = self.out_data['point_properties'].color
        opacity = self.out_data['point_properties'].opacity
        show = self.out_data['point_properties'].is_visible
        #color = self.out_data[self.active_key].color
        qcolor = QColor()
        qcolor.setRgb(*color)
        #print('color =%s' % str(color))
        palette = QPalette(self.color_edit.palette()) # make a copy of the palette
        #palette.setColor(QPalette.Active, QPalette.Base, \
                    #qcolor)
        palette.setColor(QPalette.Background, QColor('blue'))  # ButtonText
        self.color_edit.setPalette(palette)

        self.color_edit.setStyleSheet("QPushButton {"
                                      "background-color: rgb(%s, %s, %s);" % tuple(color) +
                                      #"border:1px solid rgb(255, 170, 255); "
                                      "}")

        self.all_nodes_header = QLabel("All Nodes:")
        self.point_size = QLabel("Point Size:")
        self.point_size_edit = QSpinBox(self)
        self.point_size_edit.setRange(1, 10)
        self.point_size_edit.setSingleStep(1)
        self.point_size_edit.setValue(point_size)

        self.opacity = QLabel('Opacity:')
        self.opacity_edit = QDoubleSpinBox(self)
        self.opacity_edit.setRange(0.1, 1.0)
        self.opacity_edit.setDecimals(2)
        self.opacity_edit.setSingleStep(0.05)
        self.opacity_edit.setValue(opacity)

        # show/hide
        self.checkbox_show = QCheckBox("Show")
        self.checkbox_hide = QCheckBox("Hide")
        self.checkbox_show.setChecked(show)
        self.checkbox_hide.setChecked(not show)

        #----------------------------------------------
        self.nodes_header = QLabel("Single Node:")
        self.name = QLabel("ID:")
        self.name_edit = QLineEdit('Node %i' % name)
        self.name_edit.setDisabled(True)

        self.description = QLabel("Description:")
        self.description_edit = QLineEdit(str(description))
        #self.description_edit.setDisabled(True)

        location_x = 0.1
        location_y = 0.1
        location_z = 0.1
        self.location = QLabel("Location:")
        self.location_x_edit = QDoubleSpinBox(self)
        self.location_y_edit = QDoubleSpinBox(self)
        self.location_z_edit = QDoubleSpinBox(self)
        #self.location_x_edit.setDecimals(1)
        delta_x = abs(location_x) / 100. if location_x != 0.0 else 0.1
        delta_y = abs(location_y) / 100. if location_y != 0.0 else 0.1
        delta_z = abs(location_z) / 100. if location_z != 0.0 else 0.1
        self.location_x_edit.setSingleStep(delta_x)
        self.location_y_edit.setSingleStep(delta_y)
        self.location_z_edit.setSingleStep(delta_z)
        self.location_x_edit.setValue(location_x)
        self.location_y_edit.setValue(location_y)
        self.location_z_edit.setValue(location_z)

        self.coord = QLabel("Coord:")
        self.coord_edit = QSpinBox(self)
        self.coord_edit.setRange(0, 99999999)
        #self.coord_edit.setSingleStep(1)
        self.coord_edit.setValue(0)

        self.coord_type = QLabel("Coord Type:")
        #----------------------------------------------

        # closing
        #if self._default_is_apply:
            #self.apply_button.setDisabled(True)

        self.close_button = QPushButton("Close")

        self.create_layout()
        self.set_connections()

    def update_active_key(self, index):
        name = self.active_key
        old_obj = self.out_data['points'][name]
        #self.active_key
        #self.points[self.active_key]
        old_obj[0] = str(self.description_edit.text())
        #old_obj.coord = self.description_edit.value()
        #old_obj.description = self.description_edit.value()
        #old_obj.description = self.description_edit.value()

        str_name = str(index.data().toString())
        name = int(str_name[5:])
        #i = self.keys.index(self.active_key)

        self.active_key = name
        point = self.points[self.active_key]

        #1  : ['LERoot', 0, 'R', 1.0, 2.0, 3.0],
        self.name_edit.setText(str(self.active_key))
        self.description_edit.setText(point[0])

        self.coord_edit.setValue(point[1])
        if point[2] == 'R':
            self.radio_rectangular.setChecked(True)
        elif point[2] == 'C':
            self.radio_cylindrical.setChecked(True)
        elif point[2] == 'S':
            self.radio_spherical.setChecked(True)

        self.location_x_edit.setValue(point[3])
        self.location_y_edit.setValue(point[4])
        self.location_z_edit.setValue(point[5])
        #obj = self.out_data[name]
        #point_size = obj.point_size
        #opacity = obj.opacity
        #representation = obj.representation
        #is_visible = obj.is_visible

        #self.opacity_edit.setValue(opacity)
        #self.checkbox_show.setChecked(is_visible)
        #self.checkbox_hide.setChecked(not is_visible)

    #def on_name_select(self):
        #print('on_name_select')
        #return

    def create_layout(self):
        cancel_box = QHBoxLayout()
        cancel_box.addWidget(self.close_button)

        grid1 = QGridLayout()
        grid2 = QGridLayout()

        #-----------------------------------------
        # setup
        self.radio_rectangular = QRadioButton('Rectangular')
        self.radio_cylindrical = QRadioButton('Cylindrical')
        self.radio_spherical = QRadioButton('Spherical')

        coord_type_layout = QHBoxLayout()
        coord_type_layout.addWidget(self.radio_rectangular)
        coord_type_layout.addWidget(self.radio_cylindrical)
        coord_type_layout.addWidget(self.radio_spherical)

        location_layout = QHBoxLayout()
        location_layout.addWidget(self.location_x_edit)
        location_layout.addWidget(self.location_y_edit)
        location_layout.addWidget(self.location_z_edit)

        checkboxs = QButtonGroup(self)
        checkboxs.addButton(self.checkbox_show)
        checkboxs.addButton(self.checkbox_hide)

        vbox1 = QVBoxLayout()
        vbox1.addWidget(self.checkbox_show)
        vbox1.addWidget(self.checkbox_hide)
        #vbox1.addLayout(checkboxs)

        #-----------------------------------------
        irow = 0
        grid1.addWidget(self.all_nodes_header, irow, 0)
        irow += 1

        grid1.addWidget(self.color, irow, 0)
        grid1.addWidget(self.color_edit, irow, 1)
        irow += 1

        grid1.addWidget(self.opacity, irow, 0)
        grid1.addWidget(self.opacity_edit, irow, 1)
        irow += 1

        grid1.addWidget(self.point_size, irow, 0)
        grid1.addWidget(self.point_size_edit, irow, 1)
        irow += 1

        #-----------------------------------------
        irow = 0
        grid2.addWidget(self.nodes_header, irow, 0)
        irow += 1

        grid2.addWidget(self.name, irow, 0)
        grid2.addWidget(self.name_edit, irow, 1)
        irow += 1

        grid2.addWidget(self.description, irow, 0)
        grid2.addWidget(self.description_edit, irow, 1)
        irow += 1

        #|  All Nodes:            |
        #|    Color     red       |
        #|    PointSize 3         |
        #|    Opacity   0.3       |
        #|    Show/Hide           |
        #|                        |
        #|  Name        LEwingTip |
        #|  Location    X Y Z     |
        #|  Coord       0         |
        #|  CoordType   R, C, S   |
        #|                        |
        #|   Previous     Next    |
        grid2.addWidget(self.coord, irow, 0)
        grid2.addWidget(self.coord_edit, irow, 1)
        irow += 1

        grid2.addWidget(self.coord_type, irow, 0)
        grid2.addLayout(coord_type_layout, irow, 1)
        irow += 1

        grid2.addWidget(self.location, irow, 0)
        grid2.addLayout(location_layout, irow, 1)
        irow += 1

        #------------------------------------

        vbox = QVBoxLayout()
        vbox.addLayout(grid1)
        vbox.addLayout(vbox1)
        vbox.addStretch()
        vbox.addWidget(self.table)

        vbox.addLayout(grid2)
        vbox.addStretch()
        #vbox.addWidget(self.check_apply)
        vbox.addLayout(cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        self.opacity_edit.valueChanged.connect(self.on_opacity)
        #self.connect(self.point_size, QtCore.SIGNAL('clicked()'), self.on_point_size)
        #self.point_size_edit.clicked.connect(self.on_point_size)  # need to update name??
        self.color_edit.clicked.connect(self.on_color)
        self.checkbox_show.clicked.connect(self.on_show)
        self.checkbox_hide.clicked.connect(self.on_hide)

        self.description_edit.textEdited.connect(self.on_description)  # valueChanged????
        self.coord_edit.valueChanged.connect(self.on_coord)
        self.radio_rectangular.clicked.connect(self.on_coord_type)
        self.radio_cylindrical.clicked.connect(self.on_coord_type)
        self.radio_spherical.clicked.connect(self.on_coord_type)

        self.location_x_edit.valueChanged.connect(self.on_location_x)
        self.location_y_edit.valueChanged.connect(self.on_location_y)
        self.location_z_edit.valueChanged.connect(self.on_location_z)
        #self.connect(self.check_apply, QtCore.SIGNAL('clicked()'), self.on_check_apply)
        #self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        #self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        self.close_button.clicked.connect(self.on_close)

    def on_color(self):
        obj = self.out_data['point_properties']
        rgb_color_ints = obj.color

        msg = 'Points'
        col = QColorDialog.getColor(QColor(*rgb_color_ints), self, "Choose a %s color" % msg)
        if col.isValid():
            color = col.getRgbF()[:3]
            obj.color = color
            #print('new_color =', color)
            self.color_edit.setStyleSheet("QPushButton {"
                                          "background-color: rgb(%s, %s, %s);" % tuple(obj.color) +
                                          #"border:1px solid rgb(255, 170, 255); "
                                          "}")

    def on_show(self):
        is_checked = self.checkbox_show.isChecked()
        self.out_data['point_properties'].is_visible = is_checked

    def on_hide(self):
        is_checked = self.checkbox_hide.isChecked()
        self.out_data['point_properties'].is_visible = not is_checked

    def on_point_size(self):
        point_size = self.point_size_edit.value()
        self.out_data['point_properties'].point_size = point_size

    def on_opacity(self):
        opacity = self.opacity_edit.value()
        self.out_data['point_properties'].opacity = opacity


    def on_description(self):
        #1 : ['LERoot', 0, 'R', 1.0, 2.0, 3.0],
        name = self.active_key
        description = self.description_edit.value()
        self.out_data['points'][name][0] = description

    def on_coord(self):
        #1 : ['LERoot', 0, 'R', 1.0, 2.0, 3.0],
        name = self.active_key
        coord_id = self.coord_edit.value()
        self.out_data['points'][name][1] = coord_id

    def on_coord_type(self):
        #1 : ['LERoot', 0, 'R', 1.0, 2.0, 3.0],
        name = self.active_key
        if self.radio_rectangular.isChecked():
            coord_type = 'R'
        elif self.radio_cylindrical.isChecked():
            coord_type = 'C'
        elif self.radio_spherical.isChecked():
            coord_type = 'S'
        else:
            raise NotImplementedError()
        self.out_data['points'][name][2] = coord_type

    def on_location_x(self):
        #1 : ['LERoot', 0, 'R', 1.0, 2.0, 3.0],
        name = self.active_key
        value = self.coord_edit.value()
        self.out_data['points'][name][3] = value

    def on_location_y(self):
        #1 : ['LERoot', 0, 'R', 1.0, 2.0, 3.0],
        name = self.active_key
        value = self.coord_edit.value()
        self.out_data['points'][name][4] = value

    def on_location_z(self):
        #1 : ['LERoot', 0, 'R', 1.0, 2.0, 3.0],
        name = self.active_key
        value = self.coord_edit.value()
        self.out_data['points'][name][5] = value

    def closeEvent(self, event):
        event.accept()

    #def on_default_name(self):
        #self.name_edit.setText(str(self._default_name))
        #self.name_edit.setStyleSheet(QLINEEDIT_GOOD)

    #def check_name(self, cell):
        #text = str(cell.text()).strip()
        #if len(text):
            #cell.setStyleSheet(QLINEEDIT_GOOD)
            #return text, True
        #else:
            #cell.setStyleSheet(QLINEEDIT_ERROR)
            #return None, False

    def on_validate(self):
        self.out_data['clicked_ok'] = True
        self.out_data['clicked_cancel'] = False

        old_obj = self.out_data[self.active_key]
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
            self.win_parent.on_update_gui_nodes(self.out_data)
        return passed

    def on_ok(self):
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_close(self):
        self.out_data['clicked_close'] = True
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
    green = (0, 255, 0)

    point_properties = AltGeometry(parent, 'point_properties', color=green, point_size=5, opacity=1.0)

    points = {
        1 : ['LERoot', 0, 'R', 1.0, 2.0, 3.0],
        2 : ['LETip', 42, 'S', 2.0, 3.0, 4.0],
    }
    coords = {
        0 : ['R'],
        2 : ['C'],
        42 : ['S'],
    }
    d = {
        'point_properties' : point_properties,
        'points' : points,
        'active_point' : 2,
        'coords' : coords,
    }
    main_window = EditNodeProperties(d, win_parent=None)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":  # pragma: no cover
    main()
