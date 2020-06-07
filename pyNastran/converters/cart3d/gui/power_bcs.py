"""
References:
-----------
https://wiki.python.org/moin/PyQt/Distinguishing%20between%20click%20and%20double%20click
http://www.saltycrane.com/blog/2007/12/pyqt-43-qtableview-qabstracttablemodel/
http://stackoverflow.com/questions/12152060/how-does-the-keypressevent-method-work-in-this-program
"""
#from qtpy import QtCore, QtGui
from qtpy.QtWidgets import (
    QCheckBox, QLabel, QLineEdit, QDialog, QColorDialog, QSpinBox,
    QDoubleSpinBox, QRadioButton, QHBoxLayout, QButtonGroup,
    QVBoxLayout, QGridLayout, QPushButton, QHeaderView,
)
from qtpy import QtGui
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.menus.manage_actors import Model, SingleChoiceQTableView
#from pyNastran.gui.qutils.pydialog import check_format


class BreakSurfaceMenu(QDialog):
    def __init__(self, data, win_parent=None):
        """
        +-----------------+
        | Break Surfaces  |
        +-----------------+------+
        |  EngineInlet           |
        |  EngineOutlet          |
        |                        |
        |  Name      EngineInlet |
        |  RegionMode * RegionID |
        |             * All      |
        |                        |
        |  AllowedRegions:       |
        |    Region ID      3    |
        |                        |
        |  PickMode  * All       |
        |  Pick Mode  x On/Off   |
        |  Pick Angle   20 deg   |
        |                        |
        |         Revert         |
        |     RenumberRegions    |
        |         Close          |
        +------------------------+
        """
        QDialog.__init__(self, win_parent)
        self.setWindowTitle('Break Surface')

        #default
        self.win_parent = win_parent
        self.out_data = data

        self.points = data['points']

        self.keys = sorted(self.points.keys())
        keys = self.keys
        nrows = len(keys)

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

        header = view.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        #header.setResizeMode(QtGui.QHeaderView.Stretch)
        self.table = view

        #self.representation = actor_obj.representation
        #print('rep =', self.representation)

        table = self.table
        #headers = [QtCore.QString('Groups')]

        header = table.horizontalHeader()
        header.setStretchLastSection(True)

        #----------------------------------------------
        #self._default_is_apply = False

        self.mode_header = QLabel("Mode:")

        nregions_max = 10
        pick_angle = 20.0
        region_id = 4
        all_regions = True
        self.region_id = QLabel("Region ID:")
        self.region_id_edit = QSpinBox(self)
        self.region_id_edit.setRange(1, nregions_max)
        self.region_id_edit.setSingleStep(1)
        self.region_id_edit.setValue(region_id)

        self.pick_angle = QLabel("Pick Angle:")
        self.pick_angle_edit = QDoubleSpinBox(self)
        self.pick_angle_edit.setRange(0.0, 360.0)
        self.pick_angle_edit.setDecimals(3)
        self.pick_angle_edit.setSingleStep(0.5)
        self.pick_angle_edit.setValue(pick_angle)

        # region IDs/all
        self.checkbox_region_ids = QCheckBox("Region IDs")
        self.checkbox_region_all = QCheckBox("All Regions")
        self.checkbox_region_all.setChecked(all_regions)
        self.checkbox_region_ids.setChecked(not all_regions)

        # pick mode
        self.checkbox_pick_mode = QCheckBox("Pick Mode  (Off=label)")
        self.checkbox_pick_mode.setChecked(False)

        #----------------------------------------------
        self.nodes_header = QLabel("Single Node:")
        self.name = QLabel("ID:")
        self.name_edit = QLineEdit('Node %i' % name)
        self.name_edit.setDisabled(True)

        #----------------------------------------------
        self.location_x = QLabel("X:")
        self.location_x_edit = QLineEdit('X')

        self.location_y = QLabel("Y:")
        self.location_y_edit = QLineEdit('Y')

        self.location_z = QLabel("Z:")
        self.location_z_edit = QLineEdit('Z')

        #----------------------------------------------

        # remove these...
        self.description_edit = QLineEdit('Description')
        self.coord_edit = QSpinBox()

        #----------------------------------------------

        # closing
        #if self._default_is_apply:
            #self.apply_button.setDisabled(True)

        self.close_button = QPushButton("Close")

        self.create_layout()
        #self.set_connections()

    def update_active_key(self, index):
        name = self.active_key
        old_obj = self.out_data['points'][name]
        #self.active_key
        #self.points[self.active_key]
        old_obj[0] = str(self.description_edit.text())
        #old_obj.coord = self.description_edit.value()
        #old_obj.description = self.description_edit.value()
        #old_obj.description = self.description_edit.value()

        str_name = str(index.data())
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

        self.location_x_edit.setText(str(point[3]))
        self.location_y_edit.setText(str(point[4]))
        self.location_z_edit.setText(str(point[5]))
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

        checkboxs = QButtonGroup(self)
        checkboxs.addButton(self.checkbox_region_all)
        checkboxs.addButton(self.checkbox_region_ids)

        vbox1 = QVBoxLayout()
        vbox1.addWidget(self.checkbox_region_all)
        vbox1.addWidget(self.checkbox_region_ids)
        #vbox1.addLayout(checkboxs)

        #-----------------------------------------
        irow = 0
        grid2.addWidget(self.name, irow, 0)
        grid2.addWidget(self.name_edit, irow, 1)
        irow += 1

        #grid2.addWidget(self.name, irow, 0)
        grid2.addWidget(self.description_edit, irow, 1)
        irow += 1

        grid2.addWidget(self.location_x, irow, 0)
        grid2.addWidget(self.location_x_edit, irow, 1)
        irow += 1

        grid2.addWidget(self.location_y, irow, 0)
        grid2.addWidget(self.location_y_edit, irow, 1)
        irow += 1

        grid2.addWidget(self.location_z, irow, 0)
        grid2.addWidget(self.location_z_edit, irow, 1)
        irow += 1

        #|  Name      EngineInlet |
        #|  RegionMode * RegionID |
        #|             * All      |
        #|                        |
        #|  AllowedRegions:       |
        #|    Region ID      3    |
        #|                        |
        #|  PickMode  * All       |
        #|  Pick Mode  x On/Off   |
        #|  Pick Angle   20 deg   |
        #|                        |
        #|         Revert         |
        #|     RenumberRegions    |
        #|         Close          |

        grid2.addWidget(self.region_id, irow, 0)
        grid2.addWidget(self.region_id_edit, irow, 1)
        irow += 1

        #grid2.addWidget(self.pick_mode, irow, 0)
        grid2.addWidget(self.checkbox_pick_mode, irow, 0)
        irow += 1

        grid2.addWidget(self.pick_angle, irow, 0)
        grid2.addWidget(self.pick_angle_edit, irow, 1)
        irow += 1

        #grid2.addWidget(self.pi, irow, 0)
        #grid2.addLayout(coord_type_layout, irow, 1)
        #irow += 1

        #grid2.addWidget(self.location, irow, 0)
        #grid2.addLayout(location_layout, irow, 1)
        #irow += 1

        #------------------------------------

        vbox = QVBoxLayout()
        vbox.addLayout(grid1)
        #vbox.addLayout(vbox1)
        #vbox.addStretch()
        vbox.addWidget(self.table)

        vbox.addLayout(grid2)
        vbox.addStretch()
        #vbox.addWidget(self.check_apply)
        vbox.addLayout(cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        self.opacity_edit.clicked.connect(self.on_opacity)
        self.point_size.clicked.connect(self.on_point_size)
        self.color_edit.clicked.connect(self.on_color)
        self.checkbox_show.clicked.connect(self.on_show)
        self.checkbox_hide.clicked.connect(self.on_hide)

        #self.connect(self.description_edit, QtCore.SIGNAL("valueChanged(int)"), self.on_description)
        #self.connect(self.coord_edit, QtCore.SIGNAL("valueChanged(int)"), self.on_coord)
        self.radio_rectangular.clicked.connect(self.on_coord_type)
        self.radio_cylindrical.clicked.connect(self.on_coord_type)
        self.radio_spherical.clicked.connect(self.on_coord_type)

        self.location_x_edit.clicked.connect(self.on_location_x)
        self.location_y_edit.clicked.connect(self.on_location_y)
        self.location_z_edit.clicked.connect(self.on_location_z)

        self.close_button.clicked.connect(self.on_close)

        #self.connect(self.check_apply, QtCore.SIGNAL('clicked()'), self.on_check_apply)

        #self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        #self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        self.close_button.clicked.connect(self.on_close)


    def on_color(self):
        obj = self.out_data['point_properties']
        rgb_color_ints = obj.color

        msg = 'Points'
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
        #self.name_edit.setStyleSheet("QLineEdit{background: white;}")

    #def check_float(self, cell):
        #text = cell.text()
        #try:
            #value = eval_float_from_string(text)
            #cell.setStyleSheet("QLineEdit{background: white;}")
            #return value, True
        #except ValueError:
            #cell.setStyleSheet("QLineEdit{background: red;}")
            #return None, False

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
        old_obj.point_size = self.point_size_edit.value()
        old_obj.opacity = self.opacity_edit.value()
        old_obj.is_visible = self.checkbox_show.isChecked()
        return True
        #name_value, flag0 = self.check_name(self.name_edit)
        #ox_value, flag1 = self.check_float(self.transparency_edit)
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


def main():
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)


    import sys
    from qtpy.QtWidgets import QApplication

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

    point_properties = AltGeometry(parent, 'point_properties', color=green,
                                   point_size=5, opacity=1.0)

    points = {
        1 : ['LERoot', 0, 'R', 1.0, 2.0, 3.0],
        2 : ['LETip', 42, 'S', 2.0, 3.0, 4.0],
    }
    coords = {
        0 : ['R'],
        2 : ['C'],
        42 : ['S'],
    }
    data = {
        'point_properties' : point_properties,
        'points' : points,
        'active_point' : 2,
        'coords' : coords,
    }
    main_window = BreakSurfaceMenu(data, win_parent=None)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == '__main__':   # pragma: no cover
    main()
