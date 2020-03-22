"""
References:
-----------
https://wiki.python.org/moin/PyQt/Distinguishing%20between%20click%20and%20double%20click
http://www.saltycrane.com/blog/2007/12/pyqt-43-qtableview-qabstracttablemodel/
http://stackoverflow.com/questions/12152060/how-does-the-keypressevent-method-work-in-this-program
"""
from pyNastran.gui.limits import MAX_POINT_SIZE, MAX_LINE_WIDTH
from pyNastran.gui.qt_version import qt_int as qt_version

from qtpy.QtCore import Qt#, QVariant
from qtpy import QtCore, QtGui
from qtpy.QtWidgets import (
    QLabel, QLineEdit, QPushButton, QTableView, QApplication,
    QDoubleSpinBox, QSlider, QSpinBox, QCheckBox, QHBoxLayout, QGridLayout, QVBoxLayout,
    QButtonGroup, QColorDialog, QAbstractItemView,
)

#from pyNastran.gui.menus.menu_utils import eval_float_from_string
from pyNastran.gui.utils.qt.pydialog import PyDialog
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.gui_objects.coord_properties import CoordProperties


class SingleChoiceQTableView(QTableView):
    def __init__(self, *args, **kwargs):
        self.parent2 = args[0]

        # name is not required
        self.name = None
        if 'name' in kwargs:
            self.name = kwargs['name']
            del kwargs['name']

        #super(SingleChoiceQTableView, self).__init__()
        QTableView.__init__(self, *args, **kwargs) #Use QTableView constructor
        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.setSelectionBehavior(QAbstractItemView.SelectRows)

    def get_data(self):
        return self.model().items

    def update_data(self, data):  # not needed?
        #items = self.getModel() # just the data...
        #self.model().change_data(data) # doesn't work...
        items = self.model().items
        header_labels = self.model().header_labels

        parent = self.parent2
        table_model = Model(data, header_labels, parent=parent)
        self.setModel(table_model)

    def getModel(self):  # not needed?
        model = self.model() #tableView.model()
        return model.items
        #data = []
        #for row in range(model.rowCount()):
            #data.append([])
            #for column in range(model.columnCount()):
                #index = model.index(row, column)
                ## We suppose data are strings

                #role = QtCore.Qt.DisplayRole
                #data[row].append(str(model.data(index, role).toString()))
        #return data

    #def mouseDoubleClickEvent(self, event):
        #self.last = "Double Click"
        #index = self.currentIndex()
        #self.parent2.update_active_key(index)

    #def mousePressEvent(self, event):
    def mouseReleaseEvent(self, event):
        index = self.currentIndex()
        #print('index.row() =', index.row())
        irow = index.row()
        self.selectRow(irow)
        if irow == -1:  # null case
            return
        if self.name is None:
            self.parent2.update_active_key(index)
        else:
            self.parent2.update_active_key(self.name, index)

        #index = self.currentIndex()
        #self.parent2.update_active_key(index)

    #def clicked(self, event):
        #index = self.currentIndex()
        #self.parent2.update_active_key(index)

    #def performSingleClickAction(self):
        #index = self.currentIndex()
        #self.parent2.update_active_key(index)

    #def performSingleClickAction(self):
        #if self.last == "Click":
            #self.message = "Click"
            #self.update()

    def keyPressEvent(self, event): #Reimplement the event here, in your case, do nothing
        #if event.key() == QtCore.Qt.Key_Escape:
            #self.close()
        #return
        key = event.key()
        if key == Qt.Key_Delete:
            index = self.currentIndex()
            parent = self.parent()
            irow = index.row()
            if self.name is None:
                parent.on_delete(irow)
            else:
                #print('parent =', parent, type(parent))
                #print('parent2 =', self.parent2, type(self.parent2))
                self.parent2.on_delete(self.name, irow)

            #self.parent().on_delete(index.row()) # old
            #print('pressed delete')
        elif key in [Qt.Key_Up, Qt.Key_Left]:
            index = self.currentIndex()
            nrows = len(self.getModel())
            irow = max(0, index.row() - 1)
            irow = self.selectRow(irow)
            #print('pressed up; nrows=%s' % nrows)
        elif key in [Qt.Key_Down, Qt.Key_Right]:
            index = self.currentIndex()
            nrows = len(self.getModel())
            irow = min(nrows - 1, index.row() + 1)
            irow = self.selectRow(irow)
            #print('pressed down; nrows=%s' % nrows)
        else:
            print('pressed %r' % key)

class Model(QtCore.QAbstractTableModel):

    def __init__(self, items, header_labels, parent=None, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)
        self.items = items
        self.header_labels = header_labels

    #def adding_row(index):
        ## http://stackoverflow.com/questions/13109128/pyqt-qabstracttablemodel-never-updates-when-rows-are-added
        #self.beginInsertRows(self.createIndex(0, 0), index, index)
        #print('adding ', index)

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.items)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return 1

    def data(self, index, role):
        if not index.isValid():
            #return QVariant()  # TODO: is this right?
            return
        elif role != QtCore.Qt.DisplayRole:
            #return QVariant()  # TODO: is this right?
            return

        row = index.row()
        if row < len(self.items):
            return str(self.items[row])  # TODO: is this right?
        else:
            #return QVariant()  # TODO: is this right?
            return

    def flags(self, index):
        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable #| QtCore.Qt.ItemIsEditable

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self.header_labels[section]
        return QtCore.QAbstractTableModel.headerData(self, section, orientation, role)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.close()


class EditGeometryProperties(PyDialog):
    force = True
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
        |  Active_Name    main    |
        |  Color          box     |
        |  Line_Width     2       |
        |  Point_Size     2       |
        |  Bar_Scale      2       |
        |  Opacity        0.5     |
        |  Show/Hide              |
        |                         |
        |    Apply   OK   Cancel  |
        +-------------------------+
        """
        PyDialog.__init__(self, data, win_parent)
        self.set_font_size(data['font_size'])
        del self.out_data['font_size']
        self.setWindowTitle('Edit Geometry Properties')
        self.allow_update = True

        #default
        #self.win_parent = win_parent
        #self.out_data = data

        self.keys = sorted(data.keys())
        self.keys = data.keys()
        keys = self.keys
        items = list(keys)

        #nrows = len(keys)
        active_key = 'main'
        if 'main' not in items:
            active_key = items[0]
        self.active_key = active_key

        header_labels = ['Groups']
        table_model = Model(items, header_labels, self)
        view = SingleChoiceQTableView(self) #Call your custom QTableView here
        view.setModel(table_model)

        self.table = view
        #self.opacity_edit.valueChanged.connect(self.on_opacity)
        #mListWidget, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(itemClicked(QListWidgetItem*)));
        #self.table.itemClicked.connect(self.table.mouseDoubleClickEvent)

        actor_obj = data[self.active_key]
        if isinstance(actor_obj, CoordProperties):
            opacity = 1.0
            representation = 'coord'
            show = actor_obj.is_visible
            color = None
            line_width = 0
            point_size = 0
            bar_scale = 0
            name = 'Coord'
        else:
            name = actor_obj.name
            line_width = actor_obj.line_width
            point_size = actor_obj.point_size
            bar_scale = actor_obj.bar_scale
            opacity = actor_obj.opacity
            color = actor_obj.color
            show = actor_obj.is_visible
            representation = actor_obj.representation
        self.representation = representation

        # table
        header = self.table.horizontalHeader()
        header.setStretchLastSection(True)

        self._default_is_apply = False
        self.name = QLabel("Name:")
        self.name_edit = QLineEdit(str(name))
        self.name_edit.setDisabled(True)

        self.color = QLabel("Color:")
        self.color_edit = QPushButton()
        #self.color_edit.setFlat(True)

        if color is not None:
            qcolor = QtGui.QColor()
            qcolor.setRgb(*color)
            #print('color =%s' % str(color))
            palette = QtGui.QPalette(self.color_edit.palette()) # make a copy of the palette
            #palette.setColor(QtGui.QPalette.Active, QtGui.QPalette.Base, \
                             #qcolor)
            palette.setColor(QtGui.QPalette.Background, QtGui.QColor('blue'))  # ButtonText
            self.color_edit.setPalette(palette)

            self.color_edit.setStyleSheet("QPushButton {"
                                          "background-color: rgb(%s, %s, %s);" % tuple(color) +
                                          #"border:1px solid rgb(255, 170, 255); "
                                          "}")

        self.representation_label = QLabel('Representation:')
        self.checkbox_wire = QCheckBox('Wireframe')
        self.checkbox_surf = QCheckBox('Surface/Solid')
        #print('representation = %s' % self.representation)
        #self.check_point = QCheckBox()

        self.use_slider = True
        self.is_opacity_edit_active = False
        self.is_opacity_edit_slider_active = False

        self.is_line_width_edit_active = False
        self.is_line_width_edit_slider_active = False

        self.is_point_size_edit_active = False
        self.is_point_size_edit_slider_active = False

        self.is_bar_scale_edit_active = False
        self.is_bar_scale_edit_slider_active = False

        self.opacity = QLabel("Opacity:")
        self.opacity_edit = QDoubleSpinBox(self)
        self.opacity_edit.setRange(0.1, 1.0)
        self.opacity_edit.setDecimals(1)
        self.opacity_edit.setSingleStep(0.1)
        self.opacity_edit.setValue(opacity)
        if self.use_slider:
            self.opacity_slider_edit = QSlider(QtCore.Qt.Horizontal)
            self.opacity_slider_edit.setRange(1, 10)
            self.opacity_slider_edit.setValue(opacity * 10)
            self.opacity_slider_edit.setTickInterval(1)
            self.opacity_slider_edit.setTickPosition(QSlider.TicksBelow)

        self.line_width = QLabel("Line Width:")
        self.line_width_edit = QSpinBox(self)
        self.line_width_edit.setRange(1, MAX_LINE_WIDTH)
        self.line_width_edit.setSingleStep(1)
        self.line_width_edit.setValue(line_width)
        if self.use_slider:
            self.line_width_slider_edit = QSlider(QtCore.Qt.Horizontal)
            self.line_width_slider_edit.setRange(1, MAX_LINE_WIDTH)
            self.line_width_slider_edit.setValue(line_width)
            self.line_width_slider_edit.setTickInterval(1)
            self.line_width_slider_edit.setTickPosition(QSlider.TicksBelow)

        if self.representation in ['point', 'surface']:
            self.line_width.setEnabled(False)
            self.line_width_edit.setEnabled(False)
            self.line_width_slider_edit.setEnabled(False)

        self.point_size = QLabel("Point Size:")
        self.point_size_edit = QSpinBox(self)
        self.point_size_edit.setRange(1, MAX_POINT_SIZE)
        self.point_size_edit.setSingleStep(1)
        self.point_size_edit.setValue(point_size)
        self.point_size.setVisible(False)
        self.point_size_edit.setVisible(False)
        if self.use_slider:
            self.point_size_slider_edit = QSlider(QtCore.Qt.Horizontal)
            self.point_size_slider_edit.setRange(1, MAX_POINT_SIZE)
            self.point_size_slider_edit.setValue(point_size)
            self.point_size_slider_edit.setTickInterval(1)
            self.point_size_slider_edit.setTickPosition(QSlider.TicksBelow)
            self.point_size_slider_edit.setVisible(False)

        if self.representation in ['wire', 'surface']:
            self.point_size.setEnabled(False)
            self.point_size_edit.setEnabled(False)
            if self.use_slider:
                self.point_size_slider_edit.setEnabled(False)

        self.bar_scale = QLabel("Bar Scale:")
        self.bar_scale_edit = QDoubleSpinBox(self)
        #self.bar_scale_edit.setRange(0.01, 1.0)  # was 0.1
        #self.bar_scale_edit.setRange(0.05, 5.0)
        self.bar_scale_edit.setDecimals(1)
        #self.bar_scale_edit.setSingleStep(bar_scale / 10.)
        self.bar_scale_edit.setSingleStep(0.1)
        self.bar_scale_edit.setValue(bar_scale)

        #if self.use_slider:
            #self.bar_scale_slider_edit = QSlider(QtCore.Qt.Horizontal)
            #self.bar_scale_slider_edit.setRange(1, 100)  # 1/0.05 = 100/5.0
            #self.bar_scale_slider_edit.setValue(opacity * 0.05)
            #self.bar_scale_slider_edit.setTickInterval(10)
            #self.bar_scale_slider_edit.setTickPosition(QSlider.TicksBelow)

        if self.representation != 'bar':
            self.bar_scale.setEnabled(False)
            self.bar_scale_edit.setEnabled(False)
            self.bar_scale.setVisible(False)
            self.bar_scale_edit.setVisible(False)
            #self.bar_scale_slider_edit.setVisible(False)
            #self.bar_scale_slider_edit.setEnabled(False)

        # show/hide
        self.checkbox_show = QCheckBox("Show")
        self.checkbox_hide = QCheckBox("Hide")
        self.checkbox_show.setChecked(show)
        self.checkbox_hide.setChecked(not show)

        if name == 'main':
            self.color.setEnabled(False)
            self.color_edit.setEnabled(False)
            self.point_size.setEnabled(False)
            self.point_size_edit.setEnabled(False)
            if self.use_slider:
                self.point_size_slider_edit.setEnabled(False)

        self.cancel_button = QPushButton("Close")

        self.create_layout()
        self.set_connections()

        if isinstance(actor_obj, CoordProperties):
            self.color_edit.hide()
            self.color.hide()
            self.opacity.hide()
            self.opacity_edit.hide()
            self.opacity_slider_edit.hide()
            self.line_width.hide()
            self.line_width_edit.hide()
            self.line_width_slider_edit.hide()

    def on_delete(self, irow):
        """deletes an actor based on the row number"""
        if irow == 0:  # main
            return
        nkeys = len(self.keys)
        if nkeys in [0, 1]:
            return
        name = self.keys[irow]
        nrows = nkeys - 1
        self.keys.pop(irow)

        header_labels = ['Groups']
        table_model = Model(self.keys, header_labels, self)
        self.table.setModel(table_model)

        if len(self.keys) == 0:
            self.update()
            self.set_as_null()
            return
        if irow == nrows:
            irow -= 1
        new_name = self.keys[irow]
        self.update_active_name(new_name)
        if self.is_gui:
            self.win_parent.delete_actor(name)

    def set_as_null(self):
        """sets the null case"""
        self.name.setVisible(False)
        self.name_edit.setVisible(False)
        self.color.setVisible(False)
        self.color_edit.setVisible(False)
        self.line_width.setVisible(False)
        self.line_width_edit.setVisible(False)
        self.point_size.setVisible(False)
        self.point_size_edit.setVisible(False)
        self.bar_scale.setVisible(False)
        self.bar_scale_edit.setVisible(False)
        self.opacity.setVisible(False)
        self.opacity_edit.setVisible(False)
        self.opacity_slider_edit.setVisible(False)
        self.point_size_slider_edit.setVisible(False)
        self.line_width_slider_edit.setVisible(False)
        self.checkbox_show.setVisible(False)
        self.checkbox_hide.setVisible(False)

    def on_update_geometry_properties_window(self, data):
        """Not Implemented"""
        return
        #new_keys = sorted(data.keys())
        #if self.active_key in new_keys:
            #i = new_keys.index(self.active_key)
        #else:
            #i = 0
        #self.table.update_data(new_keys)
        #self.out_data = data
        #self.update_active_key(i)

    def update_active_key(self, index):
        """
        Parameters
        ----------
        index : PyQt4.QtCore.QModelIndex
            the index of the list

        Internal Parameters
        -------------------
        name : str
            the name of obj
        obj : CoordProperties, AltGeometry
            the storage object for things like line_width, point_size, etc.
        """
        name = str(index.data())
            #print('name = %r' % name)
        #i = self.keys.index(self.active_key)
        self.update_active_name(name)

    def update_active_name(self, name):
        self.active_key = name
        self.name_edit.setText(name)
        obj = self.out_data[name]
        if isinstance(obj, CoordProperties):
            opacity = 1.0
            representation = 'coord'
            is_visible = obj.is_visible
        elif isinstance(obj, AltGeometry):
            line_width = obj.line_width
            point_size = obj.point_size
            bar_scale = obj.bar_scale
            opacity = obj.opacity
            representation = obj.representation
            is_visible = obj.is_visible

            self.color_edit.setStyleSheet("QPushButton {"
                                          "background-color: rgb(%s, %s, %s);" % tuple(obj.color) +
                                          #"border:1px solid rgb(255, 170, 255); "
                                          "}")
            self.allow_update = False
            self.force = False
            self.line_width_edit.setValue(line_width)
            self.point_size_edit.setValue(point_size)
            self.bar_scale_edit.setValue(bar_scale)
            self.force = True
            self.allow_update = True
        else:
            raise NotImplementedError(obj)

        #allowed_representations = [
            #'main', 'surface', 'coord', 'toggle', 'wire', 'point', 'bar']

        if self.representation != representation:
            self.representation = representation
            #if representation not in allowed_representations:
                #msg = 'name=%r; representation=%r is invalid\nrepresentations=%r' % (
                    #name, representation, allowed_representations)

            if self.representation == 'coord':
                self.color.setVisible(False)
                self.color_edit.setVisible(False)
                self.line_width.setVisible(False)
                self.line_width_edit.setVisible(False)
                self.point_size.setVisible(False)
                self.point_size_edit.setVisible(False)
                self.bar_scale.setVisible(False)
                self.bar_scale_edit.setVisible(False)
                self.opacity.setVisible(False)
                self.opacity_edit.setVisible(False)
                if self.use_slider:
                    self.opacity_slider_edit.setVisible(False)
                    self.point_size_slider_edit.setVisible(False)
                    self.line_width_slider_edit.setVisible(False)
                    #self.bar_scale_slider_edit.setVisible(False)
            else:
                self.color.setVisible(True)
                self.color_edit.setVisible(True)
                self.line_width.setVisible(True)
                self.line_width_edit.setVisible(True)
                self.point_size.setVisible(True)
                self.point_size_edit.setVisible(True)
                self.bar_scale.setVisible(True)
                #self.bar_scale_edit.setVisible(True)
                self.opacity.setVisible(True)
                self.opacity_edit.setVisible(True)
                if self.use_slider:
                    self.opacity_slider_edit.setVisible(True)
                    self.line_width_slider_edit.setVisible(True)
                    self.point_size_slider_edit.setVisible(True)
                    #self.bar_scale_slider_edit.setVisible(True)

                if name == 'main':
                    self.color.setEnabled(False)
                    self.color_edit.setEnabled(False)
                    self.point_size.setEnabled(False)
                    self.point_size_edit.setEnabled(False)
                    self.line_width.setEnabled(True)
                    self.line_width_edit.setEnabled(True)
                    self.bar_scale.setEnabled(False)
                    self.bar_scale_edit.setEnabled(False)
                    show_points = False
                    show_line_width = True
                    show_bar_scale = False
                    if self.use_slider:
                        self.line_width_slider_edit.setEnabled(True)
                        #self.bar_scale_slider_edit.setVisible(False)
                else:
                    self.color.setEnabled(True)
                    self.color_edit.setEnabled(True)

                    show_points = False
                    if self.representation in ['point', 'wire+point']:
                        show_points = True

                    show_line_width = False
                    if self.representation in ['wire', 'wire+point', 'wire+surf', 'bar', 'toggle']:
                        show_line_width = True

                    if representation == 'bar':
                        show_bar_scale = True
                    else:
                        show_bar_scale = False
                    #self.bar_scale_button.setVisible(show_bar_scale)
                    #self.bar_scale_edit.setSingleStep(bar_scale / 10.)
                    #if self.use_slider:
                        #self.bar_scale_slider_edit.setEnabled(False)

                self.point_size.setEnabled(show_points)
                self.point_size_edit.setEnabled(show_points)
                self.point_size.setVisible(show_points)
                self.point_size_edit.setVisible(show_points)

                self.line_width.setEnabled(show_line_width)
                self.line_width_edit.setEnabled(show_line_width)

                self.bar_scale.setEnabled(show_bar_scale)
                self.bar_scale_edit.setEnabled(show_bar_scale)
                self.bar_scale.setVisible(show_bar_scale)
                self.bar_scale_edit.setVisible(show_bar_scale)
                if self.use_slider:
                    self.point_size_slider_edit.setEnabled(show_points)
                    self.point_size_slider_edit.setVisible(show_points)
                    self.line_width_slider_edit.setEnabled(show_line_width)


            #if self.representation in ['wire', 'surface']:

        self.opacity_edit.setValue(opacity)
        #if self.use_slider:
            #self.opacity_slider_edit.setValue(opacity*10)
        self.checkbox_show.setChecked(is_visible)
        self.checkbox_hide.setChecked(not is_visible)

        passed = self.on_validate()
        #self.on_apply(force=True)  # TODO: was turned on...do I want this???
        #self.allow_update = True

    def create_layout(self):
        ok_cancel_box = QHBoxLayout()
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
        if self.use_slider:
            grid.addWidget(self.opacity_edit, irow, 2)
            grid.addWidget(self.opacity_slider_edit, irow, 1)
        else:
            grid.addWidget(self.opacity_edit, irow, 1)
        irow += 1

        grid.addWidget(self.line_width, irow, 0)
        if self.use_slider:
            grid.addWidget(self.line_width_edit, irow, 2)
            grid.addWidget(self.line_width_slider_edit, irow, 1)
        else:
            grid.addWidget(self.line_width_edit, irow, 1)
        irow += 1

        grid.addWidget(self.point_size, irow, 0)
        if self.use_slider:
            grid.addWidget(self.point_size_edit, irow, 2)
            grid.addWidget(self.point_size_slider_edit, irow, 1)
        else:
            grid.addWidget(self.point_size_edit, irow, 1)
        irow += 1

        grid.addWidget(self.bar_scale, irow, 0)
        if self.use_slider and 0:
            grid.addWidget(self.bar_scale_edit, irow, 2)
            grid.addWidget(self.bar_scale_slider_edit, irow, 1)
        else:
            grid.addWidget(self.bar_scale_edit, irow, 1)
        irow += 1

        wire_surf_checkboxes = QButtonGroup(self)
        wire_surf_checkboxes.addButton(self.checkbox_wire)
        wire_surf_checkboxes.addButton(self.checkbox_surf)

        checkboxs = QButtonGroup(self)
        checkboxs.addButton(self.checkbox_show)
        checkboxs.addButton(self.checkbox_hide)

        vbox = QVBoxLayout()
        vbox.addWidget(self.table, stretch=1)
        vbox.addLayout(grid)

        vbox1 = QVBoxLayout()
        vbox1.addWidget(self.checkbox_wire)
        vbox1.addWidget(self.checkbox_surf)

        vbox2 = QVBoxLayout()
        vbox2.addWidget(self.checkbox_show)
        vbox2.addWidget(self.checkbox_hide)

        #vbox.addLayout(vbox1)
        vbox.addLayout(vbox2)

        vbox.addStretch()
        #vbox.addWidget(self.check_apply)
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        self.opacity_edit.valueChanged.connect(self.on_opacity)
        self.line_width_edit.valueChanged.connect(self.on_line_width)
        self.point_size_edit.valueChanged.connect(self.on_point_size)
        self.bar_scale_edit.valueChanged.connect(self.on_bar_scale)

        if self.use_slider:
            self.opacity_slider_edit.valueChanged.connect(self.on_opacity_slider)
            self.line_width_slider_edit.valueChanged.connect(self.on_line_width_slider)
            self.point_size_slider_edit.valueChanged.connect(self.on_point_size_slider)
            #self.bar_scale_slider_edit.valueChanged.connect(self.on_bar_scale_slider)

        # self.connect(self.opacity_edit, QtCore.SIGNAL('clicked()'), self.on_opacity)
        # self.connect(self.line_width, QtCore.SIGNAL('clicked()'), self.on_line_width)
        # self.connect(self.point_size, QtCore.SIGNAL('clicked()'), self.on_point_size)

        self.color_edit.clicked.connect(self.on_color)
        self.checkbox_show.clicked.connect(self.on_show)
        self.checkbox_hide.clicked.connect(self.on_hide)
        self.cancel_button.clicked.connect(self.on_cancel)
        # closeEvent

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.close()

    def closeEvent(self, event):
        self.on_cancel()

    def on_color(self):
        """called when the user clicks on the color box"""
        name = self.active_key
        obj = self.out_data[name]
        rgb_color_ints = obj.color

        msg = name
        col = QColorDialog.getColor(QtGui.QColor(*rgb_color_ints), self, "Choose a %s color" % msg)
        if col.isValid():
            color_float = col.getRgbF()[:3]
            obj.color = color_float
            color_int = [int(colori * 255) for colori in color_float]
            self.color_edit.setStyleSheet("QPushButton {"
                                          "background-color: rgb(%s, %s, %s);" % tuple(color_int) +
                                          #"border:1px solid rgb(255, 170, 255); "
                                          "}")
        self.on_apply(force=self.force)
        #print(self.allow_update)

    def on_show(self):
        """shows the actor"""
        name = self.active_key
        is_checked = self.checkbox_show.isChecked()
        self.out_data[name].is_visible = is_checked
        self.on_apply(force=self.force)

    def on_hide(self):
        """hides the actor"""
        name = self.active_key
        is_checked = self.checkbox_hide.isChecked()
        self.out_data[name].is_visible = not is_checked
        self.on_apply(force=self.force)

    def on_line_width(self):
        """increases/decreases the wireframe (for solid bodies) or the bar thickness"""
        self.is_line_width_edit_active = True
        name = self.active_key
        line_width = self.line_width_edit.value()
        self.out_data[name].line_width = line_width
        if not self.is_line_width_edit_slider_active:
            if self.use_slider:
                self.line_width_slider_edit.setValue(line_width)
            self.is_line_width_edit_active = False
        self.on_apply(force=self.force)
        self.is_line_width_edit_active = False

    def on_line_width_slider(self):
        """increases/decreases the wireframe (for solid bodies) or the bar thickness"""
        self.is_line_width_edit_slider_active = True
        #name = self.active_key
        line_width = self.line_width_slider_edit.value()
        if not self.is_line_width_edit_active:
            self.line_width_edit.setValue(line_width)
        self.is_line_width_edit_slider_active = False

    def on_point_size(self):
        """increases/decreases the point size"""
        self.is_point_size_edit_active = True
        name = self.active_key
        point_size = self.point_size_edit.value()
        self.out_data[name].point_size = point_size
        if not self.is_point_size_edit_slider_active:
            if self.use_slider:
                self.point_size_slider_edit.setValue(point_size)
            self.is_point_size_edit_active = False
        self.on_apply(force=self.force)
        self.is_point_size_edit_active = False

    def on_point_size_slider(self):
        """increases/decreases the point size"""
        self.is_point_size_edit_slider_active = True
        #name = self.active_key
        point_size = self.point_size_slider_edit.value()
        if not self.is_point_size_edit_active:
            self.point_size_edit.setValue(point_size)
        self.is_point_size_edit_slider_active = False

    def on_bar_scale(self):
        """
        Vectors start at some xyz coordinate and can increase in length.
        Increases/decreases the length scale factor.
        """
        self.is_bar_scale_edit_active = True
        name = self.active_key
        float_bar_scale = self.bar_scale_edit.value()
        self.out_data[name].bar_scale = float_bar_scale
        if not self.is_bar_scale_edit_slider_active:
            #int_bar_scale = int(round(float_bar_scale * 20, 0))
            #if self.use_slider:
                #self.bar_scale_slider_edit.setValue(int_bar_scale)
            self.is_bar_scale_edit_active = False
        self.on_apply(force=self.force)
        self.is_bar_scale_edit_active = False

    def on_bar_scale_slider(self):
        """
        Vectors start at some xyz coordinate and can increase in length.
        Increases/decreases the length scale factor.
        """
        self.is_bar_scale_edit_slider_active = True
        #name = self.active_key
        int_bar_scale = self.bar_scale_slider_edit.value()
        if not self.is_bar_scale_edit_active:
            float_bar_scale = int_bar_scale / 20.
            self.bar_scale_edit.setValue(float_bar_scale)
        self.is_bar_scale_edit_slider_active = False

    def on_opacity(self):
        """
        opacity = 1.0 (solid/opaque)
        opacity = 0.0 (invisible)
        """
        self.is_opacity_edit_active = True
        name = self.active_key
        float_opacity = self.opacity_edit.value()
        self.out_data[name].opacity = float_opacity
        if not self.is_opacity_edit_slider_active:
            int_opacity = int(round(float_opacity * 10, 0))
            if self.use_slider:
                self.opacity_slider_edit.setValue(int_opacity)
            self.is_opacity_edit_active = False
        self.on_apply(force=self.force)
        self.is_opacity_edit_active = False

    def on_opacity_slider(self):
        """
        opacity = 1.0 (solid/opaque)
        opacity = 0.0 (invisible)
        """
        self.is_opacity_edit_slider_active = True
        #name = self.active_key
        int_opacity = self.opacity_slider_edit.value()
        if not self.is_opacity_edit_active:
            float_opacity = int_opacity / 10.
            self.opacity_edit.setValue(float_opacity)
        self.is_opacity_edit_slider_active = False

    def on_validate(self):
        self.out_data['clicked_ok'] = True
        self.out_data['clicked_cancel'] = False

        old_obj = self.out_data[self.active_key]
        old_obj.line_width = self.line_width_edit.value()
        old_obj.point_size = self.point_size_edit.value()
        old_obj.bar_scale = self.bar_scale_edit.value()
        old_obj.opacity = self.opacity_edit.value()
        #old_obj.color = self.color_edit
        old_obj.is_visible = self.checkbox_show.isChecked()
        return True
        #name_value, flag0 = self.check_name(self.name_edit)
        #ox_value, flag1 = check_float(self.transparency_edit)
        #if flag0 and flag1:
            #self.out_data['clicked_ok'] = True
            #return True
        #return False

    @property
    def is_gui(self):
        return hasattr(self.win_parent, 'on_update_geometry_properties')

    def on_apply(self, force=False):
        passed = self.on_validate()
        #print("passed=%s force=%s allow=%s" % (passed, force, self.allow_update))
        if (passed or force) and self.allow_update and self.is_gui:
            #print('obj = %s' % self.out_data[self.active_key])
            self.win_parent.on_update_geometry_properties(self.out_data, name=self.active_key)
        return passed

    def on_cancel(self):
        passed = self.on_apply(force=True)
        if passed:
            self.close()
            #self.destroy()


def main():  # pragma: no cover
    """gui independent way to test the program"""
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)


    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    parent = app
    red = (255, 0, 0)
    blue = (0, 0, 255)
    green = (0, 255, 0)
    purple = (255, 0, 255)

    # representation
    # * main - main mesh
    # * toggle - change with main mesh
    # * wire - always wireframe
    # * point - always points
    # * wire+point - point (vertex) and wireframe allowed
    # * surface - always surface
    # * bar - this can use bar scale
    data = {
        'font_size' : 8,
        'toggle' : AltGeometry(parent, 'toggle', color=green, line_width=3, opacity=0.2, representation='toggle'),
        'wire' : AltGeometry(parent, 'wire', color=purple, line_width=4, opacity=0.3, representation='wire'),
        'wire+point' : AltGeometry(parent, 'wire+point', color=blue, line_width=2, opacity=0.1, bar_scale=1.0, representation='wire+point'),
        'wire+surf' : AltGeometry(parent, 'wire+surf', display='surface', color=blue, line_width=2, opacity=0.1, bar_scale=1.0, representation='wire+surf'),
        'main' : AltGeometry(parent, 'main', color=red, line_width=1, opacity=0.0, representation='main'),
        'point' : AltGeometry(parent, 'point', color=blue, opacity=0.1, representation='point'),
        'surface' : AltGeometry(parent, 'surface', color=blue, opacity=0.1, representation='surface'),
    }
    main_window = EditGeometryProperties(data, win_parent=None)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":  # pragma: no cover
    main()
