# -*- coding: utf-8 -*-
from __future__ import print_function, unicode_literals
from PyQt4 import QtCore, QtGui

from numpy import setdiff1d, unique, array, hstack

from pyNastran.bdf.utils import parse_patran_syntax, parse_patran_syntax_dict
from pyNastran.bdf.cards.base_card import collapse_colon_packs


class ColorDisplay(QtGui.QWidget):
    """
    http://stackoverflow.com/questions/4624985/how-simply-display-a-qcolor-using-pyqt
    """
    def __init__(self, parent, default_color=None):
        super(ColorDisplay, self).__init__(parent)
        #if default_color is None:
        #default_color = 'red'
        #print("default color", default_color)
        self.color = default_color
        #self.color = QtGui.QColor(*default_color)
        #self.color = (0.1, 0.2, 0.3)
        #self.color = None
        #self.color = 'red'
        #self.color = QtGui.QColor((0.1, 0.2, 0.3))
        #self.update()
        self.setColor(self.color)

    def setColor(self, color):
        #print('setColor -> %s' % str(color))
        if color is not None:
            color = [int(255 * i) for i in color]
        #print('qcolor input -> %s' % str(color))
        self.color = QtGui.QColor(*color)
        self.update()

    def paintEvent(self, event=None):
        painter = QtGui.QPainter(self)
        if self.color is not None:
            painter.setBrush(QtGui.QBrush(self.color))
            painter.drawRect(self.rect())

    def getColorName(self):
        return unicode(self.color.name())


class GroupsModify(QtGui.QDialog):
    """
    +--------------------------+
    |     Groups : Modify      |
    +--------------------------+
    |                          |
    |  Name        xxx Default |
    |  Coords      xxx Default |
    |  Elements    xxx Default |
    |  Color       xxx Default |
    |  Add         xxx Add     |
    |  Remove      xxx Remove  |
    |                          |
    |     Apply OK Cancel      |
    +--------------------------+
    """
    def __init__(self, data, win_parent=None):
        self.win_parent = win_parent
        #Init the base class
        self._default_name = data['name']
        self._default_coords = data['coords']
        self._default_elements = data['elements']
        self._default_color = data['color']

        self.coords_pound = data['coords_pound']
        self.elements_pound = data['elements_pound']

        #self._default_is_discrete = data['is_discrete']

        self.out_data = data

        QtGui.QDialog.__init__(self, win_parent)
        #self.setupUi(self)
        self.setWindowTitle('Groups: Modify')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.show()

    def create_widgets(self):
        # Name
        self.name = QtGui.QLabel("Name:")
        self.name_edit = QtGui.QLineEdit(str(self._default_name).strip())
        self.name_button = QtGui.QPushButton("Default")

        # Name
        self.coords = QtGui.QLabel("Coord IDs:")
        self.coords_edit = QtGui.QLineEdit(str(self._default_coords).strip())
        self.coords_button = QtGui.QPushButton("Default")


        # elements
        self.elements = QtGui.QLabel("Element IDs:")
        self.elements_edit = QtGui.QLineEdit(str(self._default_elements).strip())
        self.elements_button = QtGui.QPushButton("Default")

        # add
        self.add = QtGui.QLabel("Add:")
        self.add_edit = QtGui.QLineEdit(str(''))
        self.add_button = QtGui.QPushButton("Add")

        # remove
        self.remove = QtGui.QLabel("Remove:")
        self.remove_edit = QtGui.QLineEdit(str(''))
        self.remove_button = QtGui.QPushButton("Remove")

        # color
        self.color = QtGui.QLabel("Color:")
        #self.color_edit = QtGui.QPushButton()
        #self.color_edit.setColor(self._default_color)
        self.color_edit = ColorDisplay(self, self._default_color)
        self.color_button = QtGui.QPushButton("Default")

        # applies a unique implicitly
        self.eids = parse_patran_syntax(str(self._default_elements), pound=self.elements_pound)
        self.cids = parse_patran_syntax(str(self._default_coords), pound=self.coords_pound)

        # continuous / discrete
        #self.checkbox_continuous = QtGui.QCheckBox("Continuous")
        #self.checkbox_discrete = QtGui.QCheckBox("Discrete")
        #self.checkbox_discrete.setChecked(self._default_is_discrete)

        # put these in a group
        #checkboxs2 = QtGui.QButtonGroup(self)
        #checkboxs2.addButton(self.checkbox_continuous)
        #checkboxs2.addButton(self.checkbox_discrete)

        # closing
        self.apply_button = QtGui.QPushButton("Apply")
        self.ok_button = QtGui.QPushButton("OK")
        self.cancel_button = QtGui.QPushButton("Cancel")

    def create_layout(self):
        grid = QtGui.QGridLayout()
        grid.addWidget(self.name, 0, 0)
        grid.addWidget(self.name_edit, 0, 1)
        grid.addWidget(self.name_button, 0, 2)

        grid.addWidget(self.coords, 1, 0)
        grid.addWidget(self.coords_edit, 1, 1)
        grid.addWidget(self.coords_button, 1, 2)

        grid.addWidget(self.elements, 2, 0)
        grid.addWidget(self.elements_edit, 2, 1)
        grid.addWidget(self.elements_button, 2, 2)

        grid.addWidget(self.color, 3, 0)
        grid.addWidget(self.color_edit, 3, 1)
        grid.addWidget(self.color_button, 3, 2)

        grid.addWidget(self.add, 4, 0)
        grid.addWidget(self.add_edit, 4, 1)
        grid.addWidget(self.add_button, 4, 2)

        grid.addWidget(self.remove, 5, 0)
        grid.addWidget(self.remove_edit, 5, 1)
        grid.addWidget(self.remove_button, 5, 2)


        ok_cancel_box = QtGui.QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)


        vbox = QtGui.QVBoxLayout()
        vbox.addLayout(grid)
        #vbox.addLayout(checkboxes)
        #vbox.addLayout(grid2)
        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)

        self.setLayout(vbox)

    def set_connections(self):
        self.connect(self.name_button, QtCore.SIGNAL('clicked()'), self.on_default_name)
        self.connect(self.coords_button, QtCore.SIGNAL('clicked()'), self.on_default_coords)
        self.connect(self.elements_button, QtCore.SIGNAL('clicked()'), self.on_default_elements)

        self.connect(self.add_button, QtCore.SIGNAL('clicked()'), self.on_add)
        self.connect(self.remove_button, QtCore.SIGNAL('clicked()'), self.on_remove)

        self.connect(self.color_edit, QtCore.SIGNAL('clicked()'), self.on_edit_color)
        #self.color_edit.clicked.connect(self.on_edit_color)


        self.connect(self.color_button, QtCore.SIGNAL('clicked()'), self.on_default_color)
        self.connect(self.apply_button, QtCore.SIGNAL('clicked()'), self.on_apply)
        self.connect(self.ok_button, QtCore.SIGNAL('clicked()'), self.on_ok)
        self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self.on_cancel)

    def closeEvent(self, event):
        event.accept()

    def on_add(self):
        adict, is_valid = self.check_patran_syntax_dict(self.add_edit)
        if not is_valid:
            #self.add_edit.setStyleSheet("QLineEdit{background: red;}")
            return

        self.eids = _add(adict, ['e', 'elem', 'element'], self.eids)
        self.cids = _add(adict, ['c', 'cid', 'coord'], self.cids)
        self._apply_cids_eids()

        self.add_edit.clear()
        self.add_edit.setStyleSheet("QLineEdit{background: white;}")

    def _apply_cids_eids(self):
        ctext = _get_collapsed_text(self.cids)
        etext = _get_collapsed_text(self.eids)

        self.coords_edit.setText(str(ctext.lstrip()))
        self.elements_edit.setText(str(etext.lstrip()))

    def on_remove(self):
        adict, is_valid = self.check_patran_syntax_dict(self.remove_edit)
        if not is_valid:
            #self.remove_edit.setStyleSheet("QLineEdit{background: red;}")
            return

        self.eids = _remove(adict, ['e', 'elem', 'element'], self.eids)
        self.cids = _remove(adict, ['c', 'cid', 'coord'], self.cids)
        self._apply_cids_eids()

        self.remove_edit.clear()
        self.remove_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_name(self):
        self.name_edit.setText(str(self._default_name))
        self.name_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_coords(self):
        self.coords_edit.setText(str(self._default_coords))
        self.coords_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_elements(self):
        self.elements_edit.setText(str(self._default_elements))
        self.elements_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_edit_color(self):
        c = [int(255 * i) for i in self.text_col]
        #print('c =', c)
        col = QtGui.QColorDialog.getColor(QtGui.QColor(*c), self, "Choose a text color")
        self.color.SetColor(col)

    def on_default_color(self):
        self.color_edit.setColor(self._default_color)
        #self.elements_edit.setStyleSheet("QLineEdit{background: white;}")

    def check_patran_syntax(self, cell, pound=None):
        text = str(cell.text())
        try:
            value = parse_patran_syntax(text, pound=pound)
            cell.setStyleSheet("QLineEdit{background: white;}")
            return value, True
        except ValueError as e:
            cell.setStyleSheet("QLineEdit{background: red;}")
            cell.setToolTip(str(e))
            return None, False

    def check_patran_syntax_dict(self, cell, pound=None):
        text = str(cell.text())
        try:
            value = parse_patran_syntax_dict(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
            cell.setToolTip('')
            return value, True
        except (ValueError, SyntaxError, KeyError) as e:
            cell.setStyleSheet("QLineEdit{background: red;}")
            cell.setToolTip(str(e))
            return None, False

    def check_float(self, cell):
        text = cell.text()
        try:
            value = float(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
            cell.setToolTip('')
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
        coords_value, flag1 = self.check_patran_syntax(self.coords_edit,
                                                       pound=self.coords_pound)
        elements_value, flag2 = self.check_patran_syntax(self.elements_edit,
                                                         pound=self.elements_pound)
        #color = self.color

        if flag0 and flag1 and flag2:
            self.out_data['name'] = name_value
            self.out_data['elements'] = elements_value
            self.out_data['coords'] = coords_value
            self.out_data['clicked_ok'] = True

            #print("name = %r" % self.name_edit.text())
            #print("min = %r" % self.min_edit.text())
            #print("max = %r" % self.max_edit.text())
            #print("format = %r" % self.format_edit.text())
            return True
        return False

    def on_apply(self):
        passed = self.on_validate()
        if passed:
            self.win_parent.on_modify_group(self.out_data)

    def on_ok(self):
        passed = self.on_validate()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.close()

def _get_collapsed_text(values):
    singles, doubles = collapse_colon_packs(values)
    text = ' '.join([str(s) for s in singles]) + ' '
    text += ' '.join([''.join([str(doublei) for doublei in double]) for double in doubles])
    return text

def _add(adict, keys, values_to_add):
    value_stack = []
    for key in keys:
        if key not in adict:
            continue
        values = adict[key]
        value_stack.append(values)
    if value_stack:
        value_stack.append(values_to_add)
        values_add = unique(hstack(value_stack))
        return values_add
    return values_to_add

def _remove(adict, keys, values_to_remove):
    value_stack = []
    for key in keys:
        if key not in adict:
            continue
        value_stack.append(adict[key])
    if value_stack:
        values_remove = unique(hstack(value_stack))
        return setdiff1d(values_to_remove, values_remove)
    return values_to_remove

def main():
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QtGui.QApplication(sys.argv)

    nice_blue = (0.1, 0.2, 0.4)
    #The Main window
    d = {
        'name' : 'asdf',
        'coords' : 0,
        'coords_pound' : 4,
        'elements_pound' : 103,
        'elements' : '1:#',
        'color' : nice_blue,
    }
    main_window = GroupsModify(d)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
