# -*- coding: utf-8 -*-
from numpy import setdiff1d, unique, array, hstack

from qtpy import QtCore
from qtpy.QtWidgets import (
    QDialog, QLineEdit, QPushButton, QGridLayout, QVBoxLayout, QHBoxLayout, QApplication,
    QColorDialog, QLabel,
)
from qtpy.QtGui import QColor

from pyNastran.bdf.utils import parse_patran_syntax, parse_patran_syntax_dict
from pyNastran.bdf.cards.collpase_card import collapse_colon_packs


class ChangeBCs(QDialog):
    """
    +--------------------------+
    |     Change BCs           |
    +--------------------------+
    |                          |
    |  Property    xxx Default |
    |  Element ID  xxx Default |
    |  Angle Tol   xxx Default |
    |                          |
    |     Apply OK Cancel      |
    +--------------------------+
    """
    def __init__(self, data, win_parent=None):
        self.win_parent = win_parent
        #Init the base class
        self._default_property = data['pid']
        self._default_elements = data['eid']
        self._default_theta = data['theta']

        self.elements_pound = data['elements_pound']
        self.eids = []
        self.cids = []

        self.out_data = data
        QDialog.__init__(self, win_parent)
        #self.setupUi(self)
        self.setWindowTitle('Groups: Modify')
        self.create_widgets()
        self.create_layout()
        #self.set_connections()
        #self.show()

    def create_widgets(self):
        # Name
        self.pid = QLabel("New Property ID:")
        self.pid_edit = QLineEdit(str(self._default_property).strip())
        self.pid_button = QPushButton("Default")

        # Name
        self.elements = QLabel("Element IDs:")
        self.elements_edit = QLineEdit(str(self._default_elements).strip())
        self.elements_button = QPushButton("Default")


        # elements
        self.theta = QLabel("Theta Neighbor Max:")
        self.theta_edit = QLineEdit(str(self._default_theta).strip())
        self.theta_button = QPushButton("Default")

        # applies a unique implicitly
        self.eids = parse_patran_syntax(str(self._default_elements), pound=self.elements_pound)
        #self.cids = parse_patran_syntax(str(self._default_coords), pound=self.coords_pound)

        # continuous / discrete
        #self.checkbox_continuous = QCheckBox("Continuous")
        #self.checkbox_discrete = QCheckBox("Discrete")
        #self.checkbox_discrete.setChecked(self._default_is_discrete)

        # put these in a group
        #checkboxs2 = QButtonGroup(self)
        #checkboxs2.addButton(self.checkbox_continuous)
        #checkboxs2.addButton(self.checkbox_discrete)

        # closing
        self.apply_button = QPushButton("Apply")
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")

    def create_layout(self):
        grid = QGridLayout()
        grid.addWidget(self.pid, 0, 0)
        grid.addWidget(self.pid_edit, 0, 1)
        grid.addWidget(self.pid_button, 0, 2)

        grid.addWidget(self.elements, 1, 0)
        grid.addWidget(self.elements_edit, 1, 1)
        grid.addWidget(self.elements_button, 1, 2)

        grid.addWidget(self.theta, 2, 0)
        grid.addWidget(self.theta_edit, 2, 1)
        grid.addWidget(self.theta_button, 2, 2)

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
        """creates the actions for the menu"""
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
        col = QColorDialog.getColor(QColor(*c), self, "Choose a text color")
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
        except ValueError as error:
            cell.setStyleSheet("QLineEdit{background: red;}")
            cell.setToolTip(str(error))
            return None, False

    def check_patran_syntax_dict(self, cell, pound=None):
        text = str(cell.text())
        try:
            value = parse_patran_syntax_dict(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
            cell.setToolTip('')
            return value, True
        except (ValueError, SyntaxError, KeyError) as error:
            cell.setStyleSheet("QLineEdit{background: red;}")
            cell.setToolTip(str(error))
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
    app = QApplication(sys.argv)

    nice_blue = (0.1, 0.2, 0.4)
    #The Main window
    d = {
        'eid' : 1,
        'pid' : 2,
        'theta' : 42.,
        'name' : 'asdf',
        'coords' : 0,
        'coords_pound' : 4,
        'elements_pound' : 103,
        'elements' : '1:#',
        'color' : nice_blue,
    }
    main_window = ChangeBCs(d)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == '__main__':   # pragma: no cover
    main()
