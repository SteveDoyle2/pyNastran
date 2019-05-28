from qtpy import QtGui
from qtpy.QtWidgets import (
    QVBoxLayout, QPushButton,
    QComboBox, QLabel, QHBoxLayout, QSpinBox,
    QGridLayout, QTextEdit, QLineEdit)

import numpy as np
from pyNastran.gui.utils.qt.pydialog import PyDialog


class ModifyMenu(PyDialog):
    def __init__(self, data, win_parent=None):
        """
        Saves the data members from data and
        performs type checks
        """
        PyDialog.__init__(self, data, win_parent)

        self._updated_preference = False

        self._default_font_size = data['font_size']

        model = data['model']
        obj = data['obj']
        variables = data['variables']
        self.obj = obj
        self.variables = variables
        #--------------------------------------------
        grid = QGridLayout()
        grid_objs = []
        for i, var in enumerate(variables):
            grid_objsi = []
            grid_objs.append(grid_objsi)
            label = QLabel(var.name + ':')
            vartype = var.vartype
            value = getattr(obj, var.var)

            grid.addWidget(label, i, 0)
            if vartype == 'lineedit':
                if isinstance(value, (list, np.ndarray)):
                    for j, valuei in enumerate(value):
                        box = QLineEdit(str(valuei))
                        grid.addWidget(box, i, j+1)
                        grid_objsi.append(box)
                    continue
                else:
                    box = QLineEdit(str(value))
            elif vartype == 'pulldown':
                #box = QLineEdit(str(value))
                pulldown = QComboBox()
                objs_dict = getattr(model, var.pulldown_objs)
                for idi, obji in sorted(objs_dict.items()):
                    pulldown.addItem('%s %i' % (obji.type, idi))
                grid.addWidget(pulldown, i, 1)
                grid_objsi.append(pulldown)
                continue

            elif vartype == 'spinner':
                box = QSpinBox()
                box.setValue(value)
            grid.addWidget(box, i, 1)
            grid_objsi.append(box)

            #self.name = name
            #self.var = var
            #self.pulldown_objs = pulldown_objs

        self.setWindowTitle('Modify %s' % obj.type)
        self.grid_objs = grid_objs

        #self.create_widgets()
        #self.create_layout()
        #self.set_connections()
        desc_comment = QHBoxLayout()
        self.desc = QLabel('Description:')
        self.comment = QTextEdit(obj.comment)
        desc_comment.addWidget(self.desc)
        desc_comment.addWidget(self.comment)

        #-----------------------------------------------------------------------
        # closing
        self.apply_button = QPushButton("Apply")
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")
        #-----------------------------------------------------------------------

        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        vbox.addLayout(desc_comment)
        #vbox.addLayout(hbox)
        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

        self.grid_objs = grid_objs
        self.on_font(self._default_font_size)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.close)

    def on_ok(self):
        is_valid = self.on_apply()
        if is_valid:
            self.close()

    def on_apply(self):
        def cast(valuei, value2):
            if isinstance(valuei, int):
                valuei2 = int(value2)
            elif isinstance(valuei, float):
                valuei2 = float(value2)
            elif isinstance(valuei, str):
                valuei2 = str(value2)
            return valuei2

        obj = self.obj
        for i, var in enumerate(self.variables):
            grid_objsi = self.grid_objs[i]
            vartype = var.vartype
            value = getattr(obj, var.var)

            if vartype == 'lineedit':
                if isinstance(value, (list, np.ndarray)):
                    for j, valuei in enumerate(value):
                        valuei2 = grid_objsi[j].text()
                        #box = QLineEdit(str(valuei))
                        #grid.addWidget(box, i, j+1)
                        value[j] = cast(valuei, valuei2)
                    continue
                else:
                    valuei2 = grid_objsi[0].text()
                    value2 = cast(value, valuei2)

            elif vartype == 'pulldown':
                #pulldown = QComboBox()
                #objs_dict = getattr(model, var.pulldown_objs)
                #for idi, obji in sorted(objs_dict.items()):
                    #pulldown.addItem('%s %i' % (obji.type, idi))
                #grid.addWidget(pulldown, i, 1)
                #grid_objsi.append(pulldown)
                continue

            elif vartype == 'spinner':
                value2 = grid_objsi[0].value()
            else:
                raise NotImplementedError(vartype)
            setattr(obj, var.var, value2)
        return True

    def on_font(self, value=None):
        """update the font for the current window"""
        if value is None:
            value = self.font_size_edit.value()
        font = QtGui.QFont()
        font.setPointSize(value)
        self.setFont(font)
