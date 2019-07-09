from qtpy import QtGui
from qtpy.QtWidgets import (
    QVBoxLayout, QPushButton,
    QComboBox, QLabel, QHBoxLayout, QSpinBox,
    QGridLayout, QTextEdit, QLineEdit)

import numpy as np
from pyNastran.gui.utils.qt.pydialog import PyDialog
from .modify_map import Var, TransposedVars


class ModifyMenu(PyDialog):
    def __init__(self, data, nastran_io, win_parent=None):
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
        self.update_function_name = data['update_function_name']
        self.nastran_io = nastran_io
        #--------------------------------------------
        i = 0
        grid_objs = {}
        grid = QGridLayout()
        for ivar, var in enumerate(variables):
            if isinstance(var, Var):
                i = add_scalar_var(model, var, obj, grid, grid_objs, i)
            elif isinstance(var, TransposedVars):
                i = add_transposed_vars(model, var, obj, grid, grid_objs, i)
            elif isinstance(var, list): # row-based list
                i = add_list_var(model, var, obj, grid, grid_objs, i)
            else:
                raise NotImplementedError(var)
            i += 1

        self.setWindowTitle('Modify %s' % obj.type)
        self.grid_objs = grid_objs

        #self.create_widgets()
        #self.create_layout()
        #self.set_connections()
        desc_comment = QHBoxLayout()
        self.desc = QLabel('Description:')
        self.comment = QTextEdit(obj.comment.replace('\n', '<br>'))
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
        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.close)

    def on_ok(self):
        is_valid = self.on_apply()
        if is_valid:
            self.close()

    def on_apply(self):
        obj = self.obj
        is_valid = False
        for i, var in enumerate(self.variables):
            if isinstance(var, list):
                for vari in var:
                    grid_objsi = self.grid_objs[vari.name]
                    try:
                        apply_scalar(vari, obj, grid_objsi) # list
                    except AttributeError:
                        print(var)
                        raise
            else:
                grid_objsi = self.grid_objs[var.name]
                try:
                    apply_scalar(var, obj, grid_objsi)
                except AttributeError:
                    print(var)
                    raise
        is_valid = True
        if is_valid:
            self.update_obj_in_gui(obj)
        return is_valid

    def update_obj_in_gui(self, obj):
        """updates the GUI object"""
        if self.nastran_io is not None:
            if self.update_function_name is not None:
                func = getattr(self.nastran_io, self.update_function_name)
                try:
                    func(obj)
                except Exception as e:
                    self.nastran_io.model.log.error(str(e))
                    return False
        return True

    def on_font(self, value=None):
        """update the font for the current window"""
        if value is None:
            value = self.font_size_edit.value()
        font = QtGui.QFont()
        font.setPointSize(value)
        self.setFont(font)

def cast(value_original, value_new):
    """
    Parameters
    ----------
    value_original : varies
        the original value
    value_new : str
        the new value

    Returns
    -------
    value_new_typed : varies
        the new value with the same type
    """
    if isinstance(value_original, int):
        value_new_typed = int(value_new)
    elif isinstance(value_original, float):
        value_new_typed = float(value_new)
    elif isinstance(value_original, str):
        value_new_typed = str(value_new)
    else:
        raise NotImplementedError('value_original=%r value_new=%r; type for value_original is not supported' % (
            value_original, value_new))
    return value_new_typed

def apply_scalar(var, obj, grid_objsi):
    vartype = var.vartype
    value = getattr(obj, var.var)

    if vartype == 'lineedit':
        if isinstance(value, (list, np.ndarray)):
            for j, valuei_original in enumerate(value):
                valuei_new = grid_objsi[j].text()

                #box = QLineEdit(str(valuei))
                #grid.addWidget(box, i, j+1)
                value[j] = cast(valuei_original, valuei_new)
            return
        else:
            valuei_new = grid_objsi[0].text()
            value2 = cast(value, valuei_new)

    elif vartype == 'pulldown':
        #pulldown = QComboBox()
        #objs_dict = getattr(model, var.pulldown_objs)
        #for idi, obji in sorted(objs_dict.items()):
            #pulldown.addItem('%s %i' % (obji.type, idi))
        #grid.addWidget(pulldown, i, 1)
        #grid_objsi.append(pulldown)
        return

    elif vartype == 'spinner':
        value2 = grid_objsi[0].value()
    elif vartype == 'lineedit_table':
        nrows = len(grid_objsi)
        ncols = len(grid_objsi[0])
        #value = getattr(obj, var.var)
        for irow, row in enumerate(grid_objsi):
            for jcol, valuei2 in enumerate(row):
                valuei_original = value[irow, jcol]
                valuei_new = valuei2.text()
                value[irow, jcol] = cast(valuei_original, valuei_new)
        return
    else:
        raise NotImplementedError(vartype)
    setattr(obj, var.var, value2)
    return

def add_scalar_var(model, var, obj, grid, grid_objs, i):
    label = QLabel(var.name + ':')
    vartype = var.vartype
    enabled = var.enabled
    value = getattr(obj, var.var)

    label.setEnabled(enabled)
    grid.addWidget(label, i, 0)
    if vartype == 'lineedit':
        if isinstance(value, (list, np.ndarray)):
            grid_objsi = []
            for j, valuei in enumerate(value):
                box = set_qlineedit(valuei)
                box.setEnabled(enabled)
                grid.addWidget(box, i, j+1)
                grid_objsi.append(box)
            grid_objs[var.name] = grid_objsi
            return i
        else:
            box = set_qlineedit(value)
            box.setEnabled(enabled)
    elif vartype == 'pulldown':
        pulldown_items = get_pulldown_items(model, var, value)

        if isinstance(value, (list, np.ndarray)):
            grid_objsi = []
            for j, valuei in enumerate(value):
                box, enabled = set_pulldown(model, var, valuei, enabled, pulldown_items)
                box.setEnabled(enabled)
                grid.addWidget(box, i, j+1)
                grid_objsi.append(box)
            grid_objs[var.name] = [grid_objsi]
            return i
        else:
            box, enabled = set_pulldown(model, var, value, enabled, pulldown_items)
            box.setEnabled(enabled)

    elif vartype == 'spinner':
        box = QSpinBox()
        box.setValue(value)
        box.setEnabled(enabled)
    elif vartype == 'lineedit_table':
        nrows, ncols = value.shape
        grid_objsi = np.zeros((nrows, ncols)).tolist()
        for irow, row in enumerate(value):
            for jcol, valuei in enumerate(row):
                box = set_qlineedit(valuei)
                grid.addWidget(box, i+irow+1, jcol)
                grid_objsi[irow][jcol] = box
        i += nrows
        grid_objs[var.name] = grid_objsi
        return i
    else:  # pragma: no cover
        raise NotImplementedError(vartype)
    grid.addWidget(box, i, 1)
    grid_objs[var.name] = [box]

    #self.name = name
    #self.var = var
    #self.pulldown_objs = pulldown_objs
    return i

def add_list_var(model, variables, obj, grid, grid_objs, i):
    #  combines multiple disjointed parameters into a single row
    # [I1, I2, I12, J]
    for j, var in enumerate(variables):
        grid_objsi = []
        grid_objs[var.name] = grid_objsi
        j2 = j + 1
        label = QLabel(var.name + ':')

        vartype = var.vartype
        enabled = var.enabled
        value = getattr(obj, var.var)

        label.setEnabled(enabled)
        grid.addWidget(label, i, j2)

        if vartype == 'lineedit':
            if isinstance(value, (list, np.ndarray)):
                raise NotImplementedError('list lineedit')
                #for j, valuei in enumerate(value):
                    #box = set_qlineedit(valuei)
                    #box.setEnabled(enabled)
                    #grid.addWidget(box, i, j+1)
                    #grid_objsi.append(box)
                #continue
            else:
                box = set_qlineedit(value)
                box.setEnabled(enabled)

        elif vartype == 'pulldown':
            pulldown_items = get_pulldown_items(model, var, value)
            if isinstance(value, list):
                raise NotImplementedError('pulldown (list)')
            else:
                box, enabled = set_pulldown(model, var, value, enabled, pulldown_items)
                box.setEnabled(enabled)

        elif vartype == 'spinner':
            box = QSpinBox()
            box.setValue(value)
            box.setEnabled(enabled)
        else:  # pragma: no cover
            raise NotImplementedError(vartype)
        grid.addWidget(box, i+1, j2)
        grid_objsi.append(box)

    #self.name = name
    #self.var = var
    #self.pulldown_objs = pulldown_objs
    i += 1
    return i

def add_transposed_vars(model, variables_tranposed, obj, grid, grid_objs, i):
    """
    combines different variables with the same dimensions into a single table,
    such as for the PCOMP

    material_id, thickness, theta, SOUT
    """
    variables = variables_tranposed.variables
    nvars = len(variables)
    nrows = 0
    #nrows = 5
    #  combines multiple disjointed parameters into a single row
    # [I1, I2, I12, J]
    for ivar, var in enumerate(variables):
        jcol = ivar

        grid_objsi = []
        grid_objs[var.name] = grid_objsi
        label = QLabel(var.name + ':')

        vartype = var.vartype
        enabled = var.enabled
        value = getattr(obj, var.var)

        label.setEnabled(enabled)
        grid.addWidget(label, i, jcol)

        if vartype == 'lineedit':
            if isinstance(value, (list, np.ndarray)):
                nrows = max(nrows, len(value))
                grid_objsi = []
                for irow, valuei in enumerate(value):
                    box = set_qlineedit(valuei)
                    box.setEnabled(enabled)
                    grid.addWidget(box, i+irow+2, jcol)
                    grid_objsi.append(box)
                grid_objs[var.name] = grid_objsi
            else:
                #box = set_qlineedit(value)
                #box.setEnabled(enabled)
                raise NotImplementedError('lineedit (scalar)')

        elif vartype == 'pulldown':
            pulldown_items = get_pulldown_items(model, var, value)
            if isinstance(value, list):
                grid_objsi = []
                nrows = max(nrows, len(value))
                for irow, valuei in enumerate(value):
                    box, enabled_actual = set_pulldown(model, var, valuei, enabled, pulldown_items)
                    box.setEnabled(enabled_actual)
                    grid.addWidget(box, i+irow+2, jcol)
                    grid_objsi.append(box)
                grid_objs[var.name] = grid_objsi
            else:
                raise NotImplementedError('pulldown (scalar)')
                #box, enabled_actual = set_pulldown(model, var, value, enabled, pulldown_items)
                #box.setEnabled(enabled_actual)

        #elif vartype == 'spinner':
            #box = QSpinBox()
            #box.setValue(value)
            #box.setEnabled(enabled)
        else:  # pragma: no cover
            raise NotImplementedError(vartype)
        #grid.addWidget(box, i+irow+1, jcol)
        #grid_objsi.append(box)

    #self.name = name
    #self.var = var
    #self.pulldown_objs = pulldown_objs
    #i += 1
    i += nrows + 1
    return i

def set_qlineedit(value):
    box = QLineEdit()
    if value is not None:
        box.setText(str(value))
    return box

def get_pulldown_items(model, var, value):
    pulldown_items = []
    if var.pulldown_allow_zero:
        pulldown_items.append('0')
    if value is None:
        return pulldown_items

    pulldown_type_limit = var.pulldown_type_limit
    pulldown_objs = var.pulldown_objs

    if isinstance(pulldown_objs, str):
        objs_dict = getattr(model, pulldown_objs)
        if pulldown_type_limit:
            for idi, obji in sorted(objs_dict.items()):
                card_type = obji.type
                if card_type in pulldown_type_limit:
                    pulldown_items.append('%s %i' % (card_type, idi))
        else:
            for idi, obji in sorted(objs_dict.items()):
                pulldown_items.append('%s %i' % (obji.type, idi))
    elif isinstance(pulldown_objs, list):
        # ['DBOX', 'TUBE', 'TUBE2', ...]
        pulldown_items = pulldown_objs
    else:  # pragma: no cover
        raise NotImplementedError(var)
    return pulldown_items

def set_pulldown(model, var, value, enabled, pulldown_items):
    pulldown = QComboBox()
    #print('var.pulldown_objs =', var, var.pulldown_objs)
    if value is None:
        return pulldown, False

    pulldown_objs = var.pulldown_objs
    if isinstance(pulldown_objs, str):
        objs_dict = getattr(model, pulldown_objs)
        assert isinstance(objs_dict, dict), objs_dict
        pulldown.addItems(pulldown_items)

        if var.pulldown_allow_zero and value == 0:
            name = '0'
        else:
            try:
                obji = objs_dict[value]
            except KeyError:
                print(var)
                values = list(objs_dict.keys())
                print(values)
                raise
            name = '%s %i' % (obji.type, value)

    elif isinstance(pulldown_objs, list):
        # ['DBOX', 'TUBE', 'TUBE2', ...]
        pulldown.addItems(pulldown_items)
        name = value
    else:  # pragma: no cover
        raise NotImplementedError(var)
    enabled = enabled and len(pulldown_items) > 1

    # pyside
    index = pulldown.findText(name)
    #pulldown.setItemText(name)
    pulldown.setCurrentIndex(index)
    return pulldown, enabled
