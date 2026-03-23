from typing import Optional, Any
from qtpy.QtWidgets import (
    QWidget, QComboBox, QLineEdit, QCheckBox,
    # QCheckBox,
    QGridLayout,
)

def create_grid_from_list(parent,
                          mylist: list[tuple[QWidget, ...]]) -> QGridLayout:
    grid = QGridLayout(parent)
    for irow, row in enumerate(mylist):
        for jcol, obj in enumerate(row):
            if isinstance(obj, str) and str == '':
                pass
            else:
                assert not isinstance(obj, str), obj
                grid.addWidget(obj, irow, jcol)
    return grid

def load_lineedits(data: dict[str, Any],
                   line_edits: list[tuple[str, int, QLineEdit]]):
    for key, index, line_edit in line_edits:
        if key not in data:
            # print(f'apply_settings: skipping key={key!r}')
            continue
        values = data[key]
        if index != -1:
            value = values[index]
        else:
            value = values

        str_value = _to_str(value)

        # print('type(value) =', type(value))
        # print(f'{key+":":<10} values={values}[{index!r}]={value!r} -> {str_value!r}')
        try:
            line_edit.setText(str_value)
        except AttributeError:  # pragma: no cover
            print(key)
            raise


def load_checkboxs(data: dict[str, Any],
                   checkboxs: list[tuple[str, QCheckBox]]) -> None:
    """
    checkboxs = [('use_rhoref', self.use_rhoref_checkbox),]
    """
    # attrs aren't stored
    for (key, checkbox) in checkboxs:
        if key not in data:
            continue
        val = data[key]
        assert isinstance(val, bool), (key, val)
        try:
            checkbox.setChecked(val)
        except AttributeError:  # pragma: no cover
            print(key)
            raise


def load_pulldowns(data: dict[str, Any],
                   pulldown_edits: list[tuple[str, QComboBox, list[str]]]) -> None:
    for key, pulldown_edit, values in pulldown_edits:
        if key not in data:
            # print(f'apply_settings: skipping key={key!r}')
            continue
        value = data[key]
        # print(f'key={key} values={values} value={value!r}')
        index = values.index(value)
        pulldown_edit.setCurrentIndex(index)


def load_min_max_lineedits(data: dict[str, Any],
                           min_max_line_edits: list[tuple[str, QLineEdit, QLineEdit]],
                           ) -> None:
    for key, line_edit_min, line_edit_max in min_max_line_edits:
        if key not in data:
            # print(f'apply_settings: skipping key={key!r}')
            continue
        values = data[key]
        value0 = _to_str(values[0])
        value1 = _to_str(values[1])
        line_edit_min.setText(value0)
        line_edit_max.setText(value1)

def _to_str(value: Optional[int | float]) -> str:
    if value is None:
        str_value = ''
    else:
        str_value = str(value)
    return str_value
