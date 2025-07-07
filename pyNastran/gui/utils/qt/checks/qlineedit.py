from __future__ import annotations
from typing import TYPE_CHECKING
import os
import numpy as np
from pyNastran.gui.utils.qt.checks.utils import (check_locale_float, is_ranged_value,
                                                 check_format_str)
if TYPE_CHECKING:  # pragma: no cover
    from qtpy.QtWidgets import QLineEdit, QSpinBox, QDoubleSpinBox

QLINEEDIT_ERROR = "QLineEdit{background: red;}"
QLINEEDIT_GOOD = "QLineEdit{background: white;}"
QLINEEDIT_DISABLED = "QLineEdit{background: lightgray;}"

QTEXTEDIT_ERROR = "QTextEdit{background: red;}"
QTEXTEDIT_GOOD = "QTextEdit{background: white;}"
QTEXTEDIT_DISABLED = "QTextEdit{background: lightgray;}"

def check_path(cell: QLineEdit) -> tuple[str, bool]:
    """verifies that the path exists"""
    path, passed = check_name_str(cell)
    if not passed:
        return '', False

    if not os.path.exists(path):
        cell.setStyleSheet(QLINEEDIT_ERROR)
        return '', False
    cell.setStyleSheet(QLINEEDIT_GOOD)
    return path, True

def check_save_path(cell: QLineEdit) -> tuple[str, bool]:
    """verifies that the path is saveable..."""
    text, passed = check_name_str(cell)
    if not passed:
        return '', False
    return text, passed

#-------------------------------------------------------------------------------
def check_int(cell: QLineEdit | QSpinBox) -> tuple[int, bool]:
    """
    Colors the cell red if the integer is invalid

    Parameters
    ----------
    cell : QLineEdit()
        a PyQt/PySide object

    Returns
    -------
    value : int / None
        int : the value as a int
        None : is_passed=False
    is_passed : bool
        is this a valid integer

    """
    text = cell.text()
    try:
        value = int(text)
        cell.setStyleSheet(QLINEEDIT_GOOD)
        return value, True
    except ValueError:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        return 0, False

def check_positive_int_or_blank(cell: QLineEdit) -> tuple[int, bool]:
    text = str(cell.text()).strip()
    if len(text) == 0:
        return -1, True
    try:
        value = int(text)
    except ValueError:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        return -1, False

    if value < 1:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        return -1, False

    cell.setStyleSheet(QLINEEDIT_GOOD)
    return value, True

#def check_float(cell):
    #text = cell.text()
    #try:
        #value = eval_float_from_string(text)
        #cell.setStyleSheet(QLINEEDIT_GOOD)
        #return value, True
    #except ValueError:
        #cell.setStyleSheet(QLINEEDIT_ERROR)
        #return None, False

def check_float(cell: QLineEdit | QSpinBox | QDoubleSpinBox) -> tuple[float, bool]:
    """
    Colors the cell red if the float is invalid

    Parameters
    ----------
    cell : QLineEdit()
        a PyQt/PySide object

    Returns
    -------
    value : float / None
        float : the value as a float
        None : is_passed=False
    is_passed : bool
        is this a valid float

    # Examples
    >>> cell = QLineEdit('3.14')
    >>> value, is_passed = check_float(cell)
    # value=3.14, is_passed=True

    >>> cell = QLineEdit('cat')
    >>> value, is_passed = check_float(cell)
    # value=0.0, is_passed=False


    """
    text = cell.text()
    value, is_valid = check_locale_float(text)
    if is_valid:
        cell.setStyleSheet(QLINEEDIT_GOOD)
        return value, True
    else:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        return np.nan, False

def check_float_ranged(cell: QLineEdit,
                       min_value=None, max_value=None,
                       min_inclusive: bool=True,
                       max_inclusive: bool=True) -> tuple[float, bool]:
    """
    Colors the cell red if the float is invalid or the value is outside
    the range [min_value, max_value].

    Parameters
    ----------
    cell : QLineEdit()
        a PyQt/PySide object
    min_value / max_value : float / None
        float : the constraint is active
        None : the constraint is inactive
    min_inclusive / max_inclusive; bool; default=True
        flips [min_value, max_value] to:
          - (min_value, max_value)
          - [min_value, max_value)
          - (min_value, max_value]

    Returns
    -------
    value : float
        float : the value as a float
        nan : is_passed=False
    is_passed : bool
        is this a valid float that meets the range requirements
    """
    value, is_passed = check_float(cell)
    if not is_passed:
        #print("failed %r" % value)
        return value, is_passed

    is_ranged = is_ranged_value(
        value, min_value=min_value, max_value=max_value,
        min_inclusive=min_inclusive, max_inclusive=max_inclusive)

    style = QLINEEDIT_GOOD
    if not is_ranged:
        value = np.nan
        style = QLINEEDIT_ERROR
    cell.setStyleSheet(style)

    return value, is_ranged

#-------------------------------------------------------------------------------
def check_name_str(cell: QLineEdit) -> tuple[str, bool]:
    """
    Verifies that the data is string-able.

    Parameters
    ----------
    cell : QLineEdit
        a QLineEdit containing a string.

    """
    cell_value = cell.text()
    try:
        text = str(cell_value).strip()
    except UnicodeEncodeError:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        return '', False

    if len(text):
        cell.setStyleSheet(QLINEEDIT_GOOD)
        return text, True
    else:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        return '', False

def check_name_length(cell: QLineEdit) -> tuple[str, bool]:
    """
    Verifies that the string has at least 1 non-whitespace character.

    Parameters
    ----------
    cell : QLineEdit
        a QLineEdit containing a string.

    """
    cell_value = cell.text()
    text = cell_value.strip()

    if len(text):
        cell.setStyleSheet(QLINEEDIT_GOOD)
        return text, True
    else:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        return '', False

def check_format(cell: QLineEdit) -> tuple[str, bool]:
    """
    Checks a QLineEdit string formatter

    Parameters
    ----------
    cell : QLineEdit
        a QLineEdit containing a string formatter like:
        {'%s', '%i', '%d', %f', '%g', '%.3f', '%e'}

    Returns
    -------
    text : str / None
        str : the validated text of the QLineEdit
        None : the format is invalid
    is_valid : bool
        The str/None flag to indicate if the string formatter is valid
    """
    text = str(cell.text())
    text2, is_valid = check_format_str(text)

    if is_valid:
        cell.setStyleSheet(QLINEEDIT_GOOD)
        return text2, True
    cell.setStyleSheet(QLINEEDIT_ERROR)
    return '', False
