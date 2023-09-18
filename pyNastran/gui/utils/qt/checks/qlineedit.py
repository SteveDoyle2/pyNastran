from __future__ import annotations
from typing import TYPE_CHECKING
import os
from pyNastran.gui.utils.qt.checks.utils import (check_locale_float, is_ranged_value,
                                                 check_format_str)
if TYPE_CHECKING:
    from qtpy.QtWidgets import QLineEdit


QLINE_EDIT_BASIC = 'QLineEdit{background: white;}'
QLINE_EDIT_ERROR = 'QLineEdit{background: red;}'

def check_path(cell: QLineEdit) -> tuple[str, bool]:
    """verifies that the path exists"""
    path, passed = check_name_str(cell)
    if not passed:
        return None, False

    if not os.path.exists(path):
        cell.setStyleSheet(QLINE_EDIT_ERROR)
        return None, False
    cell.setStyleSheet(QLINE_EDIT_BASIC)
    return path, True

def check_save_path(cell: QLineEdit) -> tuple[str, bool]:
    """verifies that the path is saveable..."""
    text, passed = check_name_str(cell)
    if not passed:
        return None, False
    return text, passed

#-------------------------------------------------------------------------------
def check_int(cell: QLineEdit) -> tuple[int, bool]:
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
        cell.setStyleSheet("QLineEdit{background: white;}")
        return value, True
    except ValueError:
        cell.setStyleSheet("QLineEdit{background: red;}")
        return None, False

def check_positive_int_or_blank(cell: QLineEdit) -> tuple[int, bool]:
    text = str(cell.text()).strip()
    if len(text) == 0:
        return None, True
    try:
        value = int(text)
    except ValueError:
        cell.setStyleSheet("QLineEdit{background: red;}")
        return None, False

    if value < 1:
        cell.setStyleSheet("QLineEdit{background: red;}")
        return None, False

    cell.setStyleSheet("QLineEdit{background: white;}")
    return value, True

#def check_float(cell):
    #text = cell.text()
    #try:
        #value = eval_float_from_string(text)
        #cell.setStyleSheet("QLineEdit{background: white;}")
        #return value, True
    #except ValueError:
        #cell.setStyleSheet("QLineEdit{background: red;}")
        #return None, False

def check_float(cell: QLineEdit) -> tuple[float, bool]:
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
    # value=None, is_passed=False


    """
    text = cell.text()
    value, is_valid = check_locale_float(text)
    if is_valid:
        cell.setStyleSheet("QLineEdit{background: white;}")
        return value, True
    else:
        cell.setStyleSheet("QLineEdit{background: red;}")
        return None, False

def check_float_ranged(cell: QLineEdit,
                       min_value=None, max_value=None,
                       min_inclusive=True, max_inclusive=True) -> tuple[str, bool]:
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
    value : float / None
        float : the value as a float
        None : is_passed=False
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

    color = 'white'
    if not is_ranged:
        value = None
        color = 'red'
    cell.setStyleSheet("QLineEdit{background: %s;}" % color)

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
        cell.setStyleSheet("QLineEdit{background: red;}")
        return None, False

    if len(text):
        cell.setStyleSheet("QLineEdit{background: white;}")
        return text, True
    else:
        cell.setStyleSheet("QLineEdit{background: red;}")
        return None, False

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
        cell.setStyleSheet("QLineEdit{background: white;}")
        return text, True
    else:
        cell.setStyleSheet("QLineEdit{background: red;}")
        return None, False

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
        cell.setStyleSheet("QLineEdit{background: white;}")
        return text2, True
    cell.setStyleSheet("QLineEdit{background: red;}")
    return None, False
