"""
defines:
 - PyDialog()

"""
from __future__ import annotations
from typing import TYPE_CHECKING
#from pyNastran.gui.qt_version import qt_version
from qtpy.QtCore import Qt
from qtpy.QtGui import QFont, QIntValidator, QDoubleValidator
from qtpy.QtWidgets import QDialog, QLineEdit

from pyNastran.bdf.utils import (
    parse_patran_syntax, parse_patran_syntax_dict)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow


def make_font(font_size, is_bold=False):
    """creates a QFont"""
    font = QFont()
    font.setPointSize(font_size)
    if is_bold:
        font.setBold(is_bold)
    return font

class QIntEdit(QLineEdit):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        qfloat_validator = QIntValidator()
        self.setValidator(qfloat_validator)

class QFloatEdit(QLineEdit):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        qfloat_validator = QDoubleValidator()
        self.setValidator(qfloat_validator)


class PyDialog(QDialog):
    """
    common class for QDialog so value checking & escape/close code
    is not repeated
    """
    def __init__(self, data, win_parent: MainWindow):
        super(PyDialog, self).__init__(win_parent)
        self.out_data = data
        self.win_parent = win_parent
        self.font_size = data['font_size']

    def set_font_size(self, font_size):
        """
        Updates the font size of all objects in the PyDialog

        Parameters
        ----------
        font_size : int
            the font size
        """
        if self.font_size == font_size:
            return
        self.font_size = font_size
        font = make_font(font_size, is_bold=False)
        self.setFont(font)

    def closeEvent(self, event):
        self.out_data['close'] = True
        event.accept()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            if hasattr(self, 'on_cancel'):
                self.on_cancel()
            else:
                self.closeEvent(event)

def check_patran_syntax(cell, pound=None):
    text = str(cell.text())
    try:
        values = parse_patran_syntax(text, pound=pound)
        cell.setStyleSheet("QLineEdit{background: white;}")
        return values, True
    except ValueError as error:
        cell.setStyleSheet("QLineEdit{background: red;}")
        cell.setToolTip(str(error))
        return None, False

def check_patran_syntax_dict(cell, pound=None):
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

def check_color(color_float):
    assert len(color_float) == 3, color_float
    assert isinstance(color_float[0], float), color_float
    color_int = [int(colori * 255) for colori in color_float]
    return color_float, color_int
