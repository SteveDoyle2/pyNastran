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
from pyNastran.gui.gui_objects.settings import FONT_SIZE_MIN, FONT_SIZE_MAX, force_ranged
from pyNastran.gui.utils.qt.checks.qlineedit import QLINEEDIT_GOOD, QLINEEDIT_ERROR

from pyNastran.bdf.utils import (
    parse_patran_syntax, parse_patran_syntax_dict)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow
    from pyNastran.gui.typing import ColorFloat


def make_font(font_size: int, is_bold=False) -> QFont:
    """creates a QFont"""
    font = QFont()
    font_size = force_ranged(font_size, min_value=FONT_SIZE_MIN, max_value=FONT_SIZE_MAX)
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

        # we set this to 0 to indicate we haven't called on_font yet
        self.font_size = 0 # data['font_size']

    def set_font_size(self, font_size: int) -> None:
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
        cell.setStyleSheet(QLINEEDIT_GOOD)
        return values, True
    except ValueError as error:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        cell.setToolTip(str(error))
        return None, False

def check_patran_syntax_dict(cell, pound=None):
    text = str(cell.text())
    try:
        value = parse_patran_syntax_dict(text)
        cell.setStyleSheet(QLINEEDIT_GOOD)
        cell.setToolTip('')
        return value, True
    except (ValueError, SyntaxError, KeyError) as error:
        cell.setStyleSheet(QLINEEDIT_ERROR)
        cell.setToolTip(str(error))
        return None, False

def check_color(color_float: ColorFloat) -> tuple[ColorFloat, Any]:
    assert len(color_float) == 3, color_float
    assert isinstance(color_float[0], float), color_float
    color_int = [int(colori * 255) for colori in color_float]
    return color_float, color_int
