"""
defines:
 - PyDialog()
"""
from __future__ import print_function

from qtpy.QtCore import Qt
from qtpy.QtGui import QFont
from qtpy.QtWidgets import QDialog

from pyNastran.bdf.utils import (
    parse_patran_syntax, parse_patran_syntax_dict)


class PyDialog(QDialog):
    """
    common class for QDialog so value checking & escape/close code
    is not repeated
    """
    def __init__(self, data, win_parent):
        super(PyDialog, self).__init__(win_parent)
        self.out_data = data
        self.win_parent = win_parent
        self.font_size = None

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
        font = QFont()
        font.setPointSize(font_size)
        self.setFont(font)

    def closeEvent(self, event):
        self.out_data['close'] = True
        event.accept()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            self.on_cancel()

def check_positive_int_or_blank(cell):
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

def check_int(cell):
    text = cell.text()
    try:
        value = int(text)
        cell.setStyleSheet("QLineEdit{background: white;}")
        return value, True
    except ValueError:
        cell.setStyleSheet("QLineEdit{background: red;}")
        return None, False

#def check_float(cell):
    #text = cell.text()
    #try:
        #value = eval_float_from_string(text)
        #cell.setStyleSheet("QLineEdit{background: white;}")
        #return value, True
    #except ValueError:
        #cell.setStyleSheet("QLineEdit{background: red;}")
        #return None, False


def check_float(cell):
    text = cell.text()
    try:
        value = float(text)
        cell.setStyleSheet("QLineEdit{background: white;}")
        return value, True
    except ValueError:
        cell.setStyleSheet("QLineEdit{background: red;}")
        return None, False

def check_format(cell):
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
    cell.setStyleSheet("QLineEdit{background: red;}")
    return None, False
