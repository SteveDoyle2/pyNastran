"""
defines:
 - PyDialog()

"""
from pyNastran.gui.qt_version import qt_version
from qtpy.QtCore import Qt
from qtpy.QtGui import QFont
from qtpy.QtWidgets import QDialog, QComboBox

from pyNastran.bdf.utils import (
    parse_patran_syntax, parse_patran_syntax_dict)

from pyNastran.gui.utils.qt.checks.qlineedit import (
    check_path, check_save_path,
    check_int, check_positive_int_or_blank,
    check_float, check_float_ranged,
    check_name_str, check_name_length, check_format, check_format_str,
)

def make_font(font_size, is_bold=False):
    """creates a QFont"""
    font = QFont()
    font.setPointSize(font_size)
    if is_bold:
        font.setBold(is_bold)
    return font

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
        font = make_font(font_size, is_bold=False)
        self.setFont(font)

    def closeEvent(self, event):
        self.out_data['close'] = True
        event.accept()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            self.on_cancel()

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

def make_combo_box(items, initial_value):
    """
    Makes a QComboBox, sets the items, and sets an initial value.

    Parameters
    ----------
    items : List[str]
        the values of the combo box
    initial_value : str
        the value to set the combo box to

    Returns
    -------
    combo_box : QComboBox
        the pulldown
    """
    assert initial_value in items, 'initial_value=%r items=%s' % (initial_value, items)
    combo_box = QComboBox()
    combo_box.addItems(items)
    set_combo_box_text(combo_box, initial_value)

    if initial_value not in items:
        msg = 'initial_value=%r is not supported in %s' % (initial_value, items)
        raise RuntimeError(msg)
    return combo_box

def set_combo_box_text(combo_box, value):
    """sets the combo_box text"""
    if qt_version == 'pyside':
        items = [combo_box.itemText(i) for i in range(combo_box.count())]
        j = items.index(value)
        combo_box.setCurrentIndex(j)
    else:
        combo_box.setCurrentText(value)

def check_color(color_float):
    assert len(color_float) == 3, color_float
    assert isinstance(color_float[0], float), color_float
    color_int = [int(colori * 255) for colori in color_float]
    return color_float, color_int
