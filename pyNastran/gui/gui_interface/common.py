from __future__ import print_function
from pyNastran.bdf.utils import (
    parse_patran_syntax, parse_patran_syntax_dict, write_patran_syntax_dict)

from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4 import QtCore, QtGui
    from PyQt4.QtGui import QDialog, QLineEdit, QFocusEvent, QFont, QPushButton
elif qt_version == 5:
    from PyQt5 import QtCore, QtGui
    from PyQt5.QtGui import QFocusEvent, QFont
    from PyQt5.QtWidgets import QDialog, QLineEdit, QPushButton

elif qt_version == 'pyside':
    from PySide import QtCore, QtGui
    from PySide.QtGui import QDialog, QLineEdit, QFocusEvent, QFont, QPushButton
else:
    raise NotImplementedError('qt_version = %r' % qt_version)

class QElementEdit(QLineEdit):
    """creates a QLineEdit that can pick element ids"""
    def __init__(self, win_parent, parent=None, *args, **kwargs):
        super(QElementEdit, self).__init__(parent, *args, **kwargs)
        self.win_parent = win_parent
        #self.focusInEvent.connect(self.on_focus)

    def focusInEvent(self, event):
        self.on_focus()
        QLineEdit.focusInEvent(self, QFocusEvent(QtCore.QEvent.FocusIn))

    def on_focus_callback(self, eids, nids):
        """the callback method for ``on_focus``"""
        eids_str = write_patran_syntax_dict({'' : eids})
        self.setText(eids_str)

    def on_focus(self):
        """called when the QElementEdit is activated"""
        self.win_parent.win_parent.on_area_pick(is_eids=True, is_nids=False,
                                                callback=self.on_focus_callback,
                                                force=True)


class QPushButtonColor(QPushButton):
    """Creates a QPushButton with a face color"""
    def __init__(self, labelcolor_int):
        QPushButton.__init__(self)

        qcolor = QtGui.QColor()
        #self.color_edit.setFlat(True)
        qcolor.setRgb(*labelcolor_int)
        palette = QtGui.QPalette(self.palette())
        palette.setColor(QtGui.QPalette.Background, QtGui.QColor('blue'))
        self.setPalette(palette)
        self.setStyleSheet(
            "QPushButton {"
            "background-color: rgb(%s, %s, %s);" % tuple(labelcolor_int) +
            #"border:1px solid rgb(255, 170, 255); "
            "}")


class PyDialog(QDialog):
    """
    common class for QDialog so value checking & escape/close code
    is not repeated
    """
    def __init__(self, data, win_parent):
        QDialog.__init__(self, win_parent)
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
        if event.key() == QtCore.Qt.Key_Escape:
            self.on_cancel()

    @staticmethod
    def check_int(cell):
        text = cell.text()
        try:
            value = int(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
            return value, True
        except ValueError:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    @staticmethod
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

    @staticmethod
    def check_float(cell):
        text = cell.text()
        try:
            value = float(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
            cell.setToolTip('')
            return value, True
        except ValueError:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    @staticmethod
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
        else:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    @staticmethod
    def check_patran_syntax(cell, pound=None):
        text = str(cell.text())
        try:
            values = parse_patran_syntax(text, pound=pound)
            cell.setStyleSheet("QLineEdit{background: white;}")
            return values, True
        except ValueError as e:
            cell.setStyleSheet("QLineEdit{background: red;}")
            cell.setToolTip(str(e))
            return None, False

    @staticmethod
    def check_patran_syntax_dict(cell, pound=None):
        text = str(cell.text())
        try:
            value = parse_patran_syntax_dict(text)
            cell.setStyleSheet("QLineEdit{background: white;}")
            cell.setToolTip('')
            return value, True
        except (ValueError, SyntaxError, KeyError) as e:
            cell.setStyleSheet("QLineEdit{background: red;}")
            cell.setToolTip(str(e))
            return None, False
