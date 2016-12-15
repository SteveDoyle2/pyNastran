
from pyNastran.bdf.utils import parse_patran_syntax, parse_patran_syntax_dict

from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4.QtGui import QDialog
elif qt_version == 5:
    from PyQt5.QtWidgets import QDialog
elif qt_version == 'pyside':
    from PySide.QtGui import QDialog
else:
    raise NotImplementedError('qt_version = %r' % qt_version)


class PyDialog(QDialog):

    def __init__(self, data, win_parent):
        QDialog.__init__(self, win_parent)
        self.out_data = data

    def closeEvent(self, event):
        self.out_data['close'] = True
        event.accept()

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.on_cancel()

    def check_float(self, cell):
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
    def check_name(cell):
        text = str(cell.text()).strip()
        if len(text):
            cell.setStyleSheet("QLineEdit{background: white;}")
            return text, True
        else:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

        if self._default_name != text:
            if self._default_name in self.out_data:
                cell.setStyleSheet("QLineEdit{background: white;}")
                return text, True
            else:
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
