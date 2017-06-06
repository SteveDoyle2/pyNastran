"""
defines:
 - ColorDisplay
"""
# -*- coding: utf-8 -*-
from __future__ import print_function, unicode_literals

from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4 import QtGui
    from PyQt4.QtGui import (
        QWidget,
        #QLabel, QLineEdit, QPushButton, QWidget, QApplication,
        #QListWidget, QGridLayout, QHBoxLayout, QVBoxLayout,
    )
elif qt_version == 5:
    from PyQt5 import QtGui
    from PyQt5.QtWidgets import (
        QWidget,
        #QLabel, QLineEdit, QPushButton, QApplication,
        #QListWidget, QGridLayout, QHBoxLayout, QVBoxLayout,
    )
elif qt_version == 'pyside':
    from PySide import QtGui
    #from PySide.QtGui import (
        #QLabel, QLineEdit, QPushButton, QWidget, QApplication,
        #QListWidget, QGridLayout, QHBoxLayout, QVBoxLayout,
    #)
else:
    raise NotImplementedError('qt_version = %r' % qt_version)


class ColorDisplay(QWidget):
    """
    http://stackoverflow.com/questions/4624985/how-simply-display-a-qcolor-using-pyqt
    """
    def __init__(self, parent, default_color=None):
        super(ColorDisplay, self).__init__(parent)
        self.color = default_color
        self.setColor(self.color)

    def setColor(self, color):
        if color is not None:
            color = [int(255 * i) for i in color]
        self.color = QtGui.QColor(*color)
        self.update()

    def paintEvent(self, event=None):
        painter = QtGui.QPainter(self)
        if self.color is not None:
            painter.setBrush(QtGui.QBrush(self.color))
            painter.drawRect(self.rect())

    def getColorName(self):
        return unicode(self.color.name())
