"""
defines:
 - ColorDisplay
"""
# -*- coding: utf-8 -*-
from __future__ import print_function, unicode_literals

from six import text_type
from qtpy import QtGui
from qtpy.QtWidgets import (
    QWidget,
)


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
        return text_type(self.color.name())
