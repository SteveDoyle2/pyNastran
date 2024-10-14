from typing import Union
from qtpy import QtGui
from qtpy.QtWidgets import QPushButton
from .version import Background
from pyNastran.gui.typing import ColorFloat, ColorInt

class QPushButtonColor(QPushButton):
    """Creates a QPushButton with a face color"""
    def __init__(self, labelcolor_int):
        QPushButton.__init__(self)

        qcolor = QtGui.QColor()
        #self.color_edit.setFlat(True)
        qcolor.setRgb(*labelcolor_int)
        palette = QtGui.QPalette(self.palette())
        palette.setColor(Background, QtGui.QColor('blue'))
        self.setPalette(palette)
        self.set_color(labelcolor_int)
        #self.setStyleSheet(
            #"QPushButton {"
            #"background-color: rgb(%s, %s, %s);" % tuple(labelcolor_int) +
            ##"border:1px solid rgb(255, 170, 255); "
            #"}")

    def set_color(self, color: ColorInt | ColorFloat) -> None:
        assert isinstance(color, (list, tuple)), color
        assert len(color) == 3, color

        if isinstance(color[0], int):
            color_int = color
        else:
            assert isinstance(color[0], float), color
            color_int = [int(255*colori) for colori in color]

        self.setStyleSheet(
            "QPushButton {"
            "background-color: rgb(%s, %s, %s);" % tuple(color_int) +
            #"border:1px solid rgb(255, 170, 255); "
            "}")
