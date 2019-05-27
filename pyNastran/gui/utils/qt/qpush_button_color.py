from qtpy import QtGui
from qtpy.QtWidgets import QPushButton


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
