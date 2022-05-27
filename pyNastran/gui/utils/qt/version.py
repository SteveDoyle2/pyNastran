from pyNastran.gui.qt_version import qt_version #qt_int
from qtpy import QtGui

if qt_version in ['pyside6', 'pyqt6']:
    Background = QtGui.QPalette.Window
else:
    Background = QtGui.QPalette.Background
