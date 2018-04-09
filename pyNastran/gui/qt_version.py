from __future__ import print_function
from qtpy import API as qt_version

if qt_version in ['pyqt', 'pyqt4', 'pyside']:
    qt_int = 4
elif qt_version in ['pyqt5', 'pyside2']:
    qt_int = 5
else:
    raise ImportError('PyQt4, PyQt5, or PySide is required; API=%r' % qt_version)

if qt_version not in ['pyqt', 'pyqt4', 'pyqt5', 'pyside',]: # 'pyside2'
    raise ImportError('PyQt4, PyQt5, or PySide is required; API=%r' % qt_version)

# required to make a pretty console
try:
    import pygments
    is_pygments = True
except ImportError:
    is_pygments = False
