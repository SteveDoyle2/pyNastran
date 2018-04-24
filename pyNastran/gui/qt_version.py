from __future__ import print_function
from qtpy import API as qt_version

if qt_version in ['pyqt', 'pyqt4']:
    qt_int = 4
    qt_version = 'pyqt4'
    from qtpy import PYQT_VERSION
elif qt_version == 'pyqt5':
    qt_int = 5
    from qtpy import PYQT_VERSION
elif qt_version == 'pyside':
    qt_int = 4
    from qtpy import PYSIDE_VERSION as PYQT_VERSION
elif qt_version == 'pyside2':
    qt_int = 5
    from qtpy import PYSIDE_VERSION as PYQT_VERSION
else:
    raise ImportError('PyQt4, PyQt5, or PySide is required; API=%r' % qt_version)

if qt_version not in ['pyqt4', 'pyqt5', 'pyside',]: # 'pyside2'
    raise ImportError('PyQt4, PyQt5, or PySide is required; API=%r' % qt_version)

# required to make a pretty console
try:
    import pygments
    is_pygments = True
except ImportError:
    is_pygments = False
