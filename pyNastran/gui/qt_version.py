"""
Figures out the "optimal" Qt version to use in a way that:

 - uses the qt version specified by the QT_API environment variable
 - picks up the already imported version
 - selects PySide, PySide2, PyQt5, PyQt4 (in that order)

"""
from __future__ import print_function
import os
import sys

API = os.environ.get('QT_API', '').lower()
if API:
    from qtpy import API as qt_version
elif 'PyQt5' in sys.modules:
    qt_version = 'pyqt5'
elif 'PySide2' in sys.modules:
    qt_version = 'pyside2'
elif 'PyQt4' in sys.modules:
    qt_version = 'pyqt4'
elif 'PySide' in sys.modules:
    qt_version = 'pyside'
else:
    found_gui = False
    try:
        import PySide  # pylint: disable=unused-import
        qt_int = 4
        qt_version = 'pyside'
        found_gui = True
    except ImportError:
        pass

    if not found_gui:
        try:
            import PySide2  # pylint: disable=unused-import
            qt_int = 5
            qt_version = 'pyside2'
            found_gui = True
        except ImportError:
            pass

    if not found_gui:
        try:
            import PyQt5  # pylint: disable=unused-import
            qt_int = 5
            qt_version = 'pyqt5'
            found_gui = True
        except ImportError:
            pass

    if not found_gui:
        try:
            import PyQt4  # pylint: disable=unused-import
            qt_int = 4
            qt_version = 'pyqt4'
            found_gui = True
        except ImportError:
            pass
    if not found_gui:
        raise ImportError('PyQt4, PyQt5, PySide, or PySide2 is required')

#if qt_version in ['pyside', 'pyside2']:
    #from qtpy import PYSIDE_VERSION as PYQT_VERSION
#elif qt_version in ['pyqt4', 'pyqt5']:
    #from qtpy import PYQT_VERSION
#else:
    #raise NotImplementedError(qt_version)

from qtpy import API as qt_version

if qt_version in ['pyqt', 'pyqt4']:
    qt_int = 4
    qt_version = 'pyqt4'
    from qtpy import PYQT_VERSION  # pylint: disable=unused-import
elif qt_version == 'pyqt5':
    qt_int = 5
    from qtpy import PYQT_VERSION  # pylint: disable=unused-import
elif qt_version == 'pyside':
    qt_int = 4
    from qtpy import PYSIDE_VERSION as PYQT_VERSION  # pylint: disable=unused-import
elif qt_version == 'pyside2':
    qt_int = 5
    from qtpy import PYSIDE_VERSION as PYQT_VERSION  # pylint: disable=unused-import
else:
    raise ImportError('PyQt4, PyQt5, PySide, or PySide2 is required; API=%r' % qt_version)

if qt_version not in ['pyqt4', 'pyqt5', 'pyside', 'pyside2']:
    raise ImportError('PyQt4, PyQt5, PySide, or PySide2 is required; API=%r' % qt_version)

# required to make a pretty console
try:
    import pygments  # pylint: disable=unused-import
    is_pygments = True
except ImportError:
    is_pygments = False
