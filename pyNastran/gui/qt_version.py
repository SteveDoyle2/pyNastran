"""
Figures out the "optimal" Qt version to use in a way that:

 - uses the qt version specified by the QT_API environment variable
 - picks up the already imported version
 - selects PySide2, PyQt5 (in that order)

"""
import os
import sys

API = os.environ.get('QT_API', '').lower()
if API:
    from qtpy import API as qt_version
elif 'PySide2' in sys.modules:
    qt_version = 'pyside2'
elif 'PyQt5' in sys.modules:
    qt_version = 'pyqt5'
elif 'PySide6' in sys.modules:
    qt_version = 'pyside6'
elif 'PyQt6' in sys.modules:
    qt_version = 'pyqt6'
else:
    found_gui = False
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
            import PySide6  # pylint: disable=unused-import
            qt_int = 6
            qt_version = 'pyside6'
            found_gui = True
        except ImportError:
            pass

    if not found_gui:
        try:
            import PyQt6  # pylint: disable=unused-import
            qt_int = 6
            qt_version = 'pyqt6'
            found_gui = True
        except ImportError:
            pass

    if not found_gui:
        raise ImportError('PyQt5/6 or PySide2/6 is required')

#if qt_version == 'pyside2':
    #from qtpy import PYSIDE_VERSION as PYQT_VERSION
#elif qt_version == 'pyqt5':
    #from qtpy import PYQT_VERSION
#else:
    #raise NotImplementedError(qt_version)

from qtpy import API as qt_version

if qt_version == 'pyqt5':
    qt_int = 5
    qt_name = 'PyQt5'
    from qtpy import PYQT_VERSION  # pylint: disable=unused-import
elif qt_version == 'pyside2':
    qt_int = 5
    qt_name = 'PySide2'
    from qtpy import PYSIDE_VERSION as PYQT_VERSION  # pylint: disable=unused-import
elif qt_version == 'pyside6':
    qt_int = 6
    qt_name = 'PySide6'
    from qtpy import PYSIDE_VERSION as PYQT_VERSION  # pylint: disable=unused-import
elif qt_version == 'pyqt6':
    qt_int = 6
    qt_name = 'PyQt6'
    from qtpy import PYQT_VERSION  # pylint: disable=unused-import
else:
    raise ImportError('PyQt5/6 or PySide2/6 is required; API=%r' % qt_version)

if qt_version not in ['pyqt5', 'pyside2', 'pyqt6', 'pyside6', ]:
    raise ImportError('PyQt5/6 or PySide2/6 is required; API=%r' % qt_version)

# required to make a pretty console
try:
    import pygments  # pylint: disable=unused-import
    is_pygments = True
except ImportError:
    is_pygments = False

QT_AGG_BACKENDS = {
    'PyQt5': 'Qt5Agg',
    'PySide2': 'Qt5Agg',
    #'PyQt6': 'QtAgg',
    'PySide6': 'QtAgg',
}
QT_AGG_BACKEND = QT_AGG_BACKENDS[qt_name]
del QT_AGG_BACKENDS
