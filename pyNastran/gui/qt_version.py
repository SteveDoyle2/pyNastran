"""
Figures out the "optimal" Qt version to use in a way that:

 - uses the qt version specified by the QT_API environment variable
 - picks up the already imported version
 - selects PySide2 or PyQt5 (in that order)

"""
import os
import sys

# there's a bug with the keyboard shortcuts K and L
ALLOW_PYSIDE6 = False
ALLOW_PYQT6 = False

#if ALLOW_PYQT6 and ALLOW_PYSIDE6:  # pragma: no cover
    #_msg = f'PyQt5/6 or PySide2/6 is required'
#elif ALLOW_PYSIDE6:
    #_msg = f'PyQt5 or PySide2/6 is required'
#if ALLOW_PYQT6:
    #_msg = f'PyQt5/6 or PySide2 is required'
#else:
_msg = f'PyQt5 or PySide2 is required'


API = os.environ.get('QT_API', '').lower()
if API:
    from qtpy import API as qt_version
elif 'PySide2' in sys.modules:
    qt_version = 'pyside2'
elif 'PyQt5' in sys.modules:
    qt_version = 'pyqt5'
elif 'PySide6' in sys.modules and ALLOW_PYSIDE6:  # pragma: no cover
    qt_version = 'pyside6'
elif 'PyQt6' in sys.modules and ALLOW_PYQT6:  # pragma: no cover
    qt_version = 'pyqt6'
else:
    found_gui = False
    try:
        import PySide2  # pylint: disable=unused-import
        qt_int = 5
        qt_version = 'pyside2'
        found_gui = True
    except ImportError:  # pragma: no cover
        pass

    if not found_gui:
        try:
            import PyQt5  # pylint: disable=unused-import
            qt_int = 5
            qt_version = 'pyqt5'
            found_gui = True
        except ImportError:  # pragma: no cover
            pass

    if not found_gui and ALLOW_PYSIDE6:  # pragma: no cover
        try:
            import PySide6  # pylint: disable=unused-import
            qt_int = 6
            qt_version = 'pyside6'
            found_gui = True
        except ImportError:  # pragma: no cover
            pass

    if not found_gui and ALLOW_PYQT6:  # pragma: no cover
        try:
            import PyQt6  # pylint: disable=unused-import
            qt_int = 6
            qt_version = 'pyqt6'
            found_gui = True
        except ImportError:
            pass

    if not found_gui:  # pragma: no cover
        raise ImportError(_msg)


from qtpy import API as qt_version

if qt_version == 'pyqt5':
    qt_int = 5
    qt_name = 'PyQt5'
    from qtpy import PYQT_VERSION  # pylint: disable=unused-import
elif qt_version == 'pyside2':
    qt_int = 5
    qt_name = 'PySide2'
    from qtpy import PYSIDE_VERSION as PYQT_VERSION  # pylint: disable=unused-import
elif ALLOW_PYSIDE6 and qt_version == 'pyside6':  # pragma: no cover
    qt_int = 6
    qt_name = 'PySide6'
    from qtpy import PYSIDE_VERSION as PYQT_VERSION  # pylint: disable=unused-import
elif ALLOW_PYQT6 and qt_version == 'pyqt6':  # pragma: no cover
    qt_int = 6
    qt_name = 'PyQt6'
    from qtpy import PYQT_VERSION  # pylint: disable=unused-import
else:  # pragma: no cover
    raise ImportError(f'{_msg}; API={qt_version!r}')

# required to make a pretty console
try:
    import pygments  # pylint: disable=unused-import
    is_pygments = True
except ModuleNotFoundError:
    is_pygments = False

QT_AGG_BACKENDS = {
    'PyQt5': 'Qt5Agg',
    'PySide2': 'Qt5Agg',
    'PyQt6': 'QtAgg',
    'PySide6': 'QtAgg',
}
QT_AGG_BACKEND = QT_AGG_BACKENDS[qt_name]
del QT_AGG_BACKENDS
