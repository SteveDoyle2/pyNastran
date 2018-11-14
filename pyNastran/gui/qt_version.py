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
    print('using default API')
    found_gui = False
    try:
        import PySide
        qt_int = 4
        qt_version = 'pyside'
        from qtpy import PYSIDE_VERSION as PYQT_VERSION
        found_gui = True
    except ImportError:
        pass

    if not found_gui:
        try:
            import PySide2
            qt_int = 5
            qt_version = 'pyside2'
            from qtpy import PYSIDE_VERSION as PYQT_VERSION
            found_gui = True
        except ImportError:
            pass

    if not found_gui:
        try:
            import PyQt5
            qt_int = 5
            qt_version = 'pyqt5'
            from qtpy import PYQT_VERSION
            found_gui = True
        except ImportError:
            pass

    if not found_gui:
        try:
            import PyQt4
            qt_int = 4
            qt_version = 'pyqt4'
            from qtpy import PYQT_VERSION
            found_gui = True
        except ImportError:
            pass
    if not found_gui:
        raise ImportError('PyQt4, PyQt5, PySide, or PySide2 is required')


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
    raise ImportError('PyQt4, PyQt5, PySide, or PySide2 is required; API=%r' % qt_version)

if qt_version not in ['pyqt4', 'pyqt5', 'pyside', 'pyside2']:
    raise ImportError('PyQt4, PyQt5, PySide, or PySide2 is required; API=%r' % qt_version)

# required to make a pretty console
try:
    import pygments
    is_pygments = True
except ImportError:
    is_pygments = False
