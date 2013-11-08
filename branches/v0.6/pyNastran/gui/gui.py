fmodes = []
try:
    import wx
    fmodes.append(0)
except ImportError:
    pass

try:
    from PySide import QtCore, QtGui
    fmodes.append(1)
except ImportError:
    pass

try:
    from PyQt4 import QtCore, QtGui
    fmodes.append(2)
except ImportError:
    pass

if __name__ == '__main__':
    if 1 in fmodes or 2 in fmodes:
        from pyNastran.gui.gui_qt import main
    elif 0 in fmodes:
        from pyNastran.gui.gui_wx import main
    else:
        raise RuntimeError('wxPython, PyQt4, or PySide must be installed')
    main()
