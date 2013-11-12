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


def main():
    if 1 in fmodes or 2 in fmodes:
        from pyNastran.gui.gui_qt import main as main2
    elif 0 in fmodes:
        from pyNastran.gui.gui_wx import main as main2
    else:
        raise RuntimeError('wxPython, PyQt4, or PySide must be installed')
    main2()


if __name__ == '__main__':
    main()
