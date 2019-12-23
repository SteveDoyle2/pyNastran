"""
Selects a matplotlib backend so you can run without a GUI/tkinter.  Supports:
 - PyQt5
 - PySide2
 - WX
 - Tkinter

"""
from pyNastran.gui import IS_DEV
if IS_DEV:
    # there is no interactive backend when testing on TravisCI
    matplotlib_backend = 'Agg'
else:
    # fails if using the terminal and PyQt/PySide & qtpy are installed
    # how do I check if there is a terminal vs just running in command line?
    #
    try:
        from pyNastran.gui.qt_version import qt_int
        matplotlib_backend = 'Qt%iAgg' % qt_int
    except ImportError:
        try:
            # hasn't been tested on a machine without a backend...
            # default matplotlib backend
            import tkinter
            matplotlib_backend = 'tkAgg'
        except ImportError:
            # no-gui backend
            matplotlib_backend = 'Agg'
