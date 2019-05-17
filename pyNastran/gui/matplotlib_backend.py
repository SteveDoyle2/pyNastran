"""
Selects a matplotlib backend so you can run without a GUI/tkinter.

Unfortunately, it currently requires qtpy to use Qt4Agg or Qt5Agg.
It's only used for interactive testing at the moment, so that can be
fixed later.

"""
from pyNastran.gui import IS_DEV
if IS_DEV:
    # there is no interactive backend when testing on TravisCI
    matplotlib_backend = 'Agg'
else:
    # fails if using the terminal and PyQt/PySide & qtpy are installed
    # how do I check if there is a terminal vs just running in command line?
    #
    # probably should support tkinter, which is the default matplotlib backend
    try:
        from pyNastran.gui.qt_version import qt_int
        matplotlib_backend = 'Qt%iAgg' % qt_int
    except ImportError:
        matplotlib_backend = 'Agg'
