"""
creates the pyNastranGUI
"""
# coding: utf-8
from __future__ import division, unicode_literals, print_function

import ctypes
# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)
import sys


# yes we're intentionally putting this here to validate the imports
# before doing lots of work
from pyNastran.gui.arg_handling import get_inputs
get_inputs()


import pyNastran
from pyNastran.gui.main_window import MainWindow


def cmd_line():
    """the setup.py entry point for ``pyNastranGUI``"""
    # this fixes the icon shown in the windows taskbar to be the custom one (not the python one)
    if sys.platform == 'win32':
        myappid = 'pynastran.pynastrangui.%s' % (pyNastran.__version__) # arbitrary string
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    from qtpy.QtWidgets import QApplication
    app = QApplication(sys.argv)

    try:
        import qtmodern.styles
    except ImportError:
        pass
    else:
        qtmodern.styles.dark(app)

    QApplication.setOrganizationName("pyNastran")
    QApplication.setOrganizationDomain(pyNastran.__website__)
    QApplication.setApplicationName("pyNastran")
    QApplication.setApplicationVersion(pyNastran.__version__)
    inputs = get_inputs()
    MainWindow(inputs)
    app.exec_()

if __name__ == '__main__':
    cmd_line()
