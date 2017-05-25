"""
creates the pyNastranGUI
"""
# coding: utf-8
from __future__ import division, unicode_literals, print_function

# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)
import sys


from pyNastran.gui.arg_handling import get_inputs
# yes we're intentionally putting this here to validate the imports
# before doing lots of work
get_inputs()

import pyNastran
from pyNastran.gui.main_window import QApplication, MainWindow


def cmd_line():
    """the setup.py entry point for ``pyNastranGUI``"""
    app = QApplication(sys.argv)
    QApplication.setOrganizationName("pyNastran")
    QApplication.setOrganizationDomain(pyNastran.__website__)
    QApplication.setApplicationName("pyNastran")
    QApplication.setApplicationVersion(pyNastran.__version__)
    inputs = get_inputs()
    MainWindow(inputs)
    sys.exit(app.exec_())

if __name__ == '__main__':
    cmd_line()
