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

    if 0:
        try:
            import qtmodern.styles
        except ImportError:
            pass
        else:
            qtmodern.styles.dark(app)

    #app.setStyle('Fusion')
    #app.setStyle('WindowsXP')

    #if 0:
        #import qtpy.QtGui as QtGui
        #import qtpy.QtCore as QtCore
        #palette = QtGui.QPalette()
        #palette.setColor(QtGui.QPalette.Window, QtGui.QColor(53,53,53))
        #palette.setColor(QtGui.QPalette.WindowText, QtCore.Qt.white)
        #palette.setColor(QtGui.QPalette.Base, QtGui.QColor(15,15,15))
        #palette.setColor(QtGui.QPalette.AlternateBase, QtGui.QColor(53,53,53))
        #palette.setColor(QtGui.QPalette.ToolTipBase, QtCore.Qt.white)
        #palette.setColor(QtGui.QPalette.ToolTipText, QtCore.Qt.white)
        #palette.setColor(QtGui.QPalette.Text, QtCore.Qt.white)
        #palette.setColor(QtGui.QPalette.Button, QtGui.QColor(53,53,53))
        #palette.setColor(QtGui.QPalette.ButtonText, QtCore.Qt.white)
        #palette.setColor(QtGui.QPalette.BrightText, QtCore.Qt.red)

        #palette.setColor(QtGui.QPalette.Highlight, QtGui.QColor(142,45,197).lighter())
        #palette.setColor(QtGui.QPalette.HighlightedText, QtCore.Qt.black)
        #app.setPalette(palette)

    if 0:
        import qtpy.QtGui as QtGui
        import qtpy.QtCore as QtCore
        from qtpy.QtGui import QPalette, QColor
        darkPalette = QtGui.QPalette()
        darkPalette.setColor(QPalette.WindowText, QColor(180, 180, 180))
        darkPalette.setColor(QPalette.Button, QColor(53, 53, 53))
        darkPalette.setColor(QPalette.Light, QColor(180, 180, 180))
        darkPalette.setColor(QPalette.Midlight, QColor(90, 90, 90))
        darkPalette.setColor(QPalette.Dark, QColor(35, 35, 35))
        darkPalette.setColor(QPalette.Text, QColor(180, 180, 180))
        darkPalette.setColor(QPalette.BrightText, QColor(180, 180, 180))
        darkPalette.setColor(QPalette.ButtonText, QColor(180, 180, 180))
        darkPalette.setColor(QPalette.Base, QColor(42, 42, 42))
        darkPalette.setColor(QPalette.Window, QColor(53, 53, 53))
        darkPalette.setColor(QPalette.Shadow, QColor(20, 20, 20))
        darkPalette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        darkPalette.setColor(QPalette.HighlightedText, QColor(180, 180, 180))
        darkPalette.setColor(QPalette.Link, QColor(56, 252, 196))
        darkPalette.setColor(QPalette.AlternateBase, QColor(66, 66, 66))
        darkPalette.setColor(QPalette.ToolTipBase, QColor(53, 53, 53))
        darkPalette.setColor(QPalette.ToolTipText, QColor(180, 180, 180))

        # disabled
        darkPalette.setColor(QPalette.Disabled, QPalette.WindowText,
                             QColor(127, 127, 127))
        darkPalette.setColor(QPalette.Disabled, QPalette.Text,
                             QColor(127, 127, 127))
        darkPalette.setColor(QPalette.Disabled, QPalette.ButtonText,
                             QColor(127, 127, 127))
        darkPalette.setColor(QPalette.Disabled, QPalette.Highlight,
                             QColor(80, 80, 80))
        darkPalette.setColor(QPalette.Disabled, QPalette.HighlightedText,
                             QColor(127, 127, 127))
        app.setPalette(darkPalette)

    QApplication.setOrganizationName("pyNastran")
    QApplication.setOrganizationDomain(pyNastran.__website__)
    QApplication.setApplicationName("pyNastran")
    QApplication.setApplicationVersion(pyNastran.__version__)
    inputs = get_inputs()
    #inputs['app'] = app
    MainWindow(inputs)
    app.exec_()

if __name__ == '__main__':
    cmd_line()
