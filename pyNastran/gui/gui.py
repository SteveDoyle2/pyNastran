"""
creates the pyNastranGUI
"""
# coding: utf-8

# we're intentionally putting this here to validate the imports
# before doing lots of work
from pyNastran.gui.arg_handling import get_inputs
get_inputs(print_inputs=True)

import sys
import ctypes
# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


import pyNastran
from pyNastran.gui.main_window import MainWindow, get_stylesheet


def cmd_line():
    """the setup.py entry point for ``pyNastranGUI``"""
    # this fixes the icon shown in the windows taskbar to be the custom one (not the python one)
    if sys.platform == 'win32':
        myappid = 'pynastran.pynastrangui.%s' % (pyNastran.__version__) # arbitrary string
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    from qtpy.QtWidgets import QApplication
    app = QApplication(sys.argv)

    if 0:  # pragma: no cover
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

    #if 1:
    stylesheet = get_stylesheet()
    if stylesheet:
        app.setStyleSheet(stylesheet)


    if 0:  # pragma: no cover
        import qtpy.QtGui as QtGui
        #import qtpy.QtCore as QtCore
        from qtpy.QtGui import QPalette, QColor
        dark_palette = QtGui.QPalette()
        dark_palette.setColor(QPalette.WindowText, QColor(180, 180, 180))
        dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.Light, QColor(180, 180, 180))
        dark_palette.setColor(QPalette.Midlight, QColor(90, 90, 90))
        dark_palette.setColor(QPalette.Dark, QColor(35, 35, 35))
        dark_palette.setColor(QPalette.Text, QColor(180, 180, 180))
        dark_palette.setColor(QPalette.BrightText, QColor(180, 180, 180))
        dark_palette.setColor(QPalette.ButtonText, QColor(180, 180, 180))
        dark_palette.setColor(QPalette.Base, QColor(42, 42, 42))
        dark_palette.setColor(QPalette.Window, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.Shadow, QColor(20, 20, 20))
        dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.HighlightedText, QColor(180, 180, 180))
        dark_palette.setColor(QPalette.Link, QColor(56, 252, 196))
        dark_palette.setColor(QPalette.AlternateBase, QColor(66, 66, 66))
        dark_palette.setColor(QPalette.ToolTipBase, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.ToolTipText, QColor(180, 180, 180))

        # disabled
        dark_palette.setColor(QPalette.Disabled, QPalette.WindowText,
                              QColor(127, 127, 127))
        dark_palette.setColor(QPalette.Disabled, QPalette.Text,
                              QColor(127, 127, 127))
        dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText,
                              QColor(127, 127, 127))
        dark_palette.setColor(QPalette.Disabled, QPalette.Highlight,
                              QColor(80, 80, 80))
        dark_palette.setColor(QPalette.Disabled, QPalette.HighlightedText,
                              QColor(127, 127, 127))
        app.setPalette(dark_palette)

    QApplication.setOrganizationName("pyNastran")
    QApplication.setOrganizationDomain(pyNastran.__website__)
    QApplication.setApplicationName("pyNastran")
    QApplication.setApplicationVersion(pyNastran.__version__)
    inputs = get_inputs(print_inputs=False)
    #inputs['app'] = app
    MainWindow(inputs)
    app.exec_()

if __name__ == '__main__':   # pragma: no cover
    cmd_line()
