"""
creates the pyNastranGUI
"""
# coding: utf-8
from __future__ import division, unicode_literals, print_function

# we're intentionally putting this here to validate the imports
# before doing lots of work
from pyNastran.gui.arg_handling import get_inputs
get_inputs()

import sys
import ctypes
# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


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

    if 0:  # pragma: no cover
        import qtpy.QtGui as QtGui
        import qtpy.QtCore as QtCore
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
    inputs = get_inputs()
    #inputs['app'] = app
    MainWindow(inputs)
    app.exec_()

def cmd_line2():  # pragma: no cover
    """Simple program that greets NAME for a total of COUNT times."""
    import argparse

    #msg = "Usage:\n"

    # INPUT format may be explicitly or implicitly defined with or
    # without an output file
    #test = ' [--test]'
    #qt = ' [--qt QT] [--plugin]'

    #msg += "  pyNastranGUI INPUT [-f FORMAT] [-o OUTPUT]\n"
    #msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    #msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    #msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)

    # You don't need to throw a -o flag
    #msg += "  pyNastranGUI INPUT OUTPUT [-f FORMAT] [-o OUTPUT]\n"
    #msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    #msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    #msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)

    dev = ''
    dev_list = []
    if not pyNastran.is_pynastrangui_exe:
        #dev = ' [--noupdate] [--test] [--qt Qt] [--plugin]'
        dev_list = ['--noupdate', '--test', '--qt', '--plugin']
        dev = ''.join([' [%s]' % devi for devi in dev_list])

    # no input/output files
    # can you ever have an OUTPUT, but no INPUT?
    usage = "  pyNastranGUI INPUT [-f FORMAT] [-o OUTPUT]\n"
    usage += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    usage += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    usage += '               [-q] [--groups] [--noupdate] [--log LOG]%s\n' % (dev)

    #parent_parser.add_argument('-g', '--geomscript', type=str, help='path to geometry script file (runs before load geometry)', action='append')
    #parent_parser.add_argument('-p', '--postscript', type=str, help='path to post script file (runs after load geometry)', action='append')
    #parent_parser.add_argument('-u', '--points_fname', type=str, help='an (nrows, 3) comma/tab/space separated list of points')
    #parent_parser.add_argument('--user_geom', type=str, help='add user specified geometry (repeatable)')

    # You don't need to throw a -o flag
    usage += "  pyNastranGUI INPUT OUTPUT [-f FORMAT] [-o OUTPUT]\n"
    usage += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    usage += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    usage += '               [-q] [--groups] [--log LOG]%s\n' % (dev)

    # no input/output files
    # can you ever have an OUTPUT, but no INPUT?
    usage += "  pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT...]\n"
    usage += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    usage += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    usage += '               [-q] [--groups] [--log LOG]%s\n' % (dev)
    #usage += '  pyNastranGUI -h | --help\n'
    usage += '  pyNastranGUI -v | --version\n'

    #msg += "\n"
    #parser = argparse.ArgumentParser(
        #prog=None, usage=None, description=None, epilog=None,
        #version=None, parents=[], formatter_class=HelpFormatter,
        #prefix_chars='-', fromfile_prefix_chars=None, argument_default=None,
        #conflict_handler='error', add_help=True)

    #usage = '[options]'
    text = (
        'Examples\n'
        '--------\n'
        '  pyNastranGUI\n'
        '  pyNastranGUI fem.bdf\n'
        '  pyNastranGUI fem.bdf fem.op2\n'
        '  pyNastranGUI --format nastran fem.dat fem.op2 -o fem2.op2\n'
    )
    import textwrap
    parent_parser = argparse.ArgumentParser(
        #prog = 'pyNastranGUI',
        #usage = usage,
        #description='A foo that bars',
        epilog="And that's how you'd foo a bar",
        #formatter_class=argparse.RawDescriptionHelpFormatter,
        #description=textwrap.dedent(text),
        #version=pyNastran.__version__,
        #add_help=False,
    )
    # positional arguments
    parent_parser.add_argument('INPUT', nargs='?', help='path to input file', type=str)
    parent_parser.add_argument('OUTPUT', nargs='?', help='path to output file', type=str)

    parent_parser.add_argument('-i', '--input', help='path to input file')
    parent_parser.add_argument('-o', '--output', help='path to output file')
    #parent_parser.add_argument('--user_geom', type=str, help='log msg')

    # double args
    parent_parser.add_argument('-f', '--format', type=str,
                               help='format type (avus, bedge, cart3d, lawgs, nastran, '
                               'openfoam_hex, openfoam_shell, openfoam_faces, panair, '
                               'stl, surf, tetgen, usm3d, ugrid, ugrid3d, #plot3d)', action='append')
    parent_parser.add_argument('-g', '--geomscript', type=str, help='path to geometry script file (runs before load geometry)', action='append')
    parent_parser.add_argument('-p', '--postscript', type=str, help='path to post script file (runs after load geometry)', action='append')
    parent_parser.add_argument('-u', '--points_fname', type=str, help='an (nrows, 3) comma/tab/space separated list of points')
    parent_parser.add_argument('--user_geom', type=str, help='add user specified geometry (repeatable)')
    parent_parser.add_argument('--log', type=str, help='{debug, info, warning, error} msg')

    # no arguments
    if dev:
        parent_parser.add_argument('--qt', type=str, help='{pyqt4, pyqt5, pyside, pyside2} msg')
        parent_parser.add_argument('--test', help='test msg', action='store_true')
        parent_parser.add_argument('--noupdate', help='noupdate msg', action='store_true')
    parent_parser.add_argument('--groups', help='enables groups', action='store_true')
    parent_parser.add_argument('--plugin', help='plugin msg', action='store_true')

    parent_parser.add_argument('-q', '--quiet', help='prints debug messages (default=True)', action='store_true')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parent_parser.add_argument('-v', '--version', action='version',
                               version=pyNastran.__version__)

    #foo_parser = argparse.ArgumentParser(parents=[parent_parser])
    #foo_parser.parse_args(['INPUT', '--format', '--output',
                           #'--geomscript', '--postscript', '--points_fname', '--user_geom',
                           #'--quiet', '--groups', '--no_update', '--log' '--help'] + dev_list)

    #msg += "  pyNastranGUI INPUT [-f FORMAT] [-o OUTPUT]\n"
    #msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    #msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    #msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)
    #parser_no_output = p

    #parser = argparse.ArgumentParser(
        #description='A foo that bars',
        #epilog="And that's how you'd foo a bar",
        #version=pyNastran.__version__,
    #)
    #parser.add_argument("square", help="display a square of a given number",
                        #type=int)
    #parser.add_argument('-v', '--version', action='version',
                        #version='%%(prog)s %s' % pyNastran.__version__)
    #parser.add_argument("-w", "--verbosity", type=int, choices=[0, 1, 2],
                        #help="increase output verbosity")
    print('a')

    mymsg = 'replacing argparse message'
    #def print_help(self, file=None):
        #if file is None:
            #file = _sys.stdout
        #self._print_message(self.format_help(), file)

    def _print_message(message, file=None):
        if message:
            if file is None:
                file = _sys.stderr
            file.write(mymsg)
    parent_parser._print_message = _print_message
    args = parent_parser.parse_args()

    print('b')
    print(args)
    import sys
    sys.exit()


if __name__ == '__main__':
    cmd_line()
