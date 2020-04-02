"""
An attempt at a user friendly Cart3d GUI
"""
# -*- coding: utf-8 -*-
import sys
import os.path

# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


from qtpy import QtCore, QtGui
from qtpy.QtWidgets import QApplication

# 3rd party
import vtk

# pyNastran
import pyNastran
from pyNastran.gui.formats import Cart3dIO#, CLASS_MAP
from pyNastran.gui.arg_handling import get_inputs
#from pyNastran.gui.qt_files.gui_qt_common import GuiQtCommon
from pyNastran.gui.gui_common import GuiCommon2


try:
    PKG_PATH = sys._MEIPASS #@UndefinedVariable
    SCRIPT_PATH = os.path.join(PKG_PATH, 'scripts')
    ICON_PATH = os.path.join(PKG_PATH, 'icons')
except:
    PKG_PATH = pyNastran.__path__[0]
    SCRIPT_PATH = os.path.join(PKG_PATH, 'gui', 'scripts')
    ICON_PATH = os.path.join(PKG_PATH, 'gui', 'icons')


class MainWindow(GuiCommon2):
    def __init__(self, inputs, **kwds):
        html_logging = True
        fmt_order = ['cart3d']

        kwds['inputs'] = inputs
        kwds['fmt_order'] = fmt_order
        kwds['html_logging'] = html_logging
        super(MainWindow, self).__init__(**kwds)
        self.base_window_title = "pyCart3d v%s"  % pyNastran.__version__

        self.build_fmts(fmt_order, stop_on_failure=True)
        logo = os.path.join(ICON_PATH, 'logo.png')
        self.logo = logo
        self.set_script_path(SCRIPT_PATH)
        self.set_icon_path(ICON_PATH)

        self.setup_gui()
        self.setup_post(inputs)

    def _cart3d_remap_bcs(self):
        pass

    def _create_cart3d_tools_and_menu_items(self):
        tools = [
            ('about_cart3d', 'About Cart3d GUI', 'tabout.png', 'CTRL+H', 'About Cart3d GUI and help on shortcuts', self.about_dialog),
            #('about', 'About Orig GUI', 'tabout.png', 'CTRL+H', 'About Nastran GUI and help on shortcuts', self.about_dialog),
        ]
        self.menu_edit = self.menubar.addMenu('Edit Cart3d')
        self.menu_help_cart3d = self.menubar.addMenu('&Help')
        self.menu_help.menuAction().setVisible(False)
        #self.file.menuAction().setVisible(False)
        #self.menu_help.

        #self.actions['about'].Disable()

        menu_items = [
            (self.menu_help_cart3d, ('about_cart3d',)),
            #(self.menu_help, ('load_geometry', 'load_results', 'script', '', 'exit')),
            #(self.menu_help2, ('load_geometry', 'load_results', 'script', '', 'exit')),
        ]
        return tools, menu_items

    def _cleanup_cart3d_tools_and_menu_items(self):
        self.menu_help.menuAction().setVisible(True)
        self.menu_help2.menuAction().setVisible(False)
        self.menu_edit.menuAction().setVisible(False)

    def about_dialog(self):
        """ Display about dialog """
        copyright = pyNastran.__pyqt_copyright__

        about = [
            'pyCart3d Qt GUI',
            '',
            'pyCart3d v%s' % pyNastran.__version__,
            copyright,
            pyNastran.__author__,
            '',
            '%s' % pyNastran.__website__,
            '',
            'Mouse',
            'Left Click - Rotate',
            'Middle Click - Pan/Recenter Rotation Point',
            'Shift + Left Click - Pan/Recenter Rotation Point',
            'Right Mouse / Wheel - Zoom',
            '',
            'Keyboard Controls',
            #'r   - reset camera view',
            #'X/x - snap to x axis',
            #'Y/y - snap to y axis',
            #'Z/z - snap to z axis',
            #'',
            #'h   - show/hide legend & info',
            'CTRL+I - take a screenshot (image)',
            'CTRL+L - cycle results',
            #'m/M    - scale up/scale down by 1.1 times',
            #'o/O    - rotate counter-clockwise/clockwise 5 degrees',
            's      - view model as a surface',
            'w      - view model as a wireframe',
            '',
            'Reload Model:  using the same filename reload the model',
        ]
        QtGui.QMessageBox.about(self, "About pyCart3d GUI", "\n".join(about))

    def on_reload(self):
        case = self.icase
        if self.format == 'usm3d':
            self.step_results_usm3d()
        else:
            self.on_load_geometry(self.infile_name, self.format)

        msg = '%s - %s - %s' % (self.format, self.infile_name, self.out_filename)
        self.window_title = msg
        self.log_command('on_reload()')
        self.cycle_results(case)

    def closeEvent(self, event):
        """
        Handling saving state before application when application is
        being closed.
        """
        settings = QtCore.QSettings()
        settings.setValue("main_WindowGeometry", self.saveGeometry())
        settings.setValue("mainWindowState", self.saveState())
        self.settings.save(settings)

        #screen_shape = QtGui.QDesktopWidget().screenGeometry()
        main_window = self.window()
        width = main_window.frameGeometry().width()
        height = main_window.frameGeometry().height()
        settings.setValue('screen_shape', (width, height))

        qpos = self.pos()
        pos = qpos.x(), qpos.y()
        settings.setValue('pos', pos)

        q_app = QApplication.instance()
        if q_app is None:
            sys.exit()
        q_app.quit()


def main():
    app = QApplication(sys.argv)
    QApplication.setOrganizationName("pyCart3d")
    QApplication.setOrganizationDomain(pyNastran.__website__)
    QApplication.setApplicationName("pyCart3d")
    QApplication.setApplicationVersion(pyNastran.__version__)

    inputs = get_inputs()
    window = MainWindow(inputs)
    sys.exit(app.exec_())

if __name__ == '__main__':  # pragma: no cover
    main()
