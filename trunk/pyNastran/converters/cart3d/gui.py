"""
An attempt at a user friendly Cart3d GUI
"""
# -*- coding: utf-8 -*-
from __future__ import division, unicode_literals, print_function
from six import string_types, iteritems
from six.moves import range

# standard library
import sys
import os.path
import traceback

# 3rd party
from numpy import ndarray, eye
import vtk
from PyQt4 import QtCore, QtGui

# pyNastran
import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.gui.formats import Cart3dIO, is_cart3d
from pyNastran.gui.arg_handling import get_inputs
from pyNastran.gui.qt_files.gui_qt_common import GuiCommon
from pyNastran.gui.gui_common import GuiCommon2


# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


if not is_cart3d:
    raise ImportError()

try:
    pkg_path = sys._MEIPASS #@UndefinedVariable
    script_path = os.path.join(pkg_path, 'scripts')
    icon_path = os.path.join(pkg_path, 'icons')
except:
    pkg_path = pyNastran.__path__[0]
    script_path = os.path.join(pkg_path, 'gui', 'scripts')
    icon_path = os.path.join(pkg_path, 'gui', 'icons')


class MainWindow(GuiCommon2, Cart3dIO):
    def __init__(self, inputs):
        html_logging = True
        GuiCommon2.__init__(self, html_logging, inputs)
        Cart3dIO.__init__(self)
        self.base_window_title = "pyCart3d v%s"  % pyNastran.__version__

        fmt_order = ['cart3d']
        self.build_fmts(fmt_order, stop_on_failure=True)
        logo = os.path.join(icon_path, 'logo.png')
        self.set_logo(logo)
        self.set_script_path(script_path)
        self.set_icon_path(icon_path)

        self.setup_gui()
        self.setup_post(inputs)

    def create_cell_picker(self):
        # cell picker
        self.cell_picker = vtk.vtkCellPicker()

    def init_cell_picker(self):
        self.is_pick = False
        self.vtk_interactor.SetPicker(self.cell_picker)
        #self.vtk_interactor.SetPicker(self.point_picker)

        def annotate_cell_picker(object, event):
            self.log_command("annotate_cell_picker()")
            picker = self.cell_picker
            if picker.GetCellId() < 0:
                #self.picker_textActor.VisibilityOff()
                pass
            else:
                worldPosition = picker.GetPickPosition()
                cell_id = picker.GetCellId()
                #ds = picker.GetDataSet()
                selPt = picker.GetSelectionPoint()
                self.log_command("annotate_picker()")
                self.log_info("worldPosition = %s" % str(worldPosition))
                self.log_info("cell_id = %s" % cell_id)
                #self.log_info("data_set = %s" % ds)
                self.log_info("selPt = %s" % str(selPt))

                #self.picker_textMapper.SetInput("(%.6f, %.6f, %.6f)"% pickPos)
                #self.picker_textActor.SetPosition(selPt[:2])
                #self.picker_textActor.VisibilityOn()

        self.cell_picker.AddObserver("EndPickEvent", annotate_cell_picker)

    def on_cell_picker(self):
        self.log_command("on_cell_picker()")
        picker = self.cell_picker
        worldPosition = picker.GetPickPosition()
        cell_id = picker.GetCellId()
        #ds = picker.GetDataSet()
        selPt = picker.GetSelectionPoint()  # get x,y pixel coordinate

        self.log_info("worldPosition = %s" % str(worldPosition))
        self.log_info("cell_id = %s" % cell_id)
        self.log_info("selPt = %s" % str(selPt))
        #self.log_info("data_set = %s" % ds)

    def about_dialog(self):
        """ Display about dialog """
        #if fmode == 1:  # PyQt
        copyright = pyNastran.__pyqt_copyright__
        #else:
            #copyright = pyNastran.__copyright__

        about = [
            'pyCart3d QT GUI',
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
        Title = self.Title
        if self.format == 'usm3d':
            self.step_results_usm3d()
        else:
            self.on_load_geometry(self.infile_name, self.format)

        msg = '%s - %s - %s' % (self.format, self.infile_name, self.out_filename)
        self.set_window_title(msg)
        self.log_command('on_reload()')
        #self.cycleResults(Title)
        for i in range(10):  #  limit on number of cycles
            if self.Title != Title:
                self.cycleResults(Title)
            else:
                break

    def closeEvent(self, event):
        """
        Handling saving state before application when application is
        being closed.
        """
        settings = QtCore.QSettings()
        settings.setValue("main_WindowGeometry", self.saveGeometry())
        settings.setValue("mainWindowState", self.saveState())
        settings.setValue("backgroundColor", self.background_col)
        QtGui.qApp.quit()

def main():
    app = QtGui.QApplication(sys.argv)
    QtGui.QApplication.setOrganizationName("pyCart3d")
    QtGui.QApplication.setOrganizationDomain(pyNastran.__website__)
    QtGui.QApplication.setApplicationName("pyCart3d")
    QtGui.QApplication.setApplicationVersion(pyNastran.__version__)

    inputs = get_inputs()
    window = MainWindow(inputs)
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
