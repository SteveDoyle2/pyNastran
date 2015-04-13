# -*- coding: utf-8 -*-
from __future__ import division, unicode_literals, print_function
from six import string_types, iteritems
from six.moves import range

# standard library
import sys
import os.path
#import traceback
#import webbrowser
#webbrowser.open("http://xkcd.com/353/")

from PyQt4 import QtCore, QtGui
print("Using PyQt4")
fmode = 1
#except ImportError:
    #try:
        #from PySide import QtCore, QtGui
        #print("Using PySide")
        #fmode = 2
    #except ImportError:
        #msg = 'Failed to import PySide or PyQt4'
        #raise ImportError(msg)
assert fmode in [1, 2]

# 3rd party
import vtk

# pyNastran
import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.gui.formats import (NastranIO, Cart3dIO, PanairIO, LaWGS_IO,
    STL_IO, TecplotIO, TetgenIO, Usm3dIO, Plot3d_io, ShabpIO,
    is_nastran, is_cart3d, is_panair, is_lawgs,
    is_shabp, is_stl, is_tecplot, is_tetgen, is_usm3d, is_plot3d)
from pyNastran.gui.arg_handling import get_inputs
from pyNastran.gui.gui_common import GuiCommon2

try:
    pkg_path = sys._MEIPASS #@UndefinedVariable
    script_path = os.path.join(pkg_path, 'scripts')
    icon_path = os.path.join(pkg_path, 'icons')
except:
    pkg_path = pyNastran.__path__[0]
    script_path = os.path.join(pkg_path, 'gui', 'scripts')
    icon_path = os.path.join(pkg_path, 'gui', 'icons')

#print('script_path = %s' % script_path)
#print('icon_path = %s' % icon_path)

# tcolorpick.png and tabout.png trefresh.png icons on LGPL license, see
# http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/color-picker-grey.png
# http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/help-hint.png
# http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/view-refresh-8.png


# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


class MainWindow(GuiCommon2, NastranIO, Cart3dIO, ShabpIO, PanairIO, LaWGS_IO, STL_IO, TetgenIO, Usm3dIO, Plot3d_io):

    def __init__(self, inputs):
        html_logging = True
        GuiCommon2.__init__(self, html_logging, inputs)

        NastranIO.__init__(self)
        Cart3dIO.__init__(self)
        PanairIO.__init__(self)
        ShabpIO.__init__(self)
        LaWGS_IO.__init__(self)
        STL_IO.__init__(self)
        TetgenIO.__init__(self)
        Usm3dIO.__init__(self)
        Plot3d_io.__init__(self)

        fmt_order = [
            'nastran', 'cart3d', 'panair', 'shabp', 'usm3d',  # results
            'lawgs', 'tetgen', 'stl',  #'plot3d',  # no results
        ]
        self.build_fmts(fmt_order, stop_on_failure=False)

        logo = os.path.join(icon_path, 'logo.png')
        self.set_logo(logo)
        self.set_script_path(script_path)
        self.set_icon_path(icon_path)

        self.setup_gui()
        self.setup_post(inputs)

    def create_cell_picker(self):
        # cell picker
        self.cell_picker = vtk.vtkCellPicker()
        #self.point_picker = vtk.vtkPointPicker()
        #self.cell_picker.SetTolerance(0.0005)

    def mousePressEvent(self, ev):
        print('press x,y = (%s, %s)' % (ev.x(), ev.y()))
        if self.is_pick:
            #self.___saveX = ev.x()
            #self.___saveY = ev.y()
            pass
        else:
            self.iren.mousePressEvent(ev)

    #def LeftButtonPressEvent(self, ev):
        #asfd

    def mouseReleaseEvent(self, ev):
        print('release x,y = (%s, %s)' % (ev.x(), ev.y()))
        if self.is_pick:
            pass
        else:
            self.iren.mousePressEvent(ev)

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
                world_position = picker.GetPickPosition()
                cell_id = picker.GetCellId()
                #ds = picker.GetDataSet()
                select_point = picker.GetSelectionPoint()
                self.log_command("annotate_picker()")
                self.log_info("world_position = %s" % str(world_position))
                self.log_info("cell_id = %s" % cell_id)
                #self.log_info("data_set = %s" % ds)
                self.log_info("selPt = %s" % str(select_point))

                #self.picker_textMapper.SetInput("(%.6f, %.6f, %.6f)"% pickPos)
                #self.picker_textActor.SetPosition(select_point[:2])
                #self.picker_textActor.VisibilityOn()

        def annotate_point_picker(object, event):
            self.log_command("annotate_point_picker()")
            picker = self.cell_picker
            if picker.GetPointId() < 0:
                #self.picker_textActor.VisibilityOff()
                pass
            else:
                world_position = picker.GetPickPosition()
                point_id = picker.GetPointId()
                #ds = picker.GetDataSet()
                select_point = picker.GetSelectionPoint()
                self.log_command("annotate_picker()")
                self.log_info("world_position = %s" % str(world_position))
                self.log_info("point_id = %s" % point_id)
                #self.log_info("data_set = %s" % ds)
                self.log_info("select_point = %s" % str(select_point))

                #self.picker_textMapper.SetInput("(%.6f, %.6f, %.6f)"% pickPos)
                #self.picker_textActor.SetPosition(select_point[:2])
                #self.picker_textActor.VisibilityOn()

        self.cell_picker.AddObserver("EndPickEvent", annotate_cell_picker)
        #self.point_picker.AddObserver("EndPickEvent", annotate_point_picker)

    def on_cell_picker(self):
        self.log_command("on_cell_picker()")
        picker = self.cell_picker
        world_position = picker.GetPickPosition()
        cell_id = picker.GetCellId()
        #ds = picker.GetDataSet()
        select_point = picker.GetSelectionPoint()  # get x,y pixel coordinate

        self.log_info("world_position = %s" % str(world_position))
        self.log_info("cell_id = %s" % cell_id)
        self.log_info("select_point = %s" % str(select_point))
        #self.log_info("data_set = %s" % ds)

    def about_dialog(self):
        """ Display about dialog """
        if fmode == 1:  # PyQt
            copyright = pyNastran.__pyqt_copyright__
        else:
            copyright = pyNastran.__copyright__

        about = [
            'pyNastran QT GUI',
            '',
            'pyNastran v%s' % pyNastran.__version__,
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
        QtGui.QMessageBox.about(self, "About pyNastran GUI", "\n".join(about))

    def on_reload(self):
        """
        Runs the reload button.

        Reload allows you to edit the input model and "reload" the data without
        having to go to the pulldown menu.  For USM3D, we dynamically load the
        latest CFD results time step, which is really handy when you're running
        a job.
        """
        Title = self.Title
        if self.format == 'usm3d':
            self.step_results_usm3d()
        else:
            self.on_load_geometry(self.infile_name, self.format)

        if self.out_filename is None:
            msg = '%s - %s' % (self.format, self.infile_name)
        else:
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
    QtGui.QApplication.setOrganizationName("pyNastran")
    QtGui.QApplication.setOrganizationDomain(pyNastran.__website__)
    QtGui.QApplication.setApplicationName("pyNastran")
    QtGui.QApplication.setApplicationVersion(pyNastran.__version__)

    inputs = get_inputs()
    window = MainWindow(inputs)
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
