# -*- coding: utf-8 -*-
# pylint: disable=C0111
from __future__ import division, unicode_literals, print_function
#from six import string_types, iteritems
from six.moves import range

# standard library
import sys
import os.path
#import traceback
#import webbrowser
#webbrowser.open("http://xkcd.com/353/")

from PyQt4 import QtCore, QtGui

import pyNastran
from pyNastran.gui.utils import check_for_newer_version

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
#from pyNastran.utils import print_bad_path
from pyNastran.gui.formats import (NastranIO, Cart3dIO, PanairIO, LaWGS_IO,
    STL_IO, TecplotIO, ShabpIO, SurfIO, UGRID_IO,
    )
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


class MainWindow(GuiCommon2, NastranIO, Cart3dIO, ShabpIO, PanairIO,
                 LaWGS_IO, STL_IO, TecplotIO, SurfIO, UGRID_IO):
    """
    MainWindow -> GuiCommon2 -> GuiCommon
    gui.py     -> gui_common -> gui_qt_common

    warp vector might finally work
    http://vtk.1045678.n5.nabble.com/How-to-get-grid-result-from-vtkWarpVector-td5727100.html

    glyphs
    http://www.itk.org/Wiki/VTK/Examples/Python/Visualization/ElevationBandsWithGlyphs

    list of VTK6 classes
    http://www.vtk.org/doc/nightly/html/annotated.html

    background grid
    http://www.vtk.org/Wiki/VTK/Examples/Python/Visualization/CubeAxesActor

    pick visible
    http://www.vtk.org/Wiki/VTK/Examples/Cxx/Filtering/ExtractVisibleCells

    plane projection
    http://www.igstk.org/Wiki/VTK/Examples/Cxx/SimpleOperations/ProjectPointPlane

    warping
    http://engronline.ee.memphis.edu/eece4731/djr_lec16.pdf

    banded filter
    http://www.igstk.org/Wiki/VTK/Examples/Cxx/VisualizationAlgorithms/BandedPolyDataContourFilter

    speeding up vtk cell loading in unstructured grids
    http://vtk.1045678.n5.nabble.com/Speed-up-cell-allocation-td5733208.html#a5733214
    """
    def __init__(self, inputs):
        html_logging = True
        fmt_order = [
            # results
            'nastran', 'cart3d', 'panair', 'shabp', 'tecplot', 'surf', 'ugrid',

            # no results
            'lawgs', 'stl',
        ]
        GuiCommon2.__init__(self, fmt_order, html_logging, inputs)

        NastranIO.__init__(self)
        Cart3dIO.__init__(self)
        LaWGS_IO.__init__(self)
        PanairIO.__init__(self)
        STL_IO.__init__(self)
        ShabpIO.__init__(self)
        SurfIO.__init__(self)
        TecplotIO.__init__(self)
        UGRID_IO.__init__(self)

        self.build_fmts(fmt_order, stop_on_failure=False)

        logo = os.path.join(icon_path, 'logo.png')
        self.set_logo(logo)
        self.set_script_path(script_path)
        self.set_icon_path(icon_path)

        self.setup_gui()
        self.setup_post(inputs)
        self._check_for_latest_version()

    def _check_for_latest_version(self):
        version_latest, version_current, is_newer = check_for_newer_version()
        if is_newer:
            url = pyNastran.__website__
            from pyNastran.gui.menus.download import DownloadWindow
            win = DownloadWindow(url, version_latest, win_parent=self)
            win.show()

    def mousePressEvent(self, ev):
        if not self.run_vtk:
            return
        #print('press x,y = (%s, %s)' % (ev.x(), ev.y()))
        if self.is_pick:
            #self.___saveX = ev.x()
            #self.___saveY = ev.y()
            pass
        else:
            self.iren.mousePressEvent(ev)

    #def LeftButtonPressEvent(self, ev):

    def mouseReleaseEvent(self, ev):
        #print('release x,y = (%s, %s)' % (ev.x(), ev.y()))
        if self.is_pick:
            pass
        else:
            self.iren.mousePressEvent(ev)

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
            'p      - pick node/element',
            's      - view model as a surface',
            'w      - view model as a wireframe',
            'f      - set rotation center (zoom out when picking',
            '         to disable clipping)',
            'e      - view model edges',
            'b      - change edge color from scalar/black',
            '',
            'Reload Model:  using the same filename, reload the model',
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
        camera = self.get_camera_data()
        Title = self.title
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
        #self.cycle_results(Title)
        for i in range(10):  #  limit on number of cycles
            if self.title != Title:
                self.cycle_results(Title)
            else:
                break
        self.on_set_camera_data(camera, show_log=False)

    def closeEvent(self, event):
        """
        Handling saving state before application when application is
        being closed.
        """
        settings = QtCore.QSettings()
        settings.setValue("main_WindowGeometry", self.saveGeometry())
        settings.setValue("mainWindowState", self.saveState())
        settings.setValue("backgroundColor", self.background_color)
        settings.setValue("textColor", self.text_color)
        settings.setValue("labelColor", self.label_color)

        #screen_shape = QtGui.QDesktopWidget().screenGeometry()
        main_window = self.window()
        width = main_window.frameGeometry().width()
        height = main_window.frameGeometry().height()
        settings.setValue('screen_shape', (width, height))

        qpos = self.pos()
        pos = qpos.x(), qpos.y()
        settings.setValue('pos', pos)
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
