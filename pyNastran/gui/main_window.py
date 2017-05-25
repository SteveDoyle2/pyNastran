"""
defines the MainWindow class
"""
# coding: utf-8
# pylint: disable=C0111
from __future__ import division, unicode_literals, print_function

# standard library
import sys
import os.path
#import traceback
#import webbrowser
#webbrowser.open("http://xkcd.com/353/")

#from six import string_types, iteritems
from six.moves import range


from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4 import QtCore, QtGui
    from PyQt4.QtGui import QApplication, QMessageBox, qApp
    print("Using PyQt4")
elif qt_version == 5:
    from PyQt5 import QtCore, QtGui
    from PyQt5.QtWidgets import QApplication, QMessageBox, qApp
    print("Using PyQt5")
elif qt_version == 'pyside':
    from PySide import QtCore, QtGui
    from PySide.QtGui import QApplication, QMessageBox, qApp
else:
    raise NotImplementedError(qt_version)

# 3rd party
import vtk

import pyNastran
from pyNastran.gui.gui_utils.utils import check_for_newer_version


# pyNastran
#from pyNastran.utils import print_bad_path
from pyNastran.gui.formats import (
    NastranIO, Cart3dIO, DegenGeomIO, PanairIO, LaWGS_IO,
    STL_IO, TecplotIO, TetgenIO, Usm3dIO, Plot3d_io, ShabpIO, ADB_IO, FastIO,
    AvusIO, SurfIO, UGRID_IO, AbaqusIO, BEdge_IO, SU2_IO,
)
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


class MainWindow(GuiCommon2, NastranIO, Cart3dIO, DegenGeomIO, ShabpIO, PanairIO,
                 LaWGS_IO, STL_IO, TetgenIO, Usm3dIO, TecplotIO, Plot3d_io, ADB_IO,
                 FastIO, AvusIO, SurfIO, UGRID_IO, AbaqusIO, BEdge_IO, SU2_IO):
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
    def __init__(self, inputs, **kwds):
        """
        inputs=None
        """
        html_logging = True
        fmt_order = [
            # results
            'nastran', 'cart3d', 'panair', 'shabp', 'usm3d',
            'tecplot', 'surf', 'ugrid',

            # no results
            'lawgs', 'stl', 'fast',
            'bedge', 'su2', 'tetgen',
        ]
        #GuiCommon2.__init__(self, fmt_order, html_logging, inputs, parent)
        kwds['inputs'] = inputs
        kwds['fmt_order'] = fmt_order
        kwds['html_logging'] = html_logging
        super(MainWindow, self).__init__(**kwds)
        #fmt_order=fmt_order, inputs=inputs,
        #html_logging=html_logging,

        if qt_version in [4, 5]:
            ADB_IO.__init__(self)
            AvusIO.__init__(self)
            BEdge_IO.__init__(self)
            NastranIO.__init__(self)
            Cart3dIO.__init__(self)
            DegenGeomIO.__init__(self)
            FastIO.__init__(self)
            LaWGS_IO.__init__(self)
            PanairIO.__init__(self)
            Plot3d_io.__init__(self)
            STL_IO.__init__(self)
            ShabpIO.__init__(self)
            SurfIO.__init__(self)
            TetgenIO.__init__(self)
            TecplotIO.__init__(self)
            Usm3dIO.__init__(self)
            UGRID_IO.__init__(self)
            AbaqusIO.__init__(self)
            SU2_IO.__init__(self)

        self.build_fmts(fmt_order, stop_on_failure=False)

        logo = os.path.join(icon_path, 'logo.png')
        self.logo = logo
        self.set_script_path(script_path)
        self.set_icon_path(icon_path)

        self.setup_gui()
        self.setup_post(inputs)
        self._check_for_latest_version(inputs['no_update'])

    def _check_for_latest_version(self, check=True):
        """
        checks the website for information regarding the latest gui version

        Looks for:
            ## pyNastran v0.7.2 has been Released (4/25/2015)
        """
        import time
        t0 = time.time()
        version_latest, version_current, is_newer = check_for_newer_version()
        if is_newer:
            url = pyNastran.__website__
            from pyNastran.gui.menus.download import DownloadWindow
            win = DownloadWindow(url, version_latest, win_parent=self)
            win.show()
        dt = time.time() - t0
        #print('dt_version_check = %.2f' % dt)

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
        copyright = pyNastran.__copyright__
        if qt_version == 'pyside':
            word = 'PySide'
            copyright_qt = pyNastran.__pyside_copyright__
        else:
            word = 'PyQt'
            copyright_qt = pyNastran.__pyqt_copyright__

        about = [
            'pyNastran %s GUI' % word,
            '',
            'pyNastran v%s' % pyNastran.__version__,
            copyright,
            copyright_qt,
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
            'CTRL+L - cycle the results forwards',
            'CTRL+K - cycle the results backwards',
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

        #message_box = QMessageBox()
        #message_box.setStyleSheet(
            #'QMessageBox {background-color: #2b5b84; color: white;}\n'
            #'QPushButton{color: white; font-size: 16px; background-color: #1d1d1d; '
            #'border-radius: 10px; padding: 10px; text-align: center;}\n'
            #' QPushButton:hover{color: #2b5b84;}')
        #message_box.setFont(self.font())
        QMessageBox.about(self, "About pyNastran GUI", "\n".join(about))

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
        case = self.icase
        if self.format == 'usm3d':
            self.step_results_usm3d()
        else:
            self.on_load_geometry(self.infile_name, self.format)

        if self.out_filename is None:
            msg = '%s - %s' % (self.format, self.infile_name)
        else:
            msg = '%s - %s - %s' % (self.format, self.infile_name, self.out_filename)
        self.window_title = msg
        self.log_command('on_reload()')
        self.cycle_results(case)
        self.on_set_camera_data(camera, show_log=False)

    def closeEvent(self, event):
        """
        Handling saving state before application when application is
        being closed.
        """
        settings = QtCore.QSettings()
        settings.setValue("mainWindowGeometry", self.saveGeometry())
        settings.setValue("mainWindowState", self.saveState())
        settings.setValue("backgroundColor", self.background_color)
        settings.setValue("textColor", self.text_color)
        settings.setValue("labelColor", self.label_color)
        settings.setValue("font_size", self.font_size)

        #screen_shape = QtGui.QDesktopWidget().screenGeometry()
        main_window = self.window()
        width = main_window.frameGeometry().width()
        height = main_window.frameGeometry().height()
        settings.setValue('screen_shape', (width, height))

        qpos = self.pos()
        pos = qpos.x(), qpos.y()
        settings.setValue('pos', pos)
        qApp.quit()
