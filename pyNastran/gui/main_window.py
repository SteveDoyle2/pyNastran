"""
defines the MainWindow class
"""
# coding: utf-8
# pylint: disable=C0111
from __future__ import division, unicode_literals, print_function

# standard library
import sys
import os.path
import imp
#import traceback
#import webbrowser
#webbrowser.open("http://xkcd.com/353/")


from pyNastran.gui.qt_version import qt_version
from qtpy import QtCore
from qtpy.QtWidgets import QMessageBox, qApp
from six.moves import urllib

# 3rd party
import vtk  # if this crashes, make sure you ran setup.py

# pyNastran
import pyNastran
from pyNastran.gui import SCRIPT_PATH, ICON_PATH
from pyNastran.gui.utils.version import check_for_newer_version
from pyNastran.gui.plugins import plugin_name_to_path
from pyNastran.gui.formats import NastranIO
from pyNastran.gui.gui_common import GuiCommon

# tcolorpick.png and tabout.png trefresh.png icons on LGPL license, see
# http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/color-picker-grey.png
# http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/help-hint.png
# http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/view-refresh-8.png


class MainWindow(GuiCommon, NastranIO):
    """
    The MainWindow class combines the base GuiCommon class with all the functionality
    with the remaining format holdout (Nastran).  It also defines which formats
    will be supported in the exe.

    MainWindow -> GuiCommon -> GuiQtCommon
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

        # these are in alphabetical order except for Nastran
        # this includes the bedge, surf, ugrid line (listed as AFLR in the gui)
        fmt_order = [
            # no results unless specified
            'nastran',  # results
            'abaqus',
            'avus',
            'bedge', 'surf', 'ugrid', 'ugrid3d', # aflr
            'cart3d',  # results
            'degen_geom',
            'fast',
            'lawgs',
            'obj',
            'openfoam_hex', 'openfoam_shell', 'openfoam_faces', # openfoam - results
            'panair',  # results
            'shabp',  # results
            'stl',
            'su2',
            'tecplot',  # results
            'tetgen',
            'usm3d',  # results
            'avl',
        ]
        #GuiCommon2.__init__(self, fmt_order, html_logging, inputs, parent)
        kwds['inputs'] = inputs
        kwds['fmt_order'] = fmt_order
        kwds['html_logging'] = html_logging
        super(MainWindow, self).__init__(**kwds)
        #fmt_order=fmt_order, inputs=inputs,
        #html_logging=html_logging,

        if qt_version in ['pyqt4', 'pyqt5', 'pyside', 'pyside2']:
            NastranIO.__init__(self)
        else:  # pragma: no cover
            raise NotImplementedError('qt_version=%r is not supported' % qt_version)

        self.build_fmts(fmt_order, stop_on_failure=False)

        self.logo = os.path.join(ICON_PATH, 'logo.png')
        self.set_script_path(SCRIPT_PATH)
        self.set_icon_path(ICON_PATH)

        is_gui = True
        if 'is_gui' in inputs:
            is_gui = inputs['is_gui']
            assert isinstance(is_gui, bool), is_gui
        self.start_logging()
        self._load_plugins()
        self.setup_gui(is_gui)
        self.setup_post(inputs)
        self._check_for_latest_version()

    def _load_plugins(self):
        """loads the plugins from pyNastran/gui/plugins.py"""
        for module_name, plugin_file, class_name in plugin_name_to_path:  # list
            if module_name in self.modules:
                raise RuntimeError('module_name=%r is already defined' % module_name)

            if not os.path.exists(plugin_file):
                # auto_wireframe is a test module and is not intended to
                # actually load unless you're testing
                if module_name != 'auto_wireframe':
                    self.log_warning('Failed to load plugin %r because %s doesnt exist' % (
                        module_name, plugin_file))
                continue

            module = imp.load_source(module_name, plugin_file)
            my_class = getattr(module, class_name)
            class_obj = my_class(self)
            self.modules[module_name] = class_obj

            # tools/checkables
            tools, checkables = class_obj.get_tools_checkables()
            self.tools += tools
            for key, is_active in checkables.items():
                self.checkables[key] = is_active

    def _check_for_latest_version(self, check=True):
        """
        checks the website for information regarding the latest gui version

        Looks for:
            ## pyNastran v0.7.2 has been Released (4/25/2015)
        """
        #import time
        #time0 = time.time()
        version_latest, unused_version_current, is_newer = check_for_newer_version()
        if is_newer and check:
            url = pyNastran.__website__
            from pyNastran.gui.menus.download import DownloadWindow
            win = DownloadWindow(url, version_latest, win_parent=self)
            win.show()
        #dt = time.time() - time0
        #print('dt_version_check = %.2f' % dt)

    def mousePressEvent(self, event):
        if not self.run_vtk:
            return
        #print('press x,y = (%s, %s)' % (ev.x(), ev.y()))
        if self.is_pick:
            #self.___saveX = ev.x()
            #self.___saveY = ev.y()
            pass
        else:
            self.vtk_interactor.mousePressEvent(event)

    #def LeftButtonPressEvent(self, ev):

    def mouseReleaseEvent(self, event):
        #print('release x,y = (%s, %s)' % (ev.x(), ev.y()))
        if self.is_pick:
            pass
        else:
            self.vtk_interactor.mousePressEvent(event)

    def open_website(self):
        """loads the pyNastran main website"""
        self._urlopen(pyNastran.__website__)

    def open_docs(self):
        """loads the pyNastran docs website"""
        url = pyNastran.__docs__
        try:
            urllib.request.urlopen(url)
        except (urllib.error.HTTPError, urllib.error.URLError):
            url = pyNastran.__docs_rtd__
        self._urlopen(url)

    def open_issue(self):
        """loads the pyNastran issue tracker"""
        self._urlopen(pyNastran.__issue__)

    def open_discussion_forum(self):
        """loads the pyNastran discussion forum website"""
        self._urlopen(pyNastran.__discussion_forum__)

    def _urlopen(self, url):
        """opens a URL"""
        import webbrowser
        if self.is_gui:
            webbrowser.open(url)

    def about_dialog(self):
        """Display about dialog"""
        copyright = pyNastran.__copyright__
        if qt_version in ['pyside', 'pyside2']:
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

        Reload allows you to edit the input model and "reload" the data
        without having to go to the pulldown menu.  If you don't like
        this behavior, implement the self.on_reload_nastran() or similar
        method for a given format.
        """
        camera = self.get_camera_data()
        unused_title = self.title
        case = self.icase

        on_reload_name = 'on_reload_%s' % self.format
        if hasattr(self, on_reload_name):
            getattr(self, on_reload_name)()  # on_reload_nastran
        else:
            self.on_load_geometry(self.infile_name, self.format, raise_error=False)

        if self.out_filename is None:
            msg = '%s - %s' % (self.format, self.infile_name)
        else:
            msg = '%s - %s - %s' % (self.format, self.infile_name, self.out_filename)
        self.window_title = msg
        self.log_command('on_reload()')
        self.cycle_results(case)
        self.on_set_camera_data(camera, show_log=False)

    def closeEvent(self, *args):
        """
        Handling saving state before application when application is
        being closed.
        """
        settings = QtCore.QSettings()
        settings.clear()
        self.settings.save(settings)

        if qApp is None:
            sys.exit()
        qApp.quit()
