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

        self._setup_supported_formats()

        logo = os.path.join(icon_path, 'logo.png')
        self.set_logo(logo)
        self.set_script_path(script_path)
        self.set_icon_path(icon_path)

        self.setup_gui()
        self.setup_post(inputs)

    def _setup_supported_formats(self):
        self.formats = {
            'cart3d' : is_cart3d,
        }
        for (name, is_on) in sorted(iteritems(self.formats)):
            if is_on:
                self.supported_formats.append(name)

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

    def on_load_geometry(self, infile_name=None, geometry_format=None):
        wildcard = ''
        is_failed = False

        if geometry_format and geometry_format.lower() not in self.supported_formats:
            is_failed = True
            if geometry_format in self.formats:
                msg = 'The import for the %r module failed.\n' % geometry_format
            else:
                msg = '%r is not a enabled format; enabled_formats=%s\n' % (geometry_format, self.supported_formats)
                msg += str("formats = %s" % str(self.formats))
            self.log_error(msg)
            return is_failed

        if infile_name:
            geometry_format = geometry_format.lower()
            print("geometry_format = %r" % geometry_format)
            if geometry_format == 'cart3d' and is_cart3d:
                has_results = True
                load_function = self.load_cart3d_geometry
            else:
                self.log_error('---invalid format=%r' % geometry_format)
                is_failed = True
                return is_failed
                raise NotImplementedError('on_load_geometry; infile_name=%r format=%r' % (infile_name, geometry_format))
            formats = [geometry_format]
            filter_index = 0
        else:
            formats = []
            load_functions = []
            has_results_list = []
            wildcard_list = []
            if is_cart3d:
                wildcard_list.append("Cart3d (*.tri; *.triq)")
                formats.append('Cart3d')
                has_results_list.append(True)
                load_functions.append(self.load_cart3d_geometry)

            wildcard = ';;'.join(wildcard_list)

            # get the filter index and filename
            if infile_name is not None and geometry_format is not None:
                filter_index = formats.index(geometry_format)
            else:
                Title = 'Choose a Geometry File to Load'
                wildcard_index, infile_name = self._create_load_file_dialog(wildcard, Title)
                if not infile_name:
                    is_failed = True
                    return is_failed # user clicked cancel
                filter_index = wildcard_list.index(wildcard_index)

            geometry_format = formats[filter_index]
            load_function = load_functions[filter_index]
            has_results = has_results_list[filter_index]

        if load_function is not None:
            self.last_dir = os.path.split(infile_name)[0]

            self.grid.Reset()
            self.grid.Modified()
            self.grid2.Reset()
            self.grid2.Modified()

            if not os.path.exists(infile_name) and geometry_format:
                msg = 'input file=%r does not exist' % infile_name
                self.log_error(msg)
                self.log_error(print_bad_path(infile_name))
                return

            if self.modelType is not None:
                # clear out old data
                name = 'clear_' + self.modelType

                # call the clear method
                try:
                    dy_method = getattr(self, 'clear_' + self.modelType)
                    dy_method()
                except:
                    print("method %r does not exist" % name)
            self.log_info("reading %s file %r" % (geometry_format, infile_name))
            try:
                has_results = load_function(infile_name, self.last_dir)
            except Exception as e:
                msg = traceback.format_exc()
                self.log_error(msg)
                raise
            self.rend.ResetCamera()

        # the model has been loaded, so we enable load_results
        if filter_index >= 0:
            self.format = formats[filter_index].lower()
            if has_results:
                enable = True
            else:
                enable = False
        else: # no file specified
            return
        self.infile_name = infile_name

        if self.out_filename is not None:
            msg = '%s - %s - %s' % (self.format, self.infile_name, self.out_filename)
        else:
            msg = '%s - %s' % (self.format, self.infile_name)
        self.set_window_title(msg)
        self.log_command("on_load_geometry(infile_name=%r, geometry_format=%r)" % (infile_name, self.format))

    def on_load_results(self, out_filename=None):
            geometry_format = self.format
            if self.format is None:
                msg ='on_load_results failed:  You need to load a file first...'
                self.log_error(msg)
                raise RuntimeError(msg)

            if out_filename in [None, False]:
                Title = 'Select a Results File for %s' % self.format
                wildcard = None
                if geometry_format == 'cart3d':
                    has_results = True
                    wildcard = "Cart3d (*.triq)"
                    load_functions = [self.load_cart3d_results]
                else:
                    msg = 'format=%r is not supported' % geometry_format
                    self.log_error(msg)
                    raise RuntimeError(msg)

                load_function = load_functions[0]
                if wildcard is None:
                    msg = 'format=%r has no method to load results' % geometry_format
                    self.log_error(msg)
                    return
                wildcard_index, out_filename = self._create_load_file_dialog(wildcard, Title)
            else:
                if geometry_format == 'cart3d':
                    load_function = self.load_cart3d_results
                else:
                    msg = 'format=%r is not supported.  Did you load a geometry model?' % geometry_format
                    self.log_error(msg)
                    raise RuntimeError(msg)

            if out_filename == '':
                return
            if not os.path.exists(out_filename):
                msg = 'result file=%r does not exist' % out_filename
                self.log_error(msg)
                return
            self.last_dir = os.path.split(out_filename)[0]
            load_function(out_filename, self.last_dir)

            self.out_filename = out_filename
            msg = '%s - %s - %s' % (self.format, self.infile_name, out_filename)
            self.set_window_title(msg)
            print("on_load_results(%r)" % out_filename)
            self.out_filename = out_filename
            self.log_command("on_load_results(%r)" % out_filename)

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
