import os
import sys
import traceback
from PyQt4 import QtCore, QtGui

import pyNastran
from pyNastran.gui.gui_common import GuiCommon2
from pyNastran.gui.arg_handling import get_inputs
from pyNastran.gui.qt_files.gui_qt_common import GuiCommon
#from pyNastran.converters.nastran.nastranIOv import NastranIO
from pyNastran.converters.stl.stl_io import STL_IO
from pyNastran.utils import print_bad_path


from pyNastran.converters.openfoam.openfoamIO import OpenFoamIO
from pyNastran.converters.ugrid.surf_io import SurfIO
from pyNastran.converters.ugrid.ugrid_io import UGRID_IO
from pyNastran.converters.aflr2.bedge_io import BEDGE_IO

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

try:
    pkg_path = sys._MEIPASS #@UndefinedVariable
    script_path = os.path.join(pkg_path, 'scripts')
    icon_path = os.path.join(pkg_path, 'icons')
except:
    pkg_path = pyNastran.__path__[0]
    script_path = os.path.join(pkg_path, 'gui', 'scripts')
    icon_path = os.path.join(pkg_path, 'gui', 'icons')


class MainWindow(GuiCommon2, STL_IO, OpenFoamIO, SurfIO, UGRID_IO, BEDGE_IO): # NastranIO,
    def __init__(self, inputs, **kwds):
        self.mesh_3d = None

        fmt_order = ['stl', 'openfoam_hex', 'openfoam_shell', 'openfoam_faces',
                     'surf', 'ugrid', 'bedge']
        html_logging = True
        kwds['inputs'] = inputs
        kwds['fmt_order'] = fmt_order
        kwds['html_logging'] = html_logging
        super(MainWindow, self).__init__(**kwds)
        STL_IO.__init__(self)
        OpenFoamIO.__init__(self)
        SurfIO.__init__(self)

    def startup(self, inputs):
        self.build_fmts(self.fmts, stop_on_failure=True)

        logo = os.path.join(icon_path, 'logo.png')

        self.logo = logo
        self.set_script_path(script_path)
        self.set_icon_path(icon_path)
        print('gui', self.supported_formats)

        self.setup_gui()
        self.menu_help2 = self.menubar.addMenu('&Help')

        self.setup_post(inputs)

    #def init_cell_picker(self):
        #pass


    def about_dialog(self):
        """ Display about dialog """
        copyright = pyNastran.__pyqt_copyright__

        about = [
            'pyNastran Qt GUI',
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
        Title = self.Title
        if 1:  # load the next file in the folder
            absdirname = os.path.abspath(os.path.dirname(self.infile_name))
            basename = os.path.basename(self.infile_name)
            base, ext = os.path.splitext(self.infile_name)
            fnames = os.listdir(absdirname)
            #print(fnames)
            fnames = [os.path.join(absdirname, fname) for fname in fnames
                      if fname.endswith(ext)]
            #print(fnames)

            fnames *= 2  # avoids dealing with index issues
            old_fname = os.path.join(absdirname, basename)
            ifile = fnames.index(old_fname)
            fname = fnames[ifile + 1]
            try:
                self.on_load_geometry(fname, self.format)
            except:
                self.infile_name = fname

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

def main():
    app = QtGui.QApplication(sys.argv)
    QtGui.QApplication.setOrganizationName("pyNastran")
    QtGui.QApplication.setOrganizationDomain(pyNastran.__website__)
    QtGui.QApplication.setApplicationName("pyNastran")
    QtGui.QApplication.setApplicationVersion(pyNastran.__version__)

    inputs = get_inputs()
    window = MainWindow(inputs)
    window.startup(inputs)
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
