import os
import sys

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from qtpy.QtGui import QMessageBox
from qtpy.QtWidgets import QApplication

import pyNastran
from pyNastran.gui.gui_common import GuiCommon2
from pyNastran.gui.arg_handling import get_inputs
#from pyNastran.gui.qt_files.gui_qt_common import GuiQtCommon
#from pyNastran.converters.nastran.nastran_io import NastranIO
from pyNastran.converters.stl.stl_io import STL_IO
from pyNastran.converters.openfoam.openfoam_io import OpenFoamIO
from pyNastran.converters.aflr.surf.surf_io import SurfIO
from pyNastran.converters.aflr.ugrid.ugrid_io import UGRID_IO
from pyNastran.converters.aflr.aflr2.bedge_io import BEdge_IO

try:
    PKG_PATH = sys._MEIPASS #@UndefinedVariable
    SCRIPT_PATH = os.path.join(PKG_PATH, 'scripts')
    ICON_PATH = os.path.join(PKG_PATH, 'icons')
except:
    PKG_PATH = pyNastran.__path__[0]
    SCRIPT_PATH = os.path.join(PKG_PATH, 'gui', 'scripts')
    ICON_PATH = os.path.join(PKG_PATH, 'gui', 'icons')


class MainWindow(GuiCommon2, STL_IO, OpenFoamIO, SurfIO, UGRID_IO, BEdge_IO): # NastranIO,
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

        logo = os.path.join(ICON_PATH, 'logo.png')

        self.logo = logo
        self.set_script_path(SCRIPT_PATH)
        self.set_icon_path(ICON_PATH)

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
        QMessageBox.about(self, "About pyNastran GUI", "\n".join(about))


    def on_reload(self):
        title = self.Title
        if 1:  # load the next file in the folder
            absdirname = os.path.abspath(os.path.dirname(self.infile_name))
            basename = os.path.basename(self.infile_name)
            ext = os.path.splitext(self.infile_name)[1]
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
        for unused_i in range(10):  #  limit on number of cycles
            if self.Title != title:
                self.cycleResults(title)
            else:
                break

def main():
    app = QApplication(sys.argv)
    QApplication.setOrganizationName("pyNastran")
    QApplication.setOrganizationDomain(pyNastran.__website__)
    QApplication.setApplicationName("pyNastran")
    QApplication.setApplicationVersion(pyNastran.__version__)

    inputs = get_inputs()
    window = MainWindow(inputs)
    window.startup(inputs)
    sys.exit(app.exec_())

if __name__ == '__main__':   # pragma: no cover
    main()
