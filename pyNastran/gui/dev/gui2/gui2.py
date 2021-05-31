import os
import sys
from typing import List, Dict, Optional, Any

#import ctypes
# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from cpylog import SimpleLogger
from cpylog.html_utils import str_to_html
import numpy as np
import vtk

import pyNastran
from qtpy import QtCore, QtGui #, API
from qtpy.QtWidgets import (
    QMainWindow, QFrame, QHBoxLayout, QAction, QMenu, QToolButton)
from qtpy.QtWidgets import QApplication
from pyNastran.gui.styles.trackball_style_camera import TrackballStyleCamera
from pyNastran.gui.menus.application_log import ApplicationLogWidget
from pyNastran.gui.menus.python_console import PythonConsoleWidget
from pyNastran.gui.gui_objects.settings import Settings

from pyNastran.gui.qt_files.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from pyNastran.gui.dev.gui2.utils import build_actions, fill_menus
from pyNastran.gui.dev.gui2.help_actions import HelpActions
from pyNastran.gui.dev.gui2.load_actions import LoadActions
#from pyNastran.gui.formats import CLASS_MAP
from pyNastran.gui.qt_files.view_actions import ViewActions
from pyNastran.gui.dev.gui2.vtk_interface import VtkInterface, ScalarBar


from pyNastran.gui.dev.gui2.format_setup import build_fmts, CLASS_MAP
PKG_PATH = pyNastran.__path__[0]

class MainWindow2(QMainWindow):
    """
    +-----------------------------------------+
    | menubar: File   Edit    View    Help    |
    +-----------------------------------------+
    | Toolbar                                 |
    +-------------------------------+---------+
    | VTK Window                    | Sidebar |
    +-------------------------------+---------+
    | Console / Logger                        |
    +-----------------------------------------+

    """
    def __init__(self):
        super().__init__()
        #self.setSize(500, 500)

        self.last_dir = ''
        self.is_gui = True
        self.dev = False
        self.debug = True
        self.html_logging = True
        self.execute_python = False
        self.performance_mode = False

        self.cases = {} # type: Dict[int, Any]
        self.form = []  # type: List[Any]
        # -----------------------------------------
        self.name = 'main'
        self.model_type = None
        self.nid_maps = {}
        self.eid_maps = {}
        self.models = {}  # type: Dict[str, Any]
        self.grid_mappers = {} # type: Dict[str, Any]
        self.main_grids = {} #  type: Dict[str, vtk.vtkUnstructuredGrid]
        self.alt_grids = {} # type: Dict[str, vtk.vtkUnstructuredGrid]
        self.geom_actors = {} # type: Dict[str, vtkLODActor]
        # -----------------------------------------
        self.settings = Settings(self)
        settings = QtCore.QSettings()
        self.settings.load(settings)

        self.actions = {}  # type: Dict[str, QAction]
        self.load_actions = LoadActions(self)
        self.view_actions = ViewActions(self)
        #self.log = SimpleLogger(level='debug', encoding='utf-8', nlevels=1, log_func=None)

        self.log = None
        self._start_logging()

        if self.html_logging is True:
            self.log_dock_widget = ApplicationLogWidget(self)
            self.log_widget = self.log_dock_widget.log_widget
            self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock_widget)
        else:
            self.log_widget = self.log

        #self.addToolBar
        self.toolbar = self.addToolBar('Show toolbar')
        self.toolbar.setObjectName('main_toolbar')
        self.toolbar.show()
        #action = QAction()
        #self.toolbar.addAction(QAction)

        self.menubar = self.menuBar()
        self._fill_menubar()
        self.statusBar().showMessage('Ready')
        self.show()

        self.format_class_map = CLASS_MAP
        fmt_order = ['cart3d', 'stl']
        self.fmts, self.supported_formats = build_fmts(
            self, self.format_class_map, fmt_order,
            self.log,
            stop_on_failure=False)
        #self.create_vtk_actors(create_rend=True)

        self.vtk_frame = QFrame()
        self.vtk_interactor = QVTKRenderWindowInteractor(parent=self.vtk_frame)
        self.set_style_as_trackball()

        self.rend = vtk.vtkRenderer()
        self.vtk_interactor.GetRenderWindow().AddRenderer(self.rend)
        self.build_vtk_frame()

        camera = self.rend.GetActiveCamera()
        if self.settings.use_parallel_projection:
            camera.ParallelProjectionOn()
        #else:
            #camera.ParallelProjectionOff()

        if self.execute_python:
            self.python_dock_widget = PythonConsoleWidget(self)
            self.python_dock_widget.setObjectName('python_console')
            self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.python_dock_widget)
        self._load_models()

    def _load_models(self) -> None:
        cart3d_filename = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\converters\cart3d\models\threePlugs.a.tri'
        #self.on_load_geometry()
        self.load_actions.on_load_geometry(
            infile_name=cart3d_filename, geometry_format='cart3d',
            name='cart3d', plot=True, raise_error=False)

        stl_filename = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\converters\stl\sphere.stl'
        self.load_actions.on_load_geometry(
            infile_name=stl_filename, geometry_format='stl',
            name='stl', plot=True, raise_error=False)

        # Render again to set the correct view
        self.render_window.Render()

    def _start_logging(self) -> None:
        if self.log is not None:
            return
        if self.html_logging is True:
            log = SimpleLogger(
                level='debug', encoding='utf-8',
                log_func=lambda w, x, y, z: self._logg_msg(w, x, y, z))
            # logging needs synchronizing, so the messages from different
            # threads would not be interleave
            self.log_mutex = QtCore.QReadWriteLock()
        else:
            log = SimpleLogger(
                level='debug', encoding='utf-8',
                #log_func=lambda x, y: print(x, y)  # no colorama
            )
        self.log = log

    def _logg_msg(self, log_type: str, filename: str, lineno: int, msg: str) -> None:
        """
        Add message to log widget trying to choose right color for it.

        Parameters
        ----------
        log_type : str
            {DEBUG, INFO, ERROR, COMMAND, WARNING} or prepend 'GUI '
        filename : str
            the active file
        lineno : int
            line number
        msg : str
            message to be displayed
        """
        if not self.html_logging:
            # standard logger
            name = '%-8s' % (log_type + ':')
            filename_n = '%s:%s' % (filename, lineno)
            msg2 = ' %-28s %s\n' % (filename_n, msg)
            print(name, msg2)
            return

        if 'DEBUG' in log_type and not self.settings.show_debug:
            return
        elif 'INFO' in log_type and not self.settings.show_info:
            return
        elif 'COMMAND' in log_type and not self.settings.show_command:
            return
        elif 'WARNING' in log_type and not self.settings.show_warning:
            return
        elif 'ERROR' in log_type and not self.settings.show_error:
            return

        if log_type in ['GUI ERROR', 'GUI COMMAND', 'GUI DEBUG', 'GUI INFO', 'GUI WARNING']:
            log_type = log_type[4:] # drop the GUI

        html_msg = str_to_html(log_type, filename, lineno, msg)

        if self.performance_mode or self.log_widget is None:
            self._log_messages.append(html_msg)
        else:
            self._log_msg(html_msg)

    def _log_msg(self, msg: str) -> None:
        """prints an HTML log message"""
        self.log_mutex.lockForWrite()
        text_cursor = self.log_widget.textCursor()
        end = text_cursor.End
        text_cursor.movePosition(end)
        text_cursor.insertHtml(msg)
        self.log_widget.ensureCursorVisible() # new message will be visible
        self.log_mutex.unlock()

    def log_info(self, msg: str) -> None:
        """ Helper funtion: log a message msg with a 'INFO:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        return self.log.simple_msg(msg, 'GUI INFO')

    def log_debug(self, msg: str) -> None:
        """ Helper funtion: log a message msg with a 'DEBUG:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        return self.log.simple_msg(msg, 'GUI DEBUG')

    def log_command(self, msg: str) -> None:
        """ Helper funtion: log a message msg with a 'COMMAND:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        return self.log.simple_msg(msg, 'GUI COMMAND')

    def log_error(self, msg: str) -> None:
        """ Helper funtion: log a message msg with a 'GUI ERROR:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        return self.log.simple_msg(msg, 'GUI ERROR')

    def log_warning(self, msg: str) -> None:
        """ Helper funtion: log a message msg with a 'WARNING:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        return self.log.simple_msg(msg, 'GUI WARNING')

    #def on_escape_null(self) -> None:
        #"""
        #The default state for Escape key is nothing.
        #"""
        #pass

    def _on_execute_python_button(self, clear=False):
        """executes the docked python console"""
        try:
            enter_data = self.python_dock_widget.enter_data
        except Exception as error:
            self.log_error(str(error))
            self.log_error('problem getting enter_data from python console')
            return
        txt = str(enter_data.toPlainText()).rstrip()
        is_passed = self._execute_python_code(txt)
        if is_passed and clear:
            enter_data.clear()

    def set_style_as_trackball(self):
        """sets the default rotation style"""
        #self._simulate_key_press('t') # change mouse style to trackball
        self.style = TrackballStyleCamera(self.vtk_interactor, self)
        self.vtk_interactor.SetInteractorStyle(self.style)

    def build_vtk_frame(self):
        """uses the vtk objects to set up the window (frame)"""
        vtk_hbox = QHBoxLayout()
        vtk_hbox.setContentsMargins(2, 2, 2, 2)

        vtk_hbox.addWidget(self.vtk_interactor)
        self.vtk_frame.setLayout(vtk_hbox)
        self.vtk_frame.setFrameStyle(QFrame.NoFrame | QFrame.Plain)
        # this is our main, 'central' widget
        self.setCentralWidget(self.vtk_frame)
        print('build_vtk_frame')

    @property
    def grid(self) -> vtk.vtkUnstructuredGrid:
        return self.main_grids[self.name]

    @property
    def iren(self) -> QVTKRenderWindowInteractor:
        return self.vtk_interactor
    @property
    def render_window(self):
        return self.vtk_interactor.GetRenderWindow()

    def render(self):
        self.vtk_interactor.GetRenderWindow().Render()

    def get_camera(self):
        return self.rend.GetActiveCamera()

    def turn_text_off(self) -> None:
        self.log.warning('turn_text_off')

    #-----------------------------------------------------------------------
    # geometry
    def set_quad_grid(self, box_name: str,
                      nodes: np.ndarray, elements: np.ndarray,
                      color: Optional[List[float]]=None,
                      line_width: float=1, opacity: float=1.) -> None:
        self.vtk_interface.set_quad_grid(box_name, nodes, elements,
                                         color=color, line_width=line_width, opacity=opacity)

    def create_global_axes(self, dim_max: float) -> None:
        self.vtk_interface.create_global_axes(dim_max)
    @property
    def scalar_bar_actor(self) -> ScalarBar:
        return self.vtk_interface.scalar_bar_actor
    # geometry
    #-----------------------------------------------------------------------
    # results post-processing
    def _finish_results_io2(self, model_name: str, form: List[Any], cases: Dict[int, Any]):
        self.form = form
        #self.cases = cases
        self.log.warning('_finish_results_io2')
    def get_new_icase(self) -> int:
        return 0
    def update_result_cases(self, cases: Dict[int, Any]) -> None:
        for case_id, case in cases.items():
            self.cases[case_id] = case
        return
    def get_form(self) -> List[Any]:
        return self.form

    #def _setup_formats(self):
        #fmt_name, _major_name, geom_wildcard, geom_func, res_wildcard, _res_func = fmt
        #from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
        #from pyNastran.converters.stl.stl_io import STL_IO
        #cart3d_class = Cart3dIO(self)
        #stl_class = STL_IO(self).get_stl_wildcard_geometry_results_functions()
        #fmts = []

        #return []

    #-----------------------------------------------------------------------
    # gui setup
    def _fill_menubar(self) -> None:
        file_actions_list = [
            'load_geometry', 'load_results', '',
            'load_custom_result', '',
            'load_csv_user_points', 'load_csv_user_geom', 'script', '', 'exit', ]

        self.vtk_interface = VtkInterface(self)
        help = HelpActions(self)
        toolbar_tools = [
            #'camera_reset', 'view',
            #'screenshot',
            'exit',
            #'min', 'max', 'map_element_fringe',
            '', # 'exit'
        ]
        menus_list = [
            ('file', '&File', file_actions_list),
            ('help', '&Help', help.actions_list),
            ('toolbar', self.toolbar, toolbar_tools),
        ]

        self.actions = self._setup_actions(help)  # type: Dict[str, QAction]
        #self.actions['pulldown'] =

        #self.combo = QtGui.QComboBox()
        #toolBar.addWidget(self.combo)
        #self.combo.insertItems(1,["One","Two","Three"])

        self.menus = fill_menus(self, menus_list, self.actions)

        #toolbar_tools = [
            #'reload', 'load_geometry', 'load_results',
            #'front_view', 'back_view', 'top_view', 'bottom_view',
            #'left_view', 'right_view',
            #'magnify', 'shrink', 'zoom',
            #'rotate_clockwise', 'rotate_cclockwise',
            #'rotation_center', 'measure_distance', 'probe_result',
            #'highlight_cell', 'highlight_node',
            #'area_pick', 'highlight_nodes_elements', 'mark_nodes_elements',
            #'wireframe', 'surface', 'edges']
        #toolbar_tools += [
            #'camera_reset', 'view', 'screenshot', 'min', 'max', 'map_element_fringe',
            #'', # 'exit'
        #]

    def _setup_actions(self, help: HelpActions) -> Dict[str, QAction]:
        icon_path = os.path.join(PKG_PATH, 'gui', 'icons')
        file_tools = [
            ('exit', '&Exit', 'texit.png', 'Ctrl+Q', 'Exit application', self.closeEvent),

            ('reload', 'Reload Model...', 'treload.png', '', 'Remove the model and reload the same geometry file', self.on_reload),
            ('load_geometry', 'Load &Geometry...', 'load_geometry.png', 'Ctrl+O', 'Loads a geometry input file', self.on_load_geometry),
            ('load_results', 'Load &Results...', 'load_results.png', 'Ctrl+R', 'Loads a results file', self.on_load_results),
            ('load_csv_user_geom', 'Load CSV User Geometry...', '', None, 'Loads custom geometry file', self.on_load_user_geom),
            ('load_csv_user_points', 'Load CSV User Points...', 'user_points.png', None, 'Loads CSV points', self.on_load_csv_points),
            ('load_custom_result', 'Load Custom Results...', '', None, 'Loads a custom results file', self.on_load_custom_results),

            ('script', 'Run Python Script...', 'python48.png', None, 'Runs pyNastranGUI in batch mode', self.on_run_script),
        ]
        checkables_set = set([])

        # setup the actions
        actions_list = file_tools + help.tools_list
        actions = build_actions(self, icon_path, actions_list, checkables_set, self.log)
        assert len(actions) > 0, actions
        return actions

    # ------------------------------------------
    # file
    def on_reload(self):
        pass

    def on_load_geometry(self):
        self.load_actions.on_load_geometry(
            infile_name=None, geometry_format=None,
            name='main', plot=True, raise_error=False)

    #def _reset_model(self, name: str) -> None:
        #self.log.info('_reset_model')

    def create_vtk_actors(self, create_rend: bool=True) -> None:
        """creates the vtk actors used by the GUI"""
        if create_rend:
            self.rend = vtk.vtkRenderer()

    @property
    def grid_selected(self):
        return self.main_grids[self.name]

    def _remove_old_geometry(self, filename: str):
        """
        >>> self.geom_actors
        {'cart3d': (vtkRenderingLODPython.vtkLODActor)000002B26C562C48,
         'stl': (vtkRenderingLODPython.vtkLODActor)000002B25024C7C8
        }
        """
        if filename in self.grid_mappers:
            #mapper = self.grid_mappers[filename]
            del self.grid_mappers[filename]
            grid = self.main_grids[filename]
            grid.FastDelete()
            del self.main_grids[filename]
            actor = self.geom_actors[filename]
            self.rend.RemoveActor(actor)
            del self.geom_actors[filename]
        #self.models = {}  # type: Dict[str, Any]
        #self.grid_mappers = {} # type: Dict[str, Any]
        #self.main_grids = {} #  type: Dict[str, vtk.vtkUnstructuredGrid]
        #self.alt_grids = {} # type: Dict[str, vtk.vtkUnstructuredGrid]
        #self.geom_actors = {} # type: Dict[str, vtkLODActor]

    def _reset_model(self, name: str) -> None:
        """resets the grids; sets up alt_grids"""
        if hasattr(self, 'main_grids') and name not in self.main_grids:
            grid = vtk.vtkUnstructuredGrid()
            grid_mapper = vtk.vtkDataSetMapper()
            grid_mapper.SetInputData(grid)

            geom_actor = vtk.vtkLODActor()
            geom_actor.DragableOff()
            geom_actor.SetMapper(grid_mapper)
            self.rend.AddActor(geom_actor)

            self.name = name
            self.models = {}
            self.main_grids[name] = grid
            self.grid_mappers[name] = grid_mapper
            self.geom_actors[name] = geom_actor
            grid.Modified()

            if 0:
                # link the current "main" to the scalar bar
                scalar_range = self.grid_selected.GetScalarRange()
                grid_mapper.ScalarVisibilityOn()
                grid_mapper.SetScalarRange(scalar_range)
                grid_mapper.SetLookupTable(self.color_function)

            #self.edge_actor = vtk.vtkLODActor()
            #self.edge_actor.DragableOff()
            #self.edge_mapper = vtk.vtkPolyDataMapper()

            # create the edges
            #self.get_edges()
        elif name in self.main_grids:
            grid = self.main_grids[name]
            grid.Reset()
            grid.Modified()
        else:
            self._setup_main_grid()

        # reset alt grids
        alt_grids = self.alt_grids
        alt_names = self.alt_grids.keys()
        for alt_name in alt_names:
            alt_grid = alt_grids[alt_name]
            alt_grid.Reset()
            alt_grid.Modified()

    def on_load_results(self):
        self.log.warning('on_load_results')
    def on_load_user_geom(self):
        self.log.warning('on_load_user_geom')
    def on_load_csv_points(self):
        self.log.warning('on_load_csv_points')
    def on_load_custom_results(self):
        self.log.warning('on_load_custom_results')
    def on_run_script(self):
        self.log.warning('on_run_script')
    # file
    # help
    # ------------------------------------------
    def _check_for_latest_version(self) -> bool:
        self.log.warning('_check_for_latest_version')
        return False
    # help
    # ------------------------------------------
    # basic functions
    @property
    def window_title(self) -> str:
        return self.getWindowTitle()

    @window_title.setter
    def window_title(self, msg: str) -> None:
        #msg2 = "%s - "  % self.base_window_title
        #msg2 += msg
        self.setWindowTitle(msg)

    def closeEvent(self, *args):
        """
        Handling saving state before application when application is
        being closed.
        """
        #settings = QtCore.QSettings()
        #settings.clear()
        #self.settings.save(settings)

        q_app = QApplication.instance()
        if q_app is None:
            sys.exit()
        q_app.quit()

def main():
    if sys.platform == 'win32':
        import ctypes
        myappid = 'pynastran.pynastrangui.%s' % (pyNastran.__version__) # arbitrary string
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    app = QApplication(sys.argv)
    QApplication.setOrganizationName('pyNastran')
    QApplication.setOrganizationDomain(pyNastran.__website__)
    QApplication.setApplicationName('pyNastran')
    QApplication.setApplicationVersion(pyNastran.__version__)

    w = MainWindow2()
    app.exec_()

if __name__ == '__main__':   # pragma: no cover
    main()
