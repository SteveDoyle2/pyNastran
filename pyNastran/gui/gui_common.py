# -*- coding: utf-8 -*-
from __future__ import division, unicode_literals, print_function
from six import string_types, iteritems
from six.moves import range

# standard library
import sys
import os.path
import datetime
import cgi #  html lib
import inspect
import traceback

import vtk
from PyQt4 import QtCore, QtGui
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from numpy import eye

import pyNastran
from pyNastran.gui.qt_files.gui_qt_common import GuiCommon
from pyNastran.gui.qt_files.qt_legend import LegendPropertiesWindow
from pyNastran.utils.log import SimpleLogger
from pyNastran.gui.ex_tree import Sidebar
from pyNastran.utils import print_bad_path


class GuiCommon2(QtGui.QMainWindow, GuiCommon):
    def __init__(self, html_logging, inputs):
        QtGui.QMainWindow.__init__(self)
        GuiCommon.__init__(self)

        self.html_logging = html_logging
        self._logo = None
        self._script_path = None
        self._icon_path = ''

        self.Title = None
        self.min_value = None
        self.max_value = None
        self.blue_to_red = False
        self._is_axes_shown = True
        self.nvalues = 9
        self.is_wireframe = False

        #-------------
        # inputs dict
        self.is_edges = inputs['is_edges']
        self.is_nodal = inputs['is_nodal']
        self.is_centroidal = inputs['is_centroidal']
        self.magnify = inputs['magnify']
        assert self.is_centroidal != self.is_nodal, "is_centroidal and is_nodal can't be the same and are set to \"%s\"" % self.is_nodal

        #self.format = ''
        debug = inputs['debug']
        self.debug = debug
        assert debug in [True, False], 'debug=%s' % debug

        #-------------
        # file
        self.format = None
        self.infile_name = None
        self.out_filename = None
        self.dirname = ''
        self.last_dir = '' # last visited directory while opening file

        #-------------
        # internal params
        self.show_info = True
        self.show_debug = True
        self.show_gui = True
        self.show_command = True
        self.coord_id = 0

        self.nCases = 0
        self.iCase = 0
        self.nNodes = 0
        self.nElements = 0

        self.base_window_title = "pyNastran v%s"  % pyNastran.__version__
        self.supported_formats = []
        self.modelType = None

        self.tools = []
        self.checkables = []

        # initializes tools/checkables
        self.set_tools()

    def set_window_title(self, msg):
        #msg2 = "%s - "  % self.base_window_title
        #msg2 += msg
        self.setWindowTitle(msg)

    def set_logo(self, logo):
        self._logo = logo

    def init_ui(self):
        """ Initialize user iterface"""
        self.resize(800, 600)
        self.statusBar().showMessage('Ready')

        # windows title and aplication icon
        self.setWindowTitle('Statusbar')
        if self._logo is not None:
            self.setWindowIcon(QtGui.QIcon(self._logo))
        self.set_window_title(self.base_window_title)

        #=========== Results widget ===================
        self.res_dock = QtGui.QDockWidget("Results", self)
        self.res_dock.setObjectName("results_obj")
        #self.res_widget = QtGui.QTextEdit()
        #self.res_widget.setReadOnly(True)
        #self.res_dock.setWidget(self.res_widget)

        self.res_widget = Sidebar(self)
        self.res_widget.clear_data()
        #self.res_widget.update_results(data)

        #self.res_widget.setWidget(sidebar)

        self.res_dock.setWidget(self.res_widget)

        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.res_dock)
        #=========== Logging widget ===================
        if self.html_logging:
            self.log_dock = QtGui.QDockWidget("Application log", self)
            self.log_dock.setObjectName("application_log")
            self.log_widget = QtGui.QTextEdit()
            self.log_widget.setReadOnly(True)
            self.log_dock.setWidget(self.log_widget)
            self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock)
        #===============================================

        self._build_menubar()

        # right sidebar
        self.res_dock.hide()

        self.build_vtk_frame()

    def load_batch_inputs(self, inputs):
        if not inputs['format']:
            return
        form = inputs['format'].lower()
        input = inputs['input']
        output = inputs['output']
        plot = True
        if output:
            plot = False

        is_failed = self.on_load_geometry(input, form, plot=plot)
        if is_failed:
            return
        if output:
            self.on_load_results(output)
        self._simulate_key_press('r')
        self.vtk_interactor.Modified()

    def set_script_path(self, script_path):
        self._script_path = script_path

    def set_icon_path(self, icon_path):
        self._icon_path = icon_path

    def set_tools(self, tools=None, checkables=None):
        if checkables is None:
            checkables = ['show_info', 'show_debug', 'show_gui', 'show_command']
        if tools is None:
            tools = [
                ('exit', '&Exit', 'texit.png', 'Ctrl+Q', 'Exit application', self.closeEvent), # QtGui.qApp.quit
                ('load_geometry', 'Load &Geometry', 'load_geometry.png', 'Ctrl+O', 'Loads a geometry input file', self.on_load_geometry),  ## @todo no picture...
                ('load_results', 'Load &Results', 'load_results.png', 'Ctrl+R', 'Loads a results file', self.on_load_results),  ## @todo no picture...
                ('back_col', 'Change background color', 'tcolorpick.png', None, 'Choose a background color', self.change_background_col),
                ('legend', 'Modify legend', 'legend.png', None, 'Set Legend', self.set_legend),
                ('axis', 'Show/Hide Axis', 'axis.png', None, 'Show/Hide Global Axis', self.on_show_hide_axes),

                ('wireframe', 'Wireframe Model', 'twireframe.png', 'w', 'Show Model as a Wireframe Model', self.on_wireframe),
                ('surface', 'Surface Model', 'tsolid.png', 's', 'Show Model as a Surface Model', self.on_surface),
                ('edges', 'Show/Hide Edges', 'tedges.png', 'e', 'Show/Hide Model Edges', self.on_flip_edges),

                ('show_info', 'Show INFO', 'show_info.png', None, 'Show "INFO" messages', self.on_show_info),
                ('show_debug', 'Show DEBUG', 'show_debug.png', None, 'Show "DEBUG" messages', self.on_show_debug),
                ('show_gui', 'Show GUI', 'show_gui.png', None, 'Show "GUI" messages', self.on_show_gui),
                ('show_command', 'Show COMMAND', 'show_command.png', None, 'Show "COMMAND" messages', self.on_show_command),

                ('magnify', 'Magnify', 'plus_zoom.png', 'M', 'Increase Magnfication', self.on_increase_magnification),
                ('shrink', 'Shrink', 'minus_zoom.png', 'm', 'Decrease Magnfication', self.on_decrease_magnification),

                ('cell_pick', 'Cell Pick', '', 'CTRL+K', 'PickTip', self.on_cell_picker),

                ('rotate_clockwise', 'Rotate Clockwise', 'tclock.png', 'o', 'Rotate Clockwise', self.on_rotate_clockwise),
                ('rotate_cclockwise', 'Rotate Counter-Clockwise', 'tcclock.png', 'O', 'Rotate Counter-Clockwise', self.on_rotate_cclockwise),

                ('scshot', 'Take a Screenshot', 'tcamera.png', 'CTRL+I', 'Take a Screenshot of current view', self.take_screenshot),
                ('about', 'About pyNastran GUI', 'tabout.png', 'CTRL+H', 'About pyCart3d GUI and help on shortcuts', self.about_dialog),
                ('creset', 'Reset camera view', 'trefresh.png', 'r', 'Reset the camera view to default', self.on_reset_camera),
                ('reload', 'Reload model', 'treload.png', 'r', 'Reload the model', self.on_reload),

                ('cycle_res', 'Cycle Results', 'cycle_results.png', 'CTRL+L', 'Changes the result case', self.cycleResults),

                ('x', 'Flips to +X Axis', 'plus_x.png', 'x', 'Flips to +X Axis', lambda: self.update_camera('+x')),
                ('y', 'Flips to +Y Axis', 'plus_y.png', 'y', 'Flips to +Y Axis', lambda: self.update_camera('+y')),
                ('z', 'Flips to +Z Axis', 'plus_z.png', 'z', 'Flips to +Z Axis', lambda: self.update_camera('+z')),

                ('X', 'Flips to -X Axis', 'minus_x.png', 'X', 'Flips to -X Axis', lambda: self.update_camera('-x')),
                ('Y', 'Flips to -Y Axis', 'minus_y.png', 'Y', 'Flips to -Y Axis', lambda: self.update_camera('-y')),
                ('Z', 'Flips to -Z Axis', 'minus_z.png', 'Z', 'Flips to -Z Axis', lambda: self.update_camera('-z')),
                ('script', 'Run Python script', 'python48.png', None, 'Runs pyCart3dGUI in batch mode', self.on_run_script),
            ]
        self.tools = tools
        self.checkables = checkables

    def add_tools(self, tools):
        self.tools += tools

    def _create_menu_items(self, actions):
        self.menu_file = self.menubar.addMenu('&File')
        self.menu_view = self.menubar.addMenu('&View')
        self.menu_window = self.menubar.addMenu('&Window')
        self.menu_help = self.menubar.addMenu('&Help')

        if self._script_path is not None and os.path.exists(self._script_path):
            scripts = [script for script in os.listdir(self._script_path) if '.py' in script]
        else:
            scripts = []

        scripts = tuple(scripts)

        #if 0:
            #print('script_path =', script_path)
            #print('scripts =', scripts)
            #self.menu_scripts = self.menubar.addMenu('&Scripts')
            #for script in scripts:
                #fname = os.path.join(script_path, script)
                #tool = (script, script, 'python48.png', None, '', lambda: self.on_run_script(fname) )
                #tools.append(tool)
        #else:
        self.menu_scripts = None

        menu_window = ['toolbar', 'reswidget']
        menu_view = ['scshot', '', 'wireframe', 'surface', 'creset', '', 'back_col', 'legend','axis', ]
        if self.html_logging:
            actions['logwidget'] = self.log_dock.toggleViewAction()
            actions['logwidget'].setStatusTip("Show/Hide application log")
            menu_view += ['', 'show_info', 'show_debug', 'show_gui', 'show_command']
            menu_window += ['logwidget']

        menu_items = [
            (self.menu_file, ('load_geometry', 'load_results', 'script', '', 'exit')),
            (self.menu_view, tuple(menu_view)),
            (self.menu_window, tuple(menu_window)),
            (self.menu_help, ('about',)),
            (self.menu_scripts, scripts),
            (self.toolbar, ('cell_pick', 'reload', 'load_geometry', 'load_results', 'cycle_res',
                            'x', 'y', 'z', 'X', 'Y', 'Z',
                            'magnify', 'shrink', 'rotate_clockwise', 'rotate_cclockwise',
                            'wireframe', 'surface', 'edges', 'creset', 'scshot', '', 'exit'))
        ]
        return menu_items

    def _build_menubar(self):
        ## toolbar
        self.toolbar = self.addToolBar('Show toolbar')
        self.toolbar.setObjectName('main_toolbar')

        ## menubar
        self.menubar = self.menuBar()

        actions = self._prepare_actions(self._icon_path, self.tools, self.checkables)
        menu_items = self._create_menu_items(actions)
        self._populate_menu(menu_items, actions)

    def _populate_menu(self, menu_items, actions):
        """populate menus and toolbar"""
        for menu, items in menu_items:
            if menu is None:
                continue
            for i in items:
                if not i:
                    menu.addSeparator()
                else:
                    if not isinstance(i, string_types):
                        raise RuntimeError('what is this...action i() = %r' % i())

                    action = actions[i] #if isinstance(i, string_types) else i()
                    menu.addAction(action)

    def _prepare_actions(self, icon_path, tools, checkables):
        """
        Prepare actions that will  be used in application in a way
        that's independent of the  menus & toolbar
        """
        actions = {}
        for (nam, txt, icon, shortcut, tip, func) in tools:
            #print("name=%s txt=%s icon=%s short=%s tip=%s func=%s" % (nam, txt, icon, short, tip, func))
            #if icon is None:
                #print("missing_icon = %r!!!" % nam)
                #icon = os.path.join(icon_path, 'no.png')

            if icon is None:
                print("missing_icon = %r!!!" % nam)
                ico = None
                #print(print_bad_path(icon))
            #elif not "/" in icon:
                #ico = QtGui.QIcon.fromTheme(icon)
            else:
                ico = QtGui.QIcon()
                pth = os.path.join(icon_path, icon)
                ico.addPixmap(QtGui.QPixmap(pth), QtGui.QIcon.Normal, QtGui.QIcon.Off)

            if nam in checkables:
                actions[nam] = QtGui.QAction(ico, txt, self, checkable=True)
                actions[nam].setChecked(True)
            else:
                actions[nam] = QtGui.QAction(ico, txt, self)

            if shortcut:
                actions[nam].setShortcut(shortcut)
            if tip:
                actions[nam].setStatusTip(tip)
            if func:
                actions[nam].triggered.connect(func)

        actions['toolbar'] = self.toolbar.toggleViewAction()
        actions['toolbar'].setStatusTip("Show/Hide application toolbar")

        actions['reswidget'] = self.res_dock.toggleViewAction()
        actions['reswidget'].setStatusTip("Show/Hide results selection")
        return actions

    def logg_msg(self, typ, msg):
        """
        Add message to log widget trying to choose right color for it.

        :param msg: message to be displayed
        """
        if not self.html_logging:
            print(typ, msg)
            return
        _fr = sys._getframe(4)  # jump to get out of the logger code
        n = _fr.f_lineno
        fn = os.path.basename(_fr.f_globals['__file__'])

        if typ == 'DEBUG' and not self.show_debug:
            return
        elif typ == 'INFO' and not self.show_info:
            return
        elif typ == 'GUI' and not self.show_gui:
            return
        elif typ == 'COMMAND' and not self.show_command:
            return

        if typ in ['GUI', 'COMMAND']:
            msg = '   fname=%-25s lineNo=%-4s   %s\n' % (fn, n, msg)

        tim = datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
        msg = cgi.escape(msg)
        #message colors
        dark_orange = '#EB9100'
        cols = {"GUI": "blue", "COMMAND":"green", "GUI ERROR":"Crimson", "DEBUG" : dark_orange}
        msg = msg.rstrip().replace('\n', '<br>')
        msg = tim + ' ' + (typ + ': ' + msg) if typ else msg
        if typ in cols:
            msg = '<font color="%s"> %s </font>' % (cols[typ], msg)

        self.log_mutex.lockForWrite()
        text_cursor = self.log_widget.textCursor()
        end = text_cursor.End
        #print("end", end)
        text_cursor.movePosition(end)
        #print(dir(text_cursor))
        text_cursor.insertHtml(msg + r"<br />")
        self.log_widget.ensureCursorVisible() # new message will be visible
        self.log_mutex.unlock()

    def log_info(self, msg):
        """ Helper funtion: log a messaage msg with a 'GUI:' prefix """
        self.log.simple_msg(msg, 'GUI')

    def log_command(self, msg):
        """ Helper funtion: log a messaage msg with a 'GUI:' prefix """
        self.log.simple_msg(msg, 'COMMAND')

    def log_error(self, msg):
        """ Helper funtion: log a messaage msg with a 'GUI:' prefix """
        self.log.simple_msg(msg, 'GUI ERROR')

    def change_background_col(self):
        """ Choose a background color """
        c = [int(255 * i) for i in self.background_col]
        col = QtGui.QColorDialog.getColor(QtGui.QColor(*c), self, "Choose a background color")
        if col.isValid():
            self.background_col = col.getRgbF()[:3]
            self.rend.SetBackground(*self.background_col)

    def create_coordinate_system(self, label='', origin=None, matrix_3x3=None, Type='xyz'):
        """
        Creates a coordinate system

        :param label:
          the coord id or other unique label (default is empty to indicate the global frame)
        :param origin:
          the origin as (3,) ndarray/list/tuple
        :param matrix_3x3:
          a standard 3x3 Nastran-style coordinate system
        :param Type:
          a string of 'xyz', 'Rtz', 'Rtp' (xyz, cylindrical, spherical)
          that changes the axis names

        .. todo::  Type is not supported ('xyz' ONLY)
        .. todo::  Can only set one coordinate system

        .. seealso::
            http://en.wikipedia.org/wiki/Homogeneous_coordinates
            http://www3.cs.stonybrook.edu/~qin/courses/graphics/camera-coordinate-system.pdf
            http://www.vtk.org/doc/nightly/html/classvtkTransform.html#ad58b847446d791391e32441b98eff151
        """
        coord_id = self.coord_id

        transform = vtk.vtkTransform()
        if origin is None and matrix_3x3 is None:
            pass
        elif origin is not None and matrix_3x3 is None:
            print('origin%s = %s' % (label, str(origin)))
            transform.Translate(*origin)
        elif matrix_3x3 is not None:  # origin can be None
            m = eye(4, dtype='float32')
            m[:3, :3] = matrix_3x3
            if origin is not None:
                m[:3, 3] = origin
            transform.SetMatrix(m.ravel())
        else:
            raise RuntimeError('unexpected coordinate system')

        axes = vtk.vtkAxesActor()
        axes.SetUserTransform(transform)

        if Type == 'xyz':
            if label:
                xlabel = 'x%s' % label
                ylabel = 'y%s' % label
                zlabel = 'z%s' % label
                axes.SetXAxisLabelText(xlabel)
                axes.SetYAxisLabelText(ylabel)
                axes.SetZAxisLabelText(zlabel)
        else:
            if Type == 'Rtz':  # cylindrical
                #x = u'R'
                #y = u'?'
                #z = 'z'
                x = 'R'
                y = 't'
                z = 'z'

            elif Type == 'Rtp':  # spherical
                xlabel = u'R'
                #ylabel = u'?'
                #z = u'?'
                x = 'R'
                y = 't'
                z = 'p'
            else:
                raise RuntimeError('invalid axis type; Type=%r' % Type)

            xlabel = '%s%s' % (x, label)
            ylabel = '%s%s' % (y, label)
            zlabel = '%s%s' % (z, label)
            axes.SetXAxisLabelText(xlabel)
            axes.SetYAxisLabelText(ylabel)
            axes.SetZAxisLabelText(zlabel)

        self.transform[coord_id] = transform
        self.axes[coord_id] = axes
        self.coord_id += 1
        self.rend.AddActor(axes)
        return self.coord_id

    def create_global_axes(self):
        self.transform = {}
        self.axes = {}
        #self.create_coordinate_system(origin=None, matrix_3x3=None, Type='Rtp')
        self.create_coordinate_system(label='', origin=None, matrix_3x3=None, Type='xyz')

    def on_show_hide_axes(self):
        # this method should handle all the coords when
        # there are more then one
        if self._is_axes_shown:
            for key, axis in iteritems(self.axes):
                axis.VisibilityOff()
        else:
            for key, axis in iteritems(self.axes):
                axis.VisibilityOn()
        self._is_axes_shown = not(self._is_axes_shown)

    def create_vtk_actors(self):
        self.rend = vtk.vtkRenderer()

        # vtk actors
        self.grid = vtk.vtkUnstructuredGrid()
        #gridResult = vtk.vtkFloatArray()

        self.grid2 = vtk.vtkUnstructuredGrid()
        #self.emptyResult = vtk.vtkFloatArray()
        #self.vectorResult = vtk.vtkFloatArray()

        # edges
        self.edgeActor = vtk.vtkActor()
        self.edgeMapper = vtk.vtkPolyDataMapper()

        self.create_cell_picker()

        # scalar bar
        self.scalarBar = vtk.vtkScalarBarActor()
        self.create_global_axes()

    def build_vtk_frame(self):
        #Frame that VTK will render on
        vtk_frame = QtGui.QFrame()
        vtk_hbox = QtGui.QHBoxLayout()
        vtk_hbox.setContentsMargins(2, 2, 2, 2)

        #Qt VTK RenderWindowInteractor
        self.vtk_interactor = QVTKRenderWindowInteractor(parent=vtk_frame)
        self.iren = self.vtk_interactor
        vtk_hbox.addWidget(self.vtk_interactor)
        vtk_frame.setLayout(vtk_hbox)
        vtk_frame.setFrameStyle(QtGui.QFrame.NoFrame | QtGui.QFrame.Plain)
        # this is our main, 'central' widget
        self.setCentralWidget(vtk_frame)

        #=============================================================
        self.vtk_interactor.GetRenderWindow().AddRenderer(self.rend)
        self.vtk_interactor.GetRenderWindow().Render()
        #self.load_nastran_geometry(None, None)
        self.textActors = {}


        #for cid, axes in iteritems(self.axes):
            #self.rend.AddActor(axes)
        self.addGeometry()
        self.addAltGeometry()
        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(*self.background_col)
        self.rend.ResetCamera()
        self._simulate_key_press('t') # change mouse style to trackball
        self.build_lookup_table()

        self.magnify = 1
        self.iText = 0
        textSize = 14 * self.magnify
        self.createText([5, 50], 'Max  ', textSize)  # text actor 0
        self.createText([5, 35], 'Min  ', textSize)  # text actor 1
        self.createText([5, 20], 'Word1', textSize)  # text actor 2
        self.createText([5, 5], 'Word2', textSize)  # text actor 3

        self.get_edges()
        if self.is_edges:
            prop = self.edgeActor.GetProperty()
            prop.EdgeVisibilityOn()
        else:
            prop = self.edgeActor.GetProperty()
            prop.EdgeVisibilityOff()

    #def _script_helper(self, python_file=False):
        #if python_file in [None, False]:
            #self.on_run_script(python_file)

    def set_legend(self):
        """
        Opens a dialog box to set:

        +--------+----------+
        |  Name  |  String  |
        +--------+----------+
        |  Min   |  Float   |
        +--------+----------+
        |  Max   |  Float   |
        +--------+----------+
        | Format | pyString |
        +--------+----------+
        """
        if not hasattr(self, 'caseKeys'):
            self.log_error('No model has been loaded.')
            return
        key = self.caseKeys[self.iCase]
        case = self.resultCases[key]
        if len(key) == 5:
            (subcaseID, resultType, vectorSize, location, data_format) = key
        elif len(key) == 6:
            (subcaseID, i, resultType, vectorSize, location, data_format) = key
        else:
            (subcaseID, i, resultType, vectorSize, location, data_format, label2) = key

        data = {
            'name' : resultType,
            'min' : case.min(),
            'max' : case.max(),
            'format' : data_format,
            'is_blue_to_red' : True,
            'is_discrete': True,
            'clicked_ok' : False,
        }
        legend = LegendPropertiesWindow(data, win_parent=self)
        legend.show()
        legend.exec_()

        if data['clicked_ok']:
            self.apply_legend(data)

    def apply_legend(self, data):
        Title = data['name']
        min_value = data['min']
        max_value = data['max']
        data_format = data['format']
        is_blue_to_red = data['is_blue_to_red']
        is_discrete = data['is_discrete']
        self.on_update_legend(Title=Title, min_value=min_value, max_value=max_value,
                              data_format=data_format,
                              is_blue_to_red=is_blue_to_red,
                              is_discrete=is_discrete)

    def on_update_legend(self, Title='Title', min_value=0., max_value=1.,
                         data_format='%.0f', is_blue_to_red=True, is_discrete=True):
        key = self.caseKeys[self.iCase]
        case = self.resultCases[key]
        if len(key) == 5:
            (subcase_id, _result_type, vector_size, location, _data_format) = key
        elif len(key) == 6:
            (subcase_id, i, _result_type, vector_size, location, _data_format) = key
        else:
            (subcase_id, i, _result_type, vector_size, location, _data_format, label2) = key

        try:
            caseName = self.iSubcaseNameMap[subcase_id]
        except KeyError:
            caseName = ('case=NA', 'label=NA')
        (subtitle, label) = caseName

        gridResult = self.build_grid_result(vector_size, location)
        norm_value, nValueSet = self.set_grid_values(gridResult, case, vector_size, min_value, max_value, is_blue_to_red=is_blue_to_red)
        self.UpdateScalarBar(Title, min_value, max_value, norm_value, data_format, is_blue_to_red=is_blue_to_red)
        self.final_grid_update(gridResult, key, subtitle, label)
        self.log_command('self.on_update_legend(Title=%r, min_value=%s, max_value=%s,\n'
                         '                      data_format=%r, is_blue_to_red=%s, is_discrete=%s)'
                         % (Title, min_value, max_value, data_format, is_blue_to_red, is_discrete))

    def on_run_script(self, python_file=False):
        print('python_file =', python_file)
        if python_file in [None, False]:
            Title = 'Choose a Python Script to Run'
            wildcard = "Python (*.py)"
            wildcard_index, infile_name = self._create_load_file_dialog(wildcard, Title)
            if not infile_name:
                is_failed = True
                return is_failed # user clicked cancel

            python_file = os.path.join(script_path, infile_name)
        execfile(python_file)
        self.log_command('self.on_run_script(%r)' % python_file)

    def on_show_info(self):
        self.show_info = not self.show_info

    def on_show_debug(self):
        self.show_debug = not self.show_debug

    def on_show_gui(self):
        self.show_gui = not self.show_gui

    def on_show_command(self):
        self.show_command = not self.show_command

    def on_reset_camera(self):
        self.log_command('on_reset_camera()')
        self._simulate_key_press('r')

    def on_surface(self):
        if self.is_wireframe:
            self.log_command('on_surface()')
            self._simulate_key_press('s')
            self.is_wireframe = False

    def on_wireframe(self):
        if not self.is_wireframe:
            self.log_command('on_wireframe()')
            self._simulate_key_press('w')
            self.is_wireframe = True

    def _update_camera(self, camera=None):
        if camera is None:
            camera = self.GetCamera()
        camera.Modified()
        self.vtk_interactor.Render()

    def zoom(self, value):
        camera = self.GetCamera()
        camera.Zoom(value)
        camera.Modified()
        self.vtk_interactor.Render()
        self.log_command('zoom(%s)' % value)

    def rotate(self, rotate_deg):
        camera = self.GetCamera()
        camera.Roll(-rotate_deg)
        camera.Modified()
        self.vtk_interactor.Render()
        self.log_command('rotate(%s)' % rotate_deg)

    def on_rotate_clockwise(self):
        self.rotate(15.0)

    def on_rotate_cclockwise(self):
        self.rotate(-15.0)

    def on_increase_magnification(self):
        self.zoom(1.1)

    def on_decrease_magnification(self):
        self.zoom(1.0/1.1)

    def on_flip_edges(self):
        self.is_edges = not self.is_edges
        self.edgeActor.SetVisibility(self.is_edges)
        #self.edgeActor.GetProperty().SetColor(0, 0, 0)  # cart3d edge color isn't black...
        self.edgeActor.Modified()
        #self.widget.Update()
        self._update_camera()
        #self.refresh()
        self.log_command('on_flip_edges()')

    def get_edges(self):
        """
        .. todo:: For some reason, the edge color is set to the parent
        surface's color instead of black
        """
        edges = vtk.vtkExtractEdges()
        if self.vtk_version[0] >= 6:
            edges.SetInputData(self.grid)
            self.edgeMapper.SetInputData(edges.GetOutput())
        else:
            edges.SetInput(self.grid)
            self.edgeMapper.SetInput(edges.GetOutput())

        self.edgeActor.SetMapper(self.edgeMapper)
        self.edgeActor.GetProperty().SetColor(0, 0, 0)

        prop = self.edgeActor.GetProperty()
        prop.SetColor(0, 0, 0)
        self.edgeActor.SetVisibility(self.is_edges)
        self.rend.AddActor(self.edgeActor)

    def createText(self, position, label, textSize=18, movable=False):
        txt = vtk.vtkTextActor()
        txt.SetInput(label)
        txtprop = txt.GetTextProperty()
        #txtprop.SetFontFamilyToArial()
        txtprop.SetFontSize(textSize)
        txtprop.SetColor(1, 1, 1)
        txt.SetDisplayPosition(*position)

        #print("dir(text) = ",dir(txt))
        txt.VisibilityOff()

        #txt.SetDisplayPosition(5,5) # bottom left
        #txt.SetDisplayPosition(5,95)
        #txt.SetPosition(0.1,0.5)

        # assign actor to the renderer
        self.rend.AddActor(txt)
        self.textActors[self.iText] = txt
        self.iText += 1

    def TurnTextOff(self):
        for (i, text) in iteritems(self.textActors):
            text.VisibilityOff()

    def TurnTextOn(self):
        for (i, text) in iteritems(self.textActors):
            text.VisibilityOn()

    def build_lookup_table(self):
        """TODO: add support for NanColors"""
        self.colorFunction = vtk.vtkColorTransferFunction()
        self.colorFunction.SetColorSpaceToHSV()
        self.colorFunction.HSVWrapOff()

        drange = [10., 20.]
        # blue - low
        # red - high
        self.colorFunction.AddRGBPoint(drange[0], 0.0, 0.0, 1.0)
        self.colorFunction.AddRGBPoint(drange[1], 1.0, 0.0, 0.0)

        self.scalarBar.SetTitle("Title1")
        self.scalarBar.SetLookupTable(self.colorFunction)
        self.scalarBar.SetOrientationToVertical()
        #print(dir(self.scalarBar))
        #print(dir(self.colorFunction))
        #self.scalarBar.SetNanColor(0., 0., 0.) # RGB color - black
        #self.scalarBar.SetNanColor(1., 1., 1., 0.) # RGBA color - white

        self.scalarBar.SetHeight(0.9)
        self.scalarBar.SetWidth(0.20)  # the width is set first
        # after the width is set, this is adjusted
        self.scalarBar.SetPosition(0.77, 0.1)

        propTitle = vtk.vtkTextProperty()
        propTitle.SetFontFamilyToArial()
        #propTitle.ItalicOff()
        propTitle.BoldOn()
        propTitle.ShadowOn()

        propLabel = vtk.vtkTextProperty()
        propLabel.BoldOff()
        propLabel.ShadowOn()

        scalar_range = self.grid.GetScalarRange()
        self.aQuadMapper.SetScalarRange(scalar_range)
        self.aQuadMapper.SetLookupTable(self.colorFunction)

        #self.scalarBar.SetTitleTextProperty(propTitle)
        #self.scalarBar.SetLabelTextProperty(propLabel)
        self.scalarBar.SetLabelFormat("%i")

        # allows 0-1 to be nice number when ranging values (gotta pick something)
        self.scalarBar.SetNumberOfLabels(11)
        self.scalarBar.SetMaximumNumberOfColors(11)

        #self.scalarBar.VisibilityOff()  # first load -> scalar bar off
        #self.scalarBar.ShadowOn()
        #self.scalarBar.RepositionableOn()
        self.rend.AddActor(self.scalarBar)
        self.scalarBar.VisibilityOff()
        #return scalarBar

    def _create_load_file_dialog(self, qt_wildcard, Title):
        # getOpenFileName return QString and we want Python string
        fname, wildcard_level = QtGui.QFileDialog.getOpenFileNameAndFilter(self, Title, self.last_dir, qt_wildcard)
        return str(wildcard_level), str(fname)

    def start_logging(self):
        if self.html_logging:
            log = SimpleLogger('debug', lambda x, y: self.logg_msg(x, y))
            # logging needs synchronizing, so the messages from different threads
            # would not be interleave
            self.log_mutex = QtCore.QReadWriteLock()
        else:
            log = SimpleLogger('debug', lambda x, y: print(x, y))
        self.log = log

    def build_fmts(self, fmt_order, stop_on_failure=False):
        fmts = []
        for fmt in fmt_order:
            if hasattr(self, 'get_%s_wildcard_geometry_results_functions' % fmt):
                func = 'get_%s_wildcard_geometry_results_functions' % fmt
                data = getattr(self, func)()
                msg = 'macro_name, geo_fmt, geo_func, res_fmt, res_func = data\n'
                msg += 'data = %s'
                assert len(data) == 5, msg % str(data)
                macro_name, geo_fmt, geo_func, res_fmt, res_func = data
                fmts.append((fmt, macro_name, geo_fmt, geo_func, res_fmt, res_func))
            else:
                if stop_on_failure:
                    func = 'get_%s_wildcard_geometry_results_functions does not exist' % fmt
                    raise RuntimeError(func)

        if len(fmts) == 0:
            RuntimeError('No formats...expected=%s' % fmt_order)
        self.fmts = fmts

        #fmts = [
            ## results
            #('nastran', 'Nastran', 'Nastran BDF (*.bdf; *.dat; *.nas)', self.load_nastran_geometry, 'Nastran OP2 (*.op2)', self.load_nastran_results),
            #('cart3d', 'Cart3d', 'Cart3d (*.tri; *.triq)', self.load_cart3d_geometry, 'Cart3d (*.triq)', self.load_cart3d_results),
            #('panair', 'Panair', 'Panair (*.inp)', self.load_panair_geometry, 'Panair (*.agps);;Panair (*.out)',  self.load_panair_results),
            #('shabp', 'S/HABP', 'Shabp (*.geo; *.mk5; *.inp)', self.load_shabp_geometry, 'Shabp (*.out)', self.load_shabp_results),
            #('usm3d', 'Usm3D', 'USM3D (*.cogsg; *.front)', self.load_usm3d_geometry, 'Usm3d (*.flo)', self.load_usm3d_results),

            ## no results
            #('lawgs', 'LaWGS', 'LaWGS (*.inp; *.wgs)', self.load_lawgs_geometry, None, None),
            #('tetgen', 'Tetgen', 'Tetgen (*.smesh)', self.load_tetgen_geometry, None, None),
            #('stl', 'STL', 'STereoLithography (*.STL)', self.load_stl_geometry, None, None),
            ##('plot3d', 'Plot3D', 'Plot3D (*.p3d; *.p3da)', self.load_plot3d_geometry, None, None),
        #]
        #self.fmts = fmts
        self.supported_formats = [fmt[0] for fmt in fmts]
        print('supported_formats = %s' % self.supported_formats)
        if len(fmts) == 0:
            raise RuntimeError('no modules were loaded...')

    def on_load_geometry(self, infile_name=None, geometry_format=None, plot=True):
        """
        Loads a baseline geometry

        :param infile_name: path to the filename (default=None -> popup)
        :param geometry_format: the geometry format for programmatic loading
        :param plot: Should the baseline geometry have results created and plotted/rendered?
                     If you're calling the on_load_results method immediately after, set it to False
                     (default=True)
        """
        wildcard = ''
        is_failed = False

        if geometry_format and geometry_format.lower() not in self.supported_formats:
            is_failed = True
            #if geometry_format in self.formats:
            msg = 'The import for the %r module failed.\n' % geometry_format
            #else:
                #msg = '%r is not a enabled format; enabled_formats=%s\n' % (geometry_format, self.supported_formats)
                #msg += str("formats = %s" % str(self.formats))
            self.log_error(msg)
            return is_failed

        if infile_name:
            geometry_format = geometry_format.lower()
            print("geometry_format = %r" % geometry_format)

            for fmt in self.fmts:
                fmt_name, _major_name, _geowild, _geofunc, _reswild, _resfunc = fmt
                if geometry_format == fmt_name:
                    load_function = _geofunc
                    if _reswild is None:
                        has_results = False
                    else:
                        has_results = True
                    break
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

            for fmt in self.fmts:
                fmt_name, _major_name, _geowild, _geofunc, _reswild, _resfunc = fmt
                formats.append(_major_name)
                wildcard_list.append(_geowild)
                load_functions.append(_geofunc)

                if _reswild is None:
                    has_results_list.append(False)
                else:
                    has_results_list.append(True)

            wildcard = ';;'.join(wildcard_list)

            # get the filter index and filename
            if infile_name is not None and geometry_format is not None:
                filter_index = formats.index(geometry_format)
            else:
                Title = 'Choose a Geometry File to Load'
                wildcard_index, infile_name = self._create_load_file_dialog(wildcard, Title)
                #print("infile_name = %r" % infile_name)
                #print("wildcard_index = %r" % wildcard_index)
                if not infile_name:
                    is_failed = True
                    return is_failed # user clicked cancel
                filter_index = wildcard_list.index(wildcard_index)

            geometry_format = formats[filter_index]
            load_function = load_functions[filter_index]
            has_results = has_results_list[filter_index]
            #return is_failed

        if load_function is not None:
            self.last_dir = os.path.split(infile_name)[0]

            self.grid.Reset()
            self.grid.Modified()
            self.grid2.Reset()
            self.grid2.Modified()
            #gridResult.Reset()
            #gridResult.Modified()

            if not os.path.exists(infile_name) and geometry_format:
                msg = 'input file=%r does not exist' % infile_name
                self.log_error(msg)
                self.log_error(print_bad_path(infile_name))
                return

            # clear out old data
            if self.modelType is not None:
                clear_name = 'clear_' + self.modelType
                try:
                    dy_method = getattr(self, clear_name)  # 'self.clear_nastran()'
                    dy_method()
                except:
                    print("method %r does not exist" % clear_name)
            self.log_info("reading %s file %r" % (geometry_format, infile_name))

            # inspect the load_geometry method to see what version it's using
            args, varargs, keywords, defaults = inspect.getargspec(load_function)
            try:
                if args[-1] == 'plot':
                    has_results = load_function(infile_name, self.last_dir, plot=plot)
                else:
                    name = load_function.__name__
                    self.log_error(str(args))
                    self.log_error("'plot' needs to be added to %r; args[-1]=%r" % (name, args[-1]))
                    has_results = load_function(infile_name, self.last_dir)
            except Exception as e:
                msg = traceback.format_exc()
                self.log_error(msg)
                #return
                raise
                #return
            #self.vtk_panel.Update()
            self.rend.ResetCamera()

        # the model has been loaded, so we enable load_results
        if filter_index >= 0:
            self.format = formats[filter_index].lower()
            if has_results:
                enable = True
            else:
                enable = False
            #self.load_results.Enable(enable)
        else: # no file specified
            return
        #print("on_load_geometry(infile_name=%r, geometry_format=None)" % infile_name)
        self.infile_name = infile_name

        self.out_filename = None
        #if self.out_filename is not None:
            #msg = '%s - %s - %s' % (self.format, self.infile_name, self.out_filename)
        #else:
        msg = '%s - %s' % (self.format, self.infile_name)
        self.set_window_title(msg)
        self.log_command("on_load_geometry(infile_name=%r, geometry_format=%r)" % (infile_name, self.format))

    def on_load_results(self, out_filename=None):
        """
        Loads a results file.  Must have called on_load_geometry first.

        :param out_filename: the path to the results file
        """
        geometry_format = self.format
        if self.format is None:
            msg = 'on_load_results failed:  You need to load a file first...'
            self.log_error(msg)
            raise RuntimeError(msg)

        if out_filename in [None, False]:
            Title = 'Select a Results File for %s' % self.format
            wildcard = None
            load_function = None

            for fmt in self.fmts:
                fmt_name, _major_name, _geowild, _geofunc, _reswild, _resfunc = fmt
                if geometry_format == fmt_name:
                    wildcard = _reswild
                    load_function = _resfunc
                    break
            else:
                msg = 'format=%r is not supported' % geometry_format
                self.log_error(msg)
                raise RuntimeError(msg)

            if wildcard is None:
                msg = 'format=%r has no method to load results' % geometry_format
                self.log_error(msg)
                return
            wildcard_index, out_filename = self._create_load_file_dialog(wildcard, Title)
        else:

            for fmt in self.fmts:
                fmt_name, _major_name, _geowild, _geofunc, _reswild, _resfunc = fmt
                print('fmt_name=%r geometry_format=%r' % (fmt_name, geometry_format))
                if fmt_name == geometry_format:
                    load_function = _resfunc
                    break
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
            #raise IOError(msg)
        self.last_dir = os.path.split(out_filename)[0]
        try:
            load_function(out_filename, self.last_dir)
        except Exception as e:
            msg = traceback.format_exc()
            self.log_error(msg)
            #return
            raise

        self.out_filename = out_filename
        msg = '%s - %s - %s' % (self.format, self.infile_name, out_filename)
        self.set_window_title(msg)
        print("on_load_results(%r)" % out_filename)
        self.out_filename = out_filename
        self.log_command("on_load_results(%r)" % out_filename)

    def setup_gui(self):
        assert self.fmts != [], 'supported_formats=%s' % self.supported_formats
        self.start_logging()
        settings = QtCore.QSettings()
        self.create_vtk_actors()

        # build GUI and restore saved application state
        self.restoreGeometry(settings.value("mainWindowGeometry").toByteArray())
        self.background_col = settings.value("backgroundColor", (0.1, 0.2, 0.4)).toPyObject()

        self.init_ui()
        self.init_cell_picker()

        self.restoreState(settings.value("mainWindowState").toByteArray())

        #-------------
        # loading
        self.show()

    def setup_post(self, inputs):
        self.load_batch_inputs(inputs)

        #-------------
        shots = inputs['shots']
        geometry_format = inputs['format']  # the active format loaded into the gui
        fname_input = inputs['input']
        fname_output = inputs['output']
        script = inputs['script']
        #-------------

        if shots is None:
            shots = []

        if shots:
        #for shot in shots:
            self.on_take_screenshot(shots)
            sys.exit('took screenshot %r' % shots)

        if script:
            self.on_run_script(script)

    def take_screenshot(self):
        """ Take a screenshot of a current view and save as a file"""
        self.on_take_screenshot(None)

    def on_take_screenshot(self, fname):
        """ Take a screenshot of a current view and save as a file"""
        if fname is None:
            filt = QtCore.QString()
            default_filename = ''

            Title = ''
            if self.Title is not None:
                Title = self.Title

            if self.out_filename is None:
                default_filename = ''
                if self.infile_name is not None:
                    base, ext = os.path.splitext(os.path.basename(self.infile_name))
                    default_filename = self.infile_name
            else:
                base, ext = os.path.splitext(os.path.basename(self.out_filename))
                default_filename = Title + '_' + base

            fname = str(QtGui.QFileDialog.getSaveFileName(self, (
                'Choose a filename and type'), default_filename, (
                    'PNG Image *.png (*.png);; '
                    'JPEG Image *.jpg *.jpeg (*.jpg, *.jpeg);; '
                    'TIFF Image *.tif *.tiff (*.tif, *.tiff);; '
                    'BMP Image *.bmp (*.bmp);; '
                    'PostScript Document *.ps (*.ps)'), filt))
            #print("fname=%r" % fname)
            if fname is None or fname == '':  # 2nd option
                return
            flt = str(filt).split()[0]
        else:
            base, ext = os.path.splitext(os.path.basename(fname))
            if ext.lower() in ['png', 'jpg', 'jpeg', 'tif', 'tiff', 'bmp', 'ps']:
                flt = ext.lower()
            else:
                flt = 'png'

        if fname:
            renderLarge = vtk.vtkRenderLargeImage()
            if self.vtk_version[0] >= 6:
                renderLarge.SetInput(self.rend)
            else:
                renderLarge.SetInput(self.rend)
            renderLarge.SetMagnification(self.magnify)

            nam, ext = os.path.splitext(fname)
            ext = ext.lower()
            for nam, exts, ob in (('PostScript', ['.ps'], vtk.vtkPostScriptWriter),
                                  ("BMP", ['.bmp'], vtk.vtkBMPWriter),
                                  ('JPG', ['.jpg', '.jpeg'], vtk.vtkJPEGWriter),
                                  ("TIFF", ['.tif', '.tiff'], vtk.vtkTIFFWriter)):
                if flt == nam:
                    fname = fname if ext in exts else fname + exts[0]
                    writer = ob()
                    break
            else:
                fname = fname if ext == '.png' else fname + '.png'
                writer = vtk.vtkPNGWriter()

            if self.vtk_version[0] >= 6:
                writer.SetInputData(renderLarge.GetOutputPort())
            else:
                writer.SetInputConnection(renderLarge.GetOutputPort())
            writer.SetFileName(fname)
            writer.Write()
            #self.log_info("Saved screenshot: " + fname)
            self.log_command('on_take_screenshot(%r)' % fname)

    def addGeometry(self):
        self.aQuadMapper = vtk.vtkDataSetMapper()
        if self.vtk_version[0] >= 6:
            self.aQuadMapper.SetInputData(self.grid)
        else:
            self.aQuadMapper.SetInput(self.grid)

        #self.warpVector = vtk.vtkWarpVector()
        #self.warpVector.SetInput(self.aQuadMapper.GetUnstructuredGridOutput())
        #aQuadMapper.SetInput(Filter.GetOutput())

        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(self.aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 0, 0)  # red
        self.rend.AddActor(geometryActor)

    def addAltGeometry(self):
        self.aQuadMapper = vtk.vtkDataSetMapper()
        if self.vtk_version[0] >= 6:
            self.aQuadMapper.SetInputData(self.grid2)
        else:
            self.aQuadMapper.SetInput(self.grid2)
        #aQuadMapper.SetInput(Filter.GetOutput())
        self.alt_geometry_actor = vtk.vtkActor()
        self.alt_geometry_actor.SetMapper(self.aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        self.alt_geometry_actor.GetProperty().SetDiffuseColor(1, 1, 0)  # green
        self.alt_geometry_actor.GetProperty().SetLineWidth(5)

        self.rend.AddActor(self.alt_geometry_actor)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()

    def on_update_scalar_bar(self, Title, min_value, max_value, data_format):
        self.Title = str(Title)
        self.min_value = float(min_value)
        self.max_value = float(max_value)
        try:
            data_format % 1
        except:
            self.log_error("failed applying the data formatter format=%r and should be of the form: '%i', '%8f', '%.2f', '%e', etc.")
            return
        self.data_format = data_format
        self.log_command('on_update_scalar_bar(%r, %r, %r')

    def ResetCamera(self):
        self.GetCamera().ResetCamera()

    def GetCamera(self):
        return self.rend.GetActiveCamera()

    def update_camera(self, code):
        camera = self.GetCamera()
        print("code =", code)
        if code == '+x':  # set x-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(1., 0., 0.)
        elif code == '-x':  # set x-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., -1.)
            camera.SetPosition(-1., 0., 0.)

        elif code == '+y':  # set y-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(0., 1., 0.)
        elif code == '-y':  # set y-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., -1.)
            camera.SetPosition(0., -1., 0.)

        elif code == '+z':  # set z-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 1., 0.)
            camera.SetPosition(0., 0., 1.)
        elif code == '-z':  # set z-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., -1., 0.)
            camera.SetPosition(0., 0., -1.)
        else:
            self.log_error('invalid camera code...%r' % code)
            return
        self._update_camera(camera)
        #print(dir(camera))
        #print(dir(self.vtk_interactor))
        #print(dir(self.rend))
        self.rend.ResetCamera()
        #self.vtk_interactor.ResetCamera()
        self.log_command('update_camera(%r)' % code)

    def _simulate_key_press(self, key):
        """
        A little hack method that simulates pressing the key for the VTK
        interactor. There is no easy way to instruct VTK to e.g. change mouse
        style to 'trackball' (as by pressing 't' key),
        (see http://public.kitware.com/pipermail/vtkusers/2011-November/119996.html)
        therefore we trick VTK to think that a key has been pressed.

        :param key: a key that VTK should be informed about, e.g. 't'
        """
        print("key = ", key)
        if key == 'f':  # change focal point
            return
        self.vtk_interactor._Iren.SetEventInformation(0, 0, 0, 0, key, 0, None)
        self.vtk_interactor._Iren.KeyPressEvent()
        self.vtk_interactor._Iren.CharEvent()

        #if key in ['y', 'z', 'X', 'Y', 'Z']:
            #self.update_camera(key)

    def _finish_results_io2(self, form, cases):
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        #print("ncases =", len(cases))
        #print("caseKeys =", self.caseKeys)

        if len(self.caseKeys) > 1:
            #print("finish_io case A")
            self.iCase = -1
            self.nCases = len(self.resultCases)  # number of keys in dictionary
        elif len(self.caseKeys) == 1:
            #print("finish_io case B")
            self.iCase = -1
            self.nCases = 1
        else:
            #print("finish_io case C")
            self.iCase = -1
            self.nCases = 0

        self.cycleResults_explicit()  # start at nCase=0
        if self.nCases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

        #data = [
        #    ('A',[]),
        #    ('B',[]),
        #    ('C',[]),
        #]

        #self.caseKeys= [(1, 'ElementID', 1, 'centroid', '%.0f'), (1, 'Region', 1, 'centroid', '%.0f')]
        data = []
        for key in self.caseKeys:
            t = (key[1], [])
            data.append(t)
        #data = self.caseKeys
        #print(data)
        self.res_widget.update_results(form)
        method = 'centroid' if self.is_centroidal else 'nodal'

        data2 = [(method, None, [])]
        self.res_widget.update_methods(data2)

    def finish_io(self, cases):
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        print("caseKeys = ", self.caseKeys)

        if len(self.resultCases) == 0:
            self.nCases = 1
            self.iCase = 0
        elif len(self.resultCases) == 1:
            self.nCases = 1
            self.iCase = 0
        else:
            self.nCases = len(self.resultCases) - 1  # number of keys in dictionary
            self.iCase = -1
        self.cycleResults()  # start at nCase=0

        if self.nCases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

    def _finish_results_io(self, cases):
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        #print("ncases =", len(cases))
        #print("caseKeys =", self.caseKeys)

        if len(self.caseKeys) > 1:
            #print("finish_io case A")
            self.iCase = -1
            self.nCases = len(self.resultCases)  # number of keys in dictionary
        elif len(self.caseKeys) == 1:
            #print("finish_io case B")
            self.iCase = -1
            self.nCases = 1
        else:
            #print("finish_io case C")
            self.iCase = -1
            self.nCases = 0

        self.cycleResults_explicit()  # start at nCase=0
        if self.nCases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

        #data = [
        #    ('A',[]),
        #    ('B',[]),
        #    ('C',[]),
        #]

        #self.caseKeys= [(1, 'ElementID', 1, 'centroid', '%.0f'), (1, 'Region', 1, 'centroid', '%.0f')]
        data = []
        for i, key in enumerate(self.caseKeys):
            t = (key[1], i, [])
            data.append(t)
            i += 1
        #data = self.caseKeys
        #print(data)
        self.res_widget.update_results(data)
        method = 'centroid' if self.is_centroidal else 'nodal'

        data2 = [(method, None, [])]
        self.res_widget.update_methods(data2)
