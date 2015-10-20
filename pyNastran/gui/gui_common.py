# -*- coding: utf-8 -*-
# pylint: disable=W0201,C0111
from __future__ import division, unicode_literals, print_function
from six import string_types, iteritems, itervalues
from six.moves import range

# standard library
import sys
import os.path
import datetime
import cgi #  html lib
import inspect
import traceback
from copy import deepcopy
from collections import OrderedDict

from PyQt4 import QtCore, QtGui
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from numpy import eye, array, zeros, loadtxt, int32
from numpy.linalg import norm

import pyNastran
from pyNastran.bdf.cards.baseCard import deprecated
from pyNastran.utils.log import SimpleLogger
from pyNastran.utils import print_bad_path

from pyNastran.gui.qt_files.gui_qt_common import GuiCommon
from pyNastran.gui.qt_files.scalar_bar import ScalarBar
from pyNastran.gui.qt_files.alt_geometry_storage import AltGeometry

from pyNastran.gui.menus.results_sidebar import Sidebar
from pyNastran.gui.menus.qt_legend import LegendPropertiesWindow
from pyNastran.gui.menus.clipping import ClippingPropertiesWindow
from pyNastran.gui.menus.camera import CameraWindow
from pyNastran.gui.menus.application_log import ApplicationLogDockWidget
from pyNastran.gui.menus.manage_actors import EditGroupProperties
from pyNastran.gui.menus.multidialog import MultiFileDialog
from pyNastran.gui.testing_methods import TestGuiCommon


class Interactor(vtk.vtkGenericRenderWindowInteractor):
    def __init__(self):
        #vtk.vtkGenericRenderWindowInteractor()
        pass

    def HighlightProp(self):
        print('highlight')


class PyNastranRenderWindowInteractor(QVTKRenderWindowInteractor):
    def __init__(self, parent=None):

        render_window = vtk.vtkRenderWindow()
        iren = Interactor()
        iren.SetRenderWindow(render_window)
        kwargs = {
            'iren' : iren,
            'rw' : render_window,
        }
        QVTKRenderWindowInteractor.__init__(self, parent=parent,
                                            iren=iren, rw=render_window)
        #self.Highlight

def loadtxt_nice(filename, delimiter=','):
    data = []
    delim = '\n\r \t' + delimiter
    with open(filename, 'r') as file_obj:
        line = file_obj.readline().strip(delim).split(delimiter)
        data.append(line)
    return array(data)

class GuiCommon2(QtGui.QMainWindow, GuiCommon, TestGuiCommon):
    def __init__(self, fmt_order, html_logging, inputs):
        QtGui.QMainWindow.__init__(self)
        GuiCommon.__init__(self)
        TestGuiCommon.__init__(self, res_widget=None)

        self.html_logging = html_logging
        self.is_testing = False
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
        # window variables
        self._legend_shown = False
        self._clipping_shown = False
        self._edit_group_properties_shown = False
        #-------------
        # inputs dict
        self.is_edges = inputs['is_edges']
        self.is_edges_black = self.is_edges
        # self.is_nodal = inputs['is_nodal']
        # self.is_centroidal = inputs['is_centroidal']
        self.magnify = inputs['magnify']

        #self.format = ''
        debug = inputs['debug']
        self.debug = debug
        assert debug in [True, False], 'debug=%s' % debug

        #-------------
        # file
        self.menu_bar_format = None
        self.format = None
        self.fmts = fmt_order
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
        self.actions = {}

        # initializes tools/checkables
        self.set_tools()

        # actor_slots
        self.text_actors = {}
        self.geometry_actors = {}
        self.alt_grids = {} #additional grids

        #geom = Geom(color, line_thickness, etc.)
        #self.geometry_properties = {
        #    'name' : Geom(),
        #}
        self.geometry_properties = {}

        self.magnify = 1
        self.iText = 0

        self.pick_state = 'node/centroid' # if self.is_centroidal else 'nodal'
        self.label_actors = {}
        self.label_ids = {}
        self.cameras = {}
        self.label_scale = 1.0 # in percent

        self.is_horizontal_scalar_bar = False
        self.scalar_bar = ScalarBar(self.is_horizontal_scalar_bar)

        self.result_cases = {}

        self.num_user_points = 0

    #def dragEnterEvent(self, e):
        #print(e)
        #print('drag event')
        #if e.mimeData().hasFormat('text/plain'):
            #e.accept()
        #else:
            #e.ignore()

    #def dropEvent(self, e):
        #print(e)
        #print('drop event')

    @property
    def scalarBar(self):
        return self.scalar_bar.scalar_bar

    @property
    def colorFunction(self):
        return self.scalar_bar.color_function

    #def get_color_function(self):
        #return self.scalar_bar.color_function

    @property
    def resultCases(self):
        return self.result_cases

    @resultCases.setter
    def resultCases(self, value):
        assert isinstance(value, dict), type(value)
        self.result_cases = value

    def set_window_title(self, msg):
        #msg2 = "%s - "  % self.base_window_title
        #msg2 += msg
        self.setWindowTitle(msg)

    def set_logo(self, logo):
        self._logo = logo

    def init_ui(self):
        """
        Initialize user iterface

        +--------------+
        | Window Title |
        +--------------+----------------+
        |  Menubar                      |
        +-------------------------------+
        |  Toolbar                      |
        +---------------------+---------+
        |                     |         |
        |                     |         |
        |                     | Results |
        |       VTK Frame     |  Dock   |
        |                     |         |
        |                     |         |
        +---------------------+---------+
        |                               |
        |      HTML Logging Dock        |
        |                               |
        +-------------------------------+
        """
        self.resize(1000, 700)
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
        #self.res_widget.update_results(data)
        #self.res_widget.setWidget(sidebar)

        self.res_dock.setWidget(self.res_widget)

        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.res_dock)
        #=========== Logging widget ===================
        if self.html_logging:
            self.log_dock = ApplicationLogDockWidget(self, execute_python=True)
            #self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock)
        #===============================================

        self._create_vtk_objects()
        self._build_menubar()

        # right sidebar
        self.res_dock.hide()
        self.build_vtk_frame()

        #compassRepresentation = vtk.vtkCompassRepresentation()
        #compassWidget = vtk.vtkCompassWidget()
        #compassWidget.SetInteractor(self.iren)
        #compassWidget.SetRepresentation(compassRepresentation)
        #compassWidget.EnabledOn()

    @property
    def dock_widget(self):
        return self.log_dock.dock_widget

    @property
    def log_widget(self):
        return self.log_dock.log_widget

    def _on_execute_python_button(self, clear=False):
        txt = str(self.log_dock.enter_data.toPlainText())
        self.log_command(txt)
        try:
            exec(txt)
        except TypeError:
            print(type(txt))
            raise
        except Exception as e:
            #traceback.print_stack()
            #traceback.print_exc(file=self.log_error)
            self.log_error(str(e))
            self.log_error(str(txt))
            return
        if clear:
            self.lock_dock.enter_data.clear()

    def load_batch_inputs(self, inputs):
        if not inputs['format']:
            return
        form = inputs['format'].lower()
        input_filename = inputs['input']
        results_filename = inputs['output']
        geom_script = inputs['geomscript']
        plot = True
        if results_filename:
            plot = False

        if geom_script is not None:
            self.on_run_script(geom_script)

        is_failed = self.on_load_geometry(input_filename, form, plot=plot)
        if is_failed:
            return
        if results_filename:
            self.on_load_results(results_filename)
        self.on_reset_camera()
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
                ('load_csv_nodal', 'Load CSV Nodal Results', '', None, 'Loads a custom nodal results file', self.on_load_nodal_results),  ## @todo no picture...
                ('load_csv_elemental', 'Load CSV Elemental Results', '', None, 'Loads a custom elemental results file', self.on_load_elemental_results),  ## @todo no picture...

                ('back_col', 'Change background color', 'tcolorpick.png', None, 'Choose a background color', self.change_background_color),
                ('label_col', 'Change label color', 'tcolorpick.png', None, 'Choose a label color', self.change_label_color),
                ('text_col', 'Change text color', 'tcolorpick.png', None, 'Choose a text color', self.change_text_color),

                ('label_clear', 'Clear current labels', '', None, 'Clear current labels', self.clear_labels),
                ('label_reset', 'Clear all labels', '', None, 'Clear all labels', self.reset_labels),

                ('legend', 'Modify legend', 'legend.png', None, 'Set Legend', self.set_legend),
                ('clipping', 'Set clipping', '', None, 'Set Clipping', self.set_clipping),
                ('axis', 'Show/Hide Axis', 'axis.png', None, 'Show/Hide Global Axis', self.on_show_hide_axes),

                ('wireframe', 'Wireframe Model', 'twireframe.png', 'w', 'Show Model as a Wireframe Model', self.on_wireframe),
                ('surface', 'Surface Model', 'tsolid.png', 's', 'Show Model as a Surface Model', self.on_surface),
                ('edges', 'Show/Hide Edges', 'tedges.png', 'e', 'Show/Hide Model Edges', self.on_flip_edges),
                ('edges_black', 'Color Edges', '', 'b', 'Set Edge Color to Color/Black', self.on_set_edge_visibility),
                ('geo_properties', 'Edit Geometry Properties', '', None, 'Change Model Color/Opacity/Line Width', self.set_actor_properties),

                ('show_info', 'Show INFO', 'show_info.png', None, 'Show "INFO" messages', self.on_show_info),
                ('show_debug', 'Show DEBUG', 'show_debug.png', None, 'Show "DEBUG" messages', self.on_show_debug),
                ('show_gui', 'Show GUI', 'show_gui.png', None, 'Show "GUI" messages', self.on_show_gui),
                ('show_command', 'Show COMMAND', 'show_command.png', None, 'Show "COMMAND" messages', self.on_show_command),

                ('magnify', 'Magnify', 'plus_zoom.png', 'M', 'Increase Magnfication', self.on_increase_magnification),
                ('shrink', 'Shrink', 'minus_zoom.png', 'm', 'Decrease Magnfication', self.on_decrease_magnification),

                #('flip_pick', 'Flip Pick', '', 'CTRL+K', 'Flips the pick state from centroidal to nodal', self.on_flip_picker),
                #('cell_pick', 'Cell Pick', '', 'c', 'Centroidal Picking', self.on_cell_picker),
                #('node_pick', 'Node Pick', '', 'n', 'Nodal Picking', self.on_node_picker),

                ('rotate_clockwise', 'Rotate Clockwise', 'tclock.png', 'o', 'Rotate Clockwise', self.on_rotate_clockwise),
                ('rotate_cclockwise', 'Rotate Counter-Clockwise', 'tcclock.png', 'O', 'Rotate Counter-Clockwise', self.on_rotate_cclockwise),

                ('scshot', 'Take a Screenshot', 'tcamera.png', 'CTRL+I', 'Take a Screenshot of current view', self.on_take_screenshot),
                ('about', 'About pyNastran GUI', 'tabout.png', 'CTRL+H', 'About pyCart3d GUI and help on shortcuts', self.about_dialog),
                ('view', 'Camera View', 'view.png', None, 'Load the camera menu', self.view_camera),
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
        if 'nastran' in self.fmts:
            tools += [
                ('caero', 'Show/Hide CAERO Panels', '', None, 'Show/Hide CAERO Panel Outlines', self.toggle_caero_panels),
                ('caero_sub', 'Toggle CAERO Subpanels', '', None, 'Show/Hide CAERO Subanel Outlines', self.toggle_caero_sub_panels),
                ('conm', 'Toggle CONMs', '', None, 'Show/Hide CONMs', self.toggle_conms),
            ]
        self.tools = tools
        self.checkables = checkables

    def deprecated(self, old_name, new_name, deprecated_version):
        deprecated(old_name, new_name, deprecated_version, levels=[-1])

    def add_tools(self, tools):
        self.deprecated('add_tools', 'removed...', '0.7')
        self.tools += tools

    def on_flip_picker(self):
        return
        # if self.pick_state == 'centroidal':
            # self.pick_state = 'nodal'

        # elif self.pick_state == 'nodal':
            # self.pick_state = 'centroidal'
        # else:
            # raise RuntimeError(self.pick_state)
        # self.log_command("on_flip_pick() # pick_state='%s'" % self.pick_state)

    def _create_menu_items(self, actions=None):
        if actions is None:
            actions = self.actions
        self.menu_file = self.menubar.addMenu('&File')
        self.menu_view = self.menubar.addMenu('&View')
        self.menu_window = self.menubar.addMenu('&Window')
        self.menu_help = self.menubar.addMenu('&Help')

        self.menu_hidden = self.menubar.addMenu('&Hidden')
        self.menu_hidden.setVisible(False)

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
                #tool = (script, script, 'python48.png', None, '',
                        #lambda: self.on_run_script(fname) )
                #tools.append(tool)
        #else:
        self.menu_scripts = None

        menu_window = ['toolbar', 'reswidget']
        menu_view = [
            'scshot', '', 'wireframe', 'surface', 'creset', '',
            'back_col', 'text_col', '',
            'label_col', 'label_clear', 'label_reset', '',
            'legend', 'geo_properties', '', 'clipping', 'axis', 'edges', 'edges_black',
        ]
        if self.html_logging:
            self.actions['logwidget'] = self.log_dock.dock_widget.toggleViewAction()
            self.actions['logwidget'].setStatusTip("Show/Hide application log")
            menu_view += ['', 'show_info', 'show_debug', 'show_gui', 'show_command']
            menu_window += ['logwidget']

        menu_items = [
            (self.menu_file, ('load_geometry', 'load_results', 'load_csv_nodal', 'load_csv_elemental', 'script', '', 'exit')),
            (self.menu_view, tuple(menu_view)),
            (self.menu_window, tuple(menu_window)),
            (self.menu_help, ('about',)),
            (self.menu_scripts, scripts),
            (self.toolbar, ('reload', 'load_geometry', 'load_results',
                            'x', 'y', 'z', 'X', 'Y', 'Z',
                            'magnify', 'shrink', 'rotate_clockwise', 'rotate_cclockwise',
                            'wireframe', 'surface', 'edges', 'creset', 'view', 'scshot', '', 'exit')),
            (self.menu_hidden, ('cycle_res',)),
            # (self.menu_scripts, ()),
            #(self._dummy_toolbar, ('cell_pick', 'node_pick'))
        ]
        return menu_items

    def _build_menubar(self):
        ## toolbar
        self.toolbar = self.addToolBar('Show toolbar')
        self.toolbar.setObjectName('main_toolbar')

        # the dummy toolbar stores actions but doesn't get shown
        # in other words, it can set shortcuts
        #self._dummy_toolbar = self.addToolBar('Dummy toolbar')
        #self._dummy_toolbar.setObjectName('dummy_toolbar')
        self.menubar = self.menuBar()

        actions = self._prepare_actions(self._icon_path, self.tools, self.checkables)
        menu_items = self._create_menu_items(actions)
        self._populate_menu(menu_items)

    def _populate_menu(self, menu_items):
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

                    action = self.actions[i] #if isinstance(i, string_types) else i()
                    menu.addAction(action)
        #self._create_plane_from_points(None)

    def _update_menu(self, menu_items):
        for menu, items in menu_items:
            menu.clear()
        self._populate_menu(menu_items)

    def _create_plane_from_points(self, points):
        origin, vx, vy, vz, x_limits, y_limits = self._fit_plane(points)

        # We create a 100 by 100 point plane to sample
        splane = vtk.vtkPlaneSource()
        plane = splane.GetOutput()

        dx = max(x_limits) - min(x_limits)
        dy = max(y_limits) - min(y_limits)
        #dx = 1.
        #dy = 3.

        # we need to offset the origin of the plane because the "origin"
        # is at the lower left corner of the plane and not the centroid
        offset = (dx * vx + dy * vy) / 2.
        origin -= offset
        splane.SetCenter(origin)

        splane.SetNormal(vz)

        # Point 1 defines the x-axis and the x-size
        # Point 2 defines the y-axis and the y-size
        splane.SetPoint1(origin + dx * vx)
        splane.SetPoint2(origin + dy * vy)

        actor = vtk.vtkLODActor()
        mapper = vtk.vtkPolyDataMapper()

        if self.vtk_version <= 5:
            mapper.SetInputData(plane)
        else:
            mapper.SetInput(plane)

        actor.GetProperty().SetColor(1., 0., 0.)
        actor.SetMapper(mapper)
        self.rend.AddActor(actor)
        splane.Update()

    def _fit_plane(self, points):
        origin = array([34.60272856552356, 16.92028913186242, 37.805958003209184])
        vx = array([1., 0., 0.])
        vy = array([0., 1., 0.])
        vz = array([0., 0., 1.])
        x_limits = [-1., 2.]
        y_limits = [0., 1.]
        return origin, vx, vy, vz, x_limits, y_limits

    def _prepare_actions(self, icon_path, tools, checkables=None):
        """
        Prepare actions that will  be used in application in a way
        that's independent of the  menus & toolbar
        """
        if checkables is None:
            checkables = []
        print('---------------------------')
        for tool in tools:
            (nam, txt, icon, shortcut, tip, func) = tool
            if nam in self.actions:
                self.log_error('trying to create a duplicate action %r' % nam)
                continue
            #print("name=%s txt=%s icon=%s shortcut=%s tip=%s func=%s"
                  #% (nam, txt, icon, shortcut, tip, func))
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
                self.actions[nam] = QtGui.QAction(ico, txt, self, checkable=True)
                self.actions[nam].setChecked(True)
            else:
                self.actions[nam] = QtGui.QAction(ico, txt, self)

            if shortcut:
                self.actions[nam].setShortcut(shortcut)
                #actions[nam].setShortcutContext(QtCore.Qt.WidgetShortcut)
            if tip:
                self.actions[nam].setStatusTip(tip)
            if func:
                self.actions[nam].triggered.connect(func)

        self.actions['toolbar'] = self.toolbar.toggleViewAction()
        self.actions['toolbar'].setStatusTip("Show/Hide application toolbar")

        self.actions['reswidget'] = self.res_dock.toggleViewAction()
        self.actions['reswidget'].setStatusTip("Show/Hide results selection")
        return self.actions

    def logg_msg(self, typ, msg):
        """
        Add message to log widget trying to choose right color for it.

        Parameters
        ----------
        msg : str
            message to be displayed
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
        cols = {
            "GUI" : "blue",
            "COMMAND" : "green",
            "GUI ERROR" : "Crimson",
            "DEBUG" : dark_orange,
            'WARNING' : "purple",
            # INFO - black
        }
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

    def change_background_color(self):
        """ Choose a background color """
        self._change_color('background', self.background_col, self.set_background_color)

    def change_label_color(self):
        """ Choose a label color """
        self._change_color('label', self.label_col, self.set_label_color)

    def change_text_color(self):
        """ Choose a text color """
        self._change_color('text', self.text_col, self.set_text_color)

    def _change_color(self, msg, rgb_color_floats, call_func):
        c = [int(255 * i) for i in rgb_color_floats]
        col = QtGui.QColorDialog.getColor(QtGui.QColor(*c), self, "Choose a %s color" % msg)
        if col.isValid():
            color = col.getRgbF()[:3]
            call_func(color)

    def set_background_color(self, color):
        """Set the background color"""
        self.background_col = color
        self.rend.SetBackground(*color)
        self.log_command('set_background_color(%s, %s, %s)' % color)

    def set_label_color(self, color):
        """Set the label color"""
        self.label_col = color
        for follower_actors in itervalues(self.label_actors):
            for follower_actor in follower_actors:
                prop = follower_actor.GetProperty()
                prop.SetColor(*color)
        self.log_command('set_label_color(%s, %s, %s)' % color)

    def set_text_color(self, color):
        """Set the text color"""
        self.text_col = color
        for text_actor in itervalues(self.text_actors):
            text_actor.GetTextProperty().SetColor(color)
        self.log_command('set_text_color(%s, %s, %s)' % color)

    def create_coordinate_system(self, label='', origin=None, matrix_3x3=None, Type='xyz'):
        """
        Creates a coordinate system

        Parameters
        ----------
        label : str
            the coord id or other unique label (default is empty to indicate the global frame)
        origin : (3, ) ndarray/list/tuple
            the origin
        matrix_3x3 : (3, 3) ndarray
            a standard Nastran-style coordinate system
        Type : str
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
        self.create_coordinate_system(label='', origin=None, matrix_3x3=None, Type='xyz')

    def create_corner_axis(self):
        axes = vtk.vtkAxesActor()
        self.corner_axis = vtk.vtkOrientationMarkerWidget()
        self.corner_axis.SetOrientationMarker(axes)
        self.corner_axis.SetInteractor(self.vtk_interactor)
        self.corner_axis.SetEnabled(1)
        self.corner_axis.InteractiveOff()

    def on_show_hide_axes(self):
        """
        show/hide axes
        """
        # this method should handle all the coords when
        # there are more then one
        if self._is_axes_shown:
            for axis in itervalues(self.axes):
                axis.VisibilityOff()
        else:
            for axis in itervalues(self.axes):
                axis.VisibilityOn()
        self._is_axes_shown = not self._is_axes_shown

    def create_vtk_actors(self):
        self.rend = vtk.vtkRenderer()

        # vtk actors
        self.grid = vtk.vtkUnstructuredGrid()
        self.grid2 = vtk.vtkUnstructuredGrid()
        #self.emptyResult = vtk.vtkFloatArray()
        #self.vectorResult = vtk.vtkFloatArray()

        # edges
        self.edgeActor = vtk.vtkLODActor()
        self.edgeMapper = vtk.vtkPolyDataMapper()

        self.create_cell_picker()

        # axes
        self.create_global_axes()

    def create_alternate_vtk_grid(self, name, color=None, line_width=5, opacity=1.0, point_size=1,
                                  bar_scale=0.0, representation=None, is_visible=True):
        self.alt_grids[name] = vtk.vtkUnstructuredGrid()
        self.geometry_properties[name] = AltGeometry(self, name, color=color,
                                                     line_width=line_width, opacity=opacity,
                                                     point_size=point_size, bar_scale=bar_scale,
                                                     representation=representation, is_visible=is_visible)

    def _create_vtk_objects(self):
        """creates some of the vtk objects"""
        #Frame that VTK will render on
        self.vtk_frame = QtGui.QFrame()

        #Qt VTK QVTKRenderWindowInteractor
        self.vtk_interactor = QVTKRenderWindowInteractor(parent=self.vtk_frame)
        #self.vtk_interactor = PyNastranRenderWindowInteractor(parent=self.vtk_frame)
        self.iren = self.vtk_interactor

    def build_vtk_frame(self):
        vtk_hbox = QtGui.QHBoxLayout()
        vtk_hbox.setContentsMargins(2, 2, 2, 2)

        vtk_hbox.addWidget(self.vtk_interactor)
        self.vtk_frame.setLayout(vtk_hbox)
        self.vtk_frame.setFrameStyle(QtGui.QFrame.NoFrame | QtGui.QFrame.Plain)
        # this is our main, 'central' widget
        self.setCentralWidget(self.vtk_frame)

        #=============================================================
        # +-----+-----+
        # |     |     |
        # |  A  |  B  |
        # |     |     |
        # +-----+-----+
        # xmin, xmax, ymin, ymax
        nframes = 1
        #nframes = 2
        if nframes == 2:
            # xmin, ymin, xmax, ymax
            frame1 = [0., 0., 0.5, 1.0]
            frame2 = [0.5, 0., 1., 1.0]
            #frames = [frame1, frame2]
            self.rend.SetViewport(*frame1)
        self.vtk_interactor.GetRenderWindow().AddRenderer(self.rend)

        if nframes == 2:
            rend = vtk.vtkRenderer()
            rend.SetViewport(*frame2)
            self.vtk_interactor.GetRenderWindow().AddRenderer(rend)

        self.vtk_interactor.GetRenderWindow().Render()
        #self.load_nastran_geometry(None, None)

        #for cid, axes in iteritems(self.axes):
            #self.rend.AddActor(axes)
        self.addGeometry()
        if nframes == 2:
            rend.AddActor(self.geom_actor)

        # initialize geometry_actors
        self.geometry_actors = {
            'main' : self.geom_actor,
        }
        #self.addAltGeometry()
        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(*self.background_col)

        self.rend.ResetCamera()
        self._simulate_key_press('t') # change mouse style to trackball
        self.build_lookup_table()

        text_size = 14 * self.magnify
        self.createText([5, 50], 'Max  ', text_size)  # text actor 0
        self.createText([5, 35], 'Min  ', text_size)  # text actor 1
        self.createText([5, 20], 'Word1', text_size)  # text actor 2
        self.createText([5, 5], 'Word2', text_size)  # text actor 3

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

    def on_run_script(self, python_file=False):
        print('python_file =', python_file)
        if python_file in [None, False]:
            Title = 'Choose a Python Script to Run'
            wildcard = "Python (*.py)"
            wildcard_index, infile_name = self._create_load_file_dialog(wildcard, Title)
            if not infile_name:
                is_failed = True
                return is_failed # user clicked cancel

            #python_file = os.path.join(script_path, infile_name)
            python_file = os.path.join(infile_name)
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
        self.vtk_interactor.Render()

    def on_surface(self):
        if self.is_wireframe:
            self.log_command('on_surface()')
            for name, actor in iteritems(self.geometry_actors):
                #if name != 'main':
                    #print('name: %s\nrep: %s' % (name, self.geometry_properties[name].representation ))
                if name == 'main' or self.geometry_properties[name].representation in ['main', 'toggle']:
                    prop = actor.GetProperty()
                    prop.SetRepresentationToSurface()
            self.is_wireframe = False
            self.vtk_interactor.Render()

    def on_wireframe(self):
        if not self.is_wireframe:
            self.log_command('on_wireframe()')
            for name, actor in iteritems(self.geometry_actors):
                #if name != 'main':
                    #print('name: %s\nrep: %s' % (name, self.geometry_properties[name].representation ))
                if name == 'main' or self.geometry_properties[name].representation in ['main', 'toggle']:
                    prop = actor.GetProperty()
                    prop.SetRepresentationToWireframe()
                #prop.SetRepresentationToPoints()
                #prop.GetPointSize()
                #prop.SetPointSize(5.0)
                #prop.ShadingOff()
            self.vtk_interactor.Render()
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
        self.zoom(1.0 / 1.1)

    def on_flip_edges(self):
        self.is_edges = not self.is_edges
        self.edgeActor.SetVisibility(self.is_edges)
        #self.edgeActor.GetProperty().SetColor(0, 0, 0)  # cart3d edge color isn't black...
        self.edgeActor.Modified()
        #self.widget.Update()
        self._update_camera()
        #self.refresh()
        self.log_command('on_flip_edges()')

    def on_set_edge_visibility(self):
        #self.edgeActor.SetVisibility(self.is_edges_black)
        self.is_edges_black = not self.is_edges_black
        if self.is_edges_black:
            prop = self.edgeActor.GetProperty()
            prop.EdgeVisibilityOn()
        else:
            prop = self.edgeActor.GetProperty()
            prop.EdgeVisibilityOff()
        self.edgeActor.Modified()
        prop.Modified()
        self.vtk_interactor.Render()

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
        self.edgeMapper.SetLookupTable(self.colorFunction)
        self.edgeMapper.SetResolveCoincidentTopologyToPolygonOffset()

        prop = self.edgeActor.GetProperty()
        prop.SetColor(0, 0, 0)
        self.edgeActor.SetVisibility(self.is_edges)
        self.rend.AddActor(self.edgeActor)

    def createText(self, position, label, text_size=18, movable=False):
        text_actor = vtk.vtkTextActor()
        text_actor.SetInput(label)
        text_prop = text_actor.GetTextProperty()
        #text_prop.SetFontFamilyToArial()
        text_prop.SetFontSize(int(text_size))
        text_prop.SetColor(self.text_col)
        text_actor.SetDisplayPosition(*position)

        text_actor.VisibilityOff()

        #txt.SetDisplayPosition(5, 5) # bottom left
        #txt.SetDisplayPosition(5, 95)
        #txt.SetPosition(0.1, 0.5)

        # assign actor to the renderer
        self.rend.AddActor(text_actor)
        self.text_actors[self.iText] = text_actor
        self.iText += 1

    def TurnTextOff(self):
        for text in itervalues(self.text_actors):
            text.VisibilityOff()

    def TurnTextOn(self):
        for text in itervalues(self.text_actors):
            text.VisibilityOn()

    def build_lookup_table(self):
        """TODO: add support for NanColors"""
        scalar_range = self.grid.GetScalarRange()
        #print('min = %s\nmax = %s' % scalar_range)
        self.aQuadMapper.SetScalarRange(scalar_range)
        self.aQuadMapper.SetLookupTable(self.colorFunction)
        self.rend.AddActor(self.scalarBar)

    def _create_load_file_dialog(self, qt_wildcard, Title):
        # getOpenFileName return QString and we want Python string
        fname, wildcard_level = QtGui.QFileDialog.getOpenFileNameAndFilter(
            self, Title, self.last_dir, qt_wildcard)
        return str(wildcard_level), str(fname)

    def _create_load_file_dialog2(self, qt_wildcard, Title):
        # getOpenFileName return QString and we want Python string
        #Title = 'Load a Tecplot Geometry/Results File'
        last_dir = ''
        #qt_wildcard = ['Tecplot Hex Binary (*.tec; *.dat)']
        dialog = MultiFileDialog()
        dialog.setWindowTitle(Title)
        dialog.setDirectory(self.last_dir)
        dialog.setFilters(qt_wildcard.split(';;'))
        if dialog.exec_() == QtGui.QDialog.Accepted:
            outfiles = dialog.selectedFiles()
            wildcard_level = dialog.selectedFilter()
            return str(wildcard_level), str(fname)
        return None, None

    def start_logging(self):
        if self.html_logging:
            log = SimpleLogger('debug', 'utf-8', lambda x, y: self.logg_msg(x, y))
            # logging needs synchronizing, so the messages from different
            # threads would not be interleave
            self.log_mutex = QtCore.QReadWriteLock()
        else:
            log = SimpleLogger('debug', 'utf-8', lambda x, y: print(x, y))
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

        self.supported_formats = [fmt[0] for fmt in fmts]
        print('supported_formats = %s' % self.supported_formats)
        if len(fmts) == 0:
            raise RuntimeError('no modules were loaded...')

    def on_load_geometry(self, infile_name=None, geometry_format=None, plot=True):
        """
        Loads a baseline geometry

        Parameters
        ----------
        infile_name : str; default=None -> popup
            path to the filename
        geometry_format : str; default=None
            the geometry format for programmatic loading
        plot : bool; default=True
            Should the baseline geometry have results created and plotted/rendered?
            If you're calling the on_load_results method immediately after, set it to False

        """
        wildcard = ''
        is_failed = False

        if geometry_format and geometry_format.lower() not in self.supported_formats:
            is_failed = True
            #if geometry_format in self.formats:
            msg = 'The import for the %r module failed.\n' % geometry_format
            #else:
            #msg += '%r is not a enabled format; enabled_formats=%s\n' % (geometry_format, self.supported_formats)
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
                #raise NotImplementedError('on_load_geometry; infile_name=%r format=%r' % (infile_name, geometry_format))
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
                title = 'Choose a Geometry File to Load'
                wildcard_index, infile_name = self._create_load_file_dialog(wildcard, title)
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

            # reset alt grids
            names = self.alt_grids.keys()
            for name in names:
                self.alt_grids[name].Reset()
                self.alt_grids[name].Modified()

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
        self.update_menu_bar()
        self.log_command("on_load_geometry(infile_name=%r, geometry_format=%r)" % (infile_name, self.format))

    def _update_menu_bar_to_format(self, fmt, method):
        self.menu_bar_format = fmt
        tools, menu_items = getattr(self, method)()
        actions = self._prepare_actions(self._icon_path, tools, self.checkables)
        self._update_menu(menu_items)

    def update_menu_bar(self):
        # the format we're switching to
        method_new = '_create_%s_tools_and_menu_items' % self.format
        method_cleanup = '_cleanup_%s_tools_and_menu_items' % self.menu_bar_format

        # the current state of the format
        #method_new = '_create_%s_tools_and_menu_items' % self.menu_bar_format
        self.menu_bar_format = 'cwo'
        if self.menu_bar_format is None:
            self._update_menu_bar_to_format(self.format, method_new)
        else:
            print('need to add %r' % method_new)
            if self.menu_bar_format != self.format:
                if hasattr(self, method_cleanup):
                #if hasattr(self, method_old):
                    self.menu_bar_format = None
                    getattr(self, method_cleanup)()

            if hasattr(self, method_new):
                self._update_menu_bar_to_format(self.format, method_new)
                    #self._update_menu_bar_to_format(self.format)
                    #actions = self._prepare_actions(self._icon_path, self.tools, self.checkables)
                    #menu_items = self._create_menu_items(actions)
                    #menu_items = self._create_menu_items()
                    #self._populate_menu(menu_items)

    def _on_load_nodal_elemental_results(self, Type, out_filename=None):
        """
        Loads a CSV/TXT results file.  Must have called on_load_geometry first.

        Parameters
        ----------
        out_filename : str / None
            the path to the results file
        """
        geometry_format = self.format
        if self.format is None:
            msg = 'on_load_results failed:  You need to load a file first...'
            self.log_error(msg)
            raise RuntimeError(msg)

        if out_filename in [None, False]:
            title = 'Select a %s Results File for %s' % (Type, self.format)
            wildcard = 'Delimited Text (*.txt; *.dat; *.csv)'
            wildcard_index, out_filename = self._create_load_file_dialog(wildcard, title)

        if out_filename == '':
            return
        if not os.path.exists(out_filename):
            msg = 'result file=%r does not exist' % out_filename
            self.log_error(msg)
            return
            #raise IOError(msg)
        # self.last_dir = os.path.split(out_filename)[0]
        try:
            self._load_csv(Type, out_filename)
        except Exception as e:
            msg = traceback.format_exc()
            self.log_error(msg)
            #return
            raise

        if 0:
            self.out_filename = out_filename
            msg = '%s - %s - %s' % (self.format, self.infile_name, out_filename)
            self.set_window_title(msg)
            self.out_filename = out_filename

        if Type == 'Nodal':
            self.log_command("_on_load_nodal_elemental_results(%r)" % out_filename)
        elif Type == 'Elemental':
            self.log_command("on_load_elemental_results(%r)" % out_filename)
        else:
            raise NotImplementedError(Type)

    def _load_csv(self, Type, out_filename):
        out_filename_short = os.path.basename(out_filename)

        ext = os.path.splitext(out_filename)[1].lower()
        if ext not in ['.csv', '.dat', '.txt']:
            raise NotImplementedError('extension=%r is not supported (use .dat, .txt, or .csv)' % ext)

        file_obj = open(out_filename, 'r')
        header_line = file_obj.readline().strip()
        if not header_line.startswith('#'):
            msg = 'Expected file of the form:\n'
            if ext in ['.dat', '.txt']:
                msg += '# var1 var2\n'
                msg += '1 2\n'
                msg += '3 4\n'
            elif ext in ['.csv']:
                # csv
                msg += '# var1, var2\n'
                msg += '1, 2\n'
                msg += '3, 4\n'
            else:
                raise NotImplementedError('extension=%r is not supported (use .dat, .txt, or .csv)' % ext)
            raise SyntaxError(msg)
        header_line = header_line.lstrip('# \t').strip()

        if ext in ['.dat', '.txt']:
            headers = header_line.split(' ')
            A = loadtxt(file_obj)
        elif ext in ['.csv']:
            headers = header_line.split(',')
            A = loadtxt(file_obj, delimiter=',')
        else:
            raise NotImplementedError('extension=%r is not supported (use .dat, .txt, or .csv)' % ext)
        file_obj.close()

        if len(A.shape) == 1:
            A = A.reshape(A.shape[0], 1)
        nrows, ncols = A.shape
        if ncols != len(headers):
            msg = 'Error loading csv/txt file\n'
            msg += 'ncols != len(headers); ncols=%s; len(headers)=%s\n'
            msg += 'headers = %s' % headers
            raise SyntaxError(msg)

        if Type == 'Nodal':
            assert nrows == self.nNodes, 'nrows=%s nnodes=%s' % (nrows, self.nNodes)
            Type2 = 'node'
        elif Type == 'Elemental':
            assert nrows == self.nElements, 'nrows=%s nelements=%s' % (nrows, self.nElements)
            Type2 = 'centroid'
        else:
            raise NotImplementedError(Type)

        formi = []
        form = self.get_form()
        icase = len(self.caseKeys)
        islot = self.caseKeys[0][0]
        for icol in range(ncols):
            datai = A[:, icol]
            header = headers[icol].strip()
            key = (islot, icase, header, 1, Type2, '%.3f')
            self.caseKeys.append(key)
            self.resultCases[key] = datai
            formi.append((header, icase, []))

            self.label_actors[header] = []
            self.label_ids[header] = set([])
            icase += 1
        form.append((out_filename_short, None, formi))

        self.nCases += ncols
        #cases[(ID, 2, 'Region', 1, 'centroid', '%i')] = regions
        self.res_widget.update_results(form)

    def on_load_nodal_results(self, out_filename=None):
        self._on_load_nodal_elemental_results('Nodal', out_filename)

    def on_load_elemental_results(self, out_filename=None):
        self._on_load_nodal_elemental_results('Elemental', out_filename)

    def on_load_results(self, out_filename=None):
        """
        Loads a results file.  Must have called on_load_geometry first.

        Parameters
        ----------
        out_filename : str / None
            the path to the results file
        """
        geometry_format = self.format
        if self.format is None:
            msg = 'on_load_results failed:  You need to load a file first...'
            self.log_error(msg)
            raise RuntimeError(msg)

        if out_filename in [None, False]:
            title = 'Select a Results File for %s' % self.format
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
            wildcard_index, out_filename = self._create_load_file_dialog(wildcard, title)
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
        nice_blue = (0.1, 0.2, 0.4)
        white = (1.0, 1.0, 1.0)
        #black = (0.0, 0.0, 0.0)
        red = (1.0, 0.0, 0.0)
        self.restoreGeometry(settings.value("mainWindowGeometry").toByteArray())
        self.background_col = settings.value("backgroundColor", nice_blue).toPyObject()
        self.label_col = settings.value("labelColor", red).toPyObject()
        self.text_col = settings.value("textColor", white).toPyObject()

        self.init_ui()
        self.init_cell_picker()

        self.restoreState(settings.value("mainWindowState").toByteArray())
        self.create_corner_axis()
        #-------------
        # loading
        self.show()

    def setup_post(self, inputs):
        self.load_batch_inputs(inputs)

        shots = inputs['shots']
        if shots is None:
            shots = []
        if shots:
        #for shot in shots:
            self.on_take_screenshot(shots)
            sys.exit('took screenshot %r' % shots)

        if inputs['user_points'] is not None:
            initial_colors = [(1.0, 0.145098039216, 1.0),
                              (0.0823529411765, 0.0823529411765, 1.0),
                              (0.0901960784314, 1.0, 0.941176470588),
                              (0.501960784314, 1.0, 0.0941176470588),
                              (1.0, 1.0, 0.117647058824),
                              (1.0, 0.662745098039, 0.113725490196)]
            for fname in inputs['user_points']:
                name = os.path.basename(fname).rsplit('.', 1)[0]
                color = initial_colors[self.num_user_points % len(initial_colors)]
                self.add_user_points(fname, name, color=color)
                self.num_user_points += 1


    def create_cell_picker(self):
        # cell picker
        self.cell_picker = vtk.vtkCellPicker()
        self.node_picker = vtk.vtkPointPicker()
        #self.cell_picker.SetTolerance(0.0005)

    def init_cell_picker(self):
        self.is_pick = False

        self.vtk_interactor.SetPicker(self.cell_picker)
        def annotate_cell_picker(object, event):
            #self.log_command("annotate_cell_picker()")
            picker = self.cell_picker
            if picker.GetCellId() < 0:
                #self.picker_textActor.VisibilityOff()
                pass
            else:
                world_position = picker.GetPickPosition()
                cell_id = picker.GetCellId()
                #ds = picker.GetDataSet()
                #select_point = picker.GetSelectionPoint()
                self.log_command("annotate_cell_picker()")
                self.log_info("XYZ Global = %s" % str(world_position))
                #self.log_info("cell_id = %s" % cell_id)
                #self.log_info("data_set = %s" % ds)
                #self.log_info("selPt = %s" % str(select_point))

                #method = 'get_result_by_cell_id()' # self.modelType
                #print('pick_state =', self.pick_state)

                duplicate_key = None

                icase = self.iCase
                key = self.caseKeys[icase]
                location = self.get_case_location(key)

                if location == 'centroid':
                    if self.pick_state == 'node/centroid':
                        duplicate_key = cell_id
                        result_name, result_value, xyz = self.get_result_by_cell_id(cell_id, world_position)
                        assert result_name in self.label_actors, result_name
                    else:
                        #cell = self.grid.GetCell(cell_id)
                        # get_nastran_centroidal_pick_state_nodal_by_xyz_cell_id()
                        method = 'get_centroidal_%s_result_pick_state_%s_by_xyz_cell_id' % (self.format, self.pick_state)
                        if hasattr(self, method):
                            methodi = getattr(self, method)
                            return_flag, value = methodi(world_position, cell_id)
                            if return_flag is True:
                                return
                        else:
                            msg = "pick_state is set to 'nodal', but the result is 'centroidal'\n"
                            msg += '  cannot find: self.%s(xyz, cell_id)' % method
                            self.log_error(msg)
                        return
                    self.log_info("%s = %s" % (result_name, result_value))
                elif location == 'node':
                    if self.pick_state == 'node/centroid':
                        result_name, result_value, node_id, xyz = self.get_result_by_xyz_cell_id(world_position, cell_id)
                        assert result_name in self.label_actors, result_name
                        assert not isinstance(xyz, int), xyz
                        duplicate_key = node_id
                    else:
                        method = 'get_nodal_%s_result_pick_state_%s_by_xyz_cell_id' % (self.format, self.pick_state)
                        if hasattr(self, method):
                            methodi = getattr(self, method)
                            return_flag, value = methodi(world_position, cell_id)
                            if return_flag is True:
                                return
                        else:
                            msg = "pick_state is set to 'centroidal', but the result is 'nodal'\n"
                            msg += '  cannot find: self.%s(xyz, cell_id)' % method
                            self.log_error(msg)
                        return
                    msg = "%s = %s" % (result_name, result_value)
                    if self.result_name in ['Node_ID', 'Node ID', 'NodeID']:
                        msg += '; xyz=(%s, %s, %s)' % tuple(xyz)
                    self.log_info(msg)
                else:
                    raise RuntimeError('invalid pick location=%r' % location)

                #print('key=%s exists=%s' % (duplicate_key, duplicate_key in self.label_ids[result_name]))
                if duplicate_key is not None and duplicate_key in self.label_ids[result_name]:
                    return
                self.label_ids[result_name].add(duplicate_key)

                #x, y, z = world_position
                x, y, z = xyz
                text = '(%.3g, %.3g, %.3g); %s' % (x, y, z, result_value)
                text = str(result_value)
                assert result_name in self.label_actors, result_name
                self._create_annotation(text, result_name, x, y, z)

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
        self.node_picker.AddObserver("EndPickEvent", annotate_point_picker)

        #self.cell_picker.AddObserver("EndPickEvent", on_cell_picker)
        #self.node_picker.AddObserver("EndPickEvent", on_node_picker)

    def _create_annotation(self, text, result_name, x, y, z):
        assert isinstance(result_name, string_types), 'result_name=%r type=%s' % (result_name, type(result_name))
        # http://nullege.com/codes/show/src%40p%40y%40pymatgen-2.9.6%40pymatgen%40vis%40structure_vtk.py/395/vtk.vtkVectorText/python
        if 1:
            source = vtk.vtkVectorText()
            source.SetText(text)

            # mappers are weird; they seem to do nothing
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())

            # the follower lets us set the position/size/color
            follower = vtk.vtkFollower()
            follower.SetMapper(mapper)
            follower.SetPosition((x, y, z))

            # 1 point = 1/72"
            # SetScale works on model scale size
            #follower.SetScale(0.5)
            follower.SetScale(self.dim_max * 0.01 * self.label_scale)

            prop = follower.GetProperty()
            prop.SetColor(self.label_col)
            #prop.SetOpacity( 0.3 );

            # we need to make sure the text rotates when the camera is changed
            camera = self.rend.GetActiveCamera()
            follower.SetCamera(camera)
        else:
            # Create a text mapper and actor to display the results of picking.
            textMapper = vtk.vtkTextMapper()
            textMapper.SetInput(text)

            tprop = textMapper.GetTextProperty()
            tprop.SetFontFamilyToArial()
            tprop.SetFontSize(10)
            tprop.BoldOn()
            tprop.ShadowOn()
            tprop.SetColor(self.label_col)

            textActor = vtk.vtkActor2D()
            textActor.GetPositionCoordinate().SetCoordinateSystemToWorld()
            textActor.SetPosition(world_position[:2])
            textActor.SetMapper(textMapper)
            follower = textActor

        # finish adding the actor
        self.rend.AddActor(follower)
        self.label_actors[result_name].append(follower)

        #self.picker_textMapper.SetInput("(%.6f, %.6f, %.6f)"% pickPos)
        #camera.GetPosition()
        #camera.GetClippingRange()
        #camera.GetFocalPoint()

    def _on_cell_picker(self, a):
        self.vtk_interactor.SetPicker(self.cell_picker)
        picker = self.cell_picker
        world_position = picker.GetPickPosition()
        cell_id = picker.GetCellId()
        select_point = picker.GetSelectionPoint()  # get x,y pixel coordinate

        self.log_info("world_position = %s" % str(world_position))
        self.log_info("cell_id = %s" % cell_id)
        self.log_info("select_point = %s" % str(select_point))

    def _on_node_picker(self, a):
        self.vtk_interactor.SetPicker(self.node_picker)
        picker = self.node_picker
        world_position = picker.GetPickPosition()
        node_id = picker.GetPointId()
        select_point = picker.GetSelectionPoint()  # get x,y pixel coordinate

        self.log_info("world_position = %s" % str(world_position))
        self.log_info("node_id = %s" % node_id)
        self.log_info("select_point = %s" % str(select_point))

    #def on_cell_picker(self):
        #self.log_command("on_cell_picker()")
        #picker = self.cell_picker
        #world_position = picker.GetPickPosition()
        #cell_id = picker.GetCellId()
        ##ds = picker.GetDataSet()
        #select_point = picker.GetSelectionPoint()  # get x,y pixel coordinate
        #self.log_info("world_position = %s" % str(world_position))
        #self.log_info("cell_id = %s" % cell_id)
        #self.log_info("select_point = %s" % str(select_point))
        #self.log_info("data_set = %s" % ds)

    def on_take_screenshot(self, fname=None, magnification=None):
        """ Take a screenshot of a current view and save as a file"""
        if fname is None or fname is False:
            filt = QtCore.QString()
            default_filename = ''

            title = ''
            if self.Title is not None:
                title = self.Title

            if self.out_filename is None:
                default_filename = ''
                if self.infile_name is not None:
                    base, ext = os.path.splitext(os.path.basename(self.infile_name))
                    default_filename = self.infile_name
                    default_filename = base + '.png'
            else:
                base, ext = os.path.splitext(os.path.basename(self.out_filename))
                default_filename = title + '_' + base + '.png'

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
            render_large = vtk.vtkRenderLargeImage()
            if self.vtk_version[0] >= 6:
                render_large.SetInput(self.rend)
            else:
                render_large.SetInput(self.rend)

            if magnification is None:
                magnify_min = 5
                magnify = self.magnify if self.magnify > magnify_min else magnify_min
            else:
                magnify = magnification

            self._update_text_size(magnify=magnify)
            render_large.SetMagnification(magnify)

            nam, ext = os.path.splitext(fname)
            ext = ext.lower()
            for nam, exts, obj in (('PostScript', ['.ps'], vtk.vtkPostScriptWriter),
                                   ("BMP", ['.bmp'], vtk.vtkBMPWriter),
                                   ('JPG', ['.jpg', '.jpeg'], vtk.vtkJPEGWriter),
                                   ("TIFF", ['.tif', '.tiff'], vtk.vtkTIFFWriter)):
                if flt == nam:
                    fname = fname if ext in exts else fname + exts[0]
                    writer = obj()
                    break
            else:
                fname = fname if ext == '.png' else fname + '.png'
                writer = vtk.vtkPNGWriter()

            if self.vtk_version[0] >= 6:
                writer.SetInputConnection(render_large.GetOutputPort())
            else:
                writer.SetInputConnection(render_large.GetOutputPort())
            writer.SetFileName(fname)
            writer.Write()
            #self.log_info("Saved screenshot: " + fname)
            self.log_command('on_take_screenshot(%r)' % fname)
            self._update_text_size(magnify=1.0)

    def _update_text_size(self, magnify=1.0):
        text_size = int(14 * magnify)
        for itext, text_actor in iteritems(self.text_actors):
            text_prop = text_actor.GetTextProperty()
            text_prop.SetFontSize(text_size)
        self.iText += 1

    def addGeometry(self):
        """
        #(N,)  for stress, x-disp
        #(N,3) for warp vectors/glyphs
        grid_result = vtk.vtkFloatArray()

        point_data = self.grid.GetPointData()
        cell_data = self.grid.GetCellData()

        self.grid.GetCellData().SetScalars(grid_result)
        self.grid.GetPointData().SetScalars(grid_result)


        self.aQuadMapper   <-input-> self.grid
        vtkDataSetMapper() <-input-> vtkUnstructuredGrid()

        self.aQuadMapper   <--map--> self.geom_actor <-add-> self.rend
        vtkDataSetMapper() <--map--> vtkActor()      <-add-> vtkRenderer()
        """
        self.aQuadMapper = vtk.vtkDataSetMapper()
        if self.vtk_version[0] >= 6:
            self.aQuadMapper.SetInputData(self.grid)
        else:
            self.aQuadMapper.SetInput(self.grid)

        if 0:
            self.warp_filter = vtk.vtkWarpVector()
            self.warp_filter.SetScaleFactor(50.0)
            self.warp_filter.SetInput(self.aQuadMapper.GetUnstructuredGridOutput())

            self.geom_filter = vtk.vtkGeometryFilter()
            self.geom_filter.SetInput(self.warp_filter.GetUnstructuredGridOutput())

            self.geom_mapper = vtk.vtkPolyDataMapper()
            self.geom_actor.setMapper(self.geom_mapper)

        if 0:
            #from vtk.numpy_interface import algorithms

            arrow = vtk.vtkArrowSource()
            arrow.PickableOff()

            self.glyph_transform = vtk.vtkTransform()
            self.glyph_transform_filter = vtk.vtkTransformPolyDataFilter()
            self.glyph_transform_filter.SetInputConnection(arrow.GetOutputPort())
            self.glyph_transform_filter.SetTransform(self.glyph_transform)

            self.glyph = vtk.vtkGlyph3D()
            #self.glyph.setInput(xxx)
            self.glyph.SetSource(self.glyph_transform_filter.GetOutput())

            self.glyph.SetVectorModeToUseVector()
            self.glyph.SetColorModeToColorByVector()
            self.glyph.SetScaleModeToScaleByVector()
            self.glyph.SetScaleFactor(1.0)

            self.append_filter = vtk.vtkAppendFilter()
            self.append_filter.AddInputConnection(self.grid.GetOutput())


        #self.warpVector = vtk.vtkWarpVector()
        #self.warpVector.SetInput(self.aQuadMapper.GetUnstructuredGridOutput())
        #aQuadMapper.SetInput(Filter.GetOutput())

        self.geom_actor = vtk.vtkLODActor()
        self.geom_actor.SetMapper(self.aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        #self.geom_actor.GetProperty().SetDiffuseColor(1, 0, 0)  # red
        self.rend.AddActor(self.geom_actor)

    def _add_alt_actors(self, grids_dict, names_to_ignore=None):
        if names_to_ignore is None:
            names_to_ignore = ['main']

        names = set(list(grids_dict.keys()))
        names_old = set(list(self.geometry_actors.keys()))
        names_old = names_old - set(names_to_ignore)
        #print('names_old1 =', names_old)

        #names_to_clear = names_old - names
        #self._remove_alt_actors(names_to_clear)
        #print('names_old2 =', names_old)
        #print('names =', names)
        for name in names:
            #print('adding %s' % name)
            grid = grids_dict[name]
            self._add_alt_geometry(grid, name)

    def _remove_alt_actors(self, names=None):
        if names is None:
            names = list(self.geometry_actors.keys())
            names.remove('main')
        for name in names:
            actor = self.geometry_actors[name]
            self.rend.RemoveActor(actor)
            del actor

    def _add_alt_geometry(self, grid, name, color=None, line_width=None, opacity=None, representation=None):
        """
        NOTE: color, line_width, opacity are ignored if name already exists
        """
        quadMapper = vtk.vtkDataSetMapper()
        if name in self.geometry_actors:
            alt_geometry_actor = self.geometry_actors[name]
            if self.vtk_version[0] >= 6:
                alt_geometry_actor.GetMapper().SetInputData(grid)
            else:
                alt_geometry_actor.GetMapper().SetInput(grid)
        else:
            if self.vtk_version[0] >= 6:
                quadMapper.SetInputData(grid)
            else:
                quadMapper.SetInput(grid)
            alt_geometry_actor = vtk.vtkActor()
            alt_geometry_actor.SetMapper(quadMapper)
            self.geometry_actors[name] = alt_geometry_actor

        #geometryActor.AddPosition(2, 0, 2)
        if name in self.geometry_properties:
            geom = self.geometry_properties[name]
        else:
            geom = AltGeometry(self, name, color=color, line_width=line_width,
                               opacity=opacity, representation=representation)
            self.geometry_properties[name] = geom

        color = geom.color_float
        opacity = geom.opacity
        point_size = geom.point_size
        representation = geom.representation
        line_width = geom.line_width
        #print('color_2014[%s] = %s' % (name, str(color)))
        assert isinstance(color[0], float), color
        assert color[0] <= 1.0, color

        prop = alt_geometry_actor.GetProperty()
        prop.SetDiffuseColor(color)
        prop.SetOpacity(opacity)
        if representation == 'point':
            prop.SetRepresentationToPoints()
            prop.SetPointSize(point_size)
        elif representation == 'surface':
            prop.SetRepresentationToSurface()
            prop.SetLineWidth(line_width)
        elif representation == 'wire':
            prop.SetRepresentationToWireframe()
            prop.SetLineWidth(line_width)

        self.rend.AddActor(alt_geometry_actor)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()

        if geom.is_visible:
            alt_geometry_actor.VisibilityOn()
        else:
            alt_geometry_actor.VisibilityOff()

        #print('current_actors = ', self.geometry_actors.keys())
        if hasattr(grid, 'Update'):
            grid.Update()
        alt_geometry_actor.Modified()

    def on_update_scalar_bar(self, title, min_value, max_value, data_format):
        self.Title = str(title)
        self.min_value = float(min_value)
        self.max_value = float(max_value)

        try:
            data_format % 1
        except:
            self.log_error("failed applying the data formatter format=%r and should be of the form: '%i', '%8f', '%.2f', '%e', etc.")
            return
        self.data_format = data_format
        self.log_command('on_update_scalar_bar(%r, %r, %r, %r)' % (title, min_value, max_value, data_format))

    def ResetCamera(self):
        self.GetCamera().ResetCamera()

    def GetCamera(self):
        return self.rend.GetActiveCamera()

    def update_camera(self, code):
        camera = self.GetCamera()
        #print("code =", code)
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
        self.rend.ResetCamera()
        self.log_command('update_camera(%r)' % code)

    def _simulate_key_press(self, key):
        """
        A little hack method that simulates pressing the key for the VTK
        interactor. There is no easy way to instruct VTK to e.g. change mouse
        style to 'trackball' (as by pressing 't' key),
        (see http://public.kitware.com/pipermail/vtkusers/2011-November/119996.html)
        therefore we trick VTK to think that a key has been pressed.

        Parameters
        ----------
        key : str
            a key that VTK should be informed about, e.g. 't'
        """
        print("key = ", key)
        if key == 'f':  # change focal point
            #print('focal_point!')
            return
        self.vtk_interactor._Iren.SetEventInformation(0, 0, 0, 0, key, 0, None)
        self.vtk_interactor._Iren.KeyPressEvent()
        self.vtk_interactor._Iren.CharEvent()

        #if key in ['y', 'z', 'X', 'Y', 'Z']:
            #self.update_camera(key)

    def _set_results(self, form, cases):
        assert len(cases) > 0, cases
        if isinstance(cases, OrderedDict):
            self.caseKeys = cases.keys()
        else:
            self.caseKeys = sorted(cases.keys())
            assert isinstance(cases, dict), type(cases)

        self.resultCases = cases

        if len(self.caseKeys) > 1:
            self.iCase = -1
            self.nCases = len(self.resultCases)  # number of keys in dictionary
        elif len(self.caseKeys) == 1:
            self.iCase = -1
            self.nCases = 1
        else:
            self.iCase = -1
            self.nCases = 0
        self.set_form(form)

    def _finish_results_io2(self, form, cases):
        self._set_results(form, cases)
        # assert len(cases) > 0, cases
        # if isinstance(cases, OrderedDict):
            # self.caseKeys = cases.keys()
        # else:
            # self.caseKeys = sorted(cases.keys())
            # assert isinstance(cases, dict), type(cases)

        self.on_update_geometry_properties(self.geometry_properties)
        # self.resultCases = cases

        #print("cases =", cases)
        #print("caseKeys =", self.caseKeys)

        self.reset_labels()
        self.cycleResults_explicit()  # start at nCase=0
        if self.nCases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

        #data = [
        #    ('A', []),
        #    ('B', []),
        #    ('C', []),
        #]

        #self.caseKeys= [(1, 'ElementID', 1, 'centroid', '%.0f'), (1, 'Region', 1, 'centroid', '%.0f')]
        data = []
        for key in self.caseKeys:
            print(key)
            if isinstance(key, int):
                obj, (i, name) = self.resultCases[key]
                t = (i, [])
            else:
                t = (key[1], [])
            data.append(t)

        self.res_widget.update_results(form)

        key = self.caseKeys[0]
        location = self.get_case_location(key)
        method = 'centroid' if location else 'nodal'

        data2 = [(method, None, [])]
        self.res_widget.update_methods(data2)

    def get_result_by_cell_id(self, cell_id, world_position):
        """should handle multiple cell_ids"""
        case_key = self.caseKeys[self.iCase]
        result_name = self.result_name
        result_values = self.resultCases[case_key][cell_id]
        cell = self.grid.GetCell(cell_id)

        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()
        cell_type = cell.GetCellType()

        if cell_type in [5, 9]:  # CTRIA3, CQUAD4
            node_xyz = zeros((nnodes, 3), dtype='float32')
            for ipoint in range(nnodes):
                point = points.GetPoint(ipoint)
                node_xyz[ipoint, :] = point
            xyz = node_xyz.mean(axis=0)
        elif cell_type in [10, 12, 13]: # CTETRA, CHEXA8, CPENTA6
            # TODO: No idea how to get the center of the face
            #       vs. a point on a face that's not exposed
            #faces = cell.GetFaces()
            #nfaces = cell.GetNumberOfFaces()
            #for iface in range(nfaces):
                #face = cell.GetFace(iface)
                #points = face.GetPoints()
            #faces
            xyz = world_position
        else:
            #self.log.error(msg)
            msg = 'cell_type=%s nnodes=%s' % (cell_type, nnodes)
            raise NotImplementedError(msg)
        return result_name, result_values, xyz

    def get_result_by_xyz_cell_id(self, node_xyz, cell_id):
        """won't handle multiple cell_ids/node_xyz"""
        case_key = self.caseKeys[self.iCase]
        result_name = self.result_name

        cell = self.grid.GetCell(cell_id)
        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()

        #node_xyz = array(node_xyz, dtype='float32')
        #point0 = array(points.GetPoint(0), dtype='float32')
        #dist_min = norm(point0 - node_xyz)
        point0 = points.GetPoint(0)
        dist_min = vtk.vtkMath.Distance2BetweenPoints(point0, node_xyz)

        point_min = point0
        imin = 0
        for ipoint in range(1, nnodes):
            #point = array(points.GetPoint(ipoint), dtype='float32')
            #dist = norm(point - node_xyz)
            point = points.GetPoint(ipoint)
            dist = vtk.vtkMath.Distance2BetweenPoints(point, node_xyz)
            if dist < dist_min:
                dist_min = dist
                imin = ipoint
                point_min = point

        node_id = cell.GetPointId(imin)
        xyz = array(point_min, dtype='float32')
        case = self.resultCases[case_key]
        if isinstance(case_key, (int, int32)):
            (obj, (i, res_name)) = case
            subcase_id = obj.subcase_id
            case = obj.get_result(i, res_name)
            result_values = case[node_id]
        else:
            result_values = case[node_id]
        assert not isinstance(xyz, int), xyz
        return result_name, result_values, node_id, xyz

    @property
    def result_name(self):
        """
        creates the self.result_name variable

        .. python ::

          if len(key) == 5:
              (subcase_id, result_type, vector_size, location, data_format) = key
          elif len(key) == 6:
              (subcase_id, j, result_type, vector_size, location, data_format) = key
          else:
              (subcase_id, j, result_type, vector_size, location, data_format, label2) = key
        """
        # case_key = (1, 'ElementID', 1, 'centroid', '%.0f')
        case_key = self.caseKeys[self.iCase]
        if isinstance(case_key, int):
            obj, (i, name) = self.resultCases[case_key]
            value = name
        else:
            if len(case_key) == 5:
                value = case_key[1]
            else:
                value = case_key[2]
        return value

    def finish_io(self, cases):
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        #print("caseKeys = ", self.caseKeys)

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

        if len(self.caseKeys) > 1:
            self.iCase = -1
            self.nCases = len(self.resultCases)  # number of keys in dictionary
        elif len(self.caseKeys) == 1:
            self.iCase = -1
            self.nCases = 1
        else:
            self.iCase = -1
            self.nCases = 0

        self.reset_labels()
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
        self.res_widget.update_results(data)
        # method = 'centroid' if self.is_centroidal else 'nodal'

        data2 = [('node/centroid', None, [])]
        self.res_widget.update_methods(data2)


    def reset_labels(self):
        """
        Wipe all labels and regenerate the key slots based on the case keys.
        This is used when changing the model.
        """
        self._remove_labels()

        # new geometry
        self.label_actors = {}
        self.label_ids = {}

        #self.caseKeys= [
            #(1, 'ElementID', 1, 'centroid', '%.0f'),
            #(1, 'Region', 1, 'centroid', '%.0f')
        #]
        for case_key in self.caseKeys:
            result_name = self.get_result_name(case_key)
            self.label_actors[result_name] = []
            self.label_ids[result_name] = set([])

    def _remove_labels(self):
        """
        Remove all labels from the current result case.
        This happens when the user explictly selects the clear label button.
        """
        if len(self.label_actors) == 0:
            self.log.warning('No actors to remove')
            return

        # existing geometry
        for result_name, actors in iteritems(self.label_actors):
            for actor in actors:
                self.rend.RemoveActor(actor)
                del actor
            self.label_actors[result_name] = []
            self.label_ids[result_name] = set([])

    def clear_labels(self):
        """
        This clears out all labels from all result cases.
        """
        if len(self.label_actors) == 0:
            self.log.warning('No actors to clear')
            return

        # existing geometry
        #case_key = self.caseKeys[self.iCase]
        result_name = self.result_name

        actors = self.label_actors[result_name]
        for actor in actors:
            self.rend.RemoveActor(actor)
            del actor
        self.label_actors[result_name] = []
        self.label_ids[result_name] = set([])

    def hide_labels(self, result_names=None, show_msg=True):
        if result_names is None:
            names = 'None)  # None -> all'
            result_names = sorted(self.label_actors.keys())
        else:
            mid = '%s,' * len(result_names)
            names = '[' + mid[:-1] + '])'

        count = 0
        for key in result_names:
            actors = self.label_actors[key]
            for actor in actors:
                actor.VisibilityOff()
                count += 1
        if count and show_msg:
            self.log_command('hide_labels(%s' % names)

    def show_labels(self, result_names=None, show_msg=True):
        if result_names is None:
            names = 'None)  # None -> all'
            result_names = sorted(self.label_actors.keys())
        else:
            mid = '%s,' * len(result_names)
            names = mid[:-1] % result_names + ')'

        count = 0
        for key in result_names:
            actors = self.label_actors[key]
            for actor in actors:
                actor.VisibilityOn()
                count += 1
        if count and show_msg:
            # yes the ) is intentionally left off because it's already been added
            self.log_command('show_labels(%s' % names)

    def update_scalar_bar(self, title, min_value, max_value, norm_value,
                        data_format, is_blue_to_red=True, is_horizontal=True, is_shown=True):
        """
        Parameters
        ----------
        title : str
            the scalar bar title
        min_value : float
            the blue value
        max_value :
            the red value
        data_format : str
            '%g','%f','%i', etc.
        is_blue_to_red : bool; default=True
            flips the order of the RGB points
        is_horizontal : bool; default=True
            makes the scalar bar horizontal
        is_shown : bool
            show the scalar bar
        """
        print("update_scalar_bar min=%s max=%s norm=%s" % (min_value, max_value, norm_value))
        self.scalar_bar.update(title, min_value, max_value, norm_value,
                               data_format, is_blue_to_red, is_horizontal,
                               is_shown=is_shown)

    #---------------------------------------------------------------------------------------
    # CAMERA MENU
    def view_camera(self):
        camera = self.rend.GetActiveCamera()
        position = camera.GetPosition()
        clip_range = camera.GetClippingRange()
        focal_point = camera.GetFocalPoint()

        data = {'cameras' : self.cameras}
        window = CameraWindow(data, win_parent=self)
        window.show()
        window.exec_()

        if data['clicked_ok']:
            self.cameras = deepcopy(data['cameras'])
            #self._apply_camera(data)
        #self.log_info('position = %s' % str(position))
        #self.log_info('clip_range = %s' % str(clip_range))
        #self.log_info('focal_point = %s' % str(focal_point))

    #def _apply_camera(self, data):
        #name = data['name']
        #self.cameras = deepcopy(data['cameras'])
        #self.on_set_camera(name)

    def on_set_camera(self, name, show_log=True):
        camera_data = self.cameras[name]
        #position, clip_range, focal_point, view_up, distance = camera_data
        self.on_set_camera_data(camera_data, show_log=show_log)

    def get_camera_data(self):
        camera = self.rend.GetActiveCamera()
        position = camera.GetPosition()
        focal_point = camera.GetFocalPoint()
        view_angle = camera.GetViewAngle()
        view_up = camera.GetViewUp()
        clip_range = camera.GetClippingRange()  # TODO: do I need this???

        parallel_scale = camera.GetParallelScale()  # TODO: do I need this???
        #parallel_proj = GetParralelProjection()
        parallel_proj = 32.
        distance = camera.GetDistance()

        # clip_range, view_up, distance
        camera_data = [position, focal_point, view_angle, view_up, clip_range, parallel_scale, parallel_proj, distance]
        return camera_data

    def set_camera_data(self, camera_data, show_log=True):
        #position, clip_range, focal_point, view_up, distance = camera_data
        position, focal_point, view_angle, view_up, clip_range, parallel_scale, parallel_proj, distance = camera_data

        camera = self.rend.GetActiveCamera()
        camera.SetPosition(position)
        camera.SetFocalPoint(focal_point)
        camera.SetViewAngle(view_angle)
        camera.SetViewUp(view_up)
        camera.SetClippingRange(clip_range)

        camera.SetParallelScale(parallel_scale)
        #parallel_proj

        camera.SetDistance(distance)

        camera.Modified()
        self.vtk_interactor.Render()
        if show_log:
            self.log_command('on_set_camera_data([%s, %s, %s, %s, %s, %s, %s, %s])'
                             % (position, focal_point, view_angle, view_up, clip_range, parallel_scale, parallel_proj, distance))

    #---------------------------------------------------------------------------------------
    # CLIPPING MENU
    def set_clipping(self):
        """
        Opens a dialog box to set:

        +--------+----------+
        |  Min   |  Float   |
        +--------+----------+
        |  Max   |  Float   |
        +--------+----------+
        """
        #if not hasattr(self, 'caseKeys'):  # TODO: maybe include...
            #self.log_error('No model has been loaded.')
            #return
        camera = self.GetCamera()
        min_clip, max_clip = camera.GetClippingRange()

        data = {
            'min' : min_clip,
            'max' : max_clip,
            'clicked_ok' : False,
            'close' : False,
        }
        if not self._clipping_shown:
            self._clipping_window = ClippingPropertiesWindow(data, win_parent=self)
            self._clipping_window.show()
            self._clipping_shown = True
            self._clipping_window.exec_()
        else:
            self._clipping_window.activateWindow()

        if data['close']:
            self._apply_clipping(data)
            del self._clipping_window
            self._clipping_shown = False
        else:
            self._clipping_window.activateWindow()

    def _apply_clipping(self, data):
        min_clip = data['min']
        max_clip = data['max']
        self.on_update_clipping(min_clip, max_clip)

    def on_update_clipping(self, min_clip=None, max_clip=None):
        camera = self.GetCamera()
        _min_clip, _max_clip = camera.GetClippingRange()
        if min_clip is None:
            min_clip = _min_clip
        if max_clip is None:
            max_clip = _max_clip
        camera.SetClippingRange(min_clip, max_clip)
        self.log_command('self.on_update_clipping(min_value=%s, max_clip=%s)'
                         % (min_clip, max_clip))

    #---------------------------------------------------------------------------------------
    # LEGEND MENU
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
        default_format = None
        default_scale = None
        if isinstance(key, (int, int32)):
            #(subcase_id, result_type, vector_size, location, data_format) = key
            (obj, (i, res_name)) = self.resultCases[key]
            #subcase_id = obj.subcase_id
            case = obj.get_result(i, res_name)
            result_type = obj.get_title(i, res_name)
            #vector_size = obj.get_vector_size(i, res_name)
            #location = obj.get_location(i, res_name)
            data_format = obj.get_data_format(i, res_name)
            scale = obj.get_scale(i, res_name)

            default_title = obj.get_default_title(i, res_name)
            default_scale = obj.get_default_scale(i, res_name)
        elif len(key) == 5:
            (subcase_id, result_type, vector_size, location, data_format) = key
            default_title = result_type
            scale = 0.0
        elif len(key) == 6:
            (subcase_id, i, result_type, vector_size, location, data_format) = key
            default_title = result_type
            scale = 0.0
        else:
            (subcase_id, i, result_type, vector_size, location, data_format, label2) = key
            default_title = result_type
            scale = 0.0

        if default_format is None:
            default_format = data_format
        if scale == 0.0:
            default_scale = 0.0
        print(key)

        data = {
            'name' : result_type,
            'min' : case.min(),
            'max' : case.max(),
            'scale' : scale,
            'format' : data_format,

            'default_title' : default_title,
            'default_scale' : default_scale,
            'default_format' : default_format,

            'is_blue_to_red' : True,
            'is_discrete': True,
            'is_horizontal': False,
            'is_shown' : True,
            'clicked_ok' : False,
            'close' : False,
        }
        if not self._legend_shown:
            self._legend_window = LegendPropertiesWindow(data, win_parent=self)
            self._legend_window.show()
            self._legend_shown = True
            self._legend_window.exec_()
        else:
            self._legend_window.activateWindow()

        if data['close']:
            self._apply_legend(data)
            self._legend_shown = False
            del self._legend_window
        else:
            self._legend_window.activateWindow()


    def _apply_legend(self, data):
        title = data['name']
        min_value = data['min']
        max_value = data['max']
        scale_value = data['scale']
        data_format = data['format']
        is_blue_to_red = data['is_blue_to_red']
        is_discrete = data['is_discrete']
        is_horizontal = data['is_horizontal']
        is_shown = data['is_shown']
        self.on_update_legend(Title=title, min_value=min_value, max_value=max_value, scale=scale_value,
                              data_format=data_format,
                              is_blue_to_red=is_blue_to_red,
                              is_discrete=is_discrete, is_horizontal=is_horizontal,
                              is_shown=is_shown)

    def on_update_legend(self, Title='Title', min_value=0., max_value=1., scale=0.0,
                         data_format='%.0f',
                         is_blue_to_red=True, is_discrete=True, is_horizontal=True,
                         is_shown=True):


        key = self.caseKeys[self.iCase]
        name_vector = None
        plot_value = self.resultCases[key] # scalar
        vector_size1 = 1
        update_3d = False
        if isinstance(key, (int, int32)):
            #(subcase_id, result_type, vector_size, location, data_format) = key
            (obj, (i, res_name)) = self.resultCases[key]
            subcase_id = obj.subcase_id
            plot_value = obj.get_plot_value(i, res_name) # vector

            result_type = obj.get_title(i, res_name)
            vector_size = obj.get_vector_size(i, res_name)
            scalar = obj.get_scalar(i, res_name)

            location = obj.get_location(i, res_name)

            #data_format = obj.get_data_format(i, res_name)
            obj.set_scale(i, res_name, scale)
            #obj.set_format(i, res_name, data_format)
            #obj.set_data_format(i, res_name, data_format)
            subtitle, label = self.get_subtitle_label(subcase_id)
            name_vector = (vector_size1, subcase_id, result_type, label, min_value, max_value, scale)
            update_3d = True
        elif len(key) == 5:
            (subcase_id, result_type, vector_size1, location, _data_format) = key
        elif len(key) == 6:
            (subcase_id, i, result_type, vector_size1, location, _data_format) = key
        else:
            (subcase_id, i, result_type, vector_size1, location, _data_format, label2) = key
        assert vector_size1 == 1, vector_size1

        if update_3d:
            self.is_horizontal_scalar_bar = is_horizontal
            self._set_case(self.result_name, self.iCase,
                           explicit=False, cycle=False, skip_click_check=True)
            return

        subtitle, label = self.get_subtitle_label(subcase_id)
        scale1 = 0.0
        name = (vector_size1, subcase_id, result_type, label, min_value, max_value, scale1)
        # if vector_size == 3:

        norm_value = float(max_value - min_value)
        # if name not in self._loaded_names:
        grid_result = self.set_grid_values(name, plot_value, vector_size1,
                                           min_value, max_value, norm_value,
                                           is_blue_to_red=is_blue_to_red)

        grid_result_vector = None
        if name_vector and 0:
            vector_size = 3
            grid_result_vector = self.set_grid_values(name_vector, plot_value, vector_size,
                                                      min_value, max_value, norm_value,
                                                      is_blue_to_red=is_blue_to_red)

        self.update_scalar_bar(Title, min_value, max_value, norm_value,
                               data_format, is_blue_to_red=is_blue_to_red,
                               is_horizontal=is_horizontal, is_shown=is_shown)

        revert_displaced = True
        self._final_grid_update(name, grid_result, None, None, None,
                                1, subcase_id, result_type, location, subtitle, label,
                                revert_displaced=revert_displaced)
        if grid_result_vector is not None:
            self._final_grid_update(name_vector, grid_result_vector, obj, i, res_name,
                                    vector_size, subcase_id, result_type, location, subtitle, label,
                                    revert_displaced=False)
            if 0:
                xyz_nominal, vector_data = obj.get_vector_result(i, res_name)
                self._update_grid(vector_data)
                self.grid.Modified()
                self.geom_actor.Modified()
                self.vtk_interactor.Render()
            #revert_displaced = False
        #self._final_grid_update(name, grid_result, None, None, None,
                                #1, subcase_id, result_type, location, subtitle, label,
                                #revert_displaced=revert_displaced)

        #self.is_horizontal_scalar_bar = is_horizontal
        self.log_command('self.on_update_legend(Title=%r, min_value=%s, max_value=%s,\n'
                         '                      data_format=%r, is_blue_to_red=%s, is_discrete=%s)'
                         % (Title, min_value, max_value, data_format, is_blue_to_red, is_discrete))

    #---------------------------------------------------------------------------------------
    # EDIT ACTOR PROPERTIES
    def set_actor_properties(self):
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
        if not len(self.geometry_properties):
            self.log_error('No secondary geometries to edit.')
            return

        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]
        #if len(key) == 5:
            #(subcase_id, result_type, vector_size, location, data_format) = key
        #elif len(key) == 6:
            #(subcase_id, i, result_type, vector_size, location, data_format) = key
        #else:
            #(subcase_id, i, result_type, vector_size, location, data_format, label2) = key

        data = deepcopy(self.geometry_properties)
        if not self._edit_group_properties_shown:
            self._edit_group_properties = EditGroupProperties(data, win_parent=self)
            self._edit_group_properties.show()
            self._edit_group_properties_shown = True
            self._edit_group_properties.exec_()
        else:
            self._edit_group_properties.activateWindow()

        if  'clicked_ok' not in data:
            self._edit_group_properties.activateWindow()

        if data['clicked_ok']:
            self.on_update_geometry_properties(data)
            self._save_geometry_properties(data)
            del self._edit_group_properties
            self._edit_group_properties_shown = False
        elif data['clicked_cancel']:
            self.on_update_geometry_properties(self.geometry_properties)
            del self._edit_group_properties
            self._edit_group_properties_shown = False

    def _save_geometry_properties(self, out_data):
        for name, group in iteritems(out_data):
            if name in ['clicked_ok', 'clicked_cancel']:
                continue

            #color2 = group.color_float
            geom_prop = self.geometry_properties[name]
            geom_prop.color = group.color
            geom_prop.line_width = group.line_width
            geom_prop.opacity = group.opacity
            geom_prop.point_size = group.point_size

    def on_update_geometry_properties(self, out_data):
        lines = []
        for name, group in iteritems(out_data):
            if name in ['clicked_ok', 'clicked_cancel']:
                continue
            changed = False
            actor = self.geometry_actors[name]
            #mapper = actor.GetMapper()
            prop = actor.GetProperty()

            color1 = prop.GetDiffuseColor()
            assert color1[1] <= 1.0, color1
            color2 = group.color_float
            #print('line2646 - name=%s color1=%s color2=%s' % (name, str(color1), str(color2)))
            #color2 = group.color

            line_width1 = prop.GetLineWidth()
            line_width2 = group.line_width
            line_width2 = max(1, line_width2)

            opacity1 = prop.GetOpacity()
            opacity2 = group.opacity
            opacity2 = max(0.1, opacity2)

            point_size1 = prop.GetPointSize()
            point_size2 = group.point_size
            point_size2 = max(1, point_size2)

            #representation = group.representation
            alt_prop = self.geometry_properties[name]
            representation = alt_prop.representation
            #is_visible1 = alt_prop.is_visible
            is_visible1 = bool(actor.GetVisibility())
            is_visible2 = group.is_visible
            #print('is_visible1=%s is_visible2=%s'  % (is_visible1, is_visible2))

            bar_scale1 = alt_prop.bar_scale
            bar_scale2 = group.bar_scale
            # bar_scale2 = max(0.0, bar_scale2)

            if color1 != color2:
                #print('color_2662[%s] = %s' % (name, str(color1)))
                assert isinstance(color1[0], float), color1
                prop.SetDiffuseColor(color2)
                changed = True
            if line_width1 != line_width2:
                line_width2 = max(1, line_width2)
                prop.SetLineWidth(line_width2)
                changed = True
            if opacity1 != opacity2:
                prop.SetOpacity(opacity2)
                changed = True
            if point_size1 != point_size2:
                prop.SetPointSize(point_size2)
                changed = True
            if bar_scale1 != bar_scale2 and bar_scale1 > 0.0:
                print('name=%s bar_scale1=%s bar_scale2=%s' % (name, bar_scale1, bar_scale2))
                self.set_bar_scale(name, bar_scale2)
            if is_visible1 != is_visible2:
                actor.SetVisibility(is_visible2)
                alt_prop.is_visible = is_visible2
                #prop.SetViPointSize(is_visible2)
                actor.Modified()
                changed = True

            if changed:
                lines.append('    %r : AltGeometry(color=(%s, %s, %s), line_width=%s, opacity=%s, point_size=%s, is_visible=%s),\n' % (
                    name, color2[0], color2[1], color2[2], line_width2, opacity2, point_size2, is_visible2))
                prop.Modified()
        self.vtk_interactor.Render()
        if lines:
            msg = 'out_data = {\n'
            msg += ''.join(lines)
            msg += '}\n'
            msg += 'self.on_update_geometry_properties(out_data)'
            self.log_command(msg)

    def set_bar_scale(self, name, bar_scale):
        print('set_bar_scale; name=%s scale=%s' % (name, bar_scale))
        bar_y = self.bar_lines[name]
        #dy = c - yaxis
        #dz = c - zaxis
        n1 = bar_y[:, :3]
        n2 = bar_y[:, 3:]
        dy = n2 - n1
        # for
        Ly = norm(dy, axis=1)
        # v = dy / Ly *  bar_scale
        # n2 = n1 + v
        print(Ly)
        nnodes = len(Ly)
        points = self.alt_grids[name].GetPoints()
        for i in range(nnodes):
            node = n1[i, :] + Ly[i] * bar_scale * dy[i, :]
            points.SetPoint(2 * i + 1, *node)
        self.alt_grids[name].Update()
        # bar_z = self.bar_lines[(name)]
        # dz = bar_z[:, :3] - bar_z[:, 3:]
        # Lz = norm(dz, axis=1)
        # bar_z[:, 3:] = bar_z[:, :3] + Lz * bar_scale


    def add_user_points(self, points_filename, name, color=None):
        if name in self.geometry_actors:
            msg = 'Name: %s is already in geometry_actors\nChoose a different name.' % name
            raise ValueError(msg)

        if color is None:
            color = (0., 1., .1)

        # create grid
        self.create_alternate_vtk_grid(name, color=color, line_width=5, opacity=1.0,
                                       point_size=1, representation='point')

        # read input file
        try:
            user_points = loadtxt(points_filename, delimiter=',')
        except ValueError:
            user_points = loadtxt_nice(points_filename, delimiter=',')
            # can't handle leading spaces?
            #raise
        npoints = user_points.shape[0]

        # allocate grid
        self.alt_grids[name].Allocate(npoints, 1000)

        # set points
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(npoints)

        for i, point in enumerate(user_points):
            points.InsertPoint(i, *point)
            elem = vtk.vtkVertex()
            elem.GetPointIds().SetId(0, i)
            self.alt_grids[name].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
        self.alt_grids[name].SetPoints(points)

        # create actor/mapper
        self._add_alt_geometry(self.alt_grids[name], name)

        # set representation to points
        self.geometry_properties[name].representation = 'point'
        actor = self.geometry_actors[name]
        prop = actor.GetProperty()
        prop.SetRepresentationToPoints()
        prop.SetPointSize(4)
