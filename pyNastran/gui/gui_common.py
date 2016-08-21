# -*- coding: utf-8 -*-
# pylint: disable=W0201,C0111
from __future__ import division, unicode_literals, print_function

# standard library
import sys
import os.path
import datetime
import cgi #  html lib
import traceback
from copy import deepcopy
from collections import OrderedDict

from six import string_types, iteritems, itervalues, PY2
from six.moves import range

import numpy as np
#from numpy import arange
#from numpy import eye, array, zeros, loadtxt
#from numpy.linalg import norm

from PyQt4 import QtCore, QtGui
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
#from vtk.util.numpy_support import numpy_to_vtk


import pyNastran
from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.utils.log import SimpleLogger
from pyNastran.utils import print_bad_path, integer_types
from pyNastran.utils.numpy_utils import loadtxt_nice

from pyNastran.gui.qt_files.gui_qt_common import GuiCommon
from pyNastran.gui.qt_files.scalar_bar import ScalarBar
from pyNastran.gui.qt_files.alt_geometry_storage import AltGeometry

from pyNastran.gui.menus.results_sidebar import Sidebar
from pyNastran.gui.menus.qt_legend import LegendPropertiesWindow
from pyNastran.gui.menus.clipping import ClippingPropertiesWindow
from pyNastran.gui.menus.camera import CameraWindow
from pyNastran.gui.menus.application_log import PythonConsoleWidget, ApplicationLogWidget
from pyNastran.gui.menus.manage_actors import EditGeometryProperties
from pyNastran.gui.menus.groups_modify import GroupsModify, Group
from pyNastran.gui.menus.modify_label_properties import ModifyLabelPropertiesMenu
from pyNastran.gui.menus.modify_picker_properties import ModifyPickerPropertiesMenu

from pyNastran.gui.testing_methods import CoordProperties
#from pyNastran.gui.menus.multidialog import MultiFileDialog
from pyNastran.gui.utils import load_csv, load_user_geom


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


class GuiCommon2(QtGui.QMainWindow, GuiCommon):
    def __init__(self, fmt_order, html_logging, inputs):
        # this will reset the background color/label color if things break
        self.reset_settings = False

        QtGui.QMainWindow.__init__(self)
        GuiCommon.__init__(self, inputs)

        self.fmts = fmt_order
        self.base_window_title = "pyNastran v%s"  % pyNastran.__version__

        # initializes tools/checkables
        self.set_tools()

        self.html_logging = html_logging
        self.execute_python = True
        self.scalar_bar = ScalarBar(self.is_horizontal_scalar_bar)
        # in,lb,s
        self.input_units = ['', '', ''] # '' means not set
        self.display_units = ['', '', '']

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
    def legend_shown(self):
        return self.scalar_bar.is_shown

    @property
    def scalarBar(self):
        return self.scalar_bar.scalar_bar

    @property
    def color_function(self):
        return self.scalar_bar.color_function

    #def get_color_function(self):
        #return self.scalar_bar.color_function

    def set_window_title(self, msg):
        #msg2 = "%s - "  % self.base_window_title
        #msg2 += msg
        self.setWindowTitle(msg)

    def set_logo(self, logo):
        """Sets the pyNastran icon path, which can be overwritten"""
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
        #self.resize(1100, 700)
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
        self.create_log_python_docks()
        #===============================================

        self.run_vtk = True
        if self.run_vtk:
            self._create_vtk_objects()
        self._build_menubar()
        #self._hide_menubar()

        # right sidebar
        self.res_dock.hide()
        if self.run_vtk:
            self.build_vtk_frame()

        #compassRepresentation = vtk.vtkCompassRepresentation()
        #compassWidget = vtk.vtkCompassWidget()
        #compassWidget.SetInteractor(self.iren)
        #compassWidget.SetRepresentation(compassRepresentation)
        #compassWidget.EnabledOn()

    def create_log_python_docks(self):
        """
        Creates the
         - HTML Log dock
         - Python Console dock
        """
        #=========== Logging widget ===================

        if self.html_logging:
            self.log_dock_widget = ApplicationLogWidget(self)
            self.log_widget = self.log_dock_widget.log_widget
            #self.log_dock_widget = QtGui.QDockWidget("Application log", self)
            #self.log_dock_widget.setObjectName("application_log")
            #self.log_widget = QtGui.QTextEdit()
            #self.log_widget.setReadOnly(True)
            #self.log_dock_widget.setWidget(self.log_widget)
            self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock_widget)

        if self.execute_python:
            self.python_dock_widget = PythonConsoleWidget(self)
            self.python_dock_widget.setObjectName("python_console")
            self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.python_dock_widget)

    def _on_execute_python_button(self, clear=False):
        """executes the docked python console"""
        txt = str(self.python_dock_widget.enter_data.toPlainText()).rstrip()
        if len(txt) == 0:
            return
        self.log_command(txt)
        try:
            exec(txt)
        except TypeError:
            print(type(txt))
            raise
        except Exception as e:
            #self.log_error(traceback.print_stack(f))
            self.log_error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)
            self.log_error(str(e))
            self.log_error(str(txt))
            raise
            return
        if clear:
            self.python_dock_widget.enter_data.clear()

    def load_batch_inputs(self, inputs):
        geom_script = inputs['geomscript']
        if geom_script is not None:
            self.on_run_script(geom_script)

        if not inputs['format']:
            return
        form = inputs['format'].lower()
        input_filenames = inputs['input']
        results_filename = inputs['output']
        plot = True
        if results_filename:
            plot = False

        #print('input_filename =', input_filename)
        if input_filenames is not None:
            for input_filename in input_filenames:
                if not os.path.exists(input_filename):
                    msg = '%s does not exist\n%s' % (
                        input_filename, print_bad_path(input_filename))
                    self.log.error(msg)
                    if self.html_logging:
                        print(msg)
                    return
            for results_filenamei in results_filename:
                #print('results_filenamei =', results_filenamei)
                if results_filenamei is not None:
                    if not os.path.exists(results_filenamei):
                        msg = '%s does not exist\n%s' % (
                            results_filenamei, print_bad_path(results_filenamei))
                        self.log.error(msg)
                        if self.html_logging:
                            print(msg)
                        return

        name = 'main'
        for i, input_filename in enumerate(input_filenames):
            if i == 0:
                name = 'main'
            else:
                name = input_filename
            print('name =', name)
            self.name = name
            #form = inputs['format'].lower()
            is_failed = self.on_load_geometry(
                infile_name=input_filename, name=name, geometry_format=form, plot=plot)
        self.name = 'main'
        print('keys =', self.nid_maps.keys())

        if is_failed:
            return
        if results_filename:
            self.on_load_results(results_filename)

        post_script = inputs['postscript']
        if post_script is not None:
            self.on_run_script(post_script)
        self.on_reset_camera()
        self.vtk_interactor.Modified()

    def set_tools(self, tools=None, checkables=None):
        """Creates the GUI tools"""
        if checkables is None:
            checkables = ['show_info', 'show_debug', 'show_gui', 'show_command']
        if tools is None:
            file_tools = [

                ('exit', '&Exit', 'texit.png', 'Ctrl+Q', 'Exit application', self.closeEvent), # QtGui.qApp.quit
                ('load_geometry', 'Load &Geometry', 'load_geometry.png', 'Ctrl+O', 'Loads a geometry input file', self.on_load_geometry),
                ('load_results', 'Load &Results', 'load_results.png', 'Ctrl+R', 'Loads a results file', self.on_load_results),

                ('load_csv_user_geom', 'Load CSV User Geometry', '', None, 'Loads custom geometry file', self.on_load_user_geom),
                ('load_csv_user_points', 'Load CSV User Points', 'user_points.png', None, 'Loads user defined points ', self.on_load_user_points),

                ('load_csv_nodal', 'Load CSV Nodal Results', '', None, 'Loads a custom nodal results file', self.on_load_nodal_results),
                ('load_csv_elemental', 'Load CSV Elemental Results', '', None, 'Loads a custom elemental results file', self.on_load_elemental_results),
                ('script', 'Run Python script', 'python48.png', None, 'Runs pyNastranGUI in batch mode', self.on_run_script),
            ]

            tools = file_tools + [
                ('back_color', 'Change background color', 'tcolorpick.png', None, 'Choose a background color', self.change_background_color),
                #('label_color', 'Change label color', 'tcolorpick.png', None, 'Choose a label color', self.change_label_color),
                ('text_color', 'Change text color', 'tcolorpick.png', None, 'Choose a text color', self.change_text_color),

                ('label_clear', 'Clear current labels', '', None, 'Clear current labels', self.clear_labels),
                ('label_modify', 'Modify label color/size', '', None, 'Edit Label Properties', self.on_set_labelsize_color),
                ('label_reset', 'Clear all labels', '', None, 'Clear all labels', self.reset_labels),

                ('picker_modify', 'Modify picker size', '', None, 'Edit Label Properties', self.on_set_picker_size),

                ('legend', 'Modify legend', 'legend.png', None, 'Set Legend', self.set_legend),
                ('clipping', 'Set clipping', '', None, 'Set Clipping', self.set_clipping),
                #('axis', 'Show/Hide Axis', 'axis.png', None, 'Show/Hide Global Axis', self.on_show_hide_axes),

                ('wireframe', 'Wireframe Model', 'twireframe.png', 'w', 'Show Model as a Wireframe Model', self.on_wireframe),
                ('surface', 'Surface Model', 'tsolid.png', 's', 'Show Model as a Surface Model', self.on_surface),
                ('geo_properties', 'Edit Geometry Properties', '', None, 'Change Model Color/Opacity/Line Width', self.edit_geometry_properties),
                ('modify_groups', 'Modify Groups', '', None, 'Create/Edit/Delete Groups', self.modify_group),

                ('create_groups_by_property_id', 'Create Groups By Property ID', '', None, 'Create Groups', self.create_groups_by_property_id),
                #('create_list', 'Create Lists through Booleans', '', None, 'Create List', self.create_list),

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

                ('screenshot', 'Take a Screenshot', 'tcamera.png', 'CTRL+I', 'Take a Screenshot of current view', self.on_take_screenshot),
                ('about', 'About pyNastran GUI', 'tabout.png', 'CTRL+H', 'About pyNastran GUI and help on shortcuts', self.about_dialog),
                ('view', 'Camera View', 'view.png', None, 'Load the camera menu', self.view_camera),
                ('camera_reset', 'Reset camera view', 'trefresh.png', 'r', 'Reset the camera view to default', self.on_reset_camera),
                ('reload', 'Reload model', 'treload.png', 'r', 'Remove the model and reload the same geometry file', self.on_reload),

                ('cycle_results', 'Cycle Results', 'cycle_results.png', 'CTRL+L', 'Changes the result case', self.cycle_results),

                ('x', 'Flips to +X Axis', 'plus_x.png', 'x', 'Flips to +X Axis', lambda: self.update_camera('+x')),
                ('y', 'Flips to +Y Axis', 'plus_y.png', 'y', 'Flips to +Y Axis', lambda: self.update_camera('+y')),
                ('z', 'Flips to +Z Axis', 'plus_z.png', 'z', 'Flips to +Z Axis', lambda: self.update_camera('+z')),

                ('X', 'Flips to -X Axis', 'minus_x.png', 'X', 'Flips to -X Axis', lambda: self.update_camera('-x')),
                ('Y', 'Flips to -Y Axis', 'minus_y.png', 'Y', 'Flips to -Y Axis', lambda: self.update_camera('-y')),
                ('Z', 'Flips to -Z Axis', 'minus_z.png', 'Z', 'Flips to -Z Axis', lambda: self.update_camera('-z')),
                ('edges', 'Show/Hide Edges', 'tedges.png', 'e', 'Show/Hide Model Edges', self.on_flip_edges),
                ('edges_black', 'Color Edges', '', 'b', 'Set Edge Color to Color/Black', self.on_set_edge_visibility),
            ]
        # print('version =', vtk.VTK_VERSION, self.vtk_version)
        #if self.vtk_version[0] < 6

        if 'nastran' in self.fmts:
            tools += [
                ('caero', 'Show/Hide CAERO Panels', '', None, 'Show/Hide CAERO Panel Outlines', self.toggle_caero_panels),
                ('caero_subpanels', 'Toggle CAERO Subpanels', '', None, 'Show/Hide CAERO Subanel Outlines', self.toggle_caero_sub_panels),
                ('conm2', 'Toggle CONM2s', '', None, 'Show/Hide CONM2s', self.toggle_conms),
            ]
        self.tools = tools
        self.checkables = checkables

    def deprecated(self, old_name, new_name, deprecated_version):
        deprecated(old_name, new_name, deprecated_version, levels=[-1])

    #def add_tools(self, tools):
        #self.deprecated('add_tools', 'removed...', '0.7')
        #self.tools += tools

    def on_flip_picker(self):
        return
        # if self.pick_state == 'centroidal':
            # self.pick_state = 'nodal'

        # elif self.pick_state == 'nodal':
            # self.pick_state = 'centroidal'
        # else:
            # raise RuntimeError(self.pick_state)
        # self.log_command("on_flip_pick() # pick_state='%s'" % self.pick_state)

    def _create_menu_items(self, actions=None, create_menu_bar=True):
        if actions is None:
            actions = self.actions

        if create_menu_bar:
            self.menu_file = self.menubar.addMenu('&File')
            self.menu_view = self.menubar.addMenu('&View')
            self.menu_window = self.menubar.addMenu('&Window')
            self.menu_help = self.menubar.addMenu('&Help')

            self.menu_hidden = self.menubar.addMenu('&Hidden')
            self.menu_hidden.menuAction().setVisible(False)

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
            'screenshot', '', 'wireframe', 'surface', 'camera_reset', '',
            'back_color', 'text_color', '',
            'label_modify', 'label_clear', 'label_reset', 'picker_modify', '',
            'legend', 'geo_properties']
        if self.is_groups:
            menu_view += ['modify_groups', 'create_groups_by_property_id']
        menu_view += [
            '', 'clipping', #'axis',
            'edges', 'edges_black',]
        if self.html_logging:
            self.actions['log_dock_widget'] = self.log_dock_widget.toggleViewAction()
            self.actions['log_dock_widget'].setStatusTip("Show/Hide application log")
            menu_view += ['', 'show_info', 'show_debug', 'show_gui', 'show_command']
            menu_window += ['log_dock_widget']
        if self.execute_python:
            self.actions['python_dock_widget'] = self.python_dock_widget.toggleViewAction()
            self.actions['python_dock_widget'].setStatusTip("Show/Hide Python Console")
            menu_window += ['python_dock_widget']

        menu_file = [
            'load_geometry', 'load_results', 'load_csv_nodal', 'load_csv_elemental',
            'load_csv_user_points', 'load_csv_user_geom', 'script', '', 'exit']
        toolbar_tools = ['reload', 'load_geometry', 'load_results',
                         'x', 'y', 'z', 'X', 'Y', 'Z',
                         'magnify', 'shrink', 'rotate_clockwise', 'rotate_cclockwise',
                         'wireframe', 'surface', 'edges']
        toolbar_tools += ['camera_reset', 'view', 'screenshot', '', 'exit']

        menu_items = []
        if create_menu_bar:
            menu_items = [
                (self.menu_file, menu_file),
                (self.menu_view, menu_view),
                (self.menu_window, menu_window),
                (self.menu_help, ('about',)),
                (self.menu_scripts, scripts),
                (self.toolbar, toolbar_tools),
                (self.menu_hidden, ('cycle_results',)),
                # (self.menu_scripts, ()),
                #(self._dummy_toolbar, ('cell_pick', 'node_pick'))
            ]
        return menu_items

    def _hide_menubar(self):
        self.toolbar.setVisible(False)
        #self.menuBar.setVisible(False)

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

    #def _create_plane_from_points(self, points):
        #origin, vx, vy, vz, x_limits, y_limits = self._fit_plane(points)

        ## We create a 100 by 100 point plane to sample
        #splane = vtk.vtkPlaneSource()
        #plane = splane.GetOutput()

        #dx = max(x_limits) - min(x_limits)
        #dy = max(y_limits) - min(y_limits)
        ##dx = 1.
        ##dy = 3.

        ## we need to offset the origin of the plane because the "origin"
        ## is at the lower left corner of the plane and not the centroid
        #offset = (dx * vx + dy * vy) / 2.
        #origin -= offset
        #splane.SetCenter(origin)

        #splane.SetNormal(vz)

        ## Point 1 defines the x-axis and the x-size
        ## Point 2 defines the y-axis and the y-size
        #splane.SetPoint1(origin + dx * vx)
        #splane.SetPoint2(origin + dy * vy)

        #actor = vtk.vtkLODActor()
        #mapper = vtk.vtkPolyDataMapper()
        ##mapper.InterpolateScalarsBeforeMappingOn()
        ##mapper.UseLookupTableScalarRangeOn()

        #if self.vtk_version <= 5:
            #mapper.SetInputData(plane)
        #else:
            #mapper.SetInput(plane)

        #actor.GetProperty().SetColor(1., 0., 0.)
        #actor.SetMapper(mapper)
        #self.rend.AddActor(actor)
        #splane.Update()

    #def _fit_plane(self, points):
        #origin = np.array([34.60272856552356, 16.92028913186242, 37.805958003209184])
        #vx = np.array([1., 0., 0.])
        #vy = np.array([0., 1., 0.])
        #vz = np.array([0., 0., 1.])
        #x_limits = [-1., 2.]
        #y_limits = [0., 1.]
        #return origin, vx, vy, vz, x_limits, y_limits

    def _prepare_actions(self, icon_path, tools, checkables=None):
        """
        Prepare actions that will  be used in application in a way
        that's independent of the  menus & toolbar
        """
        if checkables is None:
            checkables = []
        #print('---------------------------')
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
        """ Helper funtion: log a messaage msg with a 'INFO:' prefix """
        assert msg is not None, msg
        self.log.simple_msg(msg, 'INFO')

    def log_debug(self, msg):
        """ Helper funtion: log a messaage msg with a 'DEBUG:' prefix """
        assert msg is not None, msg
        self.log.simple_msg(msg, 'DEBUG')

    def log_command(self, msg):
        """ Helper funtion: log a messaage msg with a 'COMMAND:' prefix """
        assert msg is not None, msg
        self.log.simple_msg(msg, 'COMMAND')

    def log_error(self, msg):
        """ Helper funtion: log a messaage msg with a 'GUI ERROR:' prefix """
        assert msg is not None, msg
        self.log.simple_msg(msg, 'GUI ERROR')

    def log_warning(self, msg):
        """ Helper funtion: log a messaage msg with a 'WARNING:' prefix """
        assert msg is not None, msg
        self.log.simple_msg(msg, 'WARNING')

    def change_background_color(self):
        """ Choose a background color """
        self._change_color('background', self.background_color, self.set_background_color)

    def change_label_color(self):
        """ Choose a label color """
        self._change_color('label', self.label_color, self.set_label_color)

    def change_text_color(self):
        """ Choose a text color """
        self._change_color('text', self.text_color, self.set_text_color)

    def _change_color(self, msg, rgb_color_floats, call_func):
        c = [int(255 * i) for i in rgb_color_floats]
        col = QtGui.QColorDialog.getColor(QtGui.QColor(*c), self, "Choose a %s color" % msg)
        if col.isValid():
            color = col.getRgbF()[:3]
            call_func(color)

    def set_background_color(self, color):
        """
        Set the background color

        Parameters
        ----------
        color : (float, float, float)
            RGB values as floats
        """
        self.background_color = color
        self.rend.SetBackground(*color)
        self.log_command('set_background_color(%s, %s, %s)' % color)

    def set_text_color(self, color):
        """Set the text color"""
        self.text_color = color
        for text_actor in itervalues(self.text_actors):
            text_actor.GetTextProperty().SetColor(color)
        self.log_command('set_text_color(%s, %s, %s)' % color)

    def create_coordinate_system(self, dim_max, label='', origin=None, matrix_3x3=None, Type='xyz'):
        """
        Creates a coordinate system

        Parameters
        ----------
        dim_max : float
            the max model dimension; 10% of the max will be used for the coord length
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
        self.dim_max = dim_max
        scale = 0.05 * dim_max

        transform = vtk.vtkTransform()
        if origin is None and matrix_3x3 is None:
            pass
        elif origin is not None and matrix_3x3 is None:
            #print('origin%s = %s' % (label, str(origin)))
            transform.Translate(*origin)
        elif matrix_3x3 is not None:  # origin can be None
            m = np.eye(4, dtype='float32')
            m[:3, :3] = matrix_3x3
            if origin is not None:
                m[:3, 3] = origin
            transform.SetMatrix(m.ravel())
        else:
            raise RuntimeError('unexpected coordinate system')

        axes = vtk.vtkAxesActor()
        axes.DragableOff()
        axes.PickableOff()
        #axes.GetLength() # pi
        #axes.GetNormalizedShaftLength() # (0.8, 0.8, 0.8)
        #axes.GetNormalizedTipLength() # (0.2, 0.2, 0.2)
        #axes.GetOrigin() # (0., 0., 0.)
        #axes.GetScale() # (1., 1., 1.)
        #axes.GetShaftType() # 1
        #axes.GetTotalLength() # (1., 1., 1.)

        axes.SetUserTransform(transform)
        axes.SetTotalLength(scale, scale, scale)
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

        is_visible = False
        if label == '':
            label = 'Global XYZ'
            is_visible = True
        else:
            label = 'Coord %s' % label
        self.geometry_properties[label] = CoordProperties(label, Type, is_visible, scale)
        self.geometry_actors[label] = axes
        self.coord_id += 1
        self.rend.AddActor(axes)
        return self.coord_id

    def create_global_axes(self, dim_max):
        self.create_coordinate_system(dim_max, label='', origin=None, matrix_3x3=None, Type='xyz')

    def create_corner_axis(self):
        if not self.run_vtk:
            return
        axes = vtk.vtkAxesActor()
        self.corner_axis = vtk.vtkOrientationMarkerWidget()
        self.corner_axis.SetOrientationMarker(axes)
        self.corner_axis.SetInteractor(self.vtk_interactor)
        self.corner_axis.SetEnabled(1)
        self.corner_axis.InteractiveOff()

    #def on_show_hide_axes(self):
        #"""
        #show/hide axes
        #"""
        #if not self.run_vtk:
            #return
        ## this method should handle all the coords when
        ## there are more then one
        #if self._is_axes_shown:
            #for axis in itervalues(self.axes):
                #axis.VisibilityOff()
        #else:
            #for axis in itervalues(self.axes):
                #axis.VisibilityOn()
        #self._is_axes_shown = not self._is_axes_shown

    def create_vtk_actors(self):
        self.rend = vtk.vtkRenderer()

        # vtk actors
        self.grid = vtk.vtkUnstructuredGrid()
        #self.emptyResult = vtk.vtkFloatArray()
        #self.vectorResult = vtk.vtkFloatArray()

        # edges
        self.edge_actor = vtk.vtkLODActor()
        self.edge_actor.DragableOff()
        self.edge_mapper = vtk.vtkPolyDataMapper()

        self.create_cell_picker()

    def create_alternate_vtk_grid(self, name, color=None, line_width=5, opacity=1.0, point_size=1,
                                  bar_scale=0.0, representation=None, is_visible=True):
        self.alt_grids[name] = vtk.vtkUnstructuredGrid()
        self.geometry_properties[name] = AltGeometry(
            self, name, color=color,
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
        self.add_geometry()
        if nframes == 2:
            rend.AddActor(self.geom_actor)

        # initialize geometry_actors
        self.geometry_actors['main'] = self.geom_actor
        #self.geometry_actors['anti-main'] = self.not_selected_actor

        # bar scale set so you can't edit the bar scale
        white = (255, 255, 255)
        geom_props = AltGeometry(
            self, 'main', color=white, line_width=1, opacity=1.0, point_size=1,
            bar_scale=0.0, representation='main', is_visible=True)
        anti_geom_props = AltGeometry(
            self, 'anti-main', color=white, line_width=1, opacity=0.3, point_size=1,
            bar_scale=0.0, representation='main', is_visible=True)

        self.geometry_properties['main'] = geom_props
        #self.geometry_properties['anti-main'] = anti_geom_props

        #self.addAltGeometry()
        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(*self.background_color)

        self.rend.ResetCamera()
        self._simulate_key_press('t') # change mouse style to trackball
        self.build_lookup_table()

        text_size = 14
        self.create_text([5, 50], 'Max  ', text_size)  # text actor 0
        self.create_text([5, 35], 'Min  ', text_size)  # text actor 1
        self.create_text([5, 20], 'Word1', text_size)  # text actor 2
        self.create_text([5, 5], 'Word2', text_size)  # text actor 3

        self.get_edges()
        if self.is_edges:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOn()
        else:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOff()

    #def _script_helper(self, python_file=False):
        #if python_file in [None, False]:
            #self.on_run_script(python_file)

    def on_run_script(self, python_file=False):
        if python_file in [None, False]:
            title = 'Choose a Python Script to Run'
            wildcard = "Python (*.py)"
            infile_name = self._create_load_file_dialog(
                wildcard, title, self._default_python_file)[1]
            if not infile_name:
                is_failed = True
                return is_failed # user clicked cancel

            #python_file = os.path.join(script_path, infile_name)
            python_file = os.path.join(infile_name)

        print('python_file =', python_file)
        execfile(python_file)
        self._default_python_file = python_file
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
                    #print('name: %s\nrep: %s' % (
                        #name, self.geometry_properties[name].representation))
                representation = self.geometry_properties[name].representation
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
                    #print('name: %s\nrep: %s' % (
                        #name, self.geometry_properties[name].representation))
                representation = self.geometry_properties[name].representation
                if name == 'main' or representation in ['main', 'toggle']:
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
        """rotate clockwise"""
        self.rotate(15.0)

    def on_rotate_cclockwise(self):
        """rotate counter clockwise"""
        self.rotate(-15.0)

    def on_increase_magnification(self):
        """zoom in"""
        self.zoom(1.1)

    def on_decrease_magnification(self):
        """zoom out"""
        self.zoom(1.0 / 1.1)

    def on_flip_edges(self):
        """turn edges on/off"""
        self.is_edges = not self.is_edges
        self.edge_actor.SetVisibility(self.is_edges)
        #self.edge_actor.GetProperty().SetColor(0, 0, 0)  # cart3d edge color isn't black...
        self.edge_actor.Modified()
        #self.widget.Update()
        self._update_camera()
        #self.refresh()
        self.log_command('on_flip_edges()')

    def on_set_edge_visibility(self):
        #self.edge_actor.SetVisibility(self.is_edges_black)
        self.is_edges_black = not self.is_edges_black
        if self.is_edges_black:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOn()
        else:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOff()
        self.edge_actor.Modified()
        prop.Modified()
        self.vtk_interactor.Render()
        self.log_command('on_set_edge_visibility()')

    def get_edges(self):
        """Create the edge actor"""
        edges = vtk.vtkExtractEdges()
        if self.vtk_version[0] >= 6:
            # new
            edges.SetInputData(self.grid_selected)
            self.edge_mapper.SetInputConnection(edges.GetOutputPort())
        else:
            edges.SetInput(self.grid_selected)
            self.edge_mapper.SetInput(edges.GetOutput())

        self.edge_actor.SetMapper(self.edge_mapper)
        self.edge_actor.GetProperty().SetColor(0, 0, 0)
        self.edge_mapper.SetLookupTable(self.color_function)
        self.edge_mapper.SetResolveCoincidentTopologyToPolygonOffset()

        prop = self.edge_actor.GetProperty()
        prop.SetColor(0, 0, 0)
        self.edge_actor.SetVisibility(self.is_edges)
        self.rend.AddActor(self.edge_actor)

    def post_group_by_name(self, name):
        group = self.groups[name]
        self.post_group(group)
        self.group_active = name

    def post_group(self, group):
        """posts a group object"""
        eids = group.element_ids
        self.show_eids(eids)

    def get_all_eids(self):
        """get the list of all the element IDs"""
        name, result = self.get_name_result_data(0)
        if name != 'ElementID':
            name, result = self.get_name_result_data(1)
            assert name == 'ElementID', name
        return result

    def show_eids(self, eids):
        """shows the specified element IDs"""
        all_eids = self.get_all_eids()

        # remove eids that are out of range
        eids = np.intersect1d(all_eids, eids)

        # update for indices
        i = np.searchsorted(all_eids, eids)

        #eids_off = np.setdiff1d(all_eids, eids)
        #j = np.setdiff1d(all_eids, eids_off)

        self.show_ids_mask(i)

    def hide_eids(self, eids):
        """hides the specified element IDs"""
        all_eids = self.get_all_eids()

        # A-B
        eids = np.setdiff1d(all_eids, eids)

        # update for indices
        i = np.searchsorted(all_eids, eids)
        self.show_ids_mask(i)

    def create_groups_by_property_id(self):
        self._create_groups_by_name('PropertyID', 'property')
        self.log_command('create_groups_by_property_id()')

    def _create_groups_by_name(self, name, prefix):
        #eids = self.find_result_by_name('ElementID')
        #elements_pound = eids.max()
        eids = self.groups['main'].element_ids
        elements_pound = self.groups['main'].elements_pound

        result = self.find_result_by_name(name)
        ures = np.unique(result)
        for uresi in ures:
            ids = np.where(uresi == result)[0]

            name = '%s %s' % (prefix, uresi)
            element_str = ''
            group = Group(
                name, element_str, elements_pound,
                editable=True)
            group.element_ids = eids[ids]
            self.log_info('creating group=%r' % name)
            self.groups[name] = group

    def create_group_with_name(self, name, eids):
        elements_pound = self.groups['main'].elements_pound
        element_str = ''
        group = Group(
            name, element_str, elements_pound,
            editable=True)

        # TODO: make sure all the eids exist
        group.element_ids = eids
        self.log_command('create_group_with_name(%r, %r)' % (name, eids))
        self.groups[name] = group

    def find_result_by_name(self, desired_name):
        for icase in range(self.ncases):
            name, result = self.get_name_result_data(icase)
            if name == desired_name:
                return result
        raise RuntimeError('cannot find name=%r' % desired_name)

    def show_ids_mask(self, ids_to_show):
        flip_flag = True == self._show_flag
        self._update_ids_mask(ids_to_show, flip_flag, show_flag=True, render=False)
        self._update_ids_mask(ids_to_show, False, show_flag=True, render=True)
        self._show_flag = True

    def hide_ids_mask(self, ids_to_hide):
        flip_flag = False == self._show_flag
        self._update_ids_mask(ids_to_hide, flip_flag, show_flag=False, render=False)
        self._update_ids_mask(ids_to_hide, False, show_flag=False, render=True)
        self._show_flag = False

    def _update_ids_mask(self, ids_to_show, flip_flag=True, show_flag=True, render=True):
        ids = vtk.vtkIdTypeArray()
        ids.SetNumberOfComponents(1)
        #ids.SetNumberOfValues(len(ids_to_show))
        ids.Allocate(len(ids_to_show))
        for idi in ids_to_show:
            ids.InsertNextValue(idi)
        ids.Modified()

        if flip_flag:
            self.selection.RemoveAllNodes()
            self.selection_node = vtk.vtkSelectionNode()
            self.selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
            self.selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
            self.selection_node.SetSelectionList(ids)

            if not show_flag:
                self.selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
            self.selection.AddNode(self.selection_node)
        else:
            self.selection_node.SetSelectionList(ids)

        #self.grid_selected.Update() # not in vtk 6

        #ids.Update()
        #self.shown_ids.Modified()
        self.grid_selected.ShallowCopy(self.extract_selection.GetOutput())
        if 0:
            self.selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
            self.extract_selection.Update()

            #self.grid_not_selected = vtk.vtkUnstructuredGrid()
            self.grid_not_selected.ShallowCopy(self.extract_selection.GetOutput())
        self.update_all(render=render)

    def update_all(self, render=True):
        self.grid_selected.Modified()

        #selection_node.Update()
        self.selection_node.Modified()
        #selection.Update()
        self.selection.Modified()
        self.extract_selection.Update()
        self.extract_selection.Modified()

        #grid_selected.Update()
        self.grid_selected.Modified()
        #self.grid_not_selected.Modified()
        self.grid_mapper.Update()
        self.grid_mapper.Modified()
        #selected_actor.Update()
        #selected_actor.Modified()

        #right_renderer.Modified()
        #right_renderer.Update()

        self.iren.Modified()
        #interactor.Update()
        #-----------------
        self.rend.Render()
        #interactor.Start()

        self.rend.Modified()

        self.geom_actor.Modified()
        self.not_selected_actor.Modified()

        if render:
            self.vtk_interactor.Render()
            render_window = self.vtk_interactor.GetRenderWindow()
            render_window.Render()


    def _setup_element_mask(self):
        """
        starts the masking

        self.grid feeds in the geometry
        """
        ids = vtk.vtkIdTypeArray()
        ids.SetNumberOfComponents(1)

        self.selection_node = vtk.vtkSelectionNode()
        self.selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
        self.selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
        self.selection_node.SetSelectionList(ids)

        self.selection = vtk.vtkSelection()
        self.selection.AddNode(self.selection_node)

        self.extract_selection = vtk.vtkExtractSelection()
        if vtk.VTK_MAJOR_VERSION <= 5:
            self.extract_selection.SetInput(0, self.grid)
            self.extract_selection.SetInput(1, self.selection)
        else:
            self.extract_selection.SetInputData(0, self.grid)
            self.extract_selection.SetInputData(1, self.selection)
        self.extract_selection.Update()

        # In selection
        self.grid_selected = vtk.vtkUnstructuredGrid()
        self.grid_selected.ShallowCopy(self.extract_selection.GetOutput())

        #if 0:
        self.selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
        self.extract_selection.Update()

        self.grid_not_selected = vtk.vtkUnstructuredGrid()
        self.grid_not_selected.ShallowCopy(self.extract_selection.GetOutput())

        self.not_selected_mapper = vtk.vtkDataSetMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            self.not_selected_mapper.SetInputConnection(self.grid_not_selected.GetProducerPort())
        else:
            self.not_selected_mapper.SetInputData(self.grid_not_selected)

        self.not_selected_actor = vtk.vtkLODActor()
        self.not_selected_actor.DragableOff()
        self.not_selected_actor.PickableOff()
        self.not_selected_actor.SetMapper(self.not_selected_mapper)


    def create_text(self, position, label, text_size=18):
        text_actor = vtk.vtkTextActor()
        text_actor.SetInput(label)
        text_prop = text_actor.GetTextProperty()
        #text_prop.SetFontFamilyToArial()
        text_prop.SetFontSize(int(text_size))
        text_prop.SetColor(self.text_color)
        text_actor.SetDisplayPosition(*position)

        text_actor.VisibilityOff()

        #txt.SetDisplayPosition(5, 5) # bottom left
        #txt.SetDisplayPosition(5, 95)
        #txt.SetPosition(0.1, 0.5)

        # assign actor to the renderer
        self.rend.AddActor(text_actor)
        self.text_actors[self.itext] = text_actor
        self.itext += 1

    def turn_text_off(self):
        for text in itervalues(self.text_actors):
            text.VisibilityOff()

    def turn_text_on(self):
        for text in itervalues(self.text_actors):
            text.VisibilityOn()

    def build_lookup_table(self):
        """TODO: add support for NanColors"""
        scalar_range = self.grid_selected.GetScalarRange()
        #print('min = %s\nmax = %s' % scalar_range)
        self.grid_mapper.SetScalarRange(scalar_range)
        self.grid_mapper.SetLookupTable(self.color_function)
        self.rend.AddActor(self.scalarBar)

    def _create_load_file_dialog(self, qt_wildcard, title, default_filename=None):
        if default_filename is None:
            default_filename = self.last_dir
        # getOpenFileName return QString and we want Python string
        fname, wildcard_level = QtGui.QFileDialog.getOpenFileNameAndFilter(
            self, title, default_filename, qt_wildcard)
        return str(wildcard_level), str(fname)

    #def _create_load_file_dialog2(self, qt_wildcard, title):
        ## getOpenFileName return QString and we want Python string
        ##title = 'Load a Tecplot Geometry/Results File'
        #last_dir = ''
        ##qt_wildcard = ['Tecplot Hex Binary (*.tec; *.dat)']
        #dialog = MultiFileDialog()
        #dialog.setWindowTitle(title)
        #dialog.setDirectory(self.last_dir)
        #dialog.setFilters(qt_wildcard.split(';;'))
        #if dialog.exec_() == QtGui.QDialog.Accepted:
            #outfiles = dialog.selectedFiles()
            #wildcard_level = dialog.selectedFilter()
            #return str(wildcard_level), str(fname)
        #return None, None

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

    def on_load_geometry(self, infile_name=None, geometry_format=None, name='main', plot=True):
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
            #msg += '%r is not a enabled format; enabled_formats=%s\n' % (
                #geometry_format, self.supported_formats)
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
                #raise NotImplementedError('on_load_geometry; infile_name=%r format=%r' % (
                    #infile_name, geometry_format))
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
            if self.model_type is not None:
                clear_name = 'clear_' + self.model_type
                try:
                    dy_method = getattr(self, clear_name)  # 'self.clear_nastran()'
                    dy_method()
                except:
                    print("method %r does not exist" % clear_name)
            self.log_info("reading %s file %r" % (geometry_format, infile_name))

            # inspect the load_geometry method to see what version it's using
            #args, varargs, keywords, defaults = inspect.getargspec(load_function)
            try:
                #if args[-1] == 'plot':
                has_results = load_function(infile_name, self.last_dir, name=name, plot=plot)
                #else:
                    #name = load_function.__name__
                    #self.log_error(str(args))
                    #self.log_error("'plot' needs to be added to %r; "
                                   #"args[-1]=%r" % (name, args[-1]))
                    #has_results = load_function(infile_name, self.last_dir)
                    #form, cases = load_function(infile_name, self.last_dir)
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
        self.log_command("on_load_geometry(infile_name=%r, geometry_format=%r)" % (
            infile_name, self.format))

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

    def _on_load_nodal_elemental_results(self, result_type, out_filename=None):
        """
        Loads a CSV/TXT results file.  Must have called on_load_geometry first.

        Parameters
        ----------
        out_filename : str / None
            the path to the results file
        """
        # A = np.loadtxt('loadtxt_spike.txt', dtype=('float,int'))
        # dtype=[('f0', '<f8'), ('f1', '<i4')])
        # A['f0']
        # A['f1']
        geometry_format = self.format
        if self.format is None:
            msg = 'on_load_results failed:  You need to load a file first...'
            self.log_error(msg)
            raise RuntimeError(msg)

        if out_filename in [None, False]:
            title = 'Select a %s Results File for %s' % (result_type, self.format)
            wildcard = 'Delimited Text (*.txt; *.dat; *.csv)'
            out_filename = self._create_load_file_dialog(wildcard, title)[1]

        if out_filename == '':
            return
        if not os.path.exists(out_filename):
            msg = 'result file=%r does not exist' % out_filename
            self.log_error(msg)
            return
            #raise IOError(msg)
        # self.last_dir = os.path.split(out_filename)[0]
        try:
            self._load_csv(result_type, out_filename)
        except Exception as e:
            msg = traceback.format_exc()
            self.log_error(msg)
            #return
            raise

        #if 0:
            #self.out_filename = out_filename
            #msg = '%s - %s - %s' % (self.format, self.infile_name, out_filename)
            #self.set_window_title(msg)
            #self.out_filename = out_filename

        if result_type == 'Nodal':
            self.log_command("_on_load_nodal_elemental_results(%r)" % out_filename)
        elif result_type == 'Elemental':
            self.log_command("on_load_elemental_results(%r)" % out_filename)
        else:
            raise NotImplementedError(result_type)

    def _load_csv(self, result_type, out_filename):
        out_filename_short = os.path.basename(out_filename)
        A, fmt_dict, headers = load_csv(out_filename)
        #nrows, ncols, fmts
        header0 = headers[0]
        result0 = A[header0]
        nrows = result0.size

        if result_type == 'Nodal':
            assert nrows == self.nNodes, 'nrows=%s nnodes=%s' % (nrows, self.nNodes)
            result_type2 = 'node'
        elif result_type == 'Elemental':
            assert nrows == self.nElements, 'nrows=%s nelements=%s' % (nrows, self.nElements)
            result_type2 = 'centroid'
        else:
            raise NotImplementedError('result_type=%r' % result_type)

        #print('A =', A)
        formi = []
        form = self.get_form()
        icase = len(self.case_keys)
        islot = 0
        for case_key in self.case_keys:
            if isinstance(case_key, tuple):
                islot = case_key[0]
                break

        for header in headers:
            datai = A[header]
            fmti = fmt_dict[header]
            key = (islot, icase, header, 1, result_type2, fmti, '')
            self.case_keys.append(key)
            self.result_cases[key] = datai
            formi.append((header, icase, []))

            self.label_actors[header] = []
            self.label_ids[header] = set([])
            icase += 1
        form.append((out_filename_short, None, formi))

        self.ncases += len(headers)
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
            out_filename = self._create_load_file_dialog(wildcard, title)[1]
        else:

            for fmt in self.fmts:
                fmt_name, _major_name, _geowild, _geofunc, _reswild, _resfunc = fmt
                print('fmt_name=%r geometry_format=%r' % (fmt_name, geometry_format))
                if fmt_name == geometry_format:
                    load_function = _resfunc
                    break
            else:
                msg = ('format=%r is not supported.  '
                       'Did you load a geometry model?' % geometry_format)
                self.log_error(msg)
                raise RuntimeError(msg)

        if out_filename == '':
            return
        if isinstance(out_filename, string_types):
            out_filename = [out_filename]
        for out_filenamei in out_filename:
            if not os.path.exists(out_filenamei):
                msg = 'result file=%r does not exist' % out_filenamei
                self.log_error(msg)
                return
                #raise IOError(msg)
            self.last_dir = os.path.split(out_filenamei)[0]
            try:
                load_function(out_filenamei, self.last_dir)
            except Exception: #  as e
                msg = traceback.format_exc()
                self.log_error(msg)
                #return
                raise

            self.out_filename = out_filenamei
            msg = '%s - %s - %s' % (self.format, self.infile_name, out_filenamei)
            self.set_window_title(msg)
            print("on_load_results(%r)" % out_filenamei)
            self.out_filename = out_filenamei
            self.log_command("on_load_results(%r)" % out_filenamei)

    def setup_gui(self):
        assert self.fmts != [], 'supported_formats=%s' % self.supported_formats
        self.start_logging()
        settings = QtCore.QSettings()
        self.create_vtk_actors()

        # build GUI and restore saved application state
        #nice_blue = (0.1, 0.2, 0.4)
        white = (1.0, 1.0, 1.0)
        black = (0.0, 0.0, 0.0)
        #red = (1.0, 0.0, 0.0)
        grey = (119/255., 136/255., 153/255.)
        screen_shape_default = (1100, 700)
        qpos_default = self.pos()
        pos_default = qpos_default.x(), qpos_default.y()
        if PY2:
            self.restoreGeometry(settings.value("mainWindowGeometry").toByteArray())

        #self.reset_settings = False
        if self.reset_settings:
            self.background_color = grey
            self.label_color = black
            self.text_color = white
            self.resize(1100, 700)
        else:
            self.background_color = settings.value("backgroundColor", grey).toPyObject()
            self.label_color = settings.value("labelColor", black).toPyObject()
            self.text_color = settings.value("textColor", white).toPyObject()
            screen_shape = settings.value("screen_shape", screen_shape_default).toPyObject()
            #w = screen_shape.width()
            #h = screen_shape.height()
            #try:
            self.resize(screen_shape[0], screen_shape[1])
            width, height = screen_shape
            pos = settings.value("pos", pos_default).toPyObject()
            #x_pos, y_pos = pos
            #print(pos)
            #self.mapToGlobal(QtCore.QPoint(pos[0], pos[1]))
            #self.setGeometry(x_pos, y_pos, width, height)
            #except TypeError:
                #self.resize(1100, 700)


        self.init_ui()
        if self.reset_settings:
            self.res_dock.toggleViewAction()
        self.init_cell_picker()

        if PY2:
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

        self.color_order = [
            (1.0, 0.145098039216, 1.0),
            (0.0823529411765, 0.0823529411765, 1.0),
            (0.0901960784314, 1.0, 0.941176470588),
            (0.501960784314, 1.0, 0.0941176470588),
            (1.0, 1.0, 0.117647058824),
            (1.0, 0.662745098039, 0.113725490196)
        ]
        if inputs['user_points'] is not None:
            for fname in inputs['user_points']:
                self.on_load_user_points(fname)

        if inputs['user_geom'] is not None:
            for fname in inputs['user_geom']:
                self.on_load_user_geom(fname)

    def on_load_user_geom(self, csv_filename=None, name=None, color=None):
        """
        Loads a User Geometry CSV File of the form:

        #    id  x    y    z
        GRID, 1, 0.2, 0.3, 0.3
        GRID, 2, 1.2, 0.3, 0.3
        GRID, 3, 2.2, 0.3, 0.3
        GRID, 4, 5.2, 0.3, 0.3
        grid, 5, 5.2, 1.3, 2.3  # case insensitive

        #    ID, nodes
        BAR,  1, 1, 2
        TRI,  2, 1, 2, 3
        # this is a comment

        QUAD, 3, 1, 5, 3, 4
        QUAD, 4, 1, 2, 3, 4  # this is after a blank line

        #RESULT,4,CENTROID,AREA(%f),PROPERTY_ID(%i)
        # in element id sorted order: value1, value2
        #1.0, 2.0 # bar
        #1.0, 2.0 # tri
        #1.0, 2.0 # quad
        #1.0, 2.0 # quad

        #RESULT,NODE,NODEX(%f),NODEY(%f),NODEZ(%f)
        # same difference

        #RESULT,VECTOR3,GEOM,DXYZ
        # 3xN

        Parameters
        ----------
        csv_filename : str (default=None -> load a dialog)
            the path to the user geometry CSV file
        name : str (default=None -> extract from fname)
            the name for the user points
        color : (float, float, float)
            RGB values as 0.0 <= rgb <= 1.0
        """
        if csv_filename in [None, False]:
            qt_wildcard = '*.csv'
            title = 'Load User Geometry'
            csv_filename = self._create_load_file_dialog(qt_wildcard, title)[1]
            if not csv_filename:
                return

        if color is None:
            # we mod the num_user_points so we don't go outside the range
            icolor = self.num_user_points % len(self.color_order)
            color = self.color_order[icolor]
        if name is None:
            name = os.path.basename(csv_filename).rsplit('.', 1)[0]

        self._add_user_geometry(csv_filename, name, color)
        self.log_command('on_load_user_geom(%r, %r, %s)' % (
            csv_filename, name, str(color)))

    def _add_user_geometry(self, csv_filename, name, color):
        if name in self.geometry_actors:
            msg = 'Name: %s is already in geometry_actors\nChoose a different name.' % name
            raise ValueError(msg)
        if len(name) == 0:
            msg = 'Invalid Name: name=%r' % name
            raise ValueError(msg)

        point_name = name + '_point'
        geom_name = name + '_geom'


        grid_ids, xyz, bars, tris, quads = load_user_geom(csv_filename)
        nbars = len(bars)
        ntris = len(tris)
        nquads = len(quads)
        nelements = nbars + ntris + nquads
        self.create_alternate_vtk_grid(point_name, color=color, opacity=1.0,
                                       point_size=5, representation='point')

        if nelements > 0:
            nid_map = {}
            i = 0
            for nid in grid_ids:
                nid_map[nid] = i
                i += 1
            self.create_alternate_vtk_grid(geom_name, color=color, opacity=1.0,
                                           line_width=5, representation='toggle')

        # allocate
        npoints = len(grid_ids)
        self.alt_grids[point_name].Allocate(npoints, 1000)
        if nelements > 0:
            self.alt_grids[geom_name].Allocate(npoints, 1000)

        # set points
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(npoints)


        if nelements > 0:
            geom_grid = self.alt_grids[geom_name]
            for i, point in enumerate(xyz):
                points.InsertPoint(i, *point)
                elem = vtk.vtkVertex()
                elem.GetPointIds().SetId(0, i)
                self.alt_grids[point_name].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                geom_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
        else:
            for i, point in enumerate(xyz):
                points.InsertPoint(i, *point)
                elem = vtk.vtkVertex()
                elem.GetPointIds().SetId(0, i)
                self.alt_grids[point_name].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
        if nbars:
            for i, bar in enumerate(bars[:, 1:]):
                g1 = nid_map[bar[0]]
                g2 = nid_map[bar[1]]
                elem = vtk.vtkLine()
                elem.GetPointIds().SetId(0, g1)
                elem.GetPointIds().SetId(1, g2)
                geom_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

        if ntris:
            for i, tri in enumerate(tris[:, 1:]):
                g1 = nid_map[tri[0]]
                g2 = nid_map[tri[1]]
                g3 = nid_map[tri[2]]
                elem = vtk.vtkTriangle()
                elem.GetPointIds().SetId(0, g1)
                elem.GetPointIds().SetId(1, g2)
                elem.GetPointIds().SetId(2, g3)
                geom_grid.InsertNextCell(5, elem.GetPointIds())

        if nquads:
            for i, quad in enumerate(quads[:, 1:]):
                g1 = nid_map[quad[0]]
                g2 = nid_map[quad[1]]
                g3 = nid_map[quad[2]]
                g4 = nid_map[quad[3]]
                elem = vtk.vtkQuad()
                point_ids = elem.GetPointIds()
                point_ids.SetId(0, g1)
                point_ids.SetId(1, g2)
                point_ids.SetId(2, g3)
                point_ids.SetId(3, g4)
                geom_grid.InsertNextCell(9, elem.GetPointIds())

        self.alt_grids[point_name].SetPoints(points)
        if nelements > 0:
            self.alt_grids[geom_name].SetPoints(points)

        # create actor/mapper
        self._add_alt_geometry(self.alt_grids[point_name], point_name)
        if nelements > 0:
            self._add_alt_geometry(self.alt_grids[geom_name], geom_name)

        # set representation to points
        #self.geometry_properties[point_name].representation = 'point'
        #self.geometry_properties[geom_name].representation = 'toggle'
        #actor = self.geometry_actors[name]
        #prop = actor.GetProperty()
        #prop.SetRepresentationToPoints()
        #prop.SetPointSize(4)

    def on_load_user_points(self, csv_filename=None, name=None, color=None):
        """
        Loads a User Points CSV File of the form:

        1.0, 2.0, 3.0
        1.5, 2.5, 3.5

        Parameters
        ----------
        csv_filename : str (default=None -> load a dialog)
           the path to the user points CSV file
        name : str (default=None -> extract from fname)
           the name for the user points
        color : (float, float, float)
            RGB values as 0.0 <= rgb <= 1.0

        .. note:: no header line is required
        .. note:: nodes are in the global frame

        .. todo:: support changing the name
        .. todo:: support changing the color
        .. todo:: support overwriting points
        """
        if csv_filename in [None, False]:
            qt_wildcard = '*.csv'
            title = 'Load User Points'
            csv_filename = self._create_load_file_dialog(qt_wildcard, title)[1]
            if not csv_filename:
                return
        if color is None:
            # we mod the num_user_points so we don't go outside the range
            icolor = self.num_user_points % len(self.color_order)
            color = self.color_order[icolor]
        if name is None:
            sline = os.path.basename(csv_filename).rsplit('.', 1)
            name = sline[0]

        self._add_user_points(csv_filename, name, color)
        self.num_user_points += 1
        self.log_command('on_load_user_points(%r, %r, %s)' % (
            csv_filename, name, str(color)))

    def create_cell_picker(self):
        self.cell_picker = vtk.vtkCellPicker()
        self.node_picker = vtk.vtkPointPicker()
        self.cell_picker.SetTolerance(0.0005)

    def mark_node(self, nid, result_name, text):
        raise NotImplementedError()
        #i = self.node_ids.index(nid)
        #x, y, z = self.xyz_cid0[i, :]
        #self._create_annotation(text, result_name, x, y, z)

    def _cell_centroid_pick(self, cell_id, world_position):
        duplicate_key = None
        if self.pick_state == 'node/centroid':
            return_flag = False
            duplicate_key = cell_id
            result_name, result_value, xyz = self.get_result_by_cell_id(cell_id, world_position)
            assert result_name in self.label_actors, result_name
        else:
            #cell = self.grid.GetCell(cell_id)
            # get_nastran_centroidal_pick_state_nodal_by_xyz_cell_id()
            method = 'get_centroidal_%s_result_pick_state_%s_by_xyz_cell_id' % (
                self.format, self.pick_state)
            if hasattr(self, method):
                methodi = getattr(self, method)
                return_flag, value = methodi(world_position, cell_id)
                if return_flag is True:
                    return return_flag, None, None, None, None
            else:
                msg = "pick_state is set to 'nodal', but the result is 'centroidal'\n"
                msg += '  cannot find: self.%s(xyz, cell_id)' % method
                self.log_error(msg)
            return return_flag, None, None, None
        self.log_info("%s = %s" % (result_name, result_value))
        return return_flag, duplicate_key, result_value, result_name, xyz

    def _cell_node_pick(self, cell_id, world_position):
        duplicate_key = None
        if self.pick_state == 'node/centroid':
            return_flag = False
            (result_name, result_value, node_id, xyz) = self.get_result_by_xyz_cell_id(
                world_position, cell_id)
            assert result_name in self.label_actors, result_name
            assert not isinstance(xyz, int), xyz
            duplicate_key = node_id
        else:
            method = 'get_nodal_%s_result_pick_state_%s_by_xyz_cell_id' % (
                self.format, self.pick_state)
            if hasattr(self, method):
                methodi = getattr(self, method)
                return_flag, value = methodi(world_position, cell_id)
                if return_flag is True:
                    return return_flag, None, None, None, None
            else:
                msg = "pick_state is set to 'centroidal', but the result is 'nodal'\n"
                msg += '  cannot find: self.%s(xyz, cell_id)' % method
                self.log_error(msg)
            return return_flag, None, None, None
        msg = "%s = %s" % (result_name, result_value)
        if self.result_name in ['Node_ID', 'Node ID', 'NodeID']:
            msg += '; xyz=(%s, %s, %s)' % tuple(xyz)
        self.log_info(msg)
        return return_flag, duplicate_key, result_value, result_name, xyz

    def init_cell_picker(self):
        self.is_pick = False
        if not self.run_vtk:
            return

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

                #method = 'get_result_by_cell_id()' # self.model_type
                #print('pick_state =', self.pick_state)

                icase = self.icase
                key = self.case_keys[icase]
                location = self.get_case_location(key)

                if location == 'centroid':
                    out = self._cell_centroid_pick(cell_id, world_position)
                elif location == 'node':
                    out = self._cell_node_pick(cell_id, world_position)
                else:
                    raise RuntimeError('invalid pick location=%r' % location)

                return_flag, duplicate_key, result_value, result_name, xyz = out
                if return_flag is True:
                    return

                # prevent duplicate labels with the same value on the same cell
                if duplicate_key is not None and duplicate_key in self.label_ids[result_name]:
                    return
                self.label_ids[result_name].add(duplicate_key)

                if 0:
                    result_value2, xyz2 = self.convert_units(result_name, result_value, xyz)
                    result_value = result_value2
                    xyz2 = xyz
                #x, y, z = world_position
                x, y, z = xyz
                text = '(%.3g, %.3g, %.3g); %s' % (x, y, z, result_value)
                text = str(result_value)
                assert result_name in self.label_actors, result_name
                self._create_annotation(text, result_name, x, y, z)

        def annotate_point_picker(obj, event):
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
        #self.node_picker.AddObserver("EndPickEvent", on_node_picker)p

    def convert_units(self, result_name, result_value, xyz):
        #self.input_units
        #self.display_units
        return result_value, xyz

    def _create_annotation(self, text, result_name, x, y, z):
        if not isinstance(result_name, string_types):
            msg = 'result_name=%r type=%s' % (result_name, type(result_name))
            raise TypeError(msg)
        # http://nullege.com/codes/show/src%40p%40y%40pymatgen-2.9.6%40pymatgen%40vis%40structure_vtk.py/395/vtk.vtkVectorText/python

        #self.convert_units(result_name, result_value, x, y, z)
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
            follower.SetScale(self.dim_max * 0.02 * self.label_scale)

            prop = follower.GetProperty()
            prop.SetColor(self.label_color)
            #prop.SetOpacity( 0.3 );

            # we need to make sure the text rotates when the camera is changed
            camera = self.rend.GetActiveCamera()
            follower.SetCamera(camera)
        else:
            # Create a text mapper and actor to display the results of picking.
            text_mapper = vtk.vtkTextMapper()
            text_mapper.SetInput(text)

            tprop = text_mapper.GetTextProperty()
            tprop.SetFontFamilyToArial()
            tprop.SetFontSize(10)
            tprop.BoldOn()
            tprop.ShadowOn()
            tprop.SetColor(self.label_color)

            text_actor = vtk.vtkActor2D()
            text_actor.GetPositionCoordinate().SetCoordinateSystemToWorld()
            #text_actor.SetPosition(world_position[:2])
            text_actor.SetMapper(text_mapper)
            follower = text_actor

        # finish adding the actor
        self.rend.AddActor(follower)
        self.label_actors[result_name].append(follower)

        #self.picker_textMapper.SetInput("(%.6f, %.6f, %.6f)"% pickPos)
        #camera.GetPosition()
        #camera.GetClippingRange()
        #camera.GetFocalPoint()

    def _on_multi_pick(self, a):
        """
        vtkFrustumExtractor
        vtkAreaPicker
        """
        pass

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

    #def get_2d_point(self, point3d, view_matrix,
                     #projection_matrix,
                     #width, height):
        #view_projection_matrix = projection_matrix * view_matrix
        ## transform world to clipping coordinates
        #point3d = view_projection_matrix.multiply(point3d)
        #win_x = math.round(((point3d.getX() + 1) / 2.0) * width)
        ## we calculate -point3D.getY() because the screen Y axis is
        ## oriented top->down
        #win_y = math.round(((1 - point3d.getY()) / 2.0) * height)
        #return Point2D(win_x, win_y)

    #def get_3d_point(self, point2D, width, height, view_matrix, projection_matrix):
        #x = 2.0 * win_x / client_width - 1
        #y = -2.0 * win_y / client_height + 1
        #view_projection_inverse = inverse(projection_matrix * view_vatrix)
        #point3d = Point3D(x, y, 0)
        #return view_projection_inverse.multiply(point3d)

    def hide_actors(self, except_names=None):
        """
        Hide all the actors

        except_names : str, List[str], None
            list of names to exclude
        """
        if except_names is None:
            except_names = []
        elif isinstance(except_names, string_types):
            except_names = [except_names]

        # hide everything but the main grid
        for key, actor in iteritems(self.geometry_actors):
            if key not in except_names:
                actor.VisibilityOff()
            #else:
                #prop = actor.GetProperty()
                #prop.SetLineWidth(1.5)
                #prop.SetColor((0., 0., 0.))
        self.hide_axes()
        self.hide_legend()
        #self.set_background_color_to_white()

    def hide_axes(self):
        for axis in self.axes.itervalues():
            axis.VisibilityOff()
        self.corner_axis.EnabledOff()

    def show_axes(self):
        for axis in self.axes.itervalues():
            axis.VisibilityOn()
        self.corner_axis.EnabledOn()

    def hide_legend(self):
        self.scalarBar.VisibilityOff()

    def show_legend(self):
        self.scalarBar.VisibilityOn()

    def set_background_color_to_white(self):
        white = (1., 1., 1.)
        self.set_background_color(white)


    def on_take_screenshot(self, fname=None, magnify=None):
        """ Take a screenshot of a current view and save as a file"""
        if fname is None or fname is False:
            filt = QtCore.QString()
            default_filename = ''

            title = ''
            if self.title is not None:
                title = self.title

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

            if magnify is None:
                magnify_min = 1
                magnify = self.magnify if self.magnify > magnify_min else magnify_min
            else:
                magnify = magnify
            if not isinstance(magnify, integer_types):
                msg = 'magnify=%r type=%s' % (magnify, type(magnify))
                raise TypeError(msg)
            self._update_text_size(magnify=magnify)
            render_large.SetMagnification(magnify)

            # multiply linewidth by magnify
            line_widths0 = {}
            point_sizes0 = {}
            for key, geom_actor in iteritems(self.geometry_actors):
                if isinstance(geom_actor, vtk.vtkActor):
                    prop = geom_actor.GetProperty()
                    line_width0 = prop.GetLineWidth()
                    point_size0 = prop.GetPointSize()
                    line_widths0[key] = line_width0
                    point_sizes0[key] = point_size0
                    line_width = line_width0 * magnify
                    point_size = point_size0 * magnify
                    prop.SetLineWidth(line_width)
                    prop.SetPointSize(point_size)
                    prop.Modified()
                elif isinstance(geom_actor, vtk.vtkAxesActor):
                    pass
                else:
                    raise NotImplementedError(geom_actor)

            # hide corner axis
            axes_actor = self.corner_axis.GetOrientationMarker()
            axes_actor.SetVisibility(False)

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
            self.log_command('on_take_screenshot(%r, magnify=%s)' % (fname, magnify))
            self._update_text_size(magnify=1.0)

            # show corner axes
            axes_actor.SetVisibility(True)

            # set linewidth back
            for key, geom_actor in iteritems(self.geometry_actors):
                if isinstance(geom_actor, vtk.vtkActor):
                    prop = geom_actor.GetProperty()
                    prop.SetLineWidth(line_widths0[key])
                    prop.SetPointSize(point_sizes0[key])
                    prop.Modified()
                elif isinstance(geom_actor, vtk.vtkAxesActor):
                    pass
                else:
                    raise NotImplementedError(geom_actor)

    def _update_text_size(self, magnify=1.0):
        """Internal method for updating the bottom-left text when we go to take a picture"""
        text_size = int(14 * magnify)
        for text_actor in itervalues(self.text_actors):
            text_prop = text_actor.GetTextProperty()
            text_prop.SetFontSize(text_size)
        self.itext += 1  # TODO: why is this here?

    def add_geometry(self):
        """
        #(N,)  for stress, x-disp
        #(N,3) for warp vectors/glyphs
        grid_result = vtk.vtkFloatArray()

        point_data = self.grid.GetPointData()
        cell_data = self.grid.GetCellData()

        self.grid.GetCellData().SetScalars(grid_result)
        self.grid.GetPointData().SetScalars(grid_result)


        self.grid_mapper   <-input-> self.grid
        vtkDataSetMapper() <-input-> vtkUnstructuredGrid()

        self.grid_mapper   <--map--> self.geom_actor <-add-> self.rend
        vtkDataSetMapper() <--map--> vtkActor()      <-add-> vtkRenderer()
        """
        if self.is_groups:
            # solid_bending: eids 1-182
            self._setup_element_mask()
            #eids = np.arange(172)
            #eids = arange(172)
            #self.update_element_mask(eids)
        else:
            self.grid_selected = self.grid
        #print('grid_selected =', self.grid_selected)

        self.grid_mapper = vtk.vtkDataSetMapper()
        if self.vtk_version[0] <= 5:
            #self.grid_mapper.SetInput(self.grid_selected)  ## OLD
            self.grid_mapper.SetInputConnection(self.grid_selected.GetProducerPort())
        else:
            self.grid_mapper.SetInputData(self.grid_selected)


        #if 0:
            #self.warp_filter = vtk.vtkWarpVector()
            #self.warp_filter.SetScaleFactor(50.0)
            #self.warp_filter.SetInput(self.grid_mapper.GetUnstructuredGridOutput())

            #self.geom_filter = vtk.vtkGeometryFilter()
            #self.geom_filter.SetInput(self.warp_filter.GetUnstructuredGridOutput())

            #self.geom_mapper = vtk.vtkPolyDataMapper()
            #self.geom_actor.setMapper(self.geom_mapper)

        #if 0:
            ##from vtk.numpy_interface import algorithms
            #arrow = vtk.vtkArrowSource()
            #arrow.PickableOff()

            #self.glyph_transform = vtk.vtkTransform()
            #self.glyph_transform_filter = vtk.vtkTransformPolyDataFilter()
            #self.glyph_transform_filter.SetInputConnection(arrow.GetOutputPort())
            #self.glyph_transform_filter.SetTransform(self.glyph_transform)

            #self.glyph = vtk.vtkGlyph3D()
            ##self.glyph.setInput(xxx)
            #self.glyph.SetSource(self.glyph_transform_filter.GetOutput())

            #self.glyph.SetVectorModeToUseVector()
            #self.glyph.SetColorModeToColorByVector()
            #self.glyph.SetScaleModeToScaleByVector()
            #self.glyph.SetScaleFactor(1.0)

            #self.append_filter = vtk.vtkAppendFilter()
            #self.append_filter.AddInputConnection(self.grid.GetOutput())


        #self.warpVector = vtk.vtkWarpVector()
        #self.warpVector.SetInput(self.grid_mapper.GetUnstructuredGridOutput())
        #grid_mapper.SetInput(Filter.GetOutput())

        self.geom_actor = vtk.vtkLODActor()
        self.geom_actor.DragableOff()
        self.geom_actor.SetMapper(self.grid_mapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        #self.geom_actor.GetProperty().SetDiffuseColor(1, 0, 0)  # red

        #if 0:
            #id_filter = vtk.vtkIdFilter()

            #ids = np.array([1, 2, 3], dtype='int32')
            #id_array = numpy_to_vtk(
                #num_array=ids,
                #deep=True,
                #array_type=vtk.VTK_INT,
            #)

            #id_filter.SetCellIds(id_array.GetOutputPort())
            #id_filter.CellIdsOff()
            #self.grid_mapper.SetInputConnection(id_filter.GetOutputPort())
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

    def _add_alt_geometry(self, grid, name, color=None, line_width=None,
                          opacity=None, representation=None):
        """
        NOTE: color, line_width, opacity are ignored if name already exists
        """
        quad_mapper = vtk.vtkDataSetMapper()
        if name in self.geometry_actors:
            alt_geometry_actor = self.geometry_actors[name]
            if self.vtk_version[0] >= 6:
                alt_geometry_actor.GetMapper().SetInputData(grid)
            else:
                alt_geometry_actor.GetMapper().SetInput(grid)
        else:
            if self.vtk_version[0] >= 6:
                quad_mapper.SetInputData(grid)
            else:
                quad_mapper.SetInput(grid)
            alt_geometry_actor = vtk.vtkActor()
            alt_geometry_actor.SetMapper(quad_mapper)
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
        #prop.SetInterpolationToFlat()    # 0
        #prop.SetInterpolationToGouraud() # 1
        #prop.SetInterpolationToPhong()   # 2
        prop.SetDiffuseColor(color)
        prop.SetOpacity(opacity)
        #prop.Update()

        #print('prop.GetInterpolation()', prop.GetInterpolation()) # 1

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
        self.title = str(title)
        self.min_value = float(min_value)
        self.max_value = float(max_value)

        try:
            data_format % 1
        except:
            msg = ("failed applying the data formatter format=%r and "
                   "should be of the form: '%i', '%8f', '%.2f', '%e', etc.")
            self.log_error(msg)
            return
        self.data_format = data_format
        self.log_command('on_update_scalar_bar(%r, %r, %r, %r)' % (
            title, min_value, max_value, data_format))

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
        print("key_key_press = ", key)
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
            self.case_keys = list(cases.keys())
        else:
            self.case_keys = sorted(cases.keys())
            assert isinstance(cases, dict), type(cases)

        self.result_cases = cases

        if len(self.case_keys) > 1:
            self.icase = -1
            self.ncases = len(self.result_cases)  # number of keys in dictionary
        elif len(self.case_keys) == 1:
            self.icase = -1
            self.ncases = 1
        else:
            self.icase = -1
            self.ncases = 0
        self.set_form(form)

    def _finish_results_io2(self, form, cases):
        self._set_results(form, cases)
        # assert len(cases) > 0, cases
        # if isinstance(cases, OrderedDict):
            # self.case_keys = cases.keys()
        # else:
            # self.case_keys = sorted(cases.keys())
            # assert isinstance(cases, dict), type(cases)

        self.on_update_geometry_properties(self.geometry_properties, write_log=False)
        # self.result_cases = cases

        #print("cases =", cases)
        #print("case_keys =", self.case_keys)

        self.reset_labels()
        self.cycle_results_explicit()  # start at nCase=0
        if self.ncases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

        #data = [
        #    ('A', []),
        #    ('B', []),
        #    ('C', []),
        #]

        #self.case_keys = [
        #    (1, 'ElementID', 1, 'centroid', '%.0f'), (1, 'Region', 1, 'centroid', '%.0f')
        #]
        data = []
        for key in self.case_keys:
            #print(key)
            if isinstance(key, int):
                obj, (i, name) = self.result_cases[key]
                t = (i, [])
            else:
                t = (key[1], [])
            data.append(t)

        self.res_widget.update_results(form)

        key = self.case_keys[0]
        location = self.get_case_location(key)
        method = 'centroid' if location else 'nodal'

        data2 = [(method, None, [])]
        self.res_widget.update_methods(data2)

        if self.is_groups:
            #eids = np.arange(172)
            #eids = []
            #self.hide_elements_mask(eids)
            elements_pound = self.element_ids[-1]
            main_group = Group(
                'main', '', elements_pound,
                editable=False)
            anti_main_group = Group(
                'anti-main', '', elements_pound,
                editable=False)
            main_group.element_ids = self.element_ids
            self.groups['main'] = main_group
            self.groups['anti_main'] = anti_main_group
            self.post_group(main_group)
            #self.show_elements_mask(np.arange(self.nElements))

    def get_result_by_cell_id(self, cell_id, world_position):
        """TODO: should handle multiple cell_ids"""
        case_key = self.case_keys[self.icase] # int for object
        result_name = self.result_name
        case = self.result_cases[case_key]

        if isinstance(case_key, integer_types):
            (obj, (i, res_name)) = case
            subcase_id = obj.subcase_id
            case = obj.get_result(i, res_name)
        try:
            result_values = case[cell_id]
        except IndexError:
            msg = ('case[cell_id] is out of bounds; length=%s\n'
                   'result_name=%r cell_id=%r case_key=%r\n' % (
                       len(case), result_name, cell_id, case_key))
            raise IndexError(msg)

        cell = self.grid_selected.GetCell(cell_id)
        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()
        cell_type = cell.GetCellType()

        if cell_type in [5, 9, 22, 23]:  # CTRIA3, CQUAD4, CTRIA6, CQUAD8
            node_xyz = np.zeros((nnodes, 3), dtype='float32')
            for ipoint in range(nnodes):
                point = points.GetPoint(ipoint)
                node_xyz[ipoint, :] = point
            xyz = node_xyz.mean(axis=0)
        elif cell_type in [10, 12, 13]: # CTETRA4, CHEXA8, CPENTA6
            # TODO: No idea how to get the center of the face
            #       vs. a point on a face that's not exposed
            #faces = cell.GetFaces()
            #nfaces = cell.GetNumberOfFaces()
            #for iface in range(nfaces):
                #face = cell.GetFace(iface)
                #points = face.GetPoints()
            #faces
            xyz = world_position
        elif cell_type in [24]: # CTETRA10
            xyz = world_position
        elif cell_type in [3]: # CBAR
            xyz = world_position
        else:
            #self.log.error(msg)
            msg = 'cell_type=%s nnodes=%s; result_name=%s result_values=%s' % (
                cell_type, nnodes, result_name, result_values)
            self.log.error(msg)
            raise NotImplementedError(msg)
        return result_name, result_values, xyz

    def get_result_by_xyz_cell_id(self, node_xyz, cell_id):
        """won't handle multiple cell_ids/node_xyz"""
        case_key = self.case_keys[self.icase]
        result_name = self.result_name

        cell = self.grid_selected.GetCell(cell_id)
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
        xyz = np.array(point_min, dtype='float32')
        case = self.result_cases[case_key]
        if isinstance(case_key, integer_types):
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

          #if len(key) == 5:
              #(subcase_id, result_type, vector_size, location, data_format) = key
          #elif len(key) == 6:
              #(subcase_id, j, result_type, vector_size, location, data_format) = key
          else:
              (subcase_id, j, result_type, vector_size, location, data_format, label2) = key
        """
        # case_key = (1, 'ElementID', 1, 'centroid', '%.0f')
        case_key = self.case_keys[self.icase]
        if isinstance(case_key, int):
            obj, (i, name) = self.result_cases[case_key]
            value = name
        else:
            assert len(case_key) == 7, case_key
            #if len(case_key) == 5:
                #value = case_key[1]
            #else:
            value = case_key[2]
        return value

    def finish_io(self, cases):
        self.result_cases = cases
        self.case_keys = sorted(cases.keys())
        #print("case_keys = ", self.case_keys)

        if len(self.result_cases) == 0:
            self.ncases = 1
            self.icase = 0
        elif len(self.result_cases) == 1:
            self.ncases = 1
            self.icase = 0
        else:
            self.ncases = len(self.result_cases) - 1  # number of keys in dictionary
            self.icase = -1
        self.cycle_results()  # start at nCase=0

        if self.ncases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

    def _finish_results_io(self, cases):
        self.result_cases = cases
        self.case_keys = sorted(cases.keys())

        if len(self.case_keys) > 1:
            self.icase = -1
            self.ncases = len(self.result_cases)  # number of keys in dictionary
        elif len(self.case_keys) == 1:
            self.icase = -1
            self.ncases = 1
        else:
            self.icase = -1
            self.ncases = 0

        self.reset_labels()
        self.cycle_results_explicit()  # start at nCase=0
        if self.ncases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

        #data = [
        #    ('A',[]),
        #    ('B',[]),
        #    ('C',[]),
        #]

        #self.case_keys = [
        #    (1, 'ElementID', 1, 'centroid', '%.0f'), (1, 'Region', 1, 'centroid', '%.0f')
        #]
        data = []
        for i, key in enumerate(self.case_keys):
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

        #self.case_keys = [
            #(1, 'ElementID', 1, 'centroid', '%.0f'),
            #(1, 'Region', 1, 'centroid', '%.0f')
        #]
        for case_key in self.case_keys:
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
        #case_key = self.case_keys[self.icase]
        result_name = self.result_name

        actors = self.label_actors[result_name]
        for actor in actors:
            self.rend.RemoveActor(actor)
            del actor
        self.label_actors[result_name] = []
        self.label_ids[result_name] = set([])

    def resize_labels(self, result_names=None, show_msg=True):
        """
        This resizes labels for all result cases.
        TODO: not done...
        """
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
                #prop = actor.GetProperty()
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
            try:
                actors = self.label_actors[key]
            except KeyError:
                msg = 'Cant find label_actors; keys=%s' % self.label_actors.keys()
                self.log.error(msg)
                continue
            for actor in actors:
                actor.VisibilityOn()
                count += 1
        if count and show_msg:
            # yes the ) is intentionally left off because it's already been added
            self.log_command('show_labels(%s' % names)

    def update_scalar_bar(self, title, min_value, max_value, norm_value,
                          data_format,
                          nlabels=None, labelsize=None,
                          ncolors=None, colormap='jet',
                          is_low_to_high=True, is_horizontal=True,
                          is_shown=True):
        """
        Updates the Scalar Bar

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

        nlabels : int (default=None -> auto)
            the number of labels
        labelsize : int (default=None -> auto)
            the label size
        ncolors : int (default=None -> auto)
            the number of colors
        colormap : varies
            str :
                the name
            ndarray : (N, 3) float ndarry
                red-green-blue array

        is_low_to_high : bool; default=True
            flips the order of the RGB points
        is_horizontal : bool; default=True
            makes the scalar bar horizontal
        is_shown : bool
            show the scalar bar
        """
        #print("update_scalar_bar min=%s max=%s norm=%s" % (min_value, max_value, norm_value))
        self.scalar_bar.update(title, min_value, max_value, norm_value, data_format,
                               nlabels=nlabels, labelsize=labelsize,
                               ncolors=ncolors, colormap=colormap,
                               is_low_to_high=is_low_to_high, is_horizontal=is_horizontal,
                               is_shown=is_shown)

    #---------------------------------------------------------------------------------------
    # CAMERA MENU
    def view_camera(self):
        camera = self.rend.GetActiveCamera()
        #position = camera.GetPosition()
        #clip_range = camera.GetClippingRange()
        #focal_point = camera.GetFocalPoint()

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
        camera_data = [
            position, focal_point, view_angle, view_up, clip_range,
            parallel_scale, parallel_proj, distance
        ]
        return camera_data

    def on_set_camera_data(self, camera_data, show_log=True):
        """
        position : (float, float, float)
            where am I is xyz space
        focal_point : (float, float, float)
            where am I looking
        view_angle : float
            field of view (angle); perspective only?
        view_up : (float, float, float)
            up on the screen vector
        clip_range : (float, float)
            start/end distance from camera where clipping starts
        parallel_scale : float
            ???
        parallel_projection : bool (0/1)
            flag?
        distance : float
            distance to ???

        i_vector = focal_point - position
        j'_vector = view_up

        use:
           i x j' -> k
           k x i -> j
           or it's like k'
        """
        #position, clip_range, focal_point, view_up, distance = camera_data
        (position, focal_point, view_angle, view_up, clip_range,
         parallel_scale, parallel_proj, distance) = camera_data

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
            self.log_command(
                'on_set_camera_data([%s, %s, %s, %s, %s, %s, %s, %s])'
                % (position, focal_point, view_angle, view_up,
                   clip_range, parallel_scale, parallel_proj, distance))

    #---------------------------------------------------------------------------------------
    # LABEL SIZE/COLOR
    def on_set_labelsize_color(self):
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
        if not hasattr(self, 'case_keys'):
            self.log_error('No model has been loaded.')
            return

        data = {
            'size' : self.label_text_size,
            'color' : self.label_color,
            'dim_max' : self.dim_max,
            #'clicked_ok' : False,
            #'clicked_cancel' : False,
            #'close' : False,
        }
        #print(data)
        if not self._label_window_shown:
            self._label_window = ModifyLabelPropertiesMenu(data, win_parent=self)
            self._label_window.show()
            self._label_window_shown = True
            self._label_window.exec_()
        else:
            self._label_window.activateWindow()

        if 'close' not in data:
            self._label_window.activateWindow()
            return

        if data['close']:
            self._label_window_shown = False
            del self._label_window
        else:
            self._label_window.activateWindow()

    def set_labelsize_color(self, size=None, color=None):
        """
        Parameters
        ----------
        size : float
            label size
        color : (float, float, float)
            RGB values
        """
        if size is not None:
            assert isinstance(size, (int, float)), 'size=%r' % size
            self.set_labelsize(size)
        if color is not None:
            assert len(color) == 3, color
            assert isinstance(color[0], float), 'color=%r' % color
            self.set_label_color(color)

    @property
    def label_text_size(self):
        return self.dim_max * 0.02 * self.label_scale

    @label_text_size.setter
    def label_text_size(self, label_text_size):
        #self.label_text_size = self.dim_max * 0.02 * self.label_scale
        #a = b * c * d
        #d = a / bc
        self.label_scale = label_text_size / (self.dim_max * 0.02)

    def set_labelsize(self, size, render=True):
        """Updates the size of all the labels"""
        assert size >= 0., size
        self.label_text_size = size
        for result_name, follower_actors in iteritems(self.label_actors):
            for follower_actor in follower_actors:
                follower_actor.SetScale(size)
                follower_actor.Modified()
        if render:
            self.vtk_interactor.GetRenderWindow().Render()
            self.log_command('set_labelsize(%s)' % size)


    def set_label_color(self, color, render=True):
        """
        Set the label color

        Parameters
        ----------
        color : (float, float, float)
            RGB values as floats
        """
        if np.allclose(self.label_color, color):
            return
        self.label_color = color
        for follower_actors in itervalues(self.label_actors):
            for follower_actor in follower_actors:
                prop = follower_actor.GetProperty()
                prop.SetColor(*color)

        if render:
            self.vtk_interactor.GetRenderWindow().Render()
            self.log_command('set_label_color(%s, %s, %s)' % color)

    #---------------------------------------------------------------------------------------
    # PICKER MENU
    def on_set_picker_size(self):
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
        #if not hasattr(self, 'case_keys'):
            #self.log_error('No model has been loaded.')
            #return

        #print('size =', self.element_picker_size)
        size = 10.
        size = self.get_element_picker_size()
        data = {
            'size' : size,
            'dim_max' : self.dim_max,
            #'clicked_ok' : False,
            #'clicked_cancel' : False,
            #'close' : False,
        }
        #print(data)
        if not self._picker_window_shown:
            self._picker_window = ModifyPickerPropertiesMenu(data, win_parent=self)
            self._picker_window.show()
            self._picker_window_shown = True
            self._picker_window.exec_()
        else:
            self._picker_window.activateWindow()

        if 'close' not in data:
            self._picker_window.activateWindow()
            return

        if data['close']:
            self._picker_window_shown = False
            del self._picker_window
        else:
            self._picker_window.activateWindow()

    def get_element_picker_size(self):
        return self.cell_picker.GetTolerance()

    @property
    def element_picker_size(self):
        return self.get_element_picker_size()

    @element_picker_size.setter
    def element_picker_size(self, size):
        """sets the element picker size"""
        self.cell_picker.SetTolerance(size)

    def set_element_picker_size(self, size):
        """Updates the element picker size"""

        assert size >= 0., size
        self.element_picker_size = size


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
        #if not hasattr(self, 'case_keys'):  # TODO: maybe include...
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
        if not self._clipping_window_shown:
            self._clipping_window = ClippingPropertiesWindow(data, win_parent=self)
            self._clipping_window.show()
            self._clipping_window_shown = True
            self._clipping_window.exec_()
        else:
            self._clipping_window.activateWindow()

        if data['close']:
            self._apply_clipping(data)
            del self._clipping_window
            self._clipping_window_shown = False
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
        if not hasattr(self, 'case_keys') or len(self.case_keys) == 0:
            self.log_error('No model has been loaded.')
            return
        key = self.case_keys[self.icase]
        case = self.result_cases[key]
        default_format = None
        default_scale = None
        if isinstance(key, integer_types):
            #(subcase_id, result_type, vector_size, location, data_format) = key
            (obj, (i, res_name)) = self.result_cases[key]
            #subcase_id = obj.subcase_id
            case = obj.get_result(i, res_name)
            result_type = obj.get_title(i, res_name)
            nlabels, labelsize, ncolors, colormap = obj.get_nlabels_labelsize_ncolors_colormap(i, res_name)

            defaults_scalar_bar = obj.get_default_nlabels_labelsize_ncolors_colormap(i, res_name)
            default_nlabels, default_labelsize, default_ncolors, default_colormap = defaults_scalar_bar

            #vector_size = obj.get_vector_size(i, res_name)
            #location = obj.get_location(i, res_name)
            data_format = obj.get_data_format(i, res_name)
            scale = obj.get_scale(i, res_name)

            default_title = obj.get_default_title(i, res_name)
            default_scale = obj.get_default_scale(i, res_name)

            min_value, max_value = obj.get_min_max(i, res_name)
            default_min, default_max = obj.get_default_min_max(i, res_name)
        #elif len(key) == 5:
            #(subcase_id, result_type, vector_size, location, data_format) = key
            #default_title = result_type
            #scale = 0.0
        #elif len(key) == 6:
            #(subcase_id, i, result_type, vector_size, location, data_format) = key
            #default_title = result_type
            #scale = 0.0
        else:
            (subcase_id, i, result_type, vector_size, location, data_format, label2) = key
            default_title = result_type
            scale = 0.0
            default_scale = 0.0
            min_value = case.min()
            max_value = case.max()
            default_min = min_value
            default_max = max_value
            nlabels = None
            labelsize = None
            ncolors = None
            colormap = 'jet'
            default_nlabels = None
            default_labelsize = None
            default_ncolors = None
            default_colormap = 'jet'

        if default_format is None:
            default_format = data_format
        print(key)

        #if isinstance(case, ndarray):
            #if len(case.shape) == 1:
                #normi = case
            #else:
                #normi = norm(case, axis=1)
        #else:
            #raise RuntimeError('list-based results have been removed; use numpy.array')

        data = {
            'icase' : i,
            'name' : result_type,
            #'min' : normi.min(),
            #'max' : normi.max(),
            'min' : min_value,
            'max' : max_value,

            'scale' : scale,
            'format' : data_format,

            'default_min' : default_min,
            'default_max' : default_max,
            'default_title' : default_title,
            'default_scale' : default_scale,
            'default_format' : default_format,

            'default_nlabels' : default_nlabels,
            'default_labelsize' : default_labelsize,
            'default_ncolors' : default_ncolors,
            'default_colormap' : default_colormap,

            'nlabels' : nlabels,
            'labelsize' :  labelsize,
            'ncolors' : ncolors,
            'colormap' : colormap,

            'is_low_to_high' : True,
            'is_discrete': True,
            'is_horizontal': self.scalar_bar.is_horizontal,
            'is_shown' : self.scalar_bar.is_shown,
            'clicked_ok' : False,
            'close' : False,
        }
        print(data)
        if not self._legend_window_shown:
            self._legend_window = LegendPropertiesWindow(data, win_parent=self)
            self._legend_window.show()
            self._legend_window_shown = True
            self._legend_window.exec_()
        else:
            self._legend_window.activateWindow()

        if data['close']:
            if not self._legend_window._updated_legend:
                self._apply_legend(data)
            self._legend_window_shown = False
            del self._legend_window
        else:
            self._legend_window.activateWindow()

    def update_legend(self, icase, name, min_value, max_value, data_format, scale,
                      nlabels, labelsize, ncolors, colormap,
                      is_low_to_high, is_horizontal_scalar_bar):
        if not self._legend_window_shown:
            return
        self._legend_window._updated_legend = True

        key = self.case_keys[icase]
        if isinstance(key, integer_types):
            (obj, (i, name)) = self.result_cases[key]
            #subcase_id = obj.subcase_id
            #case = obj.get_result(i, name)
            #result_type = obj.get_title(i, name)
            #vector_size = obj.get_vector_size(i, name)
            #location = obj.get_location(i, name)
            #data_format = obj.get_data_format(i, name)
            #scale = obj.get_scale(i, name)
            #label2 = obj.get_header(i, name)
            default_data_format = obj.get_default_data_format(i, name)
            default_min, default_max = obj.get_default_min_max(i, name)
            default_scale = obj.get_default_scale(i, name)
            default_title = obj.get_default_title(i, name)
            out_labels = obj.get_default_nlabels_labelsize_ncolors_colormap(i, name)
            default_nlabels, default_labelsize, default_ncolors, default_colormap = out_labels
        else:
            assert len(key) == 7, key
            (subcase_id, j, result_type, vector_size, location, data_format, label2) = key
            case = self.result_cases[key]
            default_data_format = data_format
            default_min = case.min()
            default_max = case.max()
            default_scale = 0.
            default_title = result_type
            default_nlabels = None
            default_labelsize = None
            default_ncolors = None
            default_colormap = 'jet'

        assert isinstance(scale, float), 'scale=%s' % scale
        self._legend_window.update_legend(
            icase,
            name, min_value, max_value, data_format, scale,
            nlabels, labelsize,
            ncolors, colormap,
            default_title, default_min, default_max, default_data_format, default_scale,
            default_nlabels, default_labelsize,
            default_ncolors, default_colormap,
            is_low_to_high, is_horizontal_scalar_bar)
        #self.scalar_bar.set_visibility(self._legend_shown)
        #self.vtk_interactor.Render()

    def _apply_legend(self, data):
        title = data['name']
        min_value = data['min']
        max_value = data['max']
        scale_value = data['scale']
        data_format = data['format']
        is_low_to_high = data['is_low_to_high']
        is_discrete = data['is_discrete']
        is_horizontal = data['is_horizontal']
        is_shown = data['is_shown']

        nlabels = data['nlabels']
        labelsize = data['labelsize']
        ncolors = data['ncolors']
        colormap = data['colormap']

        #print('is_shown1 =', is_shown)
        self.on_update_legend(title=title, min_value=min_value, max_value=max_value,
                              scale=scale_value, data_format=data_format,
                              is_low_to_high=is_low_to_high,
                              is_discrete=is_discrete, is_horizontal=is_horizontal,
                              nlabels=nlabels, labelsize=labelsize,
                              ncolors=ncolors, colormap=colormap,
                              is_shown=is_shown)

    def on_update_legend(self, title='Title', min_value=0., max_value=1., scale=0.0,
                         data_format='%.0f',
                         is_low_to_high=True, is_discrete=True, is_horizontal=True,
                         nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                         is_shown=True):
        #print('is_shown2 =', is_shown)
        #assert is_shown == False, is_shown
        key = self.case_keys[self.icase]
        name_vector = None
        plot_value = self.result_cases[key] # scalar
        vector_size1 = 1
        update_3d = False
        if isinstance(key, integer_types):
            #(subcase_id, result_type, vector_size, location, data_format) = key
            (obj, (i, res_name)) = self.result_cases[key]
            subcase_id = obj.subcase_id
            #print('plot_value =', plot_value)

            result_type = obj.get_title(i, res_name)
            vector_size = obj.get_vector_size(i, res_name)
            if vector_size == 3:
                plot_value = obj.get_plot_value(i, res_name) # vector
                update_3d = True
                #print('setting scale=%s' % scale)
                assert isinstance(scale, float), scale
                obj.set_scale(i, res_name, scale)
            else:
                scalar_result = obj.get_scalar(i, res_name)

            location = obj.get_location(i, res_name)
            obj.set_min_max(i, res_name, min_value, max_value)
            obj.set_data_format(i, res_name, data_format)
            obj.set_nlabels_labelsize_ncolors_colormap(
                i, res_name, nlabels, labelsize, ncolors, colormap)

            #data_format = obj.get_data_format(i, res_name)
            #obj.set_format(i, res_name, data_format)
            #obj.set_data_format(i, res_name, data_format)
            subtitle, label = self.get_subtitle_label(subcase_id)
            name_vector = (vector_size1, subcase_id, result_type, label,
                           min_value, max_value, scale)

        #elif len(key) == 5:
            #(subcase_id, result_type, vector_size1, location, _data_format) = key
        #elif len(key) == 6:
            #(subcase_id, i, result_type, vector_size1, location, _data_format) = key
        else:
            (subcase_id, i, result_type, vector_size1, location, _data_format, label2) = key
            scalar_result = plot_value
        assert vector_size1 == 1, vector_size1

        #if isinstance(key, integer_types):  # vector 3
            #norm_plot_value = norm(plot_value, axis=1)
            #min_value = norm_plot_value.min()
            #max_value = norm_plot_value.max()
            #print('norm_plot_value =', norm_plot_value)

        if update_3d:
            self.is_horizontal_scalar_bar = is_horizontal
            self._set_case(self.result_name, self.icase,
                           explicit=False, cycle=False, skip_click_check=True,
                           min_value=min_value, max_value=max_value,
                           is_legend_shown=is_shown)
            return

        subtitle, label = self.get_subtitle_label(subcase_id)
        scale1 = 0.0
        # if vector_size == 3:

        name = (vector_size1, subcase_id, result_type, label, min_value, max_value, scale1)
        norm_value = float(max_value - min_value)
        # if name not in self._loaded_names:

        #if isinstance(key, integer_types):  # vector 3
            #norm_plot_value = norm(plot_value, axis=1)
            #grid_result = self.set_grid_values(name, norm_plot_value, vector_size1,
                                               #min_value, max_value, norm_value,
                                               #is_low_to_high=is_low_to_high)
        #else:
        grid_result = self.set_grid_values(name, scalar_result, vector_size1,
                                           min_value, max_value, norm_value,
                                           is_low_to_high=is_low_to_high)

        grid_result_vector = None
        #if name_vector and 0:
            #vector_size = 3
            #grid_result_vector = self.set_grid_values(name_vector, plot_value, vector_size,
                                                      #min_value, max_value, norm_value,
                                                      #is_low_to_high=is_low_to_high)

        self.update_scalar_bar(title, min_value, max_value, norm_value,
                               data_format,
                               nlabels=nlabels, labelsize=labelsize,
                               ncolors=ncolors, colormap=colormap,
                               is_low_to_high=is_low_to_high,
                               is_horizontal=is_horizontal, is_shown=is_shown)

        revert_displaced = True
        self._final_grid_update(name, grid_result, None, None, None,
                                1, subcase_id, result_type, location, subtitle, label,
                                revert_displaced=revert_displaced)
        if grid_result_vector is not None:
            self._final_grid_update(name_vector, grid_result_vector, obj, i, res_name,
                                    vector_size, subcase_id, result_type, location, subtitle, label,
                                    revert_displaced=False)
            #if 0:
                #xyz_nominal, vector_data = obj.get_vector_result(i, res_name)
                #self._update_grid(vector_data)
                #self.grid.Modified()
                #self.geom_actor.Modified()
                #self.vtk_interactor.Render()
            #revert_displaced = False
        #self._final_grid_update(name, grid_result, None, None, None,
                                #1, subcase_id, result_type, location, subtitle, label,
                                #revert_displaced=revert_displaced)

        #self.is_horizontal_scalar_bar = is_horizontal
        icase = i
        msg = ('self.on_update_legend(title=%r, min_value=%s, max_value=%s,\n'
               '                      data_format=%r, is_low_to_high=%s, is_discrete=%s,\n'
               '                      nlabels=%r, labelsize=%r, ncolors=%r, colormap=%r,\n'
               '                      is_horizontal=%r, is_shown=%r)'
               % (title, min_value, max_value, data_format, is_low_to_high, is_discrete,
                  nlabels, labelsize, ncolors, colormap, is_horizontal, is_shown))
        self.log_command(msg)
        #if is_shown:
            #pass
    #---------------------------------------------------------------------------------------
    # EDIT ACTOR PROPERTIES
    def edit_geometry_properties(self):
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
        if not hasattr(self, 'case_keys'):
            self.log_error('No model has been loaded.')
            return
        if not len(self.geometry_properties):
            self.log_error('No secondary geometries to edit.')
            return
        #print('geometry_properties.keys() =', self.geometry_properties.keys())
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]
        #if len(key) == 5:
            #(subcase_id, result_type, vector_size, location, data_format) = key
        #elif len(key) == 6:
            #(subcase_id, i, result_type, vector_size, location, data_format) = key
        #else:
            #(subcase_id, i, result_type, vector_size, location, data_format, label2) = key

        data = deepcopy(self.geometry_properties)
        if not self._edit_geometry_properties_window_shown:
            self._edit_geometry_properties = EditGeometryProperties(data, win_parent=self)
            self._edit_geometry_properties.show()
            self._edit_geometry_properties_window_shown = True
            self._edit_geometry_properties.exec_()
        else:
            self._edit_geometry_properties.activateWindow()

        if 'clicked_ok' not in data:
            self._edit_geometry_properties.activateWindow()
            return

        if data['clicked_ok']:
            self.on_update_geometry_properties(data)
            self._save_geometry_properties(data)
            del self._edit_geometry_properties
            self._edit_geometry_properties_window_shown = False
        elif data['clicked_cancel']:
            self.on_update_geometry_properties(self.geometry_properties)
            del self._edit_geometry_properties
            self._edit_geometry_properties_window_shown = False

    def _save_geometry_properties(self, out_data):
        for name, group in iteritems(out_data):
            if name in ['clicked_ok', 'clicked_cancel']:
                continue

            #color2 = group.color_float
            geom_prop = self.geometry_properties[name]
            if isinstance(geom_prop, CoordProperties):
                pass
            elif isinstance(geom_prop, AltGeometry):
                geom_prop.color = group.color
                geom_prop.line_width = group.line_width
                geom_prop.opacity = group.opacity
                geom_prop.point_size = group.point_size
            else:
                raise NotImplementedError(geom_prop)

    def modify_group(self):
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
        if not len(self.groups):  # no 'main' group
            self.log_error('No main group to create.')
            return
        print('groups.keys() =', self.groups.keys())

        data = {0 : self.groups['main']}

        i = 1
        for name, group in sorted(iteritems(self.groups)):
            if name == 'main':
                continue
            data[i] = group
            i += 1
        #data = deepcopy(self.groups)

        if not self._modify_groups_window_shown:
            self._modify_groups = GroupsModify(
                data, win_parent=self, group_active=self.group_active)
            self._modify_groups.show()
            self._modify_groups_window_shown = True
            self._modify_groups.exec_()
        else:
            self._modify_groups.activateWindow()

        if 'clicked_ok' not in data:
            self._modify_groups.activateWindow()
            return

        if data['clicked_ok']:
            self.on_update_modify_group(data)
            imain = self._modify_groups.imain
            name = self._modify_groups.keys[imain]
            self.post_group_by_name(name)
            #name =
            #self._save_geometry_properties(data)
            del self._modify_groups
            self._modify_groups_window_shown = False
        elif data['clicked_cancel']:
            self.on_update_modify_group(data)
            del self._modify_groups
            self._modify_groups_window_shown = False

    def on_update_modify_group(self, out_data):
        """
        Applies the changed groups to the different groups if
        something changed.
        """
        #self.groups = out_data
        data = {}
        for group_id, group in sorted(iteritems(out_data)):
            if not isinstance(group, Group):
                continue
            data[group.name] = group
        self.groups = data

    def on_update_geometry_properties(self, out_data, write_log=True):
        """
        Applies the changed properties to the different actors if
        something changed.

        Note that some of the values are limited.  This prevents
        points/lines from being shrunk to 0 and also the actor
        being actually "hidden" at the same time.
        """
        lines = []
        for name, group in iteritems(out_data):
            if name in ['clicked_ok', 'clicked_cancel']:
                continue
            actor = self.geometry_actors[name]
            if isinstance(actor, vtk.vtkActor):
                lines += self._update_geomtry_properties_actor(name, group, actor)
            elif isinstance(actor, vtk.vtkAxesActor):
                changed = False
                is_visible1 = bool(actor.GetVisibility())
                is_visible2 = group.is_visible
                if is_visible1 != is_visible2:
                    actor.SetVisibility(is_visible2)
                    alt_prop = self.geometry_properties[name]
                    alt_prop.is_visible = is_visible2
                    actor.Modified()
                    changed = True

                if changed:
                    lines.append('    %r : CoordProperties(is_visible=%s),\n' % (
                        name, is_visible2))
            else:
                raise NotImplementedError(actor)

        self.vtk_interactor.Render()
        if write_log and lines:
            msg = 'out_data = {\n'
            msg += ''.join(lines)
            msg += '}\n'
            msg += 'self.on_update_geometry_properties(out_data)'
            self.log_command(msg)

    def _update_geomtry_properties_actor(self, name, group, actor):
        """
        Updates an actor

        Parameters
        ----------
        name : str
            the geometry proprety to update
        group : AltGeometry()
            a storage container for all the actor's properties
        actor : vtkActor()
            the actor where the properties will be applied

        linewidth1 : int
            the active linewidth
        linewidth2 : int
            the new linewidth
        """
        lines = []
        changed = False
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

        representation = group.representation
        alt_prop = self.geometry_properties[name]
        #representation = alt_prop.representation
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
        if representation == 'bar' and bar_scale1 != bar_scale2:
            #print('name=%s rep=%r bar_scale1=%s bar_scale2=%s' % (
                #name, representation, bar_scale1, bar_scale2))
            self.set_bar_scale(name, bar_scale2)
        if is_visible1 != is_visible2:
            actor.SetVisibility(is_visible2)
            alt_prop.is_visible = is_visible2
            #prop.SetViPointSize(is_visible2)
            actor.Modified()
            changed = True

        if changed:
            lines.append('    %r : AltGeometry(self, %r, color=(%s, %s, %s), '
                         'line_width=%s, opacity=%s, point_size=%s, bar_scale=%s, '
                         'representation=%r, is_visible=%s),\n' % (
                             name, name, color2[0], color2[1], color2[2], line_width2,
                             opacity2, point_size2, bar_scale2, representation, is_visible2))
            prop.Modified()
        return lines

    def set_bar_scale(self, name, bar_scale):
        """
        Parameters
        ----------
        name : str
           the parameter to scale (e.g. TUBE_y, TUBE_z)
        bar_scale : float
           the scaling factor
        """
        #print('set_bar_scale - GuiCommon2; name=%s bar_scale=%s' % (name, bar_scale))
        if bar_scale <= 0.0:
            return
        assert bar_scale > 0.0, 'bar_scale=%r' % bar_scale

        # bar_y : (nbars, 6) float ndarray
        #     the xyz coordinates for (node1, node2) of the y/z axis of the bar
        #     xyz1 is the centroid
        #     xyz2 is the end point of the axis with a length_xyz with a bar_scale of 1.0
        bar_y = self.bar_lines[name]

        #dy = c - yaxis
        #dz = c - zaxis
        #print('bary:\n%s' % bar_y)
        xyz1 = bar_y[:, :3]
        xyz2 = bar_y[:, 3:]
        dxyz = xyz2 - xyz1

        # vectorized version of L = sqrt(dx^2 + dy^2 + dz^2)
        length_xyz = np.linalg.norm(dxyz, axis=1)
        izero = np.where(length_xyz == 0.0)[0]
        if len(izero):
            bad_eids = self.bar_eids[name][izero]
            self.log.error('The following elements have zero length...%s' % bad_eids)

        # v = dxyz / length_xyz *  bar_scale
        # xyz2 = xyz1 + v

        nnodes = len(length_xyz)
        grid = self.alt_grids[name]
        points = grid.GetPoints()
        for i in range(nnodes):
            p = points.GetPoint(2*i+1)
            #print(p)
            node = xyz1[i, :] + length_xyz[i] * bar_scale * dxyz[i, :]
            #print(p, node)
            points.SetPoint(2 * i + 1, *node)

        if hasattr(grid, 'Update'):
            #print('update....')
            grid.Update()
        grid.Modified()
        #print('update2...')

    def _add_user_points(self, points_filename, name, color):
        if name in self.geometry_actors:
            msg = 'Name: %s is already in geometry_actors\nChoose a different name.' % name
            raise ValueError(msg)
        if len(name) == 0:
            msg = 'Invalid Name: name=%r' % name
            raise ValueError(msg)

        # create grid
        self.create_alternate_vtk_grid(name, color=color, line_width=5, opacity=1.0,
                                       point_size=1, representation='point')

        assert os.path.exists(points_filename), print_bad_path(points_filename)
        # read input file
        try:
            user_points = np.loadtxt(points_filename, delimiter=',')
        except ValueError:
            user_points = loadtxt_nice(points_filename, delimiter=',')
            # can't handle leading spaces?
            #raise

        npoints = user_points.shape[0]
        if npoints == 0:
            raise RuntimeError('npoints=0 in %r' % points_filename)
        if len(user_points.shape) == 1:
            user_points = user_points.reshape(1, npoints)

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
