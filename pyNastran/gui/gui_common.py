# coding: utf-8
# pylint: disable=W0201,C0111
from __future__ import division, unicode_literals, print_function

# standard library
import sys
import os.path
import datetime
from copy import deepcopy
from collections import OrderedDict
from itertools import cycle
from math import ceil
import cgi #  html lib

from six import string_types, iteritems, itervalues
from six.moves import range

import numpy as np

from pyNastran.gui.qt_version import qt_version

from qtpy import QtCore, QtGui #, API
from qtpy.QtWidgets import (
    QMessageBox, QWidget,
    QMainWindow, QDockWidget, QFrame, QHBoxLayout, QAction)

import vtk


import pyNastran
#print('qt_version = %r' % qt_version)
if qt_version in ['pyside', 'pyqt4']:
    from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
else:
    #from vtk.qt5.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
    from pyNastran.gui.qt_files.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from pyNastran.utils.log import SimpleLogger
from pyNastran.utils import integer_types

from pyNastran.gui.qt_files.gui_qt_common import GuiCommon
from pyNastran.gui.qt_files.scalar_bar import ScalarBar

from pyNastran.gui.gui_objects.coord_properties import CoordProperties
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry

#from pyNastran.gui.styles.trackball_style_camera import TrackballStyleCamera
from pyNastran.gui.menus.menus import (
    set_legend_menu, get_legend_fringe, get_legend_disp, get_legend_vector,
    set_clipping_menu,
    set_camera_menu,
    set_preferences_menu,
    on_set_modify_groups, Group,
    Sidebar,
    ApplicationLogWidget,
    PythonConsoleWidget,
    EditGeometryProperties)

from pyNastran.gui.utils.write_gif import (
    setup_animation, update_animation_inputs, write_gif)
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_idtype


#from pyNastran.gui.menus.multidialog import MultiFileDialog
from pyNastran.gui.formats import CLASS_MAP

#class Interactor(vtk.vtkGenericRenderWindowInteractor):
    #def __init__(self):
        ##vtk.vtkGenericRenderWindowInteractor()
        #pass

    #def HighlightProp(self):
        #print('highlight')


#class PyNastranRenderWindowInteractor(QVTKRenderWindowInteractor):
    #def __init__(self, parent=None):

        #render_window = vtk.vtkRenderWindow()
        #iren = Interactor()
        #iren.SetRenderWindow(render_window)
        #kwargs = {
            #'iren' : iren,
            #'rw' : render_window,
        #}
        #QVTKRenderWindowInteractor.__init__(self, parent=parent,
                                            #iren=iren, rw=render_window)
        #self.Highlight

# http://pyqt.sourceforge.net/Docs/PyQt5/multiinheritance.html
class GuiCommon2(QMainWindow, GuiCommon):
    def __init__(self, **kwds):
        """
        fmt_order, html_logging, inputs, parent=None,
        """
        # this will reset the background color/label color if things break
        #super(QMainWindow, self).__init__(self)
        if qt_version == 'pyqt4':
            QMainWindow.__init__(self)
            GuiCommon.__init__(self, **kwds)
        elif qt_version == 'pyqt5':
            super(GuiCommon2, self).__init__(**kwds)
        elif qt_version == 'pyside':
            #super(GuiCommon2, self).__init__(**kwds) # fails

            # fails
            #QMainWindow.__init__(self)
            #GuiCommon.__init__(self, **kwds)

            #super(GuiCommon2, self).__init__(**kwds)
            #super(GuiCommon2, self).__init__(**kwds)

            #super(GuiCommon2, self).__init__(**kwds)

            QMainWindow.__init__(self)
            GuiCommon.__init__(self, **kwds)
        else:
            raise NotImplementedError(qt_version)

        self.format_class_map = CLASS_MAP
        fmt_order = kwds['fmt_order']
        inputs = kwds['inputs']

        #self.app = inputs['app']
        #del inputs['app']


        if inputs['log'] is not None:
            html_logging = False
        else:
            html_logging = kwds['html_logging']
        del kwds['html_logging']

        #if qt_version == 4:  # TODO: remove this???
            #QMainWindow.__init__(self)

        #-----------------------------------------------------------------------
        self._active_background_image = None
        self.reset_settings = False
        self.fmts = fmt_order
        self.base_window_title = "pyNastran v%s"  % pyNastran.__version__

        #defaults
        self.wildcard_delimited = 'Delimited Text (*.txt; *.dat; *.csv)'

        # initializes tools/checkables
        self.set_tools()

        self.html_logging = html_logging
        self.execute_python = True

        self.scalar_bar = ScalarBar(self.is_horizontal_scalar_bar)

        # in,lb,s
        self.input_units = ['', '', ''] # '' means not set
        self.display_units = ['', '', '']
        #self.recent_files = []

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

    def Render(self):
        #self.vtk_interactor.Render()
        self.vtk_interactor.GetRenderWindow().Render()

    @property
    def legend_shown(self):
        """determines if the legend is shown"""
        return self.scalar_bar.is_shown

    #@legend_shown.setter
    #def legend_shown(self):
        #"""show/hide the legend shown"""
        #return self.scalar_bar.is_shown

    @property
    def scalarBar(self):
        return self.scalar_bar.scalar_bar

    @property
    def color_function(self):
        return self.scalar_bar.color_function

    #def get_color_function(self):
        #return self.scalar_bar.color_function

    @property
    def logo(self):
        """Gets the pyNastran icon path, which can be overwritten"""
        return self._logo

    @logo.setter
    def logo(self, logo):
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
        self.window_title = self.base_window_title

        #=========== Results widget ===================
        self.res_dock = QDockWidget("Results", self)
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

        if self.run_vtk:
            self.build_vtk_frame()

        #compassRepresentation = vtk.vtkCompassRepresentation()
        #compassWidget = vtk.vtkCompassWidget()
        #compassWidget.SetInteractor(self.vtk_interactor)
        #compassWidget.SetRepresentation(compassRepresentation)
        #compassWidget.EnabledOn()

    def create_log_python_docks(self):
        """
        Creates the
         - HTML Log dock
         - Python Console dock
        """
        #=========== Logging widget ===================

        if self.html_logging is True:
            self.log_dock_widget = ApplicationLogWidget(self)
            self.log_widget = self.log_dock_widget.log_widget
            self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock_widget)
        else:
            self.log_widget = self.log

        if self.execute_python:
            self.python_dock_widget = PythonConsoleWidget(self)
            self.python_dock_widget.setObjectName("python_console")
            self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.python_dock_widget)

    def _on_execute_python_button(self, clear=False):
        """executes the docked python console"""
        txt = str(self.python_dock_widget.enter_data.toPlainText()).rstrip()
        is_passed = self._execute_python_code(txt)
        if is_passed and clear:
            self.python_dock_widget.enter_data.clear()

    def set_tools(self, tools=None, checkables=None):
        """Creates the GUI tools"""
        if checkables is None:
            checkables = {
                # name, is_checked
                'show_info' : True,
                'show_debug' : True,
                'show_command' : True,
                'anti_alias_0' : True,
                'anti_alias_1' : False,
                'anti_alias_2' : False,
                'anti_alias_4' : False,
                'anti_alias_8' : False,

                'rotation_center' : False,
                'measure_distance' : False,
                'probe_result' : False,
                'area_pick' : False,
                'zoom' : False,
            }

        if tools is None:
            file_tools = [

                ('exit', '&Exit', 'texit.png', 'Ctrl+Q', 'Exit application', self.closeEvent), # QtGui.qApp.quit

                ('load_geometry', 'Load &Geometry...', 'load_geometry.png', 'Ctrl+O', 'Loads a geometry input file', self.on_load_geometry),
                ('load_results', 'Load &Results...', 'load_results.png', 'Ctrl+R', 'Loads a results file', self.on_load_results),
                ('load_csv_user_geom', 'Load CSV User Geometry...', '', None, 'Loads custom geometry file', self.on_load_user_geom),
                ('load_csv_user_points', 'Load CSV User Points...', 'user_points.png', None, 'Loads CSV points', self.on_load_csv_points),
                ('load_custom_result', 'Load Custom Results...', '', None, 'Loads a custom results file', self.on_load_custom_results),

                ('script', 'Run Python Script...', 'python48.png', None, 'Runs pyNastranGUI in batch mode', self.on_run_script),
            ]

            tools = file_tools + [
                ('label_clear', 'Clear Current Labels', '', None, 'Clear current labels', self.clear_labels),
                ('label_reset', 'Clear All Labels', '', None, 'Clear all labels', self.reset_labels),

                ('legend', 'Modify Legend...', 'legend.png', None, 'Set Legend', self.set_legend),
                ('clipping', 'Set Clipping...', '', None, 'Set Clipping', self.set_clipping),
                #('axis', 'Show/Hide Axis', 'axis.png', None, 'Show/Hide Global Axis', self.on_show_hide_axes),

                ('wireframe', 'Wireframe Model', 'twireframe.png', 'w', 'Show Model as a Wireframe Model', self.on_wireframe),
                ('surface', 'Surface Model', 'tsolid.png', 's', 'Show Model as a Surface Model', self.on_surface),
                ('geo_properties', 'Edit Geometry Properties...', '', None, 'Change Model Color/Opacity/Line Width', self.edit_geometry_properties),
                ('modify_groups', 'Modify Groups...', '', None, 'Create/Edit/Delete Groups', self.on_set_modify_groups),

                ('create_groups_by_visible_result', 'Create Groups By Visible Result', '', None, 'Create Groups', self.create_groups_by_visible_result),
                ('create_groups_by_property_id', 'Create Groups By Property ID', '', None, 'Create Groups', self.create_groups_by_property_id),
                #('create_list', 'Create Lists through Booleans', '', None, 'Create List', self.create_list),

                ('show_info', 'Show INFO', 'show_info.png', None, 'Show "INFO" messages', self.on_show_info),
                ('show_debug', 'Show DEBUG', 'show_debug.png', None, 'Show "DEBUG" messages', self.on_show_debug),
                ('show_command', 'Show COMMAND', 'show_command.png', None, 'Show "COMMAND" messages', self.on_show_command),

                ('magnify', 'Magnify', 'plus_zoom.png', 'm', 'Increase Magnfication', self.on_increase_magnification),
                ('shrink', 'Shrink', 'minus_zoom.png', 'Shift+M', 'Decrease Magnfication', self.on_decrease_magnification),

                #('cell_pick', 'Cell Pick', '', 'c', 'Centroidal Picking', self.on_cell_picker),
                #('node_pick', 'Node Pick', '', 'n', 'Nodal Picking', self.on_node_picker),

                ('rotate_clockwise', 'Rotate Clockwise', 'tclock.png', 'o', 'Rotate Clockwise', self.on_rotate_clockwise),
                ('rotate_cclockwise', 'Rotate Counter-Clockwise', 'tcclock.png', 'Shift+O', 'Rotate Counter-Clockwise', self.on_rotate_cclockwise),

                ('screenshot', 'Take a Screenshot...', 'tcamera.png', 'CTRL+I', 'Take a Screenshot of current view', self.tool_actions.on_take_screenshot),
                ('about', 'About pyNastran GUI...', 'tabout.png', 'CTRL+H', 'About pyNastran GUI and help on shortcuts', self.about_dialog),
                ('view', 'Camera View', 'view.png', None, 'Load the camera menu', self.view_camera),
                ('camera_reset', 'Reset Camera View', 'trefresh.png', 'r', 'Reset the camera view to default', self.on_reset_camera),
                ('reload', 'Reload Model...', 'treload.png', '', 'Remove the model and reload the same geometry file', self.on_reload),

                ('cycle_results', 'Cycle Results', 'cycle_results.png', 'CTRL+L', 'Changes the result case', self.on_cycle_results),
                ('rcycle_results', 'Cycle Results', 'rcycle_results.png', 'CTRL+K', 'Changes the result case', self.on_rcycle_results),

                ('back_view', 'Back View', 'back.png', 'x', 'Flips to +X Axis', lambda: self.view_actions.update_camera('+x')),
                ('right_view', 'Right View', 'right.png', 'y', 'Flips to +Y Axis', lambda: self.view_actions.update_camera('+y')),
                ('top_view', 'Top View', 'top.png', 'z', 'Flips to +Z Axis', lambda: self.view_actions.update_camera('+z')),

                ('front_view', 'Front View', 'front.png', 'Shift+X', 'Flips to -X Axis', lambda: self.view_actions.update_camera('-x')),
                ('left_view', 'Left View', 'left.png', 'Shift+Y', 'Flips to -Y Axis', lambda: self.view_actions.update_camera('-y')),
                ('bottom_view', 'Bottom View', 'bottom.png', 'Shift+Z', 'Flips to -Z Axis', lambda: self.view_actions.update_camera('-z')),
                ('edges', 'Show/Hide Edges', 'tedges.png', 'e', 'Show/Hide Model Edges', self.on_flip_edges),
                ('edges_black', 'Color Edges', '', 'b', 'Set Edge Color to Color/Black', self.on_set_edge_visibility),
                ('anti_alias_0', 'Off', '', None, 'Disable Anti-Aliasing', lambda: self.on_set_anti_aliasing(0)),
                ('anti_alias_1', '1x', '', None, 'Set Anti-Aliasing to 1x', lambda: self.on_set_anti_aliasing(1)),
                ('anti_alias_2', '2x', '', None, 'Set Anti-Aliasing to 2x', lambda: self.on_set_anti_aliasing(2)),
                ('anti_alias_4', '4x', '', None, 'Set Anti-Aliasing to 4x', lambda: self.on_set_anti_aliasing(4)),
                ('anti_alias_8', '8x', '', None, 'Set Anti-Aliasing to 8x', lambda: self.on_set_anti_aliasing(8)),

                # new
                ('rotation_center', 'Set the rotation center', 'trotation_center.png', 'f', 'Pick a node for the rotation center', self.mouse_actions.on_rotation_center),

                ('measure_distance', 'Measure Distance', 'measure_distance.png', None, 'Measure the distance between two nodes', self.mouse_actions.on_measure_distance),
                ('probe_result', 'Probe', 'tprobe.png', None, 'Probe the displayed result', self.mouse_actions.on_probe_result),
                ('quick_probe_result', 'Quick Probe', '', 'p', 'Probe the displayed result', self.mouse_actions.on_quick_probe_result),
                ('zoom', 'Zoom', 'zoom.png', None, 'Zoom In', self.mouse_actions.on_zoom),
                ('font_size_increase', 'Increase Font Size', 'text_up.png', 'Ctrl+Plus', 'Increase Font Size', self.on_increase_font_size),
                ('font_size_decrease', 'Decrease Font Size', 'text_down.png', 'Ctrl+Minus', 'Decrease Font Size', self.on_decrease_font_size),
                ('set_preferences', 'Preferences...', 'preferences.png', None, 'Set GUI Preferences', self.set_preferences_menu),

                # picking
                ('area_pick', 'Area Pick', 'tarea_pick.png', None, 'Get a list of nodes/elements', self.mouse_actions.on_area_pick),
            ]

        if 'nastran' in self.fmts:
            tools += [
                ('caero', 'Show/Hide CAERO Panels', '', None, 'Show/Hide CAERO Panel Outlines', self.toggle_caero_panels),
                ('caero_subpanels', 'Toggle CAERO Subpanels', '', None, 'Show/Hide CAERO Subanel Outlines', self.toggle_caero_sub_panels),
                ('conm2', 'Toggle CONM2s', '', None, 'Show/Hide CONM2s', self.toggle_conms),
            ]
        self.tools = tools
        self.checkables = checkables

    def keyPressEvent(self, qkey_event):
        #print('qkey_event =', qkey_event.key())
        super(GuiCommon2, self).keyPressEvent(qkey_event)

    def on_increase_font_size(self):
        """used by the hidden_tools for Ctrl +"""
        self.on_set_font_size(self.settings.font_size + 1)

    def on_decrease_font_size(self):
        """used by the hidden_tools for Ctrl -"""
        self.on_set_font_size(self.settings.font_size - 1)

    def on_set_font_size(self, font_size, show_command=True):
        """changes the font size"""
        is_failed = True
        if not isinstance(font_size, int):
            self.log_error('font_size=%r must be an integer; type=%s' % (
                font_size, type(font_size)))
            return is_failed
        if font_size < 6:
            font_size = 6
        if self.settings.font_size == font_size:
            return False
        self.settings.font_size = font_size
        font = QtGui.QFont()
        font.setPointSize(self.settings.font_size)
        self.setFont(font)

        #self.toolbar.setFont(font)
        self.menu_file.setFont(font)
        self.menu_view.setFont(font)
        self.menu_window.setFont(font)
        self.menu_help.setFont(font)

        if self._legend_window_shown:
            self._legend_window.set_font_size(font_size)
        if self._clipping_window_shown:
            self._clipping_window.set_font_size(font_size)
        if self._edit_geometry_properties_window_shown:
            self._edit_geometry_properties.set_font_size(font_size)
        if self._modify_groups_window_shown:
            self._modify_groups_window.set_font_size(font_size)
        if self._preferences_window_shown:
            self._preferences_window.set_font_size(font_size)

        #self.menu_scripts.setFont(font)
        self.log_command('settings.on_set_font_size(%s)' % font_size)
        return False

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
            'set_preferences', '',
            'label_clear', 'label_reset', '',
            'legend', 'geo_properties',
            #['Anti-Aliasing', 'anti_alias_0', 'anti_alias_1', 'anti_alias_2',
            #'anti_alias_4', 'anti_alias_8',],
        ]
        if self.is_groups:
            menu_view += ['modify_groups', 'create_groups_by_property_id',
                          'create_groups_by_visible_result']
        menu_view += [
            '', 'clipping', #'axis',
            'edges', 'edges_black',]
        if self.html_logging:
            self.actions['log_dock_widget'] = self.log_dock_widget.toggleViewAction()
            self.actions['log_dock_widget'].setStatusTip("Show/Hide application log")
            menu_view += ['', 'show_info', 'show_debug', 'show_command']
            menu_window += ['log_dock_widget']
        if self.execute_python:
            self.actions['python_dock_widget'] = self.python_dock_widget.toggleViewAction()
            self.actions['python_dock_widget'].setStatusTip("Show/Hide Python Console")
            menu_window += ['python_dock_widget']

        menu_file = [
            'load_geometry', 'load_results', '',
            'load_custom_result', '',
            'load_csv_user_points', 'load_csv_user_geom', 'script', '', 'exit']
        toolbar_tools = ['reload', 'load_geometry', 'load_results',
                         'front_view', 'back_view', 'top_view', 'bottom_view', 'left_view', 'right_view',
                         'magnify', 'shrink', 'zoom',
                         'rotate_clockwise', 'rotate_cclockwise',
                         'rotation_center', 'measure_distance', 'probe_result', 'area_pick',

                         'wireframe', 'surface', 'edges']
        toolbar_tools += ['camera_reset', 'view', 'screenshot', '', 'exit']
        hidden_tools = ('cycle_results', 'rcycle_results',
                        'font_size_increase', 'font_size_decrease')

        menu_items = OrderedDict()
        if create_menu_bar:
            menu_items['file'] = (self.menu_file, menu_file)
            menu_items['view'] = (self.menu_view, menu_view)
            menu_items['main'] = (self.menu_window, menu_window)
            menu_items['help'] = (self.menu_help, ('about',))
            menu_items['scripts'] = (self.menu_scripts, scripts)
            menu_items['toolbar'] = (self.toolbar, toolbar_tools)
            menu_items['hidden'] = (self.menu_hidden, hidden_tools)
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
        assert isinstance(menu_items, dict), menu_items
        for unused_menu_name, (menu, items) in iteritems(menu_items):
            if menu is None:
                continue
            for i in items:
                if not i:
                    menu.addSeparator()
                else:
                    if isinstance(i, list):
                        sub_menu_name = i[0]
                        sub_menu = menu.addMenu(sub_menu_name)
                        for ii_count, ii in enumerate(i[1:]):
                            if not isinstance(ii, string_types):
                                raise RuntimeError('what is this...action ii() = %r' % ii())
                            action = self.actions[ii]
                            if ii_count > 0:
                                action.setChecked(False)
                            sub_menu.addAction(action)
                        continue
                    elif not isinstance(i, string_types):
                        raise RuntimeError('what is this...action i() = %r' % i())

                    try:
                        action = self.actions[i] #if isinstance(i, string_types) else i()
                    except:
                        print(self.actions.keys())
                        raise
                    menu.addAction(action)
        #self._create_plane_from_points(None)

    def _update_menu(self, menu_items):
        assert isinstance(menu_items, dict), menu_items
        for unused_name, (menu, unused_items) in iteritems(menu_items):
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

        for tool in tools:
            (name, txt, icon, shortcut, tip, func) = tool
            if name in self.actions:
                self.log_error('trying to create a duplicate action %r' % name)
                continue

            if icon is None:
                print("missing_icon = %r!!!" % name)
                ico = None
            else:
                ico = QtGui.QIcon()
                pth = os.path.join(icon_path, icon)
                ico.addPixmap(QtGui.QPixmap(pth), QtGui.QIcon.Normal, QtGui.QIcon.Off)

            if name in checkables:
                is_checked = checkables[name]
                self.actions[name] = QAction(ico, txt, self, checkable=True)
                self.actions[name].setChecked(is_checked)
            else:
                self.actions[name] = QAction(ico, txt, self)

            if shortcut:
                self.actions[name].setShortcut(shortcut)
                #actions[name].setShortcutContext(QtCore.Qt.WidgetShortcut)
            if tip:
                self.actions[name].setStatusTip(tip)
            if func:
                self.actions[name].triggered.connect(func)

        self.actions['toolbar'] = self.toolbar.toggleViewAction()
        self.actions['toolbar'].setStatusTip("Show/Hide application toolbar")

        self.actions['reswidget'] = self.res_dock.toggleViewAction()
        self.actions['reswidget'].setStatusTip("Show/Hide results selection")
        return self.actions

    def _logg_msg(self, typ, msg):
        """
        Add message to log widget trying to choose right color for it.

        Parameters
        ----------
        typ : str
            {DEBUG, INFO, ERROR, COMMAND, WARNING} or prepend 'GUI '
        msg : str
            message to be displayed
        """
        if not self.html_logging:
            print(typ, msg)
            return

        if 'DEBUG' in typ and not self.settings.show_debug:
            return
        elif 'INFO' in typ and not self.settings.show_info:
            return
        elif 'COMMAND' in typ and not self.settings.show_command:
            return

        frame = sys._getframe(4)  # jump to get out of the logger code
        lineno = frame.f_lineno
        filename = os.path.basename(frame.f_globals['__file__'])

        if typ in ['GUI ERROR', 'GUI COMMAND', 'GUI DEBUG', 'GUI INFO', 'GUI WARNING']:
            typ = typ[4:]
            msg = '   fname=%-25s:%-4s   %s\n' % (filename, lineno, msg)

        tim = datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
        msg = cgi.escape(msg)

        #message colors
        dark_orange = '#EB9100'
        colors = {
            'COMMAND' : 'blue',
            'ERROR' : 'Crimson',
            'DEBUG' : dark_orange,
            'WARNING' : 'purple',
            'INFO' : 'green',
        }
        msg = msg.rstrip().replace('\n', '<br>')
        msg = tim + ' ' + (typ + ': ' + msg) if typ else msg
        if typ in colors:
            msg = '<font color="%s"> %s </font>' % (colors[typ], msg)

        self.log_mutex.lockForWrite()
        text_cursor = self.log_widget.textCursor()
        end = text_cursor.End
        text_cursor.movePosition(end)
        text_cursor.insertHtml(msg + r"<br />")
        self.log_widget.ensureCursorVisible() # new message will be visible
        self.log_mutex.unlock()

    def log_info(self, msg):
        """ Helper funtion: log a message msg with a 'INFO:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'ERROR')
        self.log.simple_msg(msg, 'INFO')

    def log_debug(self, msg):
        """ Helper funtion: log a message msg with a 'DEBUG:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        self.log.simple_msg(msg, 'GUI DEBUG')

    def log_command(self, msg):
        """ Helper funtion: log a message msg with a 'COMMAND:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        self.log.simple_msg(msg, 'GUI COMMAND')

    def log_error(self, msg):
        """ Helper funtion: log a message msg with a 'GUI ERROR:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        self.log.simple_msg(msg, 'GUI ERROR')

    def log_warning(self, msg):
        """ Helper funtion: log a message msg with a 'WARNING:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        self.log.simple_msg(msg, 'GUI WARNING')

    def create_vtk_actors(self):
        self.rend = vtk.vtkRenderer()

        # vtk actors
        self.grid = vtk.vtkUnstructuredGrid()

        # edges
        self.edge_actor = vtk.vtkLODActor()
        self.edge_actor.DragableOff()
        self.edge_mapper = vtk.vtkPolyDataMapper()

        self.create_cell_picker()

    def _create_vtk_objects(self):
        """creates some of the vtk objects"""
        #Frame that VTK will render on
        self.vtk_frame = QFrame()

        #Qt VTK QVTKRenderWindowInteractor
        self.vtk_interactor = QVTKRenderWindowInteractor(parent=self.vtk_frame)
        #self.vtk_interactor = PyNastranRenderWindowInteractor(parent=self.vtk_frame)
        #self.set_anti_aliasing(2)

        #self._camera_event_name = 'LeftButtonPressEvent'
        self.mouse_actions.setup_mouse_buttons(mode='default')

    def on_escape_null(self):
        """
        The default state for Escape key is nothing.
        """
        pass

    def on_escape(self):
        """
        Escape key should cancel:
         - on_rotation_center

        TODO: not done...
        """
        pass

    #def remove_picker(self):
        #self.vtk_interactor.

    def set_node_picker(self):
        self.vtk_interactor.SetPicker(self.node_picker)

    def set_cell_picker(self):
        self.vtk_interactor.SetPicker(self.cell_picker)

    def set_background_image(self, image_filename='GeologicalExfoliationOfGraniteRock.jpg'):
        """adds a background image"""
        if not os.path.exists(image_filename):
            return

        fmt = os.path.splitext(image_filename)[1].lower()
        if fmt not in ['.jpg', '.jpeg', '.png', '.tif', '.tiff', '.bmp']:
            msg = 'invalid image type=%r; filename=%r' % (fmt, image_filename)
            raise NotImplementedError(msg)

        #image_reader = vtk.vtkJPEGReader()
        #image_reader = vtk.vtkPNGReader()
        #image_reader = vtk.vtkTIFFReader()
        #image_reader = vtk.vtkBMPReader()
        #image_reader = vtk.vtkPostScriptReader()  # doesn't exist?

        has_background_image = self._active_background_image is not None
        self._active_background_image = image_filename
        #if has_background_image:
            #self.image_reader.Delete()


        if fmt in ['.jpg', '.jpeg']:
            self.image_reader = vtk.vtkJPEGReader()
        elif fmt == '.png':
            self.image_reader = vtk.vtkPNGReader()
        elif fmt in ['.tif', '.tiff']:
            self.image_reader = vtk.vtkTIFFReader()
        elif fmt == '.bmp':
            self.image_reader = vtk.vtkBMPReader()
        #elif fmt == '.ps': # doesn't exist?
            #self.image_reader = vtk.vtkPostScriptReader()
        else:
            msg = 'invalid image type=%r; filename=%r' % (fmt, image_filename)
            raise NotImplementedError(msg)

        if not self.image_reader.CanReadFile(image_filename):
            print("Error reading file %s" % image_filename)
            return

        self.image_reader.SetFileName(image_filename)
        self.image_reader.Update()
        image_data = self.image_reader.GetOutput()

        if has_background_image:
            self.image_actor.SetInputData(image_data)
            self.Render()
            return

        # Create an image actor to display the image
        self.image_actor = vtk.vtkImageActor()
        self.image_actor.SetInputData(image_data)

        self.background_rend = vtk.vtkRenderer()
        self.background_rend.SetLayer(0)
        self.background_rend.InteractiveOff()
        self.background_rend.AddActor(self.image_actor)

        self.rend.SetLayer(1)
        render_window = self.vtk_interactor.GetRenderWindow()
        render_window.SetNumberOfLayers(2)

        render_window.AddRenderer(self.background_rend)

        # Set up the background camera to fill the renderer with the image
        origin = image_data.GetOrigin()
        spacing = image_data.GetSpacing()
        extent = image_data.GetExtent()

        camera = self.background_rend.GetActiveCamera()
        camera.ParallelProjectionOn()

        xc = origin[0] + 0.5*(extent[0] + extent[1]) * spacing[0]
        yc = origin[1] + 0.5*(extent[2] + extent[3]) * spacing[1]
        #xd = (extent[1] - extent[0] + 1) * spacing[0]
        yd = (extent[3] - extent[2] + 1) * spacing[1]
        d = camera.GetDistance()
        camera.SetParallelScale(0.5 * yd)
        camera.SetFocalPoint(xc, yc, 0.0)
        camera.SetPosition(xc, yc, d)

    def build_vtk_frame(self):
        vtk_hbox = QHBoxLayout()
        vtk_hbox.setContentsMargins(2, 2, 2, 2)

        vtk_hbox.addWidget(self.vtk_interactor)
        self.vtk_frame.setLayout(vtk_hbox)
        self.vtk_frame.setFrameStyle(QFrame.NoFrame | QFrame.Plain)
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

        self.set_background_image()
        self.vtk_interactor.GetRenderWindow().Render()
        #self.load_nastran_geometry(None, None)

        #for cid, axes in iteritems(self.axes):
            #self.rend.AddActor(axes)
        self.add_geometry()
        if nframes == 2:
            rend.AddActor(self.geom_actor)

        # initialize geometry_actors
        self.geometry_actors['main'] = self.geom_actor

        # bar scale set so you can't edit the bar scale
        white = (255, 255, 255)
        geom_props = AltGeometry(
            self, 'main', color=white, line_width=1, opacity=1.0, point_size=1,
            bar_scale=0.0, representation='main', is_visible=True)

        self.geometry_properties['main'] = geom_props

        #self.addAltGeometry()
        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(*self.settings.background_color)
        #self.rend.SetBackground2(*self.background_color2)

        self.rend.ResetCamera()
        self.mouse_actions.set_style_as_trackball()
        self.build_lookup_table()

        text_size = 14
        dtext_size = text_size + 1
        self.create_text([5, 5 + 3 * dtext_size], 'Max  ', text_size)  # text actor 0
        self.create_text([5, 5 + 2 * dtext_size], 'Min  ', text_size)  # text actor 1
        self.create_text([5, 5 + 1 * dtext_size], 'Word1', text_size)  # text actor 2
        self.create_text([5, 5], 'Word2', text_size)  # text actor 3

        self.get_edges()
        if self.is_edges:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOn()
        else:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOff()

    def on_show_info(self):
        """sets a flag for showing/hiding INFO messages"""
        self.settings.show_info = not self.settings.show_info

    def on_show_debug(self):
        """sets a flag for showing/hiding DEBUG messages"""
        self.settings.show_debug = not self.settings.show_debug

    def on_show_command(self):
        """sets a flag for showing/hiding COMMAND messages"""
        self.settings.show_command = not self.settings.show_command

    def on_reset_camera(self):
        self.log_command('on_reset_camera()')
        self._simulate_key_press('r')
        self.vtk_interactor.Render()

    def on_flip_edges(self):
        """turn edges on/off"""
        self.is_edges = not self.is_edges
        self.edge_actor.SetVisibility(self.is_edges)
        #self.edge_actor.GetProperty().SetColor(0, 0, 0)  # cart3d edge color isn't black...
        self.edge_actor.Modified()
        #self.widget.Update()
        #self._update_camera()
        self.Render()
        #self.refresh()
        self.log_command('on_flip_edges()')

    def on_set_edge_visibility(self):
        #self.edge_actor.SetVisibility(self.is_edges_black)
        self.is_edges_black = not self.is_edges_black
        if self.is_edges_black:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOn()
            self.edge_mapper.SetLookupTable(self.color_function_black)
        else:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOff()
            self.edge_mapper.SetLookupTable(self.color_function)
        self.edge_actor.Modified()
        prop.Modified()
        self.vtk_interactor.Render()
        self.log_command('on_set_edge_visibility()')

    def get_edges(self):
        """Create the edge actor"""
        edges = vtk.vtkExtractEdges()
        edge_mapper = self.edge_mapper
        edge_actor = self.edge_actor

        edges.SetInputData(self.grid_selected)
        edge_mapper.SetInputConnection(edges.GetOutputPort())

        edge_actor.SetMapper(edge_mapper)
        edge_actor.GetProperty().SetColor(0., 0., 0.)
        edge_mapper.SetLookupTable(self.color_function)
        edge_mapper.SetResolveCoincidentTopologyToPolygonOffset()

        prop = edge_actor.GetProperty()
        prop.SetColor(0., 0., 0.)
        edge_actor.SetVisibility(self.is_edges)
        self.rend.AddActor(edge_actor)

    def post_group_by_name(self, name):
        """posts a group with a specific name"""
        group = self.groups[name]
        self.post_group(group)
        self.group_active = name

    def post_group(self, group):
        """posts a group object"""
        eids = group.element_ids
        self.show_eids(eids)

    def get_all_eids(self):
        """get the list of all the element IDs"""
        return self.element_ids
        #name, result = self.get_name_result_data(0)
        #if name != 'ElementID':
            #name, result = self.get_name_result_data(1)
            #assert name == 'ElementID', name
        #return result

    def show_eids(self, eids):
        """shows the specified element IDs"""
        all_eids = self.get_all_eids()

        # remove eids that are out of range
        eids = np.intersect1d(all_eids, eids)

        # update for indices
        ishow = np.searchsorted(all_eids, eids)

        #eids_off = np.setdiff1d(all_eids, eids)
        #j = np.setdiff1d(all_eids, eids_off)

        self.show_ids_mask(ishow)

    def hide_eids(self, eids):
        """hides the specified element IDs"""
        all_eids = self.get_all_eids()

        # remove eids that are out of range
        eids = np.intersect1d(all_eids, eids)

        # A-B
        eids = np.setdiff1d(all_eids, eids)

        # update for indices
        ishow = np.searchsorted(all_eids, eids)
        self.show_ids_mask(ishow)

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
        """masks the specific 0-based element ids"""
        #print('ids_to_show = ', ids_to_show)
        prop = self.geom_actor.GetProperty()
        if len(ids_to_show) == self.nelements:
            #prop.BackfaceCullingOn()
            pass
        else:
            prop.BackfaceCullingOff()

        if 0:  # pragma: no cover
            self._show_ids_mask(ids_to_show)
        elif 1:
            # doesn't work for the BWB_saero.bdf
            flip_flag = True is self._show_flag
            assert self._show_flag is True, self._show_flag
            self._update_ids_mask_show(ids_to_show)
            self._show_flag = True
        elif 1:  # pragma: no cover
            # works
            flip_flag = True is self._show_flag
            assert self._show_flag is True, self._show_flag
            self._update_ids_mask_show_true(ids_to_show, flip_flag, render=False)
            self._update_ids_mask_show_true(ids_to_show, False, render=True)
            self._show_flag = True
        else:  # pragma: no cover
            # old; works; slow
            flip_flag = True is self._show_flag
            self._update_ids_mask(ids_to_show, flip_flag, show_flag=True, render=False)
            self._update_ids_mask(ids_to_show, False, show_flag=True, render=True)
            self._show_flag = True

    def hide_ids_mask(self, ids_to_hide):
        """masks the specific 0-based element ids"""
        #print('hide_ids_mask = ', hide_ids_mask)
        prop = self.geom_actor.GetProperty()
        if len(self.ids_to_hide) == 0:
            prop.BackfaceCullingOn()
        else:
            prop.BackfaceCullingOff()

        #if 0:  # pragma: no cover
        #self._hide_ids_mask(ids_to_hide)
        #else:
        # old; works; slow
        flip_flag = False is self._show_flag
        self._update_ids_mask(ids_to_hide, flip_flag, show_flag=False, render=False)
        self._update_ids_mask(ids_to_hide, False, show_flag=False, render=True)
        self._show_flag = False

    def _show_ids_mask(self, ids_to_show):
        """
        helper method for ``show_ids_mask``
        .. todo:: doesn't work
        """
        all_i = np.arange(self.nelements, dtype='int32')
        ids_to_hide = np.setdiff1d(all_i, ids_to_show)
        self._hide_ids_mask(ids_to_hide)

    def _hide_ids_mask(self, ids_to_hide):
        """
        helper method for ``hide_ids_mask``
        .. todo:: doesn't work
        """
        #print('_hide_ids_mask = ', ids_to_hide)
        ids = numpy_to_vtk_idtype(ids_to_hide)

        #self.selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)

        if 1:
            # sane; doesn't work
            self.selection_node.SetSelectionList(ids)
            ids.Modified()
            self.selection_node.Modified()
            self.selection.Modified()
            self.grid_selected.Modified()
            self.grid_selected.ShallowCopy(self.extract_selection.GetOutput())
            self.update_all(render=True)
        else:  # pragma: no cover
            # doesn't work
            self.selection.RemoveAllNodes()
            self.selection_node = vtk.vtkSelectionNode()
            self.selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
            self.selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
            #self.selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
            self.selection.AddNode(self.selection_node)
            self.selection_node.SetSelectionList(ids)

            #self.selection.RemoveAllNodes()
            #self.selection.AddNode(self.selection_node)
            self.grid_selected.ShallowCopy(self.extract_selection.GetOutput())
            self.selection_node.SetSelectionList(ids)
            self.update_all(render=True)

    def _update_ids_mask_show_false(self, ids_to_show, flip_flag=True, render=True):
        ids = numpy_to_vtk_idtype(ids_to_show)
        ids.Modified()

        if flip_flag:
            self.selection.RemoveAllNodes()
            self.selection_node = vtk.vtkSelectionNode()
            self.selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
            self.selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
            self.selection_node.SetSelectionList(ids)

            self.selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
            self.selection.AddNode(self.selection_node)
        else:
            self.selection_node.SetSelectionList(ids)

        # dumb; works
        self.grid_selected.ShallowCopy(self.extract_selection.GetOutput())
        self.update_all(render=render)

    def _update_ids_mask_show(self, ids_to_show):
        """helper method for ``show_ids_mask``"""
        ids = numpy_to_vtk_idtype(ids_to_show)
        ids.Modified()

        self.selection.RemoveAllNodes()
        self.selection_node = vtk.vtkSelectionNode()
        self.selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
        self.selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
        self.selection_node.SetSelectionList(ids)
        self.selection_node.Modified()
        self.selection.Modified()

        self.selection.AddNode(self.selection_node)

        # seems to also work
        self.extract_selection.Update()

        self.grid_selected.ShallowCopy(self.extract_selection.GetOutput())

        self.update_all(render=True)
        #if 0:
            #self.grid_selected.Modified()
            #self.vtk_interactor.Render()
            #render_window = self.vtk_interactor.GetRenderWindow()
            #render_window.Render()

    def _update_ids_mask_show_true(self, ids_to_show,
                                   flip_flag=True, render=True):  # pragma: no cover
        ids = numpy_to_vtk_idtype(ids_to_show)
        ids.Modified()

        if flip_flag:
            self.selection.RemoveAllNodes()
            self.selection_node = vtk.vtkSelectionNode()
            self.selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
            self.selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
            self.selection_node.SetSelectionList(ids)

            self.selection.AddNode(self.selection_node)
        else:
            self.selection_node.SetSelectionList(ids)

        # dumb; works
        self.grid_selected.ShallowCopy(self.extract_selection.GetOutput())
        self.update_all(render=render)

    def _update_ids_mask(self, ids_to_show, flip_flag=True, show_flag=True, render=True):
        print('flip_flag=%s show_flag=%s' % (flip_flag, show_flag))

        ids = numpy_to_vtk_idtype(ids_to_show)
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

        if 0:  # pragma: no cover
            # doesn't work...
            self.extract_selection.SetInputData(0, self.grid)
            self.extract_selection.SetInputData(1, self.selection)
        else:
            # dumb; works
            self.grid_selected.ShallowCopy(self.extract_selection.GetOutput())

        #if 0:
            #self.selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
            #self.extract_selection.Update()
        self.update_all(render=render)

    def update_all_2(self, render=True):  # pragma: no cover
        self.grid_selected.Modified()

        self.selection_node.Modified()
        self.selection.Modified()
        self.extract_selection.Update()
        self.extract_selection.Modified()

        self.grid_selected.Modified()
        self.grid_mapper.Update()
        self.grid_mapper.Modified()

        self.vtk_interactor.Modified()
        self.rend.Render()
        self.rend.Modified()

        self.geom_actor.Modified()

        if render:
            self.vtk_interactor.Render()
            render_window = self.vtk_interactor.GetRenderWindow()
            render_window.Render()

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
        self.grid_mapper.Update()
        self.grid_mapper.Modified()
        #selected_actor.Update()
        #selected_actor.Modified()

        #right_renderer.Modified()
        #right_renderer.Update()

        self.vtk_interactor.Modified()
        #interactor.Update()
        #-----------------
        self.rend.Render()
        #interactor.Start()

        self.rend.Modified()

        self.geom_actor.Modified()

        if render:
            self.vtk_interactor.Render()
            render_window = self.vtk_interactor.GetRenderWindow()
            render_window.Render()


    def _setup_element_mask(self, create_grid_selected=True):
        """
        starts the masking

        self.grid feeds in the geometry
        """
        ids = vtk.vtkIdTypeArray()
        ids.SetNumberOfComponents(1)

        # the "selection_node" is really a "selection_element_ids"
        # furthermore, it's an inverse model, so adding elements
        # hides more elements
        self.selection_node = vtk.vtkSelectionNode()
        self.selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
        self.selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
        self.selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)  # added
        self.selection_node.SetSelectionList(ids)

        self.selection = vtk.vtkSelection()
        self.selection.AddNode(self.selection_node)

        self.extract_selection = vtk.vtkExtractSelection()
        self.extract_selection.SetInputData(0, self.grid)
        self.extract_selection.SetInputData(1, self.selection)
        self.extract_selection.Update()

        # In selection
        if create_grid_selected:
            self.grid_selected = vtk.vtkUnstructuredGrid()
            self.grid_selected.ShallowCopy(self.extract_selection.GetOutput())

        #if 0:
        self.selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
        self.extract_selection.Update()

    def build_lookup_table(self):
        scalar_range = self.grid_selected.GetScalarRange()
        self.grid_mapper.SetScalarRange(scalar_range)
        self.grid_mapper.SetLookupTable(self.color_function)
        self.rend.AddActor(self.scalarBar)

    def start_logging(self):
        if self.html_logging is True:
            log = SimpleLogger('debug', 'utf-8', lambda x, y: self._logg_msg(x, y))
            # logging needs synchronizing, so the messages from different
            # threads would not be interleave
            self.log_mutex = QtCore.QReadWriteLock()
        else:
            log = SimpleLogger(
                level='debug', encoding='utf-8',
                #log_func=lambda x, y: print(x, y)  # no colorama
            )
        self.log = log

    def on_load_geometry_button(self, infile_name=None, geometry_format=None, name='main',
                                plot=True, raise_error=False):
        """action version of ``on_load_geometry``"""
        self.on_load_geometry(infile_name=infile_name, geometry_format=geometry_format,
                              name=name, plot=True, raise_error=raise_error)

    def _update_menu_bar_to_format(self, fmt, method):
        self.menu_bar_format = fmt
        tools, menu_items = getattr(self, method)()
        unused_actions = self._prepare_actions(self._icon_path, tools, self.checkables)
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

    #def _load_force(self, out_filename):
        #"""loads a deflection file"""
        #self._load_deflection_force(out_filename, is_deflection=True, is_force=False)

    def setup_gui(self):
        """
        Setup the gui

        1.  starts the logging
        2.  reapplies the settings
        3.  create pickers
        4.  create main vtk actors
        5.  shows the Qt window
        """
        assert self.fmts != [], 'supported_formats=%s' % self.supported_formats
        self.start_logging()
        settings = QtCore.QSettings()
        self.create_vtk_actors()

        # build GUI and restore saved application state
        #nice_blue = (0.1, 0.2, 0.4)
        qpos_default = self.pos()
        unused_pos_default = qpos_default.x(), qpos_default.y()

        self.reset_settings = False
        #if self.reset_settings or qt_version in [5, 'pyside']:
            #self.settings.reset_settings()
        #else:
        self.settings.load(settings)

        self.init_ui()
        if self.reset_settings:
            self.res_dock.toggleViewAction()
        self.init_cell_picker()

        unused_main_window_state = settings.value("mainWindowState")
        self.create_corner_axis()
        #-------------
        # loading
        self.show()

    def setup_post(self, inputs):
        """interface for user defined post-scripts"""
        self.load_batch_inputs(inputs)

        if inputs['user_points'] is not None:
            for fname in inputs['user_points']:
                self.on_load_user_points(fname)

        if inputs['user_geom'] is not None:
            for fname in inputs['user_geom']:
                self.on_load_user_geom(fname)
        #self.set_anti_aliasing(16)

    def create_cell_picker(self):
        """creates the vtk picker objects"""
        self.cell_picker = vtk.vtkCellPicker()
        self.node_picker = vtk.vtkPointPicker()

        self.area_picker = vtk.vtkAreaPicker()  # vtkRenderedAreaPicker?
        self.rubber_band_style = vtk.vtkInteractorStyleRubberBandPick()
        #vtk.vtkInteractorStyleRubberBand2D
        #vtk.vtkInteractorStyleRubberBand3D
        #vtk.vtkInteractorStyleRubberBandZoom
        #vtk.vtkInteractorStyleAreaSelectHover
        #vtk.vtkInteractorStyleDrawPolygon

        #vtk.vtkAngleWidget
        #vtk.vtkAngleRepresentation2D
        #vtk.vtkAngleRepresentation3D
        #vtk.vtkAnnotation
        #vtk.vtkArrowSource
        #vtk.vtkGlyph2D
        #vtk.vtkGlyph3D
        #vtk.vtkHedgeHog
        #vtk.vtkLegendBoxActor
        #vtk.vtkLegendScaleActor
        #vtk.vtkLabelPlacer

        self.cell_picker.SetTolerance(0.001)
        self.node_picker.SetTolerance(0.001)

    def mark_elements_by_different_case(self, eids, icase_result, icase_to_apply):
        """
        Marks a series of elements with custom text labels

        Parameters
        ----------
        eids : int, List[int]
            the elements to apply a message to
        icase_result : int
            the case to draw the result from
        icase_to_apply : int
            the key in label_actors to slot the result into

        TODO: fix the following
        correct   : applies to the icase_to_apply
        incorrect : applies to the icase_result

        Examples
        --------
        .. code-block::

          eids = [16563, 16564, 8916703, 16499, 16500, 8916699,
                  16565, 16566, 8916706, 16502, 16503, 8916701]
          icase_result = 22
          icase_to_apply = 25
          self.mark_elements_by_different_case(eids, icase_result, icase_to_apply)
        """
        if icase_result not in self.label_actors:
            msg = 'icase_result=%r not in label_actors=[%s]' % (
                icase_result, ', '.join(self.label_actors))
            self.log_error(msg)
            return
        if icase_to_apply not in self.label_actors:
            msg = 'icase_to_apply=%r not in label_actors=[%s]' % (
                icase_to_apply, ', '.join(self.label_actors))
            self.log_error(msg)
            return

        eids = np.unique(eids)
        unused_neids = len(eids)
        #centroids = np.zeros((neids, 3), dtype='float32')
        ieids = np.searchsorted(self.element_ids, eids)
        #print('ieids = ', ieids)

        for cell_id in ieids:
            centroid = self.cell_centroid(cell_id)
            unused_result_name, result_values, unused_xyz = self.get_result_by_cell_id(
                cell_id, centroid, icase_result)
            texti = '%s' % result_values
            xi, yi, zi = centroid
            self._create_annotation(texti, self.label_actors[icase_to_apply], xi, yi, zi)
        self.log_command('mark_elements_by_different_case(%s, %s, %s)' % (
            eids, icase_result, icase_to_apply))
        self.vtk_interactor.Render()

    def mark_nodes(self, nids, icase, text):
        """
        Marks a series of nodes with custom text labels

        Parameters
        ----------
        nids : int, List[int]
            the nodes to apply a message to
        icase : int
            the key in label_actors to slot the result into
        text : str, List[str]
            the text to display

        0 corresponds to the NodeID result
        self.mark_nodes(1, 0, 'max')
        self.mark_nodes(6, 0, 'min')
        self.mark_nodes([1, 6], 0, 'max')
        self.mark_nodes([1, 6], 0, ['max', 'min'])
        """
        if icase not in self.label_actors:
            msg = 'icase=%r not in label_actors=[%s]' % (
                icase, ', '.join(self.label_actors))
            self.log_error(msg)
            return
        i = np.searchsorted(self.node_ids, nids)
        if isinstance(text, string_types):
            text = [text] * len(i)
        else:
            assert len(text) == len(i)

        xyz = self.xyz_cid0[i, :]
        for (xi, yi, zi), texti in zip(xyz, text):
            self._create_annotation(texti, self.label_actors[icase], xi, yi, zi)
        self.vtk_interactor.Render()

    def __mark_nodes_by_result(self, nids, icases):
        """
        # mark the node 1 with the NodeID (0) result
        self.mark_nodes_by_result_case(1, 0)

        # mark the nodes 1 and 2 with the NodeID (0) result
        self.mark_nodes_by_result_case([1, 2], 0)

        # mark the nodes with the NodeID (0) and ElementID (1) result
        self.mark_nodes_by_result_case([1, 2], [0, 1])
        """
        i = np.searchsorted(self.node_ids, nids)
        if isinstance(icases, int):
            icases = [icases]

        for icase in icases:
            if icase not in self.label_actors:
                msg = 'icase=%r not in label_actors=[%s]' % (
                    icase, ', '.join(self.label_actors))
                self.log_error(msg)
                continue

            for node_id in i:
                jnid = np.where(node_id == self.node_ids)[0]
                world_position = self.xyz_cid0[jnid, :]
                out = self.get_result_by_xyz_node_id(world_position, node_id)
                _result_name, unused_result_value, node_id, node_xyz = out
                xi, yi, zi = node_xyz
                texti = 'test'
                self._create_annotation(texti, self.label_actors[icase], xi, yi, zi)
        self.vtk_interactor.Render()

    def init_cell_picker(self):
        self.is_pick = False
        if not self.run_vtk:
            return
        self.vtk_interactor.SetPicker(self.node_picker)
        self.vtk_interactor.SetPicker(self.cell_picker)
        self.mouse_actions.setup_mouse_buttons(mode='probe_result')
        self.mouse_actions.setup_mouse_buttons(mode='default')

    def convert_units(self, unused_result_name, result_value, xyz):
        #self.input_units
        #self.display_units
        return result_value, xyz

    def _create_annotation(self, text, slot, x, y, z):
        """
        Creates the actual annotation and appends it to slot

        Parameters
        ----------
        text : str
            the text to display
        x, y, z : float
            the position of the label
        slot : List[annotation]
            where to place the annotation
            self.label_actors[icase] : List[annotation]
                icase : icase
                    the key in label_actors to slot the result into
                annotation : vtkBillboardTextActor3D
                    the annotation object
            ???
        """
        if not isinstance(slot, list):
            msg = 'slot=%r type=%s' % (slot, type(slot))
            raise TypeError(msg)
        # http://nullege.com/codes/show/src%40p%40y%40pymatgen-2.9.6%40pymatgen%40vis%40structure_vtk.py/395/vtk.vtkVectorText/python

        #self.convert_units(icase, result_value, x, y, z)

        text_actor = vtk.vtkBillboardTextActor3D()
        label = text
        text_actor.SetPosition(x, y, z)
        text_actor.SetInput(label)
        text_actor.PickableOff()
        text_actor.DragableOff()
        #text_actor.SetPickable(False)

        #text_actor.SetPosition(actor.GetPosition())
        text_prop = text_actor.GetTextProperty()
        text_prop.SetFontSize(self.settings.annotation_size)
        text_prop.SetFontFamilyToArial()
        text_prop.BoldOn()
        text_prop.ShadowOn()

        text_prop.SetColor(self.settings.annotation_color)
        text_prop.SetJustificationToCentered()

        # finish adding the actor
        self.rend.AddActor(text_actor)

        #self.label_actors[icase].append(text_actor)
        slot.append(text_actor)

        #print('added label actor %r; icase=%s' % (text, icase))
        #print(self.label_actors)

        #self.picker_textMapper.SetInput("(%.6f, %.6f, %.6f)"% pickPos)
        #camera.GetPosition()
        #camera.GetClippingRange()
        #camera.GetFocalPoint()

    def _on_multi_pick(self, unused_a):
        """
        vtkFrustumExtractor
        vtkAreaPicker
        """
        pass

    def _on_cell_picker(self, unused_a):
        self.vtk_interactor.SetPicker(self.cell_picker)
        picker = self.cell_picker
        world_position = picker.GetPickPosition()
        cell_id = picker.GetCellId()
        select_point = picker.GetSelectionPoint()  # get x,y pixel coordinate

        self.log_info("world_position = %s" % str(world_position))
        self.log_info("cell_id = %s" % cell_id)
        self.log_info("select_point = %s" % str(select_point))

    def _on_node_picker(self, unused_a):
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

    def show_only(self, names):
        """
        Show these actors only

        names : str, List[str]
            names to show
            If they're hidden, show them.
            If they're shown and shouldn't be, hide them.

        ..todo :: update the GeomeryProperties
        """
        raise NotImplementedError('show_only')

    def hide_actors(self, except_names=None):
        """
        Hide all the actors

        except_names : str, List[str], None
            list of names to exclude
            None : hide all

        ..note :: If an actor is hidden and in the except_names, it will still be hidden.
        ..todo :: update the GeomeryProperties
        """
        if except_names is None:
            except_names = []
        elif isinstance(except_names, string_types):
            except_names = [except_names]

        # hide everything but the main grid
        for key, actor in iteritems(self.geometry_actors):
            if key not in except_names:
                actor.VisibilityOff()

        self.hide_axes()
        self.hide_legend()
        #self.settings.set_background_color_to_white()

    def hide_axes(self, cids=None):
        """
        ..todo :: support cids
        ..todo :: fix the coords
        """
        for axis in itervalues(self.axes):
            axis.VisibilityOff()
        self.corner_axis.EnabledOff()

    def show_axes(self, cids=None):
        """
        ..todo :: support cids
        ..todo :: fix the coords
        """
        for axis in itervalues(self.axes):
            axis.VisibilityOn()
        self.corner_axis.EnabledOn()

    def make_gif(self, gif_filename, scale, istep=None,
                 min_value=None, max_value=None,
                 animate_scale=True, animate_phase=False, animate_time=False,
                 icase=None, icase_start=None, icase_end=None, icase_delta=None,
                 time=2.0, animation_profile='0 to scale',
                 nrepeat=0, fps=30, magnify=1,
                 make_images=True, delete_images=False, make_gif=True, stop_animation=False,
                 animate_in_gui=True):
        """
        Makes an animated gif

        Parameters
        ----------
        gif_filename : str
            path to the output gif & png folder
        scale : float
            the deflection scale factor; true scale
        istep : int; default=None
            the png file number (let's you pick a subset of images)
            useful for when you press ``Step``
        stop_animation : bool; default=False
            stops the animation; don't make any images/gif
        animate_in_gui : bool; default=True
            animates the model; don't make any images/gif
            stop_animation overrides animate_in_gui
            animate_in_gui overrides make_gif

        Pick One
        --------
        animate_scale : bool; default=True
            does a deflection plot (single subcase)
        animate_phase : bool; default=False
            does a complex deflection plot (single subcase)
        animate_time : bool; default=False
            does a deflection plot (multiple subcases)

        Other
        -----
        istep : int
            the png file number (let's you pick a subset of images)
            useful for when you press ``Step``
        time : float; default=2.0
            the runtime of the gif (seconds)
        fps : int; default=30
            the frames/second

        Case Selection
        --------------
        icase : int; default=None
            None : unused
            int : the result case to plot the deflection for
                  active if animate_scale=True or animate_phase=True
        icase_start : int; default=None
            starting case id
            None : unused
            int : active if animate_time=True
        icase_end : int; default=None
            starting case id
            None : unused
            int : active if animate_time=True
        icase_delta : int; default=None
            step size
            None : unused
            int : active if animate_time=True

        Time Plot Options
        -----------------
        max_value : float; default=None
            the max value on the plot
        min_value : float; default=None
            the min value on the plot

        Options
        -------
        animation_profile : str; default='0 to scale'
            animation profile to follow
                '0 to Scale',
                '0 to Scale to 0',
                #'0 to Scale to -Scale to 0',
                '-Scale to Scale',
                '-scale to scale to -scale',
        nrepeat : int; default=0
            0 : loop infinitely
            1 : loop 1 time
            2 : loop 2 times

        Final Control Options
        ---------------------
        make_images : bool; default=True
            make the images
        delete_images : bool; default=False
            cleanup the png files at the end
        make_gif : bool; default=True
            actually make the gif at the end

        Other local variables
        ---------------------
        duration : float
           frame time (seconds)

        For one sided data
        ------------------
         - scales/phases should be one-sided
         - time should be one-sided
         - analysis_time should be one-sided
         - set onesided=True

        For two-sided data
        ------------------
         - scales/phases should be one-sided
         - time should be two-sided
         - analysis_time should be one-sided
         - set onesided=False
        """
        if stop_animation:
            return self.stop_animation()

        phases, icases, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=istep,
            animate_scale=animate_scale, animate_phase=animate_phase, animate_time=animate_time,
            icase=icase,
            icase_start=icase_start, icase_end=icase_end, icase_delta=icase_delta,
            time=time, animation_profile=animation_profile,
            fps=fps)

        parent = self

        #animate_in_gui = True
        self.stop_animation()
        if len(icases) == 1:
            pass
        elif animate_in_gui:
            class vtkAnimationCallback(object):
                """
                http://www.vtk.org/Wiki/VTK/Examples/Python/Animation
                """
                def __init__(self):
                    self.timer_count = 0
                    self.cycler = cycle(range(len(icases)))

                    self.icase0 = -1
                    self.ncases = len(icases)

                def execute(self, obj, unused_event):
                    unused_iren = obj
                    i = self.timer_count % self.ncases
                    #j = next(self.cycler)
                    unused_istep = isteps[i]
                    icase = icases[i]
                    scale = scales[i]
                    phase = phases[i]
                    if icase != self.icase0:
                        #self.cycle_results(case=icase)
                        parent.cycle_results_explicit(icase, explicit=True,
                                                      min_value=min_value, max_value=max_value)
                    try:
                        parent.update_grid_by_icase_scale_phase(icase, scale, phase=phase)
                    except AttributeError:
                        parent.log_error('Invalid Case %i' % icase)
                        parent.stop_animation()
                    self.icase0 = icase

                    parent.vtk_interactor.Render()
                    self.timer_count += 1

            # Sign up to receive TimerEvent
            callback = vtkAnimationCallback()

            observer_name = self.vtk_interactor.AddObserver('TimerEvent', callback.execute)
            self.observers['TimerEvent'] = observer_name

            # total_time not needed
            # fps
            # -> frames_per_second = 1/fps
            delay = int(1. / fps * 1000)

            # time in milliseconds
            unused_timer_id = self.vtk_interactor.CreateRepeatingTimer(delay)
            return

        is_failed = True
        try:
            is_failed = self.make_gif_helper(
                gif_filename, icases, scales,
                phases=phases, isteps=isteps,
                max_value=max_value, min_value=min_value,
                time=time, analysis_time=analysis_time, fps=fps, magnify=magnify,
                onesided=onesided, nrepeat=nrepeat,
                make_images=make_images, delete_images=delete_images, make_gif=make_gif)
        except Exception as error:
            self.log_error(str(error))
            raise
            #self.log_error(traceback.print_stack(f))
            #self.log_error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)

        if not is_failed:
            msg = (
                'make_gif(%r, %s, istep=%s,\n'
                '         min_value=%s, max_value=%s,\n'
                '         animate_scale=%s, animate_phase=%s, animate_time=%s,\n'
                '         icase=%s, icase_start=%s, icase_end=%s, icase_delta=%s,\n'
                "         time=%s, animation_profile=%r,\n"
                '         nrepeat=%s, fps=%s, magnify=%s,\n'
                '         make_images=%s, delete_images=%s, make_gif=%s, stop_animation=%s,\n'
                '         animate_in_gui=%s)\n' % (
                    gif_filename, scale, istep, min_value, max_value,
                    animate_scale, animate_phase, animate_time,
                    icase, icase_start, icase_end, icase_delta, time, animation_profile,
                    nrepeat, fps, magnify, make_images, delete_images, make_gif, stop_animation,
                    animate_in_gui)
            )
            self.log_command(msg)

        return is_failed

    def stop_animation(self):
        """removes the animation timer"""
        is_failed = False
        if 'TimerEvent' in self.observers:
            observer_name = self.observers['TimerEvent']
            self.vtk_interactor.RemoveObserver(observer_name)
            del self.observers['TimerEvent']
            self.mouse_actions.setup_mouse_buttons(mode='default', force=True)
        return is_failed

    def make_gif_helper(self, gif_filename, icases, scales, phases=None, isteps=None,
                        max_value=None, min_value=None,
                        time=2.0, analysis_time=2.0, fps=30, magnify=1,
                        onesided=True, nrepeat=0,
                        make_images=True, delete_images=False, make_gif=True):
        """
        Makes an animated gif

        Parameters
        ----------
        gif_filename : str
            path to the output gif & png folder
        icases : int / List[int]
            the result case to plot the deflection for
        scales : List[float]
            List[float] : the deflection scale factors; true scale
        phases : List[float]; default=None
            List[float] : the phase angles (degrees)
            None -> animate scale
        max_value : float; default=None
            the max value on the plot (not supported)
        min_value : float; default=None
            the min value on the plot (not supported)
        isteps : List[int]
            the png file numbers (let's you pick a subset of images)
            useful for when you press ``Step``
        time : float; default=2.0
            the runtime of the gif (seconds)
        analysis_time : float; default=2.0
            The time we actually need to simulate (seconds).
            We don't need to take extra pictures if they're just copies.
        fps : int; default=30
            the frames/second

        Options
        -------
        onesided : bool; default=True
            should the animation go up and back down
            True : the video will use images [0...N]
            False : the video will use images [0...N...0]
        nrepeat : int; default=0
            0 : loop infinitely
            1 : loop 1 time
            2 : loop 2 times

        Final Control Options
        ---------------------
        make_images : bool; default=True
            make the images
        delete_images : bool; default=False
            cleanup the png files at the end
        make_gif : bool; default=True
            actually make the gif at the end

        Other local variables
        ---------------------
        duration : float
           frame time (seconds)

        For one sided data
        ------------------
         - scales/phases should be one-sided
         - time should be one-sided
         - analysis_time should be one-sided
         - set onesided=True

        For two-sided data
        ------------------
         - scales/phases should be one-sided
         - time should be two-sided
         - analysis_time should be one-sided
         - set onesided=False
        """
        assert fps >= 1, fps
        nframes = ceil(analysis_time * fps)
        assert nframes >= 2, nframes
        unused_duration = time / nframes
        nframes = int(nframes)

        png_dirname = os.path.dirname(os.path.abspath(gif_filename))
        if not os.path.exists(png_dirname):
            os.makedirs(png_dirname)

        phases, icases, isteps, scales = update_animation_inputs(
            phases, icases, isteps, scales, analysis_time, fps)

        if gif_filename is not None:
            png_filenames = []
            fmt = gif_filename[:-4] + '_%%0%ii.png' % (len(str(nframes)))

        icase0 = -1
        is_failed = True
        if make_images:
            for istep, icase, scale, phase in zip(isteps, icases, scales, phases):
                if icase != icase0:
                    #self.cycle_results(case=icase)
                    self.cycle_results_explicit(icase, explicit=True,
                                                min_value=min_value, max_value=max_value)
                self.update_grid_by_icase_scale_phase(icase, scale, phase=phase)

                if gif_filename is not None:
                    png_filename = fmt % istep
                    self.on_take_screenshot(fname=png_filename, magnify=magnify)
                    png_filenames.append(png_filename)
        else:
            for istep in isteps:
                png_filename = fmt % istep
                png_filenames.append(png_filename)
                assert os.path.exists(png_filename), 'png_filename=%s' % png_filename

        if gif_filename is not None and png_filenames:
            is_failed = write_gif(
                gif_filename, png_filenames, time=time,
                onesided=onesided,
                nrepeat=nrepeat, delete_images=delete_images,
                make_gif=make_gif)
        return is_failed

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
            #from vtk.numpy_interface import algorithms
            #arrow = vtk.vtkArrowSource()
            #arrow.PickableOff()

            #self.glyph_transform = vtk.vtkTransform()
            #self.glyph_transform_filter = vtk.vtkTransformPolyDataFilter()
            #self.glyph_transform_filter.SetInputConnection(arrow.GetOutputPort())
            #self.glyph_transform_filter.SetTransform(self.glyph_transform)

            #self.glyph = vtk.vtkGlyph3D()
            #self.glyph.setInput(xxx)
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
        self.build_glyph()

    def build_glyph(self):
        """builds the glyph actor"""
        grid = self.grid
        glyphs = vtk.vtkGlyph3D()
        #if filter_small_forces:
            #glyphs.SetRange(0.5, 1.)

        glyphs.SetVectorModeToUseVector()
        #apply_color_to_glyph = False
        #if apply_color_to_glyph:
        #glyphs.SetScaleModeToScaleByScalar()
        glyphs.SetScaleModeToScaleByVector()
        glyphs.SetColorModeToColorByScale()
        #glyphs.SetColorModeToColorByScalar()  # super tiny
        #glyphs.SetColorModeToColorByVector()  # super tiny

        glyphs.ScalingOn()
        glyphs.ClampingOn()
        #glyphs.Update()

        glyph_source = vtk.vtkArrowSource()
        #glyph_source.InvertOn()  # flip this arrow direction
        glyphs.SetInputData(grid)


        glyphs.SetSourceConnection(glyph_source.GetOutputPort())
        #glyphs.SetScaleModeToDataScalingOff()
        #glyphs.SetScaleFactor(10.0)  # bwb
        #glyphs.SetScaleFactor(1.0)  # solid-bending
        glyph_mapper = vtk.vtkPolyDataMapper()
        glyph_mapper.SetInputConnection(glyphs.GetOutputPort())
        glyph_mapper.ScalarVisibilityOff()

        arrow_actor = vtk.vtkLODActor()
        arrow_actor.SetMapper(glyph_mapper)

        prop = arrow_actor.GetProperty()
        prop.SetColor(1., 0., 0.)
        self.rend.AddActor(arrow_actor)
        #self.grid.GetPointData().SetActiveVectors(None)
        arrow_actor.SetVisibility(False)

        self.glyph_source = glyph_source
        self.glyphs = glyphs
        self.glyph_mapper = glyph_mapper
        self.arrow_actor = arrow_actor
        #-----------------------------------------
        glyphs_centroid = vtk.vtkGlyph3D()
        glyphs_centroid.SetVectorModeToUseVector()
        glyphs_centroid.SetScaleModeToScaleByVector()
        glyphs_centroid.SetColorModeToColorByScale()
        glyphs_centroid.ScalingOn()
        glyphs_centroid.ClampingOn()
        glyphs_centroid.SetSourceConnection(glyph_source.GetOutputPort())

        glyph_mapper_centroid = vtk.vtkPolyDataMapper()
        glyph_mapper_centroid.SetInputConnection(glyphs_centroid.GetOutputPort())
        glyph_mapper_centroid.ScalarVisibilityOff()

        arrow_actor_centroid = vtk.vtkLODActor()
        arrow_actor_centroid.SetMapper(glyph_mapper_centroid)

        self.glyphs_centroid = glyphs_centroid
        self.glyph_mapper_centroid = glyph_mapper_centroid
        self.arrow_actor_centroid = arrow_actor_centroid

    def ResetCamera(self):
        self.GetCamera().ResetCamera()

    def GetCamera(self):
        return self.rend.GetActiveCamera()

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
        self.icase_disp = None
        self.icase_vector = None
        self.icase_fringe = None
        self.set_form(form)

    def _finish_results_io2(self, form, cases, reset_labels=True):
        """
        Adds results to the Sidebar

        Parameters
        ----------
        form : List[pairs]
            There are two types of pairs
            header_pair : (str, None, List[pair])
                defines a heading
                str : the sidebar label
                None : flag that there are sub-results
                List[pair] : more header/result pairs
            result_pair : (str, int, List[])
                str : the sidebar label
                int : the case id
                List[] : flag that there are no sub-results
        cases : dict[case_id] = result
            case_id : int
                the case id
            result : GuiResult
                the class that stores the result
        reset_labels : bool; default=True
            should the label actors be reset

        form = [
            'Model', None, [
                ['NodeID', 0, []],
                ['ElementID', 1, []]
                ['PropertyID', 2, []]
            ],
            'time=0.0', None, [
                ['Stress', 3, []],
                ['Displacement', 4, []]
            ],
            'time=1.0', None, [
                ['Stress', 5, []],
                ['Displacement', 6, []]
            ],
        ]
        cases = {
            0 : GuiResult(...),  # NodeID
            1 : GuiResult(...),  # ElementID
            2 : GuiResult(...),  # PropertyID
            3 : GuiResult(...),  # Stress; t=0.0
            4 : GuiResult(...),  # Displacement; t=0.0
            5 : GuiResult(...),  # Stress; t=1.0
            6 : GuiResult(...),  # Displacement; t=1.0
        }
        case_keys = [0, 1, 2, 3, 4, 5, 6]
        """
        self.turn_text_on()
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

        self.reset_labels(reset_minus1=reset_labels)
        self.cycle_results_explicit()  # start at nCase=0
        if self.ncases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

        #data = [
        #    ('A', []),
        #    ('B', []),
        #    ('C', []),
        #]

        data = []
        for key in self.case_keys:
            assert isinstance(key, integer_types), key
            unused_obj, (i, unused_name) = self.result_cases[key]
            tuple_data = (i, [])
            data.append(tuple_data)

        self.res_widget.update_results(form, self.name)

        key = self.case_keys[0]
        location = self.get_case_location(key)
        method = 'centroid' if location else 'nodal'

        data2 = [(method, None, [])]
        self.res_widget.update_methods(data2)

        if self.is_groups:
            if self.element_ids is None:
                raise RuntimeError('implement self.element_ids for this format')
            #eids = np.arange(172)
            #eids = []
            #self.hide_elements_mask(eids)
            elements_pound = self.element_ids[-1]
            main_group = Group(
                'main', '', elements_pound,
                editable=False)
            main_group.element_ids = self.element_ids
            self.groups['main'] = main_group
            self.post_group(main_group)
            #self.show_elements_mask(np.arange(self.nelements))

        for unused_module_name, module in iteritems(self.modules):
            module.post_load_geometry()

    def get_result_by_cell_id(self, cell_id, world_position, icase=None):
        """should handle multiple cell_ids"""
        if icase is None:
            icase = self.icase
        case_key = self.case_keys[icase] # int for object
        case = self.result_cases[case_key]

        (obj, (i, res_name)) = case
        unused_subcase_id = obj.subcase_id
        case = obj.get_result(i, res_name)

        try:
            result_values = case[cell_id]
        except IndexError:
            msg = ('case[cell_id] is out of bounds; length=%s\n'
                   'result_name=%r cell_id=%r case_key=%r\n' % (
                       len(case), res_name, cell_id, case_key))
            raise IndexError(msg)

        cell = self.grid_selected.GetCell(cell_id)
        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()
        cell_type = cell.GetCellType()

        if cell_type in [5, 9, 22, 23, 28]:  # CTRIA3, CQUAD4, CTRIA6, CQUAD8, CQUAD
            node_xyz = np.zeros((nnodes, 3), dtype='float32')
            for ipoint in range(nnodes):
                point = points.GetPoint(ipoint)
                node_xyz[ipoint, :] = point
            xyz = node_xyz.mean(axis=0)
        elif cell_type in [10, 12, 13, 14]: # CTETRA4, CHEXA8, CPENTA6, CPYRAM5
            # TODO: No idea how to get the center of the face
            #       vs. a point on a face that's not exposed
            #faces = cell.GetFaces()
            #nfaces = cell.GetNumberOfFaces()
            #for iface in range(nfaces):
                #face = cell.GetFace(iface)
                #points = face.GetPoints()
            #faces
            xyz = world_position
        elif cell_type in [24, 25, 26, 27]: # CTETRA10, CHEXA20, CPENTA15, CPYRAM13
            xyz = world_position
        elif cell_type in [3]: # CBAR, CBEAM, CELASx, CDAMPx, CBUSHx
            node_xyz = np.zeros((nnodes, 3), dtype='float32')
            for ipoint in range(nnodes):
                point = points.GetPoint(ipoint)
                node_xyz[ipoint, :] = point
            xyz = node_xyz.mean(axis=0)
        elif cell_type in [21]: # CBEND
            # 21-QuadraticEdge
            node_xyz = np.zeros((nnodes, 3), dtype='float32')
            for ipoint in range(nnodes):
                point = points.GetPoint(ipoint)
                node_xyz[ipoint, :] = point
            xyz = node_xyz.mean(axis=0)
        else:
            #self.log.error(msg)
            msg = 'cell_type=%s nnodes=%s; icase=%s result_values=%s' % (
                cell_type, nnodes, icase, result_values)
            self.log.error(msg)
            #VTK_LINE = 3

            #VTK_TRIANGLE = 5
            #VTK_QUADRATIC_TRIANGLE = 22

            #VTK_QUAD = 9
            #VTK_QUADRATIC_QUAD = 23

            #VTK_TETRA = 10
            #VTK_QUADRATIC_TETRA = 24

            #VTK_WEDGE = 13
            #VTK_QUADRATIC_WEDGE = 26

            #VTK_HEXAHEDRON = 12
            #VTK_QUADRATIC_HEXAHEDRON = 25

            #VTK_PYRAMID = 14
            #VTK_QUADRATIC_PYRAMID = 27
            raise NotImplementedError(msg)
        return res_name, result_values, xyz

    def cell_centroid(self, cell_id):
        """gets the cell centroid"""
        cell = self.grid_selected.GetCell(cell_id)
        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()
        centroid = np.zeros(3, dtype='float32')
        for ipoint in range(nnodes):
            point = np.array(points.GetPoint(ipoint), dtype='float32')
            centroid += point
        centroid /= nnodes
        return centroid

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
        assert isinstance(case_key, integer_types), case_key
        (obj, (i, res_name)) = case
        unused_subcase_id = obj.subcase_id
        case = obj.get_result(i, res_name)
        result_values = case[node_id]
        assert not isinstance(xyz, int), xyz
        return result_name, result_values, node_id, xyz

    @property
    def result_name(self):
        """
        creates the self.result_name variable
        """
        # case_key = (1, 'ElementID', 1, 'centroid', '%.0f')
        case_key = self.case_keys[self.icase]
        assert isinstance(case_key, integer_types), case_key
        unused_obj, (unused_i, res_name) = self.result_cases[case_key]
        return res_name

    #def finish_io(self, cases):
        #self.result_cases = cases
        #self.case_keys = sorted(cases.keys())
        ##print("case_keys = ", self.case_keys)

        #if len(self.result_cases) == 0:
            #self.ncases = 1
            #self.icase = 0
        #elif len(self.result_cases) == 1:
            #self.ncases = 1
            #self.icase = 0
        #else:
            #self.ncases = len(self.result_cases) - 1  # number of keys in dictionary
            #self.icase = -1

        #self.icase_disp = None
        #self.icase_vector = None
        #self.icase_fringe = None
        #self.cycle_results()  # start at nCase=0

        #if self.ncases:
            #self.scalarBar.VisibilityOn()
            #self.scalarBar.Modified()

    def clear_application_log(self, force=False):
        """
        Clears the application log

        Parameters
        ----------
        force : bool; default=False
            clears the dialog without asking
        """
        # popup menu
        if force:
            self.log_widget.clear()
            self.log_command('clear_application_log(force=%s)' % force)
        else:
            widget = QWidget()
            title = 'Clear Application Log'
            msg = 'Are you sure you want to clear the Application Log?'
            result = QMessageBox.question(widget, title, msg,
                                          QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if result == QMessageBox.Yes:
                self.log_widget.clear()
                self.log_command('clear_application_log(force=%s)' % force)

    def delete_actor(self, name):
        """deletes an actor and associated properties"""
        if name != 'main':
            if name in self.geometry_actors:
                actor = self.geometry_actors[name]
                self.rend.RemoveActor(actor)
                del self.geometry_actors[name]

            if name in self.geometry_properties:
                unused_prop = self.geometry_properties[name]
                del self.geometry_properties[name]
            self.Render()

    #---------------------------------------------------------------------------------------
    # CAMERA MENU
    def view_camera(self):
        set_camera_menu(self)

    #def _apply_camera(self, data):
        #name = data['name']
        #self.cameras = deepcopy(data['cameras'])
        #self.on_set_camera(name)

    #---------------------------------------------------------------------------------------
    # PICKER
    @property
    def node_picker_size(self):
        """Gets the node picker size"""
        return self.node_picker.GetTolerance()

    @node_picker_size.setter
    def node_picker_size(self, size):
        """Sets the node picker size"""
        assert size >= 0., size
        self.node_picker.SetTolerance(size)

    @property
    def element_picker_size(self):
        """Gets the element picker size"""
        return self.cell_picker.GetTolerance()

    @element_picker_size.setter
    def element_picker_size(self, size):
        """Sets the element picker size"""
        assert size >= 0., size
        self.cell_picker.SetTolerance(size)

    #---------------------------------------------------------------------------------------
    def set_preferences_menu(self):
        """
        Opens a dialog box to set:

        +--------+----------+
        |  Min   |  Float   |
        +--------+----------+
        """
        set_preferences_menu(self)

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
        set_clipping_menu(self)

    def _apply_clipping(self, data):
        min_clip = data['clipping_min']
        max_clip = data['clipping_max']
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

    def on_set_anti_aliasing(self, scale=0):
        assert isinstance(scale, int), 'scale=%r; type=%r' % (scale, type(scale))
        renwin = self.render_window
        renwin.LineSmoothingOn()
        renwin.PolygonSmoothingOn()
        renwin.PointSmoothingOn()
        renwin.SetMultiSamples(scale)
        self.vtk_interactor.Render()
        self.log_command('on_set_anti_aliasing(%r)' % (scale))

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
        set_legend_menu(self)

    def update_legend(self, icase_fringe, icase_disp, icase_vector,
                      name, min_value, max_value, data_format, scale, phase,
                      arrow_scale,
                      nlabels, labelsize, ncolors, colormap,
                      use_fringe_internal=False, use_disp_internal=False,
                      use_vector_internal=False, external_call=True):
        """
        Internal method for updating the legend

        Parameters
        ----------
        use_fringe_internal : bool; default=Falsee
            True : use the internal fringe parameters
            False : use the values that were passed in
        use_disp_internal : bool; default=Falsee
            True : use the internal of scale and phase
            False : use the values that were passed in
        use_vector_internal : bool; default=Falsee
            True : use the internal of arrow_scale
            False : use the values that were passed in
        external_call : bool; default=True
            True : allow the legend ``on_apply`` method to be called
            False : the scalar bar/displacement updating will be handled
                    manually to prevent recursion (and a crash)
        """
        if not self._legend_window_shown:
            return
        self._legend_window._updated_legend = True
        is_fringe = self._is_fringe

        out = get_legend_fringe(self, icase_fringe)
        (
            _result_type, scalar_bar, defaults_scalar_bar, data_format,
            default_format, default_title, _min_value, _max_value,
            default_min, default_max) = out

        unused_nlabels, _labelsize, _ncolors, _colormap = scalar_bar
        default_nlabels, default_labelsize, default_ncolors, default_colormap = defaults_scalar_bar
        if use_fringe_internal:
            min_value = _min_value
            max_value = _max_value
            unused_result_type = _result_type
            labelsize = _labelsize
            ncolors = _ncolors
            colormap = _colormap

        #if icase_fringe is not None:
            #key = self.case_keys[icase_fringe]
            #assert isinstance(key, integer_types), key
            #(obj, (i, name)) = self.result_cases[key]
            ##subcase_id = obj.subcase_id
            ##case = obj.get_result(i, name)
            ##result_type = obj.get_title(i, name)
            ##vector_size = obj.get_vector_size(i, name)
            ##location = obj.get_location(i, name)
            ##data_format = obj.get_data_format(i, name)
            ##scale = obj.get_scale(i, name)
            ##label2 = obj.get_header(i, name)
            #default_data_format = obj.get_default_data_format(i, name)
            #default_min, default_max = obj.get_default_min_max(i, name)
            #default_title = obj.get_default_title(i, name)
            #out_labels = obj.get_default_nlabels_labelsize_ncolors_colormap(i, name)
            #default_nlabels, default_labelsize, default_ncolors, default_colormap = out_labels
            #is_normals = obj.is_normal_result(i, name)
            #is_fringe = not is_normals

        _scale, _phase, default_scale, default_phase = get_legend_disp(
            self, icase_disp)
        #if icase_disp is not None:
            #default_scale = obj.get_default_scale(i, name)
            #default_phase = obj.get_default_phase(i, name)
        if use_disp_internal:
            scale = _scale
            phase = _phase
            #default_scale = _default_scale
            #default_phase = _default_phase


        _arrow_scale, default_arrow_scale = get_legend_vector(self, icase_vector)
        if use_vector_internal:
            arrow_scale = _arrow_scale
            #default_arrow_scale = _default_arrow_scale

        #assert isinstance(scale, float), 'scale=%s' % scale
        self._legend_window.update_legend(
            icase_fringe, icase_disp, icase_vector,
            name, min_value, max_value, data_format,
            nlabels, labelsize, ncolors, colormap, is_fringe,
            scale, phase,
            arrow_scale,

            default_title, default_min, default_max, default_format,
            default_nlabels, default_labelsize,
            default_ncolors, default_colormap,
            default_scale, default_phase,
            default_arrow_scale,
            font_size=self.settings.font_size)
        #self.scalar_bar.set_visibility(self._legend_shown)
        #self.vtk_interactor.Render()

    def _apply_legend(self, data):
        title = data['name']
        min_value = data['min']
        max_value = data['max']
        scale = data['scale']
        phase = data['phase']
        arrow_scale = data['arrow_scale']
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
                              scale=scale, phase=phase,
                              arrow_scale=arrow_scale,
                              data_format=data_format,
                              is_low_to_high=is_low_to_high,
                              is_discrete=is_discrete, is_horizontal=is_horizontal,
                              nlabels=nlabels, labelsize=labelsize,
                              ncolors=ncolors, colormap=colormap,
                              is_shown=is_shown)

    def on_update_legend(self,
                         title='Title', min_value=0., max_value=1.,
                         scale=0.0, phase=0.0,
                         arrow_scale=1.,
                         data_format='%.0f',
                         is_low_to_high=True, is_discrete=True, is_horizontal=True,
                         nlabels=None, labelsize=None, ncolors=None, colormap=None,
                         is_shown=True):
        """
        Updates the legend/model

        Parameters
        ----------
        scale : float
            displacemnt scale factor; true scale

        TODO: speed up by using existing values to skip update steps
        """
        if colormap is None:
            colormap = self.settings.colormap

        self.is_low_to_high = is_low_to_high
        self.is_horizontal_scalar_bar = is_horizontal

        #print('is_shown2 =', is_shown)
        #assert is_shown == False, is_shown
        is_normal = False
        if self.icase_fringe is not None:
            key = self.case_keys[self.icase_fringe]
            assert isinstance(key, integer_types), key
            (obj, (i, res_name)) = self.result_cases[key]
            subcase_id = obj.subcase_id

            unused_location = obj.get_location(i, res_name)
            min_value_old, max_value_old = obj.get_min_max(i, res_name)
            data_format_old = obj.get_data_format(i, res_name)
            colors_old = obj.get_nlabels_labelsize_ncolors_colormap(i, res_name)
            nlabels_old, labelsize_old, ncolors_old, colormap_old = colors_old

            update_fringe = (
                min_value != min_value_old or
                max_value != max_value_old
            )
            update_legend = (
                (
                    (nlabels, labelsize, ncolors, colormap) !=
                    (nlabels_old, labelsize_old, ncolors_old, colormap_old) or
                    data_format != data_format_old) and
                not update_fringe)

            obj.set_min_max(i, res_name, min_value, max_value)
            obj.set_data_format(i, res_name, data_format)
            obj.set_nlabels_labelsize_ncolors_colormap(
                i, res_name, nlabels, labelsize, ncolors, colormap)

            #data_format = obj.get_data_format(i, res_name)
            #obj.set_format(i, res_name, data_format)
            #obj.set_data_format(i, res_name, data_format)
            unused_subtitle, unused_label = self.get_subtitle_label(subcase_id)
            is_normal = obj.is_normal_result(i, res_name)
            #if scale != scale_old or phase != phase_old:
            #if not from_legend_menu:
            if update_fringe:
                self.on_fringe(self.icase_fringe, show_msg=False,
                               update_legend_window=False)

        if is_normal:
            return

        if self.icase_disp is not None:
            key = self.case_keys[self.icase_disp]
            assert isinstance(key, integer_types), key
            (objd, (i, res_name)) = self.result_cases[key]
            scale_old = objd.get_scale(i, res_name)
            phase_old = objd.get_phase(i, res_name)
            update_disp = scale != scale_old or phase != phase_old
            if update_disp:
                objd.set_scale(i, res_name, scale)
                objd.set_phase(i, res_name, phase)
                assert isinstance(scale, float), scale
                self.on_disp(self.icase_disp, apply_fringe=False,
                             update_legend_window=False, show_msg=False)

        if self.icase_vector is not None:
            key = self.case_keys[self.icase_vector]
            assert isinstance(key, integer_types), key
            (objv, (i, res_name)) = self.result_cases[key]
            arrow_scale_old = objv.get_scale(i, res_name)
            objv.set_scale(i, res_name, arrow_scale)
            assert isinstance(arrow_scale, float), arrow_scale
            update_vector = arrow_scale != arrow_scale_old
            if update_vector:
                self.on_vector(self.icase_vector, apply_fringe=False,
                               update_legend_window=False, show_msg=False)

        #unused_name = (vector_size1, subcase_id, result_type, label, min_value, max_value, scale1)
        #if obj.is_normal_result(i, res_name):
            #return

        if self.icase_fringe is None:
            return

        norm_value = float(max_value - min_value)
        # if name not in self._loaded_names:

        #if isinstance(key, integer_types):  # vector 3
             #norm_plot_value = norm(plot_value, axis=1)
            #grid_result = self.set_grid_values(name, norm_plot_value, vector_size1,
                                               #min_value, max_value, norm_value,
                                               #is_low_to_high=is_low_to_high)
        #else:
        if update_legend:
            self.update_scalar_bar(title, min_value, max_value, norm_value,
                                   data_format,
                                   nlabels=nlabels, labelsize=labelsize,
                                   ncolors=ncolors, colormap=colormap,
                                   is_shown=is_shown)

        msg = ('self.on_update_legend(title=%r, min_value=%s, max_value=%s,\n'
               '                      scale=%r, phase=%r,\n'
               '                      data_format=%r, is_low_to_high=%s, is_discrete=%s,\n'
               '                      nlabels=%r, labelsize=%r, ncolors=%r, colormap=%r,\n'
               '                      is_horizontal=%r, is_shown=%r)'
               % (title, min_value, max_value, scale, phase,
                  data_format, is_low_to_high, is_discrete,
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

        data = deepcopy(self.geometry_properties)
        data['font_size'] = self.settings.font_size
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

            if name not in self.geometry_properties:
                # we've deleted the actor
                continue

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

    def on_update_geometry_properties_override_dialog(self, geometry_properties):
        """
        Update the goemetry properties and overwite the options in the
        edit geometry properties dialog if it is open.

        Parameters
        -----------
        geometry_properties : dict {str : CoordProperties or AltGeometry}
            Dictionary from name to properties object. Only the names included in
            ``geometry_properties`` are modified.
        """
        if self._edit_geometry_properties_window_shown:
            # Override the output state in the edit geometry properties diaglog
            # if the button is pushed while the dialog is open. This prevent the
            # case where you close the dialog and the state reverts back to
            # before you hit the button.
            for name, prop in iteritems(geometry_properties):
                self._edit_geometry_properties.out_data[name] = prop
                if self._edit_geometry_properties.active_key == name:
                    index = self._edit_geometry_properties.table.currentIndex()
                    self._edit_geometry_properties.update_active_key(index)
        self.on_update_geometry_properties(geometry_properties)

    def on_set_modify_groups(self):
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
        on_set_modify_groups(self)

    def _apply_modify_groups(self, data):
        """called by on_set_modify_groups when apply is clicked"""
        self.on_update_modify_groups(data)
        imain = self._modify_groups_window.imain
        name = self._modify_groups_window.keys[imain]
        self.post_group_by_name(name)

    def on_update_modify_groups(self, out_data):
        """
        Applies the changed groups to the different groups if
        something changed.
        """
        #self.groups = out_data
        data = {}
        for unused_group_id, group in sorted(iteritems(out_data)):
            if not isinstance(group, Group):
                continue
            data[group.name] = group
        self.groups = data

    def on_update_geometry_properties(self, out_data, name=None, write_log=True):
        """
        Applies the changed properties to the different actors if
        something changed.

        Note that some of the values are limited.  This prevents
        points/lines from being shrunk to 0 and also the actor being
        actually "hidden" at the same time.  This prevents confusion
        when you try to show the actor and it's not visible.
        """
        lines = []
        if name is None:
            for namei, group in iteritems(out_data):
                if namei in ['clicked_ok', 'clicked_cancel']:
                    continue
                self._update_ith_geometry_properties(namei, group, lines, render=False)
        else:
            group = out_data[name]
            self._update_ith_geometry_properties(name, group, lines, render=False)

        self.vtk_interactor.Render()
        if write_log and lines:
            msg = 'out_data = {\n'
            msg += ''.join(lines)
            msg += '}\n'
            msg += 'self.on_update_geometry_properties(out_data)'
            self.log_command(msg)

    def _update_ith_geometry_properties(self, namei, group, lines, render=True):
        """updates a geometry"""
        if namei not in self.geometry_actors:
            # we've deleted the actor
            return
        actor = self.geometry_actors[namei]

        if isinstance(actor, vtk.vtkActor):
            alt_prop = self.geometry_properties[namei]
            label_actors = alt_prop.label_actors
            lines += self._update_geometry_properties_actor(namei, group, actor, label_actors)
        elif isinstance(actor, vtk.vtkAxesActor):
            changed = False
            is_visible1 = bool(actor.GetVisibility())
            is_visible2 = group.is_visible
            if is_visible1 != is_visible2:
                actor.SetVisibility(is_visible2)
                alt_prop = self.geometry_properties[namei]
                alt_prop.is_visible = is_visible2
                actor.Modified()
                changed = True

            if changed:
                lines.append('    %r : CoordProperties(is_visible=%s),\n' % (
                    namei, is_visible2))
        else:
            raise NotImplementedError(actor)
        if render:
            self.vtk_interactor.Render()

    def _update_geometry_properties_actor(self, name, group, actor, label_actors):
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
        backface_prop = actor.GetBackfaceProperty()

        if name == 'main' and backface_prop is None:
            # don't edit these
            # we're lying about the colors to make sure the
            # colors aren't reset for the Normals
            color1 = prop.GetDiffuseColor()
            color2 = color1
            assert color1[1] <= 1.0, color1
        else:
            color1 = prop.GetDiffuseColor()
            assert color1[1] <= 1.0, color1
            color2 = group.color_float
            #print('line2646 - name=%s color1=%s color2=%s' % (name, str(color1), str(color2)))
            #color2 = group.color

        opacity1 = prop.GetOpacity()
        opacity2 = group.opacity
        opacity2 = max(0.1, opacity2)

        line_width1 = prop.GetLineWidth()
        line_width2 = group.line_width
        line_width2 = max(1, line_width2)

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

        #print('name=%s color1=%s color2=%s' % (name, str(color1), str(color2)))
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
            #if backface_prop is not None:
                #backface_prop.SetOpacity(opacity2)
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
            for label_actor in label_actors:
                label_actor.SetVisibility(is_visible2)
                label_actor.Modified()
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
            #unused_point = points.GetPoint(2*i+1)
            #print(unused_point)
            node = xyz1[i, :] + length_xyz[i] * bar_scale * dxyz[i, :]
            #print(unused_point, node)
            points.SetPoint(2 * i + 1, *node)

        if hasattr(grid, 'Update'):
            #print('update....')
            grid.Update()
        grid.Modified()
        #print('update2...')
