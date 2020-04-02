# coding: utf-8
# pylint: disable=W0201,C0301
import os.path
from collections import OrderedDict
from math import ceil

import numpy as np
from cpylog import SimpleLogger

from pyNastran.gui.qt_version import qt_version

from qtpy import QtCore, QtGui #, API
from qtpy.QtWidgets import (
    QMessageBox, QWidget,
    QMainWindow, QDockWidget, QFrame, QHBoxLayout, QAction, QToolBar, QMenu, QToolButton)

import vtk


import pyNastran
#print('qt_version = %r' % qt_version)

# vtk makes poor choices regarding the selection of a backend and has no way
# to work around it
#from vtk.qt5.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from .qt_files.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from pyNastran.utils import check_path
from pyNastran.utils.numpy_utils import integer_types

#from .qt_files.gui_attributes import IS_MATPLOTLIB, IS_CUTTING_PLANE
from .qt_files.gui_vtk_common import GuiVTKCommon
from .qt_files.scalar_bar import ScalarBar

from .gui_objects.alt_geometry_storage import AltGeometry

from .menus.menus import (
    on_set_modify_groups, Group,
    Sidebar,
    ApplicationLogWidget,
    PythonConsoleWidget)

from .menus.legend.write_gif import (
    setup_animation, update_animation_inputs, write_gif, make_two_sided)
from .utils.vtk.animation_callback import AnimationCallback
from .utils.vtk.base_utils import numpy_to_vtk_idtype

try:
    from cpylog.html_utils import str_to_html
except ImportError:
    import warnings
    warnings.warn('upgrade your cpylog to v1.4')
    from .utils.html_utils import str_to_html


#from pyNastran.gui.menus.multidialog import MultiFileDialog
from pyNastran.gui.formats import CLASS_MAP


# http://pyqt.sourceforge.net/Docs/PyQt5/multiinheritance.html
class GuiCommon(QMainWindow, GuiVTKCommon):
    """this class adds in interactive/menu capability into the GUI"""
    def __init__(self, **kwds):
        """
        fmt_order, html_logging, inputs, parent=None,
        """
        # this will reset the background color/label color if things break
        #super(QMainWindow, self).__init__(self)
        if qt_version == 'pyqt5':
            super(GuiCommon, self).__init__(**kwds)
        elif qt_version == 'pyside2':
            QMainWindow.__init__(self)
            GuiVTKCommon.__init__(self, **kwds)
        else:  #: pragma: no cover
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

        self.scalar_bar = ScalarBar(self.legend_obj.is_horizontal_scalar_bar)

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

        self.res_widget = Sidebar(
            self,
            include_case_spinner=True,
            include_deflection_scale=False,
            include_vector_scale=False,
        )
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
                'show_warning' : True,
                'show_error' : True,
                'anti_alias_0' : True,
                'anti_alias_1' : False,
                'anti_alias_2' : False,
                'anti_alias_4' : False,
                'anti_alias_8' : False,

                'rotation_center' : False,
                'measure_distance' : False,
                'probe_result' : False,
                'highlight_cell' : False,
                'highlight_node' : False,
                'area_pick' : False,
                'highlight' : False,
                'zoom' : False,
            }

        if tools is None:
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

            tools = file_tools + [
                # labels
                ('label_clear', 'Clear Current Labels', '', 'CTRL+W', 'Clear current labels', self.clear_labels),
                ('label_reset', 'Clear All Labels', '', None, 'Clear all labels', self.reset_labels),

                # view
                ('wireframe', 'Wireframe Model', 'twireframe.png', 'w', 'Show Model as a Wireframe Model', self.on_wireframe),
                ('surface', 'Surface Model', 'tsolid.png', 's', 'Show Model as a Surface Model', self.on_surface),
                ('screenshot', 'Take a Screenshot...', 'tcamera.png', 'CTRL+I', 'Take a Screenshot of current view', self.tool_actions.on_take_screenshot),

                # geometry
                # Geometry:
                #  - Create
                #  - Modify
                ('geometry', 'Geometry', 'geometry.png', None, 'Geometry', self.geometry_obj.show),
                #
                # core menus
                ('legend', 'Modify Legend...', 'legend.png', 'CTRL+L', 'Set Legend', self.legend_obj.set_legend_menu),
                ('animation', 'Create Animation...', 'animation.png', 'CTRL+A', 'Create Animation', self.legend_obj.set_animation_menu),
                ('clipping', 'Set Clipping...', '', None, 'Set Clipping', self.clipping_obj.set_clipping_menu),
                ('set_preferences', 'Preferences...', 'preferences.png', 'CTRL+P', 'Set GUI Preferences', self.preferences_obj.set_preferences_menu),
                ('geo_properties', 'Edit Geometry Properties...', '', 'CTRL+E', 'Change Model Color/Opacity/Line Width', self.edit_geometry_properties_obj.edit_geometry_properties),
                ('map_element_fringe', 'Map Element Fringe', '', 'CTRL+F', 'Map Elemental Centroidal Fringe Result to Nodes', self.map_element_centroid_to_node_fringe_result),

                #('axis', 'Show/Hide Axis', 'axis.png', None, 'Show/Hide Global Axis', self.on_show_hide_axes),

                # groups
                ('modify_groups', 'Modify Groups...', '', None, 'Create/Edit/Delete Groups', self.on_set_modify_groups),
                ('create_groups_by_visible_result', 'Create Groups By Visible Result', '', None, 'Create Groups', self.create_groups_by_visible_result),
                ('create_groups_by_property_id', 'Create Groups By Property ID', '', None, 'Create Groups', self.create_groups_by_property_id),
                #('create_list', 'Create Lists through Booleans', '', None, 'Create List', self.create_list),

                # logging
                ('show_info', 'Show INFO', 'show_info.png', None, 'Show "INFO" messages', self.on_show_info),
                ('show_debug', 'Show DEBUG', 'show_debug.png', None, 'Show "DEBUG" messages', self.on_show_debug),
                ('show_command', 'Show COMMAND', 'show_command.png', None, 'Show "COMMAND" messages', self.on_show_command),
                ('show_warning', 'Show WARNING', 'show_warning.png', None, 'Show "COMMAND" messages', self.on_show_warning),
                ('show_error', 'Show ERROR', 'show_error.png', None, 'Show "COMMAND" messages', self.on_show_error),

                # zoom
                ('magnify', 'Magnify', 'plus_zoom.png', 'm', 'Increase Magnfication', self.on_increase_magnification),
                ('shrink', 'Shrink', 'minus_zoom.png', 'Shift+M', 'Decrease Magnfication', self.on_decrease_magnification),

                # rotation
                ('rotate_clockwise', 'Rotate Clockwise', 'tclock.png', 'o', 'Rotate Clockwise', self.on_rotate_clockwise),
                ('rotate_cclockwise', 'Rotate Counter-Clockwise', 'tcclock.png', 'Shift+O', 'Rotate Counter-Clockwise', self.on_rotate_cclockwise),


                #('cell_pick', 'Cell Pick', '', 'c', 'Centroidal Picking', self.on_cell_picker),
                #('node_pick', 'Node Pick', '', 'n', 'Nodal Picking', self.on_node_picker),

                # help
                ('website', 'Open pyNastran Website...', '', None, 'Open the pyNastran website', self.open_website),
                ('docs', 'Open pyNastran Docs Website...', '', None, 'Open the pyNastran documentation website', self.open_docs),
                ('report_issue', 'Report a Bug/Feature Request...', '', None, 'Open the pyNastran issue tracker', self.open_issue),
                ('discussion_forum', 'Discussion Forum Website...', '', None, 'Open the discussion forum to ask questions', self.open_discussion_forum),
                ('about', 'About pyNastran GUI...', 'tabout.png', 'CTRL+H', 'About pyNastran GUI and help on shortcuts', self.about_dialog),

                # camera
                ('view', 'Camera View', 'view.png', None, 'Load the camera menu', self.camera_obj.set_camera_menu),
                ('camera_reset', 'Reset Camera View', 'trefresh.png', 'r', 'Reset the camera view to default', self.on_reset_camera),

                # results
                ('cycle_results', 'Cycle Results', 'cycle_results.png', 'L', 'Changes the result case', self.on_cycle_results),
                ('rcycle_results', 'Cycle Results', 'rcycle_results.png', 'K', 'Changes the result case', self.on_rcycle_results),

                # view actions
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

                # mouse buttons
                ('rotation_center', 'Set the rotation center', 'trotation_center.png', 'f', 'Pick a node for the rotation center/focal point', self.mouse_actions.on_rotation_center),
                ('measure_distance', 'Measure Distance', 'measure_distance.png', None, 'Measure the distance between two nodes', self.mouse_actions.on_measure_distance),
                ('highlight_cell', 'Highlight Cell', '', None, 'Highlight a single cell', self.mouse_actions.on_highlight_cell),
                ('highlight_node', 'Highlight Node', '', None, 'Highlight a single node', self.mouse_actions.on_highlight_node),
                ('probe_result', 'Probe', 'tprobe.png', None, 'Probe the displayed result', self.mouse_actions.on_probe_result),
                ('quick_probe_result', 'Quick Probe', '', 'p', 'Probe the displayed result', self.mouse_actions.on_quick_probe_result),
                ('zoom', 'Zoom', 'zoom.png', None, 'Zoom In', self.mouse_actions.on_zoom),

                # font size
                ('font_size_increase', 'Increase Font Size', 'text_up.png', 'Ctrl+Plus', 'Increase Font Size', self.on_increase_font_size),
                ('font_size_decrease', 'Decrease Font Size', 'text_down.png', 'Ctrl+Minus', 'Decrease Font Size', self.on_decrease_font_size),

                # picking
                ('area_pick', 'Area Pick', 'tarea_pick.png', None, 'Get a list of nodes/elements', self.mouse_actions.on_area_pick),
                ('highlight', 'Highlight', 'thighlight.png', None, 'Highlight a list of nodes/elements', self.mouse_actions.on_highlight),
                ('highlight_nodes_elements', 'Highlight', 'thighlight.png', None, 'Highlight a list of nodes/elements', self.highlight_obj.set_menu),
                ('mark_nodes_elements', 'Mark', 'tmark.png', None, 'Mark a list of nodes/elements', self.mark_obj.set_menu),
            ]

        if hasattr(self, 'cutting_plane_obj'):
            tools.append(('cutting_plane', 'Cutting Plane...', 'cutting_plane.png', None, 'Create Cutting Plane', self.cutting_plane_obj.set_cutting_plane_menu))

        if 'nastran' in self.fmts:
            tools += [
                ('caero', 'Show/Hide CAERO Panels', '', None, 'Show/Hide CAERO Panel Outlines', self.toggle_caero_panels),
                ('caero_subpanels', 'Toggle CAERO Subpanels', '', None, 'Show/Hide CAERO Subanel Outlines', self.toggle_caero_sub_panels),
                ('conm2', 'Toggle CONM2s', '', None, 'Show/Hide CONM2s', self.toggle_conms),
                ('min', 'Min', '', None, 'Show/Hide Min Label', self.show_hide_min_actor),
                ('max', 'Max', '', None, 'Show/Hide Max Label', self.show_hide_max_actor),
            ]
        self.tools = tools
        self.checkables = checkables

    def keyPressEvent(self, qkey_event):
        #print('qkey_event =', qkey_event.key())
        super(GuiCommon, self).keyPressEvent(qkey_event)

    def _create_menu_bar(self, menu_bar_order=None):
        self.menu_bar_oder = menu_bar_order
        if menu_bar_order is None:
            menu_bar_order = ['menu_file', 'menu_view', 'menu_window', 'menu_help']

        for key in menu_bar_order:
            if key == 'menu_file':
                self.menu_file = self.menubar.addMenu('&File')
            elif key == 'menu_view':
                self.menu_view = self.menubar.addMenu('&View')
            elif key == 'menu_window':
                self.menu_window = self.menubar.addMenu('&Window')
            elif key == 'menu_help':
                self.menu_help = self.menubar.addMenu('&Help')
            elif isinstance(key, tuple):
                attr_name, name = key
                submenu = self.menubar.addMenu(name)
                setattr(self, attr_name, submenu)
            else:
                raise NotImplementedError(key)
        # always last
        self.menu_hidden = self.menubar.addMenu('&Hidden')
        self.menu_hidden.menuAction().setVisible(False)

    def _create_menu_items(self, actions=None, create_menu_bar=True, menu_bar_order=None):
        if actions is None:
            actions = self.actions

        if create_menu_bar:
            self._create_menu_bar(menu_bar_order=menu_bar_order)

        scripts = []
        if self._script_path is not None and os.path.exists(self._script_path):
            scripts = [script for script in os.listdir(self._script_path) if '.py' in script]
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
            'set_preferences', #'cutting_plane',
            '',
            'label_clear', 'label_reset', '',
            'legend', 'animation', 'geo_properties',
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
            menu_view += ['', 'show_info', 'show_debug', 'show_command', 'show_warning', 'show_error']
            menu_window += ['log_dock_widget']
        if self.execute_python:
            self.actions['python_dock_widget'] = self.python_dock_widget.toggleViewAction()
            self.actions['python_dock_widget'].setStatusTip("Show/Hide Python Console")
            menu_window += ['python_dock_widget']

        menu_file = [
            'load_geometry', 'load_results', '',
            'load_custom_result', '',
            'load_csv_user_points', 'load_csv_user_geom', 'script', '', 'exit']
        toolbar_tools = [
            'reload', 'load_geometry', 'load_results',
            'front_view', 'back_view', 'top_view', 'bottom_view',
            'left_view', 'right_view',
            'magnify', 'shrink', 'zoom',
            'rotate_clockwise', 'rotate_cclockwise',
            'rotation_center', 'measure_distance', 'probe_result',
            #'highlight_cell', 'highlight_node',
            'area_pick', 'highlight_nodes_elements', 'mark_nodes_elements',
            'wireframe', 'surface', 'edges']
        toolbar_tools += [
            'camera_reset', 'view', 'screenshot', 'min', 'max', 'map_element_fringe',
            '', # 'exit'
        ]
        hidden_tools = ('cycle_results', 'rcycle_results',
                        'font_size_increase', 'font_size_decrease', 'highlight')

        menu_items = OrderedDict()
        if create_menu_bar:
            menu_items['file'] = (self.menu_file, menu_file)
            menu_items['view'] = (self.menu_view, menu_view)
            menu_items['main'] = (self.menu_window, menu_window)
            menu_items['help'] = (self.menu_help, ('website', 'docs', 'report_issue', 'discussion_forum', 'about',))
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
        action_names = list(self.actions.keys())
        action_names.sort()

        #print("self.actions =", action_names)
        #for plugin in self.plugins:

        menu_items = self._create_menu_items(actions)
        self._populate_menu(menu_items)

        self.actions['show_info'].setChecked(self.settings.show_info)
        self.actions['show_debug'].setChecked(self.settings.show_debug)
        self.actions['show_command'].setChecked(self.settings.show_command)
        self.actions['show_warning'].setChecked(self.settings.show_warning)
        self.actions['show_error'].setChecked(self.settings.show_error)


    def _populate_menu(self, menu_items, actions=None):
        """populate menus and toolbar"""
        assert isinstance(menu_items, dict), menu_items

        if actions is None:
            actions = self.actions
        for unused_menu_name, (menu, items) in menu_items.items():
            if menu is None:
                continue

            for item in items:
                if not item:
                    menu.addSeparator()
                else:
                    if isinstance(item, list):
                        unused_sub_menu_name = item[0]

                        if isinstance(menu, QToolBar):
                            populate_sub_qtoolbar(menu, item, actions)
                        elif isinstance(menu, QMenu):
                            populate_sub_qmenu(menu, item, actions)
                        else:
                            raise TypeError(menu)
                        continue
                    elif not isinstance(item, str):
                        raise RuntimeError('what is this...action item() = %r' % item())

                    try:
                        action = self.actions[item] #if isinstance(item, str) else item()
                    except:
                        print(self.actions.keys())
                        raise
                    menu.addAction(action)
        #self._create_plane_from_points(None)

    def _update_menu(self, menu_items):
        assert isinstance(menu_items, dict), menu_items
        for unused_name, (menu, unused_items) in menu_items.items():
            menu.clear()
        self._populate_menu(menu_items)

    #def _create_plane_from_points(self, points):
        #origin, vx, vy, vz, x_limits, y_limits = self._fit_plane(points)

        ## We create a 100 by 100 point plane to sample
        #splane = vtk.vtkPlaneSource()
        #plane = splane.GetOutput()

        #dx = max(x_limits) - min(x_limits)
        #dy = max(y_limits) - min(y_limits)
        #dx = 1.
        #dy = 3.

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
        self._prepare_actions_helper(icon_path, tools, self.actions,
                                     checkables=checkables)

        self.actions['toolbar'] = self.toolbar.toggleViewAction()
        self.actions['toolbar'].setStatusTip("Show/Hide application toolbar")

        self.actions['reswidget'] = self.res_dock.toggleViewAction()
        self.actions['reswidget'].setStatusTip("Show/Hide results selection")
        return self.actions

    def _prepare_actions_helper(self, icon_path, tools, actions, checkables=None):
        """
        Prepare actions that will  be used in application in a way
        that's independent of the  menus & toolbar
        """
        if checkables is None:
            checkables = []

        for tool in tools:
            (name, txt, icon, shortcut, tip, func) = tool
            if name in actions:
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
                actions[name] = QAction(ico, txt, self, checkable=True)
                actions[name].setChecked(is_checked)
            else:
                actions[name] = QAction(ico, txt, self)

            if shortcut:
                actions[name].setShortcut(shortcut)
                #actions[name].setShortcutContext(QtCore.Qt.WidgetShortcut)
            if tip:
                actions[name].setStatusTip(tip)
            if func:
                actions[name].triggered.connect(func)

    def _logg_msg(self, log_type, filename, lineno, msg):
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

    def _log_msg(self, msg):
        """prints an HTML log message"""
        self.log_mutex.lockForWrite()
        text_cursor = self.log_widget.textCursor()
        end = text_cursor.End
        text_cursor.movePosition(end)
        text_cursor.insertHtml(msg)
        self.log_widget.ensureCursorVisible() # new message will be visible
        self.log_mutex.unlock()

    def log_info(self, msg):
        """ Helper funtion: log a message msg with a 'INFO:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        return self.log.simple_msg(msg, 'GUI INFO')

    def log_debug(self, msg):
        """ Helper funtion: log a message msg with a 'DEBUG:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        return self.log.simple_msg(msg, 'GUI DEBUG')

    def log_command(self, msg):
        """ Helper funtion: log a message msg with a 'COMMAND:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        return self.log.simple_msg(msg, 'GUI COMMAND')

    def log_error(self, msg):
        """ Helper funtion: log a message msg with a 'GUI ERROR:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        return self.log.simple_msg(msg, 'GUI ERROR')

    def log_warning(self, msg):
        """ Helper funtion: log a message msg with a 'WARNING:' prefix """
        if msg is None:
            msg = 'msg is None; must be a string'
            return self.log.simple_msg(msg, 'GUI ERROR')
        return self.log.simple_msg(msg, 'GUI WARNING')

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
            raise NotImplementedError('invalid image type=%r; filename=%r' % (
                fmt, image_filename))

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

        xcentroid = origin[0] + 0.5 * (extent[0] + extent[1]) * spacing[0]
        ycentroid = origin[1] + 0.5 * (extent[2] + extent[3]) * spacing[1]
        #xd = (extent[1] - extent[0] + 1) * spacing[0]
        yd = (extent[3] - extent[2] + 1) * spacing[1]
        distance = camera.GetDistance()
        camera.SetParallelScale(0.5 * yd)
        camera.SetFocalPoint(xcentroid, ycentroid, 0.0)
        camera.SetPosition(xcentroid, ycentroid, distance)

    def _create_vtk_objects(self):
        """creates some of the vtk objects"""
        #Frame that VTK will render on
        self.vtk_frame = QFrame()

        # can't build an interactor without a GUI (for testing)
        self.vtk_interactor = QVTKRenderWindowInteractor(parent=self.vtk_frame)
        #self.vtk_interactor = PyNastranRenderWindowInteractor(parent=self.vtk_frame)
        #self.set_anti_aliasing(2)

        #self._camera_event_name = 'LeftButtonPressEvent'
        self.mouse_actions.setup_mouse_buttons(mode='default')

    def build_vtk_frame(self):
        """uses the vtk objects to set up the window (frame)"""
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

        #for cid, axes in self.axes.items():
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
        self._build_vtk_frame_post()

    def on_reset_camera(self):
        self.log_command('on_reset_camera()')
        self._simulate_key_press('r')
        self.vtk_interactor.Render()

    def on_flip_edges(self):
        """turn edges on/off"""
        self.is_edges = not self.is_edges
        self.edge_actor.SetVisibility(self.is_edges)
        # cart3d edge color isn't black...
        #self.edge_actor.GetProperty().SetColor(0, 0, 0)
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

    #---------------------------------------------------------------------
    # groups

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
            # doesn't work for the bwb_saero.bdf
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
        #print('flip_flag=%s show_flag=%s' % (flip_flag, show_flag))

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

    def start_logging(self):
        if self.log is not None:
            return
        if self.html_logging is True:
            log = SimpleLogger(
                'debug', 'utf-8',
                lambda w, x, y, z: self._logg_msg(w, x, y, z))
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
                                raise_error=False):
        """action version of ``on_load_geometry``"""
        self.on_load_geometry(infile_name=infile_name, geometry_format=geometry_format,
                              name=name, plot=True, raise_error=raise_error)

    def _update_menu_bar_to_format(self, fmt, method):
        """customizes the gui to be nastran/cart3d-focused"""
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

        # TODO: what is cwo?
        self.menu_bar_format = 'cwo'
        if self.menu_bar_format is None:
            self._update_menu_bar_to_format(self.format, method_new)
        else:
            if not pyNastran.is_pynastrangui_exe:  # pragma: no cover
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

    def setup_gui(self, is_gui=True):
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
        if is_gui:
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

    def make_gif(self, gif_filename, scale, istep=None,
                 min_value=None, max_value=None,
                 animate_scale=True, animate_phase=False, animate_time=False,
                 icase_fringe=None, icase_disp=None, icase_vector=None,
                 animate_fringe=False, animate_vector=False,
                 icase_start=None, icase_end=None, icase_delta=None,
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
        icase_fringe/disp/vector : int; default=None
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
            self.stop_animation()
            return False

        is_failed = True
        try:
            out = setup_animation(
                scale, istep=istep,
                animate_scale=animate_scale, animate_phase=animate_phase, animate_time=animate_time,
                icase_fringe=icase_fringe, icase_disp=icase_disp, icase_vector=icase_vector,
                icase_start=icase_start, icase_end=icase_end, icase_delta=icase_delta,
                time=time, animation_profile=animation_profile,
                fps=fps, animate_in_gui=animate_in_gui)
        except (AssertionError, ValueError, RuntimeError, NotImplementedError) as error:
            self.log_error(str(error))
            self.stop_animation()
            return is_failed
        (phases, icases_fringe, icases_disp, icases_vector,
         isteps, scales,
         analysis_time, onesided, unused_endpoint) = out

        if animate_time:
            icase_msg = '         icase_start=%s, icase_end=%s, icase_delta=%s,\n' % (
                icase_start, icase_end, icase_delta,)
        else:
            icase_msg = (
                '         icase_fringe=%s, icase_disp=%s, icase_vector=%s, \n'
                '         animate_fringe=%s, animate_vector=%s, \n' % (
                    icase_fringe, icase_disp, icase_vector,
                    animate_fringe, animate_vector,
                ))

        #animate_in_gui = True
        self.stop_animation()
        if len(icases_disp) == 1:
            pass
        elif animate_in_gui:
            msg = (
                'make_gif(%r, %s, istep=%s,\n'
                '         min_value=%s, max_value=%s,\n'
                '         animate_scale=%s, animate_phase=%s, animate_time=%s,\n%s'
                #'         icase_fringe=%s, icase_disp=%s, icase_vector=%s, \n'
                #'         icase_start=%s, icase_end=%s, icase_delta=%s,\n'
                "         time=%s, animation_profile=%r,\n"
                '         fps=%s, stop_animation=%s, animate_in_gui=%s)\n' % (
                    gif_filename, scale, istep, min_value, max_value,
                    animate_scale, animate_phase, animate_time,
                    icase_msg,
                    #icase_fringe, icase_disp, icase_vector,
                    #icase_start, icase_end, icase_delta,
                    time, animation_profile,
                    fps, stop_animation, animate_in_gui)
            )
            self.log_command(msg)
            # onesided has no advantages for in-gui animations and creates confusion
            scales, phases, icases_fringe, icases_disp, icases_vector, isteps = make_two_sided(
                scales, phases, icases_fringe, icases_disp, icases_vector, isteps, onesided)

            self._animate_in_gui(
                min_value, max_value,
                scales, phases,
                icases_fringe, icases_disp, icases_vector,
                animate_fringe, animate_vector,
                fps)
            is_failed = False
            return is_failed

        try:
            is_failed = self.make_gif_helper(
                gif_filename, icases_fringe, icases_disp, icases_vector, scales,
                phases=phases, isteps=isteps,
                animate_fringe=animate_fringe, animate_vector=animate_vector,
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
                '         animate_scale=%s, animate_phase=%s, animate_time=%s,\n%s'
                "         time=%s, animation_profile=%r,\n"
                '         nrepeat=%s, fps=%s, magnify=%s,\n'
                '         make_images=%s, delete_images=%s, make_gif=%s, stop_animation=%s,\n'
                '         animate_in_gui=%s)\n' % (
                    gif_filename, scale, istep, min_value, max_value,
                    animate_scale, animate_phase, animate_time,
                    icase_msg,
                    time, animation_profile,
                    nrepeat, fps, magnify, make_images, delete_images, make_gif, stop_animation,
                    animate_in_gui)
            )
            self.log_command(msg)

        return is_failed

    def _animate_in_gui(self, min_value, max_value,
                        scales, phases,
                        icases_fringe, icases_disp, icases_vector,
                        animate_fringe, animate_vector,
                        fps):
        """helper method for ``make_gif``"""
        callback = AnimationCallback(self, scales, phases,
                                     icases_fringe, icases_disp, icases_vector,
                                     animate_fringe, animate_vector,
                                     min_value, max_value)

        # Sign up to receive TimerEvent
        observer_name = self.vtk_interactor.AddObserver('TimerEvent', callback.execute)
        self.observers['TimerEvent'] = observer_name

        # total_time not needed
        # fps
        # -> frames_per_second = 1/fps
        delay = int(1. / fps * 1000)

        # time in milliseconds
        unused_timer_id = self.vtk_interactor.CreateRepeatingTimer(delay)

    def stop_animation(self):
        """removes the animation timer"""
        is_failed = False
        if 'TimerEvent' in self.observers:
            observer_name = self.observers['TimerEvent']
            self.vtk_interactor.RemoveObserver(observer_name)
            del self.observers['TimerEvent']
            self.mouse_actions.setup_mouse_buttons(mode='default', force=True)
        return is_failed

    def animation_update(self, icase_fringe0, icase_disp0, icase_vector0,
                         icase_fringe, icase_disp, icase_vector, scale, phase,
                         animate_fringe, unused_animate_vector,
                         normalized_frings_scale,
                         min_value, max_value):
        """applies the animation update callback"""
        #print('icase_fringe=%r icase_fringe0=%r' % (icase_fringe, icase_fringe0))
        arrow_scale = None  # self.glyph_scale_factor * scale
        icase_vector = None
        is_legend_shown = self.scalar_bar.is_shown
        if icase_disp != icase_disp0:
            # apply the fringe
            #
            # min/max value is used only for the time plot
            # it's assumed to be a displacement result, so the fringe=displacement
            self.cycle_results_explicit(icase_disp, explicit=True,
                                        min_value=min_value, max_value=max_value)

        if icase_fringe is not None and icase_fringe != icase_fringe0:
            is_valid = self.on_fringe(icase_fringe,
                                      update_legend_window=False, show_msg=False)
            if is_legend_shown:
                # TODO: sort of a hack for the animation
                # the fringe always shows the legend, but we may not want that
                # just use whatever is active
                self.show_legend()

            if not is_valid:
                self.log_error('Invalid Fringe Case %i' % icase_fringe)
                return False

        is_valid = self.animation_update_fringe(
            icase_fringe, animate_fringe, normalized_frings_scale)
        if not is_valid:
            return is_valid

        try:
            # apply the deflection
            self.update_grid_by_icase_scale_phase(icase_disp, scale, phase=phase)
        except(AttributeError, KeyError):
            self.log_error('Invalid Displacement Case %i' % icase_disp)
            return False

        if icase_vector is not None and icase_vector != icase_vector0:
            try:
                # apply the nodal forces
                self.update_forces_by_icase_scale_phase(icase_vector, arrow_scale, phase=phase)
            except(AttributeError, KeyError):
                self.log_error('Invalid Vector Case %i' % icase_vector)
                return False
        is_valid = True
        return is_valid

    def animation_update_fringe(self, icase_fringe, animate_fringe, normalized_frings_scale):
        """helper method for ``animation_update``"""
        if animate_fringe:
            # e^(i*(theta + phase)) = sin(theta + phase) + i*cos(theta + phase)
            is_valid, data = self._update_vtk_fringe(icase_fringe, normalized_frings_scale)
            if not is_valid:
                return is_valid

            #icase = data.icase
            result_type = data.result_type
            #location = data.location
            min_value = data.min_value
            max_value = data.max_value
            #norm_value = data.norm_value
            #imin = data.imin
            #imax = data.imax
            data_format = data.data_format
            nlabels = data.nlabels
            labelsize = data.labelsize
            ncolors = data.ncolors
            colormap = data.colormap
            #subcase_id = data.subcase_id
            #subtitle = data.subtitle
            #label = data.label

            is_legend_shown = self.scalar_bar.is_shown
            self.update_scalar_bar(result_type, min_value, max_value,
                                   data_format,
                                   nlabels=nlabels, labelsize=labelsize,
                                   ncolors=ncolors, colormap=colormap,
                                   is_shown=is_legend_shown)
            #obj.get_vector_array_by_phase(i, name, )
        is_valid = True
        return is_valid

    def make_gif_helper(self, gif_filename, icases_fringe, icases_disp, icases_vector,
                        scales, phases=None, isteps=None,
                        animate_fringe=False, animate_vector=False,
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
        icases_fringe/disp/vector : int / List[int]
            the result case to plot the deflection for
        scales : List[float]
            List[float] : the deflection scale factors; true scale
        phases : List[float]; default=None
            List[float] : the phase angles (degrees)
            None -> animate scale
        max_value : float; default=None
            the max value on the plot
        min_value : float; default=None
            the min value on the plot
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

        phases, icases_fringe, icases_disp, icases_vector, isteps, scales = update_animation_inputs(
            phases, icases_fringe, icases_disp, icases_vector,
            isteps, scales, analysis_time, fps)

        if gif_filename is not None:
            png_dirname = os.path.dirname(os.path.abspath(gif_filename))
            if not os.path.exists(png_dirname):
                os.makedirs(png_dirname)

            png_filenames = []
            fmt = gif_filename[:-4] + '_%%0%ii.png' % (len(str(nframes)))

        icase_fringe0 = -1
        icase_disp0 = -1
        icase_vector0 = -1
        is_failed = True
        if make_images:
            scale_max = max(abs(scales.max()), abs(scales.min()))
            for istep, icase_fringe, icase_disp, icase_vector, scale, phase in zip(
                    isteps, icases_fringe, icases_disp, icases_vector, scales, phases):

                normalized_frings_scale = scale / scale_max
                is_valid = self.animation_update(
                    icase_fringe0, icase_disp0, icase_vector0,
                    icase_fringe, icase_disp, icase_vector,
                    scale, phase,
                    animate_fringe, animate_vector,
                    normalized_frings_scale,
                    min_value, max_value)
                if not is_valid:
                    return is_failed
                if gif_filename is not None:
                    png_filename = fmt % istep
                    self.on_take_screenshot(fname=png_filename, magnify=magnify)
                    png_filenames.append(png_filename)
        else:
            for istep in isteps:
                png_filename = fmt % istep
                png_filenames.append(png_filename)
                check_path(png_filename, 'png_filename')

        if gif_filename is not None and png_filenames:
            is_failed = write_gif(
                gif_filename, png_filenames, time=time,
                onesided=onesided,
                nrepeat=nrepeat, delete_images=delete_images,
                make_gif=make_gif)
        return is_failed

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
        #print("key_press = ", key)
        if key == 'f':  # change focal point
            #print('focal_point!')
            return
        self.vtk_interactor._Iren.SetEventInformation(0, 0, 0, 0, key, 0, None)
        self.vtk_interactor._Iren.KeyPressEvent()
        self.vtk_interactor._Iren.CharEvent()

        #if key in ['y', 'z', 'X', 'Y', 'Z']:
            #self.update_camera(key)

    def _finish_results_io2(self, model_name, form, cases, reset_labels=True):
        """
        Adds results to the Sidebar

        Parameters
        ----------
        model_name : str
            the name of the model; unused
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
            self.scalar_bar_actor.VisibilityOn()
            self.scalar_bar_actor.Modified()

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

        self.res_widget.set_case_keys(self.case_keys)
        self.res_widget.update_results(form, self.name)

        key = self.case_keys[0]
        location = self.get_case_location(key)
        method = 'centroid' if location else 'nodal'

        data2 = [(method, None, [])]
        self.res_widget.update_methods(data2)

        if self.node_ids is None:  # pragma: no cover
            raise RuntimeError('implement self.node_ids for this format')
        #if self.element_ids is None:  # pragma: no cover
            #raise RuntimeError('implement self.element_ids for this format')

        if self.is_groups:
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

        for unused_module_name, module in self.modules.items():
            module.post_load_geometry()

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
        for unused_group_id, group in sorted(out_data.items()):
            if not isinstance(group, Group):
                continue
            data[group.name] = group
        self.groups = data

def populate_sub_qmenu(menu, items, actions):
    sub_menu_name = items[0]
    sub_menu = menu.addMenu(sub_menu_name)
    for ii_count, ii in enumerate(items[1:]):
        if not isinstance(ii, str):
            raise RuntimeError('what is this...action ii() = %r' % ii())
        action = actions[ii]
        if ii_count > 0:
            action.setChecked(False)
        sub_menu.addAction(action)

def populate_sub_qtoolbar(toolbar, items, actions):
    """
    refs
    https://www.walletfox.com/course/customqtoolbutton.php
    https://stackoverflow.com/questions/9076332/qt-pyqt-how-do-i-create-a-drop-down-widget-such-as-a-qlabel-qtextbrowser-etc
    """
    sub_menu_name = items[0]
    action0 = actions[sub_menu_name]

    drop_down_menu = QMenu()

    custom_button = QToolButton()
    custom_button.setPopupMode(QToolButton.InstantPopup)
    custom_button.setMenu(drop_down_menu)
    custom_button.setDefaultAction(action0)

    toolbar.addWidget(custom_button)

    for unused_ii_count, itemi in enumerate(items[1:]):
        if not isinstance(itemi, str):
            raise RuntimeError('what is this...action ii() = %r' % itemi())
        action = actions[itemi]
        # temp
        #if ii_count > 0:
            #action.setChecked(False)
        drop_down_menu.addAction(action)  # thrown in the trash?
