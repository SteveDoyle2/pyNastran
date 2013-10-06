# -*- coding: utf-8 -*-
#try:
#    from PySide import QtCore, QtGui
#    fmode = 1
#except ImportError:
#    try:
#        from PyQt4 import QtCore, QtGui
#        fmode = 2
#    except ImportError:
#        fmode = None

from PyQt4 import QtGui, QtCore
import sys
import os.path
import cgi #  html lib
import datetime

from numpy import ndarray, amax, amin

from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import vtk

import pyNastran
from pyNastran import __version__ as version
from pyNastran.utils import print_bad_path
from pyNastran.utils.log import SimpleLogger
from nastranIO import NastranIO
pkg_path = pyNastran.__path__[0]
icon_path = os.path.join(pkg_path, 'gui', 'icons')
#image_path = os.path.join(pkg_path, 'gui_qt', 'images')

#### tcolorpick.png and tabout.png trefresh.png icons on LGPL license, see
#### http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/color-picker-grey.png
#### http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/help-hint.png
#### http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/view-refresh-8.png


from pyNastran.gui.formats import NastranIO, Cart3dIO, PanairIO, LaWGS_IO, is_nastran, is_cart3d, is_panair, is_lawgs
from pyNastran.gui.gui_inputs import get_inputs

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

class MainWindow(QtGui.QMainWindow, NastranIO, Cart3dIO, PanairIO, LaWGS_IO):
 
    def __init__(self, is_edges=False, is_nodal=False, is_centroidal=False,
                 format=None, input=None, output=None, shots=None, magnify=4, rotation=None, debug=False):
        assert debug in [True, False], 'debug=%s' % debug

        is_centroidal = not(is_nodal)
        QtGui.QMainWindow.__init__(self)
        NastranIO.__init__(self)
        Cart3dIO.__init__(self)
        PanairIO.__init__(self)
        LaWGS_IO.__init__(self)
        
        settings = QtCore.QSettings()
        self.is_nodal = is_nodal
        self.is_centroidal = is_centroidal
        self.nCases = 0
        self.iCase = 0
        self.nNodes = 0
        self.nElements = 0
        
        self.grid = vtk.vtkUnstructuredGrid()
        self.gridResult = vtk.vtkFloatArray()

        self.grid2 = vtk.vtkUnstructuredGrid()
        #self.emptyResult = vtk.vtkFloatArray()
        #self.vectorResult = vtk.vtkFloatArray()

        # edges
        self.edgeActor = vtk.vtkActor()
        self.edgeMapper = vtk.vtkPolyDataMapper()

        self.scalarBar = vtk.vtkScalarBarActor()
        
        self.format = '' # the active format loaded into the gui
        self.last_dir = '' # last visited directory while opening file

        # build GUI and restore saved application state
        self.restoreGeometry(settings.value("mainWindowGeometry").toByteArray())
        self.background_col = settings.value("backgroundColor", (0.1, 0.2, 0.4)).toPyObject()

        self.init_ui()
        self.restoreState(settings.value("mainWindowState").toByteArray())

        self.log =  SimpleLogger('debug', lambda x, y: self.logg_msg(x, y))
        # logging needs synchronizing, so the messages from different threads
        # would not be interleave
        self.log_mutex = QtCore.QReadWriteLock() 
        self.show()
        self.on_load_geometry('solid_bending.bdf', 'Nastran')
        self.on_load_results('solid_bending.op2')
        self.vtk_interactor.Modified()

    def logg_msg(self, typ, msg):
        """
        Add message to log widget trying to choose right color for it.
        @param msg message to be displayed
        """
        tim = datetime.datetime.now().strftime('[%d-%m-%Y %H:%M:%S]')
        msg = cgi.escape(msg)
        #message colors
        cols = {"GUI": "blue", "GUI ERROR":"Crimson", "DEBUG" : "Orange"}
        msg = msg.rstrip().replace('\n', '<br>')
        msg = tim + ' ' + (typ + ': ' + msg) if typ else msg
        if typ in cols:
            msg = '<font color="%s"> %s </font>' % (cols[typ], msg)
        
        self.log_mutex.lockForWrite()
        self.log_widget.textCursor().insertHtml(msg + r"<br />")
        self.log_widget.ensureCursorVisible() # new message will be visible
        self.log_mutex.unlock()

    def log_info(self, msg):
        """ Helper funtion: log a messaage msg with a 'GUI:' prefix """
        self.log.simple_msg(msg, 'GUI')

    def log_debug(self, msg):
        """ Helper funtion: log a messaage msg with a 'GUI:' prefix """
        self.log.simple_msg(msg, 'GUI ERROR')

    def change_background_col(self):
        """ Choose a background color """
        c =  [int(255 * i) for i in self.background_col]
        col = QtGui.QColorDialog.getColor(QtGui.QColor(*c), self, "Choose a background color")
        if col.isValid():
            self.background_col = col.getRgbF()[:3]
            self.rend.SetBackground(*self.background_col)

    def about_dialog(self):
        """ Display about dialog """
        about = [
            'pyNastran QT GUI',
            '',
            'pyNastran v%s' % (pyNastran.__version__),
            'Copyright '+ pyNastran.__license__,
            '  Steven Doyle 2011-2012',
            '  Marcin Gasiorek 2012',
            '',
            '%s' % (pyNastran.__website__),
            '',
            'Mouse',
            'Left Click - Rotate',
            'Middle Click - Pan/Recenter Rotation Point',
            'Shift + Left Click - Pan/Recenter Rotation Point',
            'Right Mouse - Zoom',
            '',
            'Keyboard Controls',
            'r   - reset camera view',
            'X/x - snap to x axis',
            'Y/y - snap to y axis',
            'Z/z - snap to z axis',
            '',
            #'h   - show/hide legend & info',
            'CTRL+I - take a screenshot (image)',
            #'L     - cycle op2 results',
            'm/M    - scale up/scale down by 1.1 times',
            'o/O    - rotate counter-clockwise/clockwise 5 degrees',
            's      - view model as a surface',
            'w      - view model as a wireframe',
        ]
        QtGui.QMessageBox.about(self, "About pyNastran GUI", "\n".join(about))

    def set_window_title(self, fname=None):
        msg = "pyNastran v%s"  % version
        if fname:
            msg += ' - %s' % fname
        self.setWindowTitle(msg)

    def init_ui(self):
        """ Initialize user iterface"""
        self.resize(800,600)
        self.statusBar().showMessage('Ready')

        # windows title and aplication icon
        self.setWindowTitle('Statusbar')    
        self.setWindowIcon(QtGui.QIcon(os.path.join(icon_path, 'logo.png')))
        self.set_window_title()

        ############  Results widget ##################
        self.res_dock = QtGui.QDockWidget("Results", self)
        self.res_dock.setObjectName("results_obj")
        self.res_widget = QtGui.QTextEdit()
        self.res_widget.setReadOnly(True)
        self.res_dock.setWidget(self.res_widget)
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.res_dock)
        ############  Logging widget ##################
        self.log_dock = QtGui.QDockWidget("Application log", self)
        self.log_dock.setObjectName("application_log")
        self.log_widget = QtGui.QTextEdit()
        self.log_widget.setReadOnly(True)
        self.log_dock.setWidget(self.log_widget)
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock)
        ################################################
        
        # right sidebar

        ## menubar
        self.menubar = self.menuBar()
        self.menu_file = self.menubar.addMenu('&File')
        self.menu_view = self.menubar.addMenu('&View')
        self.menu_window = self.menubar.addMenu('&Window')
        self.menu_help = self.menubar.addMenu('&Help')

        ## toolbar
        self.toolbar = self.addToolBar('Show toolbar')
        self.toolbar.setObjectName('main_toolbar')

        # prepare actions that will  be used in application
        actions = {}
        #pth = os.path.join(icon_path, 'tbdf.png')
        pth = os.path.join(icon_path, 'tbdf.png')
        #print print_bad_path(pth)
        
        # http://docs.python.org/2/library/sys.html#sys.platform
        #System  platform value
        #Linux (2.x and 3.x) 'linux2'
        #Windows 'win32'
        #Windows/Cygwin  'cygwin'
        #Mac OS X    'darwin'
        #print sys.platform
        #quit_key = 'Alt+F4' if sys.platform in ['win32', 'cygwin'] else 'Ctrl+Q'
        
        for nam, txt, icon, shortcut, tip, func in [
          ('exit', '&Exit', os.path.join(icon_path, 'texit.png'), 'Ctrl+Q', 'Exit application', QtGui.qApp.quit),
          ('open_bdf', '&Open BDF', os.path.join(icon_path, 'tbdf.png'), 'Ctrl+O', 'Loads a Geometry input file', self.on_load_geometry),
          ('open_op2', '&Open OP2', os.path.join(icon_path, 'top2.png'), None, 'Loads a OP2 results file', self.load_op2),
          ('open_f06', '&Open F06', None, None, 'Loads a F06 results file', self.load_f06),  ## @todo no picture...
          ('back_col', 'Change background color', os.path.join(icon_path, 'tcolorpick.png'), None, 'Choose a background color', self.change_background_col),

          ('wireframe', 'Wireframe Model', os.path.join(icon_path, 'twireframe.png'), 'w', 'Show Model as a Wireframe Model', lambda: self._simulate_key_press('w')),
          ('surface', 'Surface Model', os.path.join(icon_path, 'tsolid.png'), 's', 'Show Model as a Surface Model', lambda: self._simulate_key_press('s')),
          ('edges', 'Show/Hide Edges', os.path.join(icon_path, 'tedges.png'), 'e', 'Show/Hide Model Edges', lambda: self.onFlipEdges),
          ('scshot', 'Take a Screenshot', os.path.join(icon_path, 'tcamera.png'), 'CTRL+I', 'Take a Screenshot of current view', self.take_screenshot),
          ('about', 'About pyNastran GUI', os.path.join(icon_path, 'tabout.png'), 'CTRL+H', 'About pyNastran GUI and help on shortcuts', self.about_dialog),
          ('creset', 'Reset camera view', os.path.join(icon_path, 'trefresh.png'), 'r', 'Reset the camera view to default', lambda: self._simulate_key_press('r')),

          ('cycle_res', 'Cycle Results', os.path.join(icon_path, 'cycle_results.png'), 'CTRL+L', 'Changes the result case', self.cycleResults()),

          ('x', 'Flips to +X Axis', os.path.join(icon_path, '+x.png'), 'x', 'Flips to +X Axis', lambda: self.update_camera('x')),
          ('y', 'Flips to +Y Axis', os.path.join(icon_path, '+x.png'), 'y', 'Flips to +Y Axis', lambda: self.update_camera('y')),
          ('z', 'Flips to +Z Axis', os.path.join(icon_path, '+x.png'), 'z', 'Flips to +Z Axis', lambda: self.update_camera('z')),
          
          ('X', 'Flips to -X Axis', os.path.join(icon_path, '-z.png'), 'X', 'Flips to -X Axis', lambda: self.update_camera('X')),
          ('Y', 'Flips to -Y Axis', os.path.join(icon_path, '-z.png'), 'Y', 'Flips to -Y Axis', lambda: self.update_camera('Y')),
          ('Z', 'Flips to -Z Axis', os.path.join(icon_path, '-z.png'), 'Z', 'Flips to -Z Axis', lambda: self.update_camera('Z')), ]:
            #print "name=%s txt=%s icon=%s short=%s tip=%s func=%s" % (nam, txt, icon, short, tip, func)
            if icon is None:
                print "missing_icon = %r!!!" % nam
                #print print_bad_path(icon)
            #elif not "/" in icon:
                #ico = QtGui.QIcon.fromTheme(icon)
            else:
                ico = QtGui.QIcon()
                ico.addPixmap(QtGui.QPixmap(icon), QtGui.QIcon.Normal, QtGui.QIcon.Off)
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

        actions['logwidget'] = self.log_dock.toggleViewAction()
        actions['logwidget'].setStatusTip("Show/Hide application log")


        # populate menus and toolbar
        for menu, items in [(self.menu_file, ('open_bdf', 'open_op2', 'open_f06', '', 'exit')),
                           (self.menu_view,  ('scshot', '', 'wireframe', 'surface', 'creset', '', 'back_col')),
                           (self.menu_window,('toolbar', 'reswidget', 'logwidget')),
                           (self.menu_help,  ('about',)),
                           (self.toolbar, ('open_bdf', 'open_op2', 'cycle_res', 'x', 'wireframe', 'surface', 'edges', 'creset', 'scshot', '', 'exit'))]:
            for i in items:
                if not i:
                    menu.addSeparator()
                else:
                    menu.addAction(actions[i] if isinstance(i, basestring) else i())
        self.res_dock.hide()

        #Frame that VTK will render on 
        vtk_frame = QtGui.QFrame()
        vtk_hbox  = QtGui.QHBoxLayout()
        vtk_hbox.setContentsMargins(2, 2, 2, 2)

        #Qt VTK RenderWindowInteractor
        self.vtk_interactor = QVTKRenderWindowInteractor(parent = vtk_frame)
        vtk_hbox.addWidget(self.vtk_interactor)
        vtk_frame.setLayout(vtk_hbox)
        vtk_frame.setFrameStyle(QtGui.QFrame.NoFrame | QtGui.QFrame.Plain)
        # this is our main, 'central' widget
        self.setCentralWidget(vtk_frame)

        ##############################################################
        self.rend = vtk.vtkRenderer()
        self.vtk_interactor.GetRenderWindow().AddRenderer(self.rend)
        self.vtk_interactor.GetRenderWindow().Render()
        self.load_nastran_geometry(None, None)
        self.textActors = {}

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

    def onFlipEdges(self):
        self.is_edges = not(self.is_edges)
        self.edgeActor.SetVisibility(self.is_edges)
        #self.edgeActor.GetProperty().SetColor(0, 0, 0)  # cart3d edge color isn't black...
        self.edgeActor.Modified()
        #self.widget.Update()
        self.Refresh()

    def createText(self, position, label, textSize=18, movable=False):
        txt = vtk.vtkTextActor()
        txt.SetInput(label)
        txtprop = txt.GetTextProperty()
        #txtprop.SetFontFamilyToArial()
        txtprop.SetFontSize(textSize)
        txtprop.SetColor(1, 1, 1)
        txt.SetDisplayPosition(*position)

        #print "dir(text) = ",dir(txt)
        txt.VisibilityOff()

        #txt.SetDisplayPosition(5,5) # bottom left
        #txt.SetDisplayPosition(5,95)
        #txt.SetPosition(0.1,0.5)

        # assign actor to the renderer
        self.rend.AddActor(txt)
        self.textActors[self.iText] = txt
        self.iText += 1

    def TurnTextOff(self):
        for (i, text) in self.textActors.iteritems():
            text.VisibilityOff()

    def TurnTextOn(self):
        for (i, text) in self.textActors.iteritems():
            text.VisibilityOn()

    def build_lookup_table(self):
        self.colorFunction = vtk.vtkColorTransferFunction()
        self.colorFunction.SetColorSpaceToHSV()
        self.colorFunction.HSVWrapOff()

        drange = [10., 20.]
        self.colorFunction.AddRGBPoint(drange[0], 0.0, 0.0, 1.0)
        self.colorFunction.AddRGBPoint(drange[1], 1.0, 0.0, 0.0)

        self.scalarBar.SetTitle("Title1")
        self.scalarBar.SetLookupTable(self.colorFunction)
        self.scalarBar.SetOrientationToVertical()

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

    def on_load_geometry(self, infile_name=None, geometry_format=None):
        wildcard = ''
        
        if infile_name:
            if geometry_format == 'Nastran':
                has_results = True
                load_function = self.load_nastran_geometry
            elif geometry_format == 'Cart3d':
                has_results = True
                load_function = self.load_cart3d_geometry
            elif geometry_format == 'Panair':
                has_results = False
                load_function = None
            elif geometry_format == 'LaWGS':
                has_results = False
                load_function = None
            else:
                raise NotImplementedError('on_load_geometry; infile_name=%r format=%r' % (infile_name, geometry_format))
            formats = [geometry_format]
            filter_index = 0
        else:
            formats = []
            load_functions = []
            has_results_list = []
            wildcard_list = []
            if is_nastran:
                wildcard_list.append("Nastran BDF (*.bdf; *.dat; *.nas)")
                formats.append('Nastran')
                has_results_list.append(True)
                load_functions.append(self.load_nastran_geometry)
                #load_functions.append(None)
            if is_cart3d:
                wildcard_list.append("Cart3d (*.tri; *.triq)")
                formats.append('Cart3d')
                has_results_list.append(True)
                load_functions.append(self.load_cart3d_geometry)
            if is_panair:
                wildcard_list.append("Panair (*.inp)")
                formats.append('Panair')
                has_results_list.append(False)
                load_functions.append(self.load_panair_geometry)
            if is_lawgs:
                wildcard_list.append("LaWGS (*.inp)")
                formats.append('LaWGS')
                has_results_list.append(False)
                load_functions.append(None)
            wildcard = ';;'.join(wildcard_list)

            # get the filter index and filename
            if infile_name is not None and geometry_format is not None:
                filter_index = formats.index(geometry_format)
            else:
                Title = 'Choose a Geometry File to Load'
                wildcard_index, infile_name = self._create_load_file_dialog(wildcard, Title)
                #print "infile_name = %r" % infile_name
                #print "wildcard_index = %r" % wildcard_index
                if not infile_name:
                    return # user clicked cancel
                filter_index = wildcard_list.index(wildcard_index)

            geometry_format = formats[filter_index]
            load_function = load_functions[filter_index]
            has_results = has_results_list[filter_index]

        if load_function is not None:
            self.last_dir = os.path.split(infile_name)[0]
            
            #self.grid.Reset()
            self.grid.Modified()
            #self.grid2.Reset()
            self.grid2.Modified()
            #self.gridResult.Reset()
            #self.gridResult.Modified()
            
            self.log_info("reading %s file %r" % (geometry_format, infile_name))
            has_results = load_function(infile_name, self.last_dir)
            #self.vtk_panel.Update()
            self.rend.ResetCamera()

        if filter_index >= 0:
            self.format = formats[filter_index]
            if has_results:
                enable = True
            else:
                enable = False
            #self.load_results.Enable(enable)
            self.set_window_title(infile_name)
        else: # no file specified
            return
        print "on_load_results(%r)" % infile_name

    def _create_load_file_dialog(self, qt_wildcard, Title):
        # getOpenFileName return QString and we want Python string
        fname, wildcard_level = QtGui.QFileDialog.getOpenFileNameAndFilter(self, Title, self.last_dir, qt_wildcard)
        #print "d =", d # fname, wildcard_level
        #fname = str(QtGui.QFileDialog.getOpenFileName(self, Title, self.last_dir, wildcard))
        return str(wildcard_level), str(fname)

    def import_nastran_geometry(self, infile_name, geometry_format):
        """
        Load BDF file.
        """
        # getOpenFileName return QString and we want Python string
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open BDF file', self.last_dir,
        'Nastran BDF (*.bdf *.dat *.nas)'))
        
        if fname:
            self.last_dir = os.path.split(fname)[0]
            self.load_nastran_geometry(fname, self.last_dir)
            self.rend.ResetCamera()
            self.set_window_title(fname)

    def on_load_results(self, infile_name=None):
            geometry_format = self.format
            if self.format is None:
                msg ='on_load_results failed:  You need to load a file first...'
                self.log_debug(msg)
                raise RuntimeError(msg)
            
            if infile_name is None:
                Title = 'Select a Results File for %s' % self.format
                if geometry_format == 'Nastran':
                    has_results = True
                    #wildcard = "Nastran OP2 (*.op2);;Nastran PCH (*.pch);;Nastran F06 (*.f06)"
                    wildcard = "Nastran OP2 (*.op2)"
                    load_functions = [self.load_nastran_results]
                elif geometry_format == 'Cart3d':
                    has_results = True
                    wildcard = "Cart3d (*.triq)"
                    load_functions = [self.load_cart3d_results]
                elif geometry_format == 'Panair':
                    has_results = False
                    wildcard = "Panair (*.agps);;Panair (*.out)"
                    load_functions = [None]
                elif geometry_format == 'LaWGS':
                    has_results = False
                    load_functions = [None]
                else:
                    msg = 'format=%r is not supported' % geometry_format
                    self.log_debug(msg)
                    raise RuntimeError(msg)
                #scard = wildcard.split(';;')
                #n = len(load_functions)
                #wildcard = ';;'.join(scard[:n])
                load_function = load_functions[0]

                wildcard_index, infile_name = self._create_load_file_dialog(wildcard, Title)
            else:
                if geometry_format == 'Nastran':
                    load_function = self.load_nastran_results
                elif geometry_format == 'Cart3d':
                    load_function = self.load_cart3d_results
                #elif geometry_format == 'Panair':
                    #load_function = None
                #elif geometry_format == 'LaWGS':
                    #load_function = None
                else:
                    msg = 'format=%r is not supported.  Did you load a geometry model?' % geometry_format
                    self.log_debug(msg)
                    raise RuntimeError(msg)

            self.last_dir = os.path.split(infile_name)[0]
            load_function(infile_name, self.last_dir)
            print "on_load_results(%r)" % infile_name


    def load_op2(self):
        """
        Load OP2 file.
        @todo not done...
        """
        return self.on_load_results()
        
        # getOpenFileName return QString and we want Python string
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open OP2 file', self.last_dir, 'Nastran OP2 (*.op2)'))
        
        if fname:
            self.last_dir = os.path.split(fname)[0]
            self.load_nastran_results(fname, self.last_dir)
            self.rend.ResetCamera()

    def load_f06(self):
        """
        Load F06 file.
        @todo not done...
        """
        # getOpenFileName return QString and we want Python string
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open F06 file', self.last_dir, 'Nastran F06 (*.f06)'))
        
        if fname:
            self.last_dir = os.path.split(fname)[0]
            self.load_nastran_results(fname, self.last_dir)
            self.rend.ResetCamera()

    def take_screenshot(self):
        """ Take a screenshot of a current view and save as a file"""
        renderLarge = vtk.vtkRenderLargeImage()
        renderLarge.SetInput(self.rend)
        renderLarge.SetMagnification(4)

        filt = QtCore.QString() 
        fname = str(QtGui.QFileDialog.getSaveFileName(self, ('Choose a file name'
                    'and type'), '', ('PNG Image *.png (*.png);; JPEG Image '
                    '*.jpg *.jpeg (*.jpg, *.jpeg);; TIFF Image *.tif *.tiff '
                    '(*.tif, *.tiff);; BMP Image *.bmp (*.bmp);; PostScript '
                    'Document *.ps (*.ps)'), filt))

        if fname:
            flt = str(filt).split()[0]
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

            writer.SetInputConnection(renderLarge.GetOutputPort())
            writer.SetFileName(fname)
            writer.Write()
            self.log_info("Saved screenshot: " + fname)

    def closeEvent(self, event):
        """
        Handling saving state before application when application is being closed.  
        """
        settings = QtCore.QSettings()
        settings.setValue("main_WindowGeometry", self.saveGeometry())
        settings.setValue("mainWindowState", self.saveState())
        settings.setValue("backgroundColor", self.background_col)

    def _simulate_key_press(self, key):
        """
        A little hack method that simulates pressing the key for the VTK
        interactor. There is no easy way to instruct VTK to e.g. change mouse 
        style to 'trackball' (as by pressing 't' key),
        (see http://public.kitware.com/pipermail/vtkusers/2011-November/119996.html)
        therefore we trick VTK to think that a key has been pressed.
        
        :param key: a key that VTK should be informed about, e.g. 't'
        """
        print "key = ", key
        self.vtk_interactor._Iren.SetEventInformation(0, 0, 0, 0, key, 0, None)
        self.vtk_interactor._Iren.KeyPressEvent()
        self.vtk_interactor._Iren.CharEvent()

    def addGeometry(self):
        self.aQuadMapper = vtk.vtkDataSetMapper()
        self.aQuadMapper.SetInput(self.grid)
        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(self.aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 0, 0)  # red
        self.rend.AddActor(geometryActor)

    def addAltGeometry(self):
        self.aQuadMapper = vtk.vtkDataSetMapper()
        self.aQuadMapper.SetInput(self.grid2)
        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(self.aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 1, 0)  # green
        self.rend.AddActor(geometryActor)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()

    def cycleResults(self):
        print("is_nodal=%s is_centroidal=%s" % (self.is_nodal,self.is_centroidal))
        #print("nCases = %i" %(self.nCases+1))
        if self.nCases == 0:
            return

        foundCases = self.incrementCycle()
        if foundCases:
            print("incremented case")
            #self.gridResult.Reset()
            gridResult = vtk.vtkFloatArray()
            #emptyResult = vtk.vtkFloatArray()

            key = self.caseKeys[self.iCase]
            case = self.resultCases[key]
            print("len(case) = %i" % len(case))
            (subcaseID, resultType, vectorSize, location, dataFormat) = key

            if location == 'centroid' and self.is_centroidal:
                #allocationSize = vectorSize*location (where location='centroid'-> self.nElements)
                gridResult.Allocate(self.nElements, 1000)
            elif location == 'nodal' and self.is_nodal:
                #allocationSize = vectorSize*self.nNodes # (where location='node'-> self.nNodes)
                gridResult.Allocate(self.nNodes * vectorSize, 1000)
                gridResult.SetNumberOfComponents(vectorSize)
            else:
                print("***%s skipping" % location)

            #self.iSubcaseNameMap[self.isubcase] = [Subtitle,Label]
            caseName = self.iSubcaseNameMap[subcaseID]
            (subtitle, label) = caseName

            print("subcaseID=%s resultType=%s subtitle=%r label=%r"
                % (subcaseID, resultType, subtitle, label))

            if isinstance(case, ndarray):
                maxValue = amax(case)
                minValue = amin(case)
            else:                
                maxValue = case[0]
                minValue = case[0]
                for value in case:
                    maxValue = max(value, maxValue)
                    minValue = min(value, minValue)

            # flips sign to make colors go from blue -> red
            normValue = maxValue - minValue
            #print "case = ",case
            #if normValue==0.: # avoids division by 0.
            #    normValue = 1.

            valueSet = set()
            if vectorSize == 1:
                #print "minValue = ",min(case)
                for value in case:
                    #print "value", value
                    gridResult.InsertNextValue(value)
                    #if len(valueSet) < 20:
                        #valueSet.add(value)
            else:  # vectorSize=3
                pass
                #for value in case:
                #    self.gridResult.InsertNextTuple3(value)  # x,y,z

            print("max=%g min=%g norm=%g\n" % (maxValue, minValue, normValue))

            nValueSet = len(valueSet)
            if 1:
                self.textActors[0].SetInput('Max:  %g' % maxValue)  # max
                self.textActors[1].SetInput('Min:  %g' % minValue)  # min
                self.textActors[2].SetInput('Subcase=%s Subtitle: %s' % (subcaseID, subtitle))  # info

                if label:
                    self.textActors[3].SetInput('Label: %s' % label)  # info
                    self.textActors[3].VisibilityOn()
                else:
                    self.textActors[3].VisibilityOff()

            self.UpdateScalarBar(resultType, minValue, maxValue, dataFormat)
            #self.scalarBar.SetNumberOfLabels(nValueSet)
            #self.scalarBar.SetMaximumNumberOfColors(nValueSet)
            #prop = self.scalarBar.GetLabelTextProperty()
            #fontSize = prop.GetFontSize()
            #print "fontSize = ",fontSize
            #prop.SetFontSize(40)

            # TODO results can only go from centroid->node and not back to
            ## centroid
            #print dir(self.grid)
            #self.grid.Reset()
            print('gr', gridResult)
            if location == 'centroid' and self.is_centroidal:
                #self.grid.GetPointData().Reset()
                self.grid.GetCellData().SetScalars(gridResult)
                self.log_info("***plotting vector=%s skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
                self.grid.Modified()
            elif location == 'nodal' and self.is_nodal:
                self.grid.GetCellData().Reset()
                if vectorSize == 1:
                    self.log_info("***plotting vector=%s skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
                    self.grid.GetPointData().SetScalars(gridResult)
                    self.grid.Modified()
                else:
                    self.log_info("***nodal vector=%s skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
                    #pass
                    self.grid.GetPointData().SetScalars(gridResult)
                #print "***nodal skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" %(subcaseID,resultType,subtitle,label)
            else:
                self.log_info("***%s skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" % (location, subcaseID, resultType, subtitle, label))


    def incrementCycle(self):
        if self.iCase is not self.nCases:
            self.iCase += 1
        else:
            self.iCase = 0

        if len(self.caseKeys) > 0:
            #print('caseKeys =', self.caseKeys)
            key = self.caseKeys[self.iCase]
            print("key = %s" % (str(key)))
            if key[2] == 3:  # vector size=3 -> vector, skipping ???
                self.incrementCycle()
            foundCases = True
        else:
            self.log_debug("No Results found.  Many results are not supported in the GUI.\n")
            foundCases = False
        #print "next key = ",key
        return foundCases

    def UpdateScalarBar(self, Title, minValue, maxValue, dataFormat):
        """
        @param Title the scalar bar title
        @param minValue the blue value
        @param maxValue the red value
        @param dataFormat '%g','%f','%i', etc.
        """
        self.colorFunction.RemoveAllPoints()
        self.colorFunction.AddRGBPoint(minValue, 0.0, 0.0, 1.0)
        self.colorFunction.AddRGBPoint(maxValue, 1.0, 0.0, 0.0)
        #self.scalarBar.SetLookupTable(self.colorFunction)

        self.scalarBar.SetTitle(Title)
        self.scalarBar.SetLabelFormat(dataFormat)

        nValues = 11
        if (Title in ['Element_ID', 'Eids', 'Region'] and (maxValue - minValue + 1) < 11):
            nValues = int(maxValue - minValue) + 1
            #print "need to adjust axes...maxValue=%s" % maxValue
        #if dataFormat=='%.0f' and maxValue>

        self.scalarBar.SetNumberOfLabels(nValues)
        self.scalarBar.SetMaximumNumberOfColors(nValues)
        self.scalarBar.Modified()

    def ResetCamera(self):
        self.GetCamera().ResetCamera()

    def GetCamera(self):
        #print "getting camera..."
        return self.rend.GetActiveCamera()

    def update_camera(self, code):
        return
        camera = self.GetCamera()
        if code == 'x':  # set x-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(1., 0., 0.)
        elif code == 'X':  # set x-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., -1.)
            camera.SetPosition(-1., 0., 0.)

        elif code == 'y':  # set y-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(0., 1., 0.)
        elif code == 'Y':  # set y-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., -1.)
            camera.SetPosition(0., -1., 0.)

        elif code == 'z':  # set z-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 1., 0.)
            camera.SetPosition(0., 0., 1.)
        elif code == 'Z':  # set z-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., -1., 0.)
            camera.SetPosition(0., 0., -1.)
        #self.rend.Reset()

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    QtGui.QApplication.setOrganizationName("pyNastran")
    QtGui.QApplication.setOrganizationDomain(pyNastran.__website__)
    QtGui.QApplication.setApplicationName("pyNastran")
    QtGui.QApplication.setApplicationVersion(version)
    
    inputs = get_inputs
    window = MainWindow(inputs)
    sys.exit(app.exec_())