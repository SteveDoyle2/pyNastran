# -*- coding: utf-8 -*-
from PyQt4 import QtGui, QtCore
import sys
import os.path
import cgi
import datetime

from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import vtk

import pyNastran
from pyNastran.utils.log import simpleLogger
from nastranIO import NastranIO


#### tcolorpick.png and tabout.png trefresh.png icons on LGPL license, see
#### http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/color-picker-grey.png
#### http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/help-hint.png
#### http://openiconlibrary.sourceforge.net/gallery2/?./Icons/actions/view-refresh-8.png

class MainWindow(QtGui.QMainWindow, NastranIO):
 
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        settings = QtCore.QSettings()
        
        self.last_dir = '' # last visited directory while opening file
        # build GUI and restore saved application state
        self.restoreGeometry(settings.value("mainWindowGeometry").toByteArray())
        self.background_col = settings.value("backgroundColor", (0.1, 0.2, 0.4)).toPyObject()
        
        self.init_ui()
        self.restoreState(settings.value("mainWindowState").toByteArray())
        
        self.log =  simpleLogger('debug', lambda x, y: self.logg_msg(x, y))
        # logging needs synchronizing, so the messages from different threads
        # would not be interleave
        self.log_mutex = QtCore.QReadWriteLock() 
        self.show()
    
    def logg_msg(self, typ, msg):
        """
        Add message to log widget trying to choose right color for it.
        @param msg message to be displayed
        """
        tim = datetime.datetime.now().strftime('[%d-%m-%Y %H:%M:%S]')
        msg = cgi.escape(msg)
        #message colors
        cols = {"GUI": "blue", "DEBUG" : "Crimson"}
        msg = tim + ' ' + (typ + ': ' + msg) if typ else msg
        if typ in cols:
            msg = '<font color="%s"> %s </font>' % (cols[typ], msg)
        
        self.log_mutex.lockForWrite()
        self.log_widget.textCursor().insertHtml(msg + r"<br />")
        self.log_widget.ensureCursorVisible() # new message will be visible
        self.log_mutex.unlock()
    
    def log_info(self, msg):
        """ Helper funtion: log a messaage msg with a 'GUI:' prefix """
        self.log.simple_msg(msg + "\n", 'GUI')
        
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
            'Experimental pyNastran QT GUI',
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
        
    def init_ui(self):
        """ Initialize user iterface"""
        self.resize(800,600)
        self.statusBar().showMessage('Ready')
      
        # windows title and aplication icon
        self.setWindowTitle('Statusbar')    
        self.setWindowIcon(QtGui.QIcon("images/logo.png"))
        self.setWindowTitle("pyNastran experimetnal QT gui")
        
        ############  Logging widget ##################        
        self.log_dock = QtGui.QDockWidget("Application log", self)
        self.log_dock.setObjectName("application_log")
        self.log_widget = QtGui.QTextEdit()
        self.log_widget.setReadOnly(True)
        self.log_dock.setWidget(self.log_widget)
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock)
        ################################################
        
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
        for nam, txt, icon, short, tip, func in [
          ('exit', '&Exit', 'images/texit.png', 'Ctrl+Q', 'Exit application', QtGui.qApp.quit),
          ('open_bdf', '&Open BDF', 'images/topen.png', 'Ctrl+O', 'Loads a BDF input file', self.load_bdf),
          ('wireframe', 'Wireframe kodel', 'images/twireframe.png', 'w', 'Show Model as a Wireframe Model', lambda: self._simulate_key_press('w')),
          ('surface', 'Surface Model', 'images/tsolid.png', 's', 'Show Model as a Surface Model', lambda: self._simulate_key_press('s')),
          ('back_col', 'Change background color', 'images/tcolorpick.png', None, 'Choose a background color', self.change_background_col),
          ('scshot', 'Take a Screenshot', 'images/tcamera.png', 'CTRL+I', 'Take a Screenshot of current view', self.take_screenshot),
          ('about', 'About pyNastran GUI', 'images/tabout.png', 'CTRL+H', 'About pyNastran GUI and help on shortcuts', self.about_dialog),
          ('creset', 'Reset camera view', 'images/trefresh.png', 'r', 'Reset the camera view to default', lambda: self._simulate_key_press('r'))]:
            if not "/" in icon:
                ico = QtGui.QIcon.fromTheme(icon)
            else:
                ico = QtGui.QIcon()
                ico.addPixmap(QtGui.QPixmap(icon), QtGui.QIcon.Normal, QtGui.QIcon.Off)
            actions[nam] = QtGui.QAction(ico, txt, self)
            if short:
                actions[nam].setShortcut(short)
            if tip:
                actions[nam].setStatusTip(tip)
            if func:
                actions[nam].triggered.connect(func)
        
        actions['toolbar'] = self.toolbar.toggleViewAction()
        actions['toolbar'].setStatusTip("Show/Hide application toolbar")
        
        actions['logwidget'] = self.log_dock.toggleViewAction()
        actions['logwidget'].setStatusTip("Show/Hide application log")
        
      
        # populate menus and toolbar
        for menu, items in [(self.menu_file, ('open_bdf', '', 'exit')),
                           (self.menu_view,  ('scshot', '', 'wireframe', 'surface', 'creset', '', 'back_col')),
                           (self.menu_window,('toolbar', 'logwidget')),
                           (self.menu_help,  ('about',)),
                           (self.toolbar, ('open_bdf', 'wireframe', 'surface', 'creset', 'scshot', '', 'exit'))]:
            for i in items:
                if not i:
                    menu.addSeparator()
                else:
                    menu.addAction(actions[i] if isinstance(i, basestring) else i())
                

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
        self.load_nastran_geometry(None, None, False, False)

        self.addGeometry()
        self.addAltGeometry()
        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(*self.background_col)
        self.rend.ResetCamera()
        self._simulate_key_press('t') # change mouse style to trackball
        
   
    def load_bdf(self):
        """
        Load BDF file.
        """
        # getOpenFileName return QString and we want Python string
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open BDF file', self.last_dir,
        'Nastran BDF (*.bdf *.dat *.nas)'))
        
        if fname:
            self.last_dir = os.path.split(fname)[0]
            self.load_nastran_geometry(fname, self.last_dir, False, False)
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
                ("BMP", ['.bmp'], vtk.vtkBMPWriter), ('JPG', ['.jpg', '.jpeg'], 
                vtk.vtkJPEGWriter), ("TIFF", ['.tif', '.tiff'], vtk.vtkTIFFWriter)):
                                       
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
        
        @param key a key that VTK should be informed about, e.g. 't'
        """
        self.vtk_interactor._Iren.SetEventInformation(0, 0, 0, 0, key, 0, None)
        self.vtk_interactor._Iren.KeyPressEvent()
        self.vtk_interactor._Iren.CharEvent()
        
    def addGeometry(self):
        aQuadMapper = vtk.vtkDataSetMapper()
        aQuadMapper.SetInput(self.grid)
        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 0, 0)  # red
        self.rend.AddActor(geometryActor)

    def addAltGeometry(self):
        aQuadMapper = vtk.vtkDataSetMapper()
        aQuadMapper.SetInput(self.grid2)
        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 1, 0)  # green
        self.rend.AddActor(geometryActor)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()
        
    
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    QtGui.QApplication.setOrganizationName("pyNastran")
    QtGui.QApplication.setOrganizationDomain(pyNastran.__website__)
    QtGui.QApplication.setApplicationName("pyNastran")
    QtGui.QApplication.setApplicationVersion("0.6")
    
    window = MainWindow()
    sys.exit(app.exec_())
