## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
#!/usr/bin/python
# pylint: disable=C0103,C0111,W0612,R0904

import os
import wx
#import vtk
try:
    from PIL import Image
except ImportError:
    from Image import Image

import sys
#from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

#import pyNastran
import pyNastran.gui
from pyNastran.gui.guiPanel import Pan
from pyNastran.gui.gui_inputs import get_inputs
from pyNastran.utils.log import SimpleLogger

ID_SAVEAS = 803
ID_ABOUT = 3

ID_PLOT = 200

ID_SURFACE = 901
ID_WIREFRAME = 902
ID_HIDDEN = 903
ID_EDGES = 904

ID_CAMERA = 910

ID_BDF = 920
ID_OP2 = 921
ID_F06 = 922
ID_CART3D = 923
ID_LAWGS = 924
ID_PANAIR = 925

ID_EXPORT = 930


pkgPath = pyNastran.gui.__path__[0]
#print "pkgPath = %r" % pkgPath

from pyNastran.gui.formats import (NastranIO, Cart3dIO, LaWGS_IO, PanairIO,
    is_nastran, is_cart3d, is_panair, is_lawgs)


if '?' in pkgPath:
    iconPath = 'icons'
else:
    iconPath = os.path.join(pkgPath, 'icons')
#print "iconPath = %r" % iconPath

#------------------------------------------------------------------------------
# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)
#------------------------------------------------------------------------------


class AppFrame(wx.Frame):

    def __init__(self, is_edges=False, is_nodal=False, is_centroidal=False,
                 format=None, input=None, output=None, shots=None, magnify=4, rotation=None, debug=False):

        assert debug in [True, False], 'debug=%s' % debug
        if shots is None:
            shots = []

        wx.Frame.__init__(self, None, -1, size=wx.Size(800, 600),
                          title='pyNastran')
        self.log = SimpleLogger('debug')
        self.infile_name = None
        self.is_edges = is_edges
        self.is_nodal = is_nodal
        self.is_centroidal = is_centroidal
        self.magnify = magnify
        assert is_centroidal != is_nodal, "is_centroidal and is_nodal can't be the same and are set to \"%s\"" % is_nodal
        self.dirname = ''
        self.setupFrame()

        print('format=%r input=%r output=%r' % (format, input, output))
        if format is not None and format.lower() not in ['panair', 'cart3d', 'lawgs', 'nastran']:
            sys.exit('\n---invalid format=%r' % format)
        elif format and input is not None: 
            format = format.lower()
            dirname = os.path.dirname(input)
            inputbase = input

            if format=='panair' and is_panair:
                print("loading panair")
                self.frmPanel.load_panair_geometry(inputbase, dirname)
            elif format=='nastran' and is_nastran:
                print("loading nastran")
                self.frmPanel.load_nastran_geometry(inputbase, dirname)
            elif format=='cart3d' and is_cart3d:
                print("loading cart3d")
                self.frmPanel.load_cart3d_geometry(inputbase, dirname)
            elif format=='lawgs' and is_lawgs:
                print("loading lawgs")
                self.frmPanel.load_LaWGS_geometry(inputbase, dirname)
            else:
                pass
                sys.exit('\n---unsupported format=%r' % format)
            self.UpdateWindowName(input)
            self.frmPanel.Update()
        else:
            self.frmPanel.scalarBar.VisibilityOff()
            self.frmPanel.scalarBar.Modified()

        #self.frmPanel.Update()
        #if rotation:
            #if rotation == '-x':

        print 'shots = %r' % shots
        if shots:
        #for shot in shots:
            self.frmPanel.widget.onTakePicture(ID_CAMERA, shots)
            sys.exit('took screenshot %r' % shots)

    def setupFrame(self):
        """
        --------------------------------
        |        VERTICAL(VMAIN)       |
        |   -------------------------  |
        |   |                       |  |
        |   |        toolbar        |  |
        |   |                       |  |
        |   -------------------------  |
        |   |       HORIZ           |  |
        |   |         |  VERTICAL   |  |
        |   |         |             |  |
        |   |   GUI   |  sidewindow |  |
        |   |         |             |  |
        |   |         |             |  |
        |   |         |             |  |
        |   -------------------------  |
        |   |                       |  |
        |   |       statusbar       |  |
        |   |                       |  |
        |   -------------------------  |
        |                              |
        --------------------------------
        """
        # Must call before any event handler is referenced.
        self.eventsHandler = EventsHandler(self, is_nodal=self.is_nodal,
                                           is_centroidal=self.is_centroidal)

        self.frmPanel = Pan(self, is_edges=self.is_edges,
                            is_nodal=self.is_nodal,
                            is_centroidal=self.is_centroidal,
                            magnify=self.magnify,
                            log=self.log,
                            gui_parent=self,
                            size=(100, 200))

        self.buildMenuBar()
        self.buildToolBar()
        self.buildStatusBar()

        self.SetMenuBar(self.menubar)

        self.frmPanel.infile_name = self.infile_name 
        self.frmPanel.buildVTK(self.infile_name)

        windowName = self.frmPanel.getWindowName()
        self.SetTitle(windowName)
        #self.SetSize([600,600])
        #self.Centre()

        # Add them to sizer.
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.frmPanel.widget, 1, wx.EXPAND | wx.ALL, 1)

        # Add buttons in their own sizer
        if 0:
            self.redBtn = wx.Button(self.frmPanel, label='Red')
            self.greenBtn = wx.Button(self.frmPanel, label='Green')
            self.exitBtn = wx.Button(self.frmPanel, label='Exit')

            vRight = wx.BoxSizer(wx.VERTICAL)
            vRight.AddStretchSpacer()
            vRight.Add(self.greenBtn, proportion=0, flag=wx.EXPAND | wx.ALL,
                       border=5)
            vRight.Add(self.exitBtn, proportion=0, flag=wx.EXPAND | wx.ALL,
                       border=5)
            vRight.Add(self.redBtn, proportion=0, flag=wx.EXPAND | wx.ALL,
                       border=5)
            vRight.AddStretchSpacer()
            hbox.Add(vRight, 1, wx.EXPAND)

        if 0:
            tree = self.buildTree(self.frmPanel)
            vRight = wx.BoxSizer(wx.VERTICAL)
            vRight.AddStretchSpacer()
            vRight.Add(tree, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
            vRight.AddStretchSpacer()
            hbox.Add(vRight, 1, wx.EXPAND)

        # best guess at tree
        if 0:
            #self.tree  = wx.Button( self.frmPanel, label='Tree' )

            vRight = wx.BoxSizer(wx.VERTICAL)

            scroll = wx.ScrolledWindow(self, -1)
            panelRight = wx.Panel(scroll, -1)

            panelRight = wx.Panel(self, wx.EXPAND)
            tree = self.buildTree(panelRight)
            vRight.Add(tree, flag=wx.EXPAND | wx.ALL)
            #vRight.Add(tree, 1, wx.EXPAND)
            #hbox.Add(panel1, 1, wx.EXPAND)

            vRight.Add(scroll, 1, wx.EXPAND | wx.ALL)
            #panelRight.SetSizer(vRight)
            panelRight.Layout()

            hbox.Add(vRight, 1, wx.EXPAND)
            #hbox.Add(panelRight, 1, wx.EXPAND | wx.ALL)

            # SetSizer both sizers in the most senior control that has sizers in it.
            self.vMain = wx.BoxSizer(wx.VERTICAL | wx.EXPAND)
            self.vMain.Add(hbox, 1, wx.EXPAND, 5)

        #self.vMain.AddStretchSpacer()
        #self.vMain.Add(self.frmPanel.widget, 0, wx.EXPAND)
        #self.vMain.Add(self.toolbar1, 0, wx.EXPAND)
        #self.vMain.AddStretchSpacer()
        #self.vMain.Add(hbox, 0, wx.EXPAND|wx.ALL, 5)
        self.frmPanel.SetSizer(hbox)
        #self.frmPanel.SetSizer(self.vMain)
        self.frmPanel.Layout()
        #self.toolbar1.Realize()

        events = self.eventsHandler
        # Bind Controls
        #self.Bind(wx.EVT_RIGHT_DOWN, events.OnRightDown)

        # Bind View Menu
        self.Bind(wx.EVT_MENU, self.frmPanel.widget.onTakePicture, id=ID_CAMERA)
        self.Bind(wx.EVT_MENU, self.frmPanel.onSetToWireframe, id=ID_WIREFRAME)
        self.Bind(wx.EVT_MENU, self.frmPanel.onSetToSurface, id=ID_SURFACE)
        self.Bind(wx.EVT_MENU, self.frmPanel.onFlipEdges, id=ID_EDGES)

        #self.Bind(wx.EVT_MENU, self.frmPanel.onSetToFlatShading,    self.flatShading)
        #self.Bind(wx.EVT_MENU, self.frmPanel.onSetToGouraudShading, self.gouraudShading)
        #self.Bind(wx.EVT_MENU, self.frmPanel.onSetToPhongShading,   self.phongShading)

        self.Bind(wx.EVT_MENU, events.onBackgroundColor, self.bkgColorView)
        #self.Bind(wx.EVT_MENU, events.onToggleStatusBar, self.showStatusBar)
        #self.Bind(wx.EVT_MENU, events.onToggleToolBar, self.showToolBar)

        # Bind Help Menu
        self.Bind(wx.EVT_MENU, events.onAbout, id=ID_ABOUT)
    #end __init__

    def log_info(self, msg):
        print msg

    def log_debug(self, msg):
        print msg

    def buildTree(self, panel1):
        tree = wx.TreeCtrl(self, 1, wx.DefaultPosition, (-1, -1),
                           wx.TR_HIDE_ROOT | wx.TR_HAS_BUTTONS)
        root = tree.AddRoot('Programmer')
        OS = tree.AppendItem(root, 'Operating Systems')
        tree.AppendItem(OS, 'Linux')
        tree.AppendItem(OS, 'FreeBSD')
        tree.AppendItem(OS, 'OpenBSD')
        tree.AppendItem(OS, 'NetBSD')
        tree.AppendItem(OS, 'Solaris')
        pl = tree.AppendItem(root, 'Programming Languages')
        cl = tree.AppendItem(pl, 'Compiled languages')
        sl = tree.AppendItem(pl, 'Scripting languages')
        tree.AppendItem(cl, 'Java')
        tree.AppendItem(cl, 'C++')
        tree.AppendItem(cl, 'C')
        tree.AppendItem(cl, 'Pascal')
        tree.AppendItem(sl, 'Python')
        tree.AppendItem(sl, 'Ruby')
        tree.AppendItem(sl, 'Tcl')
        tree.AppendItem(sl, 'PHP')

        tk = tree.AppendItem(root, 'Toolkits')
        tree.AppendItem(tk, 'Qt')
        tree.AppendItem(tk, 'MFC')
        tree.AppendItem(tk, 'wxPython')
        tree.AppendItem(tk, 'GTK+')
        tree.AppendItem(tk, 'Swing')
        #self.Bind(wx.EVT_MENU, self.frmPanel.SetToWireframe, id=ID_WIREFRAME)
        #tree.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnSelChanged, id=1)
        return tree

    def OnSelChanged(self, event):
        item = event.GetItem()
        self.display.SetLabel(tree.GetItemText(item))

    def UpdateWindowName(self, infile_name):
        self.infile_name = infile_name
        self.frmPanel.infile_name = infile_name
        windowName = self.frmPanel.getWindowName()
        self.SetTitle(windowName)

    def buildStatusBar(self):
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetStatusText('Ready')

    def buildToolBar(self):
        events = self.eventsHandler

        #toolbar1.AddSeparator()
        #toolbar1.AddSeparator()
        #tnew  = toolbar1.AddLabelTool(wx.ID_ANY,  '', wx.Bitmap(os.path.join(iconPath,'new.png')))
        #tsave = toolbar1.AddLabelTool(ID_SAVEAS,  '', wx.Bitmap(os.path.join(iconPath,'tsave.png')))
        #tundo = toolbar1.AddLabelTool(wx.ID_UNDO, '', wx.Bitmap(os.path.join(iconPath,'tundo.png')))
        #tredo = toolbar1.AddLabelTool(wx.ID_REDO, '', wx.Bitmap(os.path.join(iconPath,'tredo.png')))

        # toolbar at top - toggles
        toolbar1 = self.CreateToolBar()
        #toolbar.AddLabelTool(self.id, '', bitmap, wx.NullBitmap, self.kind,
        #                     shortHelp=wx.MenuItem.GetLabelFromText(self.menuText),
        #             longHelp=self.helpText)

        topen = os.path.join(iconPath, 'topen.png')
        assert os.path.exists(topen), 'topen=%r' % topen

        topen = wx.Image(topen, wx.BITMAP_TYPE_ANY)
        topen = toolbar1.AddLabelTool(ID_BDF, '', wx.BitmapFromImage(topen), longHelp='Loads a BDF')

        twireframe = wx.Image(os.path.join(iconPath, 'twireframe.png'), wx.BITMAP_TYPE_ANY)
        wireframe = toolbar1.AddLabelTool(ID_WIREFRAME, '', wx.BitmapFromImage(twireframe), longHelp='Set to Wireframe Model')

        tsolid = wx.Image(os.path.join(iconPath, 'tsolid.png'), wx.BITMAP_TYPE_ANY)
        surface = toolbar1.AddLabelTool(ID_SURFACE, '', wx.BitmapFromImage(tsolid), longHelp='Set to Surface/Solid Model')

        tedges = wx.Image(os.path.join(iconPath, 'tedges.png'), wx.BITMAP_TYPE_ANY)
        edges = toolbar1.AddLabelTool(ID_EDGES, '', wx.BitmapFromImage(tedges), longHelp='Show/Hide the edges')

        tcamera = wx.Image(os.path.join(iconPath, 'tcamera.png'), wx.BITMAP_TYPE_ANY)
        camera = toolbar1.AddLabelTool(ID_CAMERA, '', wx.BitmapFromImage(tcamera), longHelp='Take a Screenshot')

        texit = wx.Image(os.path.join(iconPath, 'texit.png'), wx.BITMAP_TYPE_ANY)
        etool = toolbar1.AddLabelTool(wx.ID_EXIT, '', wx.BitmapFromImage(texit), longHelp='Exit pyNastran GUI')
        #toolbar1.EnableTool(wx.ID_REDO, False)
        toolbar1.Realize()

        self.toolbar1 = toolbar1

        # Bind File Menu
        if is_nastran:
            self.Bind(wx.EVT_TOOL, events.onLoadBDF, id=ID_BDF)
            self.Bind(wx.EVT_TOOL, events.onLoadOP2, id=ID_OP2)
        self.Bind(wx.EVT_TOOL, events.onLoadCart3d, id=ID_CART3D)
        self.Bind(wx.EVT_TOOL, events.onLoadLaWGS, id=ID_LAWGS)
        self.Bind(wx.EVT_TOOL, events.onLoadPanair, id=ID_PANAIR)
        #self.Bind(wx.EVT_TOOL, events.onExport, id=ID_EXPORT)

        self.Bind(wx.EVT_MENU, events.onExit, id=wx.ID_EXIT)
       #self.Bind(wx.EVT_TOOL, events.onExit, id=wx.ID_EXIT)

        self.Bind(wx.EVT_MENU, self.frmPanel.onSetToWireframe, id=ID_WIREFRAME)
        self.Bind(wx.EVT_MENU, self.frmPanel.onSetToSurface, id=ID_SURFACE)
        self.Bind(wx.EVT_MENU, self.frmPanel.onSetToSurface, id=ID_CAMERA)

        #self.Bind(wx.EVT_TOOL, events.onSaveAsFile, id=ID_SAVEAS)
        #self.Bind(wx.EVT_TOOL, events.onUndo, tundo)
        #self.Bind(wx.EVT_TOOL, events.onRedo, tredo)

    def buildMenuBar(self):
        events = self.eventsHandler

        menubar = wx.MenuBar()
        # --------- File Menu -------------------------------------------------
        assert os.path.exists(os.path.join(iconPath, 'topen.png'))

        fileMenu = wx.Menu()
        #fileMenu.Append(wx.ID_NEW,  '&New','does nothing')
        if is_nastran:
            loadBDF = fileMenu.Append(ID_BDF, 'Load &BDF', 'Loads a BDF Input File')
            loadOP2 = fileMenu.Append(ID_OP2, 'Load O&P2', 'Loads an OP2 Results File')
            png = wx.Image(os.path.join(iconPath, 'topen.png'), wx.BITMAP_TYPE_PNG)
            loadBDF.SetBitmap(png.ConvertToBitmap())

        if is_cart3d:
            loadCart3d = fileMenu.Append(ID_CART3D, 'Load &Cart3D', 'Loads a Cart3D Input/Results File')
        if is_lawgs:
            loadLaWGS = fileMenu.Append(ID_LAWGS, 'Load &LaWGS', 'Loads an LaWGS File')
        if is_panair:
            loadPanair = fileMenu.Append(ID_PANAIR, 'Load &Panair', 'Loads a Panair Input File')
        #export     = fileMenu.Append(ID_EXPORT,'Export to...', 'Export the Model to Another Format')
        #print "topen = ",os.path.join(iconPath,'topen.png')
        sys.stdout.flush()

        #fileMenu.Append(wx.ID_RES, 'Load OP2 &Results','Loads a OP2 - does nothing')
        #fileMenu.Append(wx.ID_SAVE, '&Save','does nothing')
        fileMenu.AppendSeparator()

        # dummy import submenu
        #imp = wx.Menu()
        #imp.Append(wx.ID_ANY, 'Import newsfeed list...')
        #imp.Append(wx.ID_ANY, 'Import bookmarks...')
        #imp.Append(wx.ID_ANY, 'Import mail...')
        #fileMenu.AppendMenu(wx.ID_ANY, 'I&mport', imp)
        exitButton = wx.MenuItem(fileMenu,
                                 wx.ID_EXIT, 'Exit', 'Exits pyNastran')
        exitButton.SetBitmap(wx.Image(os.path.join(iconPath, 'texit.png'),
                                      wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        fileMenu.AppendItem(exitButton)

        # --------- View Menu -------------------------------------------------
        # status bar at bottom - toggles
        viewMenu = wx.Menu()
        camera = viewMenu.Append(ID_CAMERA, 'Take a Screenshot', 'Take a Screenshot')
        viewMenu.AppendSeparator()
        wireframe = viewMenu.Append(ID_WIREFRAME, 'Wireframe Model', 'Show Model as a Wireframe Model')
        surface = viewMenu.Append(ID_SURFACE, 'Surface Model', 'Show Model as a Surface Model')
        #viewMenu.AppendSeparator()

        #self.flatShading    = viewMenu.Append(wx.ID_ANY, 'Flat Shading',           'Flat Shading')
        #self.gouraudShading = viewMenu.Append(wx.ID_ANY, 'Mid (Gouraud) Shading',  'Mid (Gouraud) Shading')
        #self.phongShading   = viewMenu.Append(wx.ID_ANY, 'Smooth (Phong) Shading', 'Smooth (Phong) Shading')

        # view images
        wireframe.SetBitmap(wx.Image(os.path.join(iconPath, 'twireframe.png'),
                                     wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        surface.SetBitmap(wx.Image(os.path.join(iconPath, 'tsolid.png'),
                                   wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        camera.SetBitmap(wx.Image(os.path.join(iconPath, 'tcamera.png'),
                                  wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        #wireframe = toolbar1.AddLabelTool(ID_WIREFRAME, 'Set to Wireframe Model', wx.Bitmap(os.path.join(iconPath,'twireframe.png')))
        #surface   = toolbar1.AddLabelTool(ID_SURFACE,   'Set to Surface Model',   wx.Bitmap(os.path.join(iconPath,'tsolid.png')))
        #camera    = toolbar1.AddLabelTool(ID_CAMERA,    'Take a Screenshot',      wx.Bitmap(os.path.join(iconPath,'tcamera.png')))

        viewMenu.AppendSeparator()
        self.bkgColorView = viewMenu.Append(wx.ID_ANY,
                                            'Change Background Color',
                                            'Change Background Color')
        #self.showStatusBar = viewMenu.Append(wx.ID_ANY, 'Show statusbar', 'Show Statusbar', kind=wx.ITEM_CHECK)
        #self.showToolBar   = viewMenu.Append(wx.ID_ANY, 'Show toolbar',   'Show Toolbar',   kind=wx.ITEM_CHECK)
        #viewMenu.Check(self.showStatusBar.GetId(), True)
        #viewMenu.Check(self.showToolBar.GetId(),   True)

        # --------- Plot Menu -------------------------------------------------
        #plotMenu = wx.Menu()
        #plot = plotMenu.Append(ID_PLOT, '&Plot Data','Plot Data')
        #self.Bind(wx.EVT_MENU, self.onPlot, id=ID_PLOT)

        # --------- Help / About Menu -----------------------------------------
        # help/about menu
        helpMenu = wx.Menu()
        helpMenu.Append(ID_ABOUT, '&About', 'About pyNastran')

        # menu bar
        menubar.Append(fileMenu, '&File')
        #menubar.Append(plotMenu, '&Plot')
        menubar.Append(viewMenu, '&View')
        menubar.Append(helpMenu, '&Help')
        self.menubar = menubar

    def onPlot(self, event):
        #e = Example(self)
        e = TestFrame(self)
        e.Show()
#end AppFrame class

#------------------------------------------------------------------------------


class TestFrame(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, parent, -1, "GridBagSizer Test")
        sizer = wx.GridBagSizer(hgap=5, vgap=5)

        labels = "1 2 3 4 5 6 7 8 9 0".split()
        for col in xrange(3):
            for row in xrange(3):
                bw = wx.Button(self, label=labels[row * 3 + col])
                sizer.Add(bw, pos=(row, col))

        bw = wx.Button(self, label="span 3 rows")
        sizer.Add(bw, pos=(0, 3), span=(3, 1), flag=wx.EXPAND)

        bw = wx.Button(self, label="span all columns")
        sizer.Add(bw, pos=(3, 0), span=(1, 4), flag=wx.EXPAND)

        sizer.AddGrowableCol(3)
        sizer.AddGrowableRow(3)

        self.SetSizer(sizer)
        self.Fit()


class Example(wx.Frame):

    def __init__(self, *args, **kwargs):
        super(Example, self).__init__(*args, **kwargs)

        self.InitUI()

    def InitUI(self):

        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        fitem = fileMenu.Append(wx.ID_EXIT, 'Quit', 'Quit application')
        menubar.Append(fileMenu, '&File')
        self.SetMenuBar(menubar)

        self.Bind(wx.EVT_MENU, self.OnQuit, fitem)

        self.SetSize((300, 200))
        self.SetTitle('Simple menu')
        self.Centre()
        self.Show(True)

    def onQuit(self, e):
        self.Close()


class EventsHandler(object):

    def __init__(self, parent, is_nodal=False, is_centroidal=True):
        self.parent = parent
        self.is_nodal = is_nodal
        self.is_centroidal = is_centroidal

    # File Menu

    def onLoadBDF(self, event):
        """ Open a file"""
        #print "OnOpen..."

        wildcard = "Nastran BDF (*.bdf; *.dat; *.nas)|*.bdf;*.dat;*.nas|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Nastran Input Deck to Load'
        loadFunction = self.parent.frmPanel.load_nastran_geometry
        self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)

    def onLoadOP2(self, event):
        """ Open a file"""
        #print "OnOpen..."

        if 0:
            bdf = self.parent.infile_name
            bdfBase = os.path.basename(bdf)
            #dirname = os.path.dirname(bdf)
            (op2name, op2) = os.path.splitext(bdfBase)
            op2 = os.path.join(self.parent.dirname, op2name + '.op2')

            self.parent.op2FileName = op2
            if os.path.exists(op2):
                self.parent.frmPanel.load_nastran_results(op2)
                self.parent.frmPanel.Update()
            return

        wildcard = "Nastran OP2 (*.op2)|*.op2|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Nastran Output File to Load (OP2 only)'
        loadFunction = self.parent.frmPanel.load_nastran_results
        self.createLoadFileDialog(wildcard, Title, loadFunction)

    def onLoadCart3d(self, event):
        """ Open a file"""
        #print "OnOpen..."

        #wildcard = "Cart3d (*.i.tri;*.a.tri)|(*.i.triq;*.a.triq)|" \
        # "All files (*.*)|*.*"

        wildcard = "Cart3d (*.i.tri; *.a.tri; *.triq)|*.i.tri;*.a.tri;*.triq|" \
            "All files (*.*)|*.*"
        #wildcard = "Cart3d (*.tri;)|(*.triq;)|" \
        # "All files (*.*)|*.*"

        Title = 'Choose a Cart3d Input File to Load'
        loadFunction = self.parent.frmPanel.load_cart3d_geometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)

    def onLoadLaWGS(self, event):
        """ Open a file"""
        wildcard = "LaWGS (*.wgs)|*.wgs|" \
            "All files (*.*)|*.*"

        Title = 'Choose a LaWGS File to Load'
        loadFunction = self.parent.frmPanel.load_LaWGS_geometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)

    def onLoadPanair(self, event):
        """ Open a file"""
        wildcard = "Panair (*.inp)|*.inp|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Panair Input File to Load'
        loadFunction = self.parent.frmPanel.load_panair_geometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)

    def createLoadFileDialog(self, wildcard, Title, loadFunction, updateWindowName=False):
        dlg = wx.FileDialog(None, Title, self.parent.dirname, "",
                            wildcard, wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            fileName = dlg.GetFilename()
            self.parent.dirname = dlg.GetDirectory()
            fname = os.path.join(self.parent.dirname, fileName)
            print("fname = ", fname)
            if updateWindowName:
                self.parent.UpdateWindowName(fileName)
            #self.parent.frmPanel.loadCart3dGeometry(fname,self.parent.dirname)
            loadFunction(fname, self.parent.dirname)
            self.parent.frmPanel.Update()
        dlg.Destroy()

    def onExport(self):
        if self.modelType == 'nastran':
            pass

        if self.modelType == 'cart3d':
            exportToNastran(self, fname, points, elements, regions)

        wildcard = ("Nastran OP2 (*.op2)|*.op2|"
                    "All files (*.*)|*.*")

        Title = 'Choose a Nastran Output File to Load (OP2 only)'
        loadFunction = self.parent.frmPanel.load_nastran_results
        self.createLoadFileDialog(wildcard, Title, loadFunction)

    def onExit(self, event):
        self.parent.Destroy()

    # View Menu
    def onBackgroundColor(self, event):
        dlg = wx.ColourDialog(self.parent)
        dlg.GetColourData().SetChooseFull(True)
        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetColourData()
            #print 'You selected: %s\n' % str(data.GetColour().Get())
            dlg.Destroy()
            self.SetColor(data.GetColour().Get())
            return
        dlg.Destroy()

    def SetColor(self, bkgColor):
        """
        @warning if the model is not loaded and the color is not "clicked"
        on, then rend will return None
        """
        rend = self.parent.frmPanel.widget.GetCurrentRenderer()
        if not rend:  # rend doesnt exist at first
            print 'Try again to change the color (you didnt "click" on the color; one time bug)'
            return
            rend = self.parent.frmPanel.widget.GetCurrentRenderer()

        ## bkgColor range from 0 to 255
        ## color ranges from 0 to 1
        color = [bkgColor[0] / 255., bkgColor[1] / 255., bkgColor[2] / 255.]
        rend.SetBackground(color)
        self.parent.frmPanel.widget.Render()

    def onToggleStatusBar(self, e):
        if self.parent.showStatusBar.IsChecked():
            self.parent.statusbar.Show()
        else:
            self.parent.statusbar.Hide()

    def onToggleToolBar(self, e):
        if self.parent.showToolBar.IsChecked():
            self.parent.toolbar1.Show()
        else:
            self.parent.toolbar1.Hide()

    # Help Menu
    def onAbout(self, event):
        about = [
            'pyNastran v%s' % (pyNastran.__version__),
            'Copyright %s; Steven Doyle 2011-2012\n' % (pyNastran.__license__),
            '%s' % (pyNastran.__website__),
            '',
            'Mouse',
            'Left Click - Rotate',
            'Middle Click - Pan/Recenter Rotation Point',
            'Shift + Left Click - Pan/Recenter Rotation Point',
            'Right Mouse - Zoom',
            '',
            'Keyboard Controls',
            'X/x - snap to x axis',
            'Y/y - snap to y axis',
            'Z/z - snap to z axis',
            '',
            'h   - show/hide legend & info',
            'i   - take a screenshot (image)',
            'L   - cycle op2 results',
            'm/M - scale up/scale down by 1.1 times',
            'o/O - rotate counter-clockwise/clockwise 5 degrees',
            's   - surface',
            'w   - wireframe',
        ]

        # not done
              #'',
              #'left arrow  - pan left',
              #'right arrow - pan right',
              #'up arrow    - pan up',
              #'down arrow  - pan down',
              #'',
              #'e  - show/hide edges',
              #'f   fly to rotation point',
              #'p   project point',

        dlg = wx.MessageDialog(None, '\n'.join(about), 'About',
                               wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

#end Events class

#------------------------------------------------------------------------------

def main():
    inputs = get_inputs()
    app = wx.App(redirect=False)
    appFrm = AppFrame(*inputs)
    #appFrm.Show()
    print("launching gui")
    app.MainLoop()

#end class

#==============================================================================

if __name__ == '__main__':
    main()
