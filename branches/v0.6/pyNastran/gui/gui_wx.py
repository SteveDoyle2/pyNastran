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

import pyNastran.gui
from pyNastran.gui.pan import Pan
from pyNastran.gui.arg_handling import get_inputs
from pyNastran.utils.log import SimpleLogger

# === File ===
#ID_SAVEAS = 803

ID_GEOMETRY = 201
ID_RESULTS = 202

ID_BDF = 920
ID_OP2 = 921
ID_F06 = 922
ID_CART3D = 923
ID_LAWGS = 924
ID_PANAIR = 925
ID_STL = 926
ID_TETGEN = 927
ID_USM3D = 928
ID_PLOT3D = 929

ID_EXPORT = 930

# === View ===
ID_SURFACE = 901
ID_WIREFRAME = 902
ID_HIDDEN = 903

ID_EDGES = 904
ID_CAMERA = 910

# === Help ===
ID_ABOUT = 3
#ID_PLOT = 200  # 2D grid




from pyNastran.gui.formats import (NastranIO, Cart3dIO, PanairIO, LaWGS_IO, STL_IO, TetgenIO, Usm3dIO, Plot3d_io,
    is_nastran, is_cart3d, is_panair, is_lawgs, is_stl, is_tetgen, is_usm3d, is_plot3d)
valid_formats = [
    'cart3d', 'lawgs', 'nastran', 'panair', 'plot3d', 'stl', 'tetgen', 'usm3d'
]

assert is_nastran == True, is_nastran
print('is_nastran = %r' % is_nastran)

pkgPath = pyNastran.gui.__path__[0]
#print "pkgPath = %r" % pkgPath

if '?' in pkgPath:
    iconPath = 'icons'
else:
    iconPath = os.path.join(pkgPath, 'icons')
#print("iconPath = %r" % iconPath)

#------------------------------------------------------------------------------
# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)
#------------------------------------------------------------------------------


class AppFrame(wx.Frame):

    def __init__(self, inputs):
        wx.Frame.__init__(self, None, -1, size=wx.Size(800, 600),
                          title='pyNastran')

        self.format = None
        debug = inputs['debug']
        assert debug in [True, False], 'debug=%s' % debug
        shots = inputs['shots']
        if shots is None:
            shots = []

        #self.show_info = True
        #self.show_debug = True
        #self.show_gui = True
        #self.show_command = True

        self.log = SimpleLogger('debug')
        self.infile_name = None
        self.is_edges = inputs['is_edges']
        self.is_nodal = inputs['is_nodal']
        self.is_centroidal = inputs['is_centroidal']
        self.magnify = inputs['magnify']
        assert self.is_centroidal != self.is_nodal, "is_centroidal and is_nodal can't be the same and are set to \"%s\"" % self.is_nodal
        self.dirname = ''
        self.setupFrame()
        self.load_batch_inputs(inputs)

        #self.frmPanel.Update()
        #if rotation:
            #if rotation == '-x':

        print 'shots = %r' % shots
        if shots:
        #for shot in shots:
            self.frmPanel.widget.onTakePicture(ID_CAMERA, shots)
            sys.exit('took screenshot %r' % shots)
        #print("*format =", self.format)

    def load_batch_inputs(self, inputs):
        format = inputs['format']
        input = inputs['input']
        output = inputs['output']
        print('format=%r input=%r output=%r' % (format, input, output))
        if format is not None and format.lower() not in valid_formats:
            msg = '\n---invalid format=%r' % format
            print msg
            self.frmPanel.scalarBar.VisibilityOff()
            self.frmPanel.scalarBar.Modified()
            return
            #raise IOError(msg)
            #sys.exit(msg)
        elif format and input is not None:
            format = format.lower()
            dirname = os.path.dirname(input)
            inputbase = input

            if not os.path.exists(input):
                msg = 'input=%r does not exist' % input
                print msg
                self.frmPanel.scalarBar.VisibilityOff()
                self.frmPanel.scalarBar.Modified()
                return
                #raise IOError(msg)

            if format=='panair' and is_panair:
                print("loading panair")
                self.frmPanel.load_panair_geometry(inputbase, dirname)
            elif format=='nastran' and is_nastran:
                print("loading nastran")
                self.frmPanel.load_nastran_geometry(inputbase, dirname)
                load_results = self.frmPanel.load_nastran_results
            elif format=='cart3d' and is_cart3d:
                print("loading cart3d")
                self.frmPanel.load_cart3d_geometry(inputbase, dirname)
            elif format=='lawgs' and is_lawgs:
                print("loading lawgs")
                self.frmPanel.load_LaWGS_geometry(inputbase, dirname)
            elif format=='stl' and is_stl:
                print("loading stl")
                self.frmPanel.load_stl_geometry(inputbase, dirname)
            elif format=='tetgen' and is_tetgen:
                print("loading tetgen")
                self.frmPanel.load_tetgen_geometry(inputbase, dirname)
            elif format=='usm3d' and is_usm3d:
                print("loading usm3d")
                self.frmPanel.load_usm3d_geometry(inputbase, dirname)
                load_results = self.frmPanel.load_usm3d_results
            elif format=='plot3d' and is_plot3d:
                print("loading plot3d")
                self.frmPanel.load_plot3d_geometry(inputbase, dirname)
            else:
                msg = '\n---unsupported format=%r' % format
                print msg
                return
                #raise IOError(msg)
                #sys.exit(msg)

            self.UpdateWindowName(input)

            if output:
                if os.path.exists(output):
                    print "format=%r" % format
                    print "output=%r" % output
                    load_results(output, dirname)
                else:
                    msg = 'output=%r does not exist' % output
                    #break
                    #self.frmPanel.Update()
                    raise IOError(msg)
                    #return

            self.format = format
            print("set the format=%r" % self.format)
            self.frmPanel.Update()
        else:
            self.frmPanel.scalarBar.VisibilityOff()
            self.frmPanel.scalarBar.Modified()


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
        #self.buildToolBar()
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
        #self.Bind(wx.EVT_MENU, self.frmPanel.onSetToFlatShading,    self.flatShading)
        #self.Bind(wx.EVT_MENU, self.frmPanel.onSetToGouraudShading, self.gouraudShading)
        #self.Bind(wx.EVT_MENU, self.frmPanel.onSetToPhongShading,   self.phongShading)

        #self.Bind(wx.EVT_MENU, events.onToggleStatusBar, self.showStatusBar)
        #self.Bind(wx.EVT_MENU, events.onToggleToolBar, self.showToolBar)

        # Bind Help Menu
        self.Bind(wx.EVT_MENU, events.onAbout, id=ID_ABOUT)
    #end __init__

    def log_info(self, msg):
        #if self.log_info:
        print msg

    def log_debug(self, msg):
        #if self.log_debug:
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
        #self.Bind(wx.EVT_MENU, self.frmPanel.onSetToWireframe, id=ID_WIREFRAME)
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
        return

        #tnew  = toolbar1.AddLabelTool(wx.ID_ANY,  '', wx.Bitmap(os.path.join(iconPath,'new.png')))
        #tsave = toolbar1.AddLabelTool(ID_SAVEAS,  '', wx.Bitmap(os.path.join(iconPath,'tsave.png')))
        #tundo = toolbar1.AddLabelTool(wx.ID_UNDO, '', wx.Bitmap(os.path.join(iconPath,'tundo.png')))
        #tredo = toolbar1.AddLabelTool(wx.ID_REDO, '', wx.Bitmap(os.path.join(iconPath,'tredo.png')))

        # toolbar at top - toggles
        toolbar1 = self.CreateToolBar()
        #toolbar.AddLabelTool(self.id, '', bitmap, wx.NullBitmap, self.kind,
        #                     shortHelp=wx.MenuItem.GetLabelFromText(self.menuText),
        #             longHelp=self.helpText)

        topen     = toolbar1.AddLabelTool(ID_GEOMETRY,  '', wx.BitmapFromImage(topen),      longHelp='Loads Geometry')
        wireframe = toolbar1.AddLabelTool(ID_WIREFRAME, '', wx.BitmapFromImage(twireframe), longHelp='Set to Wireframe Model')
        surface   = toolbar1.AddLabelTool(ID_SURFACE,   '', wx.BitmapFromImage(tsolid),     longHelp='Set to Surface/Solid Model')
        edges     = toolbar1.AddLabelTool(ID_EDGES,     '', wx.BitmapFromImage(tedges),     longHelp='Show/Hide the edges')
        camera    = toolbar1.AddLabelTool(ID_CAMERA,    '', wx.BitmapFromImage(tcamera),    longHelp='Take a Screenshot')
        etool     = toolbar1.AddLabelTool(wx.ID_EXIT,   '', wx.BitmapFromImage(texit),      longHelp='Exit pyNastran GUI')
        #toolbar1.EnableTool(wx.ID_REDO, False)

        # Bindings
        #self.Bind(wx.EVT_TOOL, events.onExport, id=ID_EXPORT)
        #self.Bind(wx.EVT_TOOL, events.onSaveAsFile, id=ID_SAVEAS)
        #self.Bind(wx.EVT_TOOL, events.onUndo, tundo)
        #self.Bind(wx.EVT_TOOL, events.onRedo, tredo)

    def buildMenuBar(self):
        events = self.eventsHandler

        menubar = wx.MenuBar()
        # --------- File Menu -------------------------------------------------

        formats = []
        if is_nastran: formats.append('nastran')
        if is_cart3d: formats.append('cart3d')
        if is_lawgs:  formats.append('lawgs')
        if is_panair: formats.append('panair')
        if is_plot3d: formats.append('plot3d')
        if is_stl:    formats.append('stl')
        if is_tetgen: formats.append('tetgen')
        if is_usm3d:  formats.append('usm3d')
        if formats: formats.append('')

        buttons = {
            # ---File---
            #button_name    : [ID,      selfattr, icon_path,  menu_name, tip, event]
            'load_geometry' : [ID_GEOMETRY, None, os.path.join(iconPath, 'topen.png'), 'Load &Geometry', 'Load the Geometry...', events.onLoadGeometry],
            'load_results'  : [ID_RESULTS,  None, None, 'Load &Results', 'Load the Results...', None],

            'nastran'    : [ID_BDF,     None, None, 'Load &Nastran', 'Loads a Nasran Input File',    events.onLoadBDF],
            'cart3d'     : [ID_CART3D,  None, None, 'Load &Cart3D',  'Loads a Cart3D Input/Results', events.onLoadCart3d],
            'lawgs'      : [ID_LAWGS,   None, None, 'Load &LaWGS',   'Loads an LaWGS File',          events.onLoadLaWGS],
            'panair'     : [ID_PANAIR,  None, None, 'Load &Panair',  'Loads a Panair Input File',    events.onLoadPanair],
            'plot3d'     : [ID_PLOT3D,  None, None, 'Load Pl&ot3d',  'Loads a Plot3d Input File',    events.onLoadPlot3d],
            'stl'        : [ID_STL,     None, None, 'Load &STL',     'Loads a STL Input File',       events.onLoadSTL],
            'usm3d'      : [ID_USM3D,   None, None, 'Load &Usm3d',   'Loads a Usm3d Input File',     events.onLoadUsm3d],
            'tetgen'     : [ID_TETGEN,  None, None, 'Load &Tetgen',  'Loads a Tetgen Input File',    events.onLoadTetgen],

            'exit'       : [wx.ID_EXIT, None, os.path.join(iconPath, 'texit.png'),'Exit', 'Exits pyNastran', events.onExit],

            'export'    : [ID_EXPORT,  None, None,  'Export to...' 'Export the Model to Another Format'],
            'Save'      : [wx.ID_SAVE, None, None, '&Save', 'does nothing'],
            
            # ---View---
            'camera'     : [ID_CAMERA,    None,           os.path.join(iconPath, 'tcamera.png'),    'Take a Screenshot', 'Take a Screenshot', self.frmPanel.widget.onTakePicture],
            'wireframe'  : [ID_WIREFRAME, None,           os.path.join(iconPath, 'twireframe.png'), 'Wireframe Model', 'Show Model as a Wireframe Model', self.frmPanel.onSetToWireframe],
            'edges'      : [ID_EDGES,     None,           os.path.join(iconPath, 'tedges.png'),     'Show/Hide Edges', 'Show/Hide Edges', self.frmPanel.onFlipEdges],
            'surface'    : [ID_SURFACE,   None,           os.path.join(iconPath, 'tsolid.png'),     'Surface Model', 'Show Model as a Surface Model', self.frmPanel.onSetToSurface],
            'color'      : [wx.ID_ANY,    'bkgColorView', os.path.join(iconPath, 'tcolorpick.png'), 'Change Background Color', 'Change Background Color', events.onBackgroundColor],
            
            # ---Help---
            'about' : [ID_ABOUT, None, os.path.join(iconPath, 'tabout.png'), '&About', 'About pyNastran', events.onAbout],
        }
        
        toolbar_buttons = ['load_geometry', 'camera', 'wireframe', 'surface', 'edges', '', 'exit']
        menu_order = ['&File', '&View', '&Help']
        menus = {
            '&File' : formats + ['load_geometry', 'load_results', '', 'exit'],
            '&View' : ['camera', '', 'wireframe', 'surface', 'edges', '', 'color'],
            '&Help' : ['about'],
        }
        
        self.toolbar1 = self.CreateToolBar()
        for button_name in toolbar_buttons:
            if button_name == '':
                self.toolbar1.AddSeparator()
            elif button_name is not None:
                button_args = buttons[button_name]
                (ID, selfattr, icon_pth, desc, tip, event) = button_args
                if os.path.exists(icon_pth):
                    timage = wx.Image(icon_pth, wx.BITMAP_TYPE_ANY)
                    button = self.toolbar1.AddLabelTool(ID, '', wx.BitmapFromImage(timage), longHelp=tip)
                    self.Bind(wx.EVT_TOOL, event, id=ID)
                else:
                    print "canot add %s" % desc
            else:
                print "canot add %s" % desc
        self.toolbar1.Realize()

        print('---building button menu---')
        for menu_name in menu_order:
            button_names = menus[menu_name]

            menu = wx.Menu()
            for button_name in button_names:
                if button_name == '':
                    menu.AppendSeparator()
                elif button_name in button_names:
                    button = None
                    try:
                        button_args = buttons[button_name]
                    except KeyError:
                        raise RuntimeError('menu_name=%r has not been defined.' % button_name)

                    (ID, selfattr, icon_pth, desc, tip, event) = button_args
                    if icon_pth is None:  # no icon
                        button = menu.Append(ID, desc, tip)
                    elif os.path.exists(icon_pth):
                        button = wx.MenuItem(menu, ID, desc, tip)
                        button.SetBitmap(wx.Image(icon_pth,
                                                  wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                        menu.AppendItem(button)
                    else:
                        print('***cannot find icon_pth=%r' % icon_pth)

                    if selfattr is not None:
                        print("button_name=%r; setting self.%s" % (button_name, selfattr))
                        setattr(self, selfattr, button)
                        attr = getattr(self, selfattr)
                        self.Bind(wx.EVT_MENU, event, attr)
                    else:
                        self.Bind(wx.EVT_MENU, event, id=ID)
                else:
                    print('***cannot find button=%r' % button_name)
            menubar.Append(menu, menu_name)
        print('---built button menu---')

        # --------- View Menu -------------------------------------------------
        #self.flatShading    = viewMenu.Append(wx.ID_ANY, 'Flat Shading',           'Flat Shading')
        #self.gouraudShading = viewMenu.Append(wx.ID_ANY, 'Mid (Gouraud) Shading',  'Mid (Gouraud) Shading')
        #self.phongShading   = viewMenu.Append(wx.ID_ANY, 'Smooth (Phong) Shading', 'Smooth (Phong) Shading')

        #self.showStatusBar = viewMenu.Append(wx.ID_ANY, 'Show statusbar', 'Show Statusbar', kind=wx.ITEM_CHECK)
        #self.showToolBar   = viewMenu.Append(wx.ID_ANY, 'Show toolbar',   'Show Toolbar',   kind=wx.ITEM_CHECK)
        #viewMenu.Check(self.showStatusBar.GetId(), True)
        #viewMenu.Check(self.showToolBar.GetId(),   True)

        # --------- Plot Menu -------------------------------------------------
        #plotMenu = wx.Menu()
        #plot = plotMenu.Append(ID_PLOT, '&Plot Data','Plot Data')
        #self.Bind(wx.EVT_MENU, self.onPlot, id=ID_PLOT)

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

    def onLoadGeometry(self, event):
        """ Open a file"""

        wildcard = ''
        if is_nastran:
            wildcard += "Nastran BDF (*.bdf; *.dat; *.nas)|*.bdf;*.dat;*.nas|"

        wildcard += "Cart3d (*.i.tri; *.a.tri; *.triq)|*.i.tri;*.a.tri;*.triq|"
        wildcard += "LaWGS (*.wgs)|*.wgs|"
        wildcard += "Panair (*.inp)|*.inp|"
        wildcard += "Plot3d (*.p3d; *.p3da; *.xyz)|*.p3d;*.p3da;*.xyz|"
        wildcard += "STL (*.STL)|*.STL|"
        wildcard += "Tetgen (*.smesh)|*.smesh|"
        wildcard += "Usm3D (*.cogsg; *.front)|*.cogsg;*.front|"

        Title = 'Choose an Input File to Load'
        loadFunction = self._load_geometry
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)

    def _load_geometry(self, fname, dirname):
        basename, ext = os.path.splitext(fname)
        ext = ext.lower()
        
        Format_to_load = {
            '.tri' : (self.parent.frmPanel.load_cart3d_geometry, 'cart3d'),
            '.triq' : (self.parent.frmPanel.load_cart3d_geometry, 'cart3d'),

            '.bdf' : (self.parent.frmPanel.load_nastran_geometry, 'nastran'),
            '.nas' : (self.parent.frmPanel.load_nastran_geometry, 'nastran'),
            '.dat' : (self.parent.frmPanel.load_nastran_geometry, 'nastran'),

            '.inp' : (self.parent.frmPanel.load_panair_geometry, 'panair'),

            '.p3d' : (self.parent.frmPanel.load_plot3d_geometry, 'plot3d'),
            '.p3da' : (self.parent.frmPanel.load_plot3d_geometry, 'plot3d'),
            '.xyz' : (self.parent.frmPanel.load_plot3d_geometry, 'plot3d'),

            '.wgs' : (self.parent.frmPanel.load_LaWGS_geometry, 'lawgs'),
            '.stl' : (self.parent.frmPanel.load_stl_geometry, 'stl'),
            '.smesh' : (self.parent.frmPanel.load_tetgen_geometry, 'tetgen'),

            '.cogsg' : (self.parent.frmPanel.load_usm3d_geometry, 'usm3d'),
            '.front' : (self.parent.frmPanel.load_usm3d_geometry, 'usm3d'),
        }
        if ext not in Format_to_load:
            self.Error('Cannot load %r because it has an invalid extension.' % str(fname), 'Error Loading Geometry')
            return False
        load_function, Format = Format_to_load[ext]
        load_function(fname, dirname)
        self.parent.format = Format
        #print("geo format =", Format)
        return True

    def onLoadBDF(self, event):
        """ Open a file"""
        #print "OnOpen..."

        wildcard = "Nastran BDF (*.bdf; *.dat; *.nas)|*.bdf;*.dat;*.nas|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Nastran Input Deck to Load'
        loadFunction = self.parent.frmPanel.load_nastran_geometry
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)
        if is_loaded:
            self.parent.format = 'nastran'

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
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)
        if is_loaded:
            pass

    def onLoadCart3d(self, event):
        """ Open a file"""
        wildcard = "Cart3d (*.i.tri; *.a.tri; *.triq)|*.i.tri;*.a.tri;*.triq|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Cart3d Input File to Load'
        loadFunction = self.parent.frmPanel.load_cart3d_geometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)
        if is_loaded:
            self.parent.format = 'cart3d'

    def onLoadLaWGS(self, event):
        """ Open a file"""
        wildcard = "LaWGS (*.wgs)|*.wgs|" \
            "All files (*.*)|*.*"

        Title = 'Choose a LaWGS File to Load'
        loadFunction = self.parent.frmPanel.load_LaWGS_geometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)
        if is_loaded:
            self.parent.format = 'lawgs'

    def onLoadPanair(self, event):
        """ Open a file"""
        wildcard = "Panair (*.inp)|*.inp|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Panair Input File to Load'
        loadFunction = self.parent.frmPanel.load_panair_geometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)
        if is_loaded:
            self.parent.format = 'panair'

    def onLoadSTL(self, event):
        """ Open a file"""
        wildcard = "STL (*.STL)|*.STL|" \
            "All files (*.*)|*.*"

        Title = 'Choose a STL File to Load'
        loadFunction = self.parent.frmPanel.load_stl_geometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)
        if is_loaded:
            self.parent.format = 'stl'

    def onLoadTetgen(self, event):
        """ Open a file"""
        wildcard = "Tetgen (*.smesh)|*.smesh|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Tetgen (smesh) File to Load'
        loadFunction = self.parent.frmPanel.load_tetgen_geometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)
        if is_loaded:
            self.parent.format = 'tetgen'

    def onLoadUsm3d(self, event):
        """ Open a file"""
        wildcard = "Usm3D (*.cogsg; *.front)|*.cogsg;*.front|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Usm3D (cogsg/front) Input File to Load'
        loadFunction = self.parent.frmPanel.load_usm3d_geometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)
        if is_loaded:
            self.parent.format = 'usm3d'

    def onLoadPlot3d(self, event):
        """ Open a file"""
        wildcard = "Plot3d (*.p3d; *.p3da; *.xyz)|*.p3d;*.p3da;*.xyz|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Plot3D Input File to Load'
        loadFunction = self.parent.frmPanel.load_usm3d_geometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction, updateWindowName=True)
        if is_loaded:
            self.parent.format = 'plot3d'

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
            is_loaded = True
        else:
            is_loaded = False
        dlg.Destroy()
        return is_loaded

    def Warn(self, message, caption='Warning!'):
        dlg = wx.MessageDialog(self.parent.frmPanel, message, caption, wx.OK | wx.ICON_WARNING)
        dlg.ShowModal()
        dlg.Destroy()

    def Error(self, message, caption='Error!'):
        dlg = wx.MessageDialog(self.parent.frmPanel, message, caption, wx.OK | wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()

    def onLoadResults(self, event):
        Format = self.parent.format
        if Format is None:
            self.Error('You must load an input file first.', 'No Input File Selected')
            return

        #print('res format', Format)
        Title = 'Choose a %s Output File to Load' % Format.title()
        if Format == 'nastran' and is_nastran:
            wildcard = "Nastran Results (*.op2; *.f06)|*.op2;*.f06|" \
                "All files (*.*)|*.*"
            load_function = self.parent.frmPanel.load_nastran_results
        elif Format == 'usm3d':
            wildcard = "Usm3d flow (*.flo)|*.flo|" \
                "All files (*.*)|*.*"
            load_function = self.parent.frmPanel.load_usm3d_results
        else:
            msg = 'Results for format=%r is not supported.' % Format
            #self.log.error(msg)
            self.Error(msg, 'Invalid Results Format Error')
            return

        is_loaded = self.createLoadFileDialog(wildcard, Title, load_function, updateWindowName=True)
        print("loaded the results...")
        #print("format =", Format)

    #==============================================
    def onExport(self):
        if self.modelType == 'nastran':
            pass
        elif self.modelType == 'cart3d':
            exportToNastran(self, fname, points, elements, regions)

        wildcard = ("Nastran OP2 (*.op2)|*.op2|"
                    "All files (*.*)|*.*")

        Title = 'Choose a Nastran Output File to Load (OP2 only)'
        loadFunction = self.parent.frmPanel.load_nastran_results
        is_loaded = self.createLoadFileDialog(wildcard, Title, loadFunction)

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
            'pyNastran v%s' % pyNastran.__version__,
            '%s' % pyNastran.__copyright__,
            '%s\n' % pyNastran.__author__,

            '%s' % pyNastran.__website__,
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
    inputs = get_inputs('wx')
    app = wx.App(redirect=False)
    appFrm = AppFrame(inputs)
    #appFrm.Show()
    print("launching gui")
    app.MainLoop()

#end class

#==============================================================================

if __name__ == '__main__':
    main()
