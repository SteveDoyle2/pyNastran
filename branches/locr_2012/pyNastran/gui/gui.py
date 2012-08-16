#!/usr/bin/python
# pylint: disable=C0103,C0111,W0612,R0904

import os
import wx
#import vtk
import sys
#from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

#import pyNastran
import pyNastran.gui
from pyNastran.gui.guiPanel import Pan
ID_SAVEAS = 803
ID_ABOUT = 3

ID_PLOT = 200

ID_SURFACE = 901
ID_WIREFRAME = 902
ID_HIDDEN = 903

ID_CAMERA = 910

ID_BDF = 920
ID_OP2 = 921
ID_F06 = 922
ID_CART3D = 923
ID_LAWGS = 924
ID_PANAIR = 925

ID_EXPORT = 930


pkgPath = pyNastran.gui.__path__[0]
#print "pkgPath = |%r|" %(pkgPath)

if '?' in pkgPath:
    iconPath = 'icons'
else:
    iconPath = os.path.join(pkgPath, 'icons')
#print "iconPath = |%r|" %(iconPath)

#------------------------------------------------------------------------------


class AppFrame(wx.Frame):

    def __init__(self, isEdges=False, isNodal=False, isCentroidal=False,
                 debug=False):

        wx.Frame.__init__(self, None, -1, size=wx.Size(800, 600),
                          title='pyNastran')
        self.bdfFileName = None
        self.isEdges = isEdges
        self.isNodal = isNodal
        self.isCentroidal = isCentroidal
        self.dirname = ''
        self.setupFrame()

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
        self.eventsHandler = EventsHandler(self, isNodal=self.isNodal,
                                           isCentroidal=self.isCentroidal)

        self.frmPanel = Pan(self, isEdges=self.isEdges, isNodal=self.isNodal,
                            isCentroidal=self.isCentroidal, size=(100, 200))

        self.buildMenuBar()
        self.buildToolBar()
        self.buildStatusBar()

        self.SetMenuBar(self.menubar)

        self.frmPanel.bdfFileName = self.bdfFileName
        self.frmPanel.buildVTK(self.bdfFileName)

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
        self.Bind(wx.EVT_MENU, self.frmPanel.widget.onTakePicture,
                  id=ID_CAMERA)
        self.Bind(wx.EVT_MENU, self.frmPanel.onSetToWireframe,
                  id=ID_WIREFRAME)
        self.Bind(wx.EVT_MENU, self.frmPanel.onSetToSurface,
                  id=ID_SURFACE)

        #self.Bind(wx.EVT_MENU, self.frmPanel.onSetToFlatShading,    self.flatShading)
        #self.Bind(wx.EVT_MENU, self.frmPanel.onSetToGouraudShading, self.gouraudShading)
        #self.Bind(wx.EVT_MENU, self.frmPanel.onSetToPhongShading,   self.phongShading)

        self.Bind(wx.EVT_MENU, events.onBackgroundColor, self.bkgColorView)
        #self.Bind(wx.EVT_MENU, events.onToggleStatusBar, self.showStatusBar)
        #self.Bind(wx.EVT_MENU, events.onToggleToolBar, self.showToolBar)

        # Bind Help Menu
        self.Bind(wx.EVT_MENU, events.onAbout, id=ID_ABOUT)
    #end __init__

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

    def UpdateWindowName(self, bdfFileName):
        self.bdfFileName = bdfFileName
        self.frmPanel.bdfFileName = bdfFileName
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

        topen = toolbar1.AddLabelTool(ID_BDF, '', wx.Bitmap(os.path.join(
            iconPath, 'topen.png')), longHelp='Loads a BDF')
        wireframe = toolbar1.AddLabelTool(ID_WIREFRAME, '', wx.Bitmap(os.path.join(iconPath, 'twireframe.png')), longHelp='Set to Wireframe Model')
        surface = toolbar1.AddLabelTool(ID_SURFACE, '', wx.Bitmap(os.path.join(iconPath, 'tsolid.png')), longHelp='Set to Surface/Solid Model')
        camera = toolbar1.AddLabelTool(ID_CAMERA, '', wx.Bitmap(os.path.join(
            iconPath, 'tcamera.png')), longHelp='Take a Screenshot')
        etool = toolbar1.AddLabelTool(wx.ID_EXIT, '', wx.Bitmap(os.path.join(
            iconPath, 'texit.png')), longHelp='Exit pyNastran GUI')
        #toolbar1.EnableTool(wx.ID_REDO, False)
        toolbar1.Realize()

        self.toolbar1 = toolbar1

        # Bind File Menu
        self.Bind(wx.EVT_TOOL, events.onLoadBDF, id=ID_BDF)
        self.Bind(wx.EVT_TOOL, events.onLoadOP2, id=ID_OP2)
        self.Bind(wx.EVT_TOOL, events.onLoadCart3d, id=ID_CART3D)
        self.Bind(wx.EVT_TOOL, events.onLoadLaWGS, id=ID_LAWGS)
        self.Bind(wx.EVT_TOOL, events.onLoadPanair, id=ID_PANAIR)
        #self.Bind(wx.EVT_TOOL, events.onExport,     id=ID_EXPORT)

        self.Bind(wx.EVT_MENU, events.onExit, id=wx.ID_EXIT)
       #self.Bind(wx.EVT_TOOL, events.onExit,     id=wx.ID_EXIT)

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
        fileMenu = wx.Menu()
        #fileMenu.Append(wx.ID_NEW,  '&New','does nothing')
        loadBDF = fileMenu.Append(ID_BDF,
                                  'Load &BDF',
                                  'Loads a BDF Input File')
        loadOP2 = fileMenu.Append(ID_OP2,
                                  'Load O&P2',
                                  'Loads an OP2 Results File')
        loadCart3d = fileMenu.Append(ID_CART3D,
                                     'Load &Cart3D',
                                     'Loads a Cart3D Input/Results File')
        loadLaWGS = fileMenu.Append(ID_LAWGS,
                                    'Load &LaWGS',
                                    'Loads an LaWGS File')
        loadPanair = fileMenu.Append(ID_PANAIR,
                                     'Load &Panair',
                                     'Loads a Panair Input File')
        #export     = fileMenu.Append(ID_EXPORT,'Export to...', 'Export the Model to Another Format')
        #print "topen = ",os.path.join(iconPath,'topen.png')
        sys.stdout.flush()
        assert os.path.exists(os.path.join(iconPath, 'topen.png'))
        loadBDF.SetBitmap(wx.Image(os.path.join(iconPath,
                                                'topen.png'), wx.BITMAP_TYPE_PNG).ConvertToBitmap())

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
        camera = viewMenu.Append(ID_CAMERA,
                                 'Take a Screenshot',
                                 'Take a Screenshot')
        viewMenu.AppendSeparator()
        wireframe = viewMenu.Append(ID_WIREFRAME,
                                    'Wireframe Model',
                                    'Show Model as a Wireframe Model')
        surface = viewMenu.Append(ID_SURFACE,
                                  'Surface Model',
                                  'Show Model as a Surface Model')
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

    def __init__(self, parent, isNodal=False, isCentroidal=True):
        self.parent = parent
        self.isNodal = isNodal
        self.isCentroidal = isCentroidal

    # File Menu
    def onLoadBDF(self, event):
        """ Open a file"""
        #print "OnOpen..."

        wildcard = "Nastran BDF (*.bdf; *.dat; *.nas)|*.bdf;*.dat;*.nas|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Nastran Input Deck to Load'
        loadFunction = self.parent.frmPanel.loadNastranGeometry
        self.createLoadFileDialog(wildcard, Title, loadFunction,
                                  updateWindowName=True)

    def onLoadOP2(self, event):
        """ Open a file"""
        #print "OnOpen..."

        if 0:
            bdf = self.parent.bdfFileName
            bdfBase = os.path.basename(bdf)
            #dirname = os.path.dirname(bdf)
            (op2name, op2) = os.path.splitext(bdfBase)
            op2 = os.path.join(self.parent.dirname, op2name + '.op2')

            self.parent.op2FileName = op2
            if os.path.exists(op2):
                self.parent.frmPanel.loadNastranResults(op2)
                self.parent.frmPanel.Update()
            return

        wildcard = "Nastran OP2 (*.op2)|*.op2|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Nastran Output File to Load (OP2 only)'
        loadFunction = self.parent.frmPanel.loadNastranResults
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
        loadFunction = self.parent.frmPanel.loadCart3dGeometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        self.createLoadFileDialog(wildcard, Title, loadFunction,
                                  updateWindowName=True)

    def onLoadLaWGS(self, event):
        """ Open a file"""
        wildcard = "LaWGS (*.wgs)|*.wgs|" \
            "All files (*.*)|*.*"

        Title = 'Choose a LaWGS File to Load'
        loadFunction = self.parent.frmPanel.loadLaWGSGeometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        self.createLoadFileDialog(wildcard, Title, loadFunction,
                                  updateWindowName=True)

    def onLoadPanair(self, event):
        """ Open a file"""
        wildcard = "Panair (*.inp)|*.inp|" \
            "All files (*.*)|*.*"

        Title = 'Choose a Panair Input File to Load'
        loadFunction = self.parent.frmPanel.loadPanairGeometry
        #fname = r'C:\Users\steve\Desktop\pyNastran\pyNastran\converters\cart3d\Cart3d_35000_0.825_10_0_0_0_0.i.triq'
        #dirname = ''
        #loadFunction(fname,dirname)
        self.createLoadFileDialog(wildcard, Title, loadFunction,
                                  updateWindowName=True)

    def createLoadFileDialog(self, wildcard, Title, loadFunction,
                             updateWindowName=False):
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
            loadFunction(fname, self.parent.dirname,
                         self.isNodal, self.isCentroidal)
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
        loadFunction = self.parent.frmPanel.loadNastranResults
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


def runArgParse():
    import argparse

    ver = str(pyNastran.__version__)
    parser = argparse.ArgumentParser(description='Tests to see if an OP2 will work with pyNastran.', add_help=True)  # version=ver)
    #parser.add_argument('op2FileName', metavar='op2FileName', type=str, nargs=1,
    #                   help='path to OP2 file')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', dest='quiet', action='store_true',
                       help='prints debug messages (default=True)')

    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument('-e', '--edges', dest='edges', action='store_true',
                        help='shows element edges as black lines')
    #group2.add_argument('-w','--writeBDF', dest='writeBDF', action='store_true', help='Writes the bdf to fem.bdf.out')

    group3 = parser.add_mutually_exclusive_group()
    group3.add_argument('-n', '--nodalResults', dest='isNodal',
                        action='store_true', help='plots nodal results')
    group3.add_argument('-c', '--centroidalResults', dest='isNodal',
                        action='store_false', help='plots centroidal results')

    parser.add_argument('-v', '--version', action='version', version=ver)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    #print "op2FileName = ", args.op2FileName[0]
    #print "debug       = ", not(args.quiet)

    debug = not(args.quiet)
    edges = args.edges
    isNodal = args.isNodal
    #writeBDF    = args.writeBDF
    #op2FileName = args.op2FileName[0]

    return (edges, isNodal, debug)


def main():
    isEdges = False
    isNodal = True
    #isCentroidal = True
    debug = True
    if sys.version_info < (2, 6):
        print("requires Python 2.6+ to use command line arguments...")
    else:
        if len(sys.argv) > 1:
            (isEdges, isNodal, debug) = runArgParse()
    isCentroidal = not(isNodal)

    isNodal = True
    isCentroidal = True

    app = wx.App(redirect=False)
    appFrm = AppFrame(isEdges, isNodal, isCentroidal, debug)
    #appFrm.Show()
    print("launching gui")
    app.MainLoop()

#end class

#==============================================================================

if __name__ == '__main__':
    main()
