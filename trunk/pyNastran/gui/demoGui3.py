#!/usr/bin/python

import os
import wx
import vtk
#from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

import pyNastran
from guiPanel import Pan
ID_SAVEAS = 803
ID_ABOUT = 3

ID_SURFACE   = 901
ID_WIREFRAME = 902
ID_HIDDEN    = 903

ID_CAMERA    = 910

ID_BDF = 920
ID_OP2 = 921



#------------------------------------------------------------------------------

class AppFrame( wx.Frame ) :

    def __init__( self) :
        
        wx.Frame.__init__( self, None, -1, size=wx.Size(800, 600),title='pyNastran' )
        self.bdfFileName = None
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
        self.eventsHandler = EventsHandler(self)

        self.frmPanel = Pan(self,size=(100,200))

        self.buildMenuBar()
        self.buildToolBar()
        self.buildStatusBar()
        
        self.SetMenuBar(self.menubar)

        self.frmPanel.bdfFileName = self.bdfFileName
        self.frmPanel.buildVTK(self.bdfFileName)

        windowName = self.frmPanel.getWindowName()
        self.SetTitle(windowName)
        #self.SetSize([600,600])
        self.Centre()

        # Add them to sizer.
        hbox = wx.BoxSizer( wx.HORIZONTAL )
        hbox.Add( self.frmPanel.widget, 1, wx.EXPAND|wx.ALL, 1 )

        # Add buttons in their own sizer
        if 0:
            self.redBtn   = wx.Button( self.frmPanel, label='Red' )
            self.greenBtn = wx.Button( self.frmPanel, label='Green' )
            self.exitBtn  = wx.Button( self.frmPanel, label='Exit' )

            vRight = wx.BoxSizer( wx.VERTICAL)
            vRight.AddStretchSpacer()
            vRight.Add( self.greenBtn, proportion=0, flag=wx.EXPAND|wx.ALL, border=5 )
            vRight.Add( self.exitBtn,  proportion=0, flag=wx.EXPAND|wx.ALL, border=5 )
            vRight.Add( self.redBtn,   proportion=0, flag=wx.EXPAND|wx.ALL, border=5 )
            vRight.AddStretchSpacer()
            hbox.Add(vRight,1,wx.EXPAND)

        if 0:
            tree = self.buildTree(self.frmPanel)
            vRight = wx.BoxSizer( wx.VERTICAL)
            vRight.AddStretchSpacer()
            vRight.Add( tree, proportion=1, flag=wx.EXPAND|wx.ALL, border=5 )
            vRight.AddStretchSpacer()
            hbox.Add(vRight,1,wx.EXPAND)


        # best guess at tree
        if 0:
            #self.tree  = wx.Button( self.frmPanel, label='Tree' )

            vRight = wx.BoxSizer( wx.VERTICAL)

            scroll = wx.ScrolledWindow( self, -1 )
            panelRight = wx.Panel( scroll, -1 )

            panelRight = wx.Panel(self, wx.EXPAND)
            tree = self.buildTree(panelRight)
            vRight.Add( tree,   flag=wx.EXPAND|wx.ALL)
            #vRight.Add(tree, 1, wx.EXPAND)
            #hbox.Add(panel1, 1, wx.EXPAND)
            
            vRight.Add(scroll, 1, wx.EXPAND | wx.ALL)
            #panelRight.SetSizer(vRight)
            panelRight.Layout()

            

            hbox.Add( vRight, 1, wx.EXPAND)
            #hbox.Add(panelRight, 1, wx.EXPAND | wx.ALL)

            # SetSizer both sizers in the most senior control that has sizers in it.
            self.vMain = wx.BoxSizer(wx.VERTICAL | wx.EXPAND)
            self.vMain.Add(hbox,1,wx.EXPAND,5)

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
        self.Bind(wx.EVT_MENU, self.frmPanel.widget.TakePicture, id=ID_CAMERA)
        self.Bind(wx.EVT_MENU, self.frmPanel.SetToWireframe,     id=ID_WIREFRAME)
        self.Bind(wx.EVT_MENU, self.frmPanel.SetToSurface,       id=ID_SURFACE)

        #self.Bind(wx.EVT_MENU, self.frmPanel.SetToFlatShading,    self.flatShading)
        #self.Bind(wx.EVT_MENU, self.frmPanel.SetToGouraudShading, self.gouraudShading)
        #self.Bind(wx.EVT_MENU, self.frmPanel.SetToPhongShading,   self.phongShading)

        self.Bind(wx.EVT_MENU, events.OnBackgroundColor, self.bkgColorView)
        #self.Bind(wx.EVT_MENU, events.ToggleStatusBar, self.showStatusBar)
        #self.Bind(wx.EVT_MENU, events.ToggleToolBar, self.showToolBar)
        
        # Bind Help Menu
        self.Bind(wx.EVT_MENU, events.OnAbout, id=ID_ABOUT)
    #end __init__

    def buildTree(self,panel1):
        tree = wx.TreeCtrl(self, 1, wx.DefaultPosition, (-1,-1), wx.TR_HIDE_ROOT|wx.TR_HAS_BUTTONS)
        root = tree.AddRoot('Programmer')
        os   = tree.AppendItem(root, 'Operating Systems')
        tree.AppendItem(os, 'Linux')
        tree.AppendItem(os, 'FreeBSD')
        tree.AppendItem(os, 'OpenBSD')
        tree.AppendItem(os, 'NetBSD')
        tree.AppendItem(os, 'Solaris')
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
        item =  event.GetItem()
        self.display.SetLabel(tree.GetItemText(item))

    def UpdateWindowName(self,bdfFileName):
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
        #tnew  = toolbar1.AddLabelTool(wx.ID_ANY,  '', wx.Bitmap('icons/new.png'))
        #tsave = toolbar1.AddLabelTool(ID_SAVEAS,  '', wx.Bitmap('icons/tsave.png'))
        #tundo = toolbar1.AddLabelTool(wx.ID_UNDO, '', wx.Bitmap('icons/tundo.png'))
        #tredo = toolbar1.AddLabelTool(wx.ID_REDO, '', wx.Bitmap('icons/tredo.png'))

        # toolbar at top - toggles
        toolbar1 = self.CreateToolBar()
        #toolbar.AddLabelTool(self.id, '', bitmap, wx.NullBitmap, self.kind, 
        #                     shortHelp=wx.MenuItem.GetLabelFromText(self.menuText),
        #             longHelp=self.helpText)

        topen     = toolbar1.AddLabelTool(ID_BDF,       '' ,wx.Bitmap('icons/topen.png'),     longHelp='Loads a BDF')
        wireframe = toolbar1.AddLabelTool(ID_WIREFRAME, '', wx.Bitmap('icons/twireframe.png'),longHelp='Set to Wireframe Model')
        surface   = toolbar1.AddLabelTool(ID_SURFACE,   '', wx.Bitmap('icons/tsolid.png'),    longHelp='Set to Surface/Solid Model')
        camera    = toolbar1.AddLabelTool(ID_CAMERA,    '', wx.Bitmap('icons/tcamera.png'),   longHelp='Take a Screenshot')
        etool     = toolbar1.AddLabelTool(wx.ID_EXIT,   '', wx.Bitmap('icons/texit.png'),     longHelp='Set to Surface Model')
        #toolbar1.EnableTool(wx.ID_REDO, False)
        toolbar1.Realize()

        self.toolbar1 = toolbar1


        # Bind File Menu
        self.Bind(wx.EVT_TOOL, events.OnLoadBDF,  id=ID_BDF)
        self.Bind(wx.EVT_TOOL, events.OnLoadOP2,  id=ID_OP2)

        self.Bind(wx.EVT_MENU, events.OnExit,     id=wx.ID_EXIT)
        #self.Bind(wx.EVT_TOOL, events.OnExit,     id=wx.ID_EXIT)

        self.Bind(wx.EVT_MENU, self.frmPanel.SetToWireframe, id=ID_WIREFRAME)
        self.Bind(wx.EVT_MENU, self.frmPanel.SetToSurface,   id=ID_SURFACE)
        self.Bind(wx.EVT_MENU, self.frmPanel.SetToSurface,   id=ID_CAMERA)


        #self.Bind(wx.EVT_TOOL, events.OnSaveAsFile, id=ID_SAVEAS)
        #self.Bind(wx.EVT_TOOL, events.OnUndo, tundo)
        #self.Bind(wx.EVT_TOOL, events.OnRedo, tredo)

    def buildMenuBar(self):
        events = self.eventsHandler

        menubar = wx.MenuBar()
        # file menu
        fileMenu = wx.Menu()
        #fileMenu.Append(wx.ID_NEW,  '&New','does nothing')
        loadBDF = fileMenu.Append(ID_BDF, 'Load &BDF', 'Loads a BDF Input File')
        loadOP2 = fileMenu.Append(ID_OP2, 'Load O&P2', 'Loads an OP2 Results File')
        loadBDF.SetBitmap(   wx.Image('icons/topen.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        #fileMenu.Append(wx.ID_RES, 'Load OP2 &Results','Loads a OP2 - does nothing')
        #fileMenu.Append(wx.ID_SAVE, '&Save','does nothing')
        fileMenu.AppendSeparator()

        
        # dummy import submenu
        #imp = wx.Menu()
        #imp.Append(wx.ID_ANY, 'Import newsfeed list...')
        #imp.Append(wx.ID_ANY, 'Import bookmarks...')
        #imp.Append(wx.ID_ANY, 'Import mail...')

        #fileMenu.AppendMenu(wx.ID_ANY, 'I&mport', imp)
        exitButton = wx.MenuItem(fileMenu, wx.ID_EXIT, 'Exit','Exits pyNastran')
        exitButton.SetBitmap(wx.Image('icons/texit.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        fileMenu.AppendItem(exitButton)

        # view menu
        # status bar at bottom - toggles
        viewMenu = wx.Menu()
        camera    = viewMenu.Append(ID_CAMERA, 'Take a Screenshot','Take a Screenshot')
        viewMenu.AppendSeparator()
        wireframe = viewMenu.Append(ID_WIREFRAME, 'Wireframe Model','Show Model as a Wireframe Model')
        surface   = viewMenu.Append(ID_SURFACE, 'Surface Model',  'Show Model as a Surface Model')
        #viewMenu.AppendSeparator()

        #self.flatShading    = viewMenu.Append(wx.ID_ANY, 'Flat Shading',           'Flat Shading')
        #self.gouraudShading = viewMenu.Append(wx.ID_ANY, 'Mid (Gouraud) Shading',  'Mid (Gouraud) Shading')
        #self.phongShading   = viewMenu.Append(wx.ID_ANY, 'Smooth (Phong) Shading', 'Smooth (Phong) Shading')

        # view images
        wireframe.SetBitmap(wx.Image('icons/twireframe.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        surface.SetBitmap(  wx.Image('icons/tsolid.png',     wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        camera.SetBitmap(   wx.Image('icons/tcamera.png',    wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        #wireframe = toolbar1.AddLabelTool(ID_WIREFRAME, 'Set to Wireframe Model', wx.Bitmap('icons/twireframe.png'))
        #surface   = toolbar1.AddLabelTool(ID_SURFACE,   'Set to Surface Model',   wx.Bitmap('icons/tsolid.png'))
        #camera    = toolbar1.AddLabelTool(ID_CAMERA,    'Take a Screenshot',      wx.Bitmap('icons/tcamera.png'))


        viewMenu.AppendSeparator()
        self.bkgColorView  = viewMenu.Append(wx.ID_ANY, 'Change Background Color','Change Background Color')
        #self.showStatusBar = viewMenu.Append(wx.ID_ANY, 'Show statusbar', 'Show Statusbar', kind=wx.ITEM_CHECK)
        #self.showToolBar   = viewMenu.Append(wx.ID_ANY, 'Show toolbar',   'Show Toolbar',   kind=wx.ITEM_CHECK)
        #viewMenu.Check(self.showStatusBar.GetId(), True)
        #viewMenu.Check(self.showToolBar.GetId(),   True)

        # help/about menu
        helpMenu = wx.Menu()
        helpMenu.Append(ID_ABOUT, '&About', 'About pyNastran')

        # menu bar
        menubar.Append(fileMenu, '&File')
        menubar.Append(viewMenu, '&View')
        menubar.Append(helpMenu, '&Help')
        self.menubar = menubar



#end AppFrame class

#------------------------------------------------------------------------------

class EventsHandler(object) :

    def __init__(self,parent):
        self.parent = parent

    # File Menu
    def OnLoadBDF(self, event):
        """ Open a file"""
        #print "OnOpen..."

        if 0:
            if self.parent.bdfFileName=='test_tet10.bdf':
                bdfFileName = 'box_cylindrical_coord_sys.bdf'
            else:
                bdfFileName = 'test_tet10.bdf'
            self.parent.bdfFileName = bdfFileName
            #self.parent.frmPanel.buildVTK(bdfFileName)
            self.parent.frmPanel.loadGeometry(bdfFileName)
            #self.parent.frmPanel.getWindow.Render()
            self.parent.frmPanel.Update()
            return
        
        wildcard = "Nastran BDF (*.bdf; *.dat; *.nas)|*.bdf;*.dat;*.nas|" \
         "All files (*.*)|*.*"

        dlg = wx.FileDialog(None, "Choose a Nastran Input Deck to Load ", self.parent.dirname, "", wildcard, wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            bdfFileName         = dlg.GetFilename()
            self.parent.dirname = dlg.GetDirectory()
            fname = os.path.join(self.parent.dirname, bdfFileName)
            print "fname = ",fname
            self.parent.UpdateWindowName(bdfFileName)
            self.parent.frmPanel.loadGeometry(fname,self.parent.dirname)
            self.parent.frmPanel.Update()
        dlg.Destroy()

    def OnLoadOP2(self, event):
        """ Open a file"""
        #print "OnOpen..."

        if 1:
            bdf = self.parent.bdfFileName
            bdfBase = os.path.basename(bdf)
            #dirname = os.path.dirname(bdf)
            op2name,op2 = os.path.splitext(bdfBase)
            op2 = os.path.join(self.parent.dirname,op2name+'.op2')
            
            self.parent.op2FileName = op2
            if os.path.exists(op2):
                self.parent.frmPanel.loadResults(op2)
                self.parent.frmPanel.Update()
            return
        
        wildcard = "Nastran OP2 (*.op2)|*.op2|" \
         "All files (*.*)|*.*"

        dlg = wx.FileDialog(None, "Choose a Nastran Output File to Load (OP2 only)", self.parent.dirname, "", wildcard, wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            op2FileName         = dlg.GetFilename()
            self.parent.dirname = dlg.GetDirectory()
            oname = os.path.join(self.parent.dirname, op2FileName)
            print "oname = ",oname
            #self.parent.UpdateWindowName(op2FileName)
            self.parent.frmPanel.loadResults(op2FileName)
            self.parent.frmPanel.Update()
        dlg.Destroy()

    def OnExit(self,event):
        self.parent.Destroy()

    # View Menu
    def OnBackgroundColor(self,event):
        dlg = wx.ColourDialog(self.parent)
        dlg.GetColourData().SetChooseFull(True)
        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetColourData()
            #print 'You selected: %s\n' % str(data.GetColour().Get())
            dlg.Destroy()
            self.SetColor(data.GetColour().Get())
            return
        dlg.Destroy()

    def SetColor(self,bkgColor):
        """
        @warning if the model is not loaded and the color is not "clicked"
        on, then rend will return None
        """
        rend = self.parent.frmPanel.widget.GetCurrentRenderer()
        if not rend: # rend doesnt exist at first
            print 'Try again to change the color (you didnt "click" on the color; one time bug)'
            return
            rend = self.parent.frmPanel.widget.GetCurrentRenderer()

        ## bkgColor range from 0 to 255
        ## color ranges from 0 to 1
        color = [bkgColor[0]/255.,bkgColor[1]/255.,bkgColor[2]/255.]
        rend.SetBackground(color)
        self.parent.frmPanel.widget.Render()

    def ToggleStatusBar(self, e):
        if self.parent.showStatusBar.IsChecked():
            self.parent.statusbar.Show()
        else:
            self.parent.statusbar.Hide()

    def ToggleToolBar(self, e):
        if self.parent.showToolBar.IsChecked():
            self.parent.toolbar1.Show()
        else:
            self.parent.toolbar1.Hide()

    # Help Menu
    def OnAbout(self, event):
        about = [
            'pyNastran v%s'%(pyNastran.__version__), 
            'Copyright Steven P. Doyle 2011-2012\n',
            'code.google.com/p/pynastran/',
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

def Main():
    app = wx.App( redirect=False )
    appFrm = AppFrame()
    #appFrm.Show()
    app.MainLoop()

#end class

#==============================================================================

if __name__ == '__main__' :

    Main()