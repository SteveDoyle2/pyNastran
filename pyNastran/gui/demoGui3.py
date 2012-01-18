#!/usr/bin/python

import os
import wx
import vtk
#from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

from guiPanel import Pan

#ID_OPEN = 801
ID_SAVEAS = 803
ID_ABOUT = 3

ID_SURFACE   = 901
ID_WIREFRAME = 902
ID_HIDDEN    = 903

ID_CAMERA    = 910

#------------------------------------------------------------------------------

class AppFrame( wx.Frame ) :

    def __init__( self) :
        
        wx.Frame.__init__( self, None, -1, title='pyNastran' )
        self.bdfFileName = None
        self.dirname = ''
        self.setupFrame()

    def setupFrame(self):

        # Must call before any event handler is referenced.
        self.eventsHandler = EventsHandler(self)

        self.frmPanel = Pan(self)

        self.buildMenuBar()
        self.buildToolBar()
        #self.buildToolBar2()
        self.buildStatusBar()
        
        self.SetMenuBar(self.menubar)

        self.frmPanel.bdfFileName = self.bdfFileName
        self.frmPanel.buildVTK(self.bdfFileName)

        windowName = self.frmPanel.getWindowName()
        self.SetTitle(windowName)

        # Add them to sizer.
        hbox = wx.BoxSizer( wx.HORIZONTAL )
        hbox.Add( self.frmPanel.widget, 1, wx.EXPAND|wx.ALL, 1 )

        # Add buttons in their own sizer
        if 1:
            self.redBtn   = wx.Button( self.frmPanel, label='Red' )
            self.greenBtn = wx.Button( self.frmPanel, label='Green' )
            self.exitBtn  = wx.Button( self.frmPanel, label='Exit' )

            buttonSizer = wx.BoxSizer( wx.VERTICAL )
            buttonSizer.AddStretchSpacer()
            buttonSizer.Add( self.redBtn,   proportion=0, flag=wx.EXPAND|wx.ALL, border=5 )
            buttonSizer.Add( self.greenBtn, proportion=0, flag=wx.EXPAND|wx.ALL, border=5 )
            buttonSizer.Add( self.exitBtn,  proportion=0, flag=wx.EXPAND|wx.ALL, border=5 )
            buttonSizer.AddStretchSpacer()

            hbox.Add( buttonSizer, 0, wx.EXPAND| wx.ALL, 5 )

        # SetSizer both sizers in the most senior control that has sizers in it.
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        #self.vbox.AddStretchSpacer()
        #self.vbox.Add(self.frmPanel.widget, 0, wx.EXPAND)
        #self.vbox.Add(self.toolbar1, 0, wx.EXPAND)
        self.vbox.AddStretchSpacer()
        #self.vbox.Add(hbox, 0, wx.EXPAND|wx.ALL, 5)
        self.vbox.Add(hbox)
        #self.frmPanel.SetSizer(self.vbox)
        self.frmPanel.SetSizer(hbox)
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


    def UpdateWindowName(self,bdfFileName):
        self.bdfFileName = bdfFileName
        self.frmPanel.bdfFileName = bdfFileName
        windowName = self.frmPanel.getWindowName()
        self.SetTitle(windowName)
    
    def buildStatusBar(self):
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetStatusText('Ready')

    def buildToolBar2(self):
        self.toolbar1 = wx.ToolBar(self)
        topen = self.toolbar1.AddLabelTool(wx.ID_OPEN, '', wx.Bitmap('icons/topen.png'))
        qtool = self.toolbar1.AddLabelTool(wx.ID_EXIT, '', wx.Bitmap('icons/texit.png'))
        #self.vbox.Add(self.toolbar1)

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
        topen = toolbar1.AddLabelTool(wx.ID_OPEN, '', wx.Bitmap('icons/topen.png'))
        wireframe = toolbar1.AddLabelTool(ID_WIREFRAME, 'Set to Wireframe Model', wx.Bitmap('icons/twireframe.png'))
        surface   = toolbar1.AddLabelTool(ID_SURFACE,   'Set to Surface Model',   wx.Bitmap('icons/tsolid.png'))
        camera    = toolbar1.AddLabelTool(ID_CAMERA,    'Take a Screenshot',      wx.Bitmap('icons/tcamera.png'))
        etool = toolbar1.AddLabelTool(wx.ID_EXIT,       '',                       wx.Bitmap('icons/texit.png'))
        #toolbar1.EnableTool(wx.ID_REDO, False)
        toolbar1.Realize()

        self.toolbar1 = toolbar1


        # Bind File Menu
        self.Bind(wx.EVT_TOOL, events.OnLoadBDF,  id=wx.ID_OPEN)

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
        loadBDF = fileMenu.Append(wx.ID_OPEN, '&Load BDF','Loads a BDF')
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

        dlg = wx.FileDialog(None, "Choose a file", self.parent.dirname, "", wildcard, wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            bdfFileName         = dlg.GetFilename()
            self.parent.dirname = dlg.GetDirectory()
            fname = os.path.join(self.parent.dirname, bdfFileName)
            print "fname = ",fname
            self.parent.UpdateWindowName(bdfFileName)
            self.parent.frmPanel.loadGeometry(bdfFileName)
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
            'pyNastran v0.3.0',
            'Copyright Steven P. Doyle 2011-2012\n',
            'code.google.com/p/pynastran/',
            '',
            'Controls',
              'X/x - snap to x axis',
              'Y/y - snap to y axis',
              'Z/z - snap to z axis',
              '',
              'm/M scale up/scale down',
              'o/O rotate counter-clockwise/clockwise 5 degrees',
              'w   wireframe',
              's   surface',
              'i   take a screenshot (image)',]
        
        # not done
              #'',
              #'left arrow  - pan left  (not done)',
              #'right arrow - pan right (not done)',
              #'up arrow    - pan up    (not done)',
              #'down arrow  - pan down  (not done)',
              #'',
              #'p   project point (not done)',
              #'f   fly to rotation point (not done)',

        dlg = wx.MessageDialog(None, '\n'.join(about), 'About',
                 wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

#end Events class

#------------------------------------------------------------------------------

def Main():
    app = wx.App( redirect=False )
    appFrm = AppFrame()
    appFrm.Show()
    app.MainLoop()

#end class

#==============================================================================

if __name__ == '__main__' :

    Main()