import os
import wx


#ID_OPEN = 801
ID_SAVEAS = 803
ID_ABOUT = 3


class MyPopupMenu(wx.Menu):
    def __init__(self, parent):
        super(MyPopupMenu, self).__init__()

        self.parent = parent

        mmi = wx.MenuItem(self, wx.NewId(), 'Minimize')
        self.AppendItem(mmi)
        self.Bind(wx.EVT_MENU, self.OnMinimize, mmi)

        cmi = wx.MenuItem(self, wx.NewId(), 'Close')
        self.AppendItem(cmi)
        self.Bind(wx.EVT_MENU, self.OnClose, cmi)

    def OnMinimize(self, e):
        self.parent.Iconize()

    def OnClose(self, e):
        self.parent.Close()


class Example(wx.Frame):
    def __init__(self, *args, **kwargs):
        super(Example, self).__init__(*args, **kwargs)
        self.InitUI()

    def InitUI(self):
        # max undo count
        self.count = 5

        menubar = wx.MenuBar()
        vbox = wx.BoxSizer(wx.VERTICAL)

        # file menu
        fileMenu = wx.Menu()
        fileMenu.Append(wx.ID_NEW, '&New', 'does nothing')
        fileMenu.Append(wx.ID_OPEN, '&Load BDF', 'Loads a BDF')
        fileMenu.Append(wx.ID_OPEN, 'Load OP2 &Results',
                        'Loads a OP2 - does nothing')
        fileMenu.Append(wx.ID_SAVE, '&Save', 'does nothing')
        exitButton = fileMenu.Append(
            wx.ID_EXIT, 'E&xit pyNastran', 'Exits the program')
        exitButton.SetBitmap(wx.Image(
            'icons/texit.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        #print "exitButton = ",'\n'.join(dir(exitButton))
        #print "exitButton = ",help(exitButton.SetBitmap)

        fileMenu.AppendSeparator()

        imp = wx.Menu()
        imp.Append(wx.ID_ANY, 'Import newsfeed list...')
        imp.Append(wx.ID_ANY, 'Import bookmarks...')
        imp.Append(wx.ID_ANY, 'Import mail...')

        fileMenu.AppendMenu(wx.ID_ANY, 'I&mport', imp)

        qmi = wx.MenuItem(fileMenu, wx.ID_EXIT, '&Quit\tCtrl+W')
        fileMenu.AppendItem(qmi)

        self.Bind(wx.EVT_MENU, self.OnQuit, qmi)

        # view menu
        viewMenu = wx.Menu()

        self.shst = viewMenu.Append(wx.ID_ANY, 'Show statusbar',
                                    'Show Statusbar', kind=wx.ITEM_CHECK)
        self.shtl = viewMenu.Append(wx.ID_ANY, 'Show toolbar',
                                    'Show Toolbar', kind=wx.ITEM_CHECK)

        viewMenu.Check(self.shst.GetId(), True)
        viewMenu.Check(self.shtl.GetId(), True)

        self.Bind(wx.EVT_MENU, self.ToggleStatusBar, self.shst)
        self.Bind(wx.EVT_MENU, self.ToggleToolBar, self.shtl)
        self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightDown)

        # help/about menu
        helpMenu = wx.Menu()
        helpMenu.Append(ID_ABOUT, '&About', 'About pyNastran')
        self.Bind(wx.EVT_MENU, self.OnAbout, id=ID_ABOUT)

        # menu bar
        menubar.Append(fileMenu, '&File')
        menubar.Append(viewMenu, '&View')
        menubar.Append(helpMenu, '&Help')
        self.SetMenuBar(menubar)

        # toolbar at top - toggles
        self.toolbar1 = wx.ToolBar(self)
        #self.toolbar1.AddLabelTool(wx.ID_ANY, '', wx.Bitmap('icons/new.png'))
        self.toolbar1.AddLabelTool(
            wx.ID_OPEN, '', wx.Bitmap('icons/topen.png'))
        self.toolbar1.AddLabelTool(
            ID_SAVEAS, '', wx.Bitmap('icons/tsave.png'))
        self.Bind(wx.EVT_TOOL, self.OnSaveAsFile, id=ID_SAVEAS)
        self.Bind(wx.EVT_TOOL, self.OnOpenBDF, id=wx.ID_OPEN)
        #self.toolbar1.AddSeparator()
        tundo = self.toolbar1.AddLabelTool(
            wx.ID_UNDO, '', wx.Bitmap('icons/tundo.png'))
        #self.toolbar1.AddSeparator()
        tredo = self.toolbar1.AddLabelTool(
            wx.ID_REDO, '', wx.Bitmap('icons/tredo.png'))
        self.toolbar1.EnableTool(wx.ID_REDO, False)

        self.toolbar1.Realize()

        # toolbar 2
        toolbar2 = wx.ToolBar(self)
        qtool = toolbar2.AddLabelTool(
            wx.ID_EXIT, '', wx.Bitmap('icons/texit.png'))
        toolbar2.Realize()

        vbox.Add(self.toolbar1, 0, wx.EXPAND)
        vbox.Add(toolbar2, 0, wx.EXPAND)
        self.Bind(wx.EVT_TOOL, self.OnQuit, qtool)
        self.Bind(wx.EVT_TOOL, self.OnUndo, tundo)
        self.Bind(wx.EVT_TOOL, self.OnRedo, tredo)
        self.SetSizer(vbox)

        # status bar at bottom - toggles
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetStatusText('Ready')

        self.SetSize((350, 250))
        #self.SetIcon(wx.Icon('icons/tbat.png', wx.BITMAP_TYPE_ICO))
        self.SetTitle('pyNastran')
        self.Centre()
        self.Show(True)

    def OnRightDown(self, e):
        self.PopupMenu(MyPopupMenu(self), e.GetPosition())

    def ToggleStatusBar(self, e):
        if self.shst.IsChecked():
            self.statusbar.Show()
        else:
            self.statusbar.Hide()

    def ToggleToolBar(self, e):
        if self.shtl.IsChecked():
            self.toolbar1.Show()
        else:
            self.toolbar1.Hide()

    def OnUndo(self, e):
        if self.count > 1 and self.count <= 5:
            self.count = self.count - 1

        if self.count == 1:
            self.toolbar1.EnableTool(wx.ID_UNDO, False)

        if self.count == 4:
            self.toolbar1.EnableTool(wx.ID_REDO, True)

    def OnRedo(self, e):
        if self.count < 5 and self.count >= 1:
            self.count = self.count + 1

        if self.count == 5:
            self.toolbar1.EnableTool(wx.ID_REDO, False)

        if self.count == 2:
            self.toolbar1.EnableTool(wx.ID_UNDO, True)

    def OnAbout(self, event):
        about = [
            'pyNastran v0.3.0',
            'Copyright 2011-2012\n',
            'code.google.com/p/pynastran/',
            '',
            'Controls',
            'X/x - snap to x axis',
            'Y/y - snap to axis',
            'Z/z - snap to axis',
            '',
            'left arrow  - pan left',
            'right arrow - pan right',
            'up arrow    - pan up',
            'down arrow   - pan down',
            '',
            'm/M /tscale up/scale down',
            'p   /tproject point (not done)',
            'f   /tfly to rotation point (not done)',
            'q/e /texit (to disable)',
            'o/O /trotate counter-clockwise/clockwise 5 degrees',
            'w   /twireframe',
            's   /tsurface',
            'i   /ttake a screenshot (image, not done)', ]

        dlg = wx.MessageDialog(self, '\n'.join(about), 'About',
                               wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def OnLoadBDF(self, event):
        """ Open a file"""
        #print "OnOpen..."
        self.dirname = ''
        dlg = wx.FileDialog(
            self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            fname = os.path.join(self.dirname, self.filename)
            print "fname = ", fname
            #f = open(fname, 'r')
            #self.control.SetValue(f.read())
            #print f
            #f.close()
        dlg.Destroy()

    def OnSaveAsFile(self, event):
        wcd = 'All files(*)|*|pyNastran Database (*.pndb)|*.pndb|'
        dir = os.getcwd()
        save_dlg = wx.FileDialog(self, message='Save file as...', defaultDir=dir, defaultFile='',
                                 wildcard=wcd, style=wx.SAVE | wx.OVERWRITE_PROMPT)
        if save_dlg.ShowModal() == wx.ID_OK:
            path = save_dlg.GetPath()
            try:
                print "save path = ", path
            except IOError, error:
                dlg = wx.MessageDialog(
                    self, 'Error saving file\n' + str(error))
                #dlg.ShowModal()
        save_dlg.Destroy()

    def OnQuit(self, e):
        self.Close()


def main():

    ex = wx.App(0)
    Example(None)
    #Example(title='pyNastran GUI')
    ex.MainLoop()


if __name__ == '__main__':
    main()
