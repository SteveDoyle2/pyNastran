# -*- coding: utf-8 -*-

###########################################################################
## Python code generated with wxFormBuilder (version Sep  8 2010)
## http://www.wxformbuilder.org/
##
## PLEASE DO "NOT" EDIT THIS FILE!
###########################################################################

import wx
is_wx = True

###########################################################################
## Class MyDialog3
###########################################################################

def get_string_width(name):
    #f = window.GetFont()
    #dc = wx.WindowDC(window)
    #dc.SetFont(f)
    #width, height = dc.GetTextExtent("Text to measure")
    return 100 + len(name) * 10

def bold(obj):
    f = obj.GetFont()
    f.SetWeight(wx.BOLD)
    obj.SetFont(f)

class GetChoices(wx.Dialog):

    def __init__(self, parent, args=None):
        self.args = args
        if args is None:
            args = {}

        names = args['names']
        n_names = len(names)

        widths = []
        for name in names:
            width = get_string_width(name)
            widths.append(width)
        width = max(widths)

        height = 200 + 20 * n_names
        max_height = 400
        if height > max_height:
            height = max_height

        wx.Dialog.__init__ ( self, parent,
                     id = wx.ID_ANY,
                     title = 'Select the Cases to Extract Results For',
                     pos = wx.DefaultPosition,
                     size = wx.Size(315,height),
                     style = wx.DEFAULT_DIALOG_STYLE )
        #self.build(height)


    #def build(self, height):
        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )

        bSizer1 = wx.BoxSizer( wx.VERTICAL )

        self.checkboxes = []
        #self.staticTexts= []

        #gSizer1 = wx.GridSizer(n_names + 1, 2, 0, 0 )
        gSizer1 = wx.GridSizer(n_names, 1, 0, 0 )

        staticTextA = wx.StaticText(self, wx.ID_ANY, 'Extract...', wx.DefaultPosition, wx.DefaultSize, 0)
        staticTextA.Wrap(-1)
        bold(staticTextA)

        #bSizer1.Add(staticTextA, 0, wx.ALL, 5)
        #gSizer1.Add(staticTextA, 0, wx.ALL, 5)


        self.m_scrolledWindow1 = wx.ScrolledWindow( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.EXPAND|wx.HSCROLL|wx.VSCROLL)
        self.m_scrolledWindow1.SetScrollRate( 5, 5 )

        for name in names:
            checkbox = wx.CheckBox( self.m_scrolledWindow1, wx.ID_ANY, name, wx.DefaultPosition, wx.DefaultSize, 0)
            checkbox.SetValue(True)

            gSizer1.Add(checkbox, 0, wx.ALL, 5)
            self.checkboxes.append(checkbox)

        self.m_scrolledWindow1.SetSizer( gSizer1)
        self.m_scrolledWindow1.Layout()
        gSizer1.Fit( self.m_scrolledWindow1 )

        ok_sizer = wx.StdDialogButtonSizer()

        ok_button = wx.Button(self, wx.ID_OK)
        apply_button = wx.Button(self, wx.ID_APPLY)
        cancel_button = wx.Button(self, wx.ID_CANCEL)

        apply_button.SetLabel('Select All')

        ok_sizer.AddButton(ok_button)
        ok_sizer.AddButton(apply_button)
        ok_sizer.AddButton(cancel_button)
        ok_sizer.Realize()

        bSizer1.Add(staticTextA, 0, wx.ALL, 5)
        bSizer1.Add( self.m_scrolledWindow1, 1, wx.EXPAND|wx.ALL, 1)
        bSizer1.Add(ok_sizer, 0, -1, 10)

        self.SetSizer(bSizer1)
        self.Layout()

        self.Centre(wx.BOTH)

        # Connect Events
        apply_button.Bind(wx.EVT_BUTTON, self.on_select_all)
        cancel_button.Bind(wx.EVT_BUTTON, self.on_cancel)
        ok_button.Bind(wx.EVT_BUTTON, self.on_ok)
        self.Show()

    def on_ok(self, event):
        names = self.args['names']
        names2 = []
        for name, checkbox in zip(names, self.checkboxes):
            val = checkbox.GetValue()
            if val:
                #print("name=%r val=%s" % (name, val))
                names2.append(name)
        self.args['selected'] = names2
        self.Destroy()

    def __del__(self):
        pass

    # Virtual event handlers, overide them in your derived class
    def on_select_all(self, event):
        for checkbox in self.checkboxes:
            checkbox.SetValue(True)

    def on_cancel(self, event):
        self.Destroy()

def main(names):
    ex = wx.App()

    args = {
        'names' : names,
        'selected' : None,
    }
    GetChoices(None, args)
    ex.MainLoop()
    #print("****args = ", args)
    selected = args['selected']
    print("****selected = %s" % selected)
    if selected is None:
        selected = names
    return selected

def get_choices(names=None):
    if names is None:
        names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'] * 5

    if is_wx:
        selected = main(names)
    else:
        msg = 'WxPython is supported, but PyQt and PySide are not.'
        raise NotImplementedError(msg)
    return selected

if __name__ == '__main__':  # pragma: no cover
    selected = get_choices()