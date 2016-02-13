import os

try:
    import wx
    _gui_mode = 0
except ImportError:
    try:
        from PySide import QtCore, QtGui
        _gui_mode = 1
    except ImportError:
        try:
            from PyQt4 import QtCore, QtGui
            _gui_mode = 2
        except ImportError:
            _gui_mode = None

if _gui_mode == 0:
    pass
elif _gui_mode in [1, 2]:
    class QtDialog(QtGui.QWidget):
        """Dummy GUI Object"""
        def __init__(self):
            # super(DialogDemo, self).__init__()
            QtGui.QWidget.__init__(self)

            # set the position and size of the window
            # make it really small, so we don't see this dummy window
            self.setGeometry(1, 1, 0, 0)

#----------------------------------------------------------------------
#print("f_mode =", f_mode)

def radio_pullown_dialog(Title, button_dict, nwide=3):  # pragma: no conver
    """
    buttons = [
        [header, 'checkbox', A', 'B', 'C'],
        [header2, 'pulldown', 'A', 'B'],
    ]

    header
    ------
    A     B     C
    x

    header2 v
    """
    pass

def save_file_dialog(Title, wx_wildcard, qt_wildcard, dirname=''):
    """
    creates a save file dialog in wx or PyQt4/PySide
    """
    fname = None
    if dirname == '':
        dirname = os.getcwd()

    if _gui_mode == 0: # wx
        app = wx.App(redirect=False)
        app.MainLoop()
        dlg = wx.FileDialog(None, Title, dirname, "",
                            wx_wildcard, wx.SAVE)
        app.MainLoop()

        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            fname = os.path.join(dirname, filename)

    elif _gui_mode in [1, 2]:  # PySide, PyQt4
        # checks if QApplication already exists
        app = QtGui.QApplication.instance()
        if not app: # create QApplication if it doesnt exist
            app = QtGui.QApplication([])
        form = QtDialog()
        form.show()

        fname = QtGui.QFileDialog.getSaveFileName(form, Title,
                                                  dirname, qt_wildcard)
        app.exit()
        #print("fname =", fname)
    else:
        msg = 'Could not import wx, PySide, or PyQt4.  ' \
            'Please specify the file explicitly.'
        raise ImportError(msg)
    return fname


def load_file_dialog(Title, wx_wildcard, qt_wildcard, dirname=''):
    """
    creates a load file dialog in wx or PyQt4/PySide
    """
    fname = None
    if dirname == '':
        dirname = os.getcwd()

    wildcard_level = None
    if _gui_mode == 0: # wx
        app = wx.App(redirect=False)
        app.MainLoop()
        dlg = wx.FileDialog(None, Title, dirname, "",
                            wx_wildcard, wx.OPEN)
        app.MainLoop()

        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            fname = os.path.join(dirname, filename)

    elif _gui_mode in [1, 2]:  # PySide, PyQt4
        # checks if QApplication already exists
        app = QtGui.QApplication.instance()
        if not app: # create QApplication if it doesnt exist
            app = QtGui.QApplication([])
        form = QtDialog()
        form.show()

        output = QtGui.QFileDialog.getOpenFileName(form, Title,
            dirname, qt_wildcard)
        if len(output) == 1:
            fname, wildcard_level = output
        else:
            fname = output
        app.exit()
    else:
        msg = 'Could not import wx, PySide, or PyQt4.  '\
            'Please specify the file explicitly.'
        raise ImportError(msg)
    return fname, wildcard_level


def _main():  # pragma: no conver
    """helps to test the functions"""
    wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas)|*.bdf;*.dat;*.nas|" \
        "All files (*.*)|*.*"
    wildcard_qt = "Nastran BDF (*.bdf *.dat *.pch);;All files (*)"
    stitle = 'Please select a BDF/DAT/PCH to load'
    fname = load_file_dialog(stitle, wildcard_wx, wildcard_qt)
    print("fname2 = %s" % fname)


if __name__ == '__main__':  # pragma: no conver
    _main()
