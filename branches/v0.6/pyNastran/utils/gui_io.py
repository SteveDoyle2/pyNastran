import os

try:
    import wx
    fmode = 0
except ImportError:
    try:
        from PySide import QtCore, QtGui
        fmode = 1
    except ImportError:
        try:
            from PyQt4 import QtCore, QtGui
            fmode = 2
        except ImportError:
            fmode = None

if fmode == 0:
    pass
elif fmode in [1, 2]:
    class QtDialog(QtGui.QWidget):
        def __init__(self):
            # super(DialogDemo, self).__init__()
            QtGui.QWidget.__init__(self)

            # set the position and size of the window
            # make it really small, so we don't see this dummy window
            self.setGeometry(1, 1, 0, 0)
 
#----------------------------------------------------------------------
#print "fmode =", fmode
def load_file_dialog(Title, wx_wildcard, qt_wildcard, dirname=''):
    fname = None
    if dirname == '':
        dirname = os.getcwd()

    if fmode == 0: # wx
        app = wx.App(redirect=False)
        app.MainLoop()
        dlg = wx.FileDialog(None, Title, dirname, "",
                            wx_wildcard, wx.OPEN)
        app.MainLoop()

        if dlg.ShowModal() == wx.ID_OK:
            fileName = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            fname = os.path.join(dirname, fileName)

    elif fmode in [1, 2]:  # PySide, PyQt4
        app = QtGui.QApplication.instance() # checks if QApplication already exists 
	if not app: # create QApplication if it doesnt exist 
	    app = QtGui.QApplication([])
        form = QtDialog()
        form.show()

        fname, wildcard_level = QtGui.QFileDialog.getOpenFileName(form, Title,
            dirname, qt_wildcard)
        app.exit()
        #print "fname =", fname
    else:
        raise ImportError('Could not import wx, PySide, or PyQt4.  Please specify the file explicitly.')
    return fname

if __name__ == '__main__':
    wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas)|*.bdf;*.dat;*.nas|" \
        "All files (*.*)|*.*"
    wildcard_qt = "Nastran BDF (*.bdf *.dat *.pch);;All files (*)"
    stitle = 'Please select a BDF/DAT/PCH to load'
    fname = load_file_dialog(stitle, wildcard_wx, wildcard_qt)
    print "fname2 = ", fname

