import os
import typing

#try:
    #import wx  # type: ignore
    #_gui_mode = 'wx'

try:
    from PySide import QtCore, QtGui  # type: ignore
    from PySide.QtGui import QWidget, QApplication, QFileDialog  # type: ignore
    _gui_mode = 'pyside'
except ImportError:
    try:
        from PyQt4 import QtCore, QtGui  # type: ignore
        from PyQt4.QtGui import QWidget, QApplication, QFileDialog  # type: ignore
        _gui_mode = 'pyqt'
    except ImportError:
        try:
            from PyQt5 import QtCore, QtGui  # type: ignore
            from PyQt5.QtWidgets import QWidget, QApplication, QFileDialog  # type: ignore
            _gui_mode = 'pyqt'
        except ImportError:
            try:
                import wx  # type: ignore
                _gui_mode = 'wx'
            except ImportError:
                _gui_mode = None

if _gui_mode == 'wx':
    pass
elif _gui_mode in ['pyqt', 'pyside']:
    class QtDialog(QWidget):
        """Dummy GUI Object"""
        def __init__(self):
            # super(DialogDemo, self).__init__()
            QWidget.__init__(self)

            # set the position and size of the window
            # make it really small, so we don't see this dummy window
            self.setGeometry(1, 1, 0, 0)

#----------------------------------------------------------------------
#print("f_mode =", f_mode)

def radio_pullown_dialog(Title, button_dict, nwide=3):  # pragma: no cover
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

def save_file_dialog(title, wx_wildcard, qt_wildcard, dirname=''):
    # type: (str, str, str, str) -> str
    """
    creates a save file dialog in wx or PyQt4/PySide
    """
    fname = None
    if dirname == '':
        dirname = os.getcwd()

    if _gui_mode == 'wx':
        app = wx.App(redirect=False)
        app.MainLoop()
        dlg = wx.FileDialog(None, title, dirname, "",
                            wx_wildcard, wx.SAVE)
        app.MainLoop()

        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            fname = os.path.join(dirname, filename)

    elif _gui_mode in ['pyqt', 'pyside']:
        # checks if QApplication already exists
        app = QApplication.instance()
        if not app: # create QApplication if it doesnt exist
            app = QApplication([])
        form = QtDialog()
        form.show()

        fname = QFileDialog.getSaveFileName(form, title,
                                            dirname, qt_wildcard)
        app.exit()
        #print("fname =", fname)
    else:
        msg = 'Could not import wx, PySide, or PyQt4.  ' \
            'Please specify the file explicitly.'
        raise ImportError(msg)
    return fname


def load_file_dialog(title, wx_wildcard, qt_wildcard, dirname=''):
    # type: (str, str, str, str) -> (str, str)
    """
    creates a load file dialog in wx or PyQt4/PySide
    """
    fname = None
    if dirname == '':
        dirname = os.getcwd()

    wildcard_level = None
    assert isinstance(_gui_mode, str), _gui_mode
    if _gui_mode == 'wx': # wx
        app = wx.App(redirect=False)
        app.MainLoop()
        dlg = wx.FileDialog(None, title, dirname, "",
                            wx_wildcard, wx.OPEN)
        app.MainLoop()

        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            fname = os.path.join(dirname, filename)

    elif _gui_mode in ['pyqt', 'pyside']:
        # checks if QApplication already exists
        app = QApplication.instance()
        if not app: # create QApplication if it doesnt exist
            app = QApplication([])
        form = QtDialog()
        form.show()

        output = QFileDialog.getOpenFileName(form, title,
            dirname, qt_wildcard)
        if len(output) == 1:
            fname, wildcard_level = output
        else:
            fname = output
        app.exit()
    else:
        msg = 'Could not import wx, PySide, or PyQt4/5.  '\
            'Please specify the file explicitly.' + '\ngui_mode=%r' % _gui_mode
        raise ImportError(msg)
    return fname, wildcard_level


def main():  # pragma: no cover
    """helps to test the functions"""
    wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas)|*.bdf;*.dat;*.nas|" \
        "All files (*.*)|*.*"
    wildcard_qt = "Nastran BDF (*.bdf *.dat *.pch);;All files (*)"
    stitle = 'Please select a BDF/DAT/PCH to load'
    fname = load_file_dialog(stitle, wildcard_wx, wildcard_qt)
    print("fname2 = %s" % fname)


if __name__ == '__main__':  # pragma: no cover
    main()
