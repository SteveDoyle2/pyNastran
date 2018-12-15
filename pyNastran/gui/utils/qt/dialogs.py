from pyNastran.gui.qt_version import qt_version

from qtpy.QtWidgets import QFileDialog

def open_file_dialog(self, title, default_filename,
                     file_types):
    """
    Common method for opening files


    Parameters
    ----------
    self : ???
        the main GUI; not a vbox
    title : str
        the title of the dialog
    default_filename : str
        the default directory
    file_types : str
        the wildcard
        'Nastran Geometry - Punch (*.bdf; *.dat; *.nas; *.ecd; *.pch);;All files (*)'
    """
    if qt_version == 4:
        # works in: pyqt4, pyside
        # doesn't work in: pyqt5
        fname, wildcard_level = QFileDialog.getOpenFileNameAndFilter(
            self, title, default_filename, file_types)
        return str(fname), str(wildcard_level)
    else:
        fname, flt = QFileDialog.getOpenFileName(
            self, title, default_filename, file_types)
        #flt = str(filt).strip()
    return fname, flt

def save_file_dialog(self, title, default_dirname, file_types):
    """
    Common method for saving files


    Parameters
    ----------
    self : ???
        the main GUI; not a vbox
    title : str
        the title of the dialog
    default_dirname : str
        the default directory
    file_types : str
        the wildcard
        'Nastran Geometry - Punch (*.bdf; *.dat; *.nas; *.ecd; *.pch);;All files (*)'
    """
    #asdf
    #if qt_version == 5:
    # hasn't been tested
    fname, wildcard_level = QFileDialog.getSaveFileName(self, title, default_dirname, file_types)
    return str(fname), str(wildcard_level)
    #else:
        #fname, flt = QFileDialog.getOpenFileName(
            #self, title, default_filename, file_types)
        #flt = str(filt).strip()
    #return fname, flt

def open_directory_dialog(self, title, directory=''):
    """
    Common method for selecting a directory

    Parameters
    ----------
    self : ???
        the main GUI; not a vbox
    title : str
        the title of the dialog
    directory : str
        the default directory
    """
    dirname = str(QFileDialog.getExistingDirectory(self, caption=title, directory=directory))
    return dirname
