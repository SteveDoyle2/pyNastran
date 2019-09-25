from typing import Tuple
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
    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog
    fname, flt = QFileDialog.getOpenFileName(
        self, title, default_filename, file_types, options=options)
    #flt = str(filt).strip()
    return fname, flt

def save_file_dialog(self, title: str, default_dirname: str, file_types: str) -> Tuple[str, str]:
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
    assert isinstance(title, str), 'title=%s' % title
    assert isinstance(default_dirname, str), 'default_dirname=%s' % default_dirname
    assert isinstance(file_types, str), 'file_types=%s' % file_types
    #if qt_version == 5:
    # hasn't been tested
    out = QFileDialog.getSaveFileName(parent=self, caption=title,
                                      directory=default_dirname, filter=file_types)
    fname, wildcard_level = out
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
