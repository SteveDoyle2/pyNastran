from pyNastran.gui.qt_version import qt_version

from qtpy.QtWidgets import QFileDialog, QWidget


def open_file_dialog(self, title: str, default_filename: str,
                     file_types: str) -> tuple[str, str]:
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

    Returns
    -------
    fname : str
        the selected file
    flt : str
        the filter selection, '(*.bdf; *.dat; *.nas; *.ecd; *.pch)'

    """
    if qt_version == 'pyqt6':
        fname, flt = QFileDialog.getOpenFileName(
            self, title, default_filename, file_types)
    else:
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fname, flt = QFileDialog.getOpenFileName(
            self, title, default_filename, file_types, options=options)

    #flt = str(filt).strip()
    return fname, flt

def save_file_dialog(self: QWidget, title: str, default_dirname: str,
                     file_types: str) -> tuple[str, str]:
    """
    Common method for saving files

    Parameters
    ----------
    self : QWidget
        the main GUI; not a vbox
    title : str
        the title of the dialog
    default_dirname : str
        the default directory
    file_types : str
        the wildcard
        'Nastran Geometry - Punch (*.bdf; *.dat; *.nas; *.ecd; *.pch);;All files (*)'

    Returns
    -------
    fname : str
        the saved filename (file.bdf)
    wildcard_level : str
        the extension of the file (.bdf)
    """
    assert isinstance(title, str), f'title={title!r}'
    assert isinstance(default_dirname, str), f'default_dirname={default_dirname!r}'
    assert isinstance(file_types, str), f'file_types={file_types}'

    # tested in pyside2/pyside6
    if qt_version in {'pyside6', 'pyside2'}:
        out = QFileDialog.getSaveFileName(parent=self, caption=title,
                                          dir=default_dirname, filter=file_types)
        # selected_filter; options
    else:
        # pyside2
        out = QFileDialog.getSaveFileName(parent=self, caption=title,
                                          directory=default_dirname, filter=file_types)
    fname, wildcard_level = out
    return str(fname), str(wildcard_level)
    #else:
        #fname, flt = QFileDialog.getOpenFileName(
            #self, title, default_filename, file_types)
        #flt = str(filt).strip()
    #return fname, flt

def open_directory_dialog(self, title: str, directory: str='') -> str:
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
    # tested in pyside2/pyside6
    if qt_version == 'pyside6':
        dirname = str(QFileDialog.getExistingDirectory(self, caption=title, dir=directory))
    else:
        # pyside2
        dirname = str(QFileDialog.getExistingDirectory(self, caption=title, directory=directory))
    return dirname
