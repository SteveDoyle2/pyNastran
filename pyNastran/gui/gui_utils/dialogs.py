from pyNastran.gui.qt_version import qt_version

from qtpy.QtWidgets import QFileDialog

def open_file_dialog(self, title, default_filename,
                     file_types):
    """common method for opening files"""
    if qt_version == 4:
        fname, wildcard_level = QFileDialog.getOpenFileNameAndFilter(
            self, title, default_filename, file_types)
        return str(fname), str(wildcard_level)
    else:
        fname, flt = QFileDialog.getOpenFileName(
            self, title, default_filename, file_types)
        #flt = str(filt).strip()
    return fname, flt

def open_directory_dialog(self, title):
    """common method for selecting a directory"""
    dirname = str(QFileDialog.getExistingDirectory(self, title))
    return dirname
