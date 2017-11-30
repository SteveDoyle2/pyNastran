from pyNastran.gui.qt_version import qt_version

from qtpy.QtWidgets import QFileDialog


def save_file_dialog(self, title, default_filename,
                     file_types, filt):
    """common method for saving files"""
    fname, flt = QFileDialog.getSaveFileName(
        self, title, default_filename, file_types, filt)
    return fname, flt

def open_file_dialog(self, title, default_filename,
                     file_types):
    """common method for opening files"""
    fname, flt = QFileDialog.getOpenFileName(
        self, title, default_filename, file_types)
    return fname, flt

def open_directory_dialog(self, title):
    """common method for selecting a directory"""
    dirname = str(QFileDialog.getExistingDirectory(self, title))
    return dirname
