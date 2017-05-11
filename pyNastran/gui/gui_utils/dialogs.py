from pyNastran.gui.qt_version import qt_version

if qt_version == 4:
    from PyQt4.QtGui import QFileDialog
elif qt_version == 5:
    from PyQt5.QtWidgets import QFileDialog


def save_file_dialog(self, title, default_filename,
                     file_types, filt):
    """common method for saving files"""
    if qt_version == 4:
        fname = str(QFileDialog.getSaveFileName(
            self, title, default_filename, file_types, filt))
        try:
            flt = str(filt).split()[0]
        except IndexError:
            flt = None
    else:
        fname, flt = QFileDialog.getSaveFileName(
            self, title, default_filename, file_types, filt)
        #flt = str(filt).strip()
    return fname, flt

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
