from PyQt4 import QtCore, QtGui


class MultiFileDialog(QtGui.QFileDialog):
    def __init__(self, *args, **kwargs):
        super(MultiFileDialog, self).__init__(*args, **kwargs)
        self.setOption(QtGui.QFileDialog.DontUseNativeDialog, True)
        self.setFileMode(QtGui.QFileDialog.ExistingFiles)

    def accept(self):
        QtGui.QDialog.accept(self)

class Window(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.button = QtGui.QPushButton('Test', self)
        self.button.clicked.connect(self.handleButton)
        layout = QtGui.QVBoxLayout(self)
        layout.addWidget(self.button)

    def handleButton(self):
        Title = 'Load a Tecplot Geometry/Results File'
        last_dir = ''
        qt_wildcard = ['Tecplot Hex Binary (*.tec; *.py)']
        dialog = MultiFileDialog()
        dialog.setWindowTitle(Title)
        #dialog.selectedFilter(qt_wildcard)
        dialog.setFilters(qt_wildcard)
        #dialog = MultiFileDialog(Title, last_dir, qt_wildcard)
        if dialog.exec_() == QtGui.QDialog.Accepted:
            print(dialog.selectedFiles())
            sfilter = str(dialog.selectedFilter())
            print(sfilter)

if __name__ == '__main__':  # pragma: no cover
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
