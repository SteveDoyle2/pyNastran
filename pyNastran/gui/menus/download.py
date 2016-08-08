import webbrowser
from PyQt4 import QtCore, QtGui


class ClickableQLabel(QtGui.QLabel):

    def __init(self, parent):
        QtGui.QLabel.__init__(self, parent)

    def mouseReleaseEvent(self, ev):
        self.emit(QtCore.SIGNAL('clicked()'))


class DownloadWindow(QtGui.QDialog):
    """
    +-------------------+
    | Legend Properties |
    +-----------------------+
    | Title  ______ Default |
    | Min    ______ Default |
    | Max    ______ Default |
    | Format ______ Default |
    | Scale  ______ Default |
    | Number of Colors ____ | (TODO)
    | Number of Labels ____ | (TODO)
    | Label Size       ____ | (TODO)
    |                       |
    | x Min/Max (Blue->Red) |
    | o Max/Min (Red->Blue) |
    |                       |
    | x Vertical/Horizontal |
    | x Show/Hide           |
    |                       |
    |    Apply OK Cancel    |
    +-----------------------+
    """

    def __init__(self, url, version, win_parent=None):
        self.win_parent = win_parent
        self.url = url
        self.version = version

        QtGui.QDialog.__init__(self, win_parent)
        self.setWindowTitle('pyNastran update ')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.show()

    def create_widgets(self):
        self.name = QtGui.QLabel("Version %s is now available." % self.version)
        self.link = ClickableQLabel(self.url)
        font = QtGui.QFont()
        #"Times",20,QtGui.QFont.Bold,True
        font.setUnderline(True)
        self.link.setFont(font)
        self.link.setStyleSheet("QLabel {color : blue}")

        # closing
        self.close_button = QtGui.QPushButton("Close")

    def create_layout(self):
        grid = QtGui.QGridLayout()
        grid.addWidget(self.name, 0, 0)
        grid.addWidget(self.link, 1, 0)

        close_box = QtGui.QHBoxLayout()
        close_box.addWidget(self.close_button)

        vbox = QtGui.QVBoxLayout()
        vbox.addLayout(grid)
        vbox.addStretch()
        vbox.addLayout(close_box)
        self.setLayout(vbox)

    def set_connections(self):
        self.connect(self.link, QtCore.SIGNAL('clicked()'), self.on_download)
        self.connect(self.close_button, QtCore.SIGNAL('clicked()'), self.on_cancel)
        self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.close()

    def closeEvent(self, event):
        event.accept()

    def on_download(self):
        webbrowser.open(self.url, new=0, autoraise=True)

    def on_ok(self):
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.close()


def main():
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    import pyNastran
    app = QtGui.QApplication(sys.argv)
    url = pyNastran.__website__
    version = '0.8.0'
    main_window = DownloadWindow(url, version)
    main_window.show()
    app.exec_()

if __name__ == "__main__":
    main()
