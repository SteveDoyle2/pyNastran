"""creates a popup that links to the new version of pyNastranGUI"""
import webbrowser

from pyNastran.gui.qt_version import qt_version, qt_int
from qtpy import QtCore, QtGui
from qtpy.QtWidgets import (
    QLabel, QApplication, QDialog, QGridLayout, QHBoxLayout, QVBoxLayout, QPushButton,
)

if qt_int == 4:
    class ClickableQLabel(QLabel):
        def __init(self, parent):
            QLabel.__init__(self, parent)

        def mouseReleaseEvent(self, event):
            if qt_int == 4:
                self.emit(QtCore.SIGNAL('clicked()'))
            else:
                # ????
                pass
elif qt_int == 5:
    pass
    #class ClickableQLabel(QPushButton):
        #def __init(self, text):
            #QPushButton.__init__(self, text)
            #self.setFlat(True)
else:  # pragma: no cover
    raise NotImplementedError('qt_version = %r' % qt_version)


class DownloadWindow(QDialog):
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

        QDialog.__init__(self, win_parent)
        self.setWindowTitle('pyNastran Update')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        #self.show()

    def create_widgets(self):
        self.name = QLabel("Version %s is now available." % self.version)
        if qt_int == 4:
            self.link = ClickableQLabel(self.url)
        else:
            self.link = QPushButton(self.url)
            self.link.setFlat(True)

        font = QtGui.QFont()
        #"Times",20,QtGui.QFont.Bold,True
        font.setUnderline(True)
        self.link.setFont(font)
        self.link.setStyleSheet("QLabel {color : blue}")

        # closing
        self.close_button = QPushButton("Close")

    def create_layout(self):
        grid = QGridLayout()
        grid.addWidget(self.name, 0, 0)
        grid.addWidget(self.link, 1, 0)

        close_box = QHBoxLayout()
        close_box.addWidget(self.close_button)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        vbox.addStretch()
        vbox.addLayout(close_box)
        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        if qt_int == 4:
            self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)
        self.link.clicked.connect(self.on_download)

        #self.link.linkActivated.connect(self.on_download)
        self.close_button.clicked.connect(self.on_cancel)

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


def main():  # pragma: no cover
    """test for the download menu"""
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    import pyNastran
    app = QApplication(sys.argv)
    url = pyNastran.__website__
    version = '1.0.0'
    main_window = DownloadWindow(url, version)
    main_window.show()
    app.exec_()

if __name__ == "__main__":  # pragma: no cover
    main()
