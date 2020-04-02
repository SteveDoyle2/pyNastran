from pyNastran.gui.qt_version import qt_version, is_pygments

from qtpy.QtCore import Qt
from qtpy.QtGui import QFont, QFontMetrics, QColor, QCursor
from qtpy import QtCore
from qtpy.QtWidgets import QTextEdit, QDockWidget


class HtmlLog(QTextEdit):

    def __init__(self, *args, **kwargs):
        super(HtmlLog, self).__init__(*args, **kwargs)
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.__contextMenu)

        # https://stackoverflow.com/questions/3120258/qtextedit-inserthtml-is-very-slow
        self.setReadOnly(True)
        self.setAcceptRichText(False)
        self.setContextMenuPolicy(Qt.NoContextMenu)
        #self.setOpenLinks(False)
        self.setUndoRedoEnabled(False)

    def __contextMenu(self):
        self._normalMenu = self.createStandardContextMenu()
        self._addCustomMenuItems(self._normalMenu)
        self._normalMenu.exec_(QCursor.pos())

    def _addCustomMenuItems(self, menu):
        menu.addSeparator()
        menu.addAction('Clear...', self.on_clear)

    def on_clear(self):
        """calls out to the main GUI to clear the screen"""
        self.parent().on_clear()

    #def buttonClicked(self):  # PySide2 must be <5.14.2
        #if qApp.mouseButtons() & QtCore.Qt.RightButton:
            #print(self.sender().toolTip())

    def mousePressEvent(self, event):
        """called when you click a result"""
        if event.button() == Qt.RightButton:
            pass
        else:
            QTextEdit.mousePressEvent(self, event)

    def clear(self):
        """clears out the text"""
        self.setText('')

    def __repr__(self):
        return 'HtmlLog()'

class ApplicationLogWidget(QDockWidget):
    def __init__(self, parent=None):
        QDockWidget.__init__(self, 'Application log', parent=parent)
        self.setObjectName('application_log')
        self.log_widget = HtmlLog(parent=self)
        self.setWidget(self.log_widget)

    def on_clear(self):
        parent = self.parent()
        if hasattr(parent, 'clear_application_log'):
            parent.clear_application_log()

    #def keyPressEvent(self, event):
        #key = event.key()
        #if key == Qt.Key_Delete:
            #index = self.currentIndex()
            #self.parent().on_delete(index.row())
            #print('pressed delete')


def main():  # pragma: no cover
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    from qtpy.QtWidgets import QApplication

    import sys
    import pyNastran
    app = QApplication(sys.argv)
    url = pyNastran.__website__
    #version = '1.0.0'
    main_window = ApplicationLogWidget()
    log_widget = main_window.log_widget
    msg = 'This is a message asdf'
    text_cursor = log_widget.textCursor()
    text_cursor.insertHtml(msg + r"<br />")
    #log_widget.clear()
    main_window.show()
    app.exec_()


if __name__ == '__main__':  # pragma: no cover
    main()
