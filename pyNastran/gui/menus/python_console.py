from pyNastran.gui.qt_version import qt_version, is_pygments

from qtpy.QtGui import QFont, QFontMetrics, QColor
from qtpy import QtCore
#from qtpy.QtCore import QSci
#import QScintilla
from qtpy.QtWidgets import (
    QPushButton, QTextEdit, QDockWidget,
    QVBoxLayout, QHBoxLayout, QWidget)
#import qtpy.Qsci as Qsci

QSCINTILLA_VERSION = 'N/A'
if qt_version == 'pyqt5':
    try:
        import PyQt5.Qsci as Qsci
        QSCINTILLA_VERSION = Qsci.QSCINTILLA_VERSION_STR
        IS_SCINTILLA = True
    except ImportError:
        IS_SCINTILLA = False
elif qt_version == 'pyside2':
    IS_SCINTILLA = False
else:
    raise NotImplementedError('qt_version = %r' % qt_version)

# hack for pyinstaller...if QScintilla isn't working in the exe, but it is in
# gui, enable this print block and diff the outputs
#
#import sys
#with open('modules.txt', 'w') as f:
    #mods = sorted(sys.modules.keys())
    #f.write('\n'.join(mods))


#class QSyntaxHighlighting(Qsci.QsciScintilla):
    #def __init__(self):
        #Qsci.QsciScintilla.__init__(self)
        #self.setLexer(Qsci.QsciLexerPython(self))

        #font = QFont()
        #font.setFamily('Courier')
        #self.setFont(font)


if IS_SCINTILLA:
    class SimplePythonEditorWidget(Qsci.QsciScintilla):
        ARROW_MARKER_NUM = 8

        def __init__(self, parent=None):
            """
            very similar to:
            https://stackoverflow.com/questions/40002373/qscintilla-based-text-editor-in-pyqt5-with-clickable-functions-and-variables
            """
            super(SimplePythonEditorWidget, self).__init__(parent)

            # Set the default font
            font = QFont()
            font.setFamily('Courier')
            font.setFixedPitch(True)
            font.setPointSize(10)
            self.setFont(font)
            self.setMarginsFont(font)
            self.set_font(font)

            self.setMarginLineNumbers(0, True)
            self.setMarginsBackgroundColor(QColor("#cccccc"))

            # Clickable margin 1 for showing markers
            self.setMarginSensitivity(1, True)
            self.marginClicked.connect(self.on_margin_clicked)

            self.markerDefine(Qsci.QsciScintilla.RightArrow,
                              self.ARROW_MARKER_NUM)
            self.setMarkerBackgroundColor(QColor("#ee1111"),
                                          self.ARROW_MARKER_NUM)

            # Brace matching: enable for a brace immediately before or after
            # the current position
            self.setBraceMatching(Qsci.QsciScintilla.SloppyBraceMatch)

            # Current line visible with special background color
            self.setCaretLineVisible(True)
            self.setCaretLineBackgroundColor(QColor("#ffe4e4"))

            # Set Python lexer
            # Set style for Python comments (style number 1) to a fixed-width
            # courier.
            lexer = Qsci.QsciLexerPython()
            lexer.setDefaultFont(font)
            self.setLexer(lexer)

            font_style = bytearray(str.encode("Courier"))
            self.SendScintilla(Qsci.QsciScintilla.SCI_STYLESETFONT, 1, font_style)

            # Don't want to see the horizontal scrollbar at all
            # Use raw message to Scintilla here (all messages are documented
            # here: http://www.scintilla.org/ScintillaDoc.html)
            self.SendScintilla(Qsci.QsciScintilla.SCI_SETHSCROLLBAR, 0)

            # not too small
            #self.setMinimumSize(600, 450)

        def set_font(self, font):
            # Margin 0 is used for line numbers
            fontmetrics = QFontMetrics(font)
            self.setMarginsFont(font)
            self.setMarginWidth(0, fontmetrics.width("00000") + 6)

        def on_margin_clicked(self, nmargin, nline, modifiers):
            return
            # Toggle marker for the line the margin was clicked on
            #if self.markersAtLine(nline) != 0:
                #self.markerDelete(nline, self.ARROW_MARKER_NUM)
            #else:
                #self.markerAdd(nline, self.ARROW_MARKER_NUM)

        def toPlainText(self):
            data = str(self.text())
            return data

class PythonConsoleWidget(QDockWidget):
    """
    original code pulled from:
    http://stackoverflow.com/questions/31380457/add-right-click-functionality-to-listwidget-in-pyqt4
    http://stackoverflow.com/questions/20951660/creating-a-syntax-highlighter-in-pythonpyqt4
    http://eli.thegreenplace.net/2011/04/01/sample-using-qscintilla-with-pyqt
    """
    def __init__(self, parent):
        self.parent = parent
        super(PythonConsoleWidget, self).__init__('Python Console', parent=parent)

        # I think this works by accident in qt4/5
        #super(QDockWidget, self).__init__('Python Console', parent=parent)

        self.execute_python_button = QPushButton('Execute')
        self.execute_and_clear_python_button = QPushButton('Execute and Clear')

        self.enter_data = get_code_block()

        self.setup_connections()
        self.layout()

    def layout(self):
        vbox = QVBoxLayout()
        hbox = QHBoxLayout()

        vbox.addWidget(self.enter_data)
        hbox.addWidget(self.execute_python_button)
        hbox.addWidget(self.execute_and_clear_python_button)
        vbox.addLayout(hbox)

        vbox_widget = layout_to_widget(vbox)
        self.setWidget(vbox_widget)

    def setup_connections(self):
        """sets up the callbacks"""
        self.execute_python_button.clicked.connect(self.on_execute_python_button)
        self.execute_and_clear_python_button.clicked.connect(
            self.on_execute_and_clear_python_button)
        #self.on_rig

        # TODO: enables the right click menu
        #       messes up the previous right click menu
        #self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #self.connect(self, QtCore.SIGNAL("customContextMenuRequested(QPoint)"),
        #             self.listItemRightClicked)

        #properties = self.listMenu.addAction("Properties...")
        #select_all = self.listMenu.addAction("Select All")
        #copy = self.listMenu.addAction("Copy")
        #clear = self.listMenu.addAction("Clear...")
        #properties.triggered.connect(self.on_properties)
        #select_all.triggered.connect(self.on_select_all)
        #copy.triggered.connect(self.on_copy)
        #clear.triggered.connect(self.on_clear)

        #self.connect(properties, QtCore.SIGNAL("triggered()"), self.on_properties)
        #self.connect(select_all, QtCore.SIGNAL("triggered()"), self.on_select_all)
        #self.connect(copy, QtCore.SIGNAL("triggered()"), self.on_copy)

        # we have to create a QWidget to put the console vbox into vbox_widget because
        #    self.setLayout(vbox)
        # is not supported in a QDockWidget

    def on_execute_and_clear_python_button(self):
        if self.parent is not None:
            self.parent._on_execute_python_button(clear=True)

    def on_execute_python_button(self):
        if self.parent is not None:
            self.parent._on_execute_python_button(clear=False)

    def listItemRightClicked(self, QPos):
        """
        TDOO: not done...
        """
        return
        #parentPosition = self.mapToGlobal(QtCore.QPoint(0, 0))
        #self.listMenu.move(parentPosition + QPos)
        #self.listMenu.show()

    #def on_properties(self):
        #print(1)

    #def on_select_all(self):
        #print(2)

    #def on_copy(self):
        #print(3)

    #def on_clear(self):
        #print(4)

def get_code_block():
    if is_pygments and IS_SCINTILLA:
        #self.enter_data = QSyntaxHighlighting()
        enter_data = SimplePythonEditorWidget()
    else:
        enter_data = QTextEdit()
        font = QFont()
        font.setFamily('Courier')
        enter_data.setFont(font)
    return enter_data

def layout_to_widget(layout):
    widget = QWidget()
    widget.setLayout(layout)
    return widget


def main():  # pragma: no cover
    """tests function that doesn't bleed over"""
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    from qtpy.QtWidgets import QApplication

    import sys
    #import pyNastran
    app = QApplication(sys.argv)
    main_window = PythonConsoleWidget(None)
    main_window.show()
    app.exec_()


if __name__ == '__main__':  # pragma: no cover
    main()
