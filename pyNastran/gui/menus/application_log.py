from __future__ import print_function
from pyNastran.gui.qt_version import qt_version, is_pygments

if qt_version == 4:
    from PyQt4 import QtCore
    from PyQt4.QtGui import (
        QDialog, QLabel, QLineEdit, QPushButton, QTextEdit, QDockWidget,
        QVBoxLayout, QHBoxLayout, QWidget, QFont, qApp, QMenu, QFontMetrics, QColor)
    try:
        import PyQt4.Qsci as Qsci
        is_scintilla = True
    except ImportError:
        is_scintilla = False

elif qt_version == 5:
    #from PyQt5 import QtCore, QtGui
    from PyQt5.QtGui import QFont, QFontMetrics, QColor
    from PyQt5 import QtCore
    #from PyQt5.QtCore import QSci
    #import QScintilla
    from PyQt5.QtWidgets import (
        QDialog, QLabel, QLineEdit, QPushButton, QTextEdit, QDockWidget,
        QVBoxLayout, QHBoxLayout, QWidget, qApp, QMenu)

    try:
        import PyQt5.Qsci as Qsci
        is_scintilla = True
    except ImportError:
        is_scintilla = False
    #import PyQt5.qsci
elif qt_version == 'pyside':
    from PySide import QtCore
    from PySide.QtGui import (
        QDialog, QLabel, QLineEdit, QPushButton, QTextEdit, QDockWidget,
        QVBoxLayout, QHBoxLayout, QWidget, QFont, qApp, QMenu, QFontMetrics, QColor)
    is_scintilla = False
else:
    raise NotImplementedError('qt_version = %r' % qt_version)


class HtmlLog(QTextEdit):

    def __init__(self, parent):
        QTextEdit.__init__(self, parent=parent)
        self.setReadOnly(True)

    def buttonClicked(self):
        if qApp.mouseButtons() & QtCore.Qt.RightButton:
            print(self.sender().toolTip())

class ApplicationLogWidget(QDockWidget):
    def __init__(self, parent=None):
        QDockWidget.__init__(self, 'Application log', parent=parent)
        self.setObjectName('application_log')
        self.log_widget = HtmlLog(parent=self)
        self.setWidget(self.log_widget)

#class QSyntaxHighlighting(Qsci.QsciScintilla):
    #def __init__(self):
        #Qsci.QsciScintilla.__init__(self)
        #self.setLexer(Qsci.QsciLexerPython(self))

        #font = QFont()
        #font.setFamily('Courier')
        #self.setFont(font)


if is_scintilla:
    print('using scintilla')
    class SimplePythonEditorWidget(Qsci.QsciScintilla):
        ARROW_MARKER_NUM = 8

        def __init__(self, parent=None):
            super(SimplePythonEditorWidget, self).__init__(parent)

            # Set the default font
            font = QFont()
            font.setFamily('Courier')
            font.setFixedPitch(True)
            font.setPointSize(10)
            self.setFont(font)
            self.setMarginsFont(font)

            # Margin 0 is used for line numbers
            fontmetrics = QFontMetrics(font)
            self.setMarginsFont(font)
            self.setMarginWidth(0, fontmetrics.width("00000") + 6)
            self.setMarginLineNumbers(0, True)
            self.setMarginsBackgroundColor(QColor("#cccccc"))

            # Clickable margin 1 for showing markers
            self.setMarginSensitivity(1, True)
            self.connect(self,
                         QtCore.SIGNAL('marginClicked(int, int, Qt::KeyboardModifiers)'),
                         self.on_margin_clicked)
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
            self.SendScintilla(Qsci.QsciScintilla.SCI_STYLESETFONT, 1, 'Courier')

            # Don't want to see the horizontal scrollbar at all
            # Use raw message to Scintilla here (all messages are documented
            # here: http://www.scintilla.org/ScintillaDoc.html)
            self.SendScintilla(Qsci.QsciScintilla.SCI_SETHSCROLLBAR, 0)

            # not too small
            #self.setMinimumSize(600, 450)

        def on_margin_clicked(self, nmargin, nline, modifiers):
            # Toggle marker for the line the margin was clicked on
            if self.markersAtLine(nline) != 0:
                self.markerDelete(nline, self.ARROW_MARKER_NUM)
            else:
                self.markerAdd(nline, self.ARROW_MARKER_NUM)

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
        #super(QDockWidget, self).__init__('Python Console', parent=parent) # I think this works by accident in qt4/5

        self.listMenu = QMenu()
        self.execute_python_button = QPushButton('Execute')
        self.execute_and_clear_python_button = QPushButton('Execute and Clear')

        if is_pygments and is_scintilla:
            #self.enter_data = QSyntaxHighlighting()
            self.enter_data = SimplePythonEditorWidget()
        else:
            self.enter_data = QTextEdit()
            font = QFont()
            font.setFamily('Courier')
            self.enter_data.setFont(font)
        self.setup_connections()
        self.layout()

    def layout(self):
        vbox = QVBoxLayout()
        hbox = QHBoxLayout()

        vbox.addWidget(self.enter_data)
        hbox.addWidget(self.execute_python_button)
        hbox.addWidget(self.execute_and_clear_python_button)
        vbox.addLayout(hbox)

        vbox_widget = QWidget()
        vbox_widget.setLayout(vbox)
        self.setWidget(vbox_widget)

    def setup_connections(self):
        self.execute_python_button.clicked.connect(self.on_execute_python_button)
        self.execute_and_clear_python_button.clicked.connect(
            self.on_execute_and_clear_python_button)
        #self.on_rig

        # TODO: enables the right click menu
        #       messes up the previous right click menu
        #self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #self.connect(self, QtCore.SIGNAL("customContextMenuRequested(QPoint)"),
        #             self.listItemRightClicked)

        menu_item1 = self.listMenu.addAction("Properties...")
        menu_item2 = self.listMenu.addAction("Select All...")
        menu_item3 = self.listMenu.addAction("Copy...")
        menu_item1.triggered.connect(self.menuItemClicked_1)
        menu_item2.triggered.connect(self.menuItemClicked_2)
        menu_item3.triggered.connect(self.menuItemClicked_3)
        #self.connect(menu_item1, QtCore.SIGNAL("triggered()"), self.menuItemClicked_1)
        #self.connect(menu_item2, QtCore.SIGNAL("triggered()"), self.menuItemClicked_2)
        #self.connect(menu_item3, QtCore.SIGNAL("triggered()"), self.menuItemClicked_3)

        # we have to create a QWidget to put the console vbox into vbox_widget because
        #    self.setLayout(vbox)
        # is not supported in a QDockWidget

    def on_execute_and_clear_python_button(self):
        self.parent._on_execute_python_button(clear=True)

    def on_execute_python_button(self):
        self.parent._on_execute_python_button(clear=False)

    def listItemRightClicked(self, QPos):
        """
        TDOO: not done...
        """
        return
        #parentPosition = self.mapToGlobal(QtCore.QPoint(0, 0))
        #self.listMenu.move(parentPosition + QPos)
        #self.listMenu.show()

    def menuItemClicked_1(self):
        print(1)

    def menuItemClicked_2(self):
        print(2)

    def menuItemClicked_3(self):
        print(3)
