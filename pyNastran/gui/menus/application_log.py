from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4 import QtCore
    from PyQt4.QtGui import (
        QDialog, QLabel, QLineEdit, QPushButton, QTextEdit, QDockWidget,
        QVBoxLayout, QHBoxLayout, QWidget)
elif qt_version == 5:
    #from PyQt5 import QtCore, QtGui
    from PyQt5.QtWidgets import (
        QDialog, QLabel, QLineEdit, QPushButton, QTextEdit, QDockWidget,
        QVBoxLayout, QHBoxLayout, QWidget)

class HtmlLog(QTextEdit):

    def __init__(self, parent):
        QTextEdit.__init__(self, parent=parent)
        self.setReadOnly(True)

    def buttonClicked(self):
        if QtGui.qApp.mouseButtons() & QtCore.Qt.RightButton:
            print(self.sender().toolTip())

class ApplicationLogWidget(QDockWidget):
    def __init__(self, parent=None):
        QDockWidget.__init__(self, 'Application log', parent=parent)
        self.setObjectName('application_log')
        self.log_widget = HtmlLog(parent=self)
        self.setWidget(self.log_widget)

class PythonConsoleWidget(QDockWidget):
    def __init__(self, parent):
        self.parent = parent
        super(QDockWidget, self).__init__('Python Console', parent=parent)
        vbox = QVBoxLayout()
        hbox = QHBoxLayout()

        self.enter_data = QTextEdit()
        self.execute_python_button = QPushButton('Execute')
        self.execute_and_clear_python_button = QPushButton('Execute and Clear')
        vbox.addWidget(self.enter_data)
        hbox.addWidget(self.execute_python_button)
        hbox.addWidget(self.execute_and_clear_python_button)
        vbox.addLayout(hbox)


        if qt_version == 4:
            self.connect(self.execute_python_button, QtCore.SIGNAL('clicked()'), self.on_execute_python_button)
            self.connect(self.execute_and_clear_python_button, QtCore.SIGNAL('clicked()'), self.on_execute_and_clear_python_button)
        else:
            self.execute_python_button.clicked.connect(self.on_execute_python_button)
            self.execute_and_clear_python_button.clicked.connect(self.on_execute_and_clear_python_button)

        # we have to create a QWidget to put the console vbox into vbox_widget because
        #    self.setLayout(vbox)
        # is not supported in a QDockWidget
        vbox_widget = QWidget()
        vbox_widget.setLayout(vbox)
        self.setWidget(vbox_widget)

    def on_execute_and_clear_python_button(self):
        self.parent._on_execute_python_button(clear=True)

    def on_execute_python_button(self):
        self.parent._on_execute_python_button(clear=False)
