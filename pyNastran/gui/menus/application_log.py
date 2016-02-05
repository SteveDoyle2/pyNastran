from PyQt4 import QtCore, QtGui


class ApplicationLogDockWidget(QtGui.QDockWidget):
    def __init__(self, parent, execute_python=True):
        self.parent = parent
        QtGui.QDockWidget.__init__(self, parent)
        self.dock_widget = QtGui.QDockWidget("Application log", parent)
        self.dock_widget.setObjectName("application_log")
        self.log_widget = QtGui.QTextEdit()
        self.log_widget.setReadOnly(True)

        if execute_python:
            vbox1 = QtGui.QVBoxLayout()
            vbox2 = QtGui.QVBoxLayout()
            hbox = QtGui.QHBoxLayout()

            vbox1w = QtGui.QWidget()
            vbox1w.setLayout(vbox1)

            vbox2w = QtGui.QWidget()
            vbox2w.setLayout(vbox2)

            splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
            splitter.addWidget(vbox1w)
            splitter.addWidget(vbox2w)

            self.enter_data = QtGui.QTextEdit()
            self.execute_python_button = QtGui.QPushButton("Execute")
            self.execute_and_clear_python_button = QtGui.QPushButton("Execute and Clear")
            vbox1.addWidget(self.log_widget)
            vbox2.addWidget(QtGui.QLabel('Python Console:'))
            vbox2.addWidget(self.enter_data)
            hbox.addWidget(self.execute_python_button)
            hbox.addWidget(self.execute_and_clear_python_button)
            vbox2.addLayout(hbox)

            self.connect(self.execute_python_button, QtCore.SIGNAL('clicked()'), self.on_execute_python_button)
            self.connect(self.execute_and_clear_python_button, QtCore.SIGNAL('clicked()'), self.on_execute_and_clear_python_button)
            self.dock_widget.setWidget(splitter)
        else:
            self.dock_widget.setWidget(self.log_widget)

    def on_execute_and_clear_python_button(self):
        self.parent._on_execute_python_button(clear=True)

    def on_execute_python_button(self):
        self.parent._on_execute_python_button(clear=False)

