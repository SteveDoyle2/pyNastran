from PyQt4 import QtCore, QtGui


class PythonConsoleWidget(QtGui.QDockWidget):
    def __init__(self, parent):
        self.parent = parent
        super(QtGui.QDockWidget, self).__init__('Python Console', parent=parent)
        vbox = QtGui.QVBoxLayout()
        hbox = QtGui.QHBoxLayout()

        self.enter_data = QtGui.QTextEdit()
        self.execute_python_button = QtGui.QPushButton("Execute")
        self.execute_and_clear_python_button = QtGui.QPushButton("Execute and Clear")
        vbox.addWidget(self.enter_data)
        hbox.addWidget(self.execute_python_button)
        hbox.addWidget(self.execute_and_clear_python_button)
        vbox.addLayout(hbox)


        self.connect(self.execute_python_button, QtCore.SIGNAL('clicked()'), self.on_execute_python_button)
        self.connect(self.execute_and_clear_python_button, QtCore.SIGNAL('clicked()'), self.on_execute_and_clear_python_button)

        # we have to create a QWidget to put the console vbox into vbox_widget because
        #    self.setLayout(vbox)
        # is not supported in a QDockWidget
        vbox_widget = QtGui.QWidget()
        vbox_widget.setLayout(vbox)
        self.setWidget(vbox_widget)

    def on_execute_and_clear_python_button(self):
        self.parent._on_execute_python_button(clear=True)

    def on_execute_python_button(self):
        self.parent._on_execute_python_button(clear=False)
