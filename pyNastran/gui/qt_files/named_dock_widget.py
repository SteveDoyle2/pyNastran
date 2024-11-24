from qtpy.QtWidgets import QWidget, QDockWidget


class NamedDockWidget(QDockWidget):
    def __init__(self, name: str, widget: QWidget, parent=None):
        QDockWidget.__init__(self, name, parent=parent)
        self.setObjectName(name)
        self.widgeti = widget
        self.setWidget(widget)
        #self.setFont(QFont('Arial', FONT_SIZE))
