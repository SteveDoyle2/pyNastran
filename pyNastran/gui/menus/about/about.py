"""
The preferences menu handles:
 - Font Size
 - Background Color
 - Text Color
 - Annotation Color
 - Annotation Size
 - Clipping Min
 - Clipping Max

"""
import os
import sys
import locale

from qtpy.QtCore import Qt
from qtpy import QtGui
from qtpy.QtWidgets import (
    QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
    QTabWidget, QWidget, QScrollArea, QTextEdit)

from pyNastran.gui.utils.qt.pydialog import PyDialog
#from pyNastran.gui.utils.qt.qpush_button_color import QPushButtonColor
#from pyNastran.gui.utils.qt.dialogs import save_file_dialog

CREDITS = """
pyNastran was written by Steve Doyle since 2011.  This product contains the following third party modules:

  * Python, the programming language, written by Guido van Rossum and many contributors.

  * Qt5 cross-platform GUI toolkit, developed by many contributors.

  * PyQt5 Python bindings for Qt5, by Riverbank Computing Limited.

  * Python Imaging Library, developed by Secret Labs AB and Fredrik Lundh.

  * Scintilla, a source code editor widget, written by Neil Hodgson and many contributors.

  * Docutils, tools for ReST document conversion, written by David Goodger and contributors.

  * Pygments by Georg Brandl, Armin Ronacher, Tim Hatch, and contributors.

We gratefully acknowledge the efforts of all that have contributed to these and the other open source products and tools that are used in the development of pyNastran.
""".replace('\n', '<br>')

class AboutWindow(PyDialog):
    """
    +-------------+
    | AboutWindow |
    +------------------------+
    | Origin/P1   cid  x y z |
    | P2          cid  x y z |
    | z-axis      cid  x y z |
    | tol         cid  x y z |
    |                        |
    |    Apply OK Cancel     |
    +------------------------+
    """
    def __init__(self, data, win_parent=None, show_tol=True):
        """
        Saves the data members from data and
        performs type checks
        """
        PyDialog.__init__(self, data, win_parent)

        self._updated_preference = False

        self._default_font_size = data['font_size']
        #self.out_data = data

        self.setWindowTitle('About pyNastran GUI')
        self.create_widgets(show_tol)
        self.create_layout()
        #self.set_connections()
        self.on_font(self._default_font_size)
        #self.show()

    def create_widgets(self, show_tol):
        """creates the display window"""
        # CORD2R
        #self.origin_label = QLabel("Origin:")
        #self.zaxis_label = QLabel("Z Axis:")
        #self.xz_plane_label = QLabel("XZ Plane:")

        #-----------------------------------------------------------------------
        # closing
        self.apply_button = QPushButton('Apply')
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Cancel')

    def create_layout(self):
        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        #---------------------
        version_tab, len_version = _version_tab(ok_cancel_box)
        package_tab = _package_tab(len_version)
        credits_tab = _credits_tab()
        # --------------------
        tab_widget = QTabWidget()
        tab_widget.addTab(version_tab, 'Version')
        tab_widget.addTab(package_tab, 'Packages')
        tab_widget.addTab(credits_tab, 'Credits')

        #---------------------
        vbox_outer = QVBoxLayout()
        vbox_outer.addWidget(tab_widget)
        #---------------------

        self.setLayout(vbox_outer)
        #hint = vbox.sizeHint()
        #print(hint)

        # PySide2.QtCore.QSize(516, 212)
        #hint.setHeight(hint.height() * 1.3)
        #hint.setWidth(hint.width() * 1.1)
        #self.setFixedSize(hint)

    #def set_connections(self):
        #"""creates the actions for the menu"""
        #self.method_pulldown.currentIndexChanged.connect(self.on_method)
        #self.zaxis_method_pulldown.currentIndexChanged.connect(self.on_zaxis_method)
        #self.plane_color_edit.clicked.connect(self.on_plane_color)

        #self.apply_button.clicked.connect(self.on_apply)
        #self.ok_button.clicked.connect(self.on_ok)
        #self.cancel_button.clicked.connect(self.on_cancel)
        ## closeEvent
        #return

    def on_font(self, value=None):
        """update the font for the current window"""
        if value is None:
            value = self.font_size_edit.value()
        font = QtGui.QFont()
        font.setPointSize(value)
        self.setFont(font)

    def on_ok(self):
        #passed = self.on_apply()
        #if passed:
        self.close()
        #self.destroy()

    def on_cancel(self):
        self.out_data['close'] = True
        self.close()

def get_packages(len_version=80):
    import sys
    import numpy
    import scipy
    #import matplotlib
    #import pandas
    import vtk
    from pyNastran.gui.qt_version import qt_version


    if qt_version == 'pyqt5':
        import PyQt5
        qt_name = 'PyQt5'
        _qt_version = PyQt5.__version__
    elif qt_version == 'pyside2':
        import PySide2
        qt_name = 'PySide2'
        _qt_version = PySide2.__version__
    else:
        raise NotImplementedError(qt_version)

    import importlib

    python = str(sys.version_info)
    packages = {
        'Python' : python + ' ' * (len_version - len(python) + 15),
        'numpy' : numpy.__version__,
        'scipy' : scipy.__version__,
        #'matplotlib' : matplotlib.__version__,
        #'pandas' : pandas.__version__,
        'matplotlib' : 'N/A',
        'pandas' : 'N/A',
        'vtk' : vtk.VTK_VERSION,
        #'PyQt5':,
        qt_name : _qt_version,
     }
    for name in ['matplotlib', 'pandas', 'docopt']:
        module = importlib.import_module(name, package=None)
        packages[name] = module.__version__
    return packages

def _version_tab(ok_cancel_box):
    import pyNastran
    platform = sys.platform
    localei, unused_encoding = locale.getdefaultlocale()
    try:
        os_version = str(sys.getwindowsversion())
    except:
        os_version = '???'

    import platform
    cpu = platform.processor()

    version_data = {
        'Product': 'pyNastran GUI',
        'Version': pyNastran.__version__,
        'Release Type': 'Final Release',
        'Release Date': pyNastran.__releaseDate__,
        #'Cache Directory': ,
        'OS' : f'win32 (sys.platform={platform})',
        'OS Version' : os_version,
        'CPU': cpu,
        'Memory': str(['1000 bytes']),
        'Locale': localei,
    }
    len_version = len(os_version)
    grid = grid_from_dict(version_data)

    hbox = QHBoxLayout()
    hbox.addLayout(grid)
    hbox.addStretch()

    vbox = QVBoxLayout()
    vbox.addLayout(hbox)
    vbox.addStretch()
    vbox.addLayout(ok_cancel_box)

    #---------------------
    version_tab = QWidget()
    version_tab.setLayout(vbox)

    return version_tab, len_version

def _package_tab(len_version=80):
    """makes the packages tab"""
    packages = get_packages(len_version=len_version)
    grid = grid_from_dict(packages)

    vbox = QVBoxLayout()
    vbox.addLayout(grid)
    vbox.addStretch()

    package_tab = QWidget()
    package_tab.setLayout(vbox)
    return package_tab

def grid_from_dict(mydict):
    irow = 0
    grid = QGridLayout()
    for key, valuei in mydict.items():
        label = QLabel(key + ':')
        label.setAlignment(Qt.AlignRight)

        value = QLabel(valuei)
        value.setTextInteractionFlags(Qt.TextSelectableByMouse)
        grid.addWidget(label, irow, 0)
        grid.addWidget(value, irow, 1)
        irow += 1
    return grid

#class Window(QScrollArea):
    #def __init__(self):
        #super(Window, self).__init__()
        #widget = QWidget()
        #layout = QVBoxLayout(widget)
        #layout.setAlignment(Qt.AlignTop)
        #for index in range(100):
            #layout.addWidget(QLabel('Label %02d' % index))
        #self.setWidget(widget)
        #self.setWidgetResizable(True)

def _credits_tab():
    #scroll = QScrollArea()
    #scroll.setWidget(self)
    #scroll.setWidgetResizable(True)
    ##scroll.setFixedHeight(400)
    #layout.addWidget(scroll)

    #vbox = QVBoxLayout()
    #vbox.addLayout(layout)
    #vbox.addStretch()

    scrollArea = QScrollArea()
    scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    scrollArea.setWidgetResizable(True)
    #scrollArea->setGeometry( 10, 10, 200, 200 );

    package_tab = QWidget()
    scrollArea.setWidget(package_tab)

    widget = QWidget(scrollArea)

    vbox = QVBoxLayout(widget)
    text = QTextEdit(CREDITS)
    text.setReadOnly(True)
    vbox.addWidget(text)
    #vbox.addLayout(scrollArea)

    package_tab = QWidget()
    package_tab.setLayout(vbox)
    return package_tab

def main():
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)


    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    #The Main window
    data = {
        'font_size' : 8,
        #'cids' : [0, 1, 2, 3],
        'name' : 'main',

    }
    main_window = AboutWindow(data, show_tol=True)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__":
    main()
