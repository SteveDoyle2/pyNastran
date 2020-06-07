"""The about menu credits 3rd party packages."""
import os
import sys
import locale
import platform
import importlib
import warnings
from typing import Tuple, Dict
with warnings.catch_warnings():  # avoid an imp module deprecation warning
    warnings.simplefilter("ignore")
    import setuptools

import numpy
import scipy
import vtk
import docopt
import pyNastran

import qtpy
from qtpy import QtGui
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QLabel, QPushButton, QGridLayout, QApplication, QHBoxLayout, QVBoxLayout,
    QTabWidget, QWidget, QScrollArea, QTextEdit, QMessageBox, QFrame, QGroupBox)

import cpylog
from pyNastran.gui import ICON_PATH, IS_LINUX, IS_MAC # IS_WINDOWS
from pyNastran.gui.qt_version import qt_name, PYQT_VERSION, is_pygments
from pyNastran.gui.menus.python_console import QSCINTILLA_VERSION
from pyNastran.gui.utils.qt.pydialog import PyDialog

QT = ("""
  * PyQt5 Python bindings for Qt5, by Riverbank Computing Limited.

  * Scintilla, a source code editor widget, written by Neil Hodgson and many contributors."""
    if qt_name == 'PyQt5' else """
  * PySide2 Python bindings for Qt5, by Qt for Python.""")

PYGMENTS = """
  * Pygments by Georg Brandl, Armin Ronacher, Tim Hatch, and contributors.

  * Python Imaging Library, developed by Secret Labs AB and Fredrik Lundh.
""" if is_pygments else ''

CREDITS = f"""pyNastran has been written by Steve Doyle since 2011.  This product contains the following third party modules:

  * Numpy array library, developed by many contributors.

  * Scipy scientific library, developed by many contributors.

  * Python, the programming language, written by Guido van Rossum and many contributors.

  * VTK Python bindings for Qt5, by Riverbank Computing Limited.

  * Qt5 cross-platform GUI toolkit, developed by many contributors.
{QT}
{PYGMENTS}
  * ImageIO, an animation library for writing videos, developed by many contributors.

  * WingIDE, the primary IDE used for development, by Wingware.

I gratefully acknowledge the efforts of all that have contributed to these and the other open source products and tools that are used in the development of pyNastran.
""".replace('\n', '<br>')


class AboutWindow(PyDialog):
    """
    +-------------+
    | AboutWindow |
    +-------------+
    """
    def __init__(self, data, win_parent=None, show_tol=True):
        """
        Saves the data members from data and
        performs type checks
        """
        PyDialog.__init__(self, data, win_parent)

        self._default_font_size = data['font_size']

        self.setWindowTitle('About pyNastran GUI')
        self.create_widgets(show_tol)
        self.create_layout()
        self.set_connections()
        self.on_font(self._default_font_size)

    def create_widgets(self, show_tol):
        """creates the display window"""
        #-----------------------------------------------------------------------
        # closing
        self.update_button = QPushButton('Check for Updates')
        self.ok_button = QPushButton('OK')
        #self.cancel_button = QPushButton('Cancel')

    def create_layout(self):
        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.update_button)
        ok_cancel_box.addStretch()
        ok_cancel_box.addWidget(self.ok_button)
        #ok_cancel_box.addWidget(self.cancel_button)

        #---------------------
        version_tab = _version_tab(ok_cancel_box)
        package_tab = _package_tab()
        shortcuts_tab = _shortcuts_tab()
        credits_tab = _credits_tab()
        # --------------------
        tab_widget = QTabWidget()
        tab_widget.addTab(version_tab, 'Version')
        tab_widget.addTab(package_tab, 'Packages')
        tab_widget.addTab(shortcuts_tab, 'Shortcuts')
        tab_widget.addTab(credits_tab, 'Credits')

        #---------------------
        png_filename = os.path.join(ICON_PATH, 'logo2.png')
        im = QtGui.QPixmap(png_filename)
        im_resized = im.scaled(525, 405, Qt.KeepAspectRatio)
        self.website_button = QLabel()
        self.website_button.setFrameShape(QFrame.HLine)
        self.website_button.setFrameStyle(QFrame.NoFrame)
        self.website_button.setPixmap(im_resized)
        #self.website_button.setIcon(im_resized)

        vbox_outer = QVBoxLayout()
        vbox_outer.setContentsMargins(0, 0, 0, 0)
        vbox_outer.addWidget(self.website_button)
        vbox_outer.addWidget(tab_widget)
        vbox_outer.addLayout(ok_cancel_box)
        #---------------------

        self.setLayout(vbox_outer)
        #hint = vbox.sizeHint()
        #print(hint)

        # PySide2.QtCore.QSize(516, 212)
        #hint.setHeight(hint.height() * 1.3)
        #hint.setWidth(hint.width() * 1.1)
        #self.setFixedSize(hint)

    def set_connections(self):
        #"""creates the actions for the menu"""
        #self.method_pulldown.currentIndexChanged.connect(self.on_method)
        #self.zaxis_method_pulldown.currentIndexChanged.connect(self.on_zaxis_method)
        #self.plane_color_edit.clicked.connect(self.on_plane_color)

        #self.apply_button.clicked.connect(self.on_apply)
        self.update_button.clicked.connect(self.on_update)
        self.ok_button.clicked.connect(self.on_ok)
        self.website_button.mousePressEvent = self.on_website
        #self.cancel_button.clicked.connect(self.on_cancel)

    def on_font(self, value=None):
        """update the font for the current window"""
        if value is None:
            value = self.font_size_edit.value()
        font = QtGui.QFont()
        font.setPointSize(value)
        self.setFont(font)

    def on_update(self):
        """check for a newer version"""
        is_newer = False
        if self.win_parent is not None:
            is_newer = self.win_parent._check_for_latest_version()
        if not is_newer:
            self.update_button.setDisabled(True)
            QMessageBox.about(self, 'About pyNastran GUI', 'PyNastran GUI is already up to date')

    def on_website(self, event):
        """opens the website"""
        if self.win_parent is None:
            return
        self.win_parent.open_website()

    def on_ok(self):
        """closes the window"""
        #passed = self.on_apply()
        #if passed:
        self.close()
        #self.destroy()

    def on_cancel(self):
        self.out_data['close'] = True
        self.close()

def get_packages() -> Dict[str, str]:
    """makes the packages data"""
    #python = str(sys.version_info)
    #'python_branch', 'python_revision', 'python_build', 'python_compiler', 'python_implementation',
    packages = {
        #'Python' : python + ' ' * (len_version - len(python) + 10),
        'Python': platform.python_branch(),
        #'Python revision': platform.python_revision(),
        #'Python Build': str(platform.python_build()),
        'Compiler': platform.python_compiler(),
        'Implementation': platform.python_implementation(),
        'setuptools': setuptools.__version__,
        'numpy' : numpy.__version__,
        'scipy' : scipy.__version__,
        'cpylog' : cpylog.__version__,
        'matplotlib' : 'N/A',
        'pandas' : 'N/A',
        'imageio' : 'N/A',
        'PIL' : 'N/A',
        'vtk' : vtk.VTK_VERSION,
        'qtpy' : qtpy.__version__,
        qt_name : PYQT_VERSION,
        'QScintilla2': QSCINTILLA_VERSION,
        'pygments' : 'N/A',
        'docopt-ng' : docopt.__version__,
    }
    if 'pyside' in qt_name.lower():
        del packages['QScintilla2']

    for name in ['matplotlib', 'pandas', 'imageio', 'PIL', 'pygments']:
        try:
            module = importlib.import_module(name, package=None)
        except ImportError:
            continue
        packages[name] = module.__version__
    return packages

def get_version() -> Dict[str, str]:
    """makes the version data"""
    sys_platform = sys.platform
    localei, unused_encoding = locale.getdefaultlocale()
    #try:
        #os_version = str(sys.getwindowsversion())
    #except:
        #os_version = '???'

    pmsg = [
        # 'win32_ver',
        # 'uname',
        'mac_ver', 'libc_ver',
    ]
    #if not IS_WINDOWS:
        #pmsg.remove('win32_ver')
    if not IS_LINUX:
        pmsg.remove('libc_ver')
    if not IS_MAC:
        pmsg.remove('mac_ver')

    #memory = str(sys.getsizeof(None))
    version_data = {
        'Product': 'pyNastran GUI',
        'Version': pyNastran.__version__,
        'Release Type': 'Final Release' if 'dev' not in pyNastran.__version__ else 'Developement',
        'Release Date': pyNastran.__releaseDate2__.title(),
        #'Cache Directory': ,
        'OS' : f'win32 (sys.platform={sys_platform})',
        'Platform' : platform.platform(),
        'Architecture' : str(platform.architecture()),
        'Machine' : platform.machine(),
        'Processor' : platform.processor(),

        #'Bit': bit,
        #'Memory': memory,
        'Locale': localei,
    }
    for key in pmsg:
        value = getattr(platform, key)()
        version_data[key.title()] = str(value)
    return version_data

def _version_tab(ok_cancel_box):
    """makes the version tab"""
    version_data = get_version()

    grid = grid_from_dict(version_data)

    hbox = layout_to_hlayout(grid)
    vbox = layout_to_vlayout(hbox)

    #---------------------
    version_tab = QWidget()
    version_tab.setLayout(vbox)

    return version_tab

def _package_tab():
    """makes the packages tab"""
    packages = get_packages()
    grid = grid_from_dict(packages)

    hbox = layout_to_hlayout(grid)
    vbox = layout_to_vlayout(hbox)

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

def _credits_tab():
    """creates the credits tab"""
    scroll_area = QScrollArea()
    scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    scroll_area.setWidgetResizable(True)
    scroll_widget = QWidget(scroll_area)

    widget = QWidget()
    scroll_area.setWidget(widget)

    vbox = QVBoxLayout(scroll_widget)
    text = QTextEdit(CREDITS)
    text.setReadOnly(True)
    vbox.addWidget(text)

    package_tab = QWidget()
    package_tab.setLayout(vbox)
    return package_tab

def _shortcuts_tab():
    """creates the credits tab"""
    #scroll_area = QScrollArea()
    #scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    #scroll_area.setWidgetResizable(True)
    #scroll_widget = QWidget(scroll_area)

    #widget = QWidget()
    #scroll_area.setWidget(widget)

    #vbox = QVBoxLayout(scroll_widget)
    vbox = QVBoxLayout()

    mouse_shortcuts, keyboard_shortcuts = get_shortcuts()
    mouse_grid = grid_from_dict(mouse_shortcuts)
    keyboard_grid = grid_from_dict(keyboard_shortcuts)

    mouse_group = QGroupBox('Mouse')
    mouse_group.setLayout(mouse_grid)
    mouse_layout = QVBoxLayout()
    mouse_layout.addWidget(mouse_group)
    mouse_layout2 = layout_to_hlayout(mouse_layout)

    keyboard_group = QGroupBox('Keyboard')
    keyboard_group.setLayout(keyboard_grid)
    keyboard_layout = QVBoxLayout()
    keyboard_layout.addWidget(keyboard_group)
    keyboard_layout2 = layout_to_hlayout(keyboard_layout)

    #keyboard_vbox = QGroupBox('Keyboard')
    #keyboard_vbox.addLayout(keyboard_grid)
    #text = QTextEdit(CREDITS)
    #text.setReadOnly(True)
    #vbox.addWidget(text)
    vbox.addLayout(mouse_layout2)
    vbox.addLayout(keyboard_layout2)
    vbox.addStretch()

    shortcuts_tab = QWidget()
    shortcuts_tab.setLayout(vbox)
    return shortcuts_tab

def get_shortcuts() -> Tuple[Dict[str, str], Dict[str, str]]:
    """makes the shortcuts data"""

    mouse_shortcuts = {
        'Left Click' : 'Rotate',
        'Middle Click' : 'Pan/Recenter Rotation Point',
        'Shift + Left Click' : 'Pan/Recenter Rotation Point',
        'Right Mouse / Wheel' : 'Zoom',
    }
    keyboard_shortcuts = {
        'R' : 'reset camera view',
        'Shift+X/X' : 'snap to x axis',
        'Shift+Y/Y' : 'snap to y axis',
        'Shift+Z/Z' : 'snap to z axis',
        #'',
        #'h   - show/hide legend & info',

        # shown on the menu
        #'CTRL+I - take a screenshot (image)',
        #'CTRL+W - clear the labels',
        #'CTRL+L - Legend',
        #'CTRL+A - Animation',
        #'S      - view model as a surface',
        #'W      - view model as a wireframe',

        'L' : 'cycle the results forwards',
        'K' : 'cycle the results backwards',
        'm/Shift+M' : 'scale up/scale down by 1.1 times',
        'o/Shift+O' : 'rotate counter-clockwise/clockwise 5 degrees',
        'P' : 'pick node/element',
        'F' : 'set rotation center (zoom out when picking to disable clipping)',
        'E' : 'view model edges',
        'B' : 'change edge color from scalar/black',
        #'',
        #'Reload Model:  using the same filename, reload the model',
    }
    return mouse_shortcuts, keyboard_shortcuts

def layout_to_vlayout(layout):
    vbox = QVBoxLayout()
    vbox.addLayout(layout)
    vbox.addStretch()
    return vbox

def layout_to_hlayout(layout):
    hbox = QHBoxLayout()
    hbox.addLayout(layout)
    hbox.addStretch()
    return hbox

def main():  # pragma: no cover
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
    }
    main_window = AboutWindow(data, show_tol=True)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == '__main__':   # pragma: no cover
    main()
