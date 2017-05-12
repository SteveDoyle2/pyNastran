"""
defines:
 - LegendPropertiesWindow
"""
from __future__ import print_function
import os

from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4 import QtCore#, QtGui
    from PyQt4.QtGui import (
        QApplication, QLabel, QPushButton, QLineEdit, QComboBox, QWidget, QRadioButton,
        QButtonGroup, QGridLayout, QHBoxLayout, QVBoxLayout, QFont)
elif qt_version == 5:
    #from PyQt5 import QtCore, QtGui
    from PyQt5.QtGui import QFont
    from PyQt5.QtWidgets import (
        QApplication, QLabel, QPushButton, QLineEdit, QComboBox, QWidget, QRadioButton,
        QButtonGroup, QGridLayout, QHBoxLayout, QVBoxLayout)
elif qt_version == 'pyside':
    from PySide import QtCore#, QtGui
    from PySide.QtGui import (
        QApplication, QLabel, QPushButton, QLineEdit, QComboBox, QWidget, QRadioButton,
        QButtonGroup, QGridLayout, QHBoxLayout, QVBoxLayout, QFont)
else:
    raise NotImplementedError('qt_version = %r' % qt_version)

#from pyNastran.gui.qt_files.menu_utils import eval_float_from_string
from pyNastran.gui.colormaps import colormap_keys

from pyNastran.gui.gui_interface.common import PyDialog
from pyNastran.gui.gui_interface.legend.animation import AnimationWindow


class LegendPropertiesWindow(PyDialog):
    """
    +-------------------+
    | Legend Properties |
    +-----------------------+
    | Title  ______ Default |
    | Min    ______ Default |
    | Max    ______ Default |
    | Format ______ Default |
    | Scale  ______ Default |
    | Phase  ______ Default |
    | Number of Colors ____ |
    | Number of Labels ____ |
    | Label Size       ____ | (TODO)
    | ColorMap         ____ | (TODO)
    |                       |
    | x Min/Max (Blue->Red) |
    | o Max/Min (Red->Blue) |
    |                       |
    | x Vertical/Horizontal |
    | x Show/Hide           |
    |                       |
    |        Animate        |
    |    Apply OK Cancel    |
    +-----------------------+
    """

    def __init__(self, data, win_parent=None):
        PyDialog.__init__(self, data, win_parent)

        self._updated_legend = False
        self._animation_window_shown = False

        self._icase = data['icase']
        self._default_icase = self._icase

        self._default_name = data['name']
        self._default_min = data['min']
        self._default_max = data['max']

        self._default_scale = data['default_scale']
        self._scale = data['scale']

        self._default_phase = data['default_phase']
        self._phase = data['phase']

        self._default_format = data['default_format']
        self._format = data['format']

        self._default_labelsize = data['default_labelsize']
        self._labelsize = data['labelsize']

        self._default_nlabels = data['default_nlabels']
        self._nlabels = data['nlabels']

        self._default_ncolors = data['default_ncolors']
        self._ncolors = data['ncolors']

        self._default_colormap = data['default_colormap']
        self._colormap = data['colormap']

        self._default_is_low_to_high = data['is_low_to_high']

        self._default_is_discrete = data['is_discrete']
        self._default_is_horizontal = data['is_horizontal']
        self._default_is_shown = data['is_shown']
        self._is_normals = data['is_normals']

        self._update_defaults_to_blank()

        #self.setupUi(self)
        self.setWindowTitle('Legend Properties')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        self.set_font_size(data['font_size'])

    def _update_defaults_to_blank(self):
        """Changes the default (None) to a blank string"""
        if self._default_colormap is None:
            self._default_colormap = 'jet'
        if self._default_labelsize is None:
            self._default_labelsize = ''
        if self._default_ncolors is None:
            self._default_ncolors = ''
        if self._default_nlabels is None:
            self._default_nlabels = ''

        if self._colormap is None:
            self._colormap = 'jet'
        if self._labelsize is None:
            self._labelsize = ''
        if self._ncolors is None:
            self._ncolors = ''
        if self._nlabels is None:
            self._nlabels = ''

    def update_legend(self, icase, name,
                      min_value, max_value, data_format, scale, phase,
                      nlabels, labelsize,
                      ncolors, colormap,

                      default_title, default_min_value, default_max_value,
                      default_data_format, default_scale, default_phase,
                      default_nlabels, default_labelsize,
                      default_ncolors, default_colormap,
                      is_low_to_high, is_horizontal_scalar_bar, is_normals,
                      font_size=8):
        """
        We need to update the legend if there's been a result change request
        """
        self.set_font_size(font_size)
        if icase != self._default_icase:
            self._icase = icase
            self._default_icase = icase
            self._default_name = default_title
            self._default_min = default_min_value
            self._default_max = default_max_value
            self._default_format = default_data_format
            self._default_is_low_to_high = is_low_to_high
            self._default_is_discrete = True
            self._default_is_horizontal = is_horizontal_scalar_bar
            self._default_scale = default_scale
            self._default_phase = default_phase
            self._default_nlabels = default_nlabels
            self._default_labelsize = default_labelsize
            self._default_ncolors = default_ncolors
            self._default_colormap = default_colormap
            self._is_normals = is_normals

            if colormap is None:
                colormap = 'jet'
            if labelsize is None:
                labelsize = ''
            if ncolors is None:
                ncolors = ''
            if nlabels is None:
                nlabels = ''

            self._update_defaults_to_blank()

            assert isinstance(scale, float), 'scale=%r' % scale
            assert isinstance(default_scale, float), 'default_scale=%r' % default_scale
            if self._default_scale == 0.0:
                self.scale.setEnabled(False)
                self.scale_edit.setEnabled(False)
                self.scale_button.setEnabled(False)
            else:
                self.scale.setEnabled(True)
                self.scale_edit.setEnabled(True)
                self.scale_button.setEnabled(True)

            if self._default_phase is None:
                self._phase = None
                self.phase.setEnabled(False)
                self.phase_edit.setEnabled(False)
                self.phase_button.setEnabled(False)
                self.phase_edit.setText('0.0')
                self.phase_edit.setStyleSheet("QLineEdit{background: white;}")
            else:
                self._phase = phase
                self.phase.setEnabled(True)
                self.phase_edit.setEnabled(True)
                self.phase_button.setEnabled(True)
                self.phase_edit.setText(str(phase))
                self.phase_edit.setStyleSheet("QLineEdit{background: white;}")

            #self.on_default_name()
            #self.on_default_min()
            #self.on_default_max()
            #self.on_default_format()
            #self.on_default_scale()
            # reset defaults
            self._name = name
            self.name_edit.setText(name)
            self.name_edit.setStyleSheet("QLineEdit{background: white;}")

            self.min_edit.setText(str(min_value))
            self.min_edit.setStyleSheet("QLineEdit{background: white;}")

            self.max_edit.setText(str(max_value))
            self.max_edit.setStyleSheet("QLineEdit{background: white;}")

            self.format_edit.setText(str(data_format))
            self.format_edit.setStyleSheet("QLineEdit{background: white;}")

            self._scale = scale
            self.scale_edit.setText(str(scale))
            self.scale_edit.setStyleSheet("QLineEdit{background: white;}")

            self.nlabels_edit.setText(str(nlabels))
            self.nlabels_edit.setStyleSheet("QLineEdit{background: white;}")

            self.labelsize_edit.setText(str(labelsize))
            self.labelsize_edit.setStyleSheet("QLineEdit{background: white;}")

            self.ncolors_edit.setText(str(ncolors))
            self.ncolors_edit.setStyleSheet("QLineEdit{background: white;}")

            self.colormap_edit.setCurrentIndex(colormap_keys.index(str(colormap)))

            # lots of hacking for the Normal vectors
            enable = True
            if self._is_normals:
                enable = False

            self.max.setVisible(enable)
            self.min.setVisible(enable)
            self.max_edit.setVisible(enable)
            self.min_edit.setVisible(enable)
            self.max_button.setVisible(enable)
            self.min_button.setVisible(enable)

            self.show_radio.setVisible(enable)
            self.hide_radio.setVisible(enable)
            self.low_to_high_radio.setVisible(enable)
            self.high_to_low_radio.setVisible(enable)

            self.format.setVisible(enable)
            self.format_edit.setVisible(enable)
            self.format_edit.setVisible(enable)
            self.format_button.setVisible(enable)

            self.nlabels.setVisible(enable)
            self.nlabels_edit.setVisible(enable)
            self.nlabels_button.setVisible(enable)

            self.ncolors.setVisible(enable)
            self.ncolors_edit.setVisible(enable)
            self.ncolors_button.setVisible(enable)

            self.grid2_title.setVisible(enable)
            self.vertical_radio.setVisible(enable)
            self.horizontal_radio.setVisible(enable)

            self.colormap.setVisible(enable)
            self.colormap_edit.setVisible(enable)
            self.colormap_button.setVisible(enable)

            self.on_apply()

    def create_widgets(self):
        """creates the menu objects"""
        # Name
        self.name = QLabel("Title:")
        self.name_edit = QLineEdit(str(self._default_name))
        self.name_button = QPushButton("Default")

        # Min
        self.min = QLabel("Min:")
        self.min_edit = QLineEdit(str(self._default_min))
        self.min_button = QPushButton("Default")

        # Max
        self.max = QLabel("Max:")
        self.max_edit = QLineEdit(str(self._default_max))
        self.max_button = QPushButton("Default")

        #---------------------------------------
        # Format
        self.format = QLabel("Format (e.g. %.3f, %g, %.6e):")
        self.format_edit = QLineEdit(str(self._format))
        self.format_button = QPushButton("Default")

        #---------------------------------------
        # Scale
        self.scale = QLabel("Scale:")
        self.scale_edit = QLineEdit(str(self._scale))
        self.scale_button = QPushButton("Default")
        if self._default_scale == 0.0:
            self.scale.setVisible(False)
            self.scale_edit.setVisible(False)
            self.scale_button.setVisible(False)

        # Phase
        self.phase = QLabel("Phase (deg):")
        self.phase_edit = QLineEdit(str(self._phase))
        self.phase_button = QPushButton("Default")
        if self._default_phase is None:
            self.phase.setVisible(False)
            self.phase_edit.setVisible(False)
            self.phase_button.setVisible(False)
            self.phase_edit.setText('0.0')
        #tip = QtGui.QToolTip()
        #tip.setTe
        #self.format_edit.toolTip(tip)

        #---------------------------------------
        # nlabels
        self.nlabels = QLabel("Number of Labels:")
        self.nlabels_edit = QLineEdit(str(self._nlabels))
        self.nlabels_button = QPushButton("Default")

        self.labelsize = QLabel("Label Size:")
        self.labelsize_edit = QLineEdit(str(self._labelsize))
        self.labelsize_button = QPushButton("Default")

        self.ncolors = QLabel("Number of Colors:")
        self.ncolors_edit = QLineEdit(str(self._ncolors))
        self.ncolors_button = QPushButton("Default")

        self.colormap = QLabel("Color Map:")
        self.colormap_edit = QComboBox(self)
        self.colormap_button = QPushButton("Default")
        for key in colormap_keys:
            self.colormap_edit.addItem(key)
        self.colormap_edit.setCurrentIndex(colormap_keys.index(self._colormap))

        # --------------------------------------------------------------
        # the header
        self.grid2_title = QLabel("Color Scale:")

        # red/blue or blue/red
        self.low_to_high_radio = QRadioButton('Low -> High')
        self.high_to_low_radio = QRadioButton('High -> Low')
        widget = QWidget(self)
        low_to_high_group = QButtonGroup(widget)
        low_to_high_group.addButton(self.low_to_high_radio)
        low_to_high_group.addButton(self.high_to_low_radio)
        self.low_to_high_radio.setChecked(self._default_is_low_to_high)
        self.high_to_low_radio.setChecked(not self._default_is_low_to_high)

        # horizontal / vertical
        self.horizontal_radio = QRadioButton("Horizontal")
        self.vertical_radio = QRadioButton("Vertical")
        widget = QWidget(self)
        horizontal_vertical_group = QButtonGroup(widget)
        horizontal_vertical_group.addButton(self.horizontal_radio)
        horizontal_vertical_group.addButton(self.vertical_radio)
        self.horizontal_radio.setChecked(self._default_is_horizontal)
        self.vertical_radio.setChecked(not self._default_is_horizontal)

        # on / off
        self.show_radio = QRadioButton("Show")
        self.hide_radio = QRadioButton("Hide")
        widget = QWidget(self)
        show_hide_group = QButtonGroup(widget)
        show_hide_group.addButton(self.show_radio)
        show_hide_group.addButton(self.hide_radio)
        self.show_radio.setChecked(self._default_is_shown)
        self.hide_radio.setChecked(not self._default_is_shown)

        # --------------------------------------------------------------

        if self._is_normals:
            self.max.hide()
            self.min.hide()
            self.max_edit.hide()
            self.min_edit.hide()
            self.max_button.hide()
            self.min_button.hide()

            self.format.hide()
            self.format_edit.hide()
            self.format_button.hide()

            self.nlabels.hide()
            self.nlabels_edit.hide()
            self.nlabels_button.hide()

            self.ncolors.hide()
            self.ncolors_edit.hide()
            self.ncolors_button.hide()

            self.grid2_title.hide()
            self.vertical_radio.hide()
            self.horizontal_radio.hide()
            self.show_radio.hide()
            self.hide_radio.hide()
            self.low_to_high_radio.hide()
            self.high_to_low_radio.hide()

            self.colormap.hide()
            self.colormap_edit.hide()
            self.colormap_button.hide()

        self.animate_button = QPushButton('Create Animation')

        if self._default_scale == 0.0:
            self.animate_button.setEnabled(False)
            self.animate_button.setToolTip('This must be a displacement-like result to animate')

        # closing
        self.apply_button = QPushButton("Apply")
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")

    def create_layout(self):
        """displays the menu objects"""
        grid = QGridLayout()
        grid.addWidget(self.name, 0, 0)
        grid.addWidget(self.name_edit, 0, 1)
        grid.addWidget(self.name_button, 0, 2)

        grid.addWidget(self.min, 1, 0)
        grid.addWidget(self.min_edit, 1, 1)
        grid.addWidget(self.min_button, 1, 2)

        grid.addWidget(self.max, 2, 0)
        grid.addWidget(self.max_edit, 2, 1)
        grid.addWidget(self.max_button, 2, 2)

        grid.addWidget(self.format, 3, 0)
        grid.addWidget(self.format_edit, 3, 1)
        grid.addWidget(self.format_button, 3, 2)

        grid.addWidget(self.scale, 4, 0)
        grid.addWidget(self.scale_edit, 4, 1)
        grid.addWidget(self.scale_button, 4, 2)

        grid.addWidget(self.phase, 5, 0)
        grid.addWidget(self.phase_edit, 5, 1)
        grid.addWidget(self.phase_button, 5, 2)

        grid.addWidget(self.nlabels, 6, 0)
        grid.addWidget(self.nlabels_edit, 6, 1)
        grid.addWidget(self.nlabels_button, 6, 2)

        #grid.addWidget(self.labelsize, 6, 0)
        #grid.addWidget(self.labelsize_edit, 6, 1)
        #grid.addWidget(self.labelsize_button, 6, 2)

        grid.addWidget(self.ncolors, 7, 0)
        grid.addWidget(self.ncolors_edit, 7, 1)
        grid.addWidget(self.ncolors_button, 7, 2)

        grid.addWidget(self.colormap, 8, 0)
        grid.addWidget(self.colormap_edit, 8, 1)
        grid.addWidget(self.colormap_button, 8, 2)

        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)


        grid2 = QGridLayout()
        grid2.addWidget(self.grid2_title, 0, 0)
        grid2.addWidget(self.low_to_high_radio, 1, 0)
        grid2.addWidget(self.high_to_low_radio, 2, 0)

        grid2.addWidget(self.vertical_radio, 1, 1)
        grid2.addWidget(self.horizontal_radio, 2, 1)

        grid2.addWidget(self.show_radio, 1, 2)
        grid2.addWidget(self.hide_radio, 2, 2)

        grid2.addWidget(self.animate_button, 3, 1)


        #grid2.setSpacing(0)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        #vbox.addLayout(checkboxes)
        vbox.addLayout(grid2)
        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)

        #Create central widget, add layout and set
        #central_widget = QtGui.QWidget()
        #central_widget.setLayout(vbox)
        #self.setCentralWidget(central_widget)
        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the buttons"""
        self.name_button.clicked.connect(self.on_default_name)
        self.min_button.clicked.connect(self.on_default_min)
        self.max_button.clicked.connect(self.on_default_max)
        self.format_button.clicked.connect(self.on_default_format)
        self.scale_button.clicked.connect(self.on_default_scale)
        self.phase_button.clicked.connect(self.on_default_phase)

        self.nlabels_button.clicked.connect(self.on_default_nlabels)
        self.labelsize_button.clicked.connect(self.on_default_labelsize)
        self.ncolors_button.clicked.connect(self.on_default_ncolors)
        self.colormap_button.clicked.connect(self.on_default_colormap)

        self.animate_button.clicked.connect(self.on_animate)

        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)

        if qt_version == 4:
            self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)
            #self.colormap_edit.activated[str].connect(self.onActivated)
        #else:
            # closeEvent???

    def set_font_size(self, font_size):
        """
        Updates the font size of the objects

        Parameters
        ----------
        font_size : int
            the font size
        """
        if self.font_size == font_size:
            return
        self.font_size = font_size
        font = QFont()
        font.setPointSize(font_size)
        self.setFont(font)
        self.name_edit.setFont(font)
        self.min_edit.setFont(font)
        self.max_edit.setFont(font)
        self.format_edit.setFont(font)
        self.scale_edit.setFont(font)
        self.phase_edit.setFont(font)
        self.nlabels_edit.setFont(font)
        self.labelsize_edit.setFont(font)
        self.ncolors_edit.setFont(font)

    def on_animate(self):
        name, flag0 = self.check_name(self.name_edit)
        if not flag0:
            return

        data = {
            'font_size' : self.out_data['font_size'],
            'icase' : self._icase,
            'name' : name,
            'time' : 2,
            'frames/sec' : 30,
            'resolution' : 1,
            'iframe' : 0,
            'scale' : self._scale,
            'default_scale' : self._default_scale,

            'is_scale' : self._default_phase is None,
            'phase' : self._phase,
            'default_phase' : self._default_phase,
            'dirname' : os.path.abspath(os.getcwd()),
            'clicked_ok' : False,
            'close' : False,
        }
        if not self._animation_window_shown:
            self._animation_window = AnimationWindow(data, win_parent=self)
            self._animation_window.show()
            self._animation_window_shown = True
            self._animation_window.exec_()
        else:
            self._animation_window.activateWindow()

        if data['close']:
            if not self._animation_window._updated_animation:
                #self._apply_animation(data)
                pass
            self._animation_window_shown = False
            del self._animation_window
        else:
            self._animation_window.activateWindow()

    def on_default_name(self):
        """action when user clicks 'Default' for name"""
        name = str(self._default_name)
        self.name_edit.setText(name)
        self.name_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_min(self):
        """action when user clicks 'Default' for min value"""
        self.min_edit.setText(str(self._default_min))
        self.min_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_max(self):
        """action when user clicks 'Default' for max value"""
        self.max_edit.setText(str(self._default_max))
        self.max_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_format(self):
        """action when user clicks 'Default' for the number format"""
        self.format_edit.setText(str(self._default_format))
        self.format_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_scale(self):
        """action when user clicks 'Default' for scale factor"""
        self.scale_edit.setText(str(self._default_scale))
        self.scale_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_phase(self):
        """action when user clicks 'Default' for phase angle"""
        self.phase_edit.setText(str(self._default_phase))
        self.phase_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_ncolors(self):
        """action when user clicks 'Default' for number of colors"""
        self.ncolors_edit.setText(str(self._default_ncolors))
        self.ncolors_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_colormap(self):
        """action when user clicks 'Default' for the color map"""
        self.colormap_edit.setCurrentIndex(colormap_keys.index(self._default_colormap))

    def on_default_nlabels(self):
        """action when user clicks 'Default' for number of labels"""
        self.nlabels_edit.setStyleSheet("QLineEdit{background: white;}")
        self.nlabels_edit.setText(str(self._default_nlabels))

    def on_default_labelsize(self):
        """action when user clicks 'Default' for number of labelsize"""
        self.labelsize_edit.setText(str(self._default_labelsize))
        self.labelsize_edit.setStyleSheet("QLineEdit{background: white;}")

    @staticmethod
    def check_name(cell):
        cell_value = cell.text()
        try:
            text = str(cell_value).strip()
        except UnicodeEncodeError:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

        if len(text):
            cell.setStyleSheet("QLineEdit{background: white;}")
            return text, True
        else:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    @staticmethod
    def check_colormap(cell):
        text = str(cell.text()).strip()
        if text in colormap_keys:
            cell.setStyleSheet("QLineEdit{background: white;}")
            return text, True
        else:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    def on_validate(self):
        name_value, flag0 = self.check_name(self.name_edit)
        min_value, flag1 = self.check_float(self.min_edit)
        max_value, flag2 = self.check_float(self.max_edit)
        format_value, flag3 = self.check_format(self.format_edit)
        scale, flag4 = self.check_float(self.scale_edit)
        phase, flag5 = self.check_float(self.phase_edit)

        nlabels, flag6 = self.check_positive_int_or_blank(self.nlabels_edit)
        ncolors, flag7 = self.check_positive_int_or_blank(self.ncolors_edit)
        labelsize, flag8 = self.check_positive_int_or_blank(self.labelsize_edit)
        colormap = str(self.colormap_edit.currentText())

        if all([flag0, flag1, flag2, flag3, flag4, flag5, flag6, flag7, flag8]):
            if 'i' in format_value:
                format_value = '%i'

            assert isinstance(scale, float), scale
            self.out_data['name'] = name_value
            self.out_data['min'] = min_value
            self.out_data['max'] = max_value
            self.out_data['format'] = format_value
            self.out_data['scale'] = scale
            self.out_data['phase'] = phase

            self.out_data['nlabels'] = nlabels
            self.out_data['ncolors'] = ncolors
            self.out_data['labelsize'] = labelsize
            self.out_data['colormap'] = colormap

            self.out_data['is_low_to_high'] = self.low_to_high_radio.isChecked()
            self.out_data['is_horizontal'] = self.horizontal_radio.isChecked()
            self.out_data['is_shown'] = self.show_radio.isChecked()

            self.out_data['clicked_ok'] = True
            self.out_data['close'] = True
            #print('self.out_data = ', self.out_data)
            #print("name = %r" % self.name_edit.text())
            #print("min = %r" % self.min_edit.text())
            #print("max = %r" % self.max_edit.text())
            #print("format = %r" % self.format_edit.text())
            return True
        return False

    def on_apply(self):
        passed = self.on_validate()
        if passed:
            self.win_parent._apply_legend(self.out_data)
        return passed

    def on_ok(self):
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        self.out_data['close'] = True
        self.close()


def main(): # pragma: no cover
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)


    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    #The Main window
    data1 = {
        'icase' : 1,
        'name' : 'asdf',
        'min' : 0.,
        'max' : 10,
        'scale' : 1e-12,
        'default_scale' : 1.0,

        'phase' : 0.0,
        #'default_phase' : 180.0,
        'default_phase' : None,

        'nlabels' : 11,
        'default_nlabels' : 11,

        'labelsize' : 12,
        'default_labelsize' : 12,

        'ncolors' : 13,
        'default_ncolors' : 13,

        'colormap' : 'jet',
        'default_colormap' : 'jet',

        'default_format' : '%s',
        'format' : '%g',

        'is_low_to_high': True,
        'is_discrete' : False,
        'is_horizontal' : False,
        'is_shown' : True,
    }
    main_window = LegendPropertiesWindow(data1)

    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__": # pragma: no cover
    main()
