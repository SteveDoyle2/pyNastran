"""
defines:
 - LegendPropertiesWindow

"""
import os

from qtpy import QtCore
from qtpy.QtGui import QFont
from qtpy.QtWidgets import (
    QApplication, QLabel, QPushButton, QLineEdit, QComboBox, QWidget, QRadioButton,
    QButtonGroup, QGridLayout, QHBoxLayout, QVBoxLayout)

from pyNastran.utils.numpy_utils import float_types
from pyNastran.gui.utils.colormaps import colormap_keys
from pyNastran.gui.utils.qt.pydialog import (
    PyDialog, check_float, check_format, check_name_str,
    check_positive_int_or_blank)
from pyNastran.gui.qt_version import qt_int as qt_version

ANIMATE_TOOLTIP_OFF = 'This must be a displacement-like result to animate'
ANIMATE_TOOLTIP_ON = 'Creates an scale/phase/time animation'

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

    def __init__(self, data, win_parent=None, show_animation_button=True):
        PyDialog.__init__(self, data, win_parent)
        self.is_gui = win_parent is not None

        self._updated_legend = False
        self.external_call = True

        #if win_parent is None:
        self._icase_fringe = data['icase_fringe']
        self._icase_disp = data['icase_disp']
        self._icase_vector = data['icase_vector']
        #else:
            #self._icase_fringe = data['icase']
            #self._icase_disp = data['icase']
            #self._icase_vector = data['icase']

        self._default_icase_fringe = self._icase_fringe
        self._default_icase_disp = self._icase_disp
        self._default_icase_vector = self._icase_vector
        #print('*icase_fringe=%s icase_disp=%s icase_vector=%s' % (
            #self._default_icase_fringe, self._default_icase_disp, self._default_icase_vector))

        self._default_title = data['title']
        self._default_min = data['min_value']
        self._default_max = data['max_value']

        self._default_scale = data['default_scale']
        self._scale = data['scale']

        #if win_parent is None:
        self._default_arrow_scale = data['default_arrow_scale']
        self._arrow_scale = data['arrow_scale']
        #else:
            #self._default_arrow_scale = data['default_scale']
            #self._arrow_scale = data['scale']

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
        self._is_fringe = data['is_fringe']

        self._update_defaults_to_blank()

        self.setWindowTitle('Legend Properties')
        self.create_widgets(show_animation_button=show_animation_button)
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

    def update_legend(self, icase_fringe, icase_disp, icase_vector, title,
                      min_value, max_value, data_format,
                      nlabels, labelsize, ncolors, colormap, is_fringe,
                      scale, phase,
                      arrow_scale,

                      default_title, default_min_value, default_max_value,
                      default_data_format, default_nlabels, default_labelsize,
                      default_ncolors, default_colormap,
                      default_scale, default_phase,
                      default_arrow_scale,
                      font_size=8, external_call=False):
        """
        We need to update the legend if there's been a result change request
        """
        self.external_call = external_call
        self.set_font_size(font_size)

        update_fringe = False
        update_disp = False
        update_vector = False
        #print('update_legend; fringe=%s disp=%s vector=%s' % (
            #icase_fringe, icase_disp, icase_vector))
        #print('update_legend; default: fringe=%s disp=%s vector=%s' % (
        #    self._default_icase_fringe, self._default_icase_disp, self._default_icase_vector))
        if icase_fringe != self._default_icase_fringe:
            self._icase_fringe = icase_fringe
            self._default_icase_fringe = icase_fringe
            update_fringe = True

        #is_fringe = icase_fringe is not None
        is_disp = icase_disp is not None
        is_vector = icase_vector is not None

        if icase_disp != self._default_icase_disp:
            assert isinstance(scale, float_types), 'scale=%r type=%s' % (scale, type(scale))
            #assert isinstance(default_scale, float), 'default_scale=%r' % default_scale
            self._icase_disp = icase_disp
            self._default_icase_disp = icase_disp
            self._default_scale = default_scale
            self._default_phase = default_phase
            update_disp = True

        if icase_vector != self._default_icase_vector:
            assert isinstance(arrow_scale, float_types), 'arrow_scale=%r type=%s' % (arrow_scale, type(scale))
            #assert isinstance(default_arrow_scale, float), 'default_arrow_scale=%r' % default_arrow_scale
            self._icase_vector = icase_vector
            self._default_icase_vector = icase_vector
            self._default_arrow_scale = default_arrow_scale
            update_vector = True

        #print('*update_legend; default: fringe=%s disp=%s vector=%s' % (
        #    self._default_icase_fringe, self._default_icase_disp, self._default_icase_vector))
        #print('update_fringe=%s update_disp=%s update_vector=%s' % (
        #    update_fringe, update_disp, update_vector))
        #print('is_fringe=%s is_disp=%s is_vector=%s' % (
        #    is_fringe, is_disp, is_vector))

        if update_fringe:
            self._icase_fringe = icase_fringe
            self._default_icase_fringe = icase_fringe

            self._default_title = default_title
            self._default_min = default_min_value
            self._default_max = default_max_value
            self._default_format = default_data_format
            #self._default_is_low_to_high = is_low_to_high
            self._default_is_discrete = True
            #self._default_is_horizontal = is_horizontal_scalar_bar
            self._default_nlabels = default_nlabels
            self._default_labelsize = default_labelsize
            self._default_ncolors = default_ncolors
            self._default_colormap = default_colormap
            self._is_fringe = is_fringe

            if colormap is None:
                colormap = 'jet'
            if labelsize is None:
                labelsize = ''
            if ncolors is None:
                ncolors = ''
            if nlabels is None:
                nlabels = ''

        self._update_defaults_to_blank()

        #-----------------------------------------------------------------------
        update = update_fringe or update_disp or update_vector
        if update:
            #self.scale_label.setEnabled(is_disp)
            #self.scale_edit.setEnabled(is_disp)
            #self.scale_button.setEnabled(is_disp)
            self.scale_label.setVisible(is_disp)
            self.scale_edit.setVisible(is_disp)
            self.scale_button.setVisible(is_disp)

            is_complex_disp = self._default_phase is not None
            self.phase_label.setVisible(is_complex_disp)
            self.phase_edit.setVisible(is_complex_disp)
            self.phase_button.setVisible(is_complex_disp)

            self._scale = set_cell_to_blank_if_value_is_none(self.scale_edit, scale)
            self._phase = set_cell_to_blank_if_value_is_none(self.phase_edit, phase)

        if self._default_icase_disp is None: # or self._default_icase_vector is None:
            self.animate_button.setEnabled(False)
            self.animate_button.setToolTip(ANIMATE_TOOLTIP_OFF)
        else:
            self.animate_button.setEnabled(True)
            self.animate_button.setToolTip(ANIMATE_TOOLTIP_ON)

        #-----------------------------------------------------------------------
        if update:
            #self.arrow_scale_label.setEnabled(is_vector)
            #self.arrow_scale_edit.setEnabled(is_vector)
            #self.arrow_scale_button.setEnabled(is_vector)
            self.arrow_scale_label.setVisible(is_vector)
            self.arrow_scale_edit.setVisible(is_vector)
            self.arrow_scale_button.setVisible(is_vector)
            self._arrow_scale = set_cell_to_blank_if_value_is_none(
                self.arrow_scale_edit, arrow_scale)

        #-----------------------------------------------------------------------
        if update_fringe:
            #self.on_default_title()
            #self.on_default_min()
            #self.on_default_max()
            #self.on_default_format()
            #self.on_default_scale()
            # reset defaults
            self.title_edit.setText(title)
            self.title_edit.setStyleSheet("QLineEdit{background: white;}")

            self.min_edit.setText(str(min_value))
            self.min_edit.setStyleSheet("QLineEdit{background: white;}")

            self.max_edit.setText(str(max_value))
            self.max_edit.setStyleSheet("QLineEdit{background: white;}")

            self.format_edit.setText(str(data_format))
            self.format_edit.setStyleSheet("QLineEdit{background: white;}")

            self.nlabels_edit.setText(str(nlabels))
            self.nlabels_edit.setStyleSheet("QLineEdit{background: white;}")

            self.labelsize_edit.setText(str(labelsize))
            self.labelsize_edit.setStyleSheet("QLineEdit{background: white;}")

            self.ncolors_edit.setText(str(ncolors))
            self.ncolors_edit.setStyleSheet("QLineEdit{background: white;}")

            self.colormap_edit.setCurrentIndex(colormap_keys.index(str(colormap)))

            self._set_legend_fringe(self._is_fringe)

        if update:
            self.on_apply()
            self.external_call = True
        #print('---------------------------------')

    def clear_disp(self):
        """hides dispacement blocks"""
        self._icase_disp = None
        self._default_icase_disp = None
        self.scale_label.setVisible(False)
        self.scale_edit.setVisible(False)
        self.scale_button.setVisible(False)
        self.phase_label.setVisible(False)
        self.phase_edit.setVisible(False)
        self.phase_button.setVisible(False)

    def clear_vector(self):
        """hides vector blocks"""
        self._icase_vector = None
        self._default_icase_vector = None
        self.arrow_scale_label.setVisible(False)
        self.arrow_scale_edit.setVisible(False)
        self.arrow_scale_button.setVisible(False)

    def clear(self):
        """hides fringe, displacemnt, and vector blocks"""
        self._icase_fringe = None
        self._default_icase_fringe = None
        self._set_legend_fringe(False)
        self.clear_disp()
        self.clear_vector()

    def _set_legend_fringe(self, is_fringe):
        """
        Show/hide buttons if we dont have a result.  This is used for normals.
        A result can still exist (i.e., icase_fringe is not None).
        """
        # lots of hacking for the Normal vectors
        self._is_fringe = is_fringe
        #self._default_icase_fringe = None
        enable = True
        if not is_fringe:
            enable = False

        show_title = self._icase_fringe is not None
        self.title_label.setVisible(show_title)
        self.title_edit.setVisible(show_title)
        self.title_button.setVisible(show_title)

        self.max_label.setVisible(enable)
        self.min_label.setVisible(enable)
        self.max_edit.setVisible(enable)
        self.min_edit.setVisible(enable)
        self.max_button.setVisible(enable)
        self.min_button.setVisible(enable)

        self.show_radio.setVisible(enable)
        self.hide_radio.setVisible(enable)
        self.low_to_high_radio.setVisible(enable)
        self.high_to_low_radio.setVisible(enable)

        self.format_label.setVisible(enable)
        self.format_edit.setVisible(enable)
        self.format_edit.setVisible(enable)
        self.format_button.setVisible(enable)

        self.nlabels_label.setVisible(enable)
        self.nlabels_edit.setVisible(enable)
        self.nlabels_button.setVisible(enable)

        self.ncolors_label.setVisible(enable)
        self.ncolors_edit.setVisible(enable)
        self.ncolors_button.setVisible(enable)

        self.grid2_title.setVisible(enable)
        self.vertical_radio.setVisible(enable)
        self.horizontal_radio.setVisible(enable)

        self.colormap_label.setVisible(enable)
        self.colormap_edit.setVisible(enable)
        self.colormap_button.setVisible(enable)

    def create_widgets(self, show_animation_button=True):
        """creates the menu objects"""
        # title
        self.title_label = QLabel("Title:")
        self.title_edit = QLineEdit(str(self._default_title))
        self.title_button = QPushButton("Default")

        # Min
        self.min_label = QLabel("Min:")
        self.min_edit = QLineEdit(str(self._default_min))
        self.min_button = QPushButton("Default")

        # Max
        self.max_label = QLabel("Max:")
        self.max_edit = QLineEdit(str(self._default_max))
        self.max_button = QPushButton("Default")

        #---------------------------------------
        # Format
        self.format_label = QLabel("Format (e.g. %.3f, %g, %.6e):")
        self.format_edit = QLineEdit(str(self._format))
        self.format_button = QPushButton("Default")

        #---------------------------------------
        # Scale
        self.scale_label = QLabel("True Scale:")
        self.scale_edit = QLineEdit(str(self._scale))
        self.scale_button = QPushButton("Default")
        if self._icase_disp is None:
            self.scale_label.setVisible(False)
            self.scale_edit.setVisible(False)
            self.scale_button.setVisible(False)

        # Phase
        self.phase_label = QLabel("Phase (deg):")
        self.phase_edit = QLineEdit(str(self._phase))
        self.phase_button = QPushButton("Default")
        if self._icase_disp is None or self._default_phase is None:
            self.phase_label.setVisible(False)
            self.phase_edit.setVisible(False)
            self.phase_button.setVisible(False)
            self.phase_edit.setText('0.0')

        #---------------------------------------
        self.arrow_scale_label = QLabel("Arrow Scale:")
        self.arrow_scale_edit = QLineEdit(str(self._arrow_scale))
        self.arrow_scale_button = QPushButton("Default")
        if self._icase_vector is None:
            self.arrow_scale_label.setVisible(False)
            self.arrow_scale_edit.setVisible(False)
            self.arrow_scale_button.setVisible(False)

        #tip = QtGui.QToolTip()
        #tip.setTe
        #self.format_edit.toolTip(tip)

        #---------------------------------------
        # nlabels
        self.nlabels_label = QLabel("Number of Labels:")
        self.nlabels_edit = QLineEdit(str(self._nlabels))
        self.nlabels_button = QPushButton("Default")

        self.labelsize_label = QLabel("Label Size:")
        self.labelsize_edit = QLineEdit(str(self._labelsize))
        self.labelsize_button = QPushButton("Default")

        self.ncolors_label = QLabel("Number of Colors:")
        self.ncolors_edit = QLineEdit(str(self._ncolors))
        self.ncolors_button = QPushButton("Default")

        self.colormap_label = QLabel("Color Map:")
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

        if self._icase_fringe is None:
            self.title_label.setVisible(False)
            self.title_edit.setVisible(False)
            self.title_button.setVisible(False)

        if not self._is_fringe:
            self.max_label.hide()
            self.min_label.hide()
            self.max_edit.hide()
            self.min_edit.hide()
            self.max_button.hide()
            self.min_button.hide()

            self.format_label.hide()
            self.format_edit.hide()
            self.format_button.hide()

            self.nlabels_label.hide()
            self.nlabels_edit.hide()
            self.nlabels_button.hide()

            self.ncolors_label.hide()
            self.ncolors_edit.hide()
            self.ncolors_button.hide()

            self.grid2_title.hide()
            self.vertical_radio.hide()
            self.horizontal_radio.hide()
            self.show_radio.hide()
            self.hide_radio.hide()
            self.low_to_high_radio.hide()
            self.high_to_low_radio.hide()

            self.colormap_label.hide()
            self.colormap_edit.hide()
            self.colormap_button.hide()

        self.animate_button = QPushButton('Create Animation')
        self.animate_button.setVisible(show_animation_button)
        #self.advanced_button = QPushButton('Advanced')

        if self._default_icase_disp is None: # or self._default_icase_vector is None:
            self.animate_button.setEnabled(False)
            self.animate_button.setToolTip(ANIMATE_TOOLTIP_OFF)
        else:
            self.animate_button.setEnabled(True)
            self.animate_button.setToolTip(ANIMATE_TOOLTIP_ON)

        # closing
        self.apply_button = QPushButton("Apply")
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")

    def create_layout(self):
        """displays the menu objects"""
        grid = QGridLayout()
        grid.addWidget(self.title_label, 0, 0)
        grid.addWidget(self.title_edit, 0, 1)
        grid.addWidget(self.title_button, 0, 2)

        grid.addWidget(self.min_label, 1, 0)
        grid.addWidget(self.min_edit, 1, 1)
        grid.addWidget(self.min_button, 1, 2)

        grid.addWidget(self.max_label, 2, 0)
        grid.addWidget(self.max_edit, 2, 1)
        grid.addWidget(self.max_button, 2, 2)

        grid.addWidget(self.format_label, 3, 0)
        grid.addWidget(self.format_edit, 3, 1)
        grid.addWidget(self.format_button, 3, 2)

        grid.addWidget(self.scale_label, 4, 0)
        grid.addWidget(self.scale_edit, 4, 1)
        grid.addWidget(self.scale_button, 4, 2)

        grid.addWidget(self.phase_label, 6, 0)
        grid.addWidget(self.phase_edit, 6, 1)
        grid.addWidget(self.phase_button, 6, 2)

        grid.addWidget(self.arrow_scale_label, 5, 0)
        grid.addWidget(self.arrow_scale_edit, 5, 1)
        grid.addWidget(self.arrow_scale_button, 5, 2)

        grid.addWidget(self.nlabels_label, 7, 0)
        grid.addWidget(self.nlabels_edit, 7, 1)
        grid.addWidget(self.nlabels_button, 7, 2)

        #grid.addWidget(self.labelsize_label, 6, 0)
        #grid.addWidget(self.labelsize_edit, 6, 1)
        #grid.addWidget(self.labelsize_button, 6, 2)

        grid.addWidget(self.ncolors_label, 8, 0)
        grid.addWidget(self.ncolors_edit, 8, 1)
        grid.addWidget(self.ncolors_button, 8, 2)

        grid.addWidget(self.colormap_label, 9, 0)
        grid.addWidget(self.colormap_edit, 9, 1)
        grid.addWidget(self.colormap_button, 9, 2)

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

        self.setLayout(vbox)

    def set_connections(self):
        """creates the actions for the menu"""
        self.title_button.clicked.connect(self.on_default_title)
        self.min_button.clicked.connect(self.on_default_min)
        self.max_button.clicked.connect(self.on_default_max)
        self.format_button.clicked.connect(self.on_default_format)
        self.scale_button.clicked.connect(self.on_default_scale)
        self.arrow_scale_button.clicked.connect(self.on_default_arrow_scale)
        self.phase_button.clicked.connect(self.on_default_phase)

        self.nlabels_button.clicked.connect(self.on_default_nlabels)
        self.labelsize_button.clicked.connect(self.on_default_labelsize)
        self.ncolors_button.clicked.connect(self.on_default_ncolors)
        self.colormap_button.clicked.connect(self.on_default_colormap)

        self.animate_button.clicked.connect(self.on_animate)

        self.show_radio.clicked.connect(self.on_show_hide)
        self.hide_radio.clicked.connect(self.on_show_hide)

        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)

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
        self.title_edit.setFont(font)
        self.min_edit.setFont(font)
        self.max_edit.setFont(font)
        self.format_edit.setFont(font)
        self.scale_edit.setFont(font)
        self.phase_edit.setFont(font)
        self.nlabels_edit.setFont(font)
        self.labelsize_edit.setFont(font)
        self.ncolors_edit.setFont(font)

    def on_animate(self):
        """opens the animation window"""
        title, flag0 = check_name_str(self.title_edit)
        if not flag0:
            return

        scale, flag1 = check_float(self.scale_edit)
        if not flag1:
            scale = self._scale

        data = {
            'font_size' : self.out_data['font_size'],
            'icase_fringe' : self._icase_fringe,
            'icase_disp' : self._icase_disp,
            'icase_vector' : self._icase_vector,
            'title' : title,
            'time' : 2,
            'frames/sec' : 30,
            'resolution' : 1,
            'iframe' : 0,
            'scale' : scale,
            'default_scale' : self._default_scale,

            'arrow_scale' : scale,
            'default_arrow_scale' : self._default_arrow_scale,

            'is_scale' : self._default_phase is None,
            'phase' : self._phase,
            'default_phase' : self._default_phase,
            'dirname' : os.path.abspath(os.getcwd()),
            'clicked_ok' : False,
            'close' : False,
        }
        self.win_parent.legend_obj.set_animation_window(data)

    def on_default_title(self):
        """action when user clicks 'Default' for title"""
        title = str(self._default_title)
        self.title_edit.setText(title)
        self.title_edit.setStyleSheet("QLineEdit{background: white;}")

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

    def on_default_arrow_scale(self):
        """action when user clicks 'Default' for arrow_scale factor"""
        self.arrow_scale_edit.setText(str(self._default_arrow_scale))
        self.arrow_scale_edit.setStyleSheet("QLineEdit{background: white;}")

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

    def on_show_hide(self):
        """action when user clicks the 'Show/Hide' radio button"""
        #self.colormap_edit.setCurrentIndex(colormap_keys.index(self._default_colormap))
        is_shown = self.show_radio.isChecked()
        self.nlabels_edit.setEnabled(is_shown)
        self.nlabels_button.setEnabled(is_shown)
        self.ncolors_edit.setEnabled(is_shown)
        self.ncolors_button.setEnabled(is_shown)
        self.high_to_low_radio.setEnabled(is_shown)
        self.low_to_high_radio.setEnabled(is_shown)
        self.colormap_edit.setEnabled(is_shown)
        self.colormap_button.setEnabled(is_shown)
        self.vertical_radio.setEnabled(is_shown)
        self.horizontal_radio.setEnabled(is_shown)

    def show_legend(self):
        """shows the legend"""
        self._set_legend(True)

    def hide_legend(self):
        """hides the legend"""
        self._set_legend(False)

    def _set_legend(self, is_shown):
        """shows/hides the legend"""
        self.show_radio.setChecked(is_shown)
        self.hide_radio.setChecked(not is_shown)
        #if self.is_gui:
            #self.win_parent.scalar_bar_actor.SetVisibility(is_shown)

    def on_validate(self):
        """checks to see if the ``on_apply`` method can be called"""
        show_title = self._icase_fringe is not None
        flag_title = True
        title_value = ''
        if show_title:
            title_value, flag_title = check_name_str(self.title_edit)

        flag_fringe = True
        min_value = max_value = format_value = nlabels = ncolors = labelsize = colormap = None
        if self._icase_fringe is not None:
            min_value, flag1 = check_float(self.min_edit)
            max_value, flag2 = check_float(self.max_edit)
            format_value, flag3 = check_format(self.format_edit)
            nlabels, flag4 = check_positive_int_or_blank(self.nlabels_edit)
            ncolors, flag5 = check_positive_int_or_blank(self.ncolors_edit)
            labelsize, flag6 = check_positive_int_or_blank(self.labelsize_edit)
            colormap = str(self.colormap_edit.currentText())
            if flag3 and 'i' in format_value:
                format_value = '%i'
            flag_fringe = all([flag1, flag2, flag3, flag4, flag5, flag6])

        flag_disp = True
        scale = phase = None
        if self._icase_disp is not None:
            scale, flag1 = check_float(self.scale_edit)
            phase, flag2 = check_float(self.phase_edit)
            assert isinstance(scale, float), scale
            flag_disp = all([flag1, flag2])

        flag_vector = True
        arrow_scale = None
        if self._icase_vector is not None:
            arrow_scale, flag_vector = check_float(self.arrow_scale_edit)

        if all([flag_title, flag_fringe, flag_disp, flag_vector]):
            self.out_data['title'] = title_value
            self.out_data['min_value'] = min_value
            self.out_data['max_value'] = max_value
            self.out_data['format'] = format_value
            self.out_data['scale'] = scale
            self.out_data['phase'] = phase

            self.out_data['arrow_scale'] = arrow_scale

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
            #print("title = %r" % self.title_edit.text())
            #print("min = %r" % self.min_edit.text())
            #print("max = %r" % self.max_edit.text())
            #print("format = %r" % self.format_edit.text())
            return True
        return False

    def on_apply(self):
        """applies the current values"""
        passed = self.on_validate()
        if passed and self.external_call:
            self.win_parent.legend_obj._apply_legend(self.out_data)
        self.external_call = True
        return passed

    def on_ok(self):
        """applies the current values and closes the window"""
        passed = self.on_apply()
        if passed:
            self.close()
            #self.destroy()

    def on_cancel(self):
        """closes the windows and reverts the legend"""
        self.out_data['close'] = True
        self.close()

def set_cell_to_blank_if_value_is_none(cell_edit, value):
    if value is None:
        cell_edit.setText('0.0')
    else:
        cell_edit.setText(str(value))
    cell_edit.setStyleSheet("QLineEdit{background: white;}")
    return value

#def check_colormap(cell):
    #text = str(cell.text()).strip()
    #if text in colormap_keys:
        #cell.setStyleSheet("QLineEdit{background: white;}")
        #return text, True
    #else:
        #cell.setStyleSheet("QLineEdit{background: red;}")
        #return None, False


def main(): # pragma: no cover
    """tests out the legend"""
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
        'icase_fringe' : 1,
        'icase_disp' : 62,
        'icase_vector' : 3,
        'is_fringe': True,  # False=normals or no fringe

        'font_size' : 8,
        'title' : 'asdf',
        'min_value' : 0.,
        'max_value' : 10,
        'scale' : 1e-12,
        'default_scale' : 1.0,

        'arrow_scale' : 2e-2,
        'default_arrow_scale' : 2.0,

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
