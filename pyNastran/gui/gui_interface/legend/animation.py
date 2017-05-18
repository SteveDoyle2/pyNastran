"""
defines:
 - LegendPropertiesWindow
"""
from __future__ import print_function
import os
from collections import OrderedDict

from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4.QtGui import (
        QApplication, QLabel, QPushButton, QLineEdit, QWidget, QRadioButton,
        QButtonGroup, QGridLayout, QHBoxLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox,
        QCheckBox, QGroupBox, QComboBox)
elif qt_version == 5:
    from PyQt5.QtWidgets import (
        QApplication, QLabel, QPushButton, QLineEdit, QWidget, QRadioButton,
        QButtonGroup, QGridLayout, QHBoxLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox,
        QCheckBox, QGroupBox, QComboBox)
elif qt_version == 'pyside':
    from PySide.QtGui import (
        QApplication, QLabel, QPushButton, QLineEdit, QWidget, QRadioButton,
        QButtonGroup, QGridLayout, QHBoxLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox,
        QCheckBox, QGroupBox, QComboBox)
else:
    raise NotImplementedError('qt_version = %r' % qt_version)

from pyNastran.gui.gui_interface.common import PyDialog
from pyNastran.gui.gui_utils.dialogs import open_directory_dialog, open_file_dialog
from pyNastran.gui.gui_utils.write_gif import IS_IMAGEIO


ANIMATION_PROFILES = OrderedDict()
ANIMATION_PROFILES['0 to Scale'] = [0., 1.]
#ANIMATION_PROFILES['0 to Scale to 0'] = [0., 1., 0.]
#ANIMATION_PROFILES['0 to Scale to -Scale to 0'] = [0., 1., -1., 0.]
ANIMATION_PROFILES['-Scale to Scale'] = [-1., 1.]
#ANIMATION_PROFILES['-Scale to Scale to -Scale'] = [-1., 1., -1.]
#ANIMATION_PROFILES['CSV Profile  (not supported)'] = None

class AnimationWindow(PyDialog):
    """
    +-------------------+
    | Animation         |
    +-------------------------+
    | scale   ______  Default |
    | time    ______  Default |
    |                         |
    | nframes ______  Default |
    | resolu. ______  Default |
    | Dir     ______  Browse  |
    | iFrame  ______          |
    |                         |
    | Animations:             |
    | o Scale, Phase, Time    |  # TODO: add time
    |                         |
    | x delete images         |
    | x repeat                |  # TODO: change to an integer
    | x make gif              |
    |                         |
    |      Step, RunAll       |
    |         Close           |
    +-------------------------+

    TODO: add key-frame support
    """
    def __init__(self, data, win_parent=None):
        PyDialog.__init__(self, data, win_parent)
        self.set_font_size(data['font_size'])
        self.istep = 0
        self._animate_type = 'time'

        self._updated_animation = False
        self._icase = data['icase']
        self._default_name = data['name']
        self._default_time = data['time']
        self._default_fps = data['frames/sec']
        self._default_resolution = data['resolution']

        self._scale = data['scale']
        self._default_scale = data['default_scale']
        self._default_is_scale = data['is_scale']

        self._phase = data['phase']
        self._default_phase = data['default_phase']

        self._default_dirname = data['dirname']
        self._default_gif_name = os.path.join(self._default_dirname, data['name'] + '.gif')

        self.animation_types = [
            'Animate Scale',
        ]
            #'Animate Phase',
            #'Animate Time',
            #'Animate Frequency Sweep'
        #]

        self.setWindowTitle('Animate Model')
        self.create_widgets()
        self.create_layout()
        self.set_connections()

        self.is_gui = False
        if hasattr(self.win_parent, '_updated_legend'):
            self.win_parent.is_animate_open = True
            self.is_gui = True


    def create_widgets(self):
        """creates the menu objects"""
        self.scale = QLabel("Scale:")
        self.scale_edit = QLineEdit(str(self._scale))
        self.scale_button = QPushButton("Default")
        self.scale_edit.setToolTip('Scale factor of the "deflection"')
        self.scale_button.setToolTip('Sets the scale factor of the gif to %s' % self._scale)

        self.time = QLabel("Total Time (sec):")
        self.time_edit = QDoubleSpinBox(self)
        self.time_edit.setValue(self._default_time)
        self.time_edit.setRange(0.1, 10.0)
        self.time_edit.setDecimals(2)
        self.time_edit.setSingleStep(0.1)
        self.time_button = QPushButton("Default")
        self.time_edit.setToolTip("Total time of the gif")
        self.time_button.setToolTip('Sets the total time of the gif to %.2f' % self._default_time)

        self.fps = QLabel("Frames/Second:")
        self.fps_edit = QSpinBox(self)
        self.fps_edit.setRange(10, 60)
        self.fps_edit.setSingleStep(1)
        self.fps_edit.setValue(self._default_fps)
        self.fps_button = QPushButton("Default")
        self.fps_edit.setToolTip("A higher FPS is smoother, but may not play well for large gifs")
        self.fps_button.setToolTip('Sets the FPS to %s' % self._default_fps)

        self.resolution = QLabel("Resolution Scale:")
        self.resolution_edit = QSpinBox(self)
        self.resolution_edit.setRange(1, 5)
        self.resolution_edit.setSingleStep(1)
        self.resolution_edit.setValue(self._default_resolution)
        self.resolution_button = QPushButton("Default")
        self.resolution_edit.setToolTip('Scales the window resolution by an integer factor')
        self.resolution_button.setToolTip('Sets the resolution to %s' % self._default_resolution)

        #-----------------
        # Time plot
        icase_max = 1000  # TODO: update 1000
        self.icase_start = QLabel("iCase Start:")
        self.icase_start_edit = QSpinBox(self)
        self.icase_start_edit.setRange(0, icase_max)
        self.icase_start_edit.setSingleStep(1)
        self.icase_start_edit.setValue(self._icase)
        self.icase_start_button = QPushButton("Default")

        self.icase_end = QLabel("iCase End:")
        self.icase_end_edit = QSpinBox(self)
        self.icase_end_edit.setRange(0, icase_max)
        self.icase_end_edit.setSingleStep(1)
        self.icase_end_edit.setValue(self._icase)
        self.icase_end_button = QPushButton("Default")

        self.icase_delta = QLabel("iCase Delta:")
        self.icase_delta_edit = QSpinBox(self)
        self.icase_delta_edit.setRange(1, icase_max)
        self.icase_delta_edit.setSingleStep(1)
        self.icase_delta_edit.setValue(1)
        self.icase_delta_button = QPushButton("Default")

        self.min_value = QLabel("Min Value:")
        self.min_value_edit = QLineEdit(str(0.))
        #self.min_value_edit.setRange(1, 1000)
        #self.min_value_edit.setSingleStep(1)
        #self.min_value_edit.setValue(1)
        self.min_value_button = QPushButton("Default")

        self.max_value = QLabel("Max Value:")
        self.max_value_edit = QLineEdit(str(1.))
        #self.min_value_edit.setRange(1, 1000)  # TODO: update 1000
        #self.min_value_edit.setSingleStep(1)
        #self.min_value_edit.setValue(1)
        self.max_value_button = QPushButton("Default")

        self.icase_start_edit.setToolTip('The first frame of the animation')
        self.icase_end_edit.setToolTip(
            'The last frame of the animation\n'
            'Assumes icase_start + nframes * icase_delta = icase_end')
        self.icase_delta_edit.setToolTip('The frame step size (to skip non-consecutive results)')

        self.min_value_edit.setToolTip('Min value of the legend (not supported)')
        self.max_value_edit.setToolTip('Max value of the legend (not supported)')
        #'time' : 0.,
        #'default_time' : 0,
        #'icase_start' : 10,
        #'icase_delta' : 3,
        #'min_value' : 0.,
        #'max_value' : 1000.,

        self.browse_folder = QLabel('Output Directory:')
        self.browse_folder_edit = QLineEdit(str(self._default_dirname))
        self.browse_folder_button = QPushButton('Browse')
        self.browse_folder_edit.setToolTip('Location to save the png/gif files')

        self.gif = QLabel("Gif Filename:")
        self.gif_edit = QLineEdit(str(self._default_name + '.gif'))
        self.gif_button = QPushButton('Default')
        self.gif_edit.setToolTip('Name of the gif')
        self.gif_button.setToolTip('Sets the name of the gif to %s.gif' % self._default_name)

        # scale / phase
        self.animate_scale_radio = QRadioButton("Animate Scale")
        self.animate_phase_radio = QRadioButton("Animate Phase")
        self.animate_time_radio = QRadioButton("Animate Time")
        self.animate_freq_sweeep_radio = QRadioButton("Animate Frequency Sweep")
        self.animate_scale_radio.setToolTip(
            'Animates the scale factor based on the "Animation Type"')
        self.animate_time_radio.setToolTip('Animates the time/load/mode step')

        self.animate_scale_radio.setChecked(self._default_is_scale)
        self.animate_phase_radio.setChecked(not self._default_is_scale)
        self.animate_time_radio.setChecked(False)

        msg = 'Scale : Animates the scale factor based on the "Animation Profile"\n'
        if self._default_phase is None:
            self.animate_phase_radio.setDisabled(True)
            self.animate_phase_radio.setToolTip('Animates the phase angle '
                                                '(only for complex results)')
            msg += 'Phase : Animates the phase angle (only for complex results)\n'
        else:
            self.animate_phase_radio.setToolTip("Animates the phase angle")
            msg += 'Phase : Animates the phase angle\n'
        msg += (
            'Time : Animates the time/load/mode step\n'
            'Freq Sweep : Animates a complex result across a range of frequencies '
            '(not supported)\n'
        )

        self.animate_freq_sweeep_radio.setDisabled(True)
        self.animate_freq_sweeep_radio.setToolTip(
            'Animates a complex result across a range of frequencies (not supported)')
        self.animation_type = QLabel("Animation Type:")
        animation_type = OrderedDict()
        #scale_msg = 'Scale\n'
        #phase_msg = 'Phase\n'
        #time_msg = 'Time\n'
        #animation_types = [
            #('Animate Scale', scale_msg),
            #('Animate Phase', phase_msg),
            #('Animate Time', time_msg),
            ##'Animate Frequency Sweep'
        #]

        if self._phase is not None:
            self.animation_types.append('Animate Phase')
        self.animation_types.append('Animate Time')

        self.animation_profile = QLabel("Animation Profile:")

        self.animation_profile_edit = QComboBox()
        for animation_profile in ANIMATION_PROFILES:
            self.animation_profile_edit.addItem(animation_profile)
        self.animation_profile_edit.setToolTip('The profile for a scaled GIF')

        self.animation_type_edit = QComboBox()
        # TODO: add a tooltip for each item
        for animation_type in self.animation_types:
            self.animation_type_edit.addItem(animation_type)
        #self.animation_type_edit.setToolTip('The profile for a scaled GIF')
        self.animation_type_edit.setToolTip(msg)

        self.csv_profile = QLabel("CSV profile:")
        self.csv_profile_edit = QLineEdit()
        self.csv_profile_browse_button = QPushButton('Browse')
        self.csv_profile_edit.setToolTip(
            'The path to the CSV file of (Scale1, Scale2, Scale3, ...)')

        widget = QWidget(self)
        horizontal_vertical_group = QButtonGroup(widget)
        horizontal_vertical_group.addButton(self.animate_scale_radio)
        horizontal_vertical_group.addButton(self.animate_phase_radio)
        horizontal_vertical_group.addButton(self.animate_time_radio)
        horizontal_vertical_group.addButton(self.animate_freq_sweeep_radio)

        # one / two sided
        self.onesided_radio = QRadioButton("One Sided")
        self.onesided_radio.setToolTip(
            "A one sided gif doesn't return to the starting point (e.g., 0 to 360 degrees)")
        self.twosided_radio = QRadioButton("Two Sided")
        self.twosided_radio.setToolTip(
            'A two sided gif returns to the starting point (e.g., 0 to 10 to 0)')

        if self._default_phase is None:
            self.onesided_radio.setChecked(False)
            self.twosided_radio.setChecked(True)
        else:
            self.onesided_radio.setChecked(True)
            self.twosided_radio.setChecked(False)
        widget = QWidget(self)
        horizontal_vertical_group = QButtonGroup(widget)
        horizontal_vertical_group.addButton(self.onesided_radio)
        horizontal_vertical_group.addButton(self.twosided_radio)

        # make images
        self.make_images_checkbox = QCheckBox("Make images?")
        self.make_images_checkbox.setChecked(True)

        # make images
        self.overwrite_images_checkbox = QCheckBox("Overwrite images?")
        self.overwrite_images_checkbox.setChecked(True)

        # delete images when finished
        self.delete_images_checkbox = QCheckBox("Delete images when finished?")
        self.delete_images_checkbox.setChecked(True)

        # endless loop
        self.repeat_checkbox = QCheckBox("Repeat?")
        self.repeat_checkbox.setChecked(True)
        self.repeat_checkbox.setToolTip("Repeating creates an infinitely looping gif")

        # endless loop
        self.make_gif_checkbox = QCheckBox("Make Gif?")
        if IS_IMAGEIO:
            self.make_gif_checkbox.setChecked(True)
        else:
            self.make_gif_checkbox.setChecked(False)
            self.make_gif_checkbox.setEnabled(False)
            self.make_gif_checkbox.setToolTip('imageio is not available; install it')

        # bottom buttons
        self.step_button = QPushButton("Step")
        self.run_button = QPushButton("Run All")

        #self.apply_button = QPushButton("Apply")
        #self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Close")

        #self.set_grid_time(enabled=False)
        #self.set_grid_scale(enabled=self._default_is_scale)
        if self._default_phase:
            self.on_animate_phase(force=True)
        else:
            self.on_animate_scale(force=True)

    def set_connections(self):
        """creates button actions"""
        self.scale_button.clicked.connect(self.on_default_scale)
        self.time_button.clicked.connect(self.on_default_time)

        self.fps_button.clicked.connect(self.on_default_fps)
        self.resolution_button.clicked.connect(self.on_default_resolution)
        self.browse_folder_button.clicked.connect(self.on_browse_folder)
        self.csv_profile_browse_button.clicked.connect(self.on_browse_csv)
        self.gif_button.clicked.connect(self.on_default_name)

        self.step_button.clicked.connect(self.on_step)
        self.run_button.clicked.connect(self.on_run)

        #self.animate_scale_radio.clicked.connect(self.on_animate_scale)
        #self.animate_phase_radio.clicked.connect(self.on_animate_phase)
        #self.animate_time_radio.clicked.connect(self.on_animate_time)
        self.animation_type_edit.currentIndexChanged.connect(self.on_animate)
        #self.animate_freq_sweeep_radio

        #self.apply_button.clicked.connect(self.on_apply)
        #self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)

    def on_animate(self, value):
        """
        animate pulldown

        Parameters
        ----------
        value : int
            index in animation_types
        """
        #animation_types = ['Animate Scale', 'Animate Phase', 'Animate Time',
                           #'Animate Frequency Sweep']
        animation_type = self.animation_types[value]
        if animation_type == 'Animate Scale':
            self.on_animate_scale()
        elif animation_type == 'Animate Phase':
            self.on_animate_phase()
        elif animation_type == 'Animate Time':
            self.on_animate_time()
        else:
            raise NotImplementedError('value = ', value)

    def on_animate_time(self, force=False):
        """enables the secondary input"""
        #print('on_animate_time')
        if self._animate_type == 'scale' or force:
            self.set_grid_scale(False, 'time')
        self.set_grid_time(True, 'time')
        self._animate_type = 'time'

    def on_animate_scale(self, force=False):
        """enables the secondary input"""
        #print('on_animate_scale')
        self.set_grid_scale(True, 'scale')
        if self._animate_type == 'time' or force:
            self.set_grid_time(False, 'scale')
        self._animate_type = 'scale'

    def on_animate_phase(self, force=False):
        """enables the secondary input"""
        #print('on_animate_phase')
        if self._animate_type == 'scale' or force:
            self.set_grid_scale(False, 'phase')
        if self._animate_type == 'time' or force:
            self.set_grid_time(False, 'phase')
        self._animate_type = 'phase'

    def set_grid_scale(self, enabled=True, word=''):
        """enables/disables the secondary input"""
        #print('%s-set_grid_scale; enabled = %r' % (word, enabled))
        self.animation_profile.setEnabled(enabled)
        self.animation_profile_edit.setEnabled(enabled)

        # TODO: doesn't work...
        #self.csv_profile.setEnabled(enabled)
        #self.csv_profile_edit.setEnabled(enabled)
        #self.csv_profile_button.setEnabled(enabled)


    def set_grid_time(self, enabled=True, word=''):
        """enables/disables the secondary input"""
        #print('%s-set_grid_time; enabled = %r' % (word, enabled))
        self.icase_start.setEnabled(enabled)
        self.icase_start_edit.setEnabled(enabled)
        self.icase_start_button.setEnabled(enabled)

        self.icase_end.setEnabled(enabled)
        self.icase_end_edit.setEnabled(enabled)
        self.icase_end_button.setEnabled(enabled)

        self.icase_delta.setEnabled(enabled)
        self.icase_delta_edit.setEnabled(enabled)
        self.icase_delta_button.setEnabled(enabled)

        self.min_value.setEnabled(enabled)
        self.min_value_edit.setEnabled(enabled)
        self.min_value_button.setEnabled(enabled)

        self.max_value.setEnabled(enabled)
        self.max_value_edit.setEnabled(enabled)
        self.max_value_button.setEnabled(enabled)

    def on_browse_folder(self):
        """opens a folder dialog"""
        dirname = open_directory_dialog(self, 'Select a Directory')
        if not dirname:
            return
        self.browse_folder_edit.setText(dirname)

    def on_browse_csv(self):
        """opens a file dialog"""
        default_filename = ''
        file_types = 'Delimited Text (*.txt; *.dat; *.csv)'
        dirname = open_file_dialog(self, 'Select a CSV File', default_filename, file_types)
        if not dirname:
            return
        self.csv_profile_browse_button.setText(dirname)

    def on_default_name(self):
        """sets the default gif name"""
        self.gif_edit.setText(self._default_name + '.gif')

    def on_default_scale(self):
        """sets the default displacement scale factor"""
        self.scale_edit.setText(str(self._default_scale))
        self.scale_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_time(self):
        """sets the default gif time"""
        self.time_edit.setValue(self._default_time)

    def on_default_fps(self):
        """sets the default FPS"""
        self.fps_edit.setValue(self._default_fps)

    def on_default_resolution(self):
        """sets the default image resolution scale factor"""
        self.resolution_edit.setValue(self._default_resolution)

    def create_layout(self):
        """displays the menu objects"""
        grid = QGridLayout()

        grid.addWidget(self.scale, 0, 0)
        grid.addWidget(self.scale_edit, 0, 1)
        grid.addWidget(self.scale_button, 0, 2)

        grid.addWidget(self.time, 1, 0)
        grid.addWidget(self.time_edit, 1, 1)
        grid.addWidget(self.time_button, 1, 2)

        # spacer
        spacer = QLabel('')

        grid.addWidget(self.fps, 3, 0)
        grid.addWidget(self.fps_edit, 3, 1)
        grid.addWidget(self.fps_button, 3, 2)

        grid.addWidget(self.resolution, 4, 0)
        grid.addWidget(self.resolution_edit, 4, 1)
        grid.addWidget(self.resolution_button, 4, 2)

        grid.addWidget(self.browse_folder, 5, 0)
        grid.addWidget(self.browse_folder_edit, 5, 1)
        grid.addWidget(self.browse_folder_button, 5, 2)

        grid.addWidget(self.gif, 6, 0)
        grid.addWidget(self.gif_edit, 6, 1)
        grid.addWidget(self.gif_button, 6, 2)

        grid.addWidget(self.animation_type, 7, 0)
        grid.addWidget(self.animation_type_edit, 7, 1)

        grid.addWidget(spacer, 8, 0)

        #----------
        #Time
        grid_time = QGridLayout()
        grid_time.addWidget(self.icase_start, 0, 0)
        grid_time.addWidget(self.icase_start_edit, 0, 1)
        grid_time.addWidget(self.icase_start_button, 0, 2)

        grid_time.addWidget(self.icase_end, 1, 0)
        grid_time.addWidget(self.icase_end_edit, 1, 1)
        grid_time.addWidget(self.icase_end_button, 1, 2)

        grid_time.addWidget(self.icase_delta, 2, 0)
        grid_time.addWidget(self.icase_delta_edit, 2, 1)
        grid_time.addWidget(self.icase_delta_button, 2, 2)

        #grid_time.addWidget(self.min_value, 3, 0)
        #grid_time.addWidget(self.min_value_edit, 3, 1)
        #grid_time.addWidget(self.min_value_button, 3, 2)

        #grid_time.addWidget(self.max_value, 4, 0)
        #grid_time.addWidget(self.max_value_edit, 4, 1)
        #grid_time.addWidget(self.max_value_button, 4, 2)
        grid_time.addWidget(spacer, 5, 0)

        #--------------
        grid_scale = QGridLayout()
        grid_scale.addWidget(self.animation_profile, 0, 0)
        grid_scale.addWidget(self.animation_profile_edit, 0, 1)

        #grid_scale.addWidget(self.csv_profile, 1, 0)
        #grid_scale.addWidget(self.csv_profile_edit, 1, 1)
        #grid_scale.addWidget(self.csv_profile_browse_button, 1, 2)

        self.csv_profile = QLabel("CSV profile:")
        self.csv_profile_edit = QLineEdit()
        self.csv_profile_button = QPushButton('Browse')

        #box_time = QVBoxLayout()
        # TODO: It's super annoying that the animate time box doesn't
        #       line up with the previous box
        box_scale = QGroupBox('Animate Scale')
        box_scale.setLayout(grid_scale)

        box_time = QGroupBox('Animate Time')
        box_time.setLayout(grid_time)
        #----------


        grid2 = QGridLayout()
        #grid2.addWidget(self.animate_scale_radio, 8, 0)
        #grid2.addWidget(self.animate_phase_radio, 8, 1)
        #grid2.addWidget(self.animate_time_radio, 8, 2)
        #grid2.addWidget(self.animate_freq_sweeep_radio, 8, 3)

        grid2.addWidget(self.make_images_checkbox, 10, 0)
        #grid2.addWidget(self.overwrite_images_checkbox, 10, 0)
        grid2.addWidget(self.delete_images_checkbox, 10, 1)
        grid2.addWidget(self.make_gif_checkbox, 10, 2)
        grid2.addWidget(self.repeat_checkbox, 11, 0)


        grid2.addWidget(spacer, 12, 0)
        grid_hbox = QHBoxLayout()
        grid_hbox.addWidget(spacer)
        grid_hbox.addLayout(grid2)
        grid_hbox.addWidget(spacer)

        # bottom buttons
        step_run_box = QHBoxLayout()
        step_run_box.addWidget(self.step_button)
        step_run_box.addWidget(self.run_button)

        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.cancel_button)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        vbox.addWidget(box_scale)
        vbox.addWidget(box_time)
        #vbox.addLayout(checkboxes)
        vbox.addLayout(grid_hbox)
        vbox.addStretch()
        vbox.addLayout(step_run_box)
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def on_step(self):
        """click the Step button"""
        passed, validate_out = self.on_validate()
        if passed:
            self._make_gif(validate_out, istep=self.istep)
            self.istep += 1

    def on_run(self):
        """click the Run button"""
        self.istep = 0
        passed, validate_out = self.on_validate()
        if passed:
            self._make_gif(validate_out, istep=None)
        return passed

    def _make_gif(self, validate_out, istep=None):
        """interface for making the gif"""
        scale, time, fps, magnify, output_dir, gifbase = validate_out
        if gifbase.lower().endswith('.gif'):
            gifbase = gifbase[:-4]
        gif_filename = os.path.join(output_dir, gifbase + '.gif')

        animate_scale = self.animate_scale_radio.isChecked()
        animate_phase = self.animate_phase_radio.isChecked()
        animate_time = self.animate_time_radio.isChecked()

        animate_scale = False
        animate_phase = False
        animate_time = False
        if self._animate_type == 'scale':
            animate_scale = True
        elif self._animate_type == 'phase':
            animate_phase = True
        elif self._animate_type == 'time':
            animate_time = True
        else:
            raise NotImplementedError(self._animate_type)

        make_images = self.make_images_checkbox.isChecked()
        delete_images = self.delete_images_checkbox.isChecked()
        make_gif = self.make_gif_checkbox.isChecked()
        key = str(self.animation_profile_edit.currentText())
        #profile = ANIMATION_PROFILES[key]

        #ANIMATION_PROFILES['0 to Scale'] = [0., 1.]
        #ANIMATION_PROFILES['0 to Scale to 0'] = [0., 1., 0.]
        if key == '0 to Scale':
            onesided = True
        else:
            onesided = False

        icase_start = self.icase_start_edit.value()
        icase_end = self.icase_end_edit.value()
        icase_delta = self.icase_delta_edit.value()

        #onesided = self.onesided_radio.isChecked()
        bool_repeat = self.repeat_checkbox.isChecked()  # TODO: change this to an integer
        if bool_repeat:
            nrepeat = 0
        else:
            nrepeat = 1
        #self.out_data['is_shown'] = self.show_radio.isChecked()
        icase = self._icase
        if self.is_gui:
            self.win_parent.win_parent.make_gif(
                gif_filename, scale, istep=istep,
                animate_scale=animate_scale, animate_phase=animate_phase, animate_time=animate_time,
                icase=icase, icase_start=icase_start, icase_end=icase_end, icase_delta=icase_delta,
                time=time, onesided=onesided,
                nrepeat=nrepeat, fps=fps, magnify=magnify,
                make_images=make_images, delete_images=delete_images, make_gif=make_gif,
            )

        self.out_data['clicked_ok'] = True
        self.out_data['close'] = True

    def on_validate(self):
        """checks to see if the input is valid"""
        scale, flag0 = self.check_float(self.scale_edit)
        time, flag1 = self.check_float(self.time_edit)
        fps, flag2 = self.check_float(self.fps_edit)
        magnify, flag3 = self.check_int(self.resolution_edit)
        output_dir, flag4 = self.check_path(self.browse_folder_edit)
        gifbase, flag5 = self.check_name(self.gif_edit)
        passed = all([flag0, flag1, flag2, flag3, flag4, flag5])
        return passed, (scale, time, fps, magnify, output_dir, gifbase)

    @staticmethod
    def check_name(cell):
        """verifies that the data is string-able"""
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

    def check_path(self, cell):
        """verifies that the path exists"""
        text, passed = self.check_name(cell)
        if not passed:
            return None, False

        if os.path.exists(text):
            cell.setStyleSheet("QLineEdit{background: white;}")
            return text, True
        else:
            cell.setStyleSheet("QLineEdit{background: red;}")
            return None, False

    #def on_ok(self):
        #"""click the OK button"""
        #passed = self.on_apply()
        #if passed:
            #self.win_parent._animation_window_shown = False
            #self.close()
            ##self.destroy()

    def on_cancel(self):
        """click the Cancel button"""
        self.out_data['close'] = True
        self.close()


def main(): # pragma: no cover
    """test example for AnimationWindow"""
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)


    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    #The Main window

    data2 = {
        'font_size' : 8,
        'icase' : 1,
        'name' : 'cat',
        'time' : 2,
        'frames/sec' : 30,
        'resolution' : 1,
        'iframe' : 0,
        'is_scale' : False,
        'dirname' : os.getcwd(),
        'scale' : 2.0,
        'default_scale' : 10,
        #'phase' : 0.,
        'phase' : None,
        'default_phase' : 120.,
        #'default_phase' : None,

        #'start_time' : 0.,
        #'end_time' : 0.,
        'default_time' : 0.,
        'icase_start' : 10,
        'icase_delta' : 3,
        'stress_min' : 0.,
        'stress_max' : 1000.,
    }
    main_window = AnimationWindow(data2)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__": # pragma: no cover
    main()
