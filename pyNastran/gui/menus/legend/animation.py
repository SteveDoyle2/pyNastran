"""
defines:
 - AnimationWindow

"""
import os
from collections import OrderedDict

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QApplication, QLabel, QPushButton, QLineEdit, QRadioButton,
    #QWidget, QButtonGroup,
    QGridLayout, QHBoxLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox,
    QCheckBox, QGroupBox, QComboBox, QFileDialog)
from qtpy.compat import getexistingdirectory

from pyNastran.gui.utils.qt.pydialog import (
    PyDialog, check_int, check_float, check_name_str, check_path,
    set_combo_box_text)
from pyNastran.gui.utils.qt.dialogs import open_file_dialog
from pyNastran.gui.menus.results_sidebar import ResultsWindow
from pyNastran.gui.menus.results_sidebar_utils import (
    get_cases_from_tree, #build_pruned_tree
)

from pyNastran.gui.menus.legend.write_gif import IS_IMAGEIO


ANIMATION_PROFILES = [
    #'0 to Scale',
    #'-Scale to Scale',
    '0 to Scale to 0',
    '-Scale to Scale to -Scale',
    '0 to Scale to -Scale to 0',
    'Sinusoidal: 0 to Scale to -Scale to 0',
    'Sinusoidal: Scale to -Scale to Scale',
]
#ANIMATION_PROFILES['-Scale to Scale to -Scale'] = [-1., 1., -1.]
#ANIMATION_PROFILES['CSV Profile  (not supported)'] = None

IS_TIME_FRINGE = False
HIDE_WHEN_INACTIVE = False
IS_RESULTS_SELECTOR = True
class AnimationWindow(PyDialog):
    """
    +-------------------+
    | Animation         |
    +-------------------------+
    | icase   ______          |
    | scale   ______  Default |
    | time    ______  Default |
    |                         |
    | nframes ______  Default |
    | resolu. ______  Default |
    | Dir     ______  Browse  |
    | iFrame  ______          |
    |                         |
    | Animations:             |
    | o Scale, Phase, Time    |
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
    def __init__(self, data, win_parent=None, fringe_cases=None, is_gui_parent=False):
        PyDialog.__init__(self, data, win_parent)

        # is the parent the gui?
        self.is_gui_parent = False
        self.fringe_cases = fringe_cases
        self.set_font_size(data['font_size'])
        self.istep = 0
        self._animate_type = 'time'

        self._updated_animation = False
        self._active_deformation = 0.
        self._icase_fringe = data['icase_fringe']
        self._icase_disp = data['icase_disp']
        self._icase_vector = data['icase_vector']

        self._default_title = data['title']
        self._default_time = data['time']
        self._default_fps = data['frames/sec']
        self._default_resolution = data['resolution']

        self._scale = data['scale']
        self._default_scale = data['default_scale']
        self._default_is_scale = data['is_scale']

        self._arrow_scale = data['arrow_scale']
        self._default_arrow_scale = data['default_arrow_scale']

        self._phase = data['phase']
        self._default_phase = data['default_phase']

        self._default_dirname = data['dirname']
        self._default_gif_name = os.path.join(self._default_dirname, data['title'] + '.gif')

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
        self.gui = None
        if hasattr(self.win_parent, '_updated_legend'):
            self.win_parent.is_animate_open = True
            self.is_gui = True
            self.gui = self.win_parent.win_parent

        icase_max = 1000
        if is_gui_parent:
            self.is_gui = True
            self.gui = self.win_parent
            icase_max = max(self.gui.result_cases)  # TODO: update 1000

        self.icase_fringe_edit.setRange(0, icase_max)
        self.icase_disp_edit.setRange(1, icase_max)
        self.icase_vector_edit.setRange(1, icase_max)

        self.icase_fringe_start_edit.setRange(0, icase_max)
        self.icase_fringe_end_edit.setRange(0, icase_max)
        self.icase_fringe_delta_edit.setRange(1, icase_max)

        self.icase_disp_start_edit.setRange(0, icase_max)
        self.icase_disp_end_edit.setRange(0, icase_max)
        self.icase_disp_delta_edit.setRange(1, icase_max)


        self.on_update_min_max_defaults()


    def create_widgets(self):
        """creates the menu objects"""
        self.box_scale = QGroupBox('Animate Scale')
        self.box_time = QGroupBox('Animate Time')

        #icase_max = 1000

        self.checkbox_fringe = QCheckBox('Animate')
        self.checkbox_fringe.setToolTip(
            'Animate the fringe in addition to the deflection')
        #self.checkbox_disp = QCheckBox('Animate')
        self.checkbox_fringe.setEnabled(False)

        self.icase_fringe_label = QLabel("iCase (Fringe):")
        self.icase_fringe_edit = QSpinBox(self)
        #self.icase_fringe_edit.setRange(0, icase_max)
        self.icase_fringe_edit.setSingleStep(1)
        if self._icase_fringe is not None:
            self.icase_fringe_edit.setValue(self._icase_fringe)
        self.icase_fringe_edit.setToolTip(
            'Case Number for the Scale/Phase Animation Type.\n'
            'Defaults to the result you had shown when you clicked "Create Animation".\n'
            'iCase can be seen by clicking "Apply" on a result.')

        self.icase_disp_label = QLabel("iCase (Disp):")
        self.icase_disp_edit = QSpinBox(self)
        #self.icase_disp_edit.setRange(1, icase_max)
        self.icase_disp_edit.setSingleStep(1)
        if self._icase_disp is not None:
            self.icase_disp_edit.setValue(self._icase_disp)

        self.checkbox_vector = QCheckBox('Animate')
        self.checkbox_vector.setToolTip(
            'Animate the vector in addition to the deflection')
        self.checkbox_vector.hide()
        #self.checkbox_disp = QCheckBox('Animate')

        self.icase_vector_label = QLabel("iCase (Vector):")
        self.icase_vector_edit = QSpinBox(self)
        #self.icase_vector_edit.setRange(1, icase_max)
        self.icase_vector_edit.setSingleStep(1)
        if self._icase_vector is not None:
            self.icase_vector_edit.setValue(self._icase_vector)

        self.scale_label = QLabel("True Scale:")
        self.scale_edit = QLineEdit(str(self._scale))
        self.scale_button = QPushButton("Default")
        self.scale_edit.setToolTip('Scale factor of the "deflection"')
        self.scale_button.setToolTip('Sets the scale factor of the gif to %s' % self._scale)

        self.arrow_scale_label = QLabel("Arrow Scale:")
        self.arrow_scale_edit = QLineEdit(str(self._scale))
        self.arrow_scale_button = QPushButton("Default")
        self.arrow_scale_edit.setToolTip('Scale factor of the "arrows"')
        self.arrow_scale_button.setToolTip('Sets the arrow scale factor of the gif to %s' % (
            self._arrow_scale))

        self.arrow_scale_label.setVisible(False)
        self.arrow_scale_edit.setVisible(False)
        self.arrow_scale_button.setVisible(False)

        self.time_label = QLabel("Total Time (sec):")
        self.time_edit = QDoubleSpinBox(self)
        self.time_edit.setValue(self._default_time)
        self.time_edit.setRange(0.1, 5. * 60.)
        self.time_edit.setDecimals(2)
        self.time_edit.setSingleStep(0.1)
        self.time_button = QPushButton("Default")
        self.time_edit.setToolTip("Total time of the gif")
        self.time_button.setToolTip('Sets the total time of the gif to %.2f' % self._default_time)

        self.fps_label = QLabel("Frames/Second:")
        self.fps_edit = QSpinBox(self)
        self.fps_edit.setRange(1, 60)
        self.fps_edit.setSingleStep(1)
        self.fps_edit.setValue(self._default_fps)
        self.fps_button = QPushButton("Default")
        self.fps_edit.setToolTip("A higher FPS is smoother, but may not play well for large gifs")
        self.fps_button.setToolTip('Sets the FPS to %s' % self._default_fps)

        self.resolution_label = QLabel("Resolution Scale:")
        self.resolution_edit = QSpinBox(self)
        self.resolution_edit.setRange(1, 5)
        self.resolution_edit.setSingleStep(1)
        self.resolution_edit.setValue(self._default_resolution)
        self.resolution_button = QPushButton("Default")
        self.resolution_edit.setToolTip('Scales the window resolution by an integer factor')
        self.resolution_button.setToolTip('Sets the resolution to %s' % self._default_resolution)

        #-----------------
        # Time plot
        self.fringe_label = QLabel("Fringe")

        self.icase_fringe_start_edit = QSpinBox(self)
        #self.icase_fringe_start_edit.setRange(0, icase_max)
        self.icase_fringe_start_edit.setSingleStep(1)
        self.icase_fringe_start_edit.setValue(self._icase_fringe)
        self.icase_fringe_start_button = QPushButton("Default")

        self.icase_fringe_end_edit = QSpinBox(self)
        #self.icase_fringe_end_edit.setRange(0, icase_max)
        self.icase_fringe_end_edit.setSingleStep(1)
        self.icase_fringe_end_edit.setValue(self._icase_fringe)
        self.icase_fringe_end_button = QPushButton("Default")

        self.icase_fringe_delta_edit = QSpinBox(self)
        #self.icase_fringe_delta_edit.setRange(1, icase_max)
        self.icase_fringe_delta_edit.setSingleStep(1)
        self.icase_fringe_delta_edit.setValue(1)
        self.icase_fringe_delta_button = QPushButton("Default")

        self.displacement_label = QLabel("Displacement")
        self.icase_start = QLabel("iCase Start:")
        self.icase_disp_start_edit = QSpinBox(self)
        #self.icase_disp_start_edit.setRange(0, icase_max)
        self.icase_disp_start_edit.setSingleStep(1)
        self.icase_disp_start_edit.setValue(self._icase_fringe)
        self.icase_disp_start_button = QPushButton("Default")

        self.icase_end_label = QLabel("iCase End:")
        self.icase_disp_end_edit = QSpinBox(self)
        #self.icase_disp_end_edit.setRange(0, icase_max)
        self.icase_disp_end_edit.setSingleStep(1)
        self.icase_disp_end_edit.setValue(self._icase_fringe)
        self.icase_disp_end_button = QPushButton("Default")

        self.icase_delta_label = QLabel("iCase Delta:")
        self.icase_disp_delta_edit = QSpinBox(self)
        #self.icase_disp_delta_edit.setRange(1, icase_max)
        self.icase_disp_delta_edit.setSingleStep(1)
        self.icase_disp_delta_edit.setValue(1)
        self.icase_disp_delta_button = QPushButton("Default")

        self.min_value_enable = QCheckBox()
        self.min_value_label = QLabel("Min Value:")
        self.min_value_edit = QLineEdit('')
        #self.min_value_edit.setRange(1, 1000)
        #self.min_value_edit.setSingleStep(1)
        #self.min_value_edit.setValue(1)
        self.min_value_button = QPushButton("Default")

        self.max_value_enable = QCheckBox()
        self.max_value_label = QLabel("Max Value:")
        self.max_value_edit = QLineEdit('')
        #self.min_value_edit.setRange(1, 1000)  # TODO: update 1000
        #self.min_value_edit.setSingleStep(1)
        #self.min_value_edit.setValue(1)
        self.max_value_button = QPushButton("Default")

        # TODO: enable this (uncomment) ------------------------------------------
        #self.min_value_enable.hide()
        #self.min_value.hide()
        #self.min_value_edit.hide()
        #self.min_value_button.hide()

        #self.max_value_enable.hide()
        #self.max_value.hide()
        #self.max_value_edit.hide()
        #self.max_value_button.hide()
        # TODO: enable this (uncomment) ------------------------------------------

        self.icase_disp_start_edit.setToolTip('The first frame of the animation')
        self.icase_disp_end_edit.setToolTip(
            'The last frame of the animation\n'
            'Assumes icase_start + nframes * icase_delta = icase_end')
        self.icase_disp_delta_edit.setToolTip(
            'The frame step size (to skip non-consecutive results).\n'
            'Frame skipping can be used to:\n'
            "  - skip across results that you don't want to plot\n"
            '  - adjust the FPS')

        self.min_value_edit.setToolTip('Min value of the legend (not supported)')
        self.max_value_edit.setToolTip('Max value of the legend (not supported)')
        #'time' : 0.,
        #'default_time' : 0,
        #'icase_start' : 10,
        #'icase_delta' : 3,
        #'min_value' : 0.,
        #'max_value' : 1000.,

        self.browse_folder_label = QLabel('Output Directory:')
        self.browse_folder_edit = QLineEdit(str(self._default_dirname))
        self.browse_folder_button = QPushButton('Browse...')
        self.browse_folder_edit.setToolTip('Location to save the png/gif files')

        self.gif_label = QLabel("Gif Filename:")
        self.gif_edit = QLineEdit(str(self._default_title + '.gif'))
        self.gif_button = QPushButton('Default')
        self.gif_edit.setToolTip('Name of the gif')
        self.gif_button.setToolTip('Sets the name of the gif to %s.gif' % self._default_title)
        if not IS_IMAGEIO:
            self.gif_label.setEnabled(False)
            self.gif_edit.setEnabled(False)
            self.gif_button.setEnabled(False)

        # scale / phase
        if 1: # pragma: no cover
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
        else:
            msg = 'Scale : Animates the scale factor based on the "Animation Profile"\n'
            if self._default_phase is None:
                #self.animate_phase_radio.setDisabled(True)
                #self.animate_phase_radio.setToolTip('Animates the phase angle '
                                                    #'(only for complex results)')
                msg += 'Phase : Animates the phase angle (only for complex results)\n'
            else:
                #self.animate_phase_radio.setToolTip("Animates the phase angle")
                msg += 'Phase : Animates the phase angle\n'
            msg += (
                'Time : Animates the time/load/mode step\n'
                'Freq Sweep : Animates a complex result across a range of frequencies '
                '(not supported)\n'
            )

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

        self.animation_profile_label = QLabel("Animation Profile:")

        self.animation_profile_edit = QComboBox()
        for animation_profile in ANIMATION_PROFILES:
            self.animation_profile_edit.addItem(animation_profile)
        self.animation_profile_edit.setToolTip('The profile for a scaled GIF')

        self.animation_type_edit = QComboBox()
        # TODO: add a tooltip for each item
        for animation_type in self.animation_types:
            self.animation_type_edit.addItem(animation_type)
        #self.animation_type_edit.setToolTip('The profile for a scaled GIF')
        self.animation_type_edit.setToolTip(msg.rstrip())

        self.csv_profile_label = QLabel("CSV profile:")
        self.csv_profile_edit = QLineEdit()
        self.csv_profile_browse_button = QPushButton('Browse')
        self.csv_profile_edit.setToolTip(
            'The path to the CSV file of (Scale1, Scale2, Scale3, ...)')

        #widget = QWidget(self)
        #horizontal_vertical_group = QButtonGroup(widget)
        #horizontal_vertical_group.addButton(self.animate_scale_radio)
        #horizontal_vertical_group.addButton(self.animate_phase_radio)
        #horizontal_vertical_group.addButton(self.animate_time_radio)
        #horizontal_vertical_group.addButton(self.animate_freq_sweeep_radio)

        # animate in gui
        self.animate_in_gui_checkbox = QCheckBox("Animate In GUI?")
        self.animate_in_gui_checkbox.setChecked(True)

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
        self.wipe_button = QPushButton("Wipe Deformed Shape")
        self.stop_button = QPushButton("Stop")
        self.run_button = QPushButton("Run")

        self.step_button.setToolTip('Steps through the animation (for testing)')
        self.wipe_button.setToolTip('Removes the existing "deflecton" from the animation')
        self.stop_button.setToolTip('Stops the animation')
        self.run_button.setToolTip('Creates the animation')
        self.step_button.hide()
        self.wipe_button.hide()

        self.wipe_button.setEnabled(False)
        #self.wipe_button.hide()
        self.stop_button.setEnabled(False)

        self.cancel_button = QPushButton("Close")

        #self.set_grid_time(enabled=False)
        #self.set_grid_scale(enabled=self._default_is_scale)
        if self._default_phase:
            self.on_animate_phase(force=True)
            set_combo_box_text(self.animation_type_edit, 'Animate Phase')
        else:
            self.on_animate_scale(force=True)

    def set_connections(self):
        """creates the actions for the menu"""
        self.checkbox_vector.clicked.connect(self.on_checkbox_vector)

        self.scale_button.clicked.connect(self.on_default_scale)
        self.arrow_scale_button.clicked.connect(self.on_default_arrow_scale)
        self.time_button.clicked.connect(self.on_default_time)

        self.fps_button.clicked.connect(self.on_default_fps)
        self.resolution_button.clicked.connect(self.on_default_resolution)
        self.browse_folder_button.clicked.connect(self.on_browse_folder)
        self.csv_profile_browse_button.clicked.connect(self.on_browse_csv)
        self.gif_button.clicked.connect(self.on_default_title)

        self.step_button.clicked.connect(self.on_step)
        self.wipe_button.clicked.connect(self.on_wipe)
        self.stop_button.clicked.connect(self.on_stop)
        self.run_button.clicked.connect(self.on_run)
        self.min_value_enable.clicked.connect(self.on_min_value_enable)
        self.max_value_enable.clicked.connect(self.on_max_value_enable)

        self.min_value_button.clicked.connect(self.on_min_value_default)
        self.max_value_button.clicked.connect(self.on_max_value_default)
        self.icase_disp_start_button.clicked.connect(self.on_update_min_max_defaults)

        #self.animate_scale_radio.clicked.connect(self.on_animate_scale)
        #self.animate_phase_radio.clicked.connect(self.on_animate_phase)
        #self.animate_time_radio.clicked.connect(self.on_animate_time)
        self.animation_type_edit.currentIndexChanged.connect(self.on_animate)
        #self.animate_freq_sweeep_radio

        self.cancel_button.clicked.connect(self.on_cancel)

        self.animate_in_gui_checkbox.clicked.connect(self.on_animate_in_gui)
        self.animate_in_gui_checkbox.setChecked(True)
        self.on_animate_in_gui()

    def on_checkbox_vector(self):
        is_enabled = self.checkbox_vector.isEnabled()
        is_checked = self.checkbox_vector.isChecked()
        enable_edit = is_enabled and is_checked
        self.icase_vector_label.setEnabled(is_checked)
        self.icase_vector_edit.setEnabled(enable_edit)

    def on_animate_in_gui(self):
        animate_in_gui = self.animate_in_gui_checkbox.isChecked()
        enable = not animate_in_gui
        if HIDE_WHEN_INACTIVE:
            self.make_images_checkbox.setVisible(enable)
            self.delete_images_checkbox.setVisible(enable)
            self.make_gif_checkbox.setVisible(enable)
            self.repeat_checkbox.setVisible(enable)
            self.resolution_button.setVisible(enable)
            self.resolution_label.setVisible(enable)
            self.resolution_edit.setVisible(enable)
            self.gif_label.setVisible(enable)
            self.gif_edit.setVisible(enable)
            self.gif_button.setVisible(enable)
            self.browse_folder_label.setVisible(enable)
            self.browse_folder_button.setVisible(enable)
            self.browse_folder_edit.setVisible(enable)
            self.step_button.setEnabled(enable)

        self.make_images_checkbox.setEnabled(enable)
        self.delete_images_checkbox.setEnabled(enable)
        self.repeat_checkbox.setEnabled(enable)
        self.resolution_button.setEnabled(enable)
        self.resolution_edit.setEnabled(enable)
        if IS_IMAGEIO:
            self.make_gif_checkbox.setEnabled(enable)
            self.gif_edit.setEnabled(enable)
            self.gif_button.setEnabled(enable)
        self.browse_folder_button.setEnabled(enable)
        self.browse_folder_edit.setEnabled(enable)
        self.step_button.setEnabled(enable)
        #wipe_button

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
        if HIDE_WHEN_INACTIVE:
            self.box_scale.setVisible(enabled)
            self.animation_profile_label.setVisible(enabled)
            self.animation_profile_edit.setVisible(enabled)
            #self.min_value_enable.setVisible(enabled)
            #self.max_value_enable.setVisible(enabled)

        self.animation_profile_label.setEnabled(enabled)
        self.animation_profile_edit.setEnabled(enabled)

        # TODO: doesn't work...
        #self.csv_profile.setEnabled(enabled)
        #self.csv_profile_edit.setEnabled(enabled)
        #self.csv_profile_button.setEnabled(enabled)

        self.min_value_enable.setEnabled(enabled)
        self.max_value_enable.setEnabled(enabled)
        self.on_min_value_enable()
        self.on_max_value_enable()


    def set_grid_time(self, enabled=True, word=''):
        """enables/disables the secondary input"""
        #print('%s-set_grid_time; enabled = %r' % (word, enabled))
        if HIDE_WHEN_INACTIVE:
            self.box_time.setVisible(enabled)
            self.checkbox_fringe.setVisible(not enabled)
            self.icase_fringe_label.setVisible(not enabled)
            self.icase_fringe_edit.setVisible(not enabled)
            self.icase_disp_label.setVisible(not enabled)
            self.icase_disp_edit.setVisible(not enabled)
            self.icase_vector_label.setVisible(not enabled)
            self.icase_vector_edit.setVisible(not enabled)

            #self.icase_fringe_delta_edit.setVisible(enabled)
            self.icase_disp_delta_edit.setVisible(enabled)
            self.icase_disp_delta_edit.setVisible(enabled)
            self.fps_label.setVisible(enabled)
            self.fps_edit.setVisible(enabled)
            self.fps_button.setVisible(enabled)

        self.displacement_label.setEnabled(enabled)
        self.fringe_label.setEnabled(enabled)

        self.icase_start.setEnabled(enabled)
        self.icase_disp_start_edit.setEnabled(enabled)
        self.icase_disp_start_button.setEnabled(enabled)
        self.icase_fringe_start_edit.setEnabled(enabled)
        self.icase_fringe_start_button.setEnabled(enabled)

        self.icase_end_label.setEnabled(enabled)
        self.icase_disp_end_edit.setEnabled(enabled)
        self.icase_disp_end_button.setEnabled(enabled)
        self.icase_fringe_end_edit.setEnabled(enabled)
        self.icase_fringe_end_button.setEnabled(enabled)

        self.icase_delta_label.setEnabled(enabled)
        self.icase_disp_delta_edit.setEnabled(enabled)
        self.icase_disp_delta_button.setEnabled(enabled)
        self.icase_fringe_delta_edit.setEnabled(enabled)
        self.icase_fringe_delta_button.setEnabled(enabled)

        #-----------------------------------------------------------------------
        is_min_enabled = self.min_value_enable.isChecked()
        self.min_value_label.setEnabled(is_min_enabled)
        self.min_value_edit.setEnabled(is_min_enabled)
        self.min_value_button.setEnabled(is_min_enabled)

        is_max_enabled = self.max_value_enable.isChecked()
        self.max_value_label.setEnabled(is_max_enabled)
        self.max_value_edit.setEnabled(is_max_enabled)
        self.max_value_button.setEnabled(is_max_enabled)

        self.min_value_enable.setEnabled(enabled)
        self.on_min_value_enable()
        #self.min_value.setEnabled(enabled)
        #self.min_value_edit.setEnabled(enabled)
        #self.min_value_button.setEnabled(enabled)

        self.max_value_enable.setEnabled(enabled)
        self.on_max_value_enable()
        #self.max_value.setEnabled(enabled)
        #self.max_value_edit.setEnabled(enabled)
        #self.max_value_button.setEnabled(enabled)

        self.icase_fringe_label.setEnabled(not enabled)
        self.icase_fringe_edit.setEnabled(not enabled)
        self.checkbox_fringe.setEnabled(not enabled)

        self.icase_disp_label.setEnabled(not enabled)
        self.icase_disp_edit.setEnabled(not enabled)

        self.icase_vector_label.setEnabled(not enabled)
        self.icase_vector_edit.setEnabled(not enabled)
        self.checkbox_vector.setEnabled(not enabled)
        self.on_checkbox_vector()

        self.fps_label.setEnabled(not enabled)
        self.fps_edit.setEnabled(not enabled)
        self.fps_button.setEnabled(not enabled)

    def on_min_value_enable(self):
        """
        The min edit value box is enabled when we switch to time
        and the box is checked
        """
        is_min_enabled = self.min_value_enable.isChecked() and self.min_value_enable.isEnabled()
        self.min_value_label.setEnabled(is_min_enabled)
        self.min_value_edit.setEnabled(is_min_enabled)
        self.min_value_button.setEnabled(is_min_enabled)

    def on_max_value_enable(self):
        """
        The max edit value box is enabled when we switch to time
        and the box is checked
        """
        is_max_enabled = self.max_value_enable.isChecked() and self.max_value_enable.isEnabled()
        self.max_value_label.setEnabled(is_max_enabled)
        self.max_value_edit.setEnabled(is_max_enabled)
        self.max_value_button.setEnabled(is_max_enabled)

    def on_update_min_max_defaults(self):
        """
        When the icase is changed, the min/max value default message is changed
        """
        icase = self.icase_disp_start_edit.value()
        min_value, max_value = self.get_min_max(icase)
        self.min_value_button.setToolTip('Sets the min value to %g' % min_value)
        self.max_value_button.setToolTip('Sets the max value to %g' % max_value)

    def on_min_value_default(self):
        """When min default icase is pressued, update the value"""
        icase = self.icase_disp_start_edit.value()
        min_value = self.get_min_max(icase)[0]
        self.min_value_edit.setText(str(min_value))
        self.min_value_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_max_value_default(self):
        """When max default icase is pressued, update the value"""
        icase = self.icase_disp_start_edit.value()
        max_value = self.get_min_max(icase)[1]
        self.max_value_edit.setText(str(max_value))
        self.max_value_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_browse_folder(self):
        """opens a folder dialog"""
        dirname = getexistingdirectory(parent=self, caption='Select a Directory',
                                       basedir='',
                                       options=QFileDialog.ShowDirsOnly)
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

    def on_default_title(self):
        """sets the default gif name"""
        self.gif_edit.setText(self._default_title + '.gif')

    def on_default_scale(self):
        """sets the default displacement scale factor"""
        if self.is_gui:
            icase_disp = self.icase_disp_edit.value()
            out = self.gui.legend_obj.get_legend_disp(
                icase_disp)
            unused_scale, unused_phase, default_scale, unused_default_phase = out
        else:
            default_scale = self._default_scale
        self.scale_edit.setText(str(default_scale))
        self.scale_edit.setStyleSheet("QLineEdit{background: white;}")

    def on_default_arrow_scale(self):
        """sets the default arrow scale factor"""
        if self.is_gui:
            icase_vector = self.icase_vector_edit.value()
            out = self.gui.legend_obj.get_legend_vector(icase_vector)
            unused_arrow_scale, default_arrow_scale = out
        else:
            default_arrow_scale = self._default_arrow_scale
        self.arrow_scale_edit.setText(str(default_arrow_scale))
        self.arrow_scale_edit.setStyleSheet("QLineEdit{background: white;}")

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
        irow = 0
        grid.addWidget(self.icase_fringe_label, irow, 0)
        grid.addWidget(self.icase_fringe_edit, irow, 1)
        grid.addWidget(self.checkbox_fringe, irow, 2)
        irow += 1

        grid.addWidget(self.icase_disp_label, irow, 0)
        grid.addWidget(self.icase_disp_edit, irow, 1)
        #grid.addWidget(self.checkbox_disp, irow, 2)
        irow += 1

        grid.addWidget(self.icase_vector_label, irow, 0)
        grid.addWidget(self.icase_vector_edit, irow, 1)
        grid.addWidget(self.checkbox_vector, irow, 2)
        irow += 1

        grid.addWidget(self.scale_label, irow, 0)
        grid.addWidget(self.scale_edit, irow, 1)
        grid.addWidget(self.scale_button, irow, 2)
        irow += 1

        grid.addWidget(self.arrow_scale_label, irow, 0)
        grid.addWidget(self.arrow_scale_edit, irow, 1)
        grid.addWidget(self.arrow_scale_button, irow, 2)
        irow += 1

        grid.addWidget(self.time_label, irow, 0)
        grid.addWidget(self.time_edit, irow, 1)
        grid.addWidget(self.time_button, irow, 2)
        irow += 1

        # spacer
        spacer = QLabel('')

        grid.addWidget(self.fps_label, irow, 0)
        grid.addWidget(self.fps_edit, irow, 1)
        grid.addWidget(self.fps_button, irow, 2)
        irow += 1


        grid.addWidget(self.animation_type, irow, 0)
        grid.addWidget(self.animation_type_edit, irow, 1)
        irow += 1

        grid.addWidget(spacer, irow, 0)
        irow += 1

        #----------
        #Time
        grid_time = QGridLayout()
        jrow = 0

        self.fringe_label.setAlignment(Qt.AlignCenter)
        self.displacement_label.setAlignment(Qt.AlignCenter)

        if not IS_TIME_FRINGE:
            self.fringe_label.hide()
            self.icase_fringe_delta_edit.hide()
            self.icase_fringe_start_edit.hide()
            self.icase_fringe_end_edit.hide()
            self.icase_fringe_delta_button.hide()

        grid_time.addWidget(self.displacement_label, jrow, 1)
        grid_time.addWidget(self.fringe_label, jrow, 2)
        jrow += 1

        grid_time.addWidget(self.icase_start, jrow, 0)
        grid_time.addWidget(self.icase_disp_start_edit, jrow, 1)
        grid_time.addWidget(self.icase_fringe_start_edit, jrow, 2)
        #grid_time.addWidget(self.icase_disp_start_button, jrow, 2)
        jrow += 1

        grid_time.addWidget(self.icase_end_label, jrow, 0)
        grid_time.addWidget(self.icase_disp_end_edit, jrow, 1)
        grid_time.addWidget(self.icase_fringe_end_edit, jrow, 2)
        #grid_time.addWidget(self.icase_end_button, jrow, 2)
        jrow += 1

        grid_time.addWidget(self.icase_delta_label, jrow, 0)
        grid_time.addWidget(self.icase_disp_delta_edit, jrow, 1)
        grid_time.addWidget(self.icase_fringe_delta_edit, jrow, 2)
        #grid_time.addWidget(self.icase_delta_button, jrow, 2)
        jrow += 1


        hbox_min = QHBoxLayout()
        hbox_min.addWidget(self.min_value_enable)
        hbox_min.addWidget(self.min_value_label)
        grid_time.addLayout(hbox_min, jrow, 0)
        grid_time.addWidget(self.min_value_edit, jrow, 1)
        grid_time.addWidget(self.min_value_button, jrow, 2)
        jrow += 1

        hbox_max = QHBoxLayout()
        hbox_max.addWidget(self.max_value_enable)
        hbox_max.addWidget(self.max_value_label)
        grid_time.addLayout(hbox_max, jrow, 0)
        grid_time.addWidget(self.max_value_edit, jrow, 1)
        grid_time.addWidget(self.max_value_button, jrow, 2)
        jrow += 1

        grid_time.addWidget(spacer, jrow, 0)
        jrow += 1

        #--------------
        grid_scale = QGridLayout()
        grid_scale.addWidget(self.animation_profile_label, 0, 0)
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
        self.box_scale.setLayout(grid_scale)

        self.box_time.setLayout(grid_time)
        #----------

        grid2 = QGridLayout()
        irow = 0
        #grid2.addWidget(self.animate_scale_radio, 8, 0)
        #grid2.addWidget(self.animate_phase_radio, 8, 1)
        #grid2.addWidget(self.animate_time_radio, 8, 2)
        #grid2.addWidget(self.animate_freq_sweeep_radio, 8, 3)

        grid2.addWidget(self.animate_in_gui_checkbox, irow, 0)
        irow += 1

        grid2.addWidget(self.resolution_label, irow, 0)
        grid2.addWidget(self.resolution_edit, irow, 1)
        grid2.addWidget(self.resolution_button, irow, 2)
        irow += 1

        grid2.addWidget(self.browse_folder_label, irow, 0)
        grid2.addWidget(self.browse_folder_edit, irow, 1)
        grid2.addWidget(self.browse_folder_button, irow, 2)
        irow += 1

        grid2.addWidget(self.gif_label, irow, 0)
        grid2.addWidget(self.gif_edit, irow, 1)
        grid2.addWidget(self.gif_button, irow, 2)
        irow += 1

        grid2.addWidget(self.make_images_checkbox, irow, 0)
        #grid2.addWidget(self.overwrite_images_checkbox, irow, 0)
        grid2.addWidget(self.delete_images_checkbox, irow, 1)
        grid2.addWidget(self.make_gif_checkbox, irow, 2)
        irow += 1
        grid2.addWidget(self.repeat_checkbox, irow, 0)
        irow += 1
        grid2.addWidget(spacer, irow, 0)

        grid_hbox = QHBoxLayout()
        grid_hbox.addWidget(spacer)
        grid_hbox.addLayout(grid2)
        grid_hbox.addWidget(spacer)

        # bottom buttons
        step_run_box = QHBoxLayout()
        step_run_box.addWidget(self.step_button)
        step_run_box.addWidget(self.wipe_button)
        step_run_box.addWidget(self.stop_button)
        step_run_box.addWidget(self.run_button)

        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.cancel_button)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        vbox.addWidget(self.box_scale)
        vbox.addWidget(self.box_time)
        #vbox.addLayout(checkboxes)
        vbox.addLayout(grid_hbox)
        vbox.addStretch()
        vbox.addLayout(step_run_box)
        vbox.addLayout(ok_cancel_box)

        if IS_RESULTS_SELECTOR and self.fringe_cases:
            cases = get_cases_from_tree(self.fringe_cases)
            parent = self
            name = 'main'
            data = self.fringe_cases
            choices = cases
            results_widget = ResultsWindow(parent, name, data, choices,
                                           include_clear=False, include_delete=False)
            vbox_results = QVBoxLayout()
            results_widget_label = QLabel('Results:')
            vbox_results.addWidget(results_widget_label)
            vbox_results.addWidget(results_widget)
            hbox_main = QHBoxLayout()
            hbox_main.addLayout(vbox)
            hbox_main.addLayout(vbox_results)
            self.setLayout(hbox_main)
        else:
            self.setLayout(vbox)

    def on_fringe(self, icase):
        """sets the icase fringe"""
        self.icase_fringe_edit.setValue(icase)

    def on_disp(self, icase):
        """sets the icase disp"""
        self.icase_disp_edit.setValue(icase)

    def on_vector(self, icase):
        """sets the icase vector"""
        self.icase_vector_edit.setValue(icase)

    def on_clear_results(self):
        """sink for the right click menu"""
        pass

    def on_step(self):
        """click the Step button"""
        passed, validate_out = self.on_validate()
        if passed:
            try:
                self._make_gif(validate_out, istep=self.istep)
                self.istep += 1
            except IndexError:
                self._make_gif(validate_out, istep=0)
                self.istep += 1
            self.wipe_button.setEnabled(True)

    def on_wipe(self):
        """click the Wipe button"""
        passed, validate_out = self.on_validate(wipe=True)
        if passed:
            self.istep = 0
            self._make_gif(validate_out, istep=self.istep)
            self.wipe_button.setEnabled(False)
            self.stop_button.setEnabled(False)

    def on_stop(self):
        """click the Stop button"""
        #passed, validate_out = self.on_validate()
        #if passed:
            #self._make_gif(validate_out, stop_animation=True)
        if self.is_gui:
            self.gui.stop_animation()

        self.wipe_button.setEnabled(True)
        self.stop_button.setEnabled(False)


    def on_run(self):
        """click the Run button"""
        self.istep = 0
        self.wipe_button.setEnabled(False)
        self.stop_button.setEnabled(True)

        passed, validate_out = self.on_validate()
        if passed:
            self._make_gif(validate_out, istep=None)
        return passed

    def _make_gif(self, validate_out, istep=None, stop_animation=False):
        """interface for making the gif"""
        (icase_fringe, icase_disp, icase_vector, scale, time, fps, animate_in_gui,
         magnify, output_dir, gifbase,
         min_value, max_value) = validate_out
        fps = int(fps)

        gif_filename = None
        if not stop_animation and not animate_in_gui and gifbase is not None:
            if gifbase.lower().endswith('.gif'):
                gifbase = gifbase[:-4]
            gif_filename = os.path.join(output_dir, gifbase + '.gif')

        animate_fringe = self.checkbox_fringe.isChecked()
        animate_vector = self.checkbox_vector.isChecked()

        animate_scale = self.animate_scale_radio.isChecked()
        animate_phase = self.animate_phase_radio.isChecked()
        animate_time = self.animate_time_radio.isChecked()

        if not self.checkbox_vector.isEnabled():
            icase_vector = None

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
        animation_profile = str(self.animation_profile_edit.currentText())

        icase_disp_start = self.icase_disp_start_edit.value()
        icase_disp_end = self.icase_disp_end_edit.value()
        icase_disp_delta = self.icase_disp_delta_edit.value()

        bool_repeat = self.repeat_checkbox.isChecked()  # TODO: change this to an integer
        if bool_repeat:
            nrepeat = 0
        else:
            nrepeat = 1
        #self.out_data['is_shown'] = self.show_radio.isChecked()
        #icase = self._icase

        if self.is_gui:
            self.gui.make_gif(
                gif_filename, scale, istep=istep,
                animate_scale=animate_scale, animate_phase=animate_phase, animate_time=animate_time,
                icase_fringe=icase_fringe, icase_disp=icase_disp, icase_vector=icase_vector,
                animate_fringe=animate_fringe, animate_vector=animate_vector,
                icase_start=icase_disp_start, icase_end=icase_disp_end, icase_delta=icase_disp_delta,
                time=time, animation_profile=animation_profile,
                nrepeat=nrepeat, fps=fps, magnify=magnify,
                make_images=make_images, delete_images=delete_images, make_gif=make_gif,
                stop_animation=stop_animation, animate_in_gui=animate_in_gui,
                min_value=min_value, max_value=max_value,
            )

        self.out_data['clicked_ok'] = True
        self.out_data['close'] = True

    def get_min_max(self, icase):
        if self.is_gui:
            (obj, (i, name)) = self.gui.result_cases[icase]
            min_value, max_value = obj.get_min_max(i, name)
        else:
            return 0., 1.0
        return min_value, max_value

    def on_validate(self, wipe=False):
        """checks to see if the input is valid"""
        # requires no special validation
        icase_fringe, flag0 = check_int(self.icase_fringe_edit)
        icase_disp, unused_flaga = check_int(self.icase_disp_edit)
        icase_vector, unused_flagb = check_int(self.icase_vector_edit)
        #icase_disp = self._icase_disp
        #icase_vector = self._icase_vector

        scale, flag1 = check_float(self.scale_edit)
        time, flag2 = check_float(self.time_edit)
        fps, flag3 = check_float(self.fps_edit)

        min_value = max_value = None
        flag4 = flag5 = True
        if self.min_value_edit.isEnabled():
            min_value, flag4 = check_float(self.min_value_edit)
        if self.max_value_edit.isEnabled():
            max_value, flag5 = check_float(self.max_value_edit)

        if wipe:
            animate_in_gui = False
            scale = 0.
            flag1 = True
        else:
            animate_in_gui = self.animate_in_gui_checkbox.isChecked()
            if scale == 0.0:
                self.scale_edit.setStyleSheet("QLineEdit{background: red;}")
                flag1 = False

        if animate_in_gui or wipe:
            passed = all([flag0, flag1, flag2, flag3, flag4, flag5])
            magnify, output_dir, gifbase = None, None, None
        else:
            magnify, flag6 = check_int(self.resolution_edit)
            output_dir, flag7 = check_path(self.browse_folder_edit)
            gifbase, flag8 = check_name_str(self.gif_edit)
            passed = all([flag0, flag1, flag2, flag3, flag4, flag5, flag6, flag7, flag8])
        return passed, (icase_fringe, icase_disp, icase_vector, scale, time, fps, animate_in_gui,
                        magnify, output_dir, gifbase, min_value, max_value)

    #def on_ok(self):
        #"""click the OK button"""
        #passed = self.on_apply()
        #if passed:
            #self.win_parent._animation_window_shown = False
            #self.close()
            ##self.destroy()

    def on_cancel(self):
        """click the Cancel button"""
        self.on_stop()
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

    #from pyNastran.gui.menus.legend.animation import AnimationWindow
    data2 = {
        'font_size' : 8,
        'icase_fringe' : 1,
        'icase_disp' : 2,
        'icase_vector' : 3,

        'title' : 'cat',
        'time' : 2,
        'frames/sec' : 30,
        'resolution' : 1,
        'iframe' : 0,
        'is_scale' : False,
        'dirname' : os.getcwd(),
        'scale' : 2.0,
        'default_scale' : 10,

        'arrow_scale' : 3.0,
        'default_arrow_scale' : 30,

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
    data2['phase'] = 0.  # uncomment for phase

    form = [
        [u'Geometry', None, [
            (u'NodeID', 0, []),
            (u'ElementID', 1, []),
            (u'PropertyID', 2, []),
            (u'MaterialID', 3, []),
            (u'E', 4, []),
            (u'Element Checks', None, [
                (u'ElementDim', 5, []),
                (u'Min Edge Length', 6, []),
                (u'Min Interior Angle', 7, []),
                (u'Max Interior Angle', 8, [])],
            ),],
        ],
    ]
    #[0, 1, 2, 3, 4, 5, 6, 7, 8]
    main_window = AnimationWindow(data2, fringe_cases=form)
    main_window.show()
    # Enter the main loop
    app.exec_()


if __name__ == "__main__": # pragma: no cover
    main()
