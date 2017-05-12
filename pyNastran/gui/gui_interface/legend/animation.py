"""
defines:
 - LegendPropertiesWindow
"""
from __future__ import print_function
import os
from six import integer_types
import numpy as np

from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4.QtGui import (
        QApplication, QLabel, QPushButton, QLineEdit, QWidget, QRadioButton,
        QButtonGroup, QGridLayout, QHBoxLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox, QCheckBox)
elif qt_version == 5:
    from PyQt5.QtWidgets import (
        QApplication, QLabel, QPushButton, QLineEdit, QWidget, QRadioButton,
        QButtonGroup, QGridLayout, QHBoxLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox, QCheckBox)
elif qt_version == 'pyside':
    from PySide.QtGui import (
        QApplication, QLabel, QPushButton, QLineEdit, QWidget, QRadioButton,
        QButtonGroup, QGridLayout, QHBoxLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox, QCheckBox)
else:
    raise NotImplementedError('qt_version = %r' % qt_version)

from pyNastran.gui.gui_interface.common import PyDialog
from pyNastran.gui.gui_utils.dialogs import open_directory_dialog
from pyNastran.gui.gui_utils.write_gif import IS_IMAGEIO


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

        self.setWindowTitle('Animate Model')
        self.create_widgets()
        self.create_layout()
        self.set_connections()
        if hasattr(self.win_parent, 'make_gif'):
            self.win_parent.is_animate_open = True

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
        self.icase_start = QLabel("icase_start:")
        self.icase_start_edit = QSpinBox(self)
        self.icase_start_edit.setRange(1, 5)
        self.icase_start_edit.setSingleStep(1)
        self.icase_start_edit.setValue(self._icase)
        self.icase_start_button = QPushButton("Default")

        self.icase_delta = QLabel("icase_delta:")
        self.icase_delta_edit = QSpinBox(self)
        self.icase_delta_edit.setRange(1, 1000)  # TODO: update 1000
        self.icase_delta_edit.setSingleStep(1)
        self.icase_delta_edit.setValue(1)
        self.icase_delta_button = QPushButton("Default")

        self.min_value = QLabel("Min Value:")
        self.min_value_edit = QLineEdit(str(0.))
        #self.min_value_edit.setRange(1, 1000)  # TODO: update 1000
        #self.min_value_edit.setSingleStep(1)
        #self.min_value_edit.setValue(1)
        self.min_value_button = QPushButton("Default")

        self.max_value = QLabel("Max Value:")
        self.max_value_edit = QLineEdit(str(1.))
        #self.min_value_edit.setRange(1, 1000)  # TODO: update 1000
        #self.min_value_edit.setSingleStep(1)
        #self.min_value_edit.setValue(1)
        self.max_value_button = QPushButton("Default")

        #'time' : 0.,
        #'default_time' : 0,
        #'icase_start' : 10,
        #'icase_delta' : 3,
        #'min_value' : 0.,
        #'max_value' : 1000.,


        #self.browse = QLabel("Animation File:")
        self.browse = QLabel('Output Directory:')
        self.browse_edit = QLineEdit(str(self._default_dirname))
        self.browse_button = QPushButton('Browse')
        self.browse_edit.setToolTip('Location to save the png/gif files')

        self.gif = QLabel("Gif Filename:")
        self.gif_edit = QLineEdit(str(self._default_name + '.gif'))
        self.gif_button = QPushButton('Default')
        self.gif_edit.setToolTip('Name of the gif')
        self.gif_button.setToolTip('Sets the name of the gif to %s.gif' % self._default_name)

        # scale / phase
        self.animate_scale_radio = QRadioButton("Animate Scale")
        self.animate_phase_radio = QRadioButton("Animate Phase")
        self.animate_time_radio = QRadioButton("Animate Time")
        self.animate_scale_radio.setToolTip("Animates the scale factor")
        self.animate_time_radio.setToolTip('Animates the time/load/mode step (not supported)')

        self.animate_scale_radio.setChecked(self._default_is_scale)
        self.animate_phase_radio.setChecked(not self._default_is_scale)
        self.animate_time_radio.setChecked(False)
        if self._default_phase is None:
            self.animate_phase_radio.setDisabled(True)
            self.animate_phase_radio.setToolTip('Animates the phase angle '
                                                '(only for complex results)')
        else:
            self.animate_phase_radio.setToolTip("Animates the phase angle")

        self.animate_time_radio.setDisabled(True)
        widget = QWidget(self)
        horizontal_vertical_group = QButtonGroup(widget)
        horizontal_vertical_group.addButton(self.animate_scale_radio)
        horizontal_vertical_group.addButton(self.animate_phase_radio)
        horizontal_vertical_group.addButton(self.animate_time_radio)

        # one / two sided
        self.onesided_radio = QRadioButton("One Sided")
        self.onesided_radio.setToolTip("A one sided gif doesn't return to the starting point "
                                       "(e.g., 0 to 360 degrees)")
        self.twosided_radio = QRadioButton("Two Sided")
        self.twosided_radio.setToolTip('A two sided gif returns to the starting point '
                                       '(e.g., 0 to 10 to 0)')

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

    def set_connections(self):
        """creates button actions"""
        self.scale_button.clicked.connect(self.on_default_scale)
        self.time_button.clicked.connect(self.on_default_time)

        self.fps_button.clicked.connect(self.on_default_fps)
        self.resolution_button.clicked.connect(self.on_default_resolution)
        self.browse_button.clicked.connect(self.on_browse)
        self.gif_button.clicked.connect(self.on_default_name)

        self.step_button.clicked.connect(self.on_step)
        self.run_button.clicked.connect(self.on_run)

        self.animate_scale_radio.clicked.connect(self.on_animate_time)
        self.animate_phase_radio.clicked.connect(self.on_animate_scale_phase)
        self.animate_time_radio.clicked.connect(self.on_animate_time)

        #self.apply_button.clicked.connect(self.on_apply)
        #self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)

    def on_animate_time(self):
        """enables the secondary input"""
        self.secondary_input.setEnabled(True)

    def on_animate_scale_phase(self):
        """enables the secondary input"""
        self.secondary_input.setEnabled(False)

    def on_browse(self):
        """opens a folder dialog"""
        dirname = open_directory_dialog(self, 'Select a Directory')
        if not dirname:
            return
        self.browse_edit.setText(dirname)

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
        #grid.addWidget(spacer, 2, 0)

        grid.addWidget(self.fps, 3, 0)
        grid.addWidget(self.fps_edit, 3, 1)
        grid.addWidget(self.fps_button, 3, 2)

        grid.addWidget(self.resolution, 4, 0)
        grid.addWidget(self.resolution_edit, 4, 1)
        grid.addWidget(self.resolution_button, 4, 2)

        grid.addWidget(self.browse, 5, 0)
        grid.addWidget(self.browse_edit, 5, 1)
        grid.addWidget(self.browse_button, 5, 2)

        grid.addWidget(self.gif, 6, 0)
        grid.addWidget(self.gif_edit, 6, 1)
        grid.addWidget(self.gif_button, 6, 2)

        #----------
        #Time
        if 0:
            grid.addWidget(self.icase_start, 7, 0)
            grid.addWidget(self.icase_start_edit, 7, 1)
            grid.addWidget(self.icase_start_button, 7, 2)

            grid.addWidget(self.icase_delta, 8, 0)
            grid.addWidget(self.icase_delta_edit, 8, 1)
            grid.addWidget(self.icase_delta_button, 8, 2)

            grid.addWidget(self.min_value, 9, 0)
            grid.addWidget(self.min_value_edit, 9, 1)
            grid.addWidget(self.min_value_button, 9, 2)

            grid.addWidget(self.max_value, 10, 0)
            grid.addWidget(self.max_value_edit, 10, 1)
            grid.addWidget(self.max_value_button, 10, 2)

        #----------

        grid.addWidget(spacer, 7, 0)

        grid2 = QGridLayout()
        grid2.addWidget(self.animate_scale_radio, 8, 0)
        grid2.addWidget(self.animate_phase_radio, 8, 1)
        grid2.addWidget(self.animate_time_radio, 8, 2)

        grid2.addWidget(self.twosided_radio, 9, 0)
        grid2.addWidget(self.onesided_radio, 9, 1)

        grid2.addWidget(self.repeat_checkbox, 10, 0)
        grid2.addWidget(self.delete_images_checkbox, 10, 1)
        grid2.addWidget(self.make_gif_checkbox, 10, 2)

        grid2.addWidget(spacer, 11, 0)
        grid_hbox = QHBoxLayout()
        grid_hbox.addWidget(spacer)
        grid_hbox.addLayout(grid2)
        grid_hbox.addWidget(spacer)

        #grid.addWidget(self.scale_radio, 6, 0)
        #grid.addWidget(self.phase_radio, 6, 1)
        #grid.addWidget(self.delete_images_checkbox, 6, 0)

        # bottom buttons
        step_run_box = QHBoxLayout()
        step_run_box.addWidget(self.step_button)
        step_run_box.addWidget(self.run_button)

        ok_cancel_box = QHBoxLayout()
        #ok_cancel_box.addWidget(self.apply_button)
        #ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
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

        make_images = True
        delete_images = self.delete_images_checkbox.isChecked()
        make_gif = self.make_gif_checkbox.isChecked()

        onesided = self.onesided_radio.isChecked()
        bool_repeat = self.repeat_checkbox.isChecked()  # TODO: change this to an integer
        if bool_repeat:
            nrepeat = 0
        else:
            nrepeat = 1
        #self.out_data['is_shown'] = self.show_radio.isChecked()
        icase = self._icase
        self.win_parent.win_parent.make_gif(
            gif_filename, icase, scale, istep=istep, onesided=onesided,
            animate_scale=animate_scale, animate_phase=animate_phase, animate_time=animate_time,
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
        output_dir, flag4 = self.check_path(self.browse_edit)
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

def get_analysis_time(time, onesided):
    """
    The analysis time is the time that needs to be simulated for the analysis.

    TODO: could we define time as 1/2-sided time so we can do less work?
    TODO: we could be more accurate regarding dt
          Nonesided = 5
          Ntwosided = 2 * Nonesided - 1 = 9
          Nonesided = (Ntwosided + 1) / 2

          Nframes = int(fps * t)
          Nonesided = Nframes
          Ntwosided = 2 * Nonesided - 1 = 9
          Nonesided = (Ntwosided + 1) / 2
    """
    if onesided:
        analysis_time = time / 2.
    else:
        analysis_time = time
    return analysis_time

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
        'phase' : 0.,
        'default_phase' : 120.,

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
