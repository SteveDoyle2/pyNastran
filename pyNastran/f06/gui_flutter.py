import os
import sys
import warnings
from functools import wraps
from pathlib import Path
from typing import Callable, Optional, Any

from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QRadioButton
ICON_PATH = Path('')
try:
    import json5 as json
except ModuleNotFoundError:
    warnings.warn('couldnt find json5, using json')
    import json

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

#from PyQt5.QtGui import QKeySequence
from qtpy import QtCore
from qtpy.compat import getopenfilename, getsavefilename
from qtpy.QtWidgets import (
    QLabel, QWidget,
    QApplication, QMainWindow, QMenu, QVBoxLayout, QLineEdit, QComboBox,
    QHBoxLayout, QPushButton, QGridLayout, QFileDialog,
    QAction, QTableWidget, QTableWidgetItem,
    QCheckBox, QListWidgetItem, QAbstractItemView,
    #QMessageBox, QFrame, QMenuBar, QAbstractItemView, QListView,
    QListWidget, QDockWidget,
)
from qtpy.QtGui import QColor, QIcon, QPixmap

from cpylog import SimpleLogger
from cpylog.html_utils import str_to_html
import pyNastran
from pyNastran.gui.menus.application_log import ApplicationLogWidget
from pyNastran.f06.flutter_response import FlutterResponse
from pyNastran.f06.parse_flutter import make_flutter_response

PKG_PATH = os.path.dirname(pyNastran.__file__)
F06_FILENAME = os.path.join(PKG_PATH, '..', 'models', 'aero', '2_mode_flutter', '0012_flutter.f06')

QT_RED = QColor('red')
QT_WHITE = QColor('white')


PLOT_TYPES = ['eas', 'tas', 'rho', 'q', 'alt', 'kfreq']
UNITS_IN = ['english_in', 'english_kt', 'si', 'si_mm']
UNITS_OUT = UNITS_IN

#FONT_SIZE = 12

class NamedDockWidget(QDockWidget):
    def __init__(self, name: str, widget: QWidget, parent=None):
        QDockWidget.__init__(self, name, parent=parent)
        self.setObjectName(name)
        self.widgeti = widget
        self.setWidget(widget)
        #self.setFont(QFont('Arial', FONT_SIZE))


class Action:
    def __init__(self, name: str, text: str, icon: str='', func=Callable):
        self.name = name
        self.text = text
        self.ico = icon
        self.func = func
    @property
    def icon_path(self) -> str:
        if self.ico == '':
            return self.ico
        return str(ICON_PATH / self.ico)

class LoggableGui(QMainWindow):
    def __init__(self):
        self.html_logging = False
        super().__init__()

    def _logg_msg(self, log_type: str, filename: str, lineno: int, msg: str) -> None:
        """
        Add message to log widget trying to choose right color for it.

        Parameters
        ----------
        log_type : str
            {DEBUG, INFO, ERROR, COMMAND, WARNING} or prepend 'GUI '
        filename : str
            the active file
        lineno : int
            line number
        msg : str
            message to be displayed
        """
        if not self.html_logging:
            # standard logger
            name = '%-8s' % (log_type + ':')
            filename_n = '%s:%s' % (filename, lineno)
            msg2 = ' %-28s %s\n' % (filename_n, msg)
            print(name, msg2)
            return

        # if 'DEBUG' in log_type and not self.settings.show_debug:
        #     return
        # elif 'INFO' in log_type and not self.settings.show_info:
        #     return
        # elif 'COMMAND' in log_type and not self.settings.show_command:
        #     return
        # elif 'WARNING' in log_type and not self.settings.show_warning:
        #     return
        # elif 'ERROR' in log_type and not self.settings.show_error:
        #     return

        if log_type in ['GUI ERROR', 'GUI COMMAND', 'GUI DEBUG', 'GUI INFO', 'GUI WARNING']:
            log_type = log_type[4:]  # drop the GUI

        html_msg = str_to_html(log_type, filename, lineno, msg)
        self._log_msg(html_msg)

    def _log_msg(self, msg: str) -> None:
        """prints an HTML log message"""
        self.log_mutex.lockForWrite()
        text_cursor = self.log_widget.textCursor()
        end = text_cursor.End
        text_cursor.movePosition(end)
        text_cursor.insertHtml(msg)
        self.log_widget.ensureCursorVisible() # new message will be visible
        self.log_mutex.unlock()

    def setup_logging(self):
        self.log = SimpleLogger(
            level='debug', encoding='utf-8',
            log_func=lambda w, x, y, z: self._logg_msg(w, x, y, z))
        # logging needs synchronizing, so the messages from different
        # threads would not be interleave
        self.log_mutex = QtCore.QReadWriteLock()

        self.log_dock_widget = ApplicationLogWidget(self)
        self.log_widget = self.log_dock_widget.log_widget
        #self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock_widget)
        return self.log_widget

    def on_file_exit(self):
        print('on_file_exit')
        if hasattr(self, 'on_file_save') and hasattr(self, 'save_filename'):
            self.on_file_save(self.save_filename)

    def on_file_save(self) -> None:
        if self.save_filename == '' or not os.path.exists(self.save_filename):
            self.on_file_save_as()
        else:
            self._save(self.save_filename)


class MainWindow(LoggableGui):
    def __init__(self):
        super().__init__()
        self.selected_modes = []
        self.freq_tol = -1.0
        self.save_filename = ''
        self.f06_filename = ''
        self.is_valid = False
        self.data = {}
        self.plot_type = 'eas'
        self.freq_lim = [None, None]
        self.damping_lim = [None, None]
        self.kfreq_lim = [None, None]
        self.modes = []

        self.setup_widgets()
        self.setup_layout()
        self.setup_toolbar()
        self.setup_connections()
        self._set_window_title()

    def setup_toolbar(self):
        #frame = QFrame(self)

        actions_input = {
            'file_load': Action(name='file_load', text='Load...', func=self.on_file_load, icon='folder.png'),
            'file_save': Action(name='file_save', text='Save...', func=self.on_file_save, icon='save.png'),
            'file_save_as': Action(name='file_save_as', text='Save As...', func=self.on_file_save_as),
            'file_exit': Action(name='exit', text='Exit...', icon='exit2.jpg', func=self.on_file_exit),
        }
        actions = build_actions(self, actions_input)
        file_actions = ['file_load', 'file_save', 'file_save_as', '', 'file_exit']

        self.menubar = self.menuBar()
        self.file_menu = self.menubar.addMenu('File')
        self.help_menu = self.menubar.addMenu('Help')

        help_actions = []
        menus_dict = {
            'File': (self.file_menu, file_actions),
            'Help': (self.help_menu, help_actions),
        }
        build_menus(menus_dict, actions)
        #self.file_menu.addAction(actions['file_load'])
        #self.file_menu.addAction(actions['exit'])

        #self.toolbar = self.addToolBar('Show toolbar')
        #self.toolbar.setObjectName('main_toolbar')
        self.statusbar = self.statusBar()

    def setup_modes(self):
        self.modes_widget = QListWidget(self)
        self.modes_widget.setMaximumWidth(100)
        self.modes_widget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        _set_modes_table(self.modes_widget, [0])

    # def dontcrash(func):
    #     @wraps(func)
    #     def wrapper(self, *args, **kwargs):
    #         # do something before `sum`
    #         result = func(self) # , *args, **kwargs
    #         # do something after `sum`
    #         return result
    #     return wrapper

    # @dontcrash
    def on_file_save_as(self) -> None:
        #print('on_file_save')
        title = 'Save Case File'
        qt_wildcard = '*.json'
        basedir = str(DIRNAME)
        json_filename, wildcard = getsavefilename(
            self, caption=title, basedir=basedir,
            filters=qt_wildcard,
            #options=QFileDialog.setLabelText('case.json'),
        )
        self.log.info(f'json_filename={json_filename!r} wildcard={wildcard!r}')
        self._save(json_filename)

    def _save(self, json_filename: str):
        is_valid = self.validate()
        self.log.info(f'self.data = {self.data}')
        if json_filename == '':
            return
        #print(f'json_filename={json_filename!r} wildcard={wildcard!r}')
        #print(f'self.data = {self.data}')
        with open(json_filename, 'w') as json_file:
            json.dump(self.data, json_file, indent=4)
        #print(f'fname="{fname}"')
        self.log.info(f'finished saving {json_filename!r}\n')
        self.save_filename = json_filename
        self._set_window_title()

    def _apply_settings(self, data: dict[str, Any]) -> None:
        selected_cases = data['selected_cases']

        line_edits = [
            ('output_directory', self.output_directory_edit),
        ]
        for key, line_edit in line_edits:
            value = data[key]
            line_edit.setText(value)

        pulldown_edits = [
            ('units_in', self.units_in_pulldown, UNITS_IN),
            ('units_out', self.units_out_pulldown, UNITS_OUT),
        ]
        for key, pulldown_edit, values in pulldown_edits:
            value = data[key]
            index = values.index(value)
            pulldown_edit.setCurrentIndex(index)

    # @dontcrash
    def on_file_load(self) -> None:
        #print('on_file_load')
        title = 'Load Case File'
        #dirname = os.path.dirname(CASE_DIRNAME)

        "All Files (*);;Python Files (*.py)"
        qt_wildcard = 'JSON File (*.json)'
        #dlg = QFileDialog(caption=title, directory=dirname, filter=filter_str)
        default_dirname = CASE_DIRNAME
        #if default_filename is None:
            #default_filename = self.gui.last_dir

        default_dirname_str = str(default_dirname)
        if 0:
            fname, wildcard_level = getopenfilename(
                parent=self, caption=title,
                basedir=default_dirname_str, filters=qt_wildcard,
                selectedfilter='', options=None)
        json_filename, wildcard = QFileDialog.getOpenFileName(
            self, title,
            default_dirname_str,
            filter=qt_wildcard, #options=options,
        )
        if json_filename == '':
            return
        self.save_filename = json_filename
        with open(json_filename, 'r') as json_file:
            data = json.load(json_file)
        is_valid = validate_json(data, self.log)
        self._apply_settings(data)
        self.log.info(f'finished loading {json_filename!r}')
        #return wildcard_level, fname
        self._set_window_title()

    def _set_window_title(self):
        if self.save_filename == '':
            self.setWindowTitle('Flutter Plot')
        else:
            self.setWindowTitle(f'Flutter Plot: {self.save_filename}')

    def setup_widgets(self):
        self.f06_filename_label = QLabel('F06 Filename:')
        self.f06_filename_edit = QLineEdit('', self)
        self.f06_filename_edit.setText(F06_FILENAME)
        self.f06_filename_browse = QPushButton('Browse...', self)
        # self.f06_filename_edit.setDisabled(True)
        # self.f06_filename_browse.setDisabled(True)

        self.show_points_checkbox = QCheckBox('Show Points', self)
        self.show_lines_checkbox = QCheckBox('Show Lines', self)
        self.show_points_checkbox.setChecked(True)
        self.show_lines_checkbox.setChecked(True)

        self.xlim_label = QLabel('X Limits:')
        self.xlim_edit_min = QLineEdit('0')
        self.xlim_edit_max = QLineEdit()

        self.damp_lim_label = QLabel('Damping Limits:')
        self.damp_lim_edit_min = QLineEdit('-0.3')
        self.damp_lim_edit_max = QLineEdit('0.3')

        self.freq_lim_label = QLabel('Freq Limits:')
        self.freq_lim_edit_min = QLineEdit('0')
        self.freq_lim_edit_max = QLineEdit()

        self.kfreq_lim_label = QLabel('KFreq Limits:')
        self.kfreq_lim_edit_min = QLineEdit()
        self.kfreq_lim_edit_max = QLineEdit()

        self.freq_tol_label = QLabel('Freq Tol (Hz):')
        self.freq_tol_edit = QLineEdit()

        self.subcase_label = QLabel('Subcase:')
        self.subcase_edit = QComboBox()

        units_msg = (
            "english_in: inch/s, slich/in^3\n"
            "english_ft: ft/s,   slug/ft^3\n"
            "english_kt: knots,  slug/ft^3\n"
            "si:         m/s,    kg/m^3\n"
            "si-mm:      mm/s,   Mg/mm^3\n"
        )
        self.plot_type_label = QLabel('Plot Type:')
        self.plot_type_pulldown = QComboBox(self)
        self.plot_type_pulldown.addItems(PLOT_TYPES)
        #self.plot_type_pulldown.setToolTip(units_msg)

        self.units_in_label = QLabel('Units In:')
        self.units_in_pulldown = QComboBox(self)
        self.units_in_pulldown.addItems(UNITS_IN)
        self.units_in_pulldown.setToolTip(units_msg)
        iunits_in = UNITS_IN.index('english_in')
        self.units_in_pulldown.setCurrentIndex(iunits_in)

        self.units_out_label = QLabel('Units Out:')
        self.units_out_pulldown = QComboBox()
        self.units_out_pulldown.addItems(UNITS_OUT)
        self.units_out_pulldown.setToolTip(units_msg)
        iunits_out = UNITS_IN.index('english_kt')
        self.units_out_pulldown.setCurrentIndex(iunits_out)

        self.output_directory_label = QLabel('Output Directory:')
        self.output_directory_edit = QLineEdit('', self)
        self.output_directory_browse = QPushButton('Browse...', self)
        self.output_directory_edit.setDisabled(True)
        self.output_directory_browse.setDisabled(True)

        self.f06_load_button = QPushButton('Load F06', self)
        self.ok_button = QPushButton('Run', self)
        self.setup_modes()
        self.on_plot_type()

    def on_plot_type(self) -> None:
        plot_type = self.plot_type_pulldown.currentText()

        show_kfreq = False
        if plot_type == 'kfreq':
            show_kfreq = True
        show_damp = not show_kfreq
        self.damp_lim_label.setVisible(show_damp)
        self.damp_lim_edit_min.setVisible(show_damp)
        self.damp_lim_edit_max.setVisible(show_damp)

        self.kfreq_lim_label.setVisible(show_kfreq)
        self.kfreq_lim_edit_min.setVisible(show_kfreq)
        self.kfreq_lim_edit_max.setVisible(show_kfreq)

    def setup_layout(self):
        grid = QGridLayout()
        irow = 0

        grid.addWidget(self.f06_filename_label, irow, 0)
        grid.addWidget(self.f06_filename_edit, irow, 1)
        grid.addWidget(self.f06_filename_browse, irow, 2)
        irow += 1

        grid.addWidget(self.subcase_label, irow, 0)
        grid.addWidget(self.subcase_edit, irow, 1)
        irow += 1

        grid.addWidget(self.plot_type_label, irow, 0)
        grid.addWidget(self.plot_type_pulldown, irow, 1)
        irow += 1

        grid.addWidget(self.xlim_label, irow, 0)
        grid.addWidget(self.xlim_edit_min, irow, 1)
        grid.addWidget(self.xlim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.damp_lim_label, irow, 0)
        grid.addWidget(self.damp_lim_edit_min, irow, 1)
        grid.addWidget(self.damp_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.freq_lim_label, irow, 0)
        grid.addWidget(self.freq_lim_edit_min, irow, 1)
        grid.addWidget(self.freq_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.kfreq_lim_label, irow, 0)
        grid.addWidget(self.kfreq_lim_edit_min, irow, 1)
        grid.addWidget(self.kfreq_lim_edit_max, irow, 2)
        irow += 1

        grid.addWidget(self.freq_tol_label, irow, 0)
        grid.addWidget(self.freq_tol_edit, irow, 1)
        irow += 1

        grid.addWidget(self.units_in_label, irow, 0)
        grid.addWidget(self.units_in_pulldown, irow, 1)
        irow += 1

        grid.addWidget(self.units_out_label, irow, 0)
        grid.addWidget(self.units_out_pulldown, irow, 1)
        irow += 1

        grid.addWidget(self.output_directory_label, irow, 0)
        grid.addWidget(self.output_directory_edit, irow, 1)
        grid.addWidget(self.output_directory_browse, irow, 2)
        self.output_directory_label.setDisabled(True)
        self.output_directory_edit.setDisabled(True)
        self.output_directory_browse.setDisabled(True)
        irow += 1

        radio_hbox = QHBoxLayout()
        radio_hbox.addWidget(self.show_points_checkbox)
        radio_hbox.addWidget(self.show_lines_checkbox)
        radio_hbox.addStretch(1)

        ok_cancel_hbox = QHBoxLayout()
        ok_cancel_hbox.addWidget(self.f06_load_button)
        ok_cancel_hbox.addWidget(self.ok_button)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        vbox.addLayout(radio_hbox)
        vbox.addStretch(1)
        vbox.addLayout(ok_cancel_hbox)

        #log_widget = ApplicationLogWidget(self)

        log_widget = self.setup_logging()
        self.modes_dock_widget = NamedDockWidget('Modes', self.modes_widget, self)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.modes_dock_widget)
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock_widget)

        widget = QWidget()
        widget.setLayout(vbox)
        self.setCentralWidget(widget)

    def setup_connections(self):
        self.f06_load_button.clicked.connect(self.on_load_f06)
        self.plot_type_pulldown.currentIndexChanged.connect(self.on_plot_type)
        #self.ok_button.clicked.connect(self.on_ok)
        # self.modes_widget.itemSelectionChanged.connect(self.on_modes)
        self.ok_button.clicked.connect(self.on_ok)

    # @dontcrash
    def on_load_f06(self):
        f06_filename = os.path.abspath(self.f06_filename_edit.text())
        if not os.path.exists(f06_filename):
            self.log.error(f"can't find {f06_filename}")
            return
        self.f06_filename = f06_filename
        f06_units = self.units_in_pulldown.currentText()
        out_units = self.units_out_pulldown.currentText()
        print(f06_units, out_units)
        try:
            self.response: FlutterResponse = make_flutter_response(
                f06_filename,
                f06_units=f06_units, out_units=out_units,
                make_alt=True,
                log=self.log)[1]
        except Exception as e:
            self.log.error(str(e))
            return
        self.update_modes_table(self.response.modes)

    def update_modes_table(self, modes: list[int]) -> None:
        self.modes = modes
        _set_modes_table(self.modes_widget, modes)
        self.log.info(f'modes = {self.modes}')

    def on_modes(self) -> None:
        self.validate()
        if not self.is_valid:
            return
        self.plot(self.modes)

    # @dontcrash
    def on_ok(self):
        is_valid = self.validate()
        if not is_valid:
            return

        modes = self.selected_modes
        if len(modes) == 0:
            self.log.warning(f'modes = {modes}; assuming all modes -> {self.modes}')
            modes = self.modes
            #return
        self.log.info(f'is_valid = {is_valid}\n')
        self.is_valid = True
        self.plot(modes)

    # @dontcrash
    def plot(self, modes: list[int]) -> None:
        self.log.info(f'plot; modes = {modes}\n')
        if not self.is_valid:
            self.log.warning(f'no valid\n')
            return
        plot_type = self.plot_type
        self.log.info(f'plot_type = {plot_type}\n')
        noline = not self.show_lines_checkbox.isChecked()
        nopoints = not self.show_points_checkbox.isChecked()
        # self.log.info(f'noline = {noline}\n')
        # self.log.info(f'nopoints = {nopoints}\n')

        freq_tol = self.freq_tol
        self.log.info(f'freq_tol = {freq_tol}\n')
        if noline and nopoints:
            noline = True
            nopoints = True

        xlim = self.xlim
        ylim_damping = self.damping_lim
        ylim_kfreq = None
        ylim_freq = None
        damping_limit = None

        dirname = os.path.dirname(self.f06_filename)
        basename = os.path.basename(self.f06_filename)

        base = os.path.splitext(basename)[0]
        png_filename = base + '.png'
        current_directory = os.getcwd()
        sys.stdout.flush()
        os.chdir(dirname)

        fig = plt.figure(1)
        fig.clear()
        gridspeci = gridspec.GridSpec(2, 4)
        damp_axes = fig.add_subplot(gridspeci[0, :3])
        freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        if plot_type == 'kfreq':
            # xlabel = r'KFreq [rad]; $ \omega c / (2 V)$'
            # ylabel1 = r'Structural Damping; $g = 2 \gamma $'
            # ylabel2 = 'Frequency [Hz]'
            self.response.plot_kfreq_damping2(
                fig=fig, axes1=damp_axes, axes2=freq_axes,
                modes=modes,
                #xlim=None,
                ylim1=ylim_damping, ylim2=ylim_freq,
                show=True, clear=True,
                close=False, legend=True,
                png_filename=png_filename,
            )
        else:
            try:
                self.response.plot_vg_vf(
                    fig=fig, damp_axes=damp_axes, freq_axes=freq_axes,
                    plot_type=plot_type,
                    legend=True,
                    xlim=xlim, ylim_damping=ylim_damping, ylim_freq=ylim_freq,
                    freq_tol=freq_tol,
                    modes=modes, show=True,
                    noline=noline, nopoints=nopoints,
                    vd_limit=None, damping_limit=damping_limit,
                    png_filename=png_filename,
                )
            except Exception as e:
                self.log.error(str(e))
        os.chdir(current_directory)
        #self.get_xlim()
        str(self.xlim)

    def get_xlim(self):
        default_xlim = [None, None]
        xlim_min, is_passed1 = get_float_or_none(self.xlim_edit_min)
        xlim_max, is_passed2 = get_float_or_none(self.xlim_edit_max)
        damp_lim_min, is_passed3 = get_float_or_none(self.damp_lim_edit_min)
        damp_lim_max, is_passed4 = get_float_or_none(self.damp_lim_edit_max)
        freq_tol, is_passed5 = get_float_or_none(self.freq_tol_edit)
        if is_passed5 and freq_tol is None:
            freq_tol = -1.0
        xlim = [xlim_min, xlim_max]
        damp_lim = [damp_lim_min, damp_lim_max]
        self.log.info(f'XLim = {xlim}')
        is_passed_flags = [is_passed1, is_passed2, is_passed3, is_passed4, is_passed5]
        is_passed = all(is_passed_flags)
        #print(f'is_passed_flags = {is_passed_flags}')
        print(f'freq_tol = {freq_tol}')
        return xlim, damp_lim, freq_tol, is_passed

    def get_selected_modes(self) -> list[str]:
        mode_strs = get_selected_items_flat(self.modes_widget)
        modes = [int(mode_str.split(' ')[1]) for mode_str in mode_strs]
        self.log.info(f'modes = {modes}')
        return modes

    def validate(self) -> bool:
        xlim, ydamp_lim, freq_tol, is_valid_xlim = self.get_xlim()
        subcase = 1
        selected_modes = self.get_selected_modes()
        self.selected_modes = selected_modes
        self.xlim = xlim
        self.ylim = ydamp_lim
        self.freq_tol = freq_tol

        self.plot_type = self.plot_type_pulldown.currentText()
        units_in = self.units_in_pulldown.currentText()
        units_out = self.units_out_pulldown.currentText()
        output_directory = self.output_directory_edit.text()
        self.data = {
            'modes': selected_modes,
            'xlim': xlim,
            'damp_lim': ydamp_lim,
            'subcase': subcase,
            'output_directory': output_directory,
            'units_in': units_in,
            'units_out': units_out,
        }
        is_passed = all([is_valid_xlim])
        if is_passed:
            self.xlim = xlim
            self.ylim = ydamp_lim
            #self.data = data
            # is_valid = validate_json(self.data, self.log)
            #if is_valid != is_passed:
            self.log.info(f'passed data:\n{str(self.data)}')
        else:
            self.log.error(f'failed data:\n{str(self.data)}')
        return is_passed


def get_selected_items_flat(list_widget: QListWidget) -> list[str]:
    items = list_widget.selectedItems()
    texts = []
    for item in items:
        text = item.text()
        texts.append(text)
    return texts

def build_menus(menus_dict: dict[str, tuple[QMenu, list[str]]],
                actions: dict[str, QAction]) -> None:
    for menu_name, (menu, actions_list) in menus_dict.items():
        assert isinstance(menu_name, str), menu_name
        for action_name in actions_list:
            assert isinstance(action_name, str), action_name
            if action_name == '':
                menu.addSeparator()
                continue
            #print(menu_name, action_name)
            action = actions[action_name]
            menu.addAction(action)
            #print('menu = ', menu)
            #print('action = ', action)

def build_actions(self, actions_input: dict[str, Action]) -> dict[str, QAction]:
    actions = {}
    tip = ''
    for name, action_inputi in actions_input.items():
        func = action_inputi.func
        icon_path = action_inputi.icon_path
        txt = action_inputi.text
        if icon_path and 0:
            ico = QIcon(icon_path)
            #print('icon_path = ', icon_path)
            assert os.path.exists(icon_path), icon_path
            action = QAction(txt, self)
            action.setIcon(ico)
        else:
            action = QAction(txt, self)
        #action.setText(name)
        #action.setIcon()
        #action.setCheckable()
        #action.setShortcut(QKeySequence(name))
        if tip:
            print(f'  found tip for {name}')
            action.setStatusTip(tip)
        if func:
            print(f'  found func for {name}')
            action.triggered.connect(func)
        actions[name] = action
    return actions


def validate_json(data: dict[str, Any],
                  log: SimpleLogger) -> bool:
    is_valid = True
    log.warning(f'keys = {list(data.keys())}')
    key_allowed_values = [
        ('units_in', UNITS_IN),
        ('units_out', UNITS_OUT),
        ('plot_type', PLOT_TYPES),
    ]
    for (key, allowed_values) in key_allowed_values:
        if key not in data:
            is_valid = False
            log.error(f'data[{key}] is missing; defaulting to {allowed_values[0]}')
            default_value = allowed_values[0]
            data[key] = default_value
            continue

        value = data[key]
        if value not in allowed_values:
            is_valid = False
            log.error(f'{key}={value!r} not in {allowed_values}; defaulting to {allowed_values[0]}')
            default_value = allowed_values[0]
            data[key] = default_value
    return is_valid

def get_float_or_none(line_edit: QLineEdit) -> tuple[Optional[float], bool]:
    is_passed = False
    text = line_edit.text().strip()
    if len(text) == 0:
        value = None
        is_passed = True
    else:
        try:
            value = float(text)
            is_passed = True
        except ValueError:
            value = None
            is_passed = False
    return value, is_passed


def _set_modes_table(modes_widget: QListWidget, modes: list[int]):
    modes_widget.clear()
    for imode in modes:
        mode = QListWidgetItem(f"Mode {imode}")
        mode.setSelected(True)
        modes_widget.addItem(mode)


def main():  # pragma: no cover
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    # The Main window

    main_window = MainWindow()
    main_window.show()
    # Enter the main loop
    app.exec_()


if __name__ == '__main__':
    main()
