"""
TODO: add cut, copy, paste, delete support to table
TODO: add trade study support vs. just report writing
TODO: better saving of new parameters
"""
import os
import copy
import warnings
import traceback
from pathlib import Path
from functools import wraps
from typing import Any

try:
    import json5 as json
except ModuleNotFoundError:
    warnings.warn('couldnt find json5, using json')
    import json

from pyNastran.f06.dev.flutter.utils import get_raw_json
JSON_FILENAME, USE_VTK, USE_TABS = get_raw_json(allow_vtk=False)
from pyNastran.f06.dev.flutter.trade_layout import TradeLayout
from pyNastran.f06.dev.flutter.plot_layout import PlotLayout


# from qtpy import QtCore
# from qtpy.compat import getopenfilename  # getsavefilename
# from qtpy.QtGui import QIcon, QPixmap
from qtpy.QtWidgets import (
    QWidget, QApplication,
    QAction, QTabWidget,)

import pyNastran
from pyNastran.gui.utils.qt.pydialog import make_font
from pyNastran.gui.qt_files.loggable_gui import LoggableGui

from pyNastran.f06.dev.flutter.actions_builder import Actions, Action, build_menus
from pyNastran.f06.dev.flutter.preferences_object import FlutterPreferencesObject
from pyNastran.f06.dev.flutter.preferences import (
    FLUTTER_BBOX_TO_ANCHOR_DEFAULT, LEGEND_LOC_DEFAULT,
    FONT_SIZE_DEFAULT, FLUTTER_NCOLUMNS_DEFAULT, FREQ_NDIGITS_DEFAULT, FREQ_DIVERGENCE_TOL)

from pyNastran.f06.dev.flutter.utils_qt import (
    load_checkboxs, load_lineedits, load_pulldowns, load_min_max_lineedits)
from pyNastran.f06.dev.flutter.utils import (
    validate_json, get_point_removal_str, MODE_SWITCH_METHODS,
)

PKG_PATH = Path(pyNastran.__path__[0])

AERO_PATH = PKG_PATH / '..' / 'models' / 'aero'

from pyNastran.f06.dev.flutter.vtk_data import VtkData

QLINEEDIT_WHITE = 'QLineEdit {background-color: white;}'
QLINEEDIT_RED = 'QLineEdit {background-color: red;}'
ICON_PATH = Path('')


class FlutterGui(LoggableGui):
    def __init__(self, f06_filename: str=''):
        super().__init__(html_logging=False)
        self.use_vtk = USE_VTK
        self.use_tabs = USE_TABS
        self.vtk_data = VtkData()
        self.base_settings_dict = {}

        self._export_settings_obj = FlutterPreferencesObject(self, USE_VTK, hide_vtk=True)
        self.iwindows = []

        self.divergence_legend_loc = LEGEND_LOC_DEFAULT
        self.flutter_bbox_to_anchor_x = FLUTTER_BBOX_TO_ANCHOR_DEFAULT
        self.flutter_ncolumns = FLUTTER_NCOLUMNS_DEFAULT
        self.freq_ndigits = FREQ_NDIGITS_DEFAULT
        self.freq_divergence_tol = FREQ_DIVERGENCE_TOL
        self.auto_update = True

        # Mach	Fuel	Filename
        # 0.2	0	mach0.2_fuel_0.f06
        # 0.3	40	mach0.3_fuel_40.f06
        # 0.4	100	mach0.4_fuel_100.f06
        # self.base_f06_directory = r'C:\work\code\pyNastran\models\aero\2_mode_flutter\dev'
        # self.base_f06_directory = r'C:\NASA\m4\formats\git\pyNastran\models\aero\2_mode_flutter'

        self.font_size = FONT_SIZE_DEFAULT
        self.plot_font_size = FONT_SIZE_DEFAULT
        self.show_lines = True
        self.show_points = True
        self.show_mode_number = False
        self.show_detailed_mode_info = False
        self.point_spacing = 0
        self.use_rhoref = False
        self._units_in = ''
        self._units_out = ''
        self.units_in = ''
        self.units_out = ''
        self.use_dock_widgets = self.html_logging
        self.qactions = {}
        self.nrecent_files_max = 20
        self.recent_files = []
        self.save_filename = JSON_FILENAME
        self.is_valid = False

        self.data = {}
        self.f06_filename = ''
        self.bdf_filename = ''
        self.op2_filename = ''
        self._bdf_filename_default = ''
        self._op2_filename_default = ''
        # self.subcase = 0
        self.x_plot_type = 'eas'
        self.plot_type = 'x-damp-freq'
        self.mode_switch_method = MODE_SWITCH_METHODS[0]
        self.point_removal = ''

        self.eas_lim = []
        self.tas_lim = []
        self.mach_lim = []
        self.alt_lim = []
        self.q_lim = []
        self.rho_lim = []
        self.eas_flutter_range = [None, None]
        # self.eas_diverg_range = [None, None]
        self.freq_lim = [None, None]
        # self.damping_lim = [None, None]
        self.ydamp_lim = [None, None]
        self.kfreq_lim = [None, None]
        self.eigr_lim = [None, None]
        self.eigi_lim = [None, None]
        # self.responses = {}
        # self.modes = []
        # self.selected_modes = []
        self.freq_tol = -1.0
        self.freq_tol_remove = -1.0
        self.mag_tol = -1.0
        self.damping = -1.0
        self.damping_required = -1.0
        self.damping_required_tol = -1.0
        self.vf = -1.0
        self.vl = -1.0
        self.export_to_png = True
        self.export_to_csv = False
        self.export_to_f06 = False
        self.export_to_zaero = False
        #----------------------------------------
        log_widget = self.setup_logging()
        self.plot_layout = [PlotLayout(self, ifile=0)]
        self.trade_layout = TradeLayout(self)
        self.tabs = QTabWidget()

        self.setup_layout()
        self.on_load_settings()
        if f06_filename:
            for plot_layout in self.plot_layout:
                plot_layout.f06_filename_edit.setText(f06_filename)
            # self._set_f06_default_names(f06_filename)

        self.setup_toolbar()
        self._update_recent_files_actions()
        self._set_window_title()
        self.on_font_size()
        for plot_layout in self.plot_layout:
            plot_layout.on_plot_type()
        self.setAcceptDrops(True)
        # self.on_open_new_window()
        self.show()

    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event) -> None:
        plot_layout = self.plot_layout[self.ifile]
        filenames = [url.toLocalFile() for url in event.mimeData().urls()]
        for filename in filenames:
            flower = filename.lower()
            if flower.endswith(('.f06', '.out')):
                plot_layout.f06_filename_edit.setText(filename)
            else:
                self.log.error(f'unknown extension (f06) format for {filename}')

    def setup_toolbar(self) -> None:
        # frame = QFrame(self)
        actions_dict = {
            # 'file_load': Action(name='file_load', text='Load...', func=self.on_file_load, icon='folder.png'),
            # 'file_save': Action(name='file_save', text='Save...', func=self.on_file_save, icon='save.png'),
            # 'file_save_as': Action(name='file_save_as', text='Save As...', func=self.on_file_save_as),
            'file_exit':       Action(name='exit', text='Exit...', icon='exit2.jpg', func=self.on_file_exit),
            'export_settings': Action(name='Export Settings', text='Export Settings...', icon='preferences.jpg',
                                      shortcut='Ctrl+P', func=self.on_export_settings),
        }
        actions_input = Actions(ICON_PATH, actions_dict)  # load_icon=False
        recent_files = actions_input.build_recent_file_qactions(
            self, self.recent_files, self.set_f06)
        self.qactions = actions_input.build_qactions(self)

        file_actions = [
            # 'file_load',
            # 'file_save', 'file_save_as',
            ] + recent_files + [
            'file_exit']
        view_actions = ['export_settings']

        self.menubar = self.menuBar()
        self.file_menu = self.menubar.addMenu('File')
        self.view_menu = self.menubar.addMenu('View')
        # self.help_menu = self.menubar.addMenu('Help')

        # help_actions = []
        menus_dict = {
            'File': (self.file_menu, file_actions),
            'View': (self.view_menu, view_actions),
            # 'Help': (self.help_menu, help_actions),
        }
        build_menus(menus_dict, self.qactions)
        # self.file_menu.addAction(actions['file_load'])
        # self.file_menu.addAction(actions['exit'])

        # self.toolbar = self.addToolBar('Show toolbar')
        # self.toolbar.setObjectName('main_toolbar')
        self.statusbar = self.statusBar()

    def dont_crash(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            # do something before `sum`
            try:
                result = func(self, *args, **kwargs)
            except Exception as error:
                self.log_error(str(traceback.format_exc()))
                result = None
            # do something after `sum`
            return result
        return wrapper

    def set_f06(self, ifile: int) -> None:
        f06_filename = self.recent_files[ifile]
        plot_layout = self.plot_layout[self.ifile]
        plot_layout.f06_filename_edit.setText(f06_filename)
        plot_layout.on_load_f06(None)

    def dontcrash(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            # do something before `sum`
            result = func(self) # , *args, **kwargs
            # do something after `sum`
            return result
        return wrapper

    # @dontcrash
    def on_export_settings(self) -> None:
        self._export_settings_obj.show()

    def on_file_exit(self) -> None:
        if hasattr(self, 'on_file_save') and hasattr(self, 'save_filename'):
            self.on_file_save()

    def on_file_save(self) -> None:
        if self.save_filename == '' or not os.path.exists(self.save_filename):
            self.on_file_save_as()
        else:
            # self.log.warning('on_file_save; _save')
            self._save(self.save_filename)

    def on_file_save_as(self) -> None:
        # if 0:
        #     #print('on_file_save')
        #     qt_wildcard = '*.json'
        #     basedir = str(DIRNAME)
        #     json_filename, wildcard = getsavefilename(
        #         self, caption=title, basedir=basedir,
        #         filters=qt_wildcard,
        #         #options=QFileDialog.setLabelText('data.json'),
        #     )
        #     self.log.info(f'json_filename={json_filename!r} wildcard={wildcard!r}')
        json_filename = self.save_filename
        # self.log.warning('on_file_save_as; _save')
        self._save(json_filename)

    def _save(self, json_filename: str) -> None:
        is_valid = self.validate()
        # self.log.info(f'self.data = {self.data}; is_valid={is_valid}')
        if json_filename == '' or len(self.data) == 0:
            return
        # print(f'json_filename={json_filename!r} wildcard={wildcard!r}')
        # print(f'self.data = {self.data}')

        # we are copying the OG dict so we don't modify that one
        out_data = copy.deepcopy(self.data)
        out_data['use_vtk'] = self.use_vtk
        if USE_VTK:
            out_data['vtk'] = self._vtk_window_obj.data
        out_data['preferences'] = {
            'flutter_ncolumns': self.flutter_ncolumns,
            'freq_ndigits': self.freq_ndigits,
            'freq_divergence_tol': self.freq_divergence_tol,
            'auto_update': self.auto_update,
            'flutter_bbox_to_anchor_x': self.flutter_bbox_to_anchor_x,
            'divergence_legend_loc': self.divergence_legend_loc,
            'export_to_png': self.export_to_png,
            'export_to_csv': self.export_to_csv,
            'export_to_f06': self.export_to_f06,
            'export_to_zaero': self.export_to_zaero,
            'font_size': self.font_size,
            'plot_font_size': self.plot_font_size,
        }
        trade_dict = self.trade_layout.get_dict()
        out_data.update(trade_dict)
        # out_data['trade'] = trade_dict  # TODO and export trade

        with open(json_filename, 'w') as json_file:
            json.dump(out_data, json_file, indent=4)
        # print(f'fname="{fname}"')
        self.log.info(f'finished saving {json_filename!r}\n')
        # self.log.info(f'out_data = {out_data}; is_valid={is_valid}')
        self.save_filename = json_filename
        self._set_window_title()

    def _apply_settings(self, data: dict[str, Any]) -> None:
        trade_layout = self.trade_layout
        if USE_VTK:
            self.vtk_data.apply_settings(data)
            self._vtk_window_obj.apply_settings(data)
        log = self.log
        font_size0 = self.font_size

        self.base_settings_dict = copy.deepcopy(data)
        # radios = [
        #     ('show_points', self.show_points_radio),
        # ]
        # for (key, checkbox) in radios:
        #     if key not in data:
        #         continue
        #     val = data[key]
        #     assert isinstance(val, bool), (key, val)
        #     checkbox.setChecked(val)

        # preferences = data.get('preferences', {})
        # for key, value in preferences.items():
        #     if not hasattr(self, key):
        #         self.log.error(f'load failure: skipping {key!r}={value!r} because key doesnt exist')
        #         continue
        #     self.log.error(f'load success: loaded {key!r}={value!r}')
        #     setattr(self, key, value)

        spinners = [
             # ('plot_font_size', self.plot_font_size_edit),
        ]
        for (key, spinner) in spinners:
            if key not in data:
                continue
            val = data[key]
            assert isinstance(val, int), (key, val)
            spinner.setValue(val)

        type_names = [
            (int,  ('preferences/font_size',
                    'preferences/plot_font_size',
                    'preferences/flutter_ncolumns')),
            (float, ('preferences/flutter_bbox_to_anchor_x',)),
            (str, ('preferences/divergence_legend_loc',)),
            (bool, ('preferences/export_to_png',
                    'preferences/export_to_f06',
                    'preferences/export_to_csv',
                    'preferences/export_to_zaero')),
        ]
        for value_type, keys in type_names:
            for key in keys:
                # print(f'loading key={key!r}')
                if '/' in key:
                    skey = key.split('/')
                    assert len(skey) == 2, f'key={key!r} skey={skey}'
                    key0, key1 = skey
                    if key0 not in data:
                        log.warning(f'skipping {key!r} because {key0} does not exist')
                        print(data)
                        continue
                    data0 = data[key0]
                    value = data0[key1]
                    assert hasattr(self, key1), (key, value)
                    if not isinstance(value, value_type):
                        log.warning(f'{key!r}={value!r} and is not {value_type}...skipping')
                        continue
                    # log.info(f'setting {key!r} (key1={key1!r}) -> {value!r}')
                    setattr(self, key1, value)
                else:
                    if key not in data:
                        log.warning(f'skipping {key!r}')
                        assert len(key) > 1, keys
                        continue
                    value = data[key]
                    assert hasattr(self, key), (key, value)
                    if not isinstance(value, value_type):
                        log.warning(f'{key!r}={value!r} and is not {value_type}...skipping')
                        continue
                    # log.info(f'setting {key!r} -> {value!r}')
                    setattr(self, key, value)

        checkboxs = []
        for plot_layout in self.plot_layout:
            checkboxs.extend(plot_layout.get_checkboxs())

        load_checkboxs(data, checkboxs)

        min_max_lineedits = []
        for plot_layout in self.plot_layout:
            min_max_lineedits.extend(plot_layout.get_min_max_lineedits())
        load_min_max_lineedits(data, min_max_lineedits)

        point_removal = data.get('point_removal', [])
        point_removal_str = get_point_removal_str(point_removal)
        for plot_layout in self.plot_layout:
            plot_layout.point_removal_edit.setText(point_removal_str)

        pulldown_edits = []
        line_edits = []
        for plot_layout in self.plot_layout:
            line_edits.extend(plot_layout.get_lineedits())
            pulldown_edits.extend(plot_layout.get_comboboxs())

        line_edits.extend(trade_layout.get_lineedits())
        pulldown_edits.extend(trade_layout.get_comboboxs())
        load_lineedits(data, line_edits)
        load_pulldowns(data, pulldown_edits)

        self.recent_files = []
        for fname in data['recent_files']:
            abs_path = os.path.abspath(fname)
            if abs_path not in self.recent_files:
                self.recent_files.append(abs_path)
        self.f06_filename = self.recent_files[0]
        if self.font_size != font_size0:
            self.on_set_font_size(self.font_size)

    def on_load_settings(self) -> None:
        json_filename = self.save_filename
        if not os.path.exists(json_filename):
            self.log.warning(f'unable to find {json_filename}')
            return
        # self.save_filename = json_filename
        try:
            with open(json_filename, 'r') as json_file:
                data = json.load(json_file)
            is_valid = validate_json(data, self.log)
            self._apply_settings(data)
        except Exception as e:
            self.log.error(f'failed to load {json_filename}\n{str(e)}')
            print(traceback.format_exc())
            # print(traceback.format_exception_only(e))
            # raise
            return
        self.log.info(f'finished loading {json_filename!r}')
        # return wildcard_level, fname
        self._set_window_title()

    def _set_window_title(self) -> None:
        if self.save_filename == '':
            self.setWindowTitle('Flutter Plot')
        else:
            self.setWindowTitle(f'Flutter Plot: {self.save_filename}')

    def setup_layout(self) -> None:
        tab_file = QWidget()
        tab_org = QWidget()
        iwindow = self.tabs.addTab(tab_file, 'File')
        iwindow_organize = self.tabs.addTab(tab_org, 'Organize')
        self.iwindows = [iwindow, iwindow_organize]

        tab_file.setLayout(self.plot_layout[0])
        tab_org.setLayout(self.trade_layout)
        self.setCentralWidget(self.tabs)
        tab_file.activateWindow()

    @property
    def ifile(self) -> int:
        """gets the active file index"""
        if not USE_TABS:
            return 0
        index = self.tabs.currentIndex()
        out = index - 1
        return out

    def on_new_tab(self):
        """"
        [compare, file1, +]
        """
        index = self.tabs.currentIndex()
        ifile = len(self.iwindows) - 1
        if index != ifile:
            return
        tab_plus = QWidget(self)
        iwindow_file2 = self.iwindows[-1]
        iwindow = self.tabs.addTab(tab_plus, '+')
        self.iwindows.append(iwindow)
        self.tabs.setTabText(iwindow_file2, f'File {ifile}')

    def on_font_size(self) -> None:
        # font_size = self.font_size_edit.value()
        self.on_set_font_size(self.font_size)

    def on_set_font_size(self, font_size: int) -> None:
        self.font_size = font_size
        font = make_font(font_size, is_bold=False)
        self.setFont(font)

    def validate(self) -> bool:
        return self.plot_layout[self.ifile].validate()

    def add_recent_file(self, f06_filename: str) -> None:
        path = os.path.abspath(f06_filename)
        if path in self.recent_files:
            self.recent_files.remove(path)
        self.recent_files = [path] + self.recent_files
        if len(self.recent_files) > self.nrecent_files_max:
            self.recent_files = self.recent_files[:self.nrecent_files_max]

        self._update_recent_files_actions()
        is_valid = self.validate()
        if is_valid:
            self.log.warning('add_recent_file; save')
            self._save(self.save_filename)

    def _update_recent_files_actions(self) -> None:
        current_directory = os.getcwd()
        for ifile in range(self.nrecent_files_max):
            name = f'file_{ifile}'
            try:
                fname = self.recent_files[ifile]
                show = True
            except IndexError:
                show = False

            action: QAction = self.qactions[name]
            active_text = action.text()

            if show != action.isVisible():
                action.setVisible(show)

            if show and active_text != fname:
                # only update text if it's visible
                relative_fname = os.path.relpath(fname, current_directory)
                text = fname
                if len(relative_fname) < len(fname):
                    text = relative_fname
                action.setText(text)

    def log_debug(self, msg: str) -> None:
        print(f'DEBUG: {msg}')

    def log_info(self, msg: str) -> None:
        print(f'INFO:  {msg}')

    def log_command(self, msg: str) -> None:
        print(f'COMMAND: {msg}')

    def log_warning(self, msg: str) -> None:
        print(f'WARNING: {msg}')

    def log_error(self, msg: str) -> None:
        print(f'ERROR:   {msg}')


def main(f06_filename: str='') -> None:  # pragma: no cover
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)

    # The Main window
    main_window = FlutterGui(f06_filename)
    # Enter the main loop
    app.exec_()


if __name__ == '__main__':  # pragma: no cover
    # out = split_by_pattern(['cat_dog_1.0_asdf', 'cat_dog_2.0_qwer'])
    # print(out)
    main()
