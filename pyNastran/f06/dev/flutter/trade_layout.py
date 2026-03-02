from __future__ import annotations
import os
import warnings
import traceback
from pathlib import Path
from functools import wraps
from typing import Any, TYPE_CHECKING
import natsort

from pyNastran.utils import print_bad_path
from pyNastran.utils.dev import get_files_of_type

try:
    import json5 as json
except ModuleNotFoundError:
    warnings.warn('couldnt find json5, using json')
    import json

from pyNastran.f06.dev.flutter.utils import get_raw_json
JSON_FILENAME, USE_VTK, USE_TABS = get_raw_json(allow_vtk=False)

from qtpy.QtWidgets import (
    QLabel,
    QApplication, QVBoxLayout, QComboBox,
    QPushButton, QGridLayout,
    QLineEdit, QFileDialog, QProgressBar,
)
from pyNastran.f06.dev.flutter.qtablewidgetcopy import QTableWidgetCopy
from pyNastran.f06.dev.flutter.utils_qt import create_grid_from_list
from pyNastran.f06.dev.flutter.write_report import write_report


from cpylog import SimpleLogger
import pyNastran

from pyNastran.f06.flutter_response import Limit  # FlutterResponse

PKG_PATH = Path(pyNastran.__path__[0])

AERO_PATH = PKG_PATH / '..' / 'models' / 'aero'

import pandas as pd
import tables
if TYPE_CHECKING:
    from pyNastran.f06.dev.flutter.gui_flutter_plot import FlutterGui

# QLINEEDIT_WHITE = 'QLineEdit {background-color: white;}'
# QLINEEDIT_RED = 'QLineEdit {background-color: red;}'

AXIS_VALUES = ['Equivalent Airspeed, EAS', 'Altitude',
               'Mach', 'Dynamic Pressure, Q']

class TradeLayout(QVBoxLayout):
    def dontcrash(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            # do something before `sum`
            result = func(self) # , *args, **kwargs
            # do something after `sum`
            return result
        return wrapper

    def __init__(self, parent: FlutterGui, *args, **kwargs):
        """
        Parameters
        ----------
        parent : FlutterGui
            this is a QApplication
            source for self.excel_filename, ...
        """
        super().__init__(*args, **kwargs)
        self.parent = parent
        self.log = parent.log
        #------------------------------------------------------------
        self.excel_dict: dict[str, pd.DataFrame] = {}
        self.excel_filename = 'trade_study.xlsx'
        self.base_f06_directory = ''
        self.word_filename = 'flutter_summary.docx'
        #------------------------------------------------------------

        self.excel_filename_label = QLabel('Excel Filename:', parent)
        self.excel_filename_edit = QLineEdit(parent)
        self.excel_filename_edit.setText(self.excel_filename)
        self.excel_filename_edit.setToolTip('Must click load button before selecting tab')
        self.excel_filename_browse = QPushButton('Browse...', parent)

        self.base_f06_directory_label = QLabel('Base F06 Directory:', parent)
        self.base_f06_directory_edit = QLineEdit(parent)
        # self.base_f06_directory_edit.setText(self.base_f06_directory)
        self.base_f06_directory_browse = QPushButton('Browse...', parent)
        self.base_f06_directory_load = QPushButton('Load...', parent)

        self.word_filename_label = QLabel('Word Filename:', parent)
        self.word_filename_edit = QLineEdit(parent)
        self.word_filename_edit.setText(self.word_filename)
        self.word_filename_browse = QPushButton('Browse...', parent)

        self.load_excel_button = QPushButton('Load...', parent)
        self.tab_select_pulldown = QComboBox(parent)

        self.tab_label = QLabel('Tab Name:', parent)
        self.tab_edit = QLineEdit(parent)
        self.tab_edit.setVisible(False)

        # self.tab_select_browse.setEnabled(False)
        # self.tab_select_browse.setToolTip('Must click load button before selecting tab')

        self.run_organize_button = QPushButton('Run', parent)

        # self.starting_row = QLabel('Starting Row:', parent)
        # self.starting_row = QLineEdit(parent)
        #
        # self.starting_col = QLabel('Starting Col:', parent)
        # self.starting_col = QLineEdit(parent)

        self.table_widget = QTableWidgetCopy(parent)
        self.progress_bar = QProgressBar(parent)
        self.progress_bar.setVisible(False)

        self.xaxis_label = QLabel('X-Axis', parent)
        self.yaxis_label = QLabel('Y-Axis', parent)
        self.xaxis_pulldown = QComboBox(parent)
        self.yaxis_pulldown = QComboBox(parent)
        self.xaxis_pulldown.addItems(AXIS_VALUES)
        self.yaxis_pulldown.addItems(AXIS_VALUES)

        self.config_label = QLabel('Configs', parent)
        self.config_edit = QLineEdit(parent)
        self.config_edit.setToolTip('Comma separated list of configs')
        self.xaxis_pulldown.setEnabled(False)
        self.yaxis_pulldown.setEnabled(False)
        # self.config_edit.setEnabled(False)
        #--------------------------------------
        self.setup_layout()
        self.setup_connections()

    def setup_layout(self) -> None:
        parent = self.parent

        grid_load = create_grid_from_list(parent, [
            (self.base_f06_directory_label, self.base_f06_directory_edit, self.base_f06_directory_browse, self.base_f06_directory_load),
            (self.excel_filename_label, self.excel_filename_edit, self.excel_filename_browse, self.load_excel_button),
            (self.tab_label, self.tab_select_pulldown), # self.tab_edit
        ])

        grid2 = QGridLayout(parent)
        file_row = 1
        grid2.addWidget(self.word_filename_label, file_row, 0)
        grid2.addWidget(self.word_filename_edit, file_row, 1)
        grid2.addWidget(self.word_filename_browse, file_row, 2)

        #-----------------------------------------------------------
        # grid3 = create_grid_from_list(parent, [
        #     (self.base_f06_directory_label, self.base_f06_directory_edit, self.base_f06_directory_browse, self.base_f06_directory_load),
        # ])
        grid3 = QGridLayout(parent)
        file_row = 1
        grid3.addWidget(self.xaxis_label, file_row, 0)
        grid3.addWidget(self.xaxis_pulldown, file_row, 1)
        file_row += 1
        grid3.addWidget(self.yaxis_label, file_row, 0)
        grid3.addWidget(self.yaxis_pulldown, file_row, 1)
        file_row += 1
        grid3.addWidget(self.config_label, file_row, 0)
        grid3.addWidget(self.config_edit, file_row, 1)
        file_row += 1
        #-----------------------------------------------------------

        # vbox = QVBoxLayout()
        self.addLayout(grid_load)
        # vbox.addWidget(self.load_excel_button)
        self.addWidget(self.table_widget)
        self.addStretch()
        self.addLayout(grid3)
        self.addLayout(grid2)
        self.addWidget(self.run_organize_button)
        self.addWidget(self.progress_bar)

    def setup_connections(self) -> None:
        parent = self.parent
        self.load_excel_button.clicked.connect(self.on_load_excel)
        self.tab_select_pulldown.currentIndexChanged.connect(self.on_select_excel_tab)

        # self.excel_filename_browse.clicked.connect(self.on_select_excel_tab)
        self.base_f06_directory_load.clicked.connect(self.on_base_f06_directory_load)
        self.base_f06_directory_browse.clicked.connect(self.on_base_f06_directory_browse)
        self.excel_filename_browse.clicked.connect(self.on_load_excel_file)
        self.word_filename_browse.clicked.connect(self.on_load_word_file)

        # self.base_f06_directory_browse.setEnabled(False)
        self.tab_select_pulldown.setEnabled(False)
        # self.run_organize_button.setEnabled(False)
        # self.excel_filename_browse.setEnabled(False)
        self.run_organize_button.clicked.connect(self.on_run)
        self.base_f06_directory_edit.setText(self.base_f06_directory)

    def get_settings(self) -> dict[str, Any]:
        parent = self.parent

        print(f'data = {parent.data}')
        # vl = self.vl
        # vl = self.data['vl']
        # data = {
        #     # 'bdf_filename': self.bdf_filename,
        #     # 'op2_filename': self.op2_filename,
        #     'log_scale_x': self.log_xscale_checkbox.isChecked(),
        #     'log_scale_y1': self.log_yscale1_checkbox.isChecked(),
        #     'log_scale_y2': self.log_yscale2_checkbox.isChecked(),
        #     'use_rhoref': self.use_rhoref,
        #     'show_points': self.show_points,
        #     'show_mode_number': self.show_mode_number,
        #     'show_detailed_mode_info': self.show_detailed_mode_info,
        #     'point_spacing': self.point_spacing,
        #     'show_lines': self.show_lines,
        #
        #     'recent_files': self.recent_files,
        #     'subcase': subcase,
        #     # 'modes': modes,
        #     'selected_modes': selected_modes,
        #     'x_plot_type': self.x_plot_type,
        #     'plot_type': self.plot_type,
        #     'index_lim': index_lim,
        #     'eas_lim': eas_lim,
        #     'tas_lim': tas_lim,
        #     'mach_lim': mach_lim,
        #     'alt_lim': alt_lim,
        #     'q_lim': q_lim,
        #     'rho_lim': rho_lim,
        #     'ikfreq_lim': ikfreq_lim,
        #
        #     'damp_lim': ydamp_lim,
        #     'freq_lim': freq_lim,
        #     'kfreq_lim': kfreq_lim,
        #     'eigr_lim': eigr_lim,
        #     'eigi_lim': eigi_lim,
        #     'output_directory': output_directory,
        #     'units_in': units_in,
        #     'units_out': units_out,
        #     'freq_tol': freq_tol,
        #     'freq_tol_remove': freq_tol_remove,
        #     'mag_tol': mag_tol,
        #     'vl': vl,
        #     'vf': vf,
        #     'damping': damping,
        #     'damping_required': damping_required,
        #     'damping_required_tol': damping_required_tol,
        #     'eas_flutter_range': eas_flutter_range,
        #     'point_removal': point_removal,
        #     'mode_switch_method': self.mode_switch_method,
        # }
        #
        # buggy?
        # x_plot_type = self.x_plot_type
        # if x_plot_type == 'index':
        #     xlim = self.index_lim
        # elif x_plot_type == 'eas':
        #     xlim = self.eas_lim
        # elif x_plot_type == 'tas':
        #     xlim = self.tas_lim
        # elif x_plot_type == 'mach':
        #     xlim = self.mach_lim
        # elif x_plot_type == 'alt':
        #     xlim = self.alt_lim
        # elif x_plot_type == 'q':
        #     xlim = self.q_lim
        # elif x_plot_type == 'kfreq':
        #     xlim = self.kfreq_lim
        # elif x_plot_type == 'ikfreq':
        #     xlim = self.ikfreq_lim
        # else:  # pragma: no cover
        #     log.error(f'x_plot_type={x_plot_type!r} is not supported')
        #     # raise RuntimeError(x_plot_type)
        #     xlim = (None, None)

        modes = None if len(parent.selected_modes) == 0 else parent.selected_modes
        settings = {
            'x_plot_type': 'eas',
            'nrigid_body_modes': 6,  # TODO: 6
            'f06_units': parent._units_in,
            'out_units': parent._units_out,
            'modes': modes,
            'vl_target': float(parent.data['vl']),
            'vf_target': float(parent.data['vf']),
            #'xlim_kfreq': str_limit_to_limit(self.kfreq_lim),
            #'ylim_damping': self.damping_lim,  # NO

            # self.index_lim = index_lim
            # self.tas_lim = tas_lim
            # self.mach_lim = mach_lim
            # self.alt_lim = alt_lim
            # self.q_lim = q_lim
            # self.rho_lim = rho_lim
            # self.ikfreq_lim = ikfreq_lim
            # self.ydamp_lim = ydamp_lim
            # self.kfreq_lim = kfreq_lim
            # self.freq_lim = freq_lim
            # self.eigi_lim = eigi_lim
            # self.eigr_lim = eigr_lim
            'ylim_damping': str_limit_to_limit(parent.ydamp_lim),
            'ylim_freq': str_limit_to_limit(parent.freq_lim),
            'eas_lim': str_limit_to_limit(parent.eas_lim),
            'freq_tol': float(parent.freq_tol),
            'freq_tol_remove': float(parent.freq_tol_remove),
            'damping_required': float(parent.damping_required),
            'damping_required_tol': double_or_blank(parent.damping_required_tol, default=0.0),
            'damping_limit': float(parent.damping),  # % damping
            'eas_flutter_range': str_limit_to_limit(parent.eas_flutter_range),
            'plot_font_size': parent.plot_font_size,
            'show_lines': parent.show_lines,
            'show_points': parent.show_points,
            'show_mode_number': parent.show_mode_number,
            'show_detailed_mode_info': parent.show_detailed_mode_info,
            'point_spacing': parent.point_spacing,
            'use_rhoref': parent.use_rhoref,
            'flutter_ncolumns': parent.flutter_ncolumns,
            # 'mode_switch_method': None,
            #------------------
            'divergence_legend_loc': parent.divergence_legend_loc,
            'divergence_freq_tol': parent.freq_divergence_tol,
            'flutter_bbox_to_anchor_x': parent.flutter_bbox_to_anchor_x,
            'freq_ndigits': parent.freq_ndigits,
        }
        return settings

    @dontcrash
    def on_run(self) -> None:
        is_valid = self.parent.validate()
        log = self.log
        if not is_valid:
            return

        #----------------------------------------------
        out_table = self.table_widget.get_data()
        # print(f'out_table:\n{out_table}')
        is_passed, configs = get_configs(
            out_table, self.config_edit.text(), log)
        # if len(configs) == 0:
        #     log.error('no configs...')
        #     return
        assert isinstance(configs, list), configs
        # assert isinstance(configs[0], str), configs
        #----------------------------------------------
        # configs = [config if config.strip() else for config in ]

        word_filename = self.word_filename_edit.text().strip()
        if len(word_filename) == 0:
            word_filename = 'flutter_summary.docx'
        if not word_filename.lower().endswith('.docx'):
            word_filename += 'docx'

        # docx_filename = dirname / f'{prefix}flutter.docx'
        # from pyNastran.utils import object_attributes
        # print(object_attributes(self))

        # print('settings')
        settings = self.get_settings()

        f06_filenames = out_table['Filename'].to_list()
        # print(f'f06_filenames = {f06_filenames}')
        if len(configs) == 0:
            configs = [os.path.splitext(os.path.basename(fname))[0]
                       for fname in f06_filenames]

        # print('make settings')
        print(f'settings = {settings}')

        configs, f06_filenames = remove_empty_rows(
            configs, f06_filenames, log)
        # log.info(f'f06_filenames2 = {f06_filenames}')
        # log.info(f'configs2 = {configs}')

        nfiles = len(f06_filenames)
        if nfiles == 0:
            log.error('no files found')
            return

        # Show progress bar and disable button
        self.progress_bar.setVisible(True)
        self.progress_bar.setMaximum(len(f06_filenames))
        self.progress_bar.setValue(0)
        self.run_organize_button.setEnabled(False)

        # Process files
        try:
            log.info(f'Processing {len(f06_filenames)} files...')
            write_report(
                word_filename,
                f06_filenames, configs, out_table,
                log, settings,
                progress_callback=self.update_organize_progress,
                **settings)
            log.info(f'Successfully created {word_filename}')
        except Exception as e:
            log.error(f'Failed to create Word document: {str(e)}')
            log.error(traceback.format_exc())
        finally:
            # Hide progress bar and re-enable button
            self.progress_bar.setVisible(False)
            self.run_organize_button.setEnabled(True)

    def update_organize_progress(self, current: int, total: int):
        """Update progress bar"""
        self.progress_bar.setValue(current)
        QApplication.processEvents()  # Keep GUI responsive

    def on_base_f06_directory_load(self) -> None:
        is_passed, directory = get_file_edit(
            'base_f06_directory', self.base_f06_directory_edit,
            self.log, allow_empty=True)
        if not is_passed:
            return
        if not os.path.isdir(directory):
            self.log.error(f'{directory!r} is not a directory')
            return
        filenames = get_files_of_type(directory, '.f06')
        if len(filenames) == 0:
            self.log.error(f'no .f06 files were found in {directory!r}')
            return

        # TODO: parse the filename for floats
        headers = ['Filename']
        # filenames2 = [os.path.relpath(pathi, directory) for pathi in natsort.natsorted(filenames)]
        filenames2 = list(natsort.natsorted(filenames))

        if len(filenames2) == 1:
            # can't do a trade study
            data_table = [filenames2]
        else:
            # remove the extension
            base_filenames = [os.path.splitext(os.path.basename(filename))[0] for filename in filenames2]
            data_table_initial = split_by_pattern(base_filenames)

            # define the headers
            data_row0 = data_table_initial[0]
            headers = [f'{icol+1:d}' for icol in range(len(data_row0))] + ['Filename']

            # add the filename onto the end
            data_table = [line + [filename] for line, filename in zip(data_table_initial, filenames2)]
        self.table_widget.load_table_data(headers, data_table)

    @dontcrash
    def on_load_excel(self):
        is_passed, excel_filename = get_file_edit(
            'excel_filename', self.excel_filename_edit,
            self.log)
        if not is_passed:
            return
        self.excel_dict = pd.read_excel(excel_filename, sheet_name=None)

        self.table_widget.clear()
        keys = list(self.excel_dict)
        # self.log.debug(f'df keys = {keys}')

        if len(keys):
            key0 = keys[0]
            if self.tab_select_pulldown.count() > 0:
                self.tab_select_pulldown.clear()  # seems to cause crashes...
            self.tab_select_pulldown.addItems(keys)
            self.tab_edit.setText(key0)
            self.tab_select_pulldown.setEnabled(True)
            # self.on_select_excel_tab()

    def on_base_f06_directory_browse(self):
        start_path = self.base_f06_directory_edit.text()
        directory = self._on_load_directory(
            title='Select F06 Directory', start_path=start_path)
        if directory:
            self.base_f06_directory_edit.setText(directory)
            self.base_f06_directory = directory

    def on_load_excel_file(self) -> None:
        start_path = self.excel_filename_edit.text()
        filename = self._on_load_file(
            title='Select Excel File', start_path=start_path,
            file_filter='Excel File (*.xlsx);;All Files (*)')
        if filename:
            self.excel_filename_edit.setText(filename)
            self.excel_filename = filename

    def on_load_word_file(self) -> None:
        start_path = self.word_filename_edit.text()
        filename = self._on_load_file(
            title='Select Word File', start_path=start_path,
            file_filter='Word File (*.docx);;All Files (*)')
        if filename:
            self.word_filename_edit.setText(filename)
            self.word_filename = filename

    def _on_load_file(self,
                      title: str="Select File",
                      start_path: str="",
                      file_filter: str="Text Files (*.txt);;All Files (*)") -> str:
        """
        Open a dialog to select a directory

        Parameters
        ----------
        title : str
            Dialog window title
        start_path : str
            Initial directory path

        Returns
        -------
        file_path : str
            Selected directory path, or empty string if canceled
        """
        file_path, _ = QFileDialog.getOpenFileName(
            self.parent, title, start_path, file_filter)
        if not file_path:
            return ''
        return file_path

    def _on_load_directory(self,
                           title: str="Select Directory",
                           start_path: str="") -> str:
        """
        Open a dialog to select a directory

        Parameters
        ----------
        title : str
            Dialog window title
        start_path : str
            Initial directory path

        Returns
        -------
        directory : str
            Selected directory path, or empty string if canceled
        """
        directory = QFileDialog.getExistingDirectory(
            self.parent, title, start_path,
            QFileDialog.Option.ShowDirsOnly)
        if directory:
            # print(f"Selected directory: {directory}")
            return directory
        else:
            # print("No directory selected")
            return ""

    def on_select_excel_tab(self) -> None:
        if self.tab_select_pulldown.count() == 0:
            return
        key = self.tab_select_pulldown.currentText()
        df = self.excel_dict[key]
        headers = df.columns
        data_table = df.to_numpy()
        self.table_widget.load_table_data(headers, data_table)

    def get_dict(self) -> dict[str, str]:
        excel_filename = self.excel_filename_edit.text().strip()
        base_f06_directory = self.base_f06_directory_edit.text().strip()
        word_filename = self.word_filename_edit.text().strip()

        configs = self.config_edit.text().strip(' ,')
        xaxis = self.xaxis_pulldown.currentText()
        yaxis = self.yaxis_pulldown.currentText()
        # self.word_filename = word_filename
        trade = {
            'excel_filename': excel_filename,
            'base_f06_directory': base_f06_directory,
            'word_filename': word_filename,
            'configs': configs,
            'xaxis': xaxis,
            'yaxis': yaxis,
        }
        return trade

    def get_load_settings(self, key: str,
                          prefix: str='') -> list[tuple[str, int, QLineEdit]]:
        if key == 'lineedit':
            out = [
                # should be in trade/
                (f'{prefix}excel_filename', -1, self.excel_filename_edit),
                (f'{prefix}base_f06_directory', -1, self.base_f06_directory_edit),
                (f'{prefix}word_filename', -1, self.word_filename_edit),
                (f'{prefix}configs', -1, self.config_edit)
            ]
        elif key == 'combobox':
            out = [
                ('xaxis', self.xaxis_pulldown, AXIS_VALUES),
                ('yaxis', self.yaxis_pulldown, AXIS_VALUES),
            ]
        else:  # pragma: no cover
            raise NotImplementedError(key)
        return out


def str_limit_to_limit(data: list[str | None]) -> Limit:
    data_out = []
    for value in data:
        value = value.strip()
        if len(value) == 0:
            data_out.append(None)
        else:
            data_out.append(float(value))
    return data_out

def double_or_blank(value: float | str,
                    default: float=0.0) -> float:
    if isinstance(value, str):
        value = value.strip()
        if len(value) == 0:
            return default
        value2 = float(value)
        return value2
    return float(value)


def remove_empty_rows(configs: list[str],
                      f06_filenames: list[str],
                      log: SimpleLogger) -> tuple[list[str], list[str]]:
    configs2 = []
    f06_filenames2 = []
    for config, f06_filename in zip(configs, f06_filenames):
        if not os.path.exists(f06_filename):
            log.warning(print_bad_path(f06_filename))
            continue
        configs2.append(config)
        f06_filenames2.append(f06_filename)
    configs = configs2
    f06_filenames = f06_filenames2
    return configs, f06_filenames


def split_by_pattern(strings: list[str],
                     delimiter: str='_',
                     group_common: bool=True):
    """
    Main function to split strings based on common pattern.
    """
    if not strings:
        return []

    # Split all strings
    split_lists = [s.split(delimiter) for s in strings]
    if not group_common or len(split_lists) < 2:
        return split_lists
    # split_list0 = split_lists[0]
    # print(f'split_list0 = {split_list0}')

    # Find common prefix
    common_prefix_len = 0
    min_length = min(len(parts) for parts in split_lists)

    for i in range(min_length):
        if all(parts[i] == split_lists[0][i] for parts in split_lists):
            common_prefix_len += 1
            continue
        break
    # print(f'common_prefix_len = {common_prefix_len}')
    # prefix = split_list0[:common_prefix_len]
    # print(f'prefix = {prefix}')

    # Find common suffix
    common_suffix_len = 0
    for i in range(1, min_length - common_prefix_len + 1):
        if all(parts[-i] == split_lists[0][-i] for parts in split_lists):
            common_suffix_len += 1
            continue
        break
    # print(f'common_suffix_len = {common_suffix_len}')

    # Reconstruct
    result = []
    for parts in split_lists:
        new_parts = []
        if common_prefix_len > 0:
            new_parts.append(delimiter.join(parts[:common_prefix_len]))

        middle_start = common_prefix_len
        middle_end = len(parts) - common_suffix_len if common_suffix_len > 0 else len(parts)
        new_parts.extend(parts[middle_start:middle_end])

        if common_suffix_len > 0:
            new_parts.append(delimiter.join(parts[-common_suffix_len:]))
        result.append(new_parts)
    return result


def get_configs(out_table: pd.DataFrame,
                config_text: str,
                log: SimpleLogger) -> tuple[bool, list[str]]:
    is_passed = True
    columns = out_table.columns
    config_keys = config_text.strip(', ').split(',')
    for config_key in config_keys:
        if config_key not in columns:
            log.warning(f'key={config_key!r} is not a column in the case table')

    config_headers_lower = [column.lower().strip() for column in columns]
    config_headers = [column.strip() for column in columns
                      if 'config' in config_headers_lower]

    if 'Filename' not in columns:
        log.error('missing Filename from case table')
        is_passed = False

    configs = []
    if 'config' in config_headers_lower:
        iconfig_key = config_headers_lower.index('config')
        configs = out_table[iconfig_key].to_list()
        config_headers.remove('config')
    # elif is_passed:
    #     filenames = out_table['Filename'].tolist()
    #     configs = [os.path.splitext(os.path.basename(filenames))[0]
    #                for filename in filenames]
    else:
        log.error('Missing Config from case table')
        is_passed = False

    return is_passed, configs

def get_file_edit(name: str,
                  filename_edit: QLineEdit,
                  log: SimpleLogger,
                  allow_empty: bool=False) -> tuple[bool, str]:
    is_passed = False
    filename = filename_edit.text().strip()
    if len(filename) == 0:
        if allow_empty:
            is_passed = True
            return is_passed, '.'
        else:
            log.error(f'{name}={filename!r} is empty')
            return is_passed, filename

    if not os.path.exists(filename):
        log.error(print_bad_path(filename))
        return is_passed, ''
    log.info(f'loading {name} {filename}')
    is_passed = True
    return is_passed, filename
