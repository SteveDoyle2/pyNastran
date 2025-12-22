
"""
TODO: change from dt_ms to FPS
"""
from PyQt5.QtWidgets import QPushButton, QVBoxLayout
from qtpy.QtWidgets import (
    QLabel,
    # QWidget,
    # QApplication, QMenu, QVBoxLayout, QLineEdit,
    QGridLayout,
    # QAction,
    QCheckBox, QComboBox,
    # QRadioButton,
    # QListWidgetItem, QAbstractItemView,
    # QListWidget,
    QSpinBox, QDoubleSpinBox,
)
from pyNastran.gui.utils.qt.pydialog import PyDialog, QFloatEdit, make_font, check_color
from pyNastran.gui.utils.qt.qcombobox import set_combo_box_text

FONT_SIZE_DEFAULT = 10
LEGEND_LOC_DEFAULT = 'best'
FLUTTER_BBOX_TO_ANCHOR_DEFAULT = 1.02
FLUTTER_NCOLUMNS_DEFAULT = 0
FREQ_NDIGITS_DEFAULT = 1
FREQ_DIVERGENCE_TOL = 0.05

LEGEND_LOCS = [
    'best', 'none',
    'upper right', 'upper center', 'upper left',
    'center left', 'center', 'center right',
    'lower left', 'lower center', 'lower right',
]


class FlutterPreferencesDialog(PyDialog):
    def __init__(self, data, gui_obj, win_parent=None):
        """
        Saves the data members from data and
        performs type checks
        """
        PyDialog.__init__(self, data, win_parent)
        self.gui_obj = gui_obj
        self.use_vtk = gui_obj.use_vtk
        # print(f'data = {data}')
        self.setup_widgets(data)
        self.on_font_size()
        self.setup_layout()
        self.setup_connections()
        self.setWindowTitle('Preferences')
        self.show()

    def setup_widgets(self, data: dict[str, bool]) -> None:
        self.font_size_label = QLabel('Font Size')
        self.font_size_edit = QSpinBox()
        self.font_size_edit.setValue(data['font_size'])
        # self.font_size_edit.setValue(data['font_size'])

        self.plot_font_size_label = QLabel('Plot Font Size')
        self.plot_font_size_edit = QSpinBox()
        self.plot_font_size_edit.setValue(data['plot_font_size'])

        self.auto_update_checkbox = QCheckBox('Auto_Update')
        self.auto_update_checkbox.setChecked(data['auto_update'])

        self.nphase_label = QLabel('Num Phase:')
        self.nphase_edit = QSpinBox()
        self.nphase_edit.setValue(data['nphase'])
        self.nphase_edit.setMinimum(4)
        self.nphase_edit.setMaximum(51)

        self.icase_label = QLabel('iCase:')
        self.icase_edit = QSpinBox()
        self.icase_edit.setValue(data['icase'])
        self.icase_edit.setMinimum(0)
        # self.nphase_edit.setMaximum(51)

        # self.nframes_label = QLabel('Animation Update Time (ms):')
        # self.nframes_edit = QSpinBox()
        # self.nframes_edit.setValue(data['dt_ms'])
        # self.nframes_edit.setMinimum(100)
        # self.nframes_edit.setMaximum(60)

        # Update every N milliseconds
        self.dt_ms_label = QLabel('Animation Update Time (ms):')
        self.dt_ms_edit = QSpinBox()
        self.dt_ms_edit.setValue(data['dt_ms'])
        self.dt_ms_edit.setMinimum(100)
        self.dt_ms_edit.setMaximum(5000)
        self.animate_checkbox = QCheckBox('Animate')
        self.animate_checkbox.setChecked(data['animate'])
        self.animate_checkbox.setEnabled(False)
        # self.dt_ms_edit.setMaximum(2000)

        self.divergence_legend_loc_label = QLabel('Divergence Legend Location:')
        self.divergence_legend_loc_combobox = QComboBox()
        self.divergence_legend_loc_combobox.addItems(LEGEND_LOCS)
        set_combo_box_text(self.divergence_legend_loc_combobox, data['divergence_legend_loc'])

        self.flutter_bbox_to_anchor_x_label = QLabel('Flutter BBox to Anchor:')
        self.flutter_bbox_to_anchor_x_spinner = QDoubleSpinBox()
        self.flutter_bbox_to_anchor_x_spinner.setValue(data['flutter_bbox_to_anchor_x'])
        self.flutter_bbox_to_anchor_x_spinner.setMinimum(1.0)
        self.flutter_bbox_to_anchor_x_spinner.setMaximum(2.0)
        self.flutter_bbox_to_anchor_x_spinner.setSingleStep(0.01)

        flutter_ncolumns = _get_dict_value(data, 'flutter_ncolumns', 0)
        self.flutter_ncolumns_label = QLabel('Flutter nColumns:')
        self.flutter_ncolumns_spinner = QSpinBox()
        self.flutter_ncolumns_spinner.setMinimum(0)
        self.flutter_ncolumns_spinner.setMaximum(3)
        self.flutter_ncolumns_spinner.setValue(flutter_ncolumns)

        freq_ndigits = _get_dict_value(data, 'freq_ndigits', FREQ_NDIGITS_DEFAULT)
        self.freq_ndigits_label = QLabel('Freq Digits:')
        self.freq_ndigits_spinner = QSpinBox()
        self.freq_ndigits_spinner.setMinimum(0)
        self.freq_ndigits_spinner.setMaximum(4)
        self.freq_ndigits_spinner.setValue(freq_ndigits)

        freq_diveregence_tol = _get_dict_value(data, 'freq_divergernce_tol', str(FREQ_DIVERGENCE_TOL))
        self.freq_divergence_tol_label = QLabel('Divergence Freq Tol (Hz):')
        self.freq_divergence_tol_edit = QFloatEdit('0', self)
        self.freq_divergence_tol_edit.setText(freq_diveregence_tol)
        self.freq_divergence_tol_label.setEnabled(False)
        self.freq_divergence_tol_edit.setEnabled(False)

        self.export_png_checkbox = QCheckBox('Export PNG')
        self.export_csv_checkbox = QCheckBox('Export CSV')
        self.export_f06_checkbox = QCheckBox('Export F06')
        self.export_zona_checkbox = QCheckBox('Export Zona')

        self.export_png_checkbox.setChecked(data['export_to_png'])
        self.export_csv_checkbox.setChecked(data['export_to_csv'])
        self.export_f06_checkbox.setChecked(data['export_to_f06'])
        self.export_zona_checkbox.setChecked(data['export_to_zona'])

        self.default_button = QPushButton('Default')

        if not self.use_vtk:
            objs = [self.icase_label, self.icase_edit,
                    self.nphase_label, self.nphase_edit,
                    self.dt_ms_label, self.dt_ms_edit, self.animate_checkbox]
            for obj in objs:
                obj.setVisible(False)

    def setup_layout(self) -> None:
        irow = 0
        grid = QGridLayout()
        grid.addWidget(self.auto_update_checkbox, irow, 0)
        irow += 1
        grid.addWidget(self.font_size_label, irow, 0)
        grid.addWidget(self.font_size_edit, irow, 1)
        irow += 1
        grid.addWidget(self.plot_font_size_label, irow, 0)
        grid.addWidget(self.plot_font_size_edit, irow, 1)
        irow += 1

        grid.addWidget(self.icase_label, irow, 0)
        grid.addWidget(self.icase_edit, irow, 1)
        irow += 1
        grid.addWidget(self.nphase_label, irow, 0)
        grid.addWidget(self.nphase_edit, irow, 1)
        irow += 1
        grid.addWidget(self.dt_ms_label, irow, 0)
        grid.addWidget(self.dt_ms_edit, irow, 1)
        grid.addWidget(self.animate_checkbox, irow, 2)
        irow += 1

        grid.addWidget(self.flutter_bbox_to_anchor_x_label, irow, 0)
        grid.addWidget(self.flutter_bbox_to_anchor_x_spinner, irow, 1)
        irow += 1
        grid.addWidget(self.flutter_ncolumns_label, irow, 0)
        grid.addWidget(self.flutter_ncolumns_spinner, irow, 1)
        irow += 1
        grid.addWidget(self.freq_ndigits_label, irow, 0)
        grid.addWidget(self.freq_ndigits_spinner, irow, 1)
        irow += 1

        grid.addWidget(self.freq_divergence_tol_label, irow, 0)
        grid.addWidget(self.freq_divergence_tol_edit, irow, 1)
        irow += 1
        grid.addWidget(self.divergence_legend_loc_label, irow, 0)
        grid.addWidget(self.divergence_legend_loc_combobox, irow, 1)
        irow += 1

        grid.addWidget(self.export_png_checkbox, irow, 0)
        irow += 1
        grid.addWidget(self.export_csv_checkbox, irow, 0)
        grid.addWidget(self.export_f06_checkbox, irow, 1)
        grid.addWidget(self.export_zona_checkbox, irow, 2)
        irow += 1
        grid.setColumnStretch(grid.columnCount(), 2)

        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        vbox.addWidget(self.default_button)
        # widget = QWidget()
        # widget.setLayout(vbox2)
        # self.setCentralWidget(widget)
        self.setLayout(vbox)

    def setup_connections(self) -> None:
        self.default_button.clicked.connect(self.on_default)
        self.export_csv_checkbox.clicked.connect(self.on_export_csv)
        self.export_f06_checkbox.clicked.connect(self.on_export_f06)
        self.export_png_checkbox.clicked.connect(self.on_export_png)
        self.export_zona_checkbox.clicked.connect(self.on_export_zona)
        self.font_size_edit.valueChanged.connect(self.on_font_size)
        self.plot_font_size_edit.valueChanged.connect(self.on_plot_font_size)
        #self.animate_checkbox.clicked.connect(self.on_animate)
        self.nphase_edit.valueChanged.connect(self.on_nphase)
        self.dt_ms_edit.valueChanged.connect(self.on_dt_ms)
        self.icase_edit.valueChanged.connect(self.on_icase)

        self.flutter_bbox_to_anchor_x_spinner.valueChanged.connect(self.on_flutter_bbox_to_anchor_x)
        self.flutter_ncolumns_spinner.valueChanged.connect(self.on_flutter_ncolumns)
        self.freq_ndigits_spinner.valueChanged.connect(self.on_freq_ndigits)
        self.divergence_legend_loc_combobox.currentIndexChanged.connect(self.on_divergence_legend_loc)
        self.freq_divergence_tol_edit.textChanged.connect(self.on_freq_divergence_tol)
        self.auto_update_checkbox.clicked.connect(self.on_auto_update)

    def on_font_size(self) -> None:
        font_size = self.font_size_edit.value()
        assert font_size > 0, font_size
        font = make_font(font_size, is_bold=False)
        self.setFont(font)
        self.win_parent.on_set_font_size(font_size)

    def on_plot_font_size(self) -> None:
        plot_font_size = self.plot_font_size_edit.value()
        self.win_parent.plot_font_size = plot_font_size

    def on_export_f06(self) -> None:
        self.win_parent.export_to_f06 = self.export_f06_checkbox.isChecked()

    def on_export_csv(self) -> None:
        self.win_parent.export_to_csv = self.export_csv_checkbox.isChecked()

    def on_export_png(self) -> None:
        self.win_parent.export_to_png = self.export_png_checkbox.isChecked()

    def on_export_zaero(self) -> None:
        self.win_parent.export_to_zaero = self.export_zaero_checkbox.isChecked()

    def on_default(self) -> None:
        self.font_size_edit.setValue(FONT_SIZE_DEFAULT)
        self.plot_font_size_edit.setValue(FONT_SIZE_DEFAULT)
        self.flutter_bbox_to_anchor_x_spinner.setValue(FLUTTER_BBOX_TO_ANCHOR_DEFAULT)
        self.flutter_ncolumns_spinner.setValue(FLUTTER_NCOLUMNS_DEFAULT)
        self.freq_ndigits_spinner.setValue(FREQ_NDIGITS_DEFAULT)
        set_combo_box_text(self.divergence_legend_loc_combobox, LEGEND_LOC_DEFAULT)

        self.export_png_checkbox.setChecked(True)
        self.export_f06_checkbox.setChecked(False)
        self.export_csv_checkbox.setChecked(False)
        self.export_zona_checkbox.setChecked(False)
        # self.icase_edit
        # self.nphase_edit
        # self.dt_ms_edit
        # self.animate_checkbox


    def on_dt_ms(self) -> None:
        """TODO: move this behind an apply button"""
        value = self.dt_ms_edit.value()
        self.gui_obj.on_dt_ms(value)

    def on_nphase(self) -> None:
        """TODO: move this behind an apply button"""
        value = self.nphase_edit.value()
        self.gui_obj.on_nphase(value)

    def on_icase(self) -> None:
        """TODO: move this behind an apply button"""
        value = self.icase_edit.value()
        self.gui_obj.on_icase(value)
    # 'animate': True,

    def on_flutter_ncolumns(self) -> None:
        """TODO: move this behind an apply button"""
        value = self.flutter_ncolumns_spinner.value()
        # if value == 0:
        #     value = None
        self.gui_obj.on_flutter_ncolumns(value)

    def on_freq_ndigits(self) -> None:
        """TODO: move this behind an apply button"""
        value = self.freq_ndigits_spinner.value()
        self.gui_obj.on_freq_ndigits(value)

    def on_freq_divergence_tol(self) -> None:
        value_str = self.freq_divergence_tol_edit.text()
        try:
            value = float(value_str)
        except ValueError:
            return
        self.gui_obj.on_freq_divergence_tol(value)

    def on_auto_update(self) -> None:
        value = self.auto_update_checkbox.value()
        self.gui_obj.on_auto_update(value)

    def on_flutter_bbox_to_anchor_x(self) -> None:
        """TODO: move this behind an apply button"""
        value = self.flutter_bbox_to_anchor_x_spinner.value()
        self.gui_obj.on_flutter_bbox_to_anchor_x(value)

    def on_divergence_legend_loc(self) -> None:
        """TODO: move this behind an apply button"""
        value = self.divergence_legend_loc_combobox.currentText()
        self.gui_obj.on_divergence_legend_loc(value)

    def closeEvent(self, event):
        # event.accept()
        self.gui_obj.on_close()


def _get_dict_value(data: dict, name, default):
    value = data.get(name, default)
    if value is None:
        return default
    return value
