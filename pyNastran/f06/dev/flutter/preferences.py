
"""
TODO: change from dt_ms to FPS
"""
from qtpy.QtWidgets import (
    QLabel,
    # QWidget,
    # QApplication, QMenu, QVBoxLayout, QLineEdit,
    # QHBoxLayout, QPushButton,
    QGridLayout,
    # QAction,
    QCheckBox, QComboBox,
    # QRadioButton,
    # QListWidgetItem, QAbstractItemView,
    # QListWidget,
    QSpinBox, QDoubleSpinBox,
)
from pyNastran.gui.utils.qt.pydialog import PyDialog, make_font, check_color
from pyNastran.gui.utils.qt.qcombobox import set_combo_box_text

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

        self.flutter_ncolumns_label = QLabel('Flutter nColumns:')
        self.flutter_ncolumns_spinner = QSpinBox()
        self.flutter_ncolumns_spinner.setMinimum(0)
        self.flutter_ncolumns_spinner.setMaximum(3)

        flutter_ncolumns = 0 if data['flutter_ncolumns'] is None else data['flutter_ncolumns']
        self.flutter_ncolumns_spinner.setValue(flutter_ncolumns)

        self.export_png_checkbox = QCheckBox('Export PNG')
        self.export_csv_checkbox = QCheckBox('Export CSV')
        self.export_f06_checkbox = QCheckBox('Export F06')
        self.export_zona_checkbox = QCheckBox('Export Zona')

        self.export_png_checkbox.setChecked(data['export_to_png'])
        self.export_csv_checkbox.setChecked(data['export_to_csv'])
        self.export_f06_checkbox.setChecked(data['export_to_f06'])
        self.export_zona_checkbox.setChecked(data['export_to_zona'])

        if not self.use_vtk:
            objs = [self.icase_label, self.icase_edit,
                    self.nphase_label, self.nphase_edit,
                    self.dt_ms_label, self.dt_ms_edit, self.animate_checkbox]
            for obj in objs:
                obj.setVisible(False)

    def setup_layout(self) -> None:
        irow = 0
        grid = QGridLayout()
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
        # widget = QWidget()
        # widget.setLayout(vbox2)
        # self.setCentralWidget(widget)
        self.setLayout(grid)

    def setup_connections(self) -> None:
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
        self.divergence_legend_loc_combobox.currentIndexChanged.connect(self.on_divergence_legend_loc)

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

    def on_export_zona(self) -> None:
        self.win_parent.export_to_zona = self.export_zona_checkbox.isChecked()

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
