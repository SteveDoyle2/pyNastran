from qtpy.QtWidgets import (
    QLabel,
    #QWidget,
    #QApplication, QMenu, QVBoxLayout, QLineEdit, QComboBox,
    #QHBoxLayout, QPushButton,
    QGridLayout,
    #QAction,
    QCheckBox,
    ##QRadioButton,
    #QListWidgetItem, QAbstractItemView,
    #QListWidget,
    QSpinBox,
)
from pyNastran.gui.utils.qt.pydialog import PyDialog, make_font, check_color

class FlutterPreferencesDialog(PyDialog):
    def __init__(self, data, win_parent=None):
        """
        Saves the data members from data and
        performs type checks
        """
        PyDialog.__init__(self, data, win_parent)
        #print(f'data = {data}')
        self.setup_widgets(data)
        self.on_font_size()
        self.setup_layout()
        self.setup_connections()

    def setup_widgets(self, data: dict[str, bool]) -> None:
        self.font_size_label = QLabel('Font Size')
        self.font_size_edit = QSpinBox()
        self.font_size_edit.setValue(data['font_size'])
        #self.font_size_edit.setValue(data['font_size'])

        self.plot_font_size_label = QLabel('Plot Font Size')
        self.plot_font_size_edit = QSpinBox()
        self.plot_font_size_edit.setValue(data['plot_font_size'])

        self.export_png_checkbox = QCheckBox('Export PNG')
        self.export_csv_checkbox = QCheckBox('Export CSV')
        self.export_f06_checkbox = QCheckBox('Export F06')
        self.export_zona_checkbox = QCheckBox('Export Zona')

        self.export_png_checkbox.setChecked(data['export_to_png'])
        self.export_csv_checkbox.setChecked(data['export_to_csv'])
        self.export_f06_checkbox.setChecked(data['export_to_f06'])
        self.export_zona_checkbox.setChecked(data['export_to_zona'])

    def setup_layout(self) -> None:
        irow = 0
        grid = QGridLayout()
        grid.addWidget(self.font_size_label, irow, 0)
        grid.addWidget(self.font_size_edit, irow, 1)
        irow += 1
        grid.addWidget(self.plot_font_size_label, irow, 0)
        grid.addWidget(self.plot_font_size_edit, irow, 1)
        irow += 1

        grid.addWidget(self.export_png_checkbox, irow, 0)
        irow += 1
        grid.addWidget(self.export_csv_checkbox, irow, 0)
        grid.addWidget(self.export_f06_checkbox, irow, 1)
        grid.addWidget(self.export_zona_checkbox, irow, 2)
        irow += 1
        grid.setColumnStretch(grid.columnCount(), 2)
        #widget = QWidget()
        #widget.setLayout(vbox2)
        #self.setCentralWidget(widget)
        self.setLayout(grid)

    def setup_connections(self) -> None:
        self.export_csv_checkbox.clicked.connect(self.on_export_csv)
        self.export_f06_checkbox.clicked.connect(self.on_export_f06)
        self.export_png_checkbox.clicked.connect(self.on_export_png)
        self.export_zona_checkbox.clicked.connect(self.on_export_zona)
        self.font_size_edit.valueChanged.connect(self.on_font_size)
        self.plot_font_size_edit.valueChanged.connect(self.on_plot_font_size)

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
