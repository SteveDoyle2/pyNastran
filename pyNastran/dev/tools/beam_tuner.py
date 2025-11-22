from qtpy.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QApplication, QMainWindow, QLineEdit, QGridLayout, QLabel, QWidget
from qtpy.QtWidgets import (
    QLabel, QWidget,
    QApplication, QVBoxLayout,  # QComboBox,  # QMenu, QLineEdit,
    QHBoxLayout, QPushButton, QGridLayout,
    # QAction,
    QCheckBox,
    QListWidgetItem, QAbstractItemView,
    QListWidget,  # QSpinBox, QTabWidget,  # QToolButton,
)
from PyQt5.QtGui import QFont

from pyNastran.gui.utils.qt.pydialog import QLineEdit, QFloatEdit, make_font

from pyNastran.bdf.field_writer_8 import print_float_8
from pyNastran.bdf.bdf import read_bdf, PBEAM, BDF, PSHELL, PCOMP
from typing import Any

FONT_SIZE = 12


class BeamTuner(QMainWindow):
    def __init__(self, data: dict[str, Any]):
        super().__init__()
        self.setWindowTitle('Beam Tuner')

        self.setup_properties()
        self.create_widgets()
        self.create_actions()
        self.setup_layout()

    def create_actions(self):
        self.ok_button.clicked.connect(self.on_read_model)
        # self.sidebar_widget.itemActivated.connect(self._on_update_prop_item)
        self.sidebar_widget.itemClicked.connect(self._on_update_prop_item)

    def create_widgets(self):
        bdf_filename = r'C:\work\pyNastran\models\bwb\bwb_saero.bdf'
        self.bdf_filename_label = QLabel('BDF Filename')
        self.bdf_filename_edit = QLineEdit(bdf_filename)
        self.bdf_filename_edit.setPlaceholderText('BDF Filename')
        self.ok_button = QPushButton('OK')
        self.close_button = QPushButton('Cancel')
        self.model = BDF()

    def setup_layout(self):
        grid = QGridLayout()
        irow = 1
        grid.addWidget(self.bdf_filename_label, irow, 0)
        grid.addWidget(self.bdf_filename_edit, irow, 1)

        #------------------------------
        ok_close = QHBoxLayout()
        ok_close.addWidget(self.ok_button)
        ok_close.addWidget(self.close_button)

        self.prop_grid = QGridLayout()
        self.prop_grid.addWidget(QLabel(''), 0, 0)

        #------------------------------
        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        vbox.addStretch(1)
        vbox.addLayout(self.prop_grid)
        vbox.addLayout(ok_close)
        #------------------------------
        hbox = QHBoxLayout()
        hbox.addWidget(self.sidebar_widget)
        #hbox.addStretch(1)
        hbox.addLayout(vbox)
        # self.setLayout(box)
        #------------------------------

        widget = QWidget(self)
        widget.setLayout(hbox)
        self.setCentralWidget(widget)

        font = QFont('Arial', FONT_SIZE)
        self.setFont(font)

        self.show()

    def on_read_model(self):
        bdf_filename = self.bdf_filename_edit.text()
        self.model = read_bdf(
            bdf_filename, validate=False, xref=False)

        properties = {}
        for pid, prop in sorted(self.model.properties.items()):
            line0 = prop.comment.split('\n')[0]
            comment = line0.split('$')[0].strip()
            if comment == '':
                comment = 'Prop'
            properties[pid] = (prop.type, comment)
        self._set_properties_table(self.sidebar_widget, properties)

    # def edit_property(self, pid: int) -> None:
    #     prop = self.model.properties[pid]
    #     if prop.type != 'PBEAM':
    #         return
    #     setup_pbeam(prop)

    def setup_properties(self) -> None:
        self.sidebar_widget = QListWidget(self)
        self.sidebar_widget.setMaximumWidth(200)  # was 100 when no freq
        self.sidebar_widget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        properties = {
            0: ('PSHELL', 'Prop'),
        }
        self._set_properties_table(self.sidebar_widget, properties)

    def _set_properties_table(self, sidebar_widget: QListWidget,
                              properties: dict[int, str]) -> None:
        sidebar_widget.clear()
        for pid, (pid_type, comment) in properties.items():
            # if pid == 0:
            msg = f'{pid_type} {pid}: {comment}'
            # else:
            #     msg = f'Prop {pid}; {pid_type}'
            prop_item = QListWidgetItem(msg)
            prop_item.setSelected(True)
            sidebar_widget.addItem(prop_item)

    def _on_update_prop_item(self):
        item = self.sidebar_widget.currentItem()
        #print(item)
        msg = item.text()

        type_pid = msg.split(':')[0].strip()
        pid = int(type_pid.split()[1])
        try:
            prop = self.model.properties[pid]
        except KeyError:
            print(f'couldnt find pid={pid}')
            return

        prop_type = prop.type
        if prop_type == 'PBEAM':
            func = setup_pbeam
        elif prop_type == 'PSHELL':
            func = setup_pshell
        elif prop_type == 'PCOMP':
            func = setup_pcomp
        else:
            print(f'Skipping {msg}')
            return
        print(prop_type, func)
        print('clearing grid...')
        # self.prop_grid.clear()
        clear_grid_layout(self.prop_grid)
        print('cleared grid...')
        print(func)
        func(self.prop_grid, prop)


def setup_pshell(grid: QGridLayout,
                 prop: PSHELL) -> QGridLayout:
    names = [
        ('Property ID', 'pid'),
        ('MID1', 'mid1'),
        ('MID2', 'mid2'),
        ('MID3', 'mid3'),
        ('MID4', 'mid4'),
        ('Thickness', 't',),
        ('NSM', 'nsm',),
        ('ts/t', 'tst',),
        ('12I/t^3', 'twelveIt3',),
        ('Z1', 'z1'),
        ('Z2', 'z2'),
    ]
    print(prop.get_stats())
    irow = _prop_names_to_table(grid, prop, names)
    return grid


def _prop_build_grid_array(grid: QGridLayout,
                           prop: PCOMP,
                           names: list[str],
                           names_write: list[str],
                           irow: int) -> int:
    datas = [getattr(prop, name) for name in names]
    for j, name in enumerate(names_write):
        grid.addWidget(QLabel(name), irow, j+1)
    irow += 1

    nlayers = len(datas[0])
    for ilayer in range(nlayers):
        label = QLabel(f'Layer {ilayer+1}')
        irowi = irow + ilayer
        grid.addWidget(label, irowi, 0)
        checkbox = QCheckBox()
        grid.addWidget(checkbox, irowi, 1)

    jcol = 2
    for datai in datas:
        # print('datai = ', jcol, datai)
        for irowi, val in enumerate(datai):
            irowii = irow+irowi
            # print(f'  i={irowi} j={jcol} val={val}')
            if isinstance(val, float):
                edit = QFloatEdit(_write_float(val))
            else:
                edit = QLineEdit(str(val))
            grid.addWidget(edit, irowii, jcol)
        jcol += 1

    irow += nlayers
    return irow


def _write_float(val: float) -> str:
    return print_float_8(val).strip()


def setup_pcomp(grid: QGridLayout, prop: PCOMP) -> None:
    print(prop.get_stats())
    names = [
        ('Property ID', 'pid'),
        # ('Thickness', 't',),
        ('NSM', 'nsm',),
        # ('ts/t', 'tst',),
        # ('12I/t^3', 'twelveIt3',),
        # ('TRef', 'Tref'),
        ('ge', 'ge'),
        ('sb', 'sb'),
        ('ft', 'ft'),
    ]
    # print(prop)
    irow = _prop_names_to_table(grid, prop, names)
    names2 = ['material_ids', 'thetas', 'thicknesses',  # 'souts'
              ]
    names_write = ['', 'MaterialID', 'Theta', 'Thickness',  # 'SOut'
                   ]
    _prop_build_grid_array(grid, prop, names2, names_write, irow)
    # print('done with setup_pcomp')
    return grid


def _prop_names_to_table(grid: QGridLayout,
                         prop,
                         names: list[tuple[str, str]]) -> int:
    irow = 0
    label = QLabel(str(prop.type))
    bold_font = QFont('Arial', FONT_SIZE)
    bold_font.setBold(True)
    label.setFont(bold_font)
    grid.addWidget(label, irow, 0)
    irow += 1

    for iname, namesi in enumerate(names):
        label_text, attr_name = namesi
        # print(namesi)
        label = QLabel(label_text)
        value = getattr(prop, attr_name)
        # print(label_text, attr_name, value, type(value))
        if isinstance(value, (float, int)):
            edit = QFloatEdit(_write_float(value))
        else:
            if value is None:
                value = ''
            # print(value)
            edit = QLineEdit(str(value))
        grid.addWidget(label, irow, 0)
        grid.addWidget(edit, irow, 1)
        irow += 1
        if iname == 0:
            edit.setDisabled(True)
    return irow


def setup_pbeam(grid: QGridLayout,
                prop: PBEAM) -> QGridLayout:
    names = [
        ('Area', 'Area'),
        ('I1', 'I1',),
        ('I2', 'I2',),
        ('J', 'J',),
    ]
    print(prop.get_stats())
    irow = _prop_names_to_table(grid, prop, names)
    #print('done with setup_pbeam')
    return grid


def clear_grid_layout(grid_layout):
    """Remove all widgets from the grid layout."""
    # Get a list of all items in the layout
    items = []
    for i in range(grid_layout.count()):
        items.append(grid_layout.itemAt(i))

    # Remove each item from the layout
    for item in items:
        widget = item.widget()
        if widget is not None:
            widget.setParent(None)  # Remove the widget from its parent
            widget.deleteLater()  # Schedule the widget for deletion


def main():  # pragma: no cover
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    #The Main window
    data = {
        'font_size': 9,
        'bdf_filename': '',
    }
    main_window = BeamTuner(data)
    main_window.show()
    # Enter the main loop
    app.exec_()


if __name__ == "__main__":  # pragma: no cover
    main()
