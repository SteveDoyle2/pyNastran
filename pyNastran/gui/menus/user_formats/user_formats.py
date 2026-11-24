"""
defines:
 - UserFormatsWindow

+-----------------------------------------------------------+
| User Formats                                              |
+-----------------------------------------------------------+
| File: C:/Users/me/pynastrangui_formats.json               |
|                                                           |
| Name   | Path                       | Class    | Ext      |
| cart3d | .../cart3d/cart3d_io.py    | Cart3dIO | .tri     |
|                                                           |
|  Add  Browse...  Remove                                   |
|                                                           |
|                         Apply  OK  Cancel                 |
+-----------------------------------------------------------+

"""
from __future__ import annotations
import os
from typing import Any, TYPE_CHECKING

from qtpy.QtWidgets import (
    QLabel, QPushButton, QTableWidget, QTableWidgetItem, QHeaderView,
    QHBoxLayout, QVBoxLayout, QMessageBox, QAbstractItemView)

from pyNastran.gui.utils.qt.pydialog import PyDialog
from pyNastran.gui.utils.qt.dialogs import open_file_dialog
from pyNastran.gui.user_formats import write_user_formats, UserFormats
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow

COLUMNS = ['Name', 'Path', 'Class', 'Extensions']
INAME, IPATH, ICLASS, IEXT = range(4)


class UserFormatsWindow(PyDialog):
    """
    Edits the user-defined formats file (~/pynastrangui_formats.json).

    Notes
    -----
    Apply/OK saves the json and reloads the formats in the gui, so the
    format can be used without restarting.  Comments in the json are
    not preserved by the menu.

    """
    def __init__(self, data: dict[str, Any], win_parent: MainWindow=None):
        super().__init__(data, win_parent)
        self.setWindowTitle('User Formats')
        self.filename: str = data['filename']
        user_formats: UserFormats = data['user_formats']

        self.file_label = QLabel(f'File: {self.filename}')
        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(IPATH, QHeaderView.ResizeMode.Stretch)
        self.table.setToolTip(
            'Name: format name (e.g., "pyNastranGUI -f cart3d model.tri")\n'
            'Path: a *.py file or a dotted module (pyNastran.converters.cart3d.cart3d_io)\n'
            'Class: the format class (e.g., Cart3dIO)\n'
            'Extensions: optional; comma separated (e.g., .tri, .triq)\n'
            'A user format overwrites a built-in format with the same name.')
        for name, value in user_formats.items():
            self._add_row(name, value['path'], value['class'],
                          ', '.join(value.get('extensions', [])))

        self.add_button = QPushButton('Add')
        self.browse_button = QPushButton('Browse...')
        self.remove_button = QPushButton('Remove')
        self.apply_button = QPushButton('Apply')
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Cancel')
        self.resize(800, 300)

        self.setup_layout()
        self.set_connections()
        self.set_font_size(data['font_size'])

    def setup_layout(self) -> None:
        edit_box = QHBoxLayout()
        edit_box.addWidget(self.add_button)
        edit_box.addWidget(self.browse_button)
        edit_box.addWidget(self.remove_button)
        edit_box.addStretch()

        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addStretch()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.ok_button)
        ok_cancel_box.addWidget(self.cancel_button)

        vbox = QVBoxLayout()
        vbox.addWidget(self.file_label)
        vbox.addWidget(self.table)
        vbox.addLayout(edit_box)
        vbox.addLayout(ok_cancel_box)
        self.setLayout(vbox)

    def set_connections(self) -> None:
        self.add_button.clicked.connect(self.on_add)
        self.browse_button.clicked.connect(lambda: self.on_browse())
        self.remove_button.clicked.connect(self.on_remove)
        self.apply_button.clicked.connect(self.on_apply)
        self.ok_button.clicked.connect(self.on_ok)
        self.cancel_button.clicked.connect(self.on_cancel)

    def _add_row(self, name: str='', path: str='', class_name: str='',
                 extensions: str='') -> int:
        irow = self.table.rowCount()
        self.table.insertRow(irow)
        for icol, value in enumerate([name, path, class_name, extensions]):
            self.table.setItem(irow, icol, QTableWidgetItem(value))
        return irow

    def on_add(self) -> None:
        """adds an empty row"""
        irow = self._add_row()
        self.table.setCurrentCell(irow, INAME)

    def on_browse(self, path: str='') -> None:
        """
        Picks a *.py file for the selected row (or a new row).
        The name defaults to the file name without ``_io``
        (cart3d_io.py -> cart3d).
        """
        if not path:
            irow = self.table.currentRow()
            default_dirname = ''
            if irow >= 0:
                default_dirname = os.path.dirname(self._text(irow, IPATH))
            path, unused_flt = open_file_dialog(
                self, 'Select a format module', default_dirname,
                'Python (*.py);;All files (*)')
            if not path:
                return
        irow = self.table.currentRow()
        if irow < 0:
            irow = self._add_row()
        self.table.item(irow, IPATH).setText(path)
        if not self._text(irow, INAME):
            name = os.path.splitext(os.path.basename(path))[0]
            if name.endswith('_io'):
                name = name[:-3]
            self.table.item(irow, INAME).setText(name)

    def on_remove(self) -> None:
        """removes the selected rows"""
        irows = sorted({index.row() for index in self.table.selectedIndexes()},
                       reverse=True)
        if not irows and self.table.currentRow() >= 0:
            irows = [self.table.currentRow()]
        for irow in irows:
            self.table.removeRow(irow)

    def _text(self, irow: int, icol: int) -> str:
        item = self.table.item(irow, icol)
        return '' if item is None else item.text().strip()

    def get_user_formats(self) -> tuple[UserFormats, list[str]]:
        """
        Gets the formats from the table

        Returns
        -------
        user_formats : dict[str, dict[str, Any]]
            name -> {'path': str, 'class': str, 'extensions': list[str]}
        errors : list[str]
            empty rows are skipped; partial rows are errors

        """
        user_formats: UserFormats = {}
        errors: list[str] = []
        for irow in range(self.table.rowCount()):
            name = self._text(irow, INAME).lower()
            path = self._text(irow, IPATH)
            class_name = self._text(irow, ICLASS)
            extensions = [ext for ext in self._text(irow, IEXT).replace(';', ',').split(',')
                          if ext.strip()]
            if not (name or path or class_name):
                continue
            if not name or not path or not class_name:
                errors.append(f'row {irow+1}: Name, Path, and Class are required')
                continue
            if name in user_formats:
                errors.append(f'row {irow+1}: duplicate name {name!r}')
                continue
            user_formats[name] = {'path': path, 'class': class_name,
                                  'extensions': extensions}
        return user_formats, errors

    def on_apply(self) -> bool:
        """saves the json and reloads the formats"""
        user_formats, errors = self.get_user_formats()
        if errors:
            msg = '\n'.join(errors)
            if self.win_parent is not None and not getattr(self.win_parent, 'is_gui', True):
                self.win_parent.log_error(msg)
            else:  # pragma: no cover
                QMessageBox.warning(self, 'User Formats', msg)
            return False
        write_user_formats(user_formats, self.filename)
        if self.win_parent is not None:
            self.win_parent.load_user_formats(self.filename)
        return True

    def on_ok(self) -> None:
        if self.on_apply():
            self.close()

    def on_cancel(self) -> None:
        self.close()
