import numpy as np
import pandas as pd
from qtpy.QtWidgets import (
    QInputDialog,
    QApplication, QTableWidget, QTableWidgetItem, QMenu,
    QAbstractItemView,)
from qtpy import QtCore
Qt = QtCore.Qt

class QTableWidgetCopy(QTableWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.setContextMenuPolicy(Qt.CustomContextMenu)
        # self.customContextMenuRequested.connect(self.show_context_menu)

        rename_column_support = True
        add_remove_row_support = True

        self.setVerticalScrollMode(QAbstractItemView.ScrollPerPixel)
        if rename_column_support:
            # Enable custom context menu for horizontal header
            self.horizontalHeader().setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
            self.horizontalHeader().customContextMenuRequested.connect(self.show_header_context_menu)
        if add_remove_row_support:
            # Add/remove rows
            self.verticalHeader().setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
            self.verticalHeader().customContextMenuRequested.connect(self.show_row_header_context_menu)

    def load_table_data(self,
                        headers: list[str],
                        data: np.ndarray | list[list[int | float | str]]) -> None:
        """
        headers = ['A', 'Col2', 'Col3', 'Col4']
        data = np.arange(20).reshape(4,5)
        self.load_table_data(headers, data)
        """
        # Set row and column counts
        # print(f'data = {data}')
        self.clear()
        self.setRowCount(len(data))
        self.setColumnCount(len(data[0]) if len(data) else 0)

        # Set headers (optional)
        # headers = ["Column 1", "Column 2", "Column 3"]  # Replace with your headers
        self.setHorizontalHeaderLabels(headers)

        # Populate the table
        for row_idx, row_data in enumerate(data):
            for col_idx, col_data in enumerate(row_data):
                item = QTableWidgetItem(str(col_data))
                self.setItem(row_idx, col_idx, item)

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_C and (event.modifiers() & Qt.KeyboardModifier.ControlModifier):
            self.on_copy_table(event)
        # elif event.key() == Qt.Key.Key_X and (event.modifiers() & Qt.KeyboardModifier.ControlModifier):
        #     self.on_cut_table(event)
        # elif event.key() == Qt.Key.Key_P and (event.modifiers() & Qt.KeyboardModifier.ControlModifier):
        #     self.on_paste_table(event)
        else:
            super().keyPressEvent(event)

    def on_copy_table(self, event):
        copied_cells = self.selectedIndexes()
        if not copied_cells:
            return
        # Sort to maintain order, handle rows/cols with \t and \n
        # copied_cells = sorted(copied_cells)

        copy_text = ''
        max_column = copied_cells[-1].column()
        for c in copied_cells:
            copy_text += self.item(c.row(), c.column()).text()
            if c.column() == max_column:
                copy_text += '\n'
            else:
                copy_text += '\t'
        print(f'text = {copy_text}')
        QApplication.clipboard().setText(copy_text)

    def show_context_menu(self, pos):
        if self.rowCount() == 0 or self.columnCount() == 0:
            return
        print(f'pos = {pos}')
        # Get the item at the clicked position
        item = self.itemAt(pos)
        print(f'item = {item}')
        if not item:
            return
        # Get global position to display menu correctly
        global_pos = self.mapToGlobal(pos)
        print(f'global_pos = {global_pos}')

        # Create the menu
        context_menu = QMenu(self)
        copy_action = context_menu.addAction('Copy')
        delete_action = context_menu.addAction('Delete Row')

        # Execute the menu and wait for an action selection
        action = context_menu.exec_(global_pos)

        if action == copy_action:
            # Handle copy logic
            print(f'Copying cell content: {item.text()}')
        elif action == delete_action:
            # Handle delete logic
            row = item.row()
            self.removeRow(row)
            # print(f'Deleting row {row}')

    def show_header_context_menu(self, position):
        """Show context menu when right-clicking on column header"""
        # Get the column index from the position
        column = self.horizontalHeader().logicalIndexAt(position)
        if column < 0:
            return

        # Create context menu
        menu = QMenu(self)
        rename_action = menu.addAction("Rename Column")

        # copyAction = menu.addAction('Copy')
        delete_action = menu.addAction('Delete Column')
        insert_left_action = menu.addAction('Insert Column Left')
        insert_right_action = menu.addAction('Insert Column Right')

        # Show menu and get selected action
        action = menu.exec(self.horizontalHeader().mapToGlobal(position))

        if action == rename_action:
            self.rename_column(column)
        elif action == insert_left_action:
            col = self.horizontalHeader().logicalIndexAt(position)
            self.insert_col(col)
        elif action == insert_right_action:
            col = self.horizontalHeader().logicalIndexAt(position)
            self.insert_col(col+1)
        elif action == delete_action:
            col = self.horizontalHeader().logicalIndexAt(position)
            self.delete_col(col)

    def rename_column(self, column: int):
        """Open dialog to rename column header"""
        # Get current column name
        current_name = self.horizontalHeaderItem(column).text() if self.horizontalHeaderItem(
            column) else f"Column {column}"

        # Show input dialog
        new_name, ok = QInputDialog.getText(
            self,
            "Rename Column",
            "Enter new column name:",
            text=current_name
        )

        if ok and new_name:
            self.setHorizontalHeaderItem(column, QTableWidgetItem(new_name))
    #-----------------
    # def show_col_header_context_menu(self, position):
    #     """Show context menu when right-clicking on row header (row numbers)"""
    #     # Get the col index from the position
    #     col = self.horizontalHeader().logicalIndexAt(position)
    #
    #     if col < 0:
    #         return
    #
    #     # Create context menu
    #     menu = QMenu(self)
    #     insert_left_action = menu.addAction("Insert Column Left")
    #     insert_right_action = menu.addAction("Insert Column Right")
    #     menu.addSeparator()
    #     delete_action = menu.addAction("Delete Col")
    #
    #     # Show menu and get selected action
    #     action = menu.exec(self.horizontalHeader().mapToGlobal(position))
    #
    #     if action == insert_left_action:
    #         self.insert_col(col)
    #     elif action == insert_right_action:
    #         self.insert_col(col + 1)
    #     elif action == delete_action:
    #         self.delete_col(col)

    def show_row_header_context_menu(self, position):
        """Show context menu when right-clicking on row header (row numbers)"""
        # Get the row index from the position
        row = self.verticalHeader().logicalIndexAt(position)

        if row < 0:
            return

        # Create context menu
        menu = QMenu(self)
        insert_above_action = menu.addAction("Insert Row Above")
        insert_below_action = menu.addAction("Insert Row Below")
        menu.addSeparator()
        delete_action = menu.addAction("Delete Row")

        # Show menu and get selected action
        action = menu.exec(self.verticalHeader().mapToGlobal(position))

        if action == insert_above_action:
            self.insert_row(row)
        elif action == insert_below_action:
            self.insert_row(row + 1)
        elif action == delete_action:
            self.delete_row(row)

    def insert_row(self, row: int):
        """Insert a new empty row at the specified position"""
        self.insertRow(row)

        # Initialize empty cells in the new row
        for col in range(self.columnCount()):
            self.setItem(row, col, QTableWidgetItem(""))
        # print(f'Inserted row at position {row}')

    def insert_col(self, col: int):
        """Insert a new empty col at the specified position"""
        self.insertColumn(col)

        # Initialize empty cells in the new row
        for col in range(self.columnCount()):
            self.setItem(col, col, QTableWidgetItem(""))
        # print(f'Inserted col at position {col}')

    def delete_row(self, row: int):
        """Delete the specified row"""
        if self.rowCount() > 0:
            self.removeRow(row)
            # print(f'Deleted row {row}')

    def delete_col(self, col: int):
        """Delete the specified col"""
        if self.columnCount() > 0:
            self.removeColumn(col)
            # print(f'Deleted col {col}')

    def delete_selected_rows(self):
        """Delete all selected rows (triggered by Delete key)"""
        selected_rows = set()
        for item in self.selectedItems():
            selected_rows.add(item.row())

        # Delete from highest to lowest to avoid index shifting issues
        for row in sorted(selected_rows, reverse=True):
            self.removeRow(row)

        if selected_rows:
            print(f'Deleted {len(selected_rows)} row(s)')

    def delete_selected_cols(self):
        """Delete all selected rows (triggered by Delete key)"""
        selected_cols = set()
        for item in self.selectedItems():
            selected_cols.add(item.row())

        # Delete from highest to lowest to avoid index shifting issues
        for col in sorted(selected_cols, reverse=True):
            self.removeColumn(col)

        if selected_cols:
            print(f'Deleted {len(selected_cols)} col(s)')

    def get_data(self,
                 skip_empty_rows=True,
                 strip_whitespace=True,
                 convert_numeric=False):
            """
            Extract data from QTableWidget to pandas DataFrame with advanced options

            Parameters
            ----------
            skip_empty_rows : bool
                If True, skip rows that are entirely empty (default: True)
            strip_whitespace : bool
                If True, strip leading/trailing whitespace from all cells (default: True)
            convert_numeric : bool
                If True, attempt to convert numeric strings to numbers (default: False)

            Returns
            -------
            pandas.DataFrame
                DataFrame containing the table data
            """
            rows = self.rowCount()
            cols = self.columnCount()

            # Get headers
            headers = []
            for col in range(cols):
                header = self.horizontalHeaderItem(col)
                if header is not None:
                    headers.append(header.text())
                else:
                    headers.append(f"Column_{col}")
            # print(f'headers = {headers}')

            # Get data
            data = []
            for row in range(rows):
                row_data = []
                for col in range(cols):
                    item = self.item(row, col)
                    if item is not None:
                        cell_value = item.text()
                        if strip_whitespace:
                            cell_value = cell_value.strip()
                        row_data.append(cell_value)
                    else:
                        row_data.append("")
                # print(f'row_data = {row_data}')

                # Check if row is entirely empty
                if skip_empty_rows:
                    if any(cell for cell in row_data):  # If any cell has content
                        data.append(row_data)
                else:
                    data.append(row_data)

            # Convert numeric columns if requested
            df = pd.DataFrame(data, columns=headers)
            if convert_numeric:
                df = df.apply(pd.to_numeric, errors='ignore')
            return df
