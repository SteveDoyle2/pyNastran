from __future__ import annotations
import numpy as np

from qtpy.QtWidgets import (
    QVBoxLayout, QTableWidget, QTableWidgetItem, QDialog, QHeaderView)

class ResultsDialog(QDialog):
    def __init__(self, win_parent,
                 data: np.ndarray,
                 labels: list[str],
                 title: str='Results'):
        super().__init__(win_parent)

        self.setWindowTitle(title)
        nrows, ncolumns = data.shape

        table_widget = QTableWidget(self)
        table_widget.setRowCount(nrows)
        table_widget.setColumnCount(ncolumns)
        table_widget.setHorizontalHeaderLabels(labels)
        self.table_widget = table_widget

        header = table_widget.horizontalHeader()
        for irow, row in enumerate(data):
            for jcol, value in enumerate(row):
                obj = QTableWidgetItem(str(value))
                table_widget.setItem(irow, jcol, obj)
            header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.Stretch)

        vbox = QVBoxLayout(self)
        vbox.addWidget(table_widget)
        self.setLayout(vbox)
        self.show()
