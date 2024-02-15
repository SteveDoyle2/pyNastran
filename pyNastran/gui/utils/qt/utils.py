"""small Qt utils"""
from typing import Any
from qtpy.QtWidgets import QBoxLayout, QVBoxLayout, QHBoxLayout, QGridLayout, QWidget

def add_obj_to_vbox(vbox: QVBoxLayout, widget_layout) -> None:
    if isinstance(widget_layout, QBoxLayout):
        vbox.addLayout(widget_layout)
    else:
        vbox.addWidget(widget_layout)

def create_hbox_with_widgets(widgets: list[Any]) -> QHBoxLayout:
    hbox = QHBoxLayout()
    for widget in widgets:
        hbox.addWidget(widget)
    return hbox

def add_line_widgets_to_grid(grid: QGridLayout, irow: int, widgets: list[QWidget]) -> int:
    for j, widget in enumerate(widgets):
        grid.addWidget(widget, irow, j)
    return irow + 1
