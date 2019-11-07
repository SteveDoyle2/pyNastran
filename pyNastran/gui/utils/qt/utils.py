"""small Qt utils"""
from typing import List, Any
from qtpy.QtWidgets import QBoxLayout, QHBoxLayout, QGridLayout

def add_obj_to_vbox(vbox, widget_layout):
    if isinstance(widget_layout, QBoxLayout):
        vbox.addLayout(widget_layout)
    else:
        vbox.addWidget(widget_layout)

def create_hbox_with_widgets(widgets: List[Any]):
    hbox = QHBoxLayout()
    for widget in widgets:
        hbox.addWidget(widget)
    return hbox

def add_line_widgets_to_grid(grid: QGridLayout, irow: int, widgets: List[Any]):
    for j, widget in enumerate(widgets):
        grid.addWidget(widget, irow, j)
    return irow + 1
