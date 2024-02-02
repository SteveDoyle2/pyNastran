from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from pyNastran.gui.main_window import MainWindow


class BaseGui:
    def __init__(self, gui: MainWindow):
        self.gui = gui

    @property
    def log(self):
        """links the the GUI's log"""
        return self.gui.log

    @property
    def settings(self):
        """gets the gui settings"""
        return self.gui.settings

