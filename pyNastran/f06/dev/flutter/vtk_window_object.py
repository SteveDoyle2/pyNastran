"""
defines:
 - VtkWindowObject

"""
from __future__ import annotations
from typing import TYPE_CHECKING
from qtpy.QtWidgets import QMainWindow
from pyNastran.f06.dev.flutter.vtk_window import VtkWindow

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.f06.dev.flutter.gui_flutter import FlutterGui


class VtkWindowObject:
    """defines VtkWindowObject, which is an interface to the PreferencesWindow"""
    def __init__(self, gui: FlutterGui):
        self.gui = gui
        self.window_shown = None
        self.window = None
        self.dt_ms = 150
        self.nphase = 10

    @property
    def data(self) -> dict[str, int]:
        out = {
            'dt_ms': self.dt_ms,
            'nphase': self.nphase,
            'animate': True,
        }
        return out

    def show_window(self) -> None:
        """shows the window"""
        if self.window_shown:
            self.window.show_legend()

    def hide_window(self) -> None:
        """hides the widnow"""
        if self.window_shown:
            self.window.hide_legend()

    def set_font_size(self, font_size: int) -> None:
        """sets the font size for the window"""
        if self.window_shown:
            self.window.set_font_size(font_size)
        # if self._animation_window_shown:
        #     self._animation_window.set_font_size(font_size)

    def show(self, bdf_filename: str, op2_filename: str):
        """Opens a dialog box to set"""
        gui: FlutterGui = self.gui
        # if not hasattr(gui, 'case_keys') or len(gui.case_keys) == 0:
        #     gui.log_error('No model has been loaded.')
        #     return

        data = {
            'dt_ms': self.dt_ms,
            'nphase': self.nphase,
        }
        print(f'data = {data}')
        if self.window_shown in {True, False}:
            self.window_shown = True
            self.window.set_data(data)
            print('activating...')
            self.window.activateWindow()
            self.window.show()
            print('showed')
        else:
            self.window_shown = True
            self.window = VtkWindow(
                self, gui, data, bdf_filename, op2_filename)

    # def set_data(self):
    #     asdf

    def on_close(self):
        #del self.window
        print('on close...')
        self.window_shown = False
        #del self.window
        #self.window.hide()
