"""
defines:
 - VtkWindowObject

"""
from __future__ import annotations
from typing import Any, Optional, TYPE_CHECKING
from qtpy.QtWidgets import QMainWindow
from pyNastran.f06.dev.flutter.vtk_window import VtkWindow

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.f06.dev.flutter.gui_flutter import FlutterGui

DT_MS_DEFAULT = 100
NPHASE_DEFAULT = 10
DT_MS_MIN = 100
DT_MS_MAX = 5000


class VtkWindowObject:
    """defines VtkWindowObject, which is an interface to the PreferencesWindow"""
    def __init__(self, gui: FlutterGui):
        self.gui = gui
        self.window_shown = None
        self.window = None
        self.apply_settings({})
        # self.dt_ms = DT_MS_DEFAULT
        # self.nphase = NPHASE_DEFAULT
    def apply_settings(self, data: dict[str, Any]) -> None:
        vtk_data = data.get('vtk', {})
        self.dt_ms = DT_MS_DEFAULT if 'dt_ms' not in vtk_data else int(vtk_data['dt_ms'])
        self.nphase = NPHASE_DEFAULT if 'nphase' not in vtk_data else int(vtk_data['nphase'])
        self.icase = 0 if 'icase' not in vtk_data else int(vtk_data['icase'])
        self.animate = True

    def set_preferences(self, dt_ms: Optional[int]=None,
                        nphase: Optional[int]=None,
                        icase: Optional[int]=None) -> None:
        is_updated = False
        is_timer_paused = False
        if dt_ms is not None and dt_ms != self.dt_ms:
            dt_ms = max(DT_MS_MIN, min(dt_ms, DT_MS_MAX))
            self.dt_ms = dt_ms
            if self.window_shown:
                is_updated = True
                self.window.timer.stop()  # Pause the timer
                self.window.dt_ms = self.dt_ms
                self.window.timer.setInterval(dt_ms)
        if nphase is not None and nphase != self.nphase:
            self.nphase = nphase
            if self.window_shown:
                is_updated = True
                self.window.timer.stop()
                self.window.nphase = self.nphase
        if icase is not None and icase != self.icase:
            self.icase = icase
            if self.window_shown:
                is_updated = True
                self.window.timer.stop()
                self.window.icase = self.icase

        if is_updated:
            self.window.iphase = 0
            self.window.timer.start()

    # def set_dt_ms(self, dt_ms: int) -> None:
    #     self.dt_ms = dt_ms
    #     if self.window_shown:
    #         self.window.dt_ms = self.dt_ms
    #         self.window.timer.setInterval(dt_ms)
    #         self.window.iphase = 0

    # def set_nphase(self, nphase: int) -> None:
    #     self.nphase = nphase
    #     if self.window_shown:
    #         self.window.nphase = self.nphase
    #         self.window.iphase = 0

    # def set_icase(self, icase: int) -> None:
    #     self.icase = icase
    #     if self.window_shown:
    #         self.window.icase = self.icase
    #         self.window.iphase = 0

    @property
    def data(self) -> dict[str, int]:
        out = {
            'icase': self.icase,
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
