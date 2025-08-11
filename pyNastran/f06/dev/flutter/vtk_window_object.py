"""
defines:
 - VtkWindowObject

"""
from __future__ import annotations
from typing import Any, Optional, TYPE_CHECKING
from qtpy.QtWidgets import QMainWindow
from pyNastran.utils import PathLike
from pyNastran.f06.dev.flutter.vtk_window import VtkWindow

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.f06.dev.flutter.gui_flutter import FlutterGui

DT_MS_DEFAULT = 100
NPHASE_DEFAULT = 10
DT_MS_MIN = 100
DT_MS_MAX = 5000


class VtkWindowObject:
    """defines VtkWindowObject, which is an interface to the PreferencesWindow"""
    def __init__(self, gui: FlutterGui, icon_path: PathLike):
        self.gui = gui
        self.ncase = 0
        self.icon_path = icon_path
        self.window_shown = None
        self.window = None

        self.dt_ms = DT_MS_DEFAULT
        self.nphase = NPHASE_DEFAULT
        self.icase = 0
        self.animate = True

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
                        icase: Optional[int]=None,
                        animate: Optional[bool]=None) -> None:
        """
        To avoid lag:
         - Stop the timer if something got updated
         - Update the settings
         - Restart the timer at iphase=0

        """
        is_updated = False
        apply_fringe_plot = False
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
                apply_fringe_plot = True
                self.window.timer.stop()
                self.window.icase = self.icase

        if animate is None:
            animate = self.animate
        elif animate != self.animate:
            self.animate = animate
            if self.window_shown:
                is_updated = True
                self.window.timer.stop()
                self.window.animate = self.animate

        if apply_fringe_plot:
            # avoid redoing the fringe
            self.window.apply_fringe_plot = apply_fringe_plot

        if is_updated and animate:
            # only restart the animation if the timer was stopped
            # and we're animating
            self.window.iphase = 0
            self.window.update_dphases()
            self.window.timer.start()

    @property
    def data(self) -> dict[str, int]:
        out = {
            'icase': self.icase,
            'dt_ms': self.dt_ms,
            'nphase': self.nphase,
            'animate': self.animate,
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
        #print(f'data = {data}')
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

    def reset_icase_ncase(self, icase: int, ncase: int) -> None:
        """called by VtkWindow when icase exceeds the min/max"""
        self.icase = icase
        self.gui._export_settings_obj.reset_icase_ncase(icase, ncase)

    def on_close(self):
        #del self.window
        print('on close...')
        self.window_shown = False
        #del self.window
        #self.window.hide()
