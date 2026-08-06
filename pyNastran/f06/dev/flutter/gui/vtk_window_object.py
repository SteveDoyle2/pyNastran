"""
defines:
 - VtkWindowObject

"""
from __future__ import annotations
from typing import Any, Optional, TYPE_CHECKING
from pyNastran.utils import PathLike
from pyNastran.f06.dev.flutter.vtk_data import (
    DT_MS_MIN, DT_MS_MAX,
    apply_vtk_settings, VtkData)
from pyNastran.f06.dev.flutter.vtk_window import VtkWindow

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.f06.dev.flutter.gui_flutter import FlutterGui


class VtkWindowObject:
    """defines VtkWindowObject, which is an interface to the PreferencesWindow"""
    def __init__(self, gui: FlutterGui, vtk_data: VtkData,
                 icon_path: PathLike):
        self.gui = gui
        self.ncase = 0
        self.icon_path = icon_path
        self.vtk_data = vtk_data
        self.window_shown = None
        self.window = None

        # self.dt_ms = DT_MS_DEFAULT
        # self.nphase = NPHASE_DEFAULT
        # self.icase = 0
        # self.animate = True

        self.apply_settings({})
        # self.dt_ms = DT_MS_DEFAULT
        # self.nphase = NPHASE_DEFAULT

    def apply_settings(self, data: dict[str, Any]) -> None:
        apply_vtk_settings(self, data)

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
        if dt_ms is not None and dt_ms != self.vtk_data.dt_ms:
            dt_ms = max(DT_MS_MIN, min(dt_ms, DT_MS_MAX))
            self.vtk_data.dt_ms = dt_ms
            if self.window_shown:
                is_updated = True
                self.window.timer.stop()  # Pause the timer
                self.window.timer.setInterval(dt_ms)

        if nphase is not None and nphase != self.vtk_data.nphase:
            self.vtk_data.nphase = nphase
            if self.window_shown:
                is_updated = True
                self.window.timer.stop()

        if icase is not None and icase != self.vtk_data.icase:
            self.vtk_data.icase = icase
            if self.window_shown:
                is_updated = True
                apply_fringe_plot = True
                self.window.timer.stop()

        if animate is None:
            animate = self.vtk_data.animate
        elif animate != self.vtk_data.animate:
            self.vtk_data.animate = animate
            if self.window_shown:
                is_updated = True
                self.window.timer.stop()

        if apply_fringe_plot:
            # avoid redoing the fringe
            self.window.apply_fringe_plot = apply_fringe_plot

        if is_updated and animate:
            # only restart the animation if the timer was stopped
            # and we're animating
            self.vtk_data.iphase = 0
            self.window.update_dphases()
            self.window.timer.start()

    @property
    def data(self) -> dict[str, int]:
        return self.vtk_data.to_json()
        # out = {
        #     'icase': self.icase,
        #     'dt_ms': self.dt_ms,
        #     'nphase': self.nphase,
        #     'animate': self.animate,
        # }
        # return out

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
            'dt_ms': self.vtk_data.dt_ms,
            'nphase': self.vtk_data.nphase,
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
        self.vtk_data.icase = icase
        self.gui._export_settings_obj.reset_icase_ncase(icase, ncase)

    def on_close(self):
        #del self.window
        print('on close...')
        self.window_shown = False
        #del self.window
        #self.window.hide()
