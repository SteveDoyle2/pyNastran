"""
defines:
 - PreferencesObject

"""
from __future__ import annotations
from typing import TYPE_CHECKING
from qtpy.QtWidgets import QMainWindow
from pyNastran.f06.dev.flutter.preferences import PreferencesDialog
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.f06.dev.flutter.gui_flutter import FlutterGui


class PreferencesObject:
    """defines PreferencesObject, which is an interface to the PreferencesWindow"""
    def __init__(self, gui: FlutterGui):
        self.gui = gui
        self.window_shown = None
        self.window = None

    # def show_window(self) -> None:
    #     """shows the window"""
    #     if self.window_shown:
    #         self.window.show()

    # def hide_window(self) -> None:
    #     """hides the window"""
    #     if self.window_shown:
    #         self.window.hide()

    def set_font_size(self, font_size: int) -> None:
        """sets the font size for the legend window"""
        if self.window_shown:
            self.window.set_font_size(font_size)
        # if self._animation_window_shown:
        #     self._animation_window.set_font_size(font_size)

    def show(self):
        """Opens a dialog box to set"""
        gui: FlutterGui = self.gui
        # if not hasattr(gui, 'case_keys') or len(gui.case_keys) == 0:
        #     gui.log_error('No model has been loaded.')
        #     return
        vtk_obj = gui._vtk_window_obj
        data = {
            # vtk
            'nphase': vtk_obj.nphase,
            'icase': vtk_obj.icase,
            'ncase': vtk_obj.ncase,
            'animate': vtk_obj.animate,
            'dt_ms': vtk_obj.dt_ms,

            # plotting
            'font_size': gui.font_size,
            'plot_font_size': gui.plot_font_size,
            'export_to_png': gui.export_to_png,
            'export_to_csv': gui.export_to_csv,
            'export_to_f06': gui.export_to_f06,
            'export_to_zona': gui.export_to_zona,
            'clicked_ok' : False,
            'close' : False,
        }
        if self.window_shown in {True, False}:
            self.window_shown = True
            self.window.activateWindow()
            self.window.show()
        else:
            self.window_shown = True
            self.window = PreferencesDialog(data, self, win_parent=gui)

        # if data['close']:
        #     # if not self.window._updated_legend:
        #     #     self.apply_legend(data)
        #     self.window_shown = False
        #     self.window = None
        # else:
        #     self.window.activateWindow()

    def reset_icase_ncase(self, icase: int, ncase: int) -> None:
        if self.window_shown:
            self.window.icase_edit.setValue(icase)
            self.window.icase_edit.setMaximum(ncase)

    def on_dt_ms(self, dt_ms: int) -> None:
        self.gui._vtk_window_obj.set_preferences(dt_ms=dt_ms)

    def on_nphase(self, nphase: int) -> None:
        self.gui._vtk_window_obj.set_preferences(nphase=nphase)
    def on_icase(self, icase: int) -> None:
        self.gui._vtk_window_obj.set_preferences(icase=icase)
    def on_animate(self, animate: bool) -> None:
        self.gui._vtk_window_obj.set_preferences(animate=animate)

    def on_close(self):
        #del self.window
        self.window_shown = False
        self.window.hide()
