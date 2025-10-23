"""
defines:
 - PreferencesObject

"""
from __future__ import annotations
from typing import Optional, TYPE_CHECKING
# from qtpy.QtWidgets import QMainWindow
from pyNastran.f06.dev.flutter.preferences import FlutterPreferencesDialog
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.f06.dev.flutter.gui_flutter import FlutterGui


class FlutterPreferencesObject:
    """defines PreferencesObject, which is an interface to the PreferencesWindow"""
    def __init__(self, gui: FlutterGui, use_vtk: bool):
        self.gui = gui
        self.window_shown = None
        self.window = None
        self.use_vtk = use_vtk

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
        self.gui.on_run()

    def show(self):
        """Opens a dialog box to set"""
        gui: FlutterGui = self.gui
        # if not hasattr(gui, 'case_keys') or len(gui.case_keys) == 0:
        #     gui.log_error('No model has been loaded.')
        #     return
        if hasattr(gui, '_vtk_window_obj'):
            vtk_obj = gui._vtk_window_obj
            data = {
                'nphase': vtk_obj.nphase,
                'icase': vtk_obj.icase,
                'ncase': vtk_obj.ncase,
                'animate': vtk_obj.animate,
                'dt_ms': vtk_obj.dt_ms,
            }
        else:
            data = {
                'nphase': 30,
                'icase': 0,
                'ncase': 1,
                'animate': True,
                'dt_ms': 100,
            }
        data.update({
            # plotting
            'font_size': gui.font_size,
            'plot_font_size': gui.plot_font_size,
            'export_to_png': gui.export_to_png,
            'export_to_csv': gui.export_to_csv,
            'export_to_f06': gui.export_to_f06,
            'export_to_zona': gui.export_to_zona,

            'divergence_legend_loc': gui.divergence_legend_loc,
            'flutter_bbox_to_anchor_x': gui.flutter_bbox_to_anchor_x,
            'flutter_ncolumns': gui.flutter_ncolumns,
            'freq_ndigits': gui.freq_ndigits,
            'freq_divergence_tol': gui.freq_divergence_tol,
            'auto_update': gui.auto_update,
            'clicked_ok': False,
            'close': False,
        })
        if self.window_shown in {True, False}:
            self.window_shown = True
            self.window.activateWindow()
            self.window.show()
        else:
            self.window_shown = True
            self.window = FlutterPreferencesDialog(data, self, win_parent=gui)

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

    def on_flutter_ncolumns(self, value: Optional[int]) -> None:
        self.gui.flutter_ncolumns = value
        self._on_run()

    def on_freq_ndigits(self, value: int) -> None:
        self.gui.freq_ndigits = value
        self._on_run()

    def on_freq_divergence_tol(self, value: float) -> None:
        self.gui.freq_divergence_tol = value
        self._on_run()

    def on_auto_update(self, value: bool) -> None:
        self.gui.auto_update = value
        self._on_run()

    def on_flutter_bbox_to_anchor_x(self, value: float) -> None:
        self.gui.flutter_bbox_to_anchor_x = value
        self._on_run()

    def on_divergence_legend_loc(self, value: str) -> None:
        self.gui.divergence_legend_loc = value
        self._on_run()

    def _on_run(self) -> None:
        if self.gui.auto_update:
            self.gui.on_run()

    def on_icase(self, icase: int) -> None:
        self.gui._vtk_window_obj.set_preferences(icase=icase)

    def on_animate(self, animate: bool) -> None:
        self.gui._vtk_window_obj.set_preferences(animate=animate)

    def on_close(self) -> None:
        # del self.window
        self.window_shown = False
        self.window.hide()
