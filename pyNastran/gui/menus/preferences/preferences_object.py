from __future__ import annotations
from typing import Any, TYPE_CHECKING
from pyNastran.gui.menus.preferences.preferences import PreferencesWindow
from pyNastran.gui.gui_objects.settings import (
    Settings, NASTRAN_BOOL_KEYS,
    NastranSettings, OtherSettings)
from pyNastran.gui.qt_files.base_gui import BaseGui
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow


class PreferencesObject(BaseGui):
    def __init__(self, gui: MainWindow):
        #self.gui = gui
        super().__init__(gui)
        self._preferences_window_shown = False
        self._preferences_window = None

    def set_font_size(self, font_size: int) -> None:
        """sets the font size for the preferences window"""
        if self._preferences_window_shown:
            self._preferences_window.set_font_size(font_size)

    def get_preferences_data(self) -> dict[str, Any]:
        """gets the input data for the PreferencesWindow"""
        camera = self.gui.GetCamera()
        min_clip, max_clip = camera.GetClippingRange()
        settings: Settings = self.gui.settings
        nastran_settings: NastranSettings = settings.nastran_settings
        data = {
            'picker_size' : self.gui.element_picker_size,
            'show_corner_coord' : self.gui.get_corner_axis_visiblity(),

            'min_clip' : min_clip,
            'max_clip' : max_clip,

            # not saved to the json, so they're not in add_settings_to_dict
            'nastran_real_modes_to_include': nastran_settings.real_modes_to_include,
            'nastran_complex_modes_to_include': nastran_settings.complex_modes_to_include,

            #----------------------------------------
            'clicked_ok' : False,
            'close' : False,
        }
        settings.add_settings_to_dict(data)
        if hasattr(self.gui, 'get_format_preferences'):
            data['formats'] = self.gui.get_format_preferences()
        return data

    def set_preferences_menu(self) -> None:
        """
        Opens a dialog box to set:

        +--------+----------+
        |  Max   |  Float   |
        +--------+----------+
        """
        #if not hasattr(self, 'case_keys'):  # TODO: maybe include...
            #self.log_error('No model has been loaded.')
            #return
        settings: Settings = self.gui.settings
        data = self.get_preferences_data()

        if not self._preferences_window_shown:
            self._preferences_window = PreferencesWindow(data, win_parent=self.gui)
            self._preferences_window.show()
            self._preferences_window_shown = True
            self._preferences_window.exec_()
        else:
            self._preferences_window.activateWindow()

        if data['close']:
            if not self._preferences_window._updated_preference:
                settings.on_set_font_size(data['font_size'])
            del self._preferences_window
            self._preferences_window_shown = False
        else:
            self._preferences_window.activateWindow()
