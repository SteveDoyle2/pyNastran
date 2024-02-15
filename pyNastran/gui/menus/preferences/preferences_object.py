from pyNastran.gui.menus.preferences.preferences import PreferencesWindow
from pyNastran.gui.gui_objects.settings import Settings, NASTRAN_BOOL_KEYS
from pyNastran.gui.qt_files.base_gui import BaseGui


class PreferencesObject(BaseGui):
    def __init__(self, gui):
        #self.gui = gui
        super().__init__(gui)
        self._preferences_window_shown = False
        self._preferences_window = None

    def set_font_size(self, font_size: int) -> None:
        """sets the font size for the preferences window"""
        if self._preferences_window_shown:
            self._preferences_window.set_font_size(font_size)

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

        camera = self.gui.GetCamera()
        min_clip, max_clip = camera.GetClippingRange()
        settings: Settings = self.gui.settings
        data = {
            'font_size' : settings.font_size,
            'annotation_size' : settings.annotation_size, # int
            'annotation_color' : settings.annotation_color,

            'use_startup_directory': settings.use_startup_directory,

            'use_gradient_background' : settings.use_gradient_background,
            'use_parallel_projection': settings.use_parallel_projection,
            'background_color' : settings.background_color,
            'background_color2' : settings.background_color2,

            'corner_text_size' : settings.corner_text_size,
            'corner_text_color' : settings.corner_text_color,

            'highlight_color' : settings.highlight_color,
            'highlight_opacity' : settings.highlight_opacity,
            'highlight_point_size' : settings.highlight_point_size,

            'picker_size' : self.gui.element_picker_size,
            'dim_max' : settings.dim_max,
            'coord_scale' : settings.coord_scale,
            'coord_text_scale' : settings.coord_text_scale,
            'show_corner_coord' : self.gui.get_corner_axis_visiblity(),
            'magnify' : settings.magnify,

            'min_clip' : min_clip,
            'max_clip' : max_clip,

            'clicked_ok' : False,
            'close' : False,
        }
        settings.add_model_settings_to_dict(data)

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
