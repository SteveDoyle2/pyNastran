from __future__ import print_function
from pyNastran.gui.menus.preferences.preferences import PreferencesWindow


class PreferencesObject(object):
    def __init__(self, gui):
        self.gui = gui
        self._preferences_window_shown = False
        self._preferences_window = None

    def set_font_size(self, font_size):
        """sets the font size for the preferences window"""
        if self._preferences_window_shown:
            self._preferences_window.set_font_size(font_size)

    def set_preferences_menu(self):
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
        settings = self.gui.settings
        data = {
            'font_size' : settings.font_size,
            'annotation_size' : settings.annotation_size, # int
            'annotation_color' : settings.annotation_color,

            'use_gradient_background' : settings.use_gradient_background,
            'background_color' : settings.background_color,
            'background_color2' : settings.background_color2,

            'text_size' : settings.text_size,
            'text_color' : settings.text_color,

            'highlight_color' : settings.highlight_color,
            'highlight_opacity' : settings.highlight_opacity,

            'picker_size' : self.gui.element_picker_size,
            'dim_max' : settings.dim_max,
            'coord_scale' : settings.coord_scale,
            'coord_text_scale' : settings.coord_text_scale,
            'show_corner_coord' : self.gui.get_corner_axis_visiblity(),
            'magnify' : settings.magnify,

            'min_clip' : min_clip,
            'max_clip' : max_clip,

            'nastran_is_element_quality' : settings.nastran_is_element_quality,
            'nastran_is_properties' : settings.nastran_is_properties,
            'nastran_is_bar_axes' : settings.nastran_is_bar_axes,
            'nastran_is_3d_bars' : settings.nastran_is_3d_bars,
            'nastran_is_3d_bars_update' : settings.nastran_is_3d_bars_update,
            'nastran_create_coords' : settings.nastran_create_coords,

            'clicked_ok' : False,
            'close' : False,
        }
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
