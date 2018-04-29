from __future__ import print_function
import vtk
from pyNastran.gui.menus.preferences.preferences import PreferencesWindow


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

    camera = self.GetCamera()
    min_clip, max_clip = camera.GetClippingRange()
    data = {
        'font_size' : self.settings.font_size,
        'annotation_size' : self.settings.annotation_size, # int
        'annotation_color' : self.settings.annotation_color,

        'use_gradient_background' : self.settings.use_gradient_background,
        'background_color' : self.settings.background_color,
        'background_color2' : self.settings.background_color2,

        'text_size' : self.settings.text_size,
        'text_color' : self.settings.text_color,

        'picker_size' : self.element_picker_size,
        'dim_max' : self.settings.dim_max,
        'coord_scale' : self.settings.coord_scale,
        'coord_text_scale' : self.settings.coord_text_scale,
        'magnify' : self.settings.magnify,

        'min_clip' : min_clip,
        'max_clip' : max_clip,

        'clicked_ok' : False,
        'close' : False,
    }
    if not self._preferences_window_shown:
        self._preferences_window = PreferencesWindow(data, win_parent=self)
        self._preferences_window.show()
        self._preferences_window_shown = True
        self._preferences_window.exec_()
    else:
        self._preferences_window.activateWindow()

    if data['close']:
        if not self._preferences_window._updated_preference:
            self.settings.on_set_font_size(data['font_size'])
        del self._preferences_window
        self._preferences_window_shown = False
    else:
        self._preferences_window.activateWindow()
