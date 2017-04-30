from __future__ import print_function
from pyNastran.gui.gui_interface.preferences.preferences import PreferencesWindow

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
        'font_size' : self.font_size,
        'label_size' : self.label_text_size,
        'label_color' : self.label_color,
        'background_color' : self.background_color,
        'text_color' : self.text_color,

        'picker_size' : self.element_picker_size,
        'dim_max' : self.dim_max,

        'clipping_min' : min_clip,
        'clipping_max' : max_clip,

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
            self.on_set_font_size(data['font_size'])
        del self._preferences_window
        self._preferences_window_shown = False
    else:
        self._preferences_window.activateWindow()

