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
    data = {
        'text_size' : self.font_size,
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
            self._apply_preferences(data)
        del self._preferences_window
        self._preferences_window_shown = False
    else:
        self._preferences_window.activateWindow()

