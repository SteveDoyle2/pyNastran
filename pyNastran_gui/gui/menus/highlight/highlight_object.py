from __future__ import print_function
from pyNastran.gui.menus.highlight.highlight import HighlightWindow


class HighlightObject(object):
    def __init__(self, gui):
        self.gui = gui
        self._highlight_window_shown = False
        self._highlight_window = None

    def set_font_size(self, font_size):
        """sets the font size for the preferences window"""
        if self._highlight_window_shown:
            self._highlight_window.set_font_size(font_size)

    def set_highlight_menu(self):
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

        model_name = 'main'
        data = {
            'font_size' : settings.font_size,

            'model_name' : model_name,
            'highlight_color' : settings.highlight_color,
            'highlight_opacity' : settings.highlight_opacity,

            'clicked_ok' : False,
            'close' : False,
        }
        if not self._highlight_window_shown:
            self._highlight_window = HighlightWindow(data, win_parent=self.gui)
            self._highlight_window.show()
            self._highlight_window_shown = True
            self._highlight_window.exec_()
        else:
            self._highlight_window.activateWindow()

        if data['close']:
            if not self._highlight_window._updated_highlight:
                settings.on_set_font_size(data['font_size'])
            del self._highlight_window
            self._highlight_window_shown = False
        else:
            self._highlight_window.activateWindow()
