"""
defines:
 - ClippingObject

"""
from pyNastran.gui.menus.clipping.clipping import ClippingPropertiesWindow

class ClippingObject:
    """defines ClippingObject"""
    def __init__(self, gui):
        """creates ClippingObject"""
        self.gui = gui
        self._clipping_window_shown = False
        self._clipping_window = None

    def set_font_size(self, font_size):
        """sets the font size for the edit geometry properties window"""
        if self._clipping_window_shown:
            self._clipping_window.set_font_size(font_size)

    def set_clipping_menu(self):
        """
        Opens a dialog box to set:

        +--------+----------+
        |  Min   |  Float   |
        +--------+----------+
        |  Max   |  Float   |
        +--------+----------+
        """
        #if not hasattr(self, 'case_keys'):  # TODO: maybe include...
            #self.log_error('No model has been loaded.')
            #return
        camera = self.gui.GetCamera()
        min_clip, max_clip = camera.GetClippingRange()

        data = {
            'font_size' : self.gui.settings.font_size,
            'min_clip' : min_clip,
            'max_clip' : max_clip,
            'clicked_ok' : False,
            'close' : False,
        }
        if not self._clipping_window_shown:
            self._clipping_window = ClippingPropertiesWindow(data, win_parent=self.gui)
            self._clipping_window.show()
            self._clipping_window_shown = True
            self._clipping_window.exec_()
        else:
            self._clipping_window.activateWindow()

        if data['close']:
            if not self._clipping_window._updated_clipping:
                self.apply_clipping(data)
            del self._clipping_window
            self._clipping_window_shown = False
        else:
            self._clipping_window.activateWindow()

    def apply_clipping(self, data):
        min_clip = data['min_clip']
        max_clip = data['max_clip']
        self.on_update_clipping(min_clip, max_clip)

    def on_update_clipping(self, min_clip=None, max_clip=None):
        camera = self.gui.GetCamera()
        _min_clip, _max_clip = camera.GetClippingRange()
        if min_clip is None:
            min_clip = _min_clip
        if max_clip is None:
            max_clip = _max_clip
        camera.SetClippingRange(min_clip, max_clip)
        self.gui.log_command('on_update_clipping(min_clip=%s, max_clip=%s)'
                             % (min_clip, max_clip))
