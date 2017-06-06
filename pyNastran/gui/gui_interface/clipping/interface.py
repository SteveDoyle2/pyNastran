from __future__ import print_function
from pyNastran.gui.gui_interface.clipping.clipping import ClippingPropertiesWindow

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
    camera = self.GetCamera()
    min_clip, max_clip = camera.GetClippingRange()

    data = {
        'font_size' : self.font_size,
        'clipping_min' : min_clip,
        'clipping_max' : max_clip,
        'clicked_ok' : False,
        'close' : False,
    }
    if not self._clipping_window_shown:
        self._clipping_window = ClippingPropertiesWindow(data, win_parent=self)
        self._clipping_window.show()
        self._clipping_window_shown = True
        self._clipping_window.exec_()
    else:
        self._clipping_window.activateWindow()

    if data['close']:
        if not self._clipping_window._updated_clipping:
            self._apply_clipping(data)
        del self._clipping_window
        self._clipping_window_shown = False
    else:
        self._clipping_window.activateWindow()

