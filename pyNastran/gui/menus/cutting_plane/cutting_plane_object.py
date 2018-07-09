from __future__ import print_function
from pyNastran.gui.menus.cutting_plane.cutting_plane import CuttingPlaneWindow


class CuttingPlaneObject(object):
    def __init__(self, gui):
        self.gui = gui
        self._cutting_plane_window_shown = False
        self._cutting_plane_window = None

    def set_font_size(self, font_size):
        """sets the font size for the preferences window"""
        if self._cutting_plane_window_shown:
            self._cutting_plane_window.set_font_size(font_size)

    def set_cutting_plane_menu(self):
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
        model_name = self.gui.name
        model = self.gui.models[model_name]
        if hasattr(model, 'coords'):
            cids = list(model.coords.keys())
            cids.sort()
        else:
            cids = [0]

        data = {
            'font_size' : settings.font_size,
            'cids' : cids,
            'plane_color' : (1., 0., 1.), # purple
            'name' : model_name,

            'clicked_ok' : False,
            'close' : False,
        }
        if not self._cutting_plane_window_shown:
            self._cutting_plane_window = CuttingPlaneWindow(data, win_parent=self.gui)
            self._cutting_plane_window.show()
            self._cutting_plane_window = True
            self._cutting_plane_window.exec_()
        else:
            self._cutting_plane_window.activateWindow()

        if data['close']:
            #if not self._cutting_plane_window._updated_preference:
            #    settings.on_set_font_size(data['font_size'])
            del self._cutting_plane_window
            self._cutting_plane_window_shown = False
        else:
            self._cutting_plane_window.activateWindow()
