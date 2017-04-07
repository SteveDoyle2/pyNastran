from __future__ import print_function
from pyNastran.gui.gui_interface.modify_label_properties.modify_label_properties import (
    ModifyLabelPropertiesMenu)


def on_set_labelsize_color_menu(self):
    """
    Opens a dialog box to set:

    +--------+----------+
    |  Name  |  String  |
    +--------+----------+
    |  Min   |  Float   |
    +--------+----------+
    |  Max   |  Float   |
    +--------+----------+
    | Format | pyString |
    +--------+----------+
    """
    if not hasattr(self, 'case_keys'):
        self.log_error('No model has been loaded.')
        return

    data = {
        'font_size' : self.font_size,
        'size' : self.label_text_size,
        'color' : self.label_color,
        'dim_max' : self.dim_max,
        #'clicked_ok' : False,
        #'clicked_cancel' : False,
        #'close' : False,
    }
    #print(data)
    if not self._label_window_shown:
        self._label_window = ModifyLabelPropertiesMenu(data, win_parent=self)
        self._label_window.show()
        self._label_window_shown = True
        self._label_window.exec_()
    else:
        self._label_window.activateWindow()

    if 'close' not in data:
        self._label_window.activateWindow()
        return

    if data['close']:
        self._label_window_shown = False
        del self._label_window
    else:
        self._label_window.activateWindow()
