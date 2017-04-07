from __future__ import print_function
from six import integer_types
from pyNastran.gui.gui_interface.legend.qt_legend import LegendPropertiesWindow

def set_legend_menu(self):
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
    if not hasattr(self, 'case_keys') or len(self.case_keys) == 0:
        self.log_error('No model has been loaded.')
        return
    key = self.case_keys[self.icase]
    default_format = None
    assert isinstance(key, integer_types), key
    (obj, (i, res_name)) = self.result_cases[key]
    case = obj.get_result(i, res_name)
    result_type = obj.get_title(i, res_name)
    nlabels, labelsize, ncolors, colormap = obj.get_nlabels_labelsize_ncolors_colormap(i, res_name)

    defaults_scalar_bar = obj.get_default_nlabels_labelsize_ncolors_colormap(i, res_name)
    default_nlabels, default_labelsize, default_ncolors, default_colormap = defaults_scalar_bar

    data_format = obj.get_data_format(i, res_name)
    scale = obj.get_scale(i, res_name)
    phase = obj.get_phase(i, res_name)

    default_title = obj.get_default_title(i, res_name)
    default_scale = obj.get_default_scale(i, res_name)
    default_phase = obj.get_default_phase(i, res_name)

    min_value, max_value = obj.get_min_max(i, res_name)
    default_min, default_max = obj.get_default_min_max(i, res_name)

    if default_format is None:
        default_format = data_format

    data = {
        'font_size' : self.font_size,
        'icase' : self.icase,
        'name' : result_type,
        'min' : min_value,
        'max' : max_value,

        'scale' : scale,
        'phase' : phase,
        'format' : data_format,

        'default_min' : default_min,
        'default_max' : default_max,
        'default_title' : default_title,
        'default_scale' : default_scale,
        'default_phase' : default_phase,
        'default_format' : default_format,

        'default_nlabels' : default_nlabels,
        'default_labelsize' : default_labelsize,
        'default_ncolors' : default_ncolors,
        'default_colormap' : default_colormap,

        'nlabels' : nlabels,
        'labelsize' :  labelsize,
        'ncolors' : ncolors,
        'colormap' : colormap,

        'is_low_to_high' : True,
        'is_discrete': True,
        'is_horizontal': self.scalar_bar.is_horizontal,
        'is_shown' : self.scalar_bar.is_shown,
        'is_normals' : self._is_normals,
        'clicked_ok' : False,
        'close' : False,
    }
    if not self._legend_window_shown:
        self._legend_window = LegendPropertiesWindow(data, win_parent=self)
        self._legend_window.show()
        self._legend_window_shown = True
        self._legend_window.exec_()
    else:
        self._legend_window.activateWindow()

    if data['close']:
        if not self._legend_window._updated_legend:
            self._apply_legend(data)
        self._legend_window_shown = False
        del self._legend_window
    else:
        self._legend_window.activateWindow()
