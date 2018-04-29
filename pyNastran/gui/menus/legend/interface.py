"""
defines:
 - _get_fringe(icase_fringe)
 - set_legend_menu(self)
"""
from __future__ import print_function
from six import integer_types
from pyNastran.gui.menus.legend.qt_legend import LegendPropertiesWindow

def get_legend_fringe(self, icase_fringe):
    """helper method for ``set_legend_menu``"""
    #nlabels = None
    #labelsize = None
    #ncolors = None
    #colormap = None
    result_type = None
    min_value = None
    max_value = None
    data_format = None
    scalar_bar = (None, None, None, None)
    defaults_scalar_bar = (None, None, None, None)

    #default_nlabels = None
    #default_labelsize = None
    #default_ncolors = None
    #default_colormap = None

    default_min = None
    default_max = None
    default_format = None
    default_title = None

    #title = None
    if icase_fringe is not None:
        key = self.case_keys[icase_fringe]
        assert isinstance(key, integer_types), key
        (obj, (i, res_name)) = self.result_cases[key]
        #case = obj.get_result(i, res_name)
        result_type = obj.get_title(i, res_name)

        scalar_bar = obj.get_nlabels_labelsize_ncolors_colormap(i, res_name)
        defaults_scalar_bar = obj.get_default_nlabels_labelsize_ncolors_colormap(i, res_name)

        data_format = obj.get_data_format(i, res_name)
        default_title = obj.get_default_title(i, res_name)
        min_value, max_value = obj.get_min_max(i, res_name)
        default_min, default_max = obj.get_default_min_max(i, res_name)
        default_format = obj.get_default_data_format(i, res_name)

    out = (
        result_type, scalar_bar, defaults_scalar_bar, data_format, default_format,
        default_title, min_value, max_value, default_min, default_max)
    return out


def get_legend_disp(self, icase_disp):
    """helper method for ``set_legend_menu``"""
    scale = None
    phase = None
    default_scale = None
    default_phase = None
    if icase_disp is not None:
        key = self.case_keys[icase_disp]
        (objd, (i, res_name)) = self.result_cases[key]
        scale = objd.get_scale(i, res_name)
        phase = objd.get_phase(i, res_name)
        default_scale = objd.get_default_scale(i, res_name)
        default_phase = objd.get_default_phase(i, res_name)
    return scale, phase, default_scale, default_phase

def get_legend_vector(self, icase_vector):
    """helper method for ``set_legend_menu``"""
    arrow_scale = None
    default_arrow_scale = None
    if icase_vector is not None:
        key = self.case_keys[icase_vector]
        (objv, (i, res_name)) = self.result_cases[key]
        arrow_scale = objv.get_scale(i, res_name)
        default_arrow_scale = objv.get_default_scale(i, res_name)
        #phasev = objv.get_phase(i, res_name)
        #default_phasev = objv.get_default_phase(i, res_name)
    return arrow_scale, default_arrow_scale


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

    default_format = None
    (result_type, scalar_bar, defaults_scalar_bar, data_format, default_format,
     default_title, min_value, max_value, default_min, default_max) = get_legend_fringe(
         self, self.icase_fringe)

    nlabels, labelsize, ncolors, colormap = scalar_bar
    default_nlabels, default_labelsize, default_ncolors, default_colormap = defaults_scalar_bar

    scale, phase, default_scale, default_phase = get_legend_disp(
        self, self.icase_disp)

    arrow_scale, default_arrow_scale = get_legend_vector(self, self.icase_vector)

    #arrow_scale = None
    #default_arrow_scale = None

    data = {
        'font_size' : self.settings.font_size,
        #'icase' : self.icase,
        'icase_fringe' : self.icase_fringe,  # int for normals
        'icase_disp' : self.icase_disp,
        'icase_vector' : self.icase_vector,
        'name' : result_type,
        'min_value' : min_value,
        'max_value' : max_value,

        'scale' : scale,
        'arrow_scale' : arrow_scale,
        'phase' : phase,
        'format' : data_format,

        'default_min' : default_min,
        'default_max' : default_max,
        'default_title' : default_title,
        'default_scale' : default_scale,
        'default_arrow_scale' : default_arrow_scale,
        'default_phase' : default_phase,
        'default_format' : default_format,

        'default_nlabels' : default_nlabels,
        'default_labelsize' : default_labelsize,
        'default_ncolors' : default_ncolors,
        'default_colormap' : default_colormap,

        'nlabels' : nlabels,
        'labelsize' : labelsize,
        'ncolors' : ncolors,
        'colormap' : colormap,

        'is_low_to_high' : True,
        'is_discrete': True,
        'is_horizontal': self.scalar_bar.is_horizontal,
        'is_shown' : self.scalar_bar.is_shown,
        'is_fringe' : self._is_fringe,
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
