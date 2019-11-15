"""
defines:
 - LegendObject
"""
import os
from qtpy.QtWidgets import QMainWindow
from pyNastran.gui.menus.legend.qt_legend import LegendPropertiesWindow
from pyNastran.gui.menus.legend.animation import AnimationWindow
from pyNastran.utils.numpy_utils import integer_types


class LegendObject:
    """defines LegendObject, which is an interface to the Legend Window"""
    def __init__(self, gui):
        self.gui = gui
        self._legend_window_shown = False
        self._legend_window = None
        self.is_horizontal_scalar_bar = False
        self.is_low_to_high = True

        self._animation_window_shown = False
        self._animation_window = None

    def show_legend(self):
        """shows the legend"""
        if self._legend_window_shown:
            self._legend_window.show_legend()

    def hide_legend(self):
        """hides the legend"""
        if self._legend_window_shown:
            self._legend_window.hide_legend()

    def clear_legend(self):
        """clears the legend"""
        if self._legend_window_shown:
            self._legend_window.clear()

    def _set_legend_fringe(self, is_fringe):
        if self._legend_window_shown:
            self._legend_window._set_legend_fringe(is_fringe)

    def set_font_size(self, font_size):
        """sets the font size for the legend window"""
        if self._legend_window_shown:
            self._legend_window.set_font_size(font_size)
        if self._animation_window_shown:
            self._animation_window.set_font_size(font_size)

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
        if not hasattr(self.gui, 'case_keys') or len(self.gui.case_keys) == 0:
            self.gui.log_error('No model has been loaded.')
            return

        default_format = None
        (result_type, scalar_bar, defaults_scalar_bar, data_format, default_format,
         default_title, min_value, max_value, default_min, default_max) = self.get_legend_fringe(
             self.gui.icase_fringe)

        nlabels, labelsize, ncolors, colormap = scalar_bar
        default_nlabels, default_labelsize, default_ncolors, default_colormap = defaults_scalar_bar

        scale, phase, default_scale, default_phase = self.get_legend_disp(
            self.gui.icase_disp)

        arrow_scale, default_arrow_scale = self.get_legend_vector(self.gui.icase_vector)

        #arrow_scale = None
        #default_arrow_scale = None

        data = {
            'font_size' : self.settings.font_size,
            #'icase' : self.icase,
            'icase_fringe' : self.gui.icase_fringe,  # int for normals
            'icase_disp' : self.gui.icase_disp,
            'icase_vector' : self.gui.icase_vector,
            'title' : result_type,
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
            'is_horizontal': self.gui.scalar_bar.is_horizontal,
            'is_shown' : self.gui.scalar_bar.is_shown,
            'is_fringe' : self.gui._is_fringe,
            'clicked_ok' : False,
            'close' : False,
        }
        if not isinstance(self.gui, QMainWindow): # pragma: no cover
            return # testing
        if not self._legend_window_shown:
            self._legend_window = LegendPropertiesWindow(data, win_parent=self.gui)
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

    def set_animation_menu(self):
        if not hasattr(self.gui, 'case_keys') or len(self.gui.case_keys) == 0:
            self.gui.log_error('No model has been loaded.')
            return

        default_format = None
        (result_type, scalar_bar, defaults_scalar_bar, data_format, default_format,
         default_title, min_value, max_value, default_min, default_max) = self.get_legend_fringe(
             self.gui.icase_fringe)

        #nlabels, labelsize, ncolors, colormap = scalar_bar
        #default_nlabels, default_labelsize, default_ncolors, default_colormap = defaults_scalar_bar

        scale, phase, default_scale, default_phase = self.get_legend_disp(
            self.gui.icase_disp)

        arrow_scale, default_arrow_scale = self.get_legend_vector(self.gui.icase_vector)

        data = {
            'font_size' : self.settings.font_size,
            'icase_fringe' : self.gui.icase_fringe,
            'icase_disp' : self.gui.icase_disp,
            'icase_vector' : self.gui.icase_vector,
            'title' : result_type,
            'time' : 2,
            'frames/sec' : 30,
            'resolution' : 1,
            'iframe' : 0,
            'scale' : scale,
            'default_scale' : default_scale,

            'arrow_scale' : arrow_scale,
            'default_arrow_scale' : default_arrow_scale,

            'is_scale' : default_phase is None,
            'phase' : phase,
            'default_phase' : default_phase,
            'dirname' : os.path.abspath(os.getcwd()),
            'clicked_ok' : False,
            'close' : False,
        }
        self.set_animation_window(data)

    def set_animation_window(self, data):
        if not self._animation_window_shown:
            self._animation_window = AnimationWindow(
                data, win_parent=self.gui,
                fringe_cases=self.gui.get_form(),
                is_gui_parent=True,
            )
            self._animation_window.show()
            self._animation_window_shown = True
            self._animation_window.exec_()
        else:
            self._animation_window.activateWindow()

        if data['close']:
            if not self._animation_window._updated_animation:
                #self._apply_animation(data)
                pass
            self._animation_window_shown = False
            del self._animation_window
        else:
            self._animation_window.activateWindow()

    def update_legend(self, icase_fringe, icase_disp, icase_vector,
                      name, min_value, max_value, data_format, scale, phase,
                      arrow_scale,
                      nlabels, labelsize, ncolors, colormap,
                      use_fringe_internal=False, use_disp_internal=False,
                      use_vector_internal=False, external_call=True):
        """
        Internal method for updating the legend

        Parameters
        ----------
        icase_fringe : int, None
            int : the active case being shown on the fringe (color plot)
            None : no color
        icase_disp : int, None
            int : the active displacement case
            None : no displacement
        icase_vector : int, None
            int : the active vector case
            None : no vectors
        title : str
            the title of the plot
        min_value / max_value : float, None
            float : the min/max value of the legend
            None : icase_fringe is None
        data_format : str
            a string formatter (e.g., %.4f)
        scale : float, None
            the scaling for the displacement/vector
            None : icase_disp is None
        phase : float, None
            the phase angle for the displacement/vector
            None : icase_disp is None
        arrow_scale : float, None
            the scaling for the vector
            None : icase_vector is None
        nlabels : int
            the number of legend labels
        labelsize : int
            the legend text size
        ncolors : int
            the number of colors on the legend
        colormap : str
            the selected colormap for the fringe
        use_fringe_internal : bool; default=False
            True : use the internal fringe parameters
            False : use the values that were passed in
        use_disp_internal : bool; default=False
            True : use the internal of scale and phase
            False : use the values that were passed in
        use_vector_internal : bool; default=False
            True : use the internal of arrow_scale
            False : use the values that were passed in
        external_call : bool; default=True
            True : allow the legend ``on_apply`` method to be called
            False : the scalar bar/displacement updating will be handled
                    manually to prevent recursion (and a crash)
        """
        if not self._legend_window_shown:
            return
        self._legend_window._updated_legend = True
        is_fringe = self.gui._is_fringe

        out = self.get_legend_fringe(icase_fringe)
        (
            _result_type, scalar_bar, defaults_scalar_bar, data_format,
            default_format, default_title, _min_value, _max_value,
            default_min, default_max) = out

        unused_nlabels, _labelsize, _ncolors, _colormap = scalar_bar
        default_nlabels, default_labelsize, default_ncolors, default_colormap = defaults_scalar_bar
        if use_fringe_internal:
            min_value = _min_value
            max_value = _max_value
            unused_result_type = _result_type
            labelsize = _labelsize
            ncolors = _ncolors
            colormap = _colormap

        #if icase_fringe is not None:
            #key = self.gui.case_keys[icase_fringe]
            #assert isinstance(key, integer_types), key
            #(obj, (i, name)) = self.result_cases[key]
            ##subcase_id = obj.subcase_id
            ##case = obj.get_result(i, name)
            ##result_type = obj.get_title(i, name)
            ##vector_size = obj.get_vector_size(i, name)
            ##location = obj.get_location(i, name)
            ##data_format = obj.get_data_format(i, name)
            ##scale = obj.get_scale(i, name)
            ##label2 = obj.get_header(i, name)
            #default_data_format = obj.get_default_data_format(i, name)
            #default_min, default_max = obj.get_default_min_max(i, name)
            #default_title = obj.get_default_title(i, name)
            #out_labels = obj.get_default_nlabels_labelsize_ncolors_colormap(i, name)
            #default_nlabels, default_labelsize, default_ncolors, default_colormap = out_labels
            #is_normals = obj.is_normal_result(i, name)
            #is_fringe = not is_normals

        _scale, _phase, default_scale, default_phase = self.get_legend_disp(
            icase_disp)
        #if icase_disp is not None:
            #default_scale = obj.get_default_scale(i, name)
            #default_phase = obj.get_default_phase(i, name)
        if use_disp_internal:
            scale = _scale
            phase = _phase
            #default_scale = _default_scale
            #default_phase = _default_phase


        _arrow_scale, default_arrow_scale = self.get_legend_vector(icase_vector)
        if use_vector_internal:
            arrow_scale = _arrow_scale
            #default_arrow_scale = _default_arrow_scale

        #assert isinstance(scale, float), 'scale=%s' % scale
        self._legend_window.update_legend(
            icase_fringe, icase_disp, icase_vector,
            name, min_value, max_value, data_format,
            nlabels, labelsize, ncolors, colormap, is_fringe,
            scale, phase,
            arrow_scale,

            default_title, default_min, default_max, default_format,
            default_nlabels, default_labelsize,
            default_ncolors, default_colormap,
            default_scale, default_phase,
            default_arrow_scale,
            font_size=self.settings.font_size)
        #self.scalar_bar.set_visibility(self._legend_shown)
        #self.vtk_interactor.Render()

    @property
    def settings(self):
        """gets the gui settings"""
        return self.gui.settings

    def _apply_legend(self, data):
        title = data['title']
        min_value = data['min_value']
        max_value = data['max_value']
        scale = data['scale']
        phase = data['phase']
        arrow_scale = data['arrow_scale']
        data_format = data['format']
        is_low_to_high = data['is_low_to_high']
        is_discrete = data['is_discrete']
        is_horizontal = data['is_horizontal']
        is_shown = data['is_shown']

        nlabels = data['nlabels']
        labelsize = data['labelsize']
        ncolors = data['ncolors']
        colormap = data['colormap']

        self.on_update_legend(title=title, min_value=min_value, max_value=max_value,
                              scale=scale, phase=phase,
                              arrow_scale=arrow_scale,
                              data_format=data_format,
                              is_low_to_high=is_low_to_high,
                              is_discrete=is_discrete, is_horizontal=is_horizontal,
                              nlabels=nlabels, labelsize=labelsize,
                              ncolors=ncolors, colormap=colormap,
                              is_shown=is_shown)

    def on_update_legend(self,
                         title='Title', min_value=0., max_value=1.,
                         scale=0.0, phase=0.0,
                         arrow_scale=1.,
                         data_format='%.0f',
                         is_low_to_high=True, is_discrete=True, is_horizontal=True,
                         nlabels=None, labelsize=None, ncolors=None, colormap=None,
                         is_shown=True, render=True):
        """
        Updates the legend/model

        Parameters
        ----------
        scale : float
            displacemnt scale factor; true scale

        TODO: speed up by using existing values to skip update steps
        """
        if colormap is None:
            colormap = self.settings.colormap

        is_shown_old = self.gui.scalar_bar.is_shown
        is_horizontal_old = self.is_horizontal_scalar_bar
        is_low_to_high_old = self.is_low_to_high

        self.is_low_to_high = is_low_to_high
        self.is_horizontal_scalar_bar = is_horizontal

        #print('is_shown2 =', is_shown)
        #assert is_shown == False, is_shown
        is_normal = False
        update_legend = False
        location = 'centroid'
        if self.gui.icase_fringe is not None:
            key = self.gui.case_keys[self.gui.icase_fringe]
            assert isinstance(key, integer_types), key
            (obj, (i, res_name)) = self.gui.result_cases[key]
            subcase_id = obj.subcase_id

            location = obj.get_location(i, res_name)
            min_value_old, max_value_old = obj.get_min_max(i, res_name)
            data_format_old = obj.get_data_format(i, res_name)
            colors_old = obj.get_nlabels_labelsize_ncolors_colormap(i, res_name)
            nlabels_old, labelsize_old, ncolors_old, colormap_old = colors_old

            update_fringe = (
                min_value != min_value_old or
                max_value != max_value_old
            )
            update_legend = (
                (
                    (nlabels, labelsize, ncolors, colormap) !=
                    (nlabels_old, labelsize_old, ncolors_old, colormap_old) or
                    data_format != data_format_old or
                    is_shown != is_shown_old or
                    is_horizontal != is_horizontal_old or
                    is_low_to_high != is_low_to_high_old) and
                not update_fringe)

            obj.set_min_max(i, res_name, min_value, max_value)
            obj.set_data_format(i, res_name, data_format)
            obj.set_nlabels_labelsize_ncolors_colormap(
                i, res_name, nlabels, labelsize, ncolors, colormap)

            #data_format = obj.get_data_format(i, res_name)
            #obj.set_format(i, res_name, data_format)
            #obj.set_data_format(i, res_name, data_format)
            unused_subtitle, unused_label = self.gui.get_subtitle_label(subcase_id)
            is_normal = obj.is_normal_result(i, res_name)
            #if scale != scale_old or phase != phase_old:
            #if not from_legend_menu:
            if update_fringe:
                self.gui.on_fringe(self.gui.icase_fringe, show_msg=False,
                                   update_legend_window=False)

        if is_normal:
            return

        if self.gui.icase_disp is not None:
            key = self.gui.case_keys[self.gui.icase_disp]
            assert isinstance(key, integer_types), key
            (objd, (i, res_name)) = self.gui.result_cases[key]
            scale_old = objd.get_scale(i, res_name)
            phase_old = objd.get_phase(i, res_name)
            update_disp = scale != scale_old or phase != phase_old
            if update_disp:
                objd.set_scale(i, res_name, scale)
                objd.set_phase(i, res_name, phase)
                assert isinstance(scale, float), scale
                self.gui.on_disp(self.gui.icase_disp, apply_fringe=False,
                                 update_legend_window=False, show_msg=False)

        if self.gui.icase_vector is not None:
            key = self.gui.case_keys[self.gui.icase_vector]
            assert isinstance(key, integer_types), key
            (objv, (i, res_name)) = self.gui.result_cases[key]
            arrow_scale_old = objv.get_scale(i, res_name)
            objv.set_scale(i, res_name, arrow_scale)
            assert isinstance(arrow_scale, float), arrow_scale
            update_vector = arrow_scale != arrow_scale_old
            if update_vector:
                self.gui.on_vector(self.gui.icase_vector, apply_fringe=False,
                                   update_legend_window=False, show_msg=False)

        #unused_name = (vector_size1, subcase_id, result_type, label, min_value, max_value, scale1)
        #if obj.is_normal_result(i, res_name):
            #return

        if self.gui.icase_fringe is None:
            return

        #norm_value = float(max_value - min_value)
        # if name not in self._loaded_names:

        #if isinstance(key, integer_types):  # vector 3
             #norm_plot_value = norm(plot_value, axis=1)
            #grid_result = self.set_grid_values(name, norm_plot_value, vector_size1,
                                               #is_low_to_high=is_low_to_high)
        #else:
        if update_legend:
            self.gui.update_scalar_bar(title, min_value, max_value,
                                       data_format,
                                       nlabels=nlabels, labelsize=labelsize,
                                       ncolors=ncolors, colormap=colormap,
                                       is_shown=is_shown)
            self.gui.update_contour_filter(nlabels, location, min_value, max_value)
        if render:
            self.gui.Render()

        msg = (
            f'self.on_update_legend(title={title!r}, min_value={min_value}, max_value={max_value},\n'
            f'                      scale={scale}, phase={phase},\n'
            f'                      data_format={data_format!r}, is_low_to_high={is_low_to_high}, '
            f'is_discrete={is_discrete},\n'
            f'                      nlabels={nlabels}, labelsize={labelsize}, '
            f'ncolors={ncolors}, colormap={colormap!r},\n'
            f'                      is_horizontal={is_horizontal}, is_shown={is_shown})'
        )
        self.gui.log_command(msg)
        #if is_shown:
            #pass

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
            key = self.gui.case_keys[icase_fringe]
            assert isinstance(key, integer_types), key
            (obj, (i, res_name)) = self.gui.result_cases[key]
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
            key = self.gui.case_keys[icase_disp]
            (objd, (i, res_name)) = self.gui.result_cases[key]
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
            key = self.gui.case_keys[icase_vector]
            (objv, (i, res_name)) = self.gui.result_cases[key]
            arrow_scale = objv.get_scale(i, res_name)
            default_arrow_scale = objv.get_default_scale(i, res_name)
            #phasev = objv.get_phase(i, res_name)
            #default_phasev = objv.get_default_phase(i, res_name)
        return arrow_scale, default_arrow_scale
