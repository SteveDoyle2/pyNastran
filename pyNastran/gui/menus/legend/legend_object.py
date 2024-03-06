"""
defines:
 - LegendObject

"""
from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
import os
from qtpy.QtWidgets import QMainWindow
from pyNastran.gui.menus.legend.qt_legend import LegendPropertiesWindow
from pyNastran.gui.menus.legend.animation import AnimationWindow
from pyNastran.gui.qt_files.base_gui import BaseGui
from pyNastran.utils.numpy_utils import integer_types
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow


class LegendObject(BaseGui):
    """defines LegendObject, which is an interface to the Legend Window"""
    def __init__(self, gui):
        super().__init__(gui)
        #self.gui = gui
        self._legend_window_shown = False
        self._legend_window = None
        self.is_low_to_high = True

        self._animation_window_shown = False
        self._animation_window = None

    def show_legend(self) -> None:
        """shows the legend"""
        if self._legend_window_shown:
            self._legend_window.show_legend()

    def hide_legend(self) -> None:
        """hides the legend"""
        if self._legend_window_shown:
            self._legend_window.hide_legend()

    def clear_legend(self) -> None:
        """clears the legend"""
        if self._legend_window_shown:
            self._legend_window.clear()

    def _set_legend_fringe(self, is_fringe: int) -> None:
        if self._legend_window_shown:
            self._legend_window._set_legend_fringe(is_fringe)

    def set_font_size(self, font_size: int) -> None:
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
        gui = self.gui
        if not hasattr(gui, 'case_keys') or len(gui.case_keys) == 0:
            gui.log_error('No model has been loaded.')
            return

        default_format = None
        (is_method_array,
         legend_title, scalar_bar, defaults_scalar_bar,
         data_format, default_format,
         default_title, min_value, max_value,
         default_min, default_max) = self.get_legend_fringe(gui.icase_fringe)

        nlabels, labelsize, ncolors, colormap = scalar_bar
        default_nlabels, default_labelsize, default_ncolors, default_colormap = defaults_scalar_bar

        scale, phase, default_scale, default_phase = self.get_legend_disp(
            gui.icase_disp)

        arrow_scale, default_arrow_scale = self.get_legend_vector(gui.icase_vector)

        #arrow_scale = None
        #default_arrow_scale = None


        data = {
            'font_size' : self.settings.font_size,
            #'icase' : self.icase,
            'icase_fringe' : gui.icase_fringe,  # int for normals
            'icase_disp' : gui.icase_disp,
            'icase_vector' : gui.icase_vector,
            'title' : legend_title,
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
            'is_horizontal': gui.scalar_bar.is_horizontal,
            'is_shown' : gui.scalar_bar.is_shown,
            'is_fringe' : gui._is_fringe,
            'clicked_ok' : False,
            'close' : False,
        }
        gui: MainWindow = self.gui
        if not isinstance(gui, QMainWindow): # pragma: no cover
            return # testing
        if not self._legend_window_shown:
            self._legend_window = LegendPropertiesWindow(data, win_parent=gui)
            self._legend_window.show()
            self._legend_window_shown = True
            self._legend_window.exec_()
        else:
            self._legend_window.activateWindow()

        if data['close']:
            if not self._legend_window._updated_legend:
                self.apply_legend(data)
            self._legend_window_shown = False
            del self._legend_window
        else:
            self._legend_window.activateWindow()

    def set_animation_menu(self) -> None:
        gui: QMainWindow = self.gui
        if not hasattr(gui, 'case_keys') or len(gui.case_keys) == 0:
            gui.log_error('No model has been loaded.')
            return

        default_format = None
        (is_method_array,
         legend_title, scalar_bar, defaults_scalar_bar,
         data_format, default_format,
         default_title, min_value, max_value,
         default_min, default_max) = self.get_legend_fringe(gui.icase_fringe)

        # set the title as the default
        if is_method_array and legend_title == default_title:
            legend_title = ''

        #nlabels, labelsize, ncolors, colormap = scalar_bar
        #default_nlabels, default_labelsize, default_ncolors, default_colormap = defaults_scalar_bar

        scale, phase, default_scale, default_phase = self.get_legend_disp(
            gui.icase_disp)

        arrow_scale, default_arrow_scale = self.get_legend_vector(gui.icase_vector)

        data = {
            'font_size' : self.settings.font_size,
            'icase_fringe' : gui.icase_fringe,
            'icase_disp' : gui.icase_disp,
            'icase_vector' : gui.icase_vector,
            'title' : legend_title,
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

    def set_animation_window(self, data) -> None:
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

    def update_legend(self,
                      icase_fringe: Optional[int],
                      icase_disp: Optional[int],
                      icase_vector: Optional[int],
                      title: str,
                      min_value: Optional[float],
                      max_value: Optional[float],
                      data_format: str,
                      scale: Optional[float],
                      phase: Optional[float],
                      arrow_scale: Optional[float],
                      nlabels: int, labelsize: int,
                      ncolors: int, colormap: str,
                      use_fringe_internal: bool=False,
                      use_disp_internal: bool=False,
                      use_vector_internal: bool=False,
                      external_call: bool=True) -> None:
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
        if not self._legend_window_shown:  # pragma: no cover
            return
        gui: MainWindow = self.gui
        self._legend_window._updated_legend = True
        is_fringe = gui._is_fringe

        out = self.get_legend_fringe(icase_fringe)
        (
            is_method_array,
            _legend_title, scalar_bar, defaults_scalar_bar,
            data_format, default_format, default_title,
            _min_value, _max_value,
            default_min, default_max) = out
        if is_method_array and title == default_title:
            title = ''

        unused_nlabels, _labelsize, _ncolors, _colormap = scalar_bar
        default_nlabels, default_labelsize, default_ncolors, default_colormap = defaults_scalar_bar
        if use_fringe_internal:
            min_value = _min_value
            max_value = _max_value
            unused_legend_title = _legend_title
            labelsize = _labelsize
            ncolors = _ncolors
            colormap = _colormap

        #if icase_fringe is not None:
            #key = gui.case_keys[icase_fringe]
            #assert isinstance(key, integer_types), key
            #(obj, (i, name)) = self.result_cases[key]
            ##subcase_id = obj.subcase_id
            ##fringe, case = obj.get_fringe_vector_result(i, name)
            ##legend_title = obj.get_legend_title(i, name)
            ##vector_size = obj.get_vector_size(i, name)
            ##location = obj.get_location(i, name)
            ##data_format = obj.get_data_format(i, name)
            ##scale = obj.get_scale(i, name)
            ##label2 = obj.get_annotation(i, name)
            #default_data_format = obj.get_default_data_format(i, name)
            #default_min, default_max = obj.get_default_min_max(i, name)
            #default_title = obj.get_default_legend_title(i, name)
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
            title, min_value, max_value, data_format,
            nlabels, labelsize, ncolors, colormap, is_fringe,
            scale, phase,
            arrow_scale,

            default_title, default_min, default_max, default_format,
            default_nlabels, default_labelsize,
            default_ncolors, default_colormap,
            default_scale, default_phase,
            default_arrow_scale,
            is_method_array=is_method_array,
            font_size=self.settings.font_size)
        #self.scalar_bar.set_visibility(self._legend_shown)
        #self.vtk_interactor.Render()

    def apply_legend(self, data: dict[str, Any]) -> None:
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

        self.on_update_legend(
            title=title,
            min_value=min_value, max_value=max_value,
            scale=scale, phase=phase,
            arrow_scale=arrow_scale,
            data_format=data_format,
            is_low_to_high=is_low_to_high,
            is_discrete=is_discrete,
            is_horizontal=is_horizontal,
            nlabels=nlabels, labelsize=labelsize,
            ncolors=ncolors, colormap=colormap,
            is_shown=is_shown)

    def on_update_legend(self,
                         title: str='Title',
                         min_value: float=0.,
                         max_value: float=1.,
                         scale: float=0.0,
                         phase: float=0.0,
                         arrow_scale: float=1.,
                         data_format: str='%.0f',
                         is_low_to_high: bool=True,
                         is_discrete: bool=True,
                         is_horizontal: bool=True,
                         nlabels=None, labelsize=None,
                         ncolors=None, colormap=None,
                         is_shown: bool=True,
                         render: bool=True) -> None:
        """
        Updates the legend/model

        Parameters
        ----------
        scale : float
            displacemnt scale factor; true scale

        TODO: speed up by using existing values to skip update steps
        """
        nlabels = None if nlabels == -1 else nlabels
        labelsize = None if labelsize == -1 else labelsize
        ncolors = None if ncolors == -1 else ncolors
        #colormap
        gui: MainWindow = self.gui
        if colormap is None:
            colormap = self.settings.colormap

        is_shown_old = gui.scalar_bar.is_shown
        is_horizontal_old = gui.settings.is_horizontal_scalar_bar
        is_low_to_high_old = self.is_low_to_high

        self.is_low_to_high = is_low_to_high

        #print('is_shown2 =', is_shown)
        #assert is_shown == False, is_shown
        is_normal = False
        update_legend = False
        location = 'centroid'

        if gui.icase_disp is not None:
            key = gui.case_keys[gui.icase_disp]
            assert isinstance(key, integer_types), key
            (objd, (i, res_name)) = gui.result_cases[key]
            scale_old = objd.get_scale(i, res_name)
            phase_old = objd.get_phase(i, res_name)
            update_disp = scale != scale_old or phase != phase_old
            if update_disp:
                objd.set_scale(i, res_name, scale)
                objd.set_phase(i, res_name, phase)
                assert isinstance(scale, float), scale
                gui.on_disp(gui.icase_disp, apply_fringe=False,
                            update_legend_window=False, show_msg=False)

        if gui.icase_vector is not None:
            key = gui.case_keys[gui.icase_vector]
            assert isinstance(key, integer_types), key
            (objv, (i, res_name)) = gui.result_cases[key]
            #arrow_scale_old = objv.get_arrow_scale(i, res_name)
            #objv.set_arrow_scale(i, res_name, arrow_scale)
            arrow_scale_old = objv.get_scale(i, res_name)
            objv.set_scale(i, res_name, arrow_scale)
            assert isinstance(arrow_scale, float), arrow_scale
            update_vector = arrow_scale != arrow_scale_old
            if update_vector:
                gui.on_vector(gui.icase_vector, apply_fringe=False,
                              update_legend_window=False, show_msg=False)

        if gui.icase_fringe is not None:
            key = gui.case_keys[gui.icase_fringe]
            assert isinstance(key, integer_types), key
            (obj, (i, res_name)) = gui.result_cases[key]
            subcase_id = obj.subcase_id

            location = obj.get_location(i, res_name)
            min_value_old, max_value_old = obj.get_min_max(i, res_name)
            data_format_old = obj.get_data_format(i, res_name)
            colors_old = obj.get_nlabels_labelsize_ncolors_colormap(i, res_name)
            nlabels_old, labelsize_old, ncolors_old, colormap_old = colors_old
            title_old = obj.get_legend_title(i, res_name)

            update_fringe = (
                min_value != min_value_old or
                max_value != max_value_old
            )
            update_legend = (
                (
                    (nlabels, labelsize, ncolors, colormap) !=
                    (nlabels_old, labelsize_old, ncolors_old, colormap_old) or
                    title != title_old or
                    data_format != data_format_old or
                    is_shown != is_shown_old or
                    is_horizontal != is_horizontal_old or
                    is_low_to_high != is_low_to_high_old) and
                not update_fringe)
            gui.settings.is_horizontal_scalar_bar = is_horizontal

            try:
                obj.set_min_max(i, res_name, min_value, max_value)
            except TypeError:
                gui.log_error(f'Error setting min/max; i={i} res_name={res_name} '
                              f'min_value={min_value} max_value={max_value}\nobj={str(obj)}')

            ## TODO: Allows you to break the NodeID / GuiResult
            ##       but lets the fancy tables to use the default.
            #if title == '':
            obj.set_legend_title(i, res_name, title)

            # That legend_title is not what is shown.
            # It's setting the default flag, so overwrite it :)
            title = obj.get_legend_title(i, res_name)

            obj.set_data_format(i, res_name, data_format)
            obj.set_nlabels_labelsize_ncolors_colormap(
                i, res_name, nlabels, labelsize, ncolors, colormap)
            is_normal = obj.is_normal_result(i, res_name)

            #data_format = obj.get_data_format(i, res_name)
            #obj.set_format(i, res_name, data_format)
            #obj.set_data_format(i, res_name, data_format)
            unused_subtitle, unused_label = gui.get_subtitle_label(subcase_id)
            #if scale != scale_old or phase != phase_old:
            #if not from_legend_menu:
            if update_fringe:
                gui.on_fringe(gui.icase_fringe, show_msg=False,
                              update_legend_window=False)

        if is_normal:
            return

        #unused_name = (vector_size1, subcase_id, legend_title, label,
                       #min_value, max_value, scale1)
        #if obj.is_normal_result(i, res_name):
            #return

        if gui.icase_fringe is None:
            return

        #norm_value = float(max_value - min_value)
        # if name not in self._loaded_names:

        #if isinstance(key, integer_types):  # vector 3
             #norm_plot_value = norm(plot_value, axis=1)
            #grid_result = self.set_vtk_fringe(name, norm_plot_value, vector_size1,
                                              #is_low_to_high=is_low_to_high)
        #else:
        if update_legend:
            gui.update_scalar_bar(title, min_value, max_value,
                                  data_format,
                                  nlabels=nlabels, labelsize=labelsize,
                                  ncolors=ncolors, colormap=colormap,
                                  is_horizontal=is_horizontal,
                                  is_shown=is_shown)
            gui.update_contour_filter(nlabels, location, min_value, max_value)
        if render:
            gui.Render()

        msg = (
            f'self.on_update_legend(title={title!r}, min_value={min_value}, max_value={max_value},\n'
            f'                      scale={scale}, phase={phase},\n'
            f'                      data_format={data_format!r}, is_low_to_high={is_low_to_high}, '
            f'is_discrete={is_discrete},\n'
            f'                      nlabels={nlabels}, labelsize={labelsize}, '
            f'ncolors={ncolors}, colormap={colormap!r},\n'
            f'                      is_horizontal={is_horizontal}, is_shown={is_shown})'
        )
        gui.log_command(msg)
        #if is_shown:
            #pass

    def get_legend_fringe(self, icase_fringe: Optional[int]) -> tuple:
        """helper method for ``set_legend_menu``"""
        #nlabels = None
        #labelsize = None
        #ncolors = None
        #colormap = None
        is_method_array = False
        legend_title = None
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
        gui = self.gui
        if icase_fringe is not None:
            key = gui.case_keys[icase_fringe]
            assert isinstance(key, integer_types), key
            (obj, (i, res_name)) = gui.result_cases[key]
            is_method_array = obj.is_method_array
            #fringe, case = obj.get_fringe_vector_result(i, res_name)
            legend_title = obj.get_legend_title(i, res_name)

            scalar_bar = obj.get_nlabels_labelsize_ncolors_colormap(i, res_name)
            defaults_scalar_bar = obj.get_default_nlabels_labelsize_ncolors_colormap(i, res_name)
            unused_method = obj.get_methods(i, res_name)[0]

            data_format = obj.get_data_format(i, res_name)
            default_title = obj.get_default_legend_title(i, res_name)
            min_value, max_value = obj.get_min_max(i, res_name)
            default_min, default_max = obj.get_default_min_max(i, res_name)
            default_format = obj.get_default_data_format(i, res_name)

        out = (
            is_method_array,
            legend_title,
            scalar_bar, defaults_scalar_bar,
            data_format, default_format,
            default_title, min_value, max_value,
            default_min, default_max)
        return out

    def get_legend_disp(self, icase_disp: Optional[int],
                        ) -> tuple[Optional[float], Optional[float],
                                   Optional[float], Optional[float]]:
        """helper method for ``set_legend_menu``"""
        scale = None
        phase = None
        default_scale = None
        default_phase = None
        gui: MainWindow = self.gui
        if icase_disp is not None:
            key = gui.case_keys[icase_disp]
            (objd, (i, res_name)) = gui.result_cases[key]
            scale = objd.get_scale(i, res_name)
            phase = objd.get_phase(i, res_name)
            default_scale = objd.get_default_scale(i, res_name)
            default_phase = objd.get_default_phase(i, res_name)
        return scale, phase, default_scale, default_phase

    def get_legend_vector(self, icase_vector: Optional[int],
                          ) -> tuple[Optional[float], Optional[float]]:
        """helper method for ``set_legend_menu``"""
        arrow_scale = None
        default_arrow_scale = None
        gui: MainWindow = self.gui
        if icase_vector is not None:
            key = gui.case_keys[icase_vector]
            (objv, (i, res_name)) = gui.result_cases[key]
            arrow_scale = objv.get_scale(i, res_name)
            default_arrow_scale = objv.get_default_scale(i, res_name)
            #phasev = objv.get_phase(i, res_name)
            #default_phasev = objv.get_default_phase(i, res_name)
        return arrow_scale, default_arrow_scale
