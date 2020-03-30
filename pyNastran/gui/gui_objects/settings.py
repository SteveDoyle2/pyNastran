"""
defines:
 - Settings(parent)
   - reset_settings(self)
   - load(self, settings)
   - save(self, settings)
   - on_increase_text_size(self)
   - on_decrease_font_size(self)
   - on_set_font_size(self, font_size, show_command=True)
   - set_annotation_size_color(self, size=None, color=None)
   - set_annotation_size(self, size, render=True)
   - set_annotation_color(self, color, render=True)
   - set_background_color_to_white(self)
   - set_background_color(self, color)
   - set_text_color(self, color)
   - update_text_size(self, magnify=1.0)

 - repr_settings(settings)

"""
from __future__ import annotations
from typing import List, Dict, Any, TYPE_CHECKING
import numpy as np
from qtpy import QtGui

from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.gui_objects.coord_properties import CoordProperties
from pyNastran.gui.gui_objects.utils import get_setting
from pyNastran.utils import object_attributes
if TYPE_CHECKING:  # pragma: no cover
    import vtk


BLACK = (0.0, 0.0, 0.0)
WHITE = (1., 1., 1.)
GREY = (119/255., 136/255., 153/255.)
ORANGE = (229/255., 92/255., 0.)
HIGHLIGHT_OPACITY = 0.9
HIGHLIGHT_POINT_SIZE = 10.
HIGHLIGHT_LINE_THICKNESS = 5.
ANNOTATION_SIZE = 18
FONT_SIZE = 8
TEXT_SIZE = 14
COORD_SCALE = 0.05  # in percent of max dimension
COORD_TEXT_SCALE = 0.5 # percent of nominal

NASTRAN_BOOL_KEYS = [
    'nastran_create_coords',
    'nastran_is_properties',
    'nastran_is_element_quality',
    'nastran_is_bar_axes',
    'nastran_is_3d_bars', 'nastran_is_3d_bars_update',
    'nastran_is_shell_mcids', 'nastran_is_update_conm2',

    'nastran_stress', 'nastran_plate_stress', 'nastran_composite_plate_stress',
    'nastran_strain', 'nastran_plate_strain', 'nastran_composite_plate_strain',
    'nastran_rod_stress', 'nastran_bar_stress', 'nastran_beam_stress',
    'nastran_rod_strain', 'nastran_bar_strain', 'nastran_beam_strain',
    'nastran_spring_stress', 'nastran_solid_stress',
    'nastran_spring_strain', 'nastran_solid_strain',

    'nastran_force',
    'nastran_bar_force', 'nastran_beam_force', 'nastran_plate_force',
    'nastran_spring_force', 'nastran_gap_force', 'nastran_cbush_force',
]

class Settings:
    """storage class for various settings"""
    def __init__(self, parent):
        """
        Creates the Settings object

        Parameters
        ----------
        parent : MainWindow()
            used by the class to access the MainWindow
        """
        self.parent = parent

        # rgb tuple
        self.use_gradient_background = True
        self.background_color = GREY
        self.background_color2 = GREY

        # TODO: what is an annotation color?
        self.annotation_color = BLACK

        # text in the lower left corner
        self.text_size = TEXT_SIZE
        self.text_color = BLACK

        # used for highlight actors
        self.highlight_color = ORANGE
        self.highlight_opacity = HIGHLIGHT_OPACITY
        self.highlight_point_size = HIGHLIGHT_POINT_SIZE
        self.highlight_line_thickness = HIGHLIGHT_LINE_THICKNESS

        self.show_info = True
        self.show_debug = True
        self.show_command = True
        self.show_warning = True
        self.show_error = True

        # int
        self.annotation_size = ANNOTATION_SIZE
        self.font_size = FONT_SIZE
        self.magnify = 5

        # floats
        self.coord_scale = COORD_SCALE
        self.coord_text_scale = COORD_TEXT_SCALE
        self.coord_linewidth = 2.0

        # string
        self.colormap = 'jet' # 'viridis'

        # not stored
        self.dim_max = 1.0
        #self.annotation_scale = 1.0

        self.nastran_is_element_quality = True
        self.nastran_is_properties = True
        self.nastran_is_3d_bars = True
        self.nastran_is_3d_bars_update = True
        self.nastran_create_coords = True
        self.nastran_is_bar_axes = True
        self.nastran_is_shell_mcids = True
        self.nastran_is_update_conm2 = True

        self.nastran_stress = True
        self.nastran_spring_stress = True
        self.nastran_rod_stress = True
        self.nastran_bar_stress = True
        self.nastran_beam_stress = True
        self.nastran_plate_stress = True
        self.nastran_composite_plate_stress = True
        self.nastran_solid_stress = True

        self.nastran_strain = True
        self.nastran_spring_strain = True
        self.nastran_rod_strain = True
        self.nastran_bar_strain = True
        self.nastran_beam_strain = True
        self.nastran_plate_strain = True
        self.nastran_composite_plate_strain = True
        self.nastran_solid_strain = True

        self.nastran_force = True
        self.nastran_spring_force = True
        self.nastran_cbush_force = True
        self.nastran_gap_force = True
        self.nastran_bar_force = True
        self.nastran_beam_force = True
        self.nastran_plate_force = True


    def reset_settings(self):
        """helper method for ``setup_gui``"""
        # rgb tuple
        self.use_gradient_background = True
        self.background_color = GREY
        self.background_color2 = GREY

        self.annotation_size = ANNOTATION_SIZE
        self.annotation_color = BLACK

        self.text_size = TEXT_SIZE
        self.text_color = BLACK

        self.highlight_color = ORANGE
        self.highlight_opacity = HIGHLIGHT_OPACITY
        self.highlight_point_size = HIGHLIGHT_POINT_SIZE
        self.highlight_line_thickness = HIGHLIGHT_LINE_THICKNESS

        self.show_info = True
        self.show_debug = True
        self.show_command = True
        self.show_warning = True
        self.show_error = True

        # int
        self.font_size = FONT_SIZE
        self.magnify = 5

        # float
        self.coord_scale = COORD_SCALE
        self.coord_text_scale = COORD_TEXT_SCALE
        self.coord_linewidth = 2.0

        # string
        self.colormap = 'jet' # 'viridis'

        self.parent.resize(1100, 700)

        # not stored
        self.dim_max = 1.0
        #self.annotation_scale = 1.0

        self.nastran_is_element_quality = True
        self.nastran_is_properties = True
        self.nastran_is_3d_bars = True
        self.nastran_is_3d_bars_update = True
        self.nastran_create_coords = True
        self.nastran_is_bar_axes = True
        self.nastran_is_shell_mcids = True
        self.nastran_is_update_conm2 = True

        self.nastran_stress = True
        self.nastran_spring_stress = True
        self.nastran_rod_stress = True
        self.nastran_bar_stress = True
        self.nastran_beam_stress = True
        self.nastran_plate_stress = True
        self.nastran_composite_plate_stress = True
        self.nastran_solid_stress = True

        self.nastran_strain = True
        self.nastran_spring_strain = True
        self.nastran_rod_strain = True
        self.nastran_bar_strain = True
        self.nastran_beam_strain = True
        self.nastran_plate_strain = True
        self.nastran_composite_plate_strain = True
        self.nastran_solid_strain = True

        self.nastran_force = True
        self.nastran_spring_force = True
        self.nastran_cbush_force = True
        self.nastran_gap_force = True
        self.nastran_bar_force = True
        self.nastran_beam_force = True
        self.nastran_plate_force = True

    def load(self, settings):
        """helper method for ``setup_gui``"""
        #red = (1.0, 0.0, 0.0)
        screen_shape_default = (1100, 700)

        setting_keys = [str(key) for key in settings.childKeys()]

        # sets the window size/position
        main_window_geometry = get_setting(
            settings, setting_keys, ['main_window_geometry', 'mainWindowGeometry'], None)
        if main_window_geometry is not None:
            self.parent.restoreGeometry(main_window_geometry)

        # this is the gui font
        self._set_setting(settings, setting_keys, ['font_size'], self.font_size, auto_type=int)

        # the info/debug/gui/command preferences
        self._set_setting(settings, setting_keys, ['show_info'], self.show_info,
                          True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_debug'], self.show_debug,
                          True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_command'], self.show_command,
                          True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_warning'], self.show_warning,
                          True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_error'], self.show_error,
                          True, auto_type=bool)

        # the vtk panel background color
        self._set_setting(settings, setting_keys, ['use_gradient_background'],
                          False, auto_type=bool)
        self._set_setting(settings, setting_keys, ['background_color', 'backgroundColor'],
                          GREY, auto_type=float)
        self._set_setting(settings, setting_keys, ['background_color2'], GREY, auto_type=float)

        # scales the coordinate systems
        self._set_setting(settings, setting_keys, ['coord_scale'], COORD_SCALE, auto_type=float)
        self._set_setting(settings, setting_keys, ['coord_text_scale'], COORD_TEXT_SCALE, auto_type=float)

        # this is for the 3d annotation
        self._set_setting(settings, setting_keys, ['annotation_color', 'labelColor'],
                          BLACK, auto_type=float)
        self._set_setting(settings, setting_keys, ['annotation_size'], ANNOTATION_SIZE, auto_type=int) # int
        if isinstance(self.annotation_size, float):
            # throw the float in the trash as it's from an old version of vtk
            self.annotation_size = ANNOTATION_SIZE
        elif isinstance(self.annotation_size, int):
            pass
        else:
            print('annotation_size = ', self.annotation_size)

        self._set_setting(settings, setting_keys, ['magnify'], self.magnify, auto_type=int)

        # this is the text in the lower left corner
        self._set_setting(settings, setting_keys, ['text_color', 'textColor'],
                          BLACK, auto_type=float)
        self._set_setting(settings, setting_keys, ['text_size'], TEXT_SIZE, auto_type=int)

        # highlight
        self._set_setting(settings, setting_keys, ['highlight_color'],
                          ORANGE, auto_type=float)
        self._set_setting(settings, setting_keys, ['highlight_opacity'],
                          HIGHLIGHT_OPACITY, auto_type=float)
        self._set_setting(settings, setting_keys, ['highlight_point_size'],
                          HIGHLIGHT_POINT_SIZE, auto_type=float)
        self._set_setting(settings, setting_keys, ['highlight_line_thickness'],
                          HIGHLIGHT_LINE_THICKNESS, auto_type=float)
        #self._set_setting(settings, setting_keys, ['highlight_style'],
                          #HIGHLIGHT_OPACITY, auto_type=float)

        # default colormap for legend
        self._set_setting(settings, setting_keys, ['colormap'],
                          'jet')


        # general gui sizing
        screen_shape = self._set_setting(settings, setting_keys, ['screen_shape'],
                                         screen_shape_default, save=False, auto_type=int)

        #try:
            #screen_shape = settings.value("screen_shape", screen_shape_default)
        #except (TypeError, AttributeError):
            #screen_shape = screen_shape_default

        #if 'recent_files' in setting_keys:
        try:
            self.parent.recent_files = settings.value("recent_files", self.recent_files)
        except (TypeError, AttributeError):
            pass

        for key in NASTRAN_BOOL_KEYS:
            default = getattr(self, key)
            self._set_setting(settings, setting_keys, [key],
                              default, save=True, auto_type=bool)

        #w = screen_shape.width()
        #h = screen_shape.height()
        #try:
        if screen_shape:
            self.parent.resize(screen_shape[0], screen_shape[1])
            #width, height = screen_shape

        font = QtGui.QFont()
        font.setPointSize(self.font_size)
        self.parent.setFont(font)

        #if 0:
            #pos_default = 0, 0
            #pos = settings.value("pos", pos_default)
            #x_pos, y_pos = pos
            #print(pos)
            #self.mapToGlobal(QtCore.QPoint(pos[0], pos[1]))
            #y_pos = pos_default[0]
            #self.parent.setGeometry(x_pos, y_pos, width, height)
        #except TypeError:
            #self.resize(1100, 700)
        is_loaded = True
        return is_loaded

    def _set_setting(self, settings, setting_keys: List[str],
                     setting_names: List[str], default: Any,
                     save: bool=True, auto_type=None) -> Any:
        """
        helper method for ``reapply_settings``
        """
        set_name = setting_names[0]
        value = get_setting(settings, setting_keys, setting_names, default,
                            auto_type=auto_type)
        if save:
            setattr(self, set_name, value)
        return value

    def save(self, settings, is_testing: bool=False) -> None:
        """saves the settings"""
        if not is_testing:
            settings.setValue('main_window_geometry', self.parent.saveGeometry())
            settings.setValue('mainWindowState', self.parent.saveState())

        # rgb tuple
        settings.setValue('use_gradient_background', self.use_gradient_background)
        settings.setValue('background_color', self.background_color)
        settings.setValue('background_color2', self.background_color2)
        settings.setValue('annotation_color', self.annotation_color)
        settings.setValue('text_color', self.text_color)

        settings.setValue('highlight_color', self.highlight_color)
        settings.setValue('highlight_opacity', self.highlight_opacity)

        settings.setValue('show_info', self.show_info)
        settings.setValue('show_debug', self.show_debug)
        settings.setValue('show_command', self.show_command)
        settings.setValue('show_warning', self.show_warning)
        settings.setValue('show_error', self.show_error)

        # int
        settings.setValue('font_size', self.font_size)
        settings.setValue('annotation_size', self.annotation_size)
        settings.setValue('magnify', self.magnify)

        # float
        settings.setValue('text_size', self.text_size)
        settings.setValue('coord_scale', self.coord_scale)
        settings.setValue('coord_text_scale', self.coord_text_scale)

        # str
        settings.setValue('colormap', self.colormap)

        # format-specific
        for key in NASTRAN_BOOL_KEYS:
            value = getattr(self, key)
            settings.setValue(key, value)


        #screen_shape = QtGui.QDesktopWidget().screenGeometry()
        if not is_testing:
            main_window = self.parent.window()
            width = main_window.frameGeometry().width()
            height = main_window.frameGeometry().height()
            settings.setValue('screen_shape', (width, height))

            qpos = self.parent.pos()
            pos = qpos.x(), qpos.y()
            settings.setValue('pos', pos)

    #---------------------------------------------------------------------------
    # FONT SIZE
    def on_increase_font_size(self):
        """shrinks the overall GUI font size"""
        self.on_set_font_size(self.font_size + 1)

    def on_decrease_font_size(self):
        """shrinks the overall GUI font size"""
        self.on_set_font_size(self.font_size - 1)

    def on_set_font_size(self, font_size, show_command=True):
        """updates the GUI font size"""
        return self.parent.on_set_font_size(font_size, show_command=show_command)

    #---------------------------------------------------------------------------
    # ANNOTATION SIZE/COLOR
    def set_annotation_size_color(self, size=None, color=None):
        """
        Parameters
        ----------
        size : float
            annotation size
        color : (float, float, float)
            RGB values
        """
        if size is not None:
            assert isinstance(size, int), 'size=%r' % size
            self.set_annotation_size(size)
        if color is not None:
            assert len(color) == 3, color
            assert isinstance(color[0], float), 'color=%r' % color
            self.set_annotation_color(color)

    def set_annotation_size(self, size, render=True):
        """Updates the size of all the annotations"""
        assert size >= 0, size
        assert isinstance(size, int), size
        if self.annotation_size == size:
            return
        self.annotation_size = size

        # min/max
        for actor in self.parent.min_max_actors:
            actor.GetTextProperty().SetFontSize(size)
            actor.Modified()

        # case attached annotations (typical)
        for follower_actors in self.parent.label_actors.values():
            for follower_actor in follower_actors:
                follower_actor.GetTextProperty().SetFontSize(size)
                follower_actor.Modified()

        # geometry property attached annotations (e.g., flaps)
        for obj in self.parent.geometry_properties.values():
            if isinstance(obj, CoordProperties):
                continue
            elif isinstance(obj, AltGeometry):
                pass
            else:
                raise NotImplementedError(obj)

            follower_actors = obj.label_actors
            for follower_actor in follower_actors:
                follower_actor.GetTextProperty().SetFontSize(size)
                follower_actor.Modified()

        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()
            self.parent.log_command('settings.set_annotation_size(%s)' % size)

    def set_coord_scale(self, coord_scale, render=True):
        """sets the coordinate system size"""
        self.coord_scale = coord_scale
        self.update_coord_scale(coord_scale, render=render)

    def set_coord_text_scale(self, coord_text_scale, render=True):
        """sets the coordinate system text size"""
        self.coord_text_scale = coord_text_scale
        self.update_coord_text_scale(coord_text_scale, render=render)

    def update_coord_scale(self, coord_scale=None, coord_text_scale=None,
                           linewidth=None, render=True):
        """internal method for updating the coordinate system size"""
        if coord_scale is None:
            coord_scale = self.coord_scale
        #if coord_text_scale:
            #self.update_coord_text_scale(coord_text_scale=coord_text_scale, render=False)

        dim_max = self.dim_max
        scale = coord_scale * dim_max

        for unused_coord_id, axes in self.parent.axes.items():
            axes.SetTotalLength(scale, scale, scale) # was coord_scale
            #axes.SetScale(magnify, magnify, magnify)
            #if linewidth:
                #xaxis = axes.GetXAxisShaftProperty()
                #yaxis = axes.GetXAxisShaftProperty()
                #zaxis = axes.GetXAxisShaftProperty()
                #lw = xaxis.GetLineWidth()  #  1.0
                #xaxis.SetLineWidth(linewidth)
                #yaxis.SetLineWidth(linewidth)
                #zaxis.SetLineWidth(linewidth)
            #print(f'coord_scale coord_id={unused_coord_id} scale={scale} lw={linewidth}')

        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()

    def scale_coord(self, magnify: float, render=True):
        """internal method for scaling the coordinate system size"""
        for unused_coord_id, axes in self.parent.axes.items():
            axes.SetScale(magnify)
        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()

    def update_coord_text_scale(self, coord_text_scale=None, render=True):
        """internal method for updating the coordinate system size"""
        if coord_text_scale is None:
            coord_text_scale = self.coord_text_scale

        update_axes_text_size(self.parent.axes, coord_text_scale,
                              width=1.0, height=0.25)
        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()

    def set_annotation_color(self, color, render=True):
        """
        Set the annotation color

        Parameters
        ----------
        color : (float, float, float)
            RGB values as floats
        """
        if np.allclose(self.annotation_color, color):
            return
        self.annotation_color = color

        # min/max
        for min_max_actor in self.parent.min_max_actors:
            #print(dir(min_max_actor))
            prop = min_max_actor.GetProperty()
            prop.SetColor(*color)

        # case attached annotations (typical)
        for follower_actors in self.parent.label_actors.values():
            for follower_actor in follower_actors:
                prop = follower_actor.GetProperty()
                prop.SetColor(*color)

        # geometry property attached annotations (e.g., flaps)
        for obj in self.parent.geometry_properties.values():
            if isinstance(obj, CoordProperties):
                continue
            elif isinstance(obj, AltGeometry):
                pass
            else:
                raise NotImplementedError(obj)

            follower_actors = obj.label_actors
            for follower_actor in follower_actors:
                prop = follower_actor.GetProperty()
                prop.SetColor(*color)

        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()
            self.parent.log_command('settings.set_annotation_color(%s, %s, %s)' % color)

    #---------------------------------------------------------------------------
    def set_background_color_to_white(self, render=True):
        """sets the background color to white; used by gif writing?"""
        self.set_gradient_background(use_gradient_background=False, render=False)
        self.set_background_color(WHITE, render=render)

    def set_gradient_background(self, use_gradient_background=False, render=True):
        """enables/diables the gradient background"""
        self.use_gradient_background = use_gradient_background
        self.parent.rend.SetGradientBackground(self.use_gradient_background)
        if render:
            self.parent.vtk_interactor.Render()

    def set_background_color(self, color, render=True):
        """
        Set the background color

        Parameters
        ----------
        color : (float, float, float)
            RGB values as floats
        """
        self.background_color = color
        self.parent.rend.SetBackground(*color)
        if render:
            self.parent.vtk_interactor.Render()
        self.parent.log_command('settings.set_background_color(%s, %s, %s)' % color)

    def set_background_color2(self, color, render=True):
        """
        Set the background color

        Parameters
        ----------
        color : (float, float, float)
            RGB values as floats
        """
        self.background_color2 = color
        self.parent.rend.SetBackground2(*color)
        if render:
            self.parent.vtk_interactor.Render()
        self.parent.log_command('settings.set_background_color2(%s, %s, %s)' % color)

    def set_highlight_color(self, color):
        """
        Set the highlight color

        Parameters
        ----------
        color : (float, float, float)
            RGB values as floats
        """
        self.highlight_color = color
        self.parent.log_command('settings.set_highlight_color(%s, %s, %s)' % color)

    def set_highlight_opacity(self, opacity):
        """
        Set the highlight opacity

        Parameters
        ----------
        opacity : float
            0.0 : invisible
            1.0 : solid
        """
        self.highlight_opacity = opacity
        self.parent.log_command('settings.set_highlight_opacity(%s)' % opacity)

    #---------------------------------------------------------------------------
    # TEXT ACTORS - used for lower left notes

    def set_text_color(self, color, render=True):
        """
        Set the text color

        Parameters
        ----------
        color : (float, float, float)
            RGB values as floats
        """
        self.text_color = color
        for text_actor in self.parent.text_actors.values():
            text_actor.GetTextProperty().SetColor(color)
        if render:
            self.parent.vtk_interactor.Render()
        self.parent.log_command('settings.set_text_color(%s, %s, %s)' % color)

    def set_text_size(self, text_size, render=True):
        """
        Set the text color

        Parameters
        ----------
        text_size : int
            the lower left text size (typical 14)
        """
        i = 0
        dtext_size = text_size + 1
        self.text_size = text_size
        for text_actor in self.parent.text_actors.values():
            text_prop = text_actor.GetTextProperty()
            text_prop.SetFontSize(text_size)

            position = [5, 5 + i * dtext_size]
            text_actor.SetDisplayPosition(*position)
            i += 1
        if render:
            self.parent.vtk_interactor.Render()
        self.parent.log_command('settings.set_text_size(%s)' % text_size)

    def update_text_size(self, magnify=1.0):
        """Internal method for updating the bottom-left text when we go to take a picture"""
        text_size = int(14 * magnify)
        for text_actor in self.parent.text_actors.values():
            text_prop = text_actor.GetTextProperty()
            text_prop.SetFontSize(text_size)

    def set_magnify(self, magnify=5):
        """sets the screenshot magnification factor (int)"""
        self.magnify = magnify

    def __repr__(self):
        msg = '<Settings>\n'
        for key in object_attributes(self, mode='public', keys_to_skip=['parent']):
            value = getattr(self, key)
            if isinstance(value, tuple):
                value = str(value)
            msg += '  %r = %r\n' % (key, value)
        return msg

def update_axes_text_size(axes: Dict[int, vtk.vtkAxes],
                          coord_text_scale: float,
                          width: float=1.0, height: float=0.25):
    """updates the coordinate system text size"""
    # width doesn't set the width
    # it being very large (old=0.1) makes the width constraint inactive

    for unused_coord_id, axis in axes.items():
        #print(f'coord_text_scale coord_id={unused_coord_id} coord_text_scale={coord_text_scale}')
        texts = [
            axis.GetXAxisCaptionActor2D(),
            axis.GetYAxisCaptionActor2D(),
            axis.GetZAxisCaptionActor2D(),
        ]
        # this doesn't set the width
        # this being very large (old=0.1) makes the width constraint inactive
        for text in texts:
            text.SetWidth(coord_text_scale * width)
            text.SetHeight(coord_text_scale * height)

def isfloat(value):
    """is the value floatable"""
    try:
        float(value)
        return True
    except ValueError:
        return False

def repr_settings(settings):
    """works on a QSettings, not a Settings"""
    msg = 'QSettings:\n'
    for key in sorted(settings.allKeys()):
        value = settings.value(key)
        msg += '    %r : %r\n' % (key, value)
    return msg
