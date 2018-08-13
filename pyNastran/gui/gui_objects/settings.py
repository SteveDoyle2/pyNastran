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
from __future__ import print_function
from six import itervalues, iteritems # PY3
import numpy as np
from qtpy import QtGui

from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.gui_objects.coord_properties import CoordProperties
from pyNastran.gui.gui_objects.utils import get_setting

BLACK = (0.0, 0.0, 0.0)
WHITE = (1., 1., 1.)
GREY = (119/255., 136/255., 153/255.)


class Settings(object):
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
        self.annotation_color = BLACK

        self.text_size = 14
        self.text_color = BLACK

        self.show_info = True
        self.show_debug = True
        self.show_command = True
        self.show_warning = True
        self.show_error = True

        # int
        self.annotation_size = 18
        self.font_size = 8
        self.magnify = 5

        # floats
        self.coord_scale = 0.05  # in percent of max dimension
        self.coord_text_scale = 0.5 # percent of nominal

        # string
        self.colormap = 'jet' # 'viridis'

        # not stored
        self.dim_max = 1.0
        #self.annotation_scale = 1.0

    def reset_settings(self):
        """helper method for ``setup_gui``"""
        # rgb tuple
        self.use_gradient_background = True
        self.background_color = GREY
        self.background_color2 = GREY

        self.annotation_color = BLACK
        self.text_color = BLACK

        self.show_info = True
        self.show_debug = True
        self.show_command = True
        self.show_warning = True
        self.show_error = True

        # int
        self.text_size = 14
        self.annotation_size = 18
        self.font_size = 8
        self.magnify = 5

        # float
        self.coord_scale = 0.05
        self.coord_text_scale = 0.5

        # string
        self.colormap = 'jet' # 'viridis'

        self.parent.resize(1100, 700)

        # not stored
        self.dim_max = 1.0
        #self.annotation_scale = 1.0

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
        self._set_setting(settings, setting_keys, ['show_info'], self.show_info, True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_debug'], self.show_debug, True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_command'], self.show_command, True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_warning'], self.show_warning, True, auto_type=bool)
        self._set_setting(settings, setting_keys, ['show_error'], self.show_error, True, auto_type=bool)

        # the vtk panel background color
        self._set_setting(settings, setting_keys, ['use_gradient_background'],
                          False, auto_type=bool)
        self._set_setting(settings, setting_keys, ['background_color', 'backgroundColor'],
                          GREY, auto_type=float)
        self._set_setting(settings, setting_keys, ['background_color2'], GREY, auto_type=float)

        # scales the coordinate systems
        self._set_setting(settings, setting_keys, ['coord_scale'], self.coord_scale, auto_type=float)
        self._set_setting(settings, setting_keys, ['coord_text_scale'], self.coord_text_scale, auto_type=float)

        # this is for the 3d annotation
        self._set_setting(settings, setting_keys, ['annotation_color', 'labelColor'],
                          BLACK, auto_type=float)
        self._set_setting(settings, setting_keys, ['annotation_size'], 18, auto_type=int) # int
        if isinstance(self.annotation_size, float):
            # throw the float in the trash as it's from an old version of vtk
            self.annotation_size = 18
        elif isinstance(self.annotation_size, int):
            pass
        else:
            print('annotation_size = ', self.annotation_size)

        self._set_setting(settings, setting_keys, ['magnify'], self.magnify, auto_type=int)

        # this is the text in the lower left corner
        self._set_setting(settings, setting_keys, ['text_color', 'textColor'],
                          BLACK, auto_type=float)
        self._set_setting(settings, setting_keys, ['text_size'], 14, auto_type=int)

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

        #w = screen_shape.width()
        #h = screen_shape.height()
        #try:
        if screen_shape:
            self.parent.resize(screen_shape[0], screen_shape[1])
            #width, height = screen_shape

        font = QtGui.QFont()
        font.setPointSize(self.font_size)
        self.parent.setFont(font)

        #if 0 and PY3:
            #pos_default = 0, 0
            #pos = settings.value("pos", pos_default)
            #x_pos, y_pos = pos
            #print(pos)
            #self.mapToGlobal(QtCore.QPoint(pos[0], pos[1]))
            #y_pos = pos_default[0]
            #self.parent.setGeometry(x_pos, y_pos, width, height)
        #except TypeError:
            #self.resize(1100, 700)

    def _set_setting(self, settings, setting_keys, setting_names, default,
                     save=True, auto_type=None):
        """
        helper method for ``reapply_settings``
        """
        set_name = setting_names[0]
        value = get_setting(settings, setting_keys, setting_names, default,
                            auto_type=auto_type)
        if save:
            setattr(self, set_name, value)
        return value

    def save(self, settings):
        """saves the settings"""
        settings.setValue('main_window_geometry', self.parent.saveGeometry())
        settings.setValue('mainWindowState', self.parent.saveState())

        # rgb tuple
        settings.setValue('use_gradient_background', self.use_gradient_background)
        settings.setValue('background_color', self.background_color)
        settings.setValue('background_color2', self.background_color2)
        settings.setValue('annotation_color', self.annotation_color)
        settings.setValue('text_color', self.text_color)

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
        settings.setValue('coord_scale', self.coord_scale)
        settings.setValue('coord_text_scale', self.coord_text_scale)

        # str
        settings.setValue('colormap', self.colormap)

        #screen_shape = QtGui.QDesktopWidget().screenGeometry()
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

        # case attached annotations (typical)
        for follower_actors in itervalues(self.parent.label_actors):
            for follower_actor in follower_actors:
                follower_actor.GetTextProperty().SetFontSize(size)
                follower_actor.Modified()

        # geometry property attached annotations (e.g., flaps)
        for obj in itervalues(self.parent.geometry_properties):
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

    def update_coord_scale(self, coord_scale=None, render=True):
        """internal method for updating the coordinate system size"""
        if coord_scale is None:
            coord_scale = self.coord_scale
        dim_max = self.dim_max
        scale = coord_scale * dim_max

        for unused_coord_id, axes in iteritems(self.parent.axes):
            axes.SetTotalLength(scale, scale, scale)
        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()

    def update_coord_text_scale(self, coord_text_scale=None, render=True):
        """internal method for updating the coordinate system size"""
        if coord_text_scale is None:
            coord_text_scale = self.coord_text_scale

        for unused_coord_id, axes in iteritems(self.parent.axes):
            texts = [
                axes.GetXAxisCaptionActor2D(),
                axes.GetYAxisCaptionActor2D(),
                axes.GetZAxisCaptionActor2D(),
            ]
            # this doesn't set the width
            # this being very large (old=0.1) makes the width constraint inactive
            width = 1.0
            height = 0.25
            for text in texts:
                text.SetWidth(coord_text_scale * width)
                text.SetHeight(coord_text_scale * height)

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

        # case attached annotations (typical)
        for follower_actors in itervalues(self.parent.label_actors):
            for follower_actor in follower_actors:
                prop = follower_actor.GetProperty()
                prop.SetColor(*color)

        # geometry property attached annotations (e.g., flaps)
        for obj in itervalues(self.parent.geometry_properties):
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
        for text_actor in itervalues(self.parent.text_actors):
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
        for text_actor in itervalues(self.parent.text_actors):
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
        for text_actor in itervalues(self.parent.text_actors):
            text_prop = text_actor.GetTextProperty()
            text_prop.SetFontSize(text_size)

    def set_magnify(self, magnify=5):
        """sets the screenshot magnification factor (int)"""
        self.magnify = magnify

    def __repr__(self):
        return '<Settings>'


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
