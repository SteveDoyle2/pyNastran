from __future__ import print_function
from six import itervalues, PY3
import numpy as np
import vtk
from qtpy import QtGui

from pyNastran.gui.qt_files.alt_geometry_storage import AltGeometry
from pyNastran.gui.qt_files.coord_properties import CoordProperties

BLACK = (0.0, 0.0, 0.0)
WHITE = (1., 1., 1.)
GREY = (119/255., 136/255., 153/255.)
USE_ANNOTATION_INT = int(vtk.VTK_VERSION[0]) >= 7

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
        self.background_color = GREY
        self.annotation_color = BLACK
        self.text_color = BLACK

        # int
        self.annotation_size_int = 18
        self.font_size = 8

        # not stored
        self.dim_max = 1.0
        self.annotation_scale = 1.0

    def reset_settings(self):
        """helper method for ``setup_gui``"""
        # rgb tuple
        self.background_color = GREY
        self.annotation_color = BLACK
        self.text_color = BLACK

        # int
        self.annotation_size_int = 18
        self.font_size = 8

        self.parent.resize(1100, 700)

        # not stored
        self.dim_max = 1.0
        self.annotation_scale = 1.0

    def load(self, settings):
        """helper method for ``setup_gui``"""
        #red = (1.0, 0.0, 0.0)
        screen_shape_default = (1100, 700)

        setting_keys = [str(key) for key in settings.childKeys()]

        # sets the window size/position
        main_window_geometry = self._get_setting(
            settings, setting_keys, ['main_window_geometry', 'mainWindowGeometry'], None)
        if main_window_geometry is not None:
            self.parent.restoreGeometry(main_window_geometry)

        # this is the gui font
        self._set_setting(settings, setting_keys, ['font_size'], self.font_size)

        # the vtk panel background color
        self._set_setting(settings, setting_keys, ['background_color', 'backgroundColor'], GREY)

        # this is for the 3d annotation
        self._set_setting(settings, setting_keys, ['annotation_color', 'labelColor'], BLACK)
        self._set_setting(settings, setting_keys, ['annotation_size_int'], 18)

        # this is the text in the lower left corner
        self._set_setting(settings, setting_keys, ['text_color', 'textColor'], BLACK)
        screen_shape = self._set_setting(settings, setting_keys, ['screen_shape'],
                                         screen_shape_default, save=False)

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
        self.parent.resize(screen_shape[0], screen_shape[1])
        width, height = screen_shape

        font = QtGui.QFont()
        font.setPointSize(self.font_size)
        #self.app.setFont(font)
        self.parent.setFont(font)

        if 0 and PY3:
            pos_default = 0, 0
            pos = settings.value("pos", pos_default)
            x_pos, y_pos = pos
            #print(pos)
            #self.mapToGlobal(QtCore.QPoint(pos[0], pos[1]))
            y_pos = pos_default[0]
            self.parent.setGeometry(x_pos, y_pos, width, height)
        #except TypeError:
            #self.resize(1100, 700)

    def _get_setting(self, settings, setting_keys, setting_names, default):
        """
        helper method for ``reapply_settings``

        does this, but for a variable number of input names, but one output name:
            screen_shape = settings.value("screen_shape", screen_shape_default)

        If the registry key is not defined, the default is used.
        """
        set_name = setting_names[0]
        pull_name = None
        for key in setting_names:
            if key in setting_keys:
                pull_name = key
                break
        if pull_name is None:
            value = default
        else:
            value = settings.value(pull_name, default)
        return value

    def _set_setting(self, settings, setting_keys, setting_names, default, save=True):
        """
        helper method for ``reapply_settings``
        """
        set_name = setting_names[0]
        value = self._get_setting(settings, setting_keys, setting_names, default)
        if save:
            setattr(self.parent, set_name, value)
        return value

    def save(self, settings):
        """saves the settings"""
        settings.setValue('main_window_geometry', self.parent.saveGeometry())
        settings.setValue('mainWindowState', self.parent.saveState())

        # rgb tuple
        settings.setValue('background_color', self.background_color)
        settings.setValue('annotation_color', self.annotation_color)
        settings.setValue('text_color', self.text_color)

        # int
        settings.setValue('font_size', self.font_size)
        settings.setValue('annotation_size_int', self.annotation_size_int)

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
    def on_increase_text_size(self):
        self.on_set_font_size(self.font_size + 1)

    def on_decrease_text_size(self):
        self.on_set_font_size(self.font_size - 1)

    def on_set_font_size(self, font_size, show_command=True):
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
            assert isinstance(size, (int, float)), 'size=%r' % size
            self.set_annotation_size(size)
        if color is not None:
            assert len(color) == 3, color
            assert isinstance(color[0], float), 'color=%r' % color
            self.set_annotation_color(color)

    @property
    def annotation_text_size(self):
        return self.dim_max * 0.02 * self.annotation_scale

    @annotation_text_size.setter
    def annotation_text_size(self, annotation_text_size):
        # annotation_text_size = dim_max * 0.02 * annotation_scale
        #a = b * c * d
        #d = a / bc
        self.annotation_scale = annotation_text_size / (self.dim_max * 0.02)

    def set_annotation_size(self, size, render=True):
        """Updates the size of all the annotations"""
        assert size >= 0., size
        self.annotation_text_size = size

        if USE_ANNOTATION_INT:
            size = int(size)
            assert size > 0, size

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

        else:
            # case attached annotations (typical)
            for follower_actors in itervalues(self.parent.label_actors):
                for follower_actor in follower_actors:
                    follower_actor.SetScale(size)
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
                    follower_actor.SetScale(size)
                    follower_actor.Modified()

        if render:
            self.parent.vtk_interactor.GetRenderWindow().Render()
            self.parent.log_command('settings.set_annotation_size(%s)' % size)

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
    def set_background_color_to_white(self):
        """sets the background color to white; used by gif writing?"""
        self.set_background_color(WHITE)

    def set_background_color(self, color):
        """
        Set the background color

        Parameters
        ----------
        color : (float, float, float)
            RGB values as floats
        """
        self.background_color = color
        self.parent.rend.SetBackground(*color)
        self.parent.vtk_interactor.Render()
        self.parent.log_command('settings.set_background_color(%s, %s, %s)' % color)

    #---------------------------------------------------------------------------
    # TEXT ACTORS - used for lower left notes

    def set_text_color(self, color):
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
        self.parent.vtk_interactor.Render()
        self.parent.log_command('settings.set_text_color(%s, %s, %s)' % color)

    def update_text_size(self, magnify=1.0):
        """Internal method for updating the bottom-left text when we go to take a picture"""
        text_size = int(14 * magnify)
        for text_actor in itervalues(self.parent.text_actors):
            text_prop = text_actor.GetTextProperty()
            text_prop.SetFontSize(text_size)

def repr_settings(settings):
    """works on a QSettings, not a Settings"""
    msg = 'QSettings:\n'
    for key in sorted(settings.allKeys()):
        value = settings.value(key)
        msg += '    %r : %r\n' % (key, value)
    return msg
