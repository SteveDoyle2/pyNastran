from pyNastran.gui.qt_files.base_gui import BaseGui
from pyNastran.gui.menus.highlight.highlight import HighlightWindow
from pyNastran.gui.gui_objects.settings import Settings


class HighlightObject(BaseGui):
    def __init__(self, gui):
        #self.gui = gui
        super().__init__(gui)
        self._highlight_window_shown = False
        self._highlight_window = None

    def set_font_size(self, font_size):
        """sets the font size for the preferences window"""
        if self._highlight_window_shown:
            self._highlight_window.set_font_size(font_size)

    def set_menu(self):
        """
        Opens a dialog box to set:

        +--------+----------+
        |  Max   |  Float   |
        +--------+----------+
        """
        #if not hasattr(self, 'case_keys'):  # TODO: maybe include...
            #self.log_error('No model has been loaded.')
            #return

        settings = self.gui.settings
        if not hasattr(self.gui, 'case_keys'):
            self.gui.log_error('No model has been loaded.')
            return

        model_name = 'main'
        data = {
            'font_size' : settings.font_size,
            'model_name' : model_name,
            'clicked_ok' : False,
            'close' : False,
        }
        gui = self.gui
        nodes = gui.get_node_ids(model_name=model_name)
        elements = gui.get_element_ids(model_name=model_name)
        if nodes is None or elements is None:
            gui.log.error('No model was found')
            return

        if not self._highlight_window_shown:
            self._highlight_window = HighlightWindow(data, win_parent=self.gui, menu_type='highlight')
            self._highlight_window.show()
            self._highlight_window_shown = True
            self._highlight_window.exec_()
        else:
            self._highlight_window.activateWindow()

        if data['close']:
            if not self._highlight_window._updated_window:
                settings.on_set_font_size(data['font_size'])
            del self._highlight_window
            self._highlight_window_shown = False
        else:
            self._highlight_window.activateWindow()


class MarkObject(BaseGui):
    def __init__(self, gui):
        #self.gui = gui
        super().__init__(gui)
        self._mark_window_shown = False
        self._mark_window = None

    def set_font_size(self, font_size: int) -> None:
        """sets the font size for the preferences window"""
        if self._mark_window_shown:
            self._mark_window.set_font_size(font_size)

    def set_menu(self):
        """
        Opens a dialog box to set:

        +--------+----------+
        |  Max   |  Float   |
        +--------+----------+
        """
        #if not hasattr(self, 'case_keys'):  # TODO: maybe include...
            #self.log_error('No model has been loaded.')
            #return

        settings: Settings = self.gui.settings
        if not hasattr(self.gui, 'case_keys'):
            self.gui.log_error('No model has been loaded.')
            return

        model_name = 'main'
        data = {
            'font_size' : settings.font_size,
            'model_name' : model_name,
            'clicked_ok' : False,
            'close' : False,
        }
        gui = self.gui
        nodes = gui.get_node_ids(model_name=model_name)
        elements = gui.get_element_ids(model_name=model_name)
        if nodes is None or elements is None:
            gui.log.error('No model was found')
            return

        if not self._mark_window_shown:
            self._mark_window = HighlightWindow(data, win_parent=self.gui, menu_type='mark')
            self._mark_window.show()
            self._mark_window_shown = True
            self._mark_window.exec_()
        else:
            self._mark_window.activateWindow()

        if data['close']:
            if not self._mark_window._updated_window:
                settings.on_set_font_size(data['font_size'])
            del self._mark_window
            self._mark_window_shown = False
        else:
            self._mark_window.activateWindow()
