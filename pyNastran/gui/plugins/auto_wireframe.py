from pyNastran.gui.plugins.module import Module

class AutoWireframe(Module):
    def __init__(self, gui):
        Module.__init__(self, gui)
        self._create_menu_items = self.gui._create_menu_items
        self.gui.create_menu_items = self._create_menu_items

    #def get_tools_checkables(self):
        #return [], {}

    def post_load_geometry(self):
        self.gui.on_wireframe()

    def _create_menu_items(self, actions=None, create_menu_bar=True):
        menu_items = self._create_menu_items()
        #_file, _view, _window, _help, _scripts, _toolbar, _hidden = menu_items
        #menu_items[0]
        #menu_items
        #    (self.menu_file, menu_file),
        #    (self.menu_view, menu_view),
        #    (self.menu_window, menu_window),
        #    (self.menu_help, ('about',)),
        #    (self.menu_scripts, scripts),
        #    (self.toolbar, toolbar_tools),
        #    (self.menu_hidden, hidden_tools),
        #    # (self.menu_scripts, ()),
        #    #(self._dummy_toolbar, ('cell_pick', 'node_pick'))
        #]
