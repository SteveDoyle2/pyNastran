class Module:
    def __init__(self, gui):
        self.gui = gui

    def get_tools_checkables(self):
        """
        Gets the module's tools and checkables
        Returns
        -------
        tools : List[tool]
            tool : List[name, menu_name path, command, desc, func)
            ('load_geometry', 'Load &Geometry...', 'load_geometry.png',
             'Ctrl+O', 'Loads a geometry input file', self.on_load_geometry)
        checkables : Dict[name]=bool
            is the action checkable

    Examples
    --------
        tools = [
            ('load_bdf', 'Load BDF...', 'load_bdf.png', 'Ctrl+B', 'Loads a file', self.on_load_bdf),
        ]
        ]
        checkables = {
            'load_bdf' : False,
        }
        """
        return [], {}
