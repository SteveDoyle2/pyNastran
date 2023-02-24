class BaseGui:
    def __init__(self, gui):
        self.gui = gui

    @property
    def log(self):
        """links the the GUI's log"""
        return self.gui.log

    @property
    def settings(self):
        """gets the gui settings"""
        return self.gui.settings

