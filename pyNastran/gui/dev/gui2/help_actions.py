from __future__ import annotations
import urllib.request
import urllib
import urllib.error
import webbrowser
from typing import TYPE_CHECKING

import pyNastran
from pyNastran.gui.menus.about.about import AboutWindow

if TYPE_CHECKING:  # pragma: no cover
    from qtpy.QtWidgets import QMainWindow


class HelpActions:
    def __init__(self, gui: QMainWindow):
        self.gui = gui

        self.tools_list = [
            # help
            ('website', 'Open pyNastran Website...', '', None, 'Open the pyNastran website', self.open_website),
            ('docs', 'Open pyNastran Docs Website...', '', None, 'Open the pyNastran documentation website', self.open_docs),
            ('report_issue', 'Report a Bug/Feature Request...', '', None, 'Open the pyNastran issue tracker', self.open_issue),
            ('discussion_forum', 'Discussion Forum Website...', '', None, 'Open the discussion forum to ask questions', self.open_discussion_forum),
            ('about', 'About pyNastran GUI...', 'tabout.png', 'CTRL+H', 'About pyNastran GUI and help on shortcuts', self.about_dialog),
        ]
        self.actions_list = ['website', 'docs', 'report_issue', 'discussion_forum', 'about',]

    # help
    def open_website(self) -> None:
        """opens the website"""
        self._openbrowswer(pyNastran.__website__)

    def open_docs(self) -> None:
        url = pyNastran.__docs__
        try:
            urllib.request.urlopen(url)
        except (urllib.error.HTTPError, urllib.error.URLError):
            return
        self._openbrowswer(url)

    def open_issue(self) -> None:
        """loads the pyNastran issue tracker"""
        self._openbrowswer(pyNastran.__issue__)

    def open_discussion_forum(self) -> None:
        """loads the pyNastran discussion forum website"""
        self._openbrowswer(pyNastran.__discussion_forum__)

    def about_dialog(self) -> None:
        """Display about dialog"""
        data = {
            'font_size': self.gui.settings.font_size,
            #'font_size': self.gui.settings.font_size,
        }
        win = AboutWindow(data, win_parent=self.gui, show_tol=True)
        win.show()

    def _openbrowswer(self, url: str) -> None:
        """opens a URL"""
        if self.gui.is_gui:
            webbrowser.open(url)
