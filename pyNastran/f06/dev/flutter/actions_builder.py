from pathlib import Path
from functools import partial
from typing import Callable
from qtpy.QtWidgets import QAction


class Action:
    def __init__(self, name: str, text: str, icon: str='',
                 shortcut: str='',
                 func=Callable, show: bool=True):
        self.name = name
        self.text = text
        self.ico = icon
        self.shortcut = shortcut
        self.func = func
        self.show = show

    def __repr__(self) -> str:
        return f'Action(name={self.name}, text={self.text}'

    @property
    def icon_path(self) -> str:
        return self.ico


class Actions:
    def __init__(self, dirname: str | Path, actions: dict[str, Action]):
        self.dirname = str(dirname)
        if self.dirname == '.':
            self.dirname = ''
        assert self.dirname == '', self.dirname
        self.actions = actions

    def build_qactions(self, win_parent) -> dict[str, QAction]:
        qactions = {}
        tip = ''
        dirname = self.dirname
        load_icon = (dirname != '')
        for name, action_inputi in self.actions.items():
            assert isinstance(action_inputi, Action), action_inputi
            func = action_inputi.func
            icon_path = action_inputi.icon_path
            txt = action_inputi.text
            show = action_inputi.show
            shortcut = action_inputi.shortcut

            if icon_path and load_icon:
                ico = QIcon(icon_path)
                #print('icon_path = ', icon_path)
                assert os.path.exists(icon_path), icon_path
                action = QAction(txt, win_parent)
                action.setIcon(ico)
            else:
                action = QAction(txt, win_parent)

            if not show:
                action.setVisible(show)
            if shortcut:
                action.setShortcut(shortcut)
            #action.setText(name)
            #action.setIcon()
            #action.setCheckable()
            #action.setShortcut(QKeySequence(name))
            if tip:
                #print(f'  found tip for {name}')
                action.setStatusTip(tip)
            if func:
                #print(f'  found func for {name}')
                action.triggered.connect(func)
            qactions[name] = action
        return qactions

    def build_recent_file_qactions(self, win_parent,
                                   recent_filenames: list[str],
                                   set_func) -> list[str]:
        recent_files = []
        nfiles = len(recent_filenames)
        for ifile in range(win_parent.nrecent_files_max):
            name = f'file_{ifile}'
            pth = name
            func = partial(set_func, ifile)
            show = (ifile <= nfiles)
            self.actions[name] = Action(
                name=name, text=pth, func=func, show=show)
            recent_files.append(name)

        if len(recent_files):
            recent_files.append('') # dash line
        return recent_files
