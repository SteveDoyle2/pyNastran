from typing import Callable


class Action:
    def __init__(self, name: str, text: str, icon: str='',
                 shortcut: str='',
                 tip: str='',
                 func=Callable, show: bool=True):
        self.name = name
        self.text = text
        self.ico = icon
        self.shortcut = shortcut
        self.tip = tip
        self.func = func
        self.show = show

    def __repr__(self) -> str:
        return f'Action(name={self.name}, text={self.text}'

    @property
    def icon_path(self) -> str:
        return self.ico
