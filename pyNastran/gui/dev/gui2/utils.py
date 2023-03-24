from __future__ import annotations
import os
from typing import Tuple, List, Dict, Set, Callable, Any, TYPE_CHECKING

from qtpy import QtGui
from qtpy.QtWidgets import QAction, QToolBar, QMenu

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from qtpy.QtWidgets import QMainWindow

def build_actions(self: QMainWindow,
                  base_actions: dict[str, QAction],
                  icon_path: str,
                  tools_list: list[tuple[str, str, str, str, str, Callable[Any]]],
                  checkables_set: Set[str],
                  log: SimpleLogger) -> dict[str, Any]:
    checkables = {}

    actions = base_actions
    for tool in tools_list:
        (name, txt, icon, shortcut, tip, func) = tool
        if name in actions:
            log.error('trying to create a duplicate action %r' % name)
            continue

        if icon is None:
            msg = f'missing_icon = {name!r}!!!'
            print(msg)
            self.log.warning(msg)
            ico = None
        else:
            ico = QtGui.QIcon()
            pth = os.path.join(icon_path, icon)
            if icon.endswith('.png') and not os.path.exists(pth):
                self.log.warning(str((name, pth)))

            ico.addPixmap(QtGui.QPixmap(pth), QtGui.QIcon.Normal, QtGui.QIcon.Off)

        if name in checkables:
            is_checked = checkables[name]
            actions[name] = QAction(ico, txt, self, checkable=True)
            actions[name].setChecked(is_checked)
        else:
            actions[name] = QAction(ico, txt, self)

        if shortcut:
            actions[name].setShortcut(shortcut)
            #actions[name].setShortcutContext(QtCore.Qt.WidgetShortcut)
        if tip:
            actions[name].setStatusTip(tip)
        if func:
            actions[name].triggered.connect(func)
    return actions


def fill_menus(self: QMainWindow,
               menus_list: list[tuple[str, str, list[str]]],
               actions: dict[str, QAction],
               allow_missing_actions: bool=False) -> None:
    assert len(self.actions) > 0, self.actions
    menus = {}
    for name, header, actions_to_add in menus_list:
        for i, item in enumerate(menus_list):
            if item != '':
                break
        actions_to_add = actions_to_add[i:]
        if len(actions_to_add) == 0:
            continue

        #file_menu = self.menubar.addMenu('&Help')
        if isinstance(header, str):
            menu = self.menubar.addMenu(header)
            assert isinstance(menu, QMenu), menu
        elif isinstance(header, QToolBar):
            menu = header
        else:
            raise TypeError(header)

        for action_name in actions_to_add:
            if action_name == '':
                menu.addSeparator()
                continue
            try:
                action = self.actions[action_name]
            except KeyError:
                if not allow_missing_actions:
                    raise
                self.log.warning(f'missing action {action_name}')
            menu.addAction(action)
        menus[name] = menu
    return menus
