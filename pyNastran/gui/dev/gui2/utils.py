from __future__ import annotations
import os
from typing import TYPE_CHECKING

from qtpy import QtGui
from qtpy.QtWidgets import QAction

if TYPE_CHECKING:
    from cpylog import SimpleLogger
    from qtpy.QtWidgets import QMainWindow

def build_actions(self: QMainWindow,
                  icon_path: str,
                  tools_list: List[Tuple[str, str, str, str, str, Any]],
                  checkables_set: Set[str],
                  log: SimpleLogger) -> Dict[str, Any]:
    checkables = {}

    actions = {}
    for tool in tools_list:
        (name, txt, icon, shortcut, tip, func) = tool
        if name in actions:
            log.error('trying to create a duplicate action %r' % name)
            continue

        if icon is None:
            print(f'missing_icon = {name!r}!!!')
            ico = None
        else:
            ico = QtGui.QIcon()
            pth = os.path.join(icon_path, icon)
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
               menus_list: List[Tuple[str, str, List[str]]],
               actions: Dict[str, QAction],
               ) -> None:
    assert len(self.actions) > 0, self.actions
    menus = {}
    for name, header, actions_to_add in menus_list:
        #file_menu = self.menubar.addMenu('&Help')
        menu = self.menubar.addMenu(header)
        for action_name in actions_to_add:
            if action_name == '':
                menu.addSeparator()
                continue
            action = self.actions[action_name]
            menu.addAction(action)
        menus[name] = menu
    return menus
