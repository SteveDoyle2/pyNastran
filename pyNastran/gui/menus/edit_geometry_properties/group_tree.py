"""
Defines the grouped tree widget used by the ``Edit Geometry Properties`` menu.

``GroupTreeView`` is a drop-in replacement for the flat ``SingleChoiceQTableView``
that adds:

 - collapsible group rows ('Coords', 'Bar Axes', ...) built by ``group_names.py``
 - extended selection, so ctrl/shift-click picks an arbitrary set of actors

Both mechanisms feed the same output: ``selected_names()`` returns the flat list
of actor names the edits should apply to.  Selecting a group row implicitly
selects all of its members.

"""
from __future__ import annotations
from typing import Any, Callable, Optional

from qtpy.QtCore import Qt
from qtpy import QtCore, QtGui
from qtpy.QtWidgets import QTreeView, QAbstractItemView

from pyNastran.gui.menus.edit_geometry_properties.group_names import (
    group_geometry_names, PRIMARY_NAME)

# stashed on each QStandardItem so we can recover the actor name(s) from a row
NAME_ROLE = Qt.UserRole + 1
IS_GROUP_ROLE = Qt.UserRole + 2


class GroupTreeView(QTreeView):
    """a QTreeView of actor names, bucketed into collapsible groups"""

    def __init__(self, parent, data: dict[str, Any],
                 group_func: Optional[Callable]=None,
                 group_rules: Optional[list]=None,
                 name: Optional[str]=None):
        QTreeView.__init__(self, parent)
        self.parent2 = parent
        self.name = name
        self.group_func = group_func
        self.group_rules = group_rules

        self.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setUniformRowHeights(True)
        self.setAlternatingRowColors(True)
        self.setExpandsOnDoubleClick(True)
        self.setEditTriggers(QAbstractItemView.NoEditTriggers)

        self.set_data(data)

    # ------------------------------------------------------------------
    # model building
    # ------------------------------------------------------------------
    def set_data(self, data: dict[str, Any]) -> None:
        """(re)builds the tree from a geometry_properties-like dict"""
        model = QtGui.QStandardItemModel(self)
        model.setHorizontalHeaderLabels(['Groups'])
        root = model.invisibleRootItem()

        groups = group_geometry_names(data, group_func=self.group_func,
                                      rules=self.group_rules)
        for group_name, names in groups:
            if group_name is None:
                for namei in names:
                    root.appendRow(_make_item(namei, namei, is_group=False))
                continue

            label = f'{group_name}  ({len(names)})'
            group_item = _make_item(label, names, is_group=True)
            font = group_item.font()
            font.setBold(True)
            group_item.setFont(font)
            for namei in names:
                group_item.appendRow(_make_item(namei, namei, is_group=False))
            root.appendRow(group_item)

        self.setModel(model)
        self.expandAll()
        header = self.header()
        header.setStretchLastSection(True)

        selection_model = self.selectionModel()
        if selection_model is not None:
            selection_model.selectionChanged.connect(self.on_selection_changed)

    # ------------------------------------------------------------------
    # selection
    # ------------------------------------------------------------------
    def selected_names(self) -> list[str]:
        """
        the flat list of actor names the edits apply to

        Selecting a group row expands to all of its members.  Order follows
        the tree, and duplicates (group + one of its children both selected)
        are collapsed.

        """
        model = self.model()
        names: list[str] = []
        seen = set()
        for index in self.selectedIndexes():
            if index.column() != 0:
                continue
            item = model.itemFromIndex(index)
            if item is None:
                continue
            payload = item.data(NAME_ROLE)
            namesi = payload if isinstance(payload, list) else [payload]
            for namei in namesi:
                if namei is not None and namei not in seen:
                    seen.add(namei)
                    names.append(namei)
        return names

    def is_group_selected(self) -> bool:
        """True if any selected row is a group header"""
        model = self.model()
        for index in self.selectedIndexes():
            if index.column() != 0:
                continue
            item = model.itemFromIndex(index)
            if item is not None and item.data(IS_GROUP_ROLE):
                return True
        return False

    def select_name(self, name: str) -> None:
        """selects the row for a single actor name"""
        model = self.model()
        matches = model.match(
            model.index(0, 0), NAME_ROLE, name, 1,
            Qt.MatchExactly | Qt.MatchRecursive)
        if not matches:
            return
        index = matches[0]
        self.setCurrentIndex(index)
        self.selectionModel().select(
            index,
            QtCore.QItemSelectionModel.ClearAndSelect | QtCore.QItemSelectionModel.Rows)
        self.scrollTo(index)

    def on_selection_changed(self, selected, deselected) -> None:
        names = self.selected_names()
        if not names:
            return
        self.parent2.update_active_names(names, is_group=self.is_group_selected())

    # ------------------------------------------------------------------
    # events
    # ------------------------------------------------------------------
    def keyPressEvent(self, event) -> None:
        key = event.key()
        if key == Qt.Key_Delete:
            names = self.selected_names()
            # never delete the main mesh
            names = [namei for namei in names if namei != PRIMARY_NAME]
            if names:
                self.parent2.on_delete_names(names)
            return
        QTreeView.keyPressEvent(self, event)


def _make_item(label: str, payload, is_group: bool) -> QtGui.QStandardItem:
    item = QtGui.QStandardItem(str(label))
    item.setEditable(False)
    item.setData(payload, NAME_ROLE)
    item.setData(is_group, IS_GROUP_ROLE)
    return item
