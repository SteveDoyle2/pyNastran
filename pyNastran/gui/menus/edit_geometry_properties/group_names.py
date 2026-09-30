"""
Generic name-grouping machinery for the ``Edit Geometry Properties`` menu.

The GUI stores every actor in a flat ``dict[str, AltGeometry | CoordProperties]``
keyed by name, so a model with 40 coordinate systems produces 40 rows that must
each be toggled by hand.  This module collapses those flat names into groups
so the menu can show a tree.

This file is **format agnostic** and deliberately contains no knowledge of
Nastran/Cart3d/UGRID naming.  There are two ways a name gets a group, checked
in this order:

1. an explicit ``group`` attribute on the properties object, e.g.::

       AltGeometry(self, 'BAR_y', group='Bar Axes', representation='bar')

2. a list of ``(group_name, regex)`` rules supplied by the format's IO class

The attribute always wins, so the IO layer can be made precise over time
without anyone having to touch the regexes.  Set ``group=''`` on an object to
explicitly force it to stay ungrouped.

A format contributes its rules by defining a
``get_<format>_geometry_group_rules()`` method on its IO class, following the
same optional-hook convention as ``_create_<format>_tools_and_menu_items()``.
See ``pyNastran/converters/nastran/gui/geometry_groups.py`` for the reference
implementation.  Formats that define no hook simply get an ungrouped flat
list, exactly as before.

"""
from __future__ import annotations
import re
from typing import Any, Callable, Optional

#: name of the optional hook each format's IO class may define
GROUP_RULES_FUNC = 'get_%s_geometry_group_rules'

#: (group_name, regex) pairs; the first matching rule wins, so order matters.
#: Empty by design - rules belong to the format-specific converters.
GROUP_RULES: list[tuple[str, str]] = []

#: names that are bookkeeping flags in ``out_data``, not actors
NON_ACTOR_KEYS = frozenset({'clicked_ok', 'clicked_cancel', 'font_size'})

#: 'main' is the primary mesh; it should always sort first and never be grouped
PRIMARY_NAME = 'main'

#: compiled-regex cache, keyed by the rule tuple; ``group_of`` is called once
#: per actor name, so this avoids recompiling the whole table each time
_COMPILED_CACHE: dict[tuple, list[tuple[str, re.Pattern]]] = {}


def _compile(rules: Optional[list[tuple[str, str]]]) -> list[tuple[str, re.Pattern]]:
    if rules is None:
        rules = GROUP_RULES
    if not rules:
        return []
    key = tuple(rules)
    compiled = _COMPILED_CACHE.get(key)
    if compiled is None:
        compiled = [(group_name, re.compile(pattern))
                    for group_name, pattern in rules]
        _COMPILED_CACHE[key] = compiled
    return compiled


def group_of(name: str,
             obj: Any=None,
             rules: Optional[list[tuple[str, str]]]=None) -> Optional[str]:
    """
    Determines the group a single actor belongs to.

    Parameters
    ----------
    name : str
        the actor name (the ``geometry_properties`` key)
    obj : AltGeometry | CoordProperties | None
        the properties object; if it has a ``group`` attribute, that wins
    rules : list[(str, str)] | None
        overrides ``GROUP_RULES``

    Returns
    -------
    group : str | None
        the group name, or None if the actor is ungrouped

    """
    if name == PRIMARY_NAME:
        return None

    # an explicit attribute always beats the regexes; group='' forces ungrouped
    group = getattr(obj, 'group', None)
    if group is not None:
        group = str(group).strip()
        return group if group else None

    for group_name, regex in _compile(rules):
        if regex.search(name):
            return group_name
    return None


def get_group_rules_for_format(gui) -> list[tuple[str, str]]:
    """
    Asks the active format's IO class for its grouping rules.

    Mirrors the optional-hook convention used by
    ``_create_<format>_tools_and_menu_items`` in ``gui_common.py``: the
    method is looked up by name and skipped if absent, so a format that
    doesn't care about grouping needs no changes at all.

    Nastran mixes its IO into MainWindow, while other formats live in
    ``gui.load_actions.model_objs``, so both are checked.

    Parameters
    ----------
    gui : MainWindow | None
        the gui; ``gui.format`` names the active format

    Returns
    -------
    rules : list[(str, str)]
        ``(group_name, regex)`` pairs; empty if the format has no hook

    """
    fmt = getattr(gui, 'format', None)
    if not fmt:
        return []

    func_name = GROUP_RULES_FUNC % fmt

    # nastran is a mixin on MainWindow; everything else is a model_obj
    owners = [gui]
    model_objs = getattr(getattr(gui, 'load_actions', None), 'model_objs', None)
    if model_objs:
        owners.extend(model_objs.values())

    for owner in owners:
        if hasattr(owner, func_name):
            try:
                return list(getattr(owner, func_name)())
            except Exception as error:  # pragma: no cover
                log = getattr(gui, 'log', None)
                if log is not None:
                    log.error(f'{func_name} failed: {error}')
                return []
    return []


def natural_key(name: str) -> list:
    """sorts 'Coord 2' before 'Coord 10' instead of after it"""
    return [int(token) if token.isdigit() else token.lower()
            for token in re.split(r'(\d+)', name)]


def group_geometry_names(
        data: dict[str, Any],
        group_func: Optional[Callable[[str, Any], Optional[str]]]=None,
        rules: Optional[list[tuple[str, str]]]=None,
        min_group_size: int=2) -> list[tuple[Optional[str], list[str]]]:
    """
    Buckets the flat ``geometry_properties`` dict into display groups.

    Parameters
    ----------
    data : dict[str, AltGeometry | CoordProperties]
        the ``out_data`` dict; bookkeeping keys are ignored
    group_func : callable(name, obj) -> str | None; default=None
        overrides the default attribute-then-regex lookup entirely
    rules : list[(str, str)] | None
        overrides ``GROUP_RULES`` (ignored if ``group_func`` is given)
    min_group_size : int; default=2
        groups with fewer members than this are dissolved back into
        standalone rows; a group of one saves nobody any clicks

    Returns
    -------
    groups : list[(str | None, list[str])]
        ordered ``(group_name, member_names)`` pairs.  A group name of None
        means the members are standalone top-level rows.  ``'main'`` is always
        first.

    """
    if group_func is None:
        def group_func(name, obj):
            return group_of(name, obj, rules=rules)

    # group -> members, plus the order each group was first seen
    members: dict[Optional[str], list[str]] = {}
    order: list[Optional[str]] = []
    standalone: list[str] = []

    for name in data:
        if name in NON_ACTOR_KEYS:
            continue
        group = group_func(name, data[name])
        if group is None:
            # keep each ungrouped name as its own slot so ordering is stable
            standalone.append(name)
            if name not in order:
                order.append(name)
            continue
        if group not in members:
            members[group] = []
            order.append(group)
        members[group].append(name)

    standalone_set = set(standalone)
    groups: list[tuple[Optional[str], list[str]]] = []
    for key in order:
        if key in standalone_set and key not in members:
            groups.append((None, [key]))
            continue
        names = sorted(members[key], key=natural_key)
        if len(names) < min_group_size:
            # dissolve: show the lone member as a plain row
            groups.extend((None, [name]) for name in names)
        else:
            groups.append((key, names))

    # the main mesh always goes on top
    groups.sort(key=lambda pair: pair[1] != [PRIMARY_NAME])
    return groups
