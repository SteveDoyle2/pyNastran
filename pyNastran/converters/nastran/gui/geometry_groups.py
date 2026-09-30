"""
Nastran-specific grouping rules for the ``Edit Geometry Properties`` menu.

``NastranIO`` creates one actor per coordinate system, per spline, per SPC set,
etc., which makes the menu's flat list unusable on a real model.  These rules
bucket those names into collapsible groups so a whole family can be
hidden/shown/recolored at once.

The rules are consumed by ``pyNastran/gui/menus/edit_geometry_properties/
group_names.py`` via the ``get_nastran_geometry_group_rules`` hook, following
the same optional-hook convention as ``_create_nastran_tools_and_menu_items``.

Each entry is ``(group_name, regex)`` and the FIRST match wins, so order
matters; put the specific patterns above the general ones.  The regexes are
applied with ``re.search`` against the raw actor name.

An actor can always opt out of these rules by carrying an explicit ``group``
attribute (``group=''`` forces it to stay ungrouped), which takes priority
over everything here.

"""
from __future__ import annotations

#: (group_name, regex) pairs for the names produced by NastranIO.
#: The comment on each line is the name-producing site in the converter.
NASTRAN_GROUP_RULES: list[tuple[str, str]] = [
    # tool_actions.py: 'Global XYZ' / 'Coord 1' / 'Coord 10'
    ('Coords', r'^(Coord\s+\S+|Global XYZ)$'),

    # nastran_io.py: bar_type + '_y'/'_z' -> 'BAR_y', 'TUBE_z', 'I_y', 'HAT1_z'
    ('Bar Axes', r'^[A-Z][A-Z0-9_]*_[yz]$'),

    # 'SPC=100', 'SPC=100: Subcases=1, 2', 'spc_id=100'
    ('SPCs', r'^(SPC=|spc_id=)'),

    # 'MPC=3_lines', 'mpc_id=3_dependent'
    ('MPCs', r'^(MPC=|mpc_id=)'),

    # 'SUPORT', 'suport1_id=5'
    ('SUPORTs', r'^(SUPORT$|suport1_id=)'),

    # 'rigid_dependent', 'rigid_independent', 'rigid_lines'
    ('Rigid Elements', r'^rigid_(dependent|independent|lines)$'),

    # 'caero', 'caero_boxes', 'caero_control_surfaces'
    ('CAEROs', r'^caero'),

    # nastran_io.py:893 lowercases the card name, so these are
    # 'spline1_600000_boxes' / 'spline1_600000_structure_points'.
    # Matched case-insensitively anyway, plus the 'all_spline_points' summary
    # actor that belongs with them.
    ('Splines', r'(?i)^(spline\d*_\d+_|all_spline_points$)'),

    # 'ELEV_control_surface'
    ('Control Surfaces', r'_control_surface$'),

    # 'MONPNT1: WING xyz', 'MONPNT1: WING GRIDs; cid=0'
    ('MONPNTs', r'^MONPNT\d*\s*:'),

    # 'element coord', 'material coord', 'mcid ply=3'
    ('Element Coords', r'^(element coord|material coord|mcid ply=)'),

    ('Microphones', r'^MICPNT-'),
    ('AML Regions', r'^AMLREG-'),

    # shear-moment-torque / cutting plane: 'smt_plane', 'smt_vector', ...
    ('Cutting Plane', r'^smt_'),
]


class NastranGeometryGroups:
    """mixin supplying the Edit Geometry Properties grouping rules"""

    def get_nastran_geometry_group_rules(self) -> list[tuple[str, str]]:
        """
        Groups the actors NastranIO creates (coords, splines, SPCs, ...).

        Returns
        -------
        rules : list[(str, str)]
            ``(group_name, regex)`` pairs, first match wins

        """
        return NASTRAN_GROUP_RULES
