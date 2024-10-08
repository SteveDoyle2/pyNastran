Tool = tuple[str, str, str, str]

CONTROLS: dict[str, tuple[str, str]] = {
    'reset': ('R', 'reset camera view'),
    'xsnap': ('Shift+X/X', 'snap to x axis'),
    'ysnap': ('Shift+Y/Y', 'snap to y axis'),
    'zsnap': ('Shift+Z/Z', 'snap to z axis'),
    # '',
    # 'h   - show/hide legend & info'),

    # shown on the menu
    # 'CTRL+I - take a screenshot (image)'),
    # 'CTRL+W - clear the labels'),
    # 'CTRL+L - Legend'),
    # 'CTRL+A - Animation'),
    # 'S      - view model as a surface'),
    # 'W      - view model as a wireframe'),

    'fwd_cycle': ('L', 'cycle the results forwards'),
    'rev_cycle': ('K', 'cycle the results backwards'),
    'scale': ('m/Shift+M', 'scale up/scale down by 1.1 times'),
    'rotate': ('o/Shift+O', 'rotate counter-clockwise/clockwise 5 degrees'),
    'pick': ('P', 'pick node/element'),
    'edge1': ('E', 'view model edges'),
    'edge2': ('B', 'change edge color from scalar/black'),
    # '',
    # 'Reload Model:  using the same filename, reload the model',
    'rotation_center': ('f', 'set rotation center (zoom out when picking to disable clipping)'),
}

BASE_TOOLS = {
#BASE_FILE_TOOLS: dict[str, Tool] = {
    # flag,  label,   picture,     shortcut, tooltip
    'exit': ('&Exit', 'texit.png', 'Ctrl+Q', 'Exit application'),

    'reload': (               'Reload Model...',           'treload.png',       '',       'Remove the model and reload the same geometry file'),
    'load_geometry': (        'Load &Geometry...',         'load_geometry.png', 'Ctrl+O', 'Loads a geometry input file'),
    'load_results': (        'Load &Results...',          'load_results.png',  'Ctrl+R', 'Loads a results file'),
    'load_csv_user_geom': (   'Load CSV User Geometry...', '',                  '',       'Loads custom geometry file'),
    'load_csv_user_points': ( 'Load CSV User Points...',   'user_points.png',   '',       'Loads CSV points'),
    'load_custom_result': (   'Load Custom Results...',    '',                  '',       'Loads a custom results file'),

    'save_vtk': ( 'Export VTK...',        '',             '', 'Export a VTK file'),
    'script': (  'Run Python Script...', 'python48.png', '', 'Runs pyNastranGUI in batch mode'),
#}
#
#
#BASE_TOOLS: dict[str, Tool] = {
    # labels
    'label_clear': ('Clear Current Labels', '', 'CTRL+W', 'Clear current labels'),
    'label_reset': ('Clear All Labels',     '', '',       'Clear all labels'),

    # view
    'wireframe': ('Wireframe Model',      'twireframe.png', 'w',      'Show Model as a Wireframe Model'),
    'surface': ('Surface Model',        'tsolid.png',     's',      'Show Model as a Surface Model'),
    'screenshot': ('Take a Screenshot...', 'tcamera.png',    'CTRL+I', 'Take a Screenshot of current view'),

    # geometry
    # Geometry:
    #  - Create
    #  - Modify
    'geometry': ('Geometry', 'geometry.png', '', 'Geometry'),
    #
    # core menus
    'legend':             ('Modify Legend...',           'legend.png',      'CTRL+L', 'Set Legend'),
    'animation':          ('Create Animation...',        'animation.png',   'CTRL+A', 'Create Animation'),
    'clipping':           ('Set Clipping...',            '',                '',       'Set Clipping'),
    'set_preferences':    ('Preferences...',             'preferences.png', 'CTRL+P', 'Set GUI Preferences'),
    'geo_properties':     ('Edit Geometry Properties...', '',               'CTRL+E', 'Change Model Color/Opacity/Line Width'),
    'map_element_fringe': ('Map Element Fringe',          '',               'CTRL+F', 'Map Elemental Centroidal Fringe Result to Nodes'),

    #('axis', 'Show/Hide Axis', 'axis.png', None, 'Show/Hide Global Axis', self.on_show_hide_axes),

    # groups
    'modify_groups': (                     'Modify Groups...',                '', 'CTRL+G', 'Create/Edit/Delete Groups'),
    'create_groups_by_visible_result': (   'Create Groups By Visible Result', '', '',       'Create Groups'),
    'create_groups_by_property_id': (      'Create Groups By Property ID',    '', '',       'Create Groups'),
    'create_groups_by_model_group': (      'Create Groups By Model Group',    '', '',       'Create Groups'),
    #('create_list', 'Create Lists through Booleans', '', None, 'Create List', self.create_list),

    # logging
    'show_info': (    'Show INFO',    'show_info.png',    '', 'Show "INFO" messages'),
    'show_debug': (   'Show DEBUG',   'show_debug.png',   '', 'Show "DEBUG" messages'),
    'show_command': ( 'Show COMMAND', 'show_command.png', '', 'Show "COMMAND" messages'),
    'show_warning': ( 'Show WARNING', 'show_warning.png', '', 'Show "COMMAND" messages'),
    'show_error': (   'Show ERROR',   'show_error.png',   '', 'Show "COMMAND" messages'),

    # zoom
    'magnify': ( 'Zoom In', 'plus_zoom.png',  'm',       'Increase Magnfication'),
    'shrink': ( 'Zoom Out', 'minus_zoom.png', 'Shift+M', 'Decrease Magnfication'),

    # rotation
    'rotate_clockwise': (  'Rotate Clockwise',         'tclock.png',  'o',       'Rotate Clockwise'),
    'rotate_cclockwise': ( 'Rotate Counter-Clockwise', 'tcclock.png', 'Shift+O', 'Rotate Counter-Clockwise'),


    #('cell_pick', 'Cell Pick', '', 'c', 'Centroidal Picking', self.on_cell_picker),
    #('node_pick', 'Node Pick', '', 'n', 'Nodal Picking', self.on_node_picker),

    # help
    'website': (          'Open pyNastran Website...',       '',           '',       'Open the pyNastran website'),
    'docs': (             'Open pyNastran Docs Website...',  '',           '',       'Open the pyNastran documentation website'),
    'report_issue': (     'Report a Bug/Feature Request...', '',           '',       'Open the pyNastran issue tracker'),
    'discussion_forum': ( 'Discussion Forum Website...',     '',           '',       'Open the discussion forum to ask questions'),
    'about': (            'About pyNastran GUI...',          'tabout.png', 'CTRL+H', 'About pyNastran GUI and help on shortcuts'),

    # camera
    'view': (         'Camera View',       'view.png',     '',  'Load the camera menu'),
    'camera_reset': ( 'Reset Camera View', 'trefresh.png', 'r', 'Reset the camera view to default'),

    # results
    'cycle_results': (  'Cycle Results', 'cycle_results.png',  'L', 'Changes the result case'),
    'rcycle_results': ( 'Cycle Results', 'rcycle_results.png', 'k', 'Changes the result case'),

    # view actions
    'back_view': (   'Back View',   'back.png',   'x',       'Flips to +X Axis'),
    'right_view': (  'Right View',  'right.png',  'y',       'Flips to +Y Axis'),
    'top_view': (    'Top View',    'top.png',    'z',       'Flips to +Z Axis'),
    'front_view': (  'Front View',  'front.png',  'Shift+X', 'Flips to -X Axis'),
    'left_view': (   'Left View',   'left.png',   'Shift+Y', 'Flips to -Y Axis'),
    'bottom_view': ( 'Bottom View', 'bottom.png', 'Shift+Z', 'Flips to -Z Axis'),


    'edges': (       'Show/Hide Edges',   'tedges.png',       'e', 'Show/Hide Model Edges'),
    'edges_black': ( 'Color Edges Black', 'tedges_color.png', 'b', 'Set Edge Color to Color/Black'),
    'anti_alias_0': ( 'Off',              '',                 '',  'Disable Anti-Aliasing'),
    'anti_alias_1': ( '1x',               '',                 '',  'Set Anti-Aliasing to 1x'),
    'anti_alias_2': ( '2x',               '',                 '',  'Set Anti-Aliasing to 2x'),
    'anti_alias_4': ( '4x',               '',                 '',  'Set Anti-Aliasing to 4x'),
    'anti_alias_8': ( '8x',               '',                 '',  'Set Anti-Aliasing to 8x'),

    # mouse buttons
    'rotation_center': (  'Set the Rotation Center', 'trotation_center.png', 'f', 'Pick a node for the rotation center/focal point'),
    'measure_distance': ( 'Measure Distance',        'measure_distance.png', '',  'Measure the distance between two nodes'),
    'highlight_cell': (   'Highlight Cell',          '',                     '',  'Highlight a single cell'),
    'highlight_node': (   'Highlight Node',          '',                     '',  'Highlight a single node'),

    # name, gui_name, png, shortcut, desc, func
    'probe_result': (       'Probe the model and mark it with the value of a node/element', 'tprobe.png', '',  'Probe the displayed result'),
    'quick_probe_result': ( 'Quick Probe',                                                  '',           'p', 'Probe the displayed result'),

    #Probe all the results at the given location (slow!)
    'probe_result_all': (       'Probe All Results', 'tprobe_all.png', '',  'Probe results for all cases'),
    'quick_probe_result_all': ( 'Quick Probe All',   '',               'a', 'Probe all cases'),

    'zoom': ( 'Zoom', 'zoom.png', '', 'Zoom In'),

    # font size
    'font_size_increase': ( 'Increase Font Size', 'text_up.png',   'Ctrl+Plus', 'Increase Font Size'),
    'font_size_decrease': ( 'Decrease Font Size', 'text_down.png', 'Ctrl+Minus', 'Decrease Font Size'),

    # picking
    'area_pick': (                'Area Pick', 'tarea_pick.png', '', 'Get a list of nodes/elements'),
    'highlight': (                'Highlight', 'thighlight.png', '', 'Highlight a list of nodes/elements'),
    'highlight_nodes_elements': ( 'Highlight', 'thighlight.png', '', 'Highlight a list of nodes/elements'),
    'mark_nodes_elements': (      'Mark',      'tmark.png',      '', 'Mark a list of nodes/elements'),
}
