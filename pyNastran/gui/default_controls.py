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

BASE_TOOL_SHORTCUTS = {
#BASE_FILE_TOOLS: dict[str, str] = {
    # shortcut
    'exit': 'Ctrl+Q',

    #'reload': (               'Reload Model...',           'treload.png',       '',       'Remove the model and reload the same geometry file'),
    'load_geometry': 'Ctrl+O',
    'load_results': 'Ctrl+R',
    #'load_csv_user_geom': (   'Load CSV User Geometry...', '',                  '',       'Loads custom geometry file'),
    #'load_csv_user_points': ( 'Load CSV User Points...',   'user_points.png',   '',       'Loads CSV points'),
    #'load_custom_result': (   'Load Custom Results...',    '',                  '',       'Loads a custom results file'),

    #'save_vtk': ( 'Export VTK...',        '',             '', 'Export a VTK file'),
    #'script': (  'Run Python Script...', 'python48.png', '', 'Runs pyNastranGUI in batch mode'),
#}
#
#
#BASE_TOOLS: dict[str, str] = {
    # labels
    'label_clear': 'CTRL+W',
    #'label_reset': ('Clear All Labels',     '', '',       'Clear all labels'),

    # view
    'wireframe': 'w',
    'surface': 's',
    'screenshot': 'CTRL+I',

    # geometry
    # Geometry:
    #  - Create
    #  - Modify
    #'geometry': ('Geometry', 'geometry.png', '', 'Geometry'),
    #
    # core menus
    'legend':             'CTRL+L',
    'animation':          'CTRL+A',
    #'clipping':           ('Set Clipping...',            '',                '',       'Set Clipping'),
    'set_preferences':    'CTRL+P',
    'geo_properties':     'CTRL+E',
    #'map_element_fringe': 'CTRL+F',

    #('axis', 'Show/Hide Axis', 'axis.png', None, 'Show/Hide Global Axis', self.on_show_hide_axes),

    # groups
    'modify_groups': 'CTRL+G',
    #create_groups_by_visible_result': (   'Create Groups By Visible Result', '', '',       'Create Groups'),
    #'create_groups_by_property_id': (      'Create Groups By Property ID',    '', '',       'Create Groups'),
    #'create_groups_by_model_group': (      'Create Groups By Model Group',    '', '',       'Create Groups'),
    #('create_list', 'Create Lists through Booleans', '', None, 'Create List', self.create_list),

    # logging
    #'show_info': (    'Show INFO',    'show_info.png',    '', 'Show "INFO" messages'),
    #'show_debug': (   'Show DEBUG',   'show_debug.png',   '', 'Show "DEBUG" messages'),
    #'show_command': ( 'Show COMMAND', 'show_command.png', '', 'Show "COMMAND" messages'),
    #'show_warning': ( 'Show WARNING', 'show_warning.png', '', 'Show "COMMAND" messages'),
    #'show_error': (   'Show ERROR',   'show_error.png',   '', 'Show "COMMAND" messages'),

    # zoom
    'magnify': 'm',
    'shrink': 'Shift+M',

    # rotation
    'rotate_clockwise': 'o',
    'rotate_cclockwise': 'Shift+O',


    #('cell_pick', 'Cell Pick', '', 'c', 'Centroidal Picking', self.on_cell_picker),
    #('node_pick', 'Node Pick', '', 'n', 'Nodal Picking', self.on_node_picker),

    # help
    #'website': (          'Open pyNastran Website...',       '',           '',       'Open the pyNastran website'),
    #'docs': (             'Open pyNastran Docs Website...',  '',           '',       'Open the pyNastran documentation website'),
    #'report_issue': (     'Report a Bug/Feature Request...', '',           '',       'Open the pyNastran issue tracker'),
    #'discussion_forum': ( 'Discussion Forum Website...',     '',           '',       'Open the discussion forum to ask questions'),
    'about': 'CTRL+H',

    # camera
    #'view': (         'Camera View',       'view.png',     '',  'Load the camera menu'),
    'camera_reset': 'r',

    # results
    'cycle_results': 'L',
    'rcycle_results': 'k',

    # view actions
    'back_view': 'x',
    'right_view': 'y',
    'top_view': 'z',
    'front_view': 'Shift+X',
    'left_view': 'Shift+Y',
    'bottom_view': 'Shift+Z',

    'edges': 'e',
    'edges_black': 'b',
    #'anti_alias_0': ( 'Off',              '',                 '',  'Disable Anti-Aliasing'),
    #'anti_alias_1': ( '1x',               '',                 '',  'Set Anti-Aliasing to 1x'),
    #'anti_alias_2': ( '2x',               '',                 '',  'Set Anti-Aliasing to 2x'),
    #'anti_alias_4': ( '4x',               '',                 '',  'Set Anti-Aliasing to 4x'),
    #'anti_alias_8': ( '8x',               '',                 '',  'Set Anti-Aliasing to 8x'),

    # mouse buttons
    'rotation_center': 'f',
    #'measure_distance': ( 'Measure Distance',        'measure_distance.png', '',  'Measure the distance between two nodes'),
    #'highlight_cell': (   'Highlight Cell',          '',                     '',  'Highlight a single cell'),
    #'highlight_node': (   'Highlight Node',          '',                     '',  'Highlight a single node'),

    # name, gui_name, png, shortcut, desc, func
    #'probe_result': (       'Probe the model and mark it with the value of a node/element', 'tprobe.png', '',  'Probe the displayed result'),
    'quick_probe_result': 'p',

    #Probe all the results at the given location (slow!)
   # 'probe_result_all': (       'Probe All Results', 'tprobe_all.png', '',  'Probe results for all cases'),
    'quick_probe_result_all': 'a',

    #'zoom': ( 'Zoom', 'zoom.png', '', 'Zoom In'),

    # font size
    'font_size_increase': 'Ctrl+Plus',
    'font_size_decrease': 'Ctrl+Minus',

    # picking
    #'area_pick': (                'Area Pick', 'tarea_pick.png', '', 'Get a list of nodes/elements'),
    #'highlight': (                'Highlight', 'thighlight.png', '', 'Highlight a list of nodes/elements'),
    #'highlight_nodes_elements': ( 'Highlight', 'thighlight.png', '', 'Highlight a list of nodes/elements'),
    #'mark_nodes_elements': (      'Mark',      'tmark.png',      '', 'Mark a list of nodes/elements'),
}
