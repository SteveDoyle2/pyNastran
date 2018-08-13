# plugins should be defined in reverse order
# the last plugin will get the final callback
#
# the path for auto_wireframe is intended to not be correct (so it doesn't load),
# but is a demonstration of how to use a plugin
#
#

plugin_name_to_path = [
    ('auto_wireframe', r'C:\pyNastran\pyNastran\gui\plugins\auto_wireframe.py', 'AutoWireframe'),
    #'aero_panels'
]

# consider adding order control based on the function being called
