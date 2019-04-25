# plugins should be defined in reverse order
# the last plugin will get the final callback
#
# the path for auto_wireframe is intended to not be correct (so it doesn't load),
# but is a demonstration of how to use a plugin
#
#
import os
#PKG_PATH = pyNastran.__file__
PLUGIN_DIR = os.path.dirname(__file__)
plugin_name_to_path = [
    #('auto_wireframe', os.path.join(PLUGIN_DIR, 'auto_wireframe.py'), 'AutoWireframe'),
    #('rfs_viewer', os.path.join(PLUGIN_DIR, 'rfs', 'rfs_viewer.py'), 'RFSViewer'),
    #'aero_panels'
]

# consider adding order control based on the function being called
