"""requires qtpy"""
import os
import unittest
from qtpy import QtGui


from pyNastran.gui.menus.test.test_gui_menu import UsesQApplication
from pyNastran.gui.gui import MainWindow

import pyNastran
PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')
PLUGIN_DIR = os.path.dirname(__file__)


class TestGUI(UsesQApplication):

    #def setUp(self):
        #If you override setup, tearDown, make sure
        #to have a super call
        #super(MainWindow, self).setUp()

    #def tearDown(self):
        #super(MainWindow, self).tearDown()

    def test_real_gui(self):
        inputs = {
            'debug' : True,
            'is_groups' : True,
            'is_gui' : False,
            #'log' : 'debug',
            'log' : None,
            'geomscript' : None,
            'postscript' : None,
            'format' : None,
            'user_points' : None,
            'user_geom' : None,
        }
        gui = MainWindow(inputs)
        #gui.html_logging = False
        gui.is_gui = False
        gui.open_docs()
        gui.open_issue()
        gui.open_discussion_forum()

        # default plugins
        gui._load_plugins()

        plugin_name_to_path_bad = [
            ('bad_module1', os.path.join(PLUGIN_DIR, 'spike_module.py'), 'SpikeModule_Bad'),
            ('bad_module2', 'spike_module_doesnt_exist.py', 'SpikeModule'),
            #('rfs_viewer', os.path.join(PLUGIN_DIR, 'rfs', 'rfs_viewer.py'), 'RFSViewer'),
        ]

        plugin_name_to_path_good = [
            ('spike_module', os.path.join(PLUGIN_DIR, 'spike_module.py'), 'SpikeModule'),
        ]
        # load bad plugins and a good plugin - #3
        gui._load_plugins(plugin_name_to_path_bad)
        gui._load_plugins(plugin_name_to_path_good)

        bdf_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes.dat')
        #op2_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes_m2.op2')
        gui.on_load_geometry_button(
            infile_name=bdf_filename, geometry_format='nastran', name='main', raise_error=False)

        # TODO: fix this...
        #with self.assertRaises(AttributeError):
        gui.clear_application_log(force=True)
        gui.on_reset_camera()

        gui.view_actions.on_show_hide_max_actor()
        gui.view_actions.on_show_hide_min_actor()
        gui.view_actions.on_show_hide_max_actor()
        gui.view_actions.on_show_hide_min_actor()

        gui._cycle_results()
        gui.closeEvent()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
