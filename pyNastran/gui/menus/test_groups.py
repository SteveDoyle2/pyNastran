import os
import unittest
from qtpy import QtGui


from pyNastran.gui.menus.test_menu import UsesQApplication
from pyNastran.gui.gui import MainWindow

import pyNastran
PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


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

        bdf_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes.dat')
        #op2_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes_m2.op2')
        gui.on_load_geometry_button(
            infile_name=bdf_filename, geometry_format='nastran', name='main', raise_error=False)

        # TODO: fix this...
        #with self.assertRaises(AttributeError):
        gui.clear_application_log(force=True)
        gui.on_reset_camera()

        gui.show_hide_max_actor(render=True)
        gui.show_hide_min_actor(render=True)
        gui.show_hide_max_actor(render=False)
        gui.show_hide_min_actor(render=False)

        gui._cycle_results()
        gui.closeEvent()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
