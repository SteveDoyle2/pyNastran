import os

from pyNastran.gui.testing_methods import add_dummy_gui_functions
from pyNastran.converters.shabp.shabp_io import ShabpIO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'shabp')

import unittest

class ShabpGUI(unittest.TestCase):

    def test_shabp_geometry(self):
        return
        #geometry_filename = os.path.join(model_path, 'M100.inp')
        #agps_filename = os.path.join(model_path, 'agps')
        #out_filename = os.path.join(model_path, 'panair.out')
        dirname = None

        test = ShabpIO()
        test.is_nodal = True
        test.is_centroidal = True

        test = ShabpIO()
        test.is_nodal = True
        test.is_centroidal = False
        add_dummy_gui_functions(test)


        test.load_shabp_geometry('models/NAC6.INP', '')
        #test.load_nastran_geometry(geometry_filename, None)
        #test.load_shabp_geometry(geometry_filename, dirname)

    def test_shabp_results(self):
        pass
        #geometry_filename = os.path.join(model_path, 'M100.inp')
        #agps_filename = os.path.join(model_path, 'agps')
        #out_filename = os.path.join(model_path, 'panair.out')
        #dirname = None

        #test = ShabpIO()
        #test.is_nodal = True
        #test.is_centroidal = True

        #add_dummy_gui_functions(test)

        #test.load_panair_geometry(geometry_filename, dirname)
        #test.load_panair_results(agps_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

