import os

from pyNastran.gui.testing_methods import add_dummy_gui_functions
from pyNastran.converters.panair.panairIO import PanairIO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'panair','M100')

import unittest

class PanairGUI(unittest.TestCase):

    def test_m100_geometry(self):
        geometry_filename = os.path.join(model_path, 'M100.inp')
        agps_filename = os.path.join(model_path, 'agps')
        out_filename = os.path.join(model_path, 'panair.out')
        dirname = None

        test = PanairIO()
        test.is_nodal = True
        test.is_centroidal = True

        add_dummy_gui_functions(test)

        #test.load_nastran_geometry(geometry_filename, None)
        test.load_panair_geometry(geometry_filename, dirname)

    def test_m100_results(self):
        geometry_filename = os.path.join(model_path, 'M100.inp')
        agps_filename = os.path.join(model_path, 'agps')
        out_filename = os.path.join(model_path, 'panair.out')
        dirname = None

        test = PanairIO()
        test.is_nodal = True
        test.is_centroidal = True

        add_dummy_gui_functions(test)

        test.load_panair_geometry(geometry_filename, dirname)
        test.load_panair_results(agps_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

