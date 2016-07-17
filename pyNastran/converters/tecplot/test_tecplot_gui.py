import os

from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.tecplot.tecplot_io import TecplotIO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tecplot')

import unittest

class TecplotGUI(TecplotIO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        TecplotIO.__init__(self)

    def _remove_old_cart3d_geometry(self, tecplot_filename):
        pass

class TestTecplotGUI(unittest.TestCase):

    def test_tecplot_geometry(self):
        dirname = None
        test = TecplotGUI()
        tecplot_filename = os.path.join(model_path, 'models', 'ascii', 'point_fetri_2d_02.dat')
        test.load_tecplot_geometry(tecplot_filename, '')
        #test.load_nastran_geometry(geometry_filename, None)
        #test.load_shabp_geometry(geometry_filename, dirname)

    def test_tecplot_results(self):
        pass
        #test.load_panair_geometry(geometry_filename, dirname)
        #test.load_panair_results(agps_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

