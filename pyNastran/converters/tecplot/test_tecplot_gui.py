import os
import unittest

from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.tecplot.tecplot_io import TecplotIO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tecplot')


class TecplotGUI(TecplotIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = TecplotIO(self)
        self.build_fmts(['tecplot'], stop_on_failure=True)

    def _remove_old_cart3d_geometry(self, tecplot_filename):
        pass

class TestTecplotGUI(unittest.TestCase):

    def test_tecplot_geometry(self):
        test = TecplotGUI()
        tecplot_filename = os.path.join(model_path, 'models', 'ascii', 'point_fetri_2d_02.dat')
        test.model.load_tecplot_geometry(tecplot_filename, '')
        #test.model.load_nastran_geometry(geometry_filename)
        #test.model.load_shabp_geometry(geometry_filename)

    def test_tecplot_results(self):
        pass
        #test.model.load_panair_geometry(geometry_filename)
        #test.model.load_panair_results(agps_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

