import os
import unittest
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.tecplot.tecplot_io import TecplotIO
from pyNastran.converters.aflr.ugrid.ugrid3d_to_tecplot import (
    ugrid3d_to_tecplot_filename, read_ugrid)


PKG_PATH = pyNastran.__path__[0]
TECPLOT_PATH = os.path.join(PKG_PATH, 'converters', 'tecplot')
UGRID_PATH = os.path.join(PKG_PATH, 'converters', 'aflr', 'ugrid', 'models')


class TecplotGUI(TecplotIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = TecplotIO(self)
        self.build_fmts(['tecplot'], stop_on_failure=True)

class TestTecplotGUI(unittest.TestCase):

    def test_tecplot_geometry(self):
        test = TecplotGUI()
        tecplot_filename = os.path.join(TECPLOT_PATH, 'models', 'ascii', 'point_fetri_2d_02.dat')
        #test.model.load_tecplot_geometry(tecplot_filename, '')
        test.on_load_geometry(tecplot_filename, geometry_format='tecplot', raise_error=True)
        #test.model.load_nastran_geometry(geometry_filename)
        #test.model.load_shabp_geometry(geometry_filename)

    #def test_tecplot_results(self):
        #pass
        #test.model.load_panair_geometry(geometry_filename)
        #test.model.load_panair_results(agps_filename)

    def test_tecplot_box(self):
        """simple UGRID3D box model"""
        ugrid_filename = os.path.join(UGRID_PATH, 'box.b8.ugrid')
        #log = get_logger(level='warning')
        log = None
        tecplot_filename2 = os.path.join(TECPLOT_PATH, 'box.plt')

        ugrid_model = read_ugrid(ugrid_filename, log=log)
        tecplot = ugrid3d_to_tecplot_filename(ugrid_filename, tecplot_filename2)
        tecplot.write_tecplot(tecplot_filename2, res_types=None,
                              is_points=True,
                              adjust_nids=True)

        test = TecplotGUI()
        tecplot_filename = os.path.join(TECPLOT_PATH, 'models', 'ascii', 'point_fetri_2d_02.dat')
        test.on_load_geometry(tecplot_filename, geometry_format='tecplot', raise_error=True)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

