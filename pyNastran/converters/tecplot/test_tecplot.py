import os

from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.tecplot.tecplot import read_tecplot
from pyNastran.converters.tecplot.tecplot import read_tecplot

import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tecplot', 'models')

import unittest

class TestTecplot(unittest.TestCase):

    def test_tecplot_geometry(self):
        tecplot_filename1 = os.path.join(model_path, 'ascii', 'point_fetri_2d_02.dat')
        tecplot_filename2 = os.path.join(model_path, 'ascii', 'point_fetri_2d_02.dat_out')

        tecplot = read_tecplot(tecplot_filename1)
        #tecplot.write_tecplot(tecplot_filename2, res_types=None,
                              #is_points=True, adjust_nids=True)
        #os.remove(tecplot_filename2)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

