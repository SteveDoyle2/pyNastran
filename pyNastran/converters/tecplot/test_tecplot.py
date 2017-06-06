import os

from pyNastran.converters.tecplot.tecplot import read_tecplot
from pyNastran.converters.tecplot.tecplot_to_nastran import tecplot_to_nastran_filename
from pyNastran.converters.nastran.nastran_to_tecplot import nastran_to_tecplot, nastran_to_tecplot_filename
from pyNastran.utils.log import get_logger
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tecplot', 'models')
nastran_path = os.path.join(pkg_path, '..', 'models')

import unittest

class TestTecplot(unittest.TestCase):

    def test_tecplot_01(self):
        log = get_logger(level='warning')
        tecplot_filename1 = os.path.join(model_path, 'ascii', 'point_fetri_2d_02.dat')
        tecplot_filename2 = os.path.join(model_path, 'ascii', 'point_fetri_2d_02.dat_out')

        tecplot = read_tecplot(tecplot_filename1, log=log)
        #tecplot.write_tecplot(tecplot_filename2, res_types=None,
                              #is_points=True, adjust_nids=True)
        #os.remove(tecplot_filename2)

    def test_tecplot_02(self):
        log = get_logger(level='warning')
        nastran_filename1 = os.path.join(nastran_path, 'solid_bending', 'solid_bending.bdf')
        nastran_filename2 = os.path.join(nastran_path, 'solid_bending', 'solid_bending2.bdf')
        tecplot_filename = os.path.join(nastran_path, 'solid_bending', 'solid_bending.plt')
        tecplot = nastran_to_tecplot_filename(nastran_filename1, tecplot_filename, log=log)
        #tecplot.write_tecplot(tecplot_filename)
        #tecplot_to_nastran_filename(tecplot_filename, nastran_filename2)
        #os.remove(nastran_filename2)
        #os.remove(tecplot_filename)

    def _test_tecplot_02(self):
        nastran_filename1 = os.path.join(nastran_path, 'solid_bending', 'solid_bending.bdf')
        nastran_filename2 = os.path.join(nastran_path, 'solid_bending', 'solid_bending2.bdf')
        tecplot_filename = os.path.join(nastran_path, 'solid_bending', 'solid_bending.plt')
        tecplot = nastran_to_tecplot_filename(nastran_filename1, tecplot_filename, log=log)
        tecplot_to_nastran_filename(tecplot_filename, nastran_filename2, log=log)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

