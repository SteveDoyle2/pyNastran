import os
import unittest
from cpylog import get_logger

import pyNastran
from pyNastran.bdf.bdf import read_bdf
from pyNastran.converters.tecplot.tecplot import read_tecplot
from pyNastran.converters.tecplot.tecplot_to_nastran import tecplot_to_nastran_filename
from pyNastran.converters.nastran.nastran_to_tecplot import (
    nastran_to_tecplot, nastran_to_tecplot_filename)

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'tecplot', 'models')
NASTRAN_MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class TestTecplot(unittest.TestCase):

    def test_tecplot_01(self):
        log = get_logger(level='warning')
        tecplot_filename1 = os.path.join(MODEL_PATH, 'ascii', 'point_fetri_2d_02.dat')
        #tecplot_filename2 = os.path.join(MODEL_PATH, 'ascii', 'point_fetri_2d_02.dat_out')

        tecplot = read_tecplot(tecplot_filename1, log=log)
        #tecplot.write_tecplot(tecplot_filename2, res_types=None,
                              #is_points=True, adjust_nids=True)
        #os.remove(tecplot_filename2)

    def test_tecplot_02(self):
        log = get_logger(level='warning')
        nastran_filename1 = os.path.join(NASTRAN_MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        nastran_filename2 = os.path.join(NASTRAN_MODEL_PATH, 'solid_bending', 'solid_bending2.bdf')
        tecplot_filename = os.path.join(NASTRAN_MODEL_PATH, 'solid_bending', 'solid_bending.plt')
        tecplot = nastran_to_tecplot_filename(nastran_filename1, tecplot_filename, log=log)
        #tecplot.write_tecplot(tecplot_filename)
        tecplot_to_nastran_filename(tecplot_filename, nastran_filename2, log=log)
        #os.remove(nastran_filename2)
        #os.remove(tecplot_filename)

        bdf_model = read_bdf(nastran_filename1, log=log)
        unused_tecplot = nastran_to_tecplot(bdf_model)

    def test_tecplot_03(self):
        log = get_logger(level='warning')
        nastran_filename = os.path.join(NASTRAN_MODEL_PATH, 'elements', 'static_elements.bdf')
        tecplot_filename = os.path.join(NASTRAN_MODEL_PATH, 'elements', 'static_elements.plt')
        tecplot = nastran_to_tecplot_filename(nastran_filename, tecplot_filename, log=log)
        #tecplot2 = read_tecplot(tecplot_filename)

        bdf_model = read_bdf(nastran_filename, log=log)
        with self.assertRaises(RuntimeError):
            unused_tecplot = nastran_to_tecplot(bdf_model)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
