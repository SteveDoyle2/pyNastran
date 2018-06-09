"""
Defines TestUGrid
"""
from __future__ import print_function
import os
import unittest

import pyNastran
from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid
from pyNastran.converters.aflr.ugrid.ugrid3d_to_nastran import ugrid3d_to_nastran
from pyNastran.converters.aflr.ugrid.ugrid3d_to_tecplot import (
    ugrid_to_tecplot, ugrid3d_to_tecplot_filename, read_ugrid)
from pyNastran.utils.log import get_logger

PKG_PATH = pyNastran.__path__[0]
UGRID_PATH = os.path.join(PKG_PATH, 'converters', 'aflr', 'ugrid', 'models')
TECPLOT_PATH = os.path.join(PKG_PATH, 'converters', 'tecplot', 'models')
NASTRAN_PATH = os.path.join(PKG_PATH, '..', 'models')

class TestUgrid(unittest.TestCase):
    """runs ugrid2d/3d tests"""

    def test_ugrid_01(self):
        """tests solid_bending.bdf"""
        nastran_filename1 = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending.bdf')
        ugrid_filename = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending.b8.ugrid')
        log = get_logger(level='warning')

        nastran_to_ugrid(nastran_filename1, ugrid_filename_out=ugrid_filename,
                         properties=None, check_shells=False, check_solids=True, log=log)
        assert os.path.exists(ugrid_filename), ugrid_filename

        nastran_filename2 = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending2.bdf')
        model = ugrid3d_to_nastran(ugrid_filename, nastran_filename2,
                                   include_shells=True, include_solids=True,
                                   convert_pyram_to_penta=False,
                                   encoding=None, size=16,
                                   is_double=False, log=log)
        model.skin_solids()
        assert os.path.exists(nastran_filename2), nastran_filename2

        #tecplot_filename1 = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending.plt')
        #ugrid3d_to_tecplot_filename(model, tecplot_filename1)
        #assert os.path.exists(tecplot_filename1), tecplot_filename1

        tecplot_filename2 = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending2.plt')
        tecplot = ugrid_to_tecplot(model)
        tecplot.write_tecplot(tecplot_filename2, res_types=None,
                              is_points=True,
                              adjust_nids=True)
        assert os.path.exists(tecplot_filename2), tecplot_filename2

    def test_ugrid3d_gui_box(self):
        """simple UGRID3D box model"""
        ugrid_filename = os.path.join(UGRID_PATH, 'box.b8.ugrid')
        log = get_logger(level='warning')
        tecplot_filename2 = os.path.join(UGRID_PATH, 'box.plt')

        ugrid_model = read_ugrid(ugrid_filename, log=log)
        tecplot = ugrid3d_to_tecplot_filename(ugrid_filename, tecplot_filename2)
        tecplot = ugrid_to_tecplot(ugrid_filename)
        tecplot = ugrid_to_tecplot(ugrid_model)
        tecplot.write_tecplot(tecplot_filename2, res_types=None,
                              is_points=True,
                              adjust_nids=True)
        assert os.path.exists(tecplot_filename2), tecplot_filename2


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
