"""
Defines TestUGrid
"""
from __future__ import print_function
import os
import unittest

import pyNastran
from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid
from pyNastran.converters.aflr.ugrid.ugrid3d_to_nastran import ugrid3d_to_nastran
from pyNastran.converters.aflr.ugrid.ugrid3d_to_tecplot import ugrid_to_tecplot
from pyNastran.utils.log import get_logger

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tecplot', 'models')
nastran_path = os.path.join(pkg_path, '..', 'models')

class TestUgrid(unittest.TestCase):
    """runs ugrid2d/3d tests"""
    def test_ugrid_01(self):
        """tests solid_bending.bdf"""
        nastran_filename1 = os.path.join(nastran_path, 'solid_bending', 'solid_bending.bdf')
        ugrid_filename = os.path.join(nastran_path, 'solid_bending', 'solid_bending.b8.ugrid')
        log = get_logger(level='warning')

        nastran_to_ugrid(nastran_filename1, ugrid_filename_out=ugrid_filename,
                         properties=None, check_shells=False, check_solids=True, log=log)
        assert os.path.exists(ugrid_filename), ugrid_filename

        nastran_filename2 = os.path.join(nastran_path, 'solid_bending', 'solid_bending2.bdf')
        model = ugrid3d_to_nastran(ugrid_filename, nastran_filename2,
                                   include_shells=True, include_solids=True,
                                   convert_pyram_to_penta=False,
                                   encoding=None, size=16,
                                   is_double=False, log=log)
        assert os.path.exists(nastran_filename2), nastran_filename2

        #tecplot_filename1 = os.path.join(nastran_path, 'solid_bending', 'solid_bending.plt')
        #ugrid3d_to_tecplot_filename(model, tecplot_filename1)
        #assert os.path.exists(tecplot_filename1), tecplot_filename1

        tecplot_filename2 = os.path.join(nastran_path, 'solid_bending', 'solid_bending2.plt')
        tecplot = ugrid_to_tecplot(model)
        tecplot.write_tecplot(tecplot_filename2, res_types=None,
                              is_points=True,
                              adjust_nids=True)
        assert os.path.exists(tecplot_filename2), tecplot_filename2


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
