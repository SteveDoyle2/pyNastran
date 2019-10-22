"""Defines TestUGrid"""
import os
import unittest
import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.bdf.bdf import read_bdf
from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid
from pyNastran.converters.nastran.nastran_to_ugrid3d import merge_ugrid3d_and_bdf_to_ugrid3d_filename
from pyNastran.converters.aflr.ugrid.ugrid3d_to_nastran import ugrid3d_to_nastran
from pyNastran.converters.aflr.ugrid.ugrid3d_to_tecplot import (
    ugrid_to_tecplot, ugrid3d_to_tecplot_filename, read_ugrid)
from pyNastran.converters.format_converter import cmd_line_format_converter

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

        unused_ugrid_model = nastran_to_ugrid(
            nastran_filename1, ugrid_filename_out=ugrid_filename,
            properties=None, check_shells=False, check_solids=True, log=log)
        assert os.path.exists(ugrid_filename), ugrid_filename

        nastran_filename2 = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending2.bdf')
        ugrid_model = ugrid3d_to_nastran(
            ugrid_filename, nastran_filename2,
            include_shells=True, include_solids=True,
            convert_pyram_to_penta=False,
            encoding=None, size=16,
            is_double=False, log=log)
        argv = ['format_converter', 'ugrid', ugrid_filename,
                'nastran', 'shell_solid_bending.bdf']
        cmd_line_format_converter(argv=argv, quiet=True)

        nastran_filename3 = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending3.bdf')
        tris, quads = ugrid_model.skin_solids()
        ugrid_model.tris = tris
        ugrid_model.quads = quads
        ugrid_model.pids = np.ones(len(tris) + len(quads))

        ugrid_model.write_bdf(nastran_filename3)

        unused_bdf_model = read_bdf(nastran_filename3, log=log)
        #print(bdf_model.get_bdf_stats())
        assert os.path.exists(nastran_filename3), nastran_filename3

        #tecplot_filename1 = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending.plt')
        #ugrid3d_to_tecplot_filename(model, tecplot_filename1)
        #assert os.path.exists(tecplot_filename1), tecplot_filename1

        tecplot_filename2 = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending2.plt')
        tecplot, unused_zone = ugrid_to_tecplot(ugrid_model, log=log)
        tecplot.write_tecplot(tecplot_filename2, res_types=None,
                              adjust_nids=True)
        assert os.path.exists(tecplot_filename2), tecplot_filename2

        ugrid_filename_out = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending.b8.ugrid_out')
        pshell_pids_to_remove = []
        merge_ugrid3d_and_bdf_to_ugrid3d_filename(
            ugrid_filename, nastran_filename3, ugrid_filename_out,
            pshell_pids_to_remove,
            update_equivalence=True, tol=0.01, log=log)
        assert os.path.exists(ugrid_filename_out), ugrid_filename_out
        os.remove(nastran_filename2)
        os.remove(nastran_filename3)
        os.remove(tecplot_filename2)
        os.remove(ugrid_filename)
        os.remove(ugrid_filename_out)
        os.remove('shell_solid_bending.bdf')

    def test_ugrid3d_gui_box(self):
        """simple UGRID3D box model"""
        ugrid_filename = os.path.join(UGRID_PATH, 'box.b8.ugrid')
        log = get_logger(level='warning')
        tecplot_filename2 = os.path.join(UGRID_PATH, 'box.plt')
        tecplot_filename3 = os.path.join(UGRID_PATH, 'slice.plt')

        ugrid_model = read_ugrid(ugrid_filename, log=log)
        tecplot = ugrid3d_to_tecplot_filename(ugrid_filename, tecplot_filename2, log=log)
        tecplot, unused_zone = ugrid_to_tecplot(ugrid_filename, log=log)
        tecplot, unused_zone = ugrid_to_tecplot(ugrid_model, log=log)
        tecplot.write_tecplot(tecplot_filename2, res_types=None,
                              adjust_nids=True)
        assert os.path.exists(tecplot_filename2), tecplot_filename2

        argv = ['format_converter', 'ugrid', ugrid_filename, 'tecplot', tecplot_filename3, '-z 0.0']
        cmd_line_format_converter(argv=argv, quiet=True)
        os.remove(tecplot_filename2)
        os.remove(tecplot_filename3)

    def test_ugrid3d_convert(self):
        argv = ['format_converter', 'afrl', 'junk.b8.ugrid', 'stl', 'cart3d.stl']
        with self.assertRaises(NotImplementedError):
            cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'ugrid', 'junk.b8.ugrid', 'abaqus', 'cart3d.stl']
        with self.assertRaises(AssertionError):
            cmd_line_format_converter(argv=argv, quiet=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
