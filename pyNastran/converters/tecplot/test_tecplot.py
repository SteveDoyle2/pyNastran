import os
import unittest
from cpylog import get_logger

import pyNastran
from pyNastran.bdf.bdf import read_bdf
from pyNastran.converters.tecplot.tecplot import read_tecplot, split_headers
from pyNastran.converters.tecplot.tecplot_to_nastran import tecplot_to_nastran_filename
from pyNastran.converters.tecplot.tecplot_to_cart3d import tecplot_to_cart3d_filename
from pyNastran.converters.nastran.nastran_to_tecplot import (
    nastran_to_tecplot, nastran_to_tecplot_filename)
from pyNastran.converters.format_converter import cmd_line_format_converter

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'tecplot', 'models')
NASTRAN_MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class TestTecplot(unittest.TestCase):

    def test_split_headers(self):
        headers = [
            ('Zone I=    17, J=    17, K=     1, F=POINT', 4)
        ]
        for header, nheaders in headers:
            headers = split_headers(header)
            assert len(headers) == nheaders, headers

    def test_tecplot_ascii_models(self):
        tecplot_filenames = [
            'ascii/3dgeom.dat', #  good; multi-zone, geometry
            'ascii/block_febrick_3d.dat', # 3d unstructured block; good
            #'ascii/block_fetet_3d.dat', # bad; no decimal values
            'ascii/channel.dat', # 2d structured point; good
            #'ascii/cylinder_slice.dat', # 3d structured point; block 2 has poor formatting
            #'ascii/cylindrical.dat',  # 3d structured empty lines; good
            'ascii/ell.dat', # 2d; good
            'ascii/humanoid_quad.dat', # good
            'ascii/humanoid_tri.dat', # good

            'ascii/movie.dat',  # csv -> good
            'ascii/multzn2d.dat',  #  2d structured; good
            'ascii/plane_slice.dat',  # 2d structured multi-line; good
            #'ascii/point_febrick_3d_02.dat',  # difficult header and funny write bug; bad
            'ascii/point_fequad_2d.dat',  # 2d; good
            'ascii/point_fetet_3d.dat',  # good
            'ascii/point_fetri_2d_01.dat',  # good
            'ascii/point_fetri_2d_02.dat',  # good
            'ascii/point_fetri_2d_03.dat',  # good

            #'ascii/simp3dbk.dat',  # 3d structured block - bad
            'ascii/simp3dpt.dat', #  good
            #'ascii/simpscat.dat', #  bad -> text
            #'ascii/simpxy.dat',  # no xyz; it's a plot -> bad
            #'ascii/simpxy2.dat',  # no xyz; it's a plot -> bad
            'ascii/tiny.dat',  # good
        ]
        log = get_logger(log=None, level='warning', encoding='utf-8')
        junk_plt = os.path.join(MODEL_PATH, 'junk.plt')
        for fname in tecplot_filenames:
            tecplot_filename = os.path.join(MODEL_PATH, fname)
            #print(fname)
            log.info('read %r' % fname)
            model = read_tecplot(tecplot_filename, log=log)
            str(model)
            model.write_tecplot(junk_plt, res_types=None, adjust_nids=True)
        os.remove(junk_plt)

    def test_tecplot_01(self):
        """CTRIA3 elements"""
        log = get_logger(level='warning')
        tecplot_filename1 = os.path.join(MODEL_PATH, 'ascii', 'point_fetri_2d_02.dat')
        #tecplot_filename2 = os.path.join(MODEL_PATH, 'ascii', 'point_fetri_2d_02.dat_out')

        unused_tecplot = read_tecplot(tecplot_filename1, log=log)
        #tecplot.write_tecplot(tecplot_filename2, res_types=None,
                              #is_points=True, adjust_nids=True)
        #os.remove(tecplot_filename2)

        #tecplot_to_cart3d_filename()
        argv = ['format_converter', 'tecplot', tecplot_filename1, 'cart3d', 'cart3d.tri']
        cmd_line_format_converter(argv=argv, quiet=True)
        os.remove('cart3d.tri')

        argv = ['format_converter', 'tecplot', tecplot_filename1, 'stl', 'cart3d.stl']
        cmd_line_format_converter(argv=argv, quiet=True)
        os.remove('cart3d.stl')

    def test_tecplot_02(self):
        """CTETRA10 elements"""
        log = get_logger(level='warning')
        nastran_filename1 = os.path.join(NASTRAN_MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        nastran_filename2 = os.path.join(NASTRAN_MODEL_PATH, 'solid_bending', 'solid_bending2.bdf')
        tecplot_filename = os.path.join(NASTRAN_MODEL_PATH, 'solid_bending', 'solid_bending.plt')
        tecplot_filename2 = os.path.join(NASTRAN_MODEL_PATH, 'solid_bending', 'solid_bending2.plt')
        unused_tecplot = nastran_to_tecplot_filename(nastran_filename1, tecplot_filename, log=log)
        #tecplot.write_tecplot(tecplot_filename)
        tecplot_to_nastran_filename(tecplot_filename, nastran_filename2, log=log)
        #os.remove(nastran_filename2)
        #os.remove(tecplot_filename)

        bdf_model = read_bdf(nastran_filename1, log=log)
        unused_tecplot = nastran_to_tecplot(bdf_model)

        argv = ['format_converter', 'tecplot', tecplot_filename, 'tecplot', tecplot_filename2]
        cmd_line_format_converter(argv=argv, quiet=True)

    def test_tecplot_03(self):
        log = get_logger(level='warning')
        nastran_filename = os.path.join(NASTRAN_MODEL_PATH, 'elements', 'static_elements.bdf')
        tecplot_filename = os.path.join(NASTRAN_MODEL_PATH, 'elements', 'static_elements.plt')
        unused_tecplot = nastran_to_tecplot_filename(nastran_filename, tecplot_filename, log=log)
        #tecplot2 = read_tecplot(tecplot_filename)

        bdf_model = read_bdf(nastran_filename, log=log)
        with self.assertRaises(RuntimeError):
            unused_tecplot = nastran_to_tecplot(bdf_model)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
