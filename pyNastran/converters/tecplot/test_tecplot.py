import os
from io import StringIO
import unittest
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import read_bdf
from pyNastran.converters.tecplot.tecplot import read_tecplot
from pyNastran.converters.tecplot.read_ascii import split_headers
from pyNastran.converters.tecplot.tecplot_to_nastran import tecplot_to_nastran_filename
from pyNastran.converters.tecplot.tecplot_to_cart3d import tecplot_to_cart3d_filename
from pyNastran.converters.nastran.nastran_to_tecplot import (
    nastran_to_tecplot, nastran_to_tecplot_filename)
from pyNastran.converters.tecplot.utils import merge_tecplot
from pyNastran.converters.format_converter import cmd_line_format_converter

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'tecplot', 'models')
NASTRAN_MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class TestTecplot(unittest.TestCase):

    def test_split_headers(self):
        log = SimpleLogger(level='info', encoding='utf-8')
        headers = [
            ('Zone I=    17, J=    17, K=     1, F=POINT', 4)
        ]
        for header, nheaders in headers:
            headers = split_headers(header, log)
            assert len(headers) == nheaders, headers

    def test_tecplot_binary_models(self):
        tecplot_filenames = [
            'binary/point_febrick_3d_02.plt',  # works
        ]
        log = SimpleLogger(level='warning', encoding='utf-8')
        binary_plt = os.path.join(MODEL_PATH, 'ascii.plt')
        ascii_plt = os.path.join(MODEL_PATH, 'binary.plt')
        for fname in tecplot_filenames:
            tecplot_filename = os.path.join(MODEL_PATH, fname)
            #print(fname)
            log.info(f'read {fname!r}')
            model = read_tecplot(tecplot_filename, log=log)
            str(model)
            model.write_tecplot_ascii(ascii_plt)
            model.write_tecplot_binary(binary_plt)
            read_tecplot(ascii_plt, log=log)
            read_tecplot(binary_plt, log=log)
        #os.remove(junk_plt)

    def test_tecplot_ascii_structured(self):
        tecplot_filenames = [
            'ascii/3dgeom.dat', #  good; multi-zone, geometry
            'ascii/channel.dat', # 2d structured point; good
            #'ascii/cylinder_slice.dat', # 3d structured point; block 2 has poor formatting
            #'ascii/cylindrical.dat',  # 3d structured empty lines; good

            'ascii/movie.dat',  # csv -> good
            'ascii/multzn2d.dat',  #  2d structured; good
            'ascii/plane_slice.dat',  # 2d structured multi-line; good
            #'ascii/simp3dbk.dat',  # 3d structured block - bad
            'ascii/simp3dpt.dat', #  good
        ]
        log = SimpleLogger(level='warning', encoding='utf-8')
        ascii_plt = os.path.join(MODEL_PATH, 'ascii.plt')
        binary_plt = os.path.join(MODEL_PATH, 'binary.plt')
        #junk_plt = os.path.join(MODEL_PATH, 'junk.plt')
        for fname in tecplot_filenames:
            tecplot_filename = os.path.join(MODEL_PATH, fname)
            #print(fname)
            log.info(f'read {fname}')
            model = read_tecplot(tecplot_filename, log=log)
            str(model)
            assert model.zones[0].is_structured, f'{fname} {model.zones[0]}'
            model.write_tecplot_ascii(ascii_plt)
            #model.write_tecplot_binary(binary_plt)
            #read_tecplot(ascii_plt, log=log)
            #read_tecplot(binary_plt, log=log)
            #model.write_tecplot(junk_plt, res_types=None, adjust_nids=True)
        #os.remove(junk_plt)

    def test_merge_tecplot(self):
        """tests merge_tecplot"""
        tecplot_filenames = [
            #os.path.join(MODEL_PATH, 'ascii/humanoid_quad.dat'),
            #os.path.join(MODEL_PATH, 'ascii/humanoid_quad.dat'),
            os.path.join(MODEL_PATH, 'ascii/humanoid_tri.dat'),
            os.path.join(MODEL_PATH, 'ascii/humanoid_tri.dat'),
        ]
        tecplot_filename_out = os.path.join(MODEL_PATH, 'junk.plt')
        log = SimpleLogger(level='debug', encoding='utf-8')
        args = tecplot_filenames + ['-b', '-o', tecplot_filename_out]
        merge_tecplot(argv=args)

        args = tecplot_filenames + ['-v']
        with self.assertRaises(SystemExit):
            merge_tecplot(argv=args)

    def test_tecplot_ascii_unstructured(self):
        tecplot_filenames = [
            'ascii/block_febrick_3d.dat', # 3d unstructured block; good
            'ascii/block_fetet_3d.dat', # bad; no decimal values
            'ascii/humanoid_quad.dat', # good
            'ascii/humanoid_tri.dat', # good
            'ascii/ell.dat', # 2d; good

            #'ascii/point_febrick_3d_02.dat',  # difficult header and funny write bug; bad
            'ascii/point_fetet_3d.dat',  # good
            'ascii/point_fetri_2d_01.dat',  # good
            'ascii/point_fetri_2d_02.dat',  # good
            'ascii/point_fetri_2d_03.dat',  # good

            #'ascii/simpscat.dat', #  bad -> text
            #'ascii/simpxy.dat',  # no xyz; it's a plot -> bad
            #'ascii/simpxy2.dat',  # no xyz; it's a plot -> bad

            'ascii/point_fequad_2d.dat',  # 2d; good
            'ascii/tiny.dat',  # good
        ]
        log = SimpleLogger(level='warning', encoding='utf-8')
        ascii_plt = os.path.join(MODEL_PATH, 'ascii.plt')
        binary_plt = os.path.join(MODEL_PATH, 'binary.plt')
        for fname in tecplot_filenames:
            tecplot_filename = os.path.join(MODEL_PATH, fname)
            #print(fname)
            log.info(f'read {fname}')
            model = read_tecplot(tecplot_filename, log=log)
            str(model)
            assert model.zones[0].is_unstructured, f'{fname} {model.zones[0]}'
            model.write_tecplot_ascii(ascii_plt)
            model.write_tecplot_binary(binary_plt)
            read_tecplot(ascii_plt, log=log)
            read_tecplot(binary_plt, log=log)
            #model.write_tecplot(junk_plt, res_types=None, adjust_nids=True)
        #os.remove(junk_plt)

    #def test_tecplot_ascii_unstructured_no_binary(self):
        #tecplot_filenames = [
        #]
        #log = SimpleLogger(level='warning', encoding='utf-8')
        #ascii_plt = os.path.join(MODEL_PATH, 'ascii.plt')
        #binary_plt = os.path.join(MODEL_PATH, 'binary.plt')
        #for fname in tecplot_filenames:
            #tecplot_filename = os.path.join(MODEL_PATH, fname)
            ##print(fname)
            #log.info(f'read {fname}')
            #model = read_tecplot(tecplot_filename, log=log)
            #str(model)
            #assert model.zones[0].is_unstructured, f'{fname} {model.zones[0]}'
            #model.write_tecplot_ascii(ascii_plt)
            #model.write_tecplot_binary(binary_plt)
            #read_tecplot(ascii_plt, log=log)
            #with self.assertRaises(Exception):
            #read_tecplot(binary_plt, log=log)
            #model.write_tecplot(junk_plt, res_types=None, adjust_nids=True)
        #os.remove(junk_plt)

    def test_tecplot_ctria3(self):
        """CTRIA3 elements"""
        log = SimpleLogger(level='warning')
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

    def test_tecplot_ctetra(self):
        """CTETRA10 elements"""
        log = SimpleLogger(level='warning')
        model_path = os.path.join(NASTRAN_MODEL_PATH, 'solid_bending')
        nastran_filename1 = os.path.join(model_path, 'solid_bending.bdf')
        nastran_filename2 = os.path.join(model_path, 'solid_bending2.bdf')
        tecplot_filename = os.path.join(model_path, 'solid_bending.plt')
        tecplot_filename2 = os.path.join(model_path, 'solid_bending2.plt')
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
        log = SimpleLogger(level='warning')
        nastran_filename = os.path.join(NASTRAN_MODEL_PATH, 'elements', 'static_elements.bdf')
        tecplot_filename = os.path.join(NASTRAN_MODEL_PATH, 'elements', 'static_elements.plt')
        unused_tecplot = nastran_to_tecplot_filename(nastran_filename, tecplot_filename, log=log)
        #tecplot2 = read_tecplot(tecplot_filename)

        bdf_model = read_bdf(nastran_filename, log=log)
        with self.assertRaises(RuntimeError):
            unused_tecplot = nastran_to_tecplot(bdf_model)

    def test_tecplot_360_point_viscous(self):
        lines = [
            'TITLE     = "tecplot geometry and solution file"',
            'VARIABLES = "x"',
            '"y"',
            '"z"',
            '"cp"',
            '"cf_x"',
            '"cf_y"',
            '"cf_z"',
            'ZONE T="\"boundary 1 nose-fuselage\""',
            ' STRANDID=1001, SOLUTIONTIME=5000',
            ' Nodes=4, Elements=1, ZONETYPE=FEQuadrilateral',
            ' DATAPACKING=POINT',
            ' DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )',
            '0. 0. 0. 1. 2. 3. 4.',
            '0. 1. 0. 11. 12. 13. 14.',
            '0. 0. 1. 21. 22. 23. 24.',
            '0. 1. 1. 31. 32. 33. 34.',
            '1 2 4 3',

            'ZONE T="\"boundary 3 wing-tail\""',
            ' STRANDID=1003, SOLUTIONTIME=5000',
            ' Nodes=4, Elements=1, ZONETYPE=FEQuadrilateral',
            ' DATAPACKING=POINT',
            ' DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )',
            '1. 0. 0. 1. 2. 3. 4.',
            '1. 1. 0. 11. 12. 13. 14.',
            '1. 0. 1. 21. 22. 23. 24.',
            '1. 1. 1. 31. 32. 33. 34.',
            '1 2 4 3',

            'ZONE T="\"boundary 4 fairing\""',
            ' STRANDID=1003, SOLUTIONTIME=5000',
            ' Nodes=3, Elements=1, ZONETYPE=FETriangle',
            ' DATAPACKING=POINT',
            ' DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )',
            '1. 0. 0. 1. 2. 3. 4.',
            '1. 1. 0. 11. 12. 13. 14.',
            '1. 0. 1. 21. 22. 23. 24.',
            '1 2 3',
        ]
        tecplot_file = StringIO()
        tecplot_file.write('\n'.join(lines))
        tecplot_file.seek(0)
        log = SimpleLogger(level='warning')
        model = read_tecplot(tecplot_file, use_cols=None, dtype=None,
                             filetype='ascii', log=log, debug=False)
        nodes, tris, quads, tets, hexas, zone_ids, names = model.stack_geometry()
        assert nodes.shape == (11, 3)
        assert tris.shape == (1, 3)
        assert quads.shape == (2, 4)
        assert len(tets) == 0
        assert len(hexas) == 0
        cart3d_filename = os.path.join(MODEL_PATH, 'junk.tri')
        cart3d_model = tecplot_to_cart3d_filename(
            model, cart3d_filename, remove_degenerate_tris=False,
            log=log, debug=True)
        assert cart3d_model.nodes.shape == (11, 3), cart3d_model.nodes.shape
        assert cart3d_model.elements.shape == (5, 3), cart3d_model.elements.shape

        cart3d_model = tecplot_to_cart3d_filename(
            model, cart3d_filename, remove_degenerate_tris=True,
            log=log, debug=True)
        assert cart3d_model.nodes.shape == (11, 3), cart3d_model.nodes.shape
        assert cart3d_model.elements.shape == (3, 3), cart3d_model.elements.shape

    def test_tecplot_360_point_euler(self):
        lines = [
            'TITLE     = "Solution mapped to surface triangulation"',
            'VARIABLES = "x"',
            '"y"',
            '"z"',
            '"Cp"',
            'DATASETAUXDATA Common.AngleOfAttack="5.000000"',
            'DATASETAUXDATA Common.DensityVar="5"',
            'DATASETAUXDATA Common.GasConstant="0.7142857143"',
            'DATASETAUXDATA Common.PressureVar="9"',
            'DATASETAUXDATA Common.ReferenceMachNumber="0.800000"',
            'DATASETAUXDATA Common.UVar="6"',
            'DATASETAUXDATA Common.VectorVarsAreVelocity="TRUE"',
            'DATASETAUXDATA Common.VVar="7"',
            'DATASETAUXDATA Common.WVar="8"',
            'ZONE T="Surface"',
            'STRANDID=0, SOLUTIONTIME=0',
            'Nodes=4, Elements=1, ZONETYPE=FETriangle',
            'DATAPACKING=POINT',
            'DT=(DOUBLE DOUBLE DOUBLE DOUBLE )',
            '0. 0. 0. 1.',
            '0. 1. 0. 2.',
            '0. 0. 1. 3.',
            '0. 1. 1. 4.',
            '1 2 3',
        ]
        file = StringIO()
        file.write('\n'.join(lines))
        file.seek(0)
        log = SimpleLogger(level='warning')
        model = read_tecplot(file, use_cols=None, dtype=None,
                             filetype='ascii', log=log, debug=False)
        nodes, tris, quads, tets, hexas, zone_ids, names = model.stack_geometry()
        assert nodes.shape == (4, 3)
        assert tris.shape == (1, 3)
        assert len(quads) == 0
        assert len(tets) == 0
        assert len(hexas) == 0

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
