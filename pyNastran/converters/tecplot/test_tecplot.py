import os
from io import StringIO
from pathlib import Path
import unittest

import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import read_bdf
from pyNastran.converters.tecplot.tecplot import read_tecplot
from pyNastran.converters.tecplot.read_ascii import split_headers, parse_varlocation
from pyNastran.converters.tecplot.tecplot_to_nastran import tecplot_to_nastran_filename
from pyNastran.converters.tecplot.tecplot_to_cart3d import tecplot_to_cart3d_filename
from pyNastran.converters.nastran.nastran_to_tecplot import (
    nastran_to_tecplot, nastran_to_tecplot_filename)
from pyNastran.converters.tecplot.utils import merge_tecplot
from pyNastran.converters.format_converter import cmd_line_format_converter

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / 'converters' / 'tecplot' / 'models'
NASTRAN_MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class TestTecplot(unittest.TestCase):

    def test_mode1(self):
        tecplot_filename = (PKG_PATH / 'bdf' / 'cards' / 'aero' /
                            'examples' / 'flutter' / 'case1' / 'mode1.dat')
        read_tecplot(tecplot_filename, filetype='ascii')

    def test_split_headers(self):
        log = SimpleLogger(level='info', encoding='utf-8')
        headers = [
            ('Zone I=    17, J=    17, K=     1, F=POINT', 4)
        ]
        for header, nheaders in headers:
            headers = split_headers(header, log)
            assert len(headers) == nheaders, headers

    def test_split_headers_varlocation(self):
        """commas inside VARLOCATION=(...) don't split the header"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        header = ('ZONE T="wing", n=6, e=3, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK,'
                  'VARLOCATION=([1, 2, 3]=NODAL,[4-6]=CELLCENTERED)')
        headers = split_headers(header, log)
        assert len(headers) == 6, headers
        assert headers[-1] == 'VARLOCATION=([1, 2, 3]=NODAL,[4-6]=CELLCENTERED)', headers

    def test_parse_varlocation(self):
        """
        Tests:
         - explicit lists, ranges, and unlisted (default NODAL) variables
         - bad locations / out of range ids / unparseable values raise
        """
        is_cell = parse_varlocation('([1, 2, 3]=NODAL,[4, 5, 6]=CELLCENTERED)', 6)
        assert is_cell.tolist() == [False, False, False, True, True, True], is_cell
        is_cell = parse_varlocation('([3-5,7]=CELLCENTERED)', 7)
        assert is_cell.tolist() == [False, False, True, True, True, False, True], is_cell
        is_cell = parse_varlocation('( [4] = cellcentered )', 4)
        assert is_cell.tolist() == [False, False, False, True], is_cell

        with self.assertRaises(ValueError):
            parse_varlocation('([4]=FACE)', 4)
        with self.assertRaises(ValueError):
            parse_varlocation('([5]=CELLCENTERED)', 4)
        with self.assertRaises(ValueError):
            parse_varlocation('CELLCENTERED', 4)

    def test_tecplot_ascii_varlocation(self):
        """
        Reads an FE BLOCK zone with NODAL xyz and CELLCENTERED results.

        Tests (exact; values are short decimals):
         - nodal xyz -> zone_data; cell-centered results -> element_data
         - values split across lines / mixed with commas
         - a 2nd zone that reuses the VARIABLES
         - cell-centered before nodal raises NotImplementedError
        """
        log = SimpleLogger(level='warning', encoding='utf-8')
        lines = (
            'TITLE = "varlocation"\n'
            'VARIABLES = "X", "Y", "Z", "Cp", "Region"\n'
            'ZONE T="tris", N=4, E=2, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK,\n'
            'VARLOCATION=([1-3]=NODAL,[4,5]=CELLCENTERED)\n'
            '0.0 1.0\n1.0 0.0\n'
            '0.0 0.0 1.0 1.0\n'
            '0.0, 0.0, 0.0, 0.0\n'
            '-0.5 0.25\n'
            '3.0 4.0\n'
            '1 2 3\n'
            '1 3 4\n'
            'ZONE T="quad", N=4, E=1, ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK,\n'
            'VARLOCATION=([1, 2, 3]=NODAL,[4, 5]=CELLCENTERED)\n'
            '0.0 1.0 1.0 0.0\n'
            '0.0 0.0 1.0 1.0\n'
            '2.0 2.0 2.0 2.0\n'
            '0.75\n'
            '7.0\n'
            '1 2 3 4\n'
        )
        model = read_tecplot(StringIO(lines), filetype='ascii', log=log)
        assert len(model.zones) == 2, model.zones
        tri_zone, quad_zone = model.zones
        assert tri_zone.variables == ['X', 'Y', 'Z', 'Cp', 'Region'], tri_zone.variables
        assert np.array_equal(tri_zone.xyz, [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]]), tri_zone.xyz
        assert np.array_equal(tri_zone.element_data, [[-0.5, 3.0], [0.25, 4.0]]), tri_zone.element_data
        assert np.array_equal(tri_zone.tri_elements, [[0, 1, 2], [0, 2, 3]]), tri_zone.tri_elements
        assert np.array_equal(quad_zone.xyz[:, 2], [2., 2., 2., 2.]), quad_zone.xyz
        assert np.array_equal(quad_zone.element_data, [[0.75, 7.0]]), quad_zone.element_data
        assert np.array_equal(quad_zone.quad_elements, [[0, 1, 2, 3]]), quad_zone.quad_elements

        bad_lines = lines.replace('[1-3]=NODAL,[4,5]=CELLCENTERED', '[1]=CELLCENTERED')
        with self.assertRaises(NotImplementedError):
            read_tecplot(StringIO(bad_lines), filetype='ascii', log=log)

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
            model.is2d
            model.is3d
            model.result_variables
            model.write_tecplot_ascii(ascii_plt)
            model.write_tecplot_binary(binary_plt)
            read_tecplot(ascii_plt, log=log)
            read_tecplot(binary_plt, log=log)
        #os.remove(junk_plt)

    def test_tecplot_ascii_structured(self):
        tecplot_filenames = [
            'ascii/3dgeom.dat', #  good; multi-zone, geometry
            'ascii/channel.dat', # 2d structured point; good

            #    'Zone I=    17, J=    17K=     1, F=POINT'
            # -> 'Zone I=    17, J=    17,K=     1, F=POINT'
            #'ascii/cylinder_slice.dat', # 3d structured point; block 2 has poor formatting
            'ascii/cylindrical.dat',  # 3d structured empty lines at end; good

            'ascii/movie.dat',  # csv -> good
            'ascii/multzn2d.dat',  #  2d structured; good
            'ascii/plane_slice.dat',  # 2d structured multi-line; good
            'ascii/simp3dbk.dat',  # 3d structured block - good
            'ascii/simp3dpt.dat', #  3d structured point - good
        ]
        log = SimpleLogger(level='warning', encoding='utf-8')
        #log = SimpleLogger(level='debug', encoding='utf-8')
        ascii_plt = os.path.join(MODEL_PATH, 'ascii.plt')
        binary_plt = os.path.join(MODEL_PATH, 'binary.plt')
        #junk_plt = os.path.join(MODEL_PATH, 'junk.plt')
        for fname in tecplot_filenames:
            tecplot_filename = os.path.join(MODEL_PATH, fname)
            #print(fname)
            log.info(f'read {fname}')
            model = read_tecplot(tecplot_filename, log=log)
            str(model)
            model.is2d
            model.is3d
            model.result_variables
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
        log = SimpleLogger(level='warning', encoding='utf-8')
        args = tecplot_filenames + ['-b', '-o', tecplot_filename_out]
        merge_tecplot(argv=args, log=log)

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
            model.is2d
            model.is3d
            model.result_variables
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

    def test_tecplot_ctria3_split(self):
        """CTRIA3 elements"""
        log = SimpleLogger(level='debug')

        lines = [
            'TITLE     = "Solution mapped to surface triangulation"',
            'VARIABLES = "x"',
            '"y"',
            '"z"',
            '"Cp"',
            'ZONE T="Surface"',
            'Nodes=4, Elements=1, ZONETYPE=FETriangle',
            'DATAPACKING=POINT',
            '0. 0. 0. 1.',
            '0. 1. 0. 2.',
            '0. 0. 1. 3.',
            '0. 1. 1. 4.',  # unused node
            '1 2 3',
        ]
        tecplot_filename1 = StringIO()
        tecplot_filename1.write('\n'.join(lines))
        tecplot_filename1.seek(0)
        tecplot_filename1 = os.path.join(MODEL_PATH, 'ascii', 'point_fetri_2d_02.dat')
        #tecplot_filename1 = os.path.join(MODEL_PATH, 'ascii', 'humanoid_tri.dat')

        tecplot_filename2 = os.path.join(MODEL_PATH, 'ascii', 'tria_split.plt')
        model = read_tecplot(tecplot_filename1, log=log)
        #type_limit = ('tri')
        model.demote_elements()
        model.split_elements(ntri_nodes=3)
        model.write_tecplot_binary(tecplot_filename2)
        model.is2d
        model.is3d
        model.result_variables

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
        assert len(tets) == 0, tets
        assert len(hexas) == 0, hexas
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
