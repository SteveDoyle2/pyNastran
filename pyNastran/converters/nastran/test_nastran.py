"""tests the Nastran converters"""
import os
import unittest
import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf

from pyNastran.converters.format_converter import cmd_line_format_converter
from pyNastran.converters.nastran.nastran_to_cart3d import nastran_to_cart3d, nastran_to_cart3d_filename
from pyNastran.converters.nastran.nastran_to_stl import nastran_to_stl
from pyNastran.converters.nastran.nastran_to_surf import nastran_to_surf, clear_out_solids
from pyNastran.converters.nastran.nastran_to_tecplot import nastran_to_tecplot, nastran_to_tecplot_filename
from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid
from pyNastran.converters.aflr.ugrid.ugrid_reader import read_ugrid
from pyNastran.converters.cart3d.cart3d import read_cart3d
from pyNastran.bdf.mesh_utils.skin_solid_elements import write_skin_solid_faces

import pyNastran.converters.nastran.nastran_to_ugrid3d
from pyNastran.converters.tecplot.tecplot_to_nastran import nastran_tables_to_tecplot_filenames

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')
DIRNAME = os.path.dirname(__file__)

class FakeCase:
    def __init__(self, times: np.ndarray):
        self._times = times
        self.headers = ['a', 'b']
        self.title = 'title'
        self.subtitle = 'subtitle'
        ntimes = 4
        nresults = len(self.headers)
        self.data = np.zeros((ntimes, 36, 2))


class TestNastran(unittest.TestCase):
    def test_nastran_to_cart3d_se2(self):
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 2, [1, 2, 3, 4])
        model.add_ctria3(2, 2, [2, 3, 4])

        model.add_cquadr(3, 2, [1, 2, 3, 4])
        model.add_ctriar(4, 2, [2, 3, 4])
        model.add_pshell(2, mid1=100, t=0.1)
        model.add_mat1(100, 3.0e7, None, 0.3)
        model.add_cbar(101, 200, [3, 2], [1., 1., 1.], None)
        model.add_pbarl(200, 100, 'ROD', [1.])
        model.card_count = {
            'GRID': 4,
            'CQUAD4': 1,
            'CQUADR': 1,
            'CTRIA3': 1,
            'CTRIAR': 1,
            'CBAR': 1,
        }
        cart3d = nastran_to_cart3d(model)
        cart3d.flip_model()
        assert len(cart3d.nodes) == 4
        assert len(cart3d.elements) == 6
        cart3d.remove_elements([0], remove_associated_nodes=True)
        assert len(cart3d.nodes) == 4
        assert len(cart3d.elements) == 5
        cart3d.keep_elements([2], remove_associated_nodes=True)
        assert len(cart3d.nodes) == 3
        assert len(cart3d.elements) == 1

        bdf_filename = os.path.join(DIRNAME, 'nastran_to_cart3d.bdf')
        model.write_bdf(bdf_filename)
        cart3d_filename = os.path.join(DIRNAME, 'nastran_to_cart3d.tri')
        nastran_to_cart3d_filename(bdf_filename, cart3d_filename)

    def test_nastran_to_cart3d(self):
        model = BDF(debug=False)
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        model.add_grid(4, [1., 1., 0.])
        model.add_grid(5, [0., 1., 0.])
        model.add_cquad4(10, 2, [2, 3, 4, 5])
        model.add_ctria3(12, 2, [2, 3, 4])

        model.add_cquadr(100, 2, [2, 3, 4, 5])
        model.add_ctriar(121, 2, [2, 3, 4])
        model.add_pshell(2, mid1=100, t=0.1)
        model.add_mat1(100, 3.0e7, None, 0.3)
        model.add_cbar(101, 200, [3, 2], [1., 1., 1.], None)
        model.add_pbarl(200, 100, 'ROD', [1.])
        model.card_count = {
            'GRID': 4,
            'CQUAD4': 1,
            'CQUADR': 1,
            'CTRIA3': 1,
            'CTRIAR': 1,
            'CBAR': 1,
        }
        cart3d = nastran_to_cart3d(model)
        cart3d.flip_model()
        assert len(cart3d.nodes) == 4
        assert len(cart3d.elements) == 6
        cart3d.remove_elements([0], remove_associated_nodes=True)
        assert len(cart3d.nodes) == 4
        assert len(cart3d.elements) == 5
        cart3d.keep_elements([2], remove_associated_nodes=True)
        assert len(cart3d.nodes) == 3
        assert len(cart3d.elements) == 1

        bdf_filename = os.path.join(DIRNAME, 'nastran_to_cart3d.bdf')
        model.write_bdf(bdf_filename)
        cart3d_filename = os.path.join(DIRNAME, 'nastran_to_cart3d.tri')
        nastran_to_cart3d_filename(bdf_filename, cart3d_filename)

    def test_nastran_to_tecplot(self):
        """tests a large number of elements and results in SOL 101"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
        tecplot_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.plt')
        tecplot_filename2 = os.path.join(MODEL_PATH, 'elements', 'static_elements2.plt')
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = read_bdf(bdf_filename, log=log)
        with self.assertRaises(RuntimeError):
            nastran_to_tecplot(model)
        nastran_to_tecplot_filename(bdf_filename, tecplot_filename, log=log)

        argv = ['format_converter', 'nastran', bdf_filename, 'tecplot', tecplot_filename2]
        with self.assertRaises(RuntimeError):
            cmd_line_format_converter(argv=argv, quiet=True)

    def test_nastran_to_ugrid_01(self):
        bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')

        size = 8
        debug = False
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = read_bdf(bdf_filename, log=log, debug=debug)
        #log = model.log
        #model.get_element_faces()
        skin_bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin.bdf')
        write_skin_solid_faces(model, skin_bdf_filename, write_solids=True,
                               write_shells=True,
                               size=size, is_double=False, encoding=None)

        bdf_model = read_bdf(skin_bdf_filename, log=log, debug=debug)
        ugrid_filename_out = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin.b8.ugrid')
        ugrid_filename_out2 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin2.b8.ugrid')
        nastran_to_ugrid(bdf_model, ugrid_filename_out, properties=None,
                         check_shells=True, check_solids=True)
        ugrid = read_ugrid(ugrid_filename_out, encoding=None, log=log,
                           debug=debug)

        skin_bdf_filename2 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin2.bdf')
        skin_cart3d_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin2.tri')
        skin_cart3d_filename3 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin3.tri')
        skin_stl_filename3 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin3.stl')

        #msg += "  format_converter nastran   <INPUT> <format2> <OUTPUT> [-o <OP2>] --no_xref\n"
        #msg += "  format_converter <format1> <INPUT> tecplot   <OUTPUT> [-r RESTYPE...] [-b] [--block] [-x <X>] [-y <Y>] [-z <Z>] [--scale SCALE]\n"
        #msg += "  format_converter <format1> <INPUT> stl       <OUTPUT> [-b]  [--scale SCALE]\n"
        #msg += "  format_converter cart3d    <INPUT> <format2> <OUTPUT> [-b]  [--scale SCALE]\n"
        #msg += "  format_converter <format1> <INPUT> <format2> <OUTPUT> [--scale SCALE]\n"
        argv = ['format_converter', 'nastran', bdf_filename, 'ugrid', ugrid_filename_out2]
        with self.assertRaises(RuntimeError):
            cmd_line_format_converter(argv=argv, quiet=True)

        #argv = ['format_converter', 'nastran', bdf_filename, 'cart3d', skin_cart3d_filename3]
        #cmd_line_format_converter(argv=argv)

        #argv = ['format_converter', 'nastran', bdf_filename, 'stl', skin_stl_filename3]
        #cmd_line_format_converter(argv=argv)

        ugrid.write_bdf(skin_bdf_filename2, include_shells=True, include_solids=True,
                        convert_pyram_to_penta=True, encoding=None,
                        size=size, is_double=False)
        read_bdf(skin_bdf_filename2, log=log, debug=debug)

        with self.assertRaises(NotImplementedError):
            nastran_to_cart3d_filename(skin_bdf_filename2, skin_cart3d_filename)

        ugrid.write_bdf(skin_bdf_filename2, include_shells=True, include_solids=False,
                        convert_pyram_to_penta=True, encoding=None,
                        size=size, is_double=False)

        nastran_to_cart3d_filename(skin_bdf_filename2, skin_cart3d_filename)
        read_cart3d(skin_cart3d_filename, log=log)

        os.remove(ugrid_filename_out)
        os.remove(skin_bdf_filename)
        os.remove(skin_bdf_filename2)
        os.remove(skin_cart3d_filename)

    def test_nastran_to_stl(self):
        """tests nastran_to_stl"""
        bdf_filename = os.path.join(MODEL_PATH, 'plate', 'plate.bdf')
        stl_filename = os.path.join(MODEL_PATH, 'plate', 'plate.stl')
        log = SimpleLogger(level='warning', encoding='utf-8')
        nastran_to_stl(bdf_filename, stl_filename, is_binary=False, log=log)

    def test_nastran_to_tecplot_case(self):
        bdf_filename = os.path.join(MODEL_PATH, 'plate', 'plate.bdf')
        log = SimpleLogger(level='warning', encoding='utf-8')
        bdf_model = read_bdf(bdf_filename, log=log, debug=False)

        times = np.arange(2)
        case = FakeCase(times)
        tecplot_filename_base = 'cat%d'
        nastran_tables_to_tecplot_filenames(tecplot_filename_base, bdf_model, case,
                                            variables=None, ivars=None)

    def test_format_converter(self):
        """tests nastran_to_stl"""
        bdf_filename = os.path.join(MODEL_PATH, 'plate', 'plate.bdf')
        bdf_filename2 = os.path.join(MODEL_PATH, 'plate', 'plate2.bdf')

        stl_filename = os.path.join(MODEL_PATH, 'plate', 'plate.stl')
        ugrid_filename = os.path.join(MODEL_PATH, 'plate', 'plate.b8.ugrid')
        cart3d_filename = os.path.join(MODEL_PATH, 'plate', 'plate.tri')
        tecplot_filename = os.path.join(MODEL_PATH, 'plate', 'plate.plt')

        argv = ['format_converter', 'nastran', bdf_filename, 'stl', stl_filename]
        cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'nastran', bdf_filename, 'tecplot', tecplot_filename]
        cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'nastran', bdf_filename, 'ugrid', ugrid_filename]
        with self.assertRaises(RuntimeError):
            cmd_line_format_converter(argv=argv, quiet=True)


        #argv = ['format_converter', 'nastran', bdf_filename, 'cart3d', cart3d_filename]
        #cmd_line_format_converter(argv=argv, quiet=True)
        #os.remove(stl_filename)
        #os.remove(cart3d_filename)
        #os.remove(tecplot_filename)
        # -------------------------
        tecplot_filename2 = os.path.join(MODEL_PATH, 'plate', 'plate2.plt')
        argv = ['format_converter', 'stl', stl_filename, 'nastran', bdf_filename2]
        cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'stl', stl_filename, 'tecplot', tecplot_filename2]
        with self.assertRaises(AssertionError):
            cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'stl', tecplot_filename, 'ugrid', ugrid_filename]
        with self.assertRaises(AssertionError):
            cmd_line_format_converter(argv=argv, quiet=True)

        os.remove(bdf_filename2)
        #os.remove(stl_filename)
        #os.remove(cart3d_filename)
        #os.remove(tecplot_filename)
        # -------------------------
        argv = ['format_converter', 'tecplot', tecplot_filename, 'nastran', bdf_filename2]
        cmd_line_format_converter(argv=argv, quiet=True)

        #argv = ['format_converter', 'tecplot', tecplot_filename, 'stl', stl_filename]
        #cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'tecplot', tecplot_filename, 'ugrid', ugrid_filename]
        with self.assertRaises(AssertionError):
            cmd_line_format_converter(argv=argv, quiet=True)

        os.remove(bdf_filename2)
        os.remove(stl_filename)
        #os.remove(cart3d_filename)
        os.remove(tecplot_filename)

    def test_clear_out_solids(self):
        """tests clear_out_solids"""
        deck = (
            "$ pyNastran: punch=True\n"
            "GRID,1\n"
            "GRID,2\n"
            "GRID,3\n"
            "GRID,4\n"
            "GRID,5\n"
            "GRID,6\n"
            "GRID,7\n"
            "GRID,8\n"

            "GRID,9\n"
            "GRID,10\n"
            "GRID,11\n"
            "GRID,12\n"
            "CHEXA,1,1, 5,6,7,8,9,10,\n"
            ",7,8\n"
            "CQUAD4,2,200, 1,2,3,4\n"
            # doesn't work
            #"CHEXA,1,1, 1,2,3,4,5,6,\n"
            #",7,8\n"
            #"CQUAD4,2,200, 8,9,10,11\n"
            "PSHELL,200,1000,0.1\n"
            "PSOLID,100,1000\n"
            "MAT1,1000,3.0e7,,0.3\n"
        )

        bdf_filename = 'deck.bdf'
        bdf_clean_filename = 'clean.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(deck)

        log = SimpleLogger(level='warning', encoding='utf-8')
        model = read_bdf(bdf_filename, xref=False, log=log)
        clear_out_solids(model, bdf_clean_filename, renumber=True,
                         equivalence=False, equivalence_tol=0.01)

        model = read_bdf(bdf_clean_filename, log=log)
        assert len(model.nodes) == 4, len(model.nodes)
        assert len(model.elements) == 1, len(model.elements)
        os.remove(bdf_filename)
        os.remove(bdf_clean_filename)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
