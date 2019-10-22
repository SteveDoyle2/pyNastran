"""tests non-gui related Cart3d class/interface"""
import os
import unittest
from numpy import array_equal, allclose
from cpylog import get_logger

import pyNastran
from pyNastran.converters.cart3d.cart3d import read_cart3d, comp2tri
from pyNastran.converters.cart3d.cart3d_to_nastran import cart3d_to_nastran_filename, cart3d_to_nastran_model
from pyNastran.converters.cart3d.cart3d_to_stl import cart3d_to_stl_filename
from pyNastran.converters.cart3d.cart3d_to_tecplot import cart3d_to_tecplot
from pyNastran.converters.cart3d.input_c3d_reader import read_input_c3d
import pyNastran.converters.cart3d.input_cntl_reader
from pyNastran.converters.format_converter import cmd_line_format_converter


PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'cart3d', 'models')

class TestCart3d(unittest.TestCase):

    def test_cart3d_io_01(self):
        """geometry"""
        lines = (
            "7 6\n"
            "0.000000 0.000000 0.000000\n"
            "1.000000 0.000000 0.000000\n"
            "2.000000 0.000000 0.000000\n"
            "1.000000 1.000000 0.000000\n"
            "2.000000 1.000000 0.000000\n"
            "1.000000 -1.000000 0.000000\n"
            "2.000000 -1.000000 0.000000\n"
            "1 4 2\n"
            "2 4 5\n"
            "2 5 3\n"
            "2 6 1\n"
            "5 6 2\n"
            "5 5 2\n"
            "1\n"
            "2\n"
            "3\n"
            "2\n"
            "4\n"
            "6\n"
        )
        log = get_logger(level='warning', encoding='utf-8')
        infile_name = os.path.join(MODEL_PATH, 'flat_full.tri')
        out_name = os.path.join(MODEL_PATH, 'combined.tri')
        with open(infile_name, 'w') as f:
            f.write(lines)

        model = comp2tri([infile_name, infile_name], out_name, is_binary=False, float_fmt='%6.7f',
                         log=log)
        read_cart3d(out_name, log=log, debug=False)

        cart3d = read_cart3d(infile_name, log=log, debug=False)
        assert len(cart3d.points) == 7, 'npoints=%s' % len(cart3d.points)
        assert len(cart3d.elements) == 6, 'nelements=%s' % len(cart3d.elements)
        assert len(cart3d.regions) == 6, 'nregions=%s' % len(cart3d.regions)
        assert len(cart3d.loads) == 0, 'nloads=%s' % len(cart3d.loads)
        os.remove(out_name)

        argv = ['format_converter', 'cart3d', infile_name, 'cart3d', out_name, '-b']
        cmd_line_format_converter(argv=argv, quiet=True)

        os.remove(infile_name)
        os.remove(out_name)

    def test_cart3d_io_02(self):
        """geometry + results"""
        lines = (
            "5 3 6\n"
            "0. 0. 0.\n"
            "1. 0. 0.\n"
            "2. 0. 0.\n"
            "1. 1. 0.\n"
            "2. 1. 0.\n"
            "1 4 2\n"
            "2 4 5\n"
            "2 5 3\n"
            "1\n"
            "2\n"
            "3\n"
            "1.\n"
            "1. 1. 1. 1. 1.\n"
            "2.\n"
            "2. 2. 2. 2. 2.\n"
            "3.\n"
            "3. 3. 3. 3. 3.\n"
            "4.\n"
            "4. 4. 4. 4. 4.\n"
            "5.\n"
            "5. 5. 5. 5. 5.\n"
        )
        cart3d_filename = os.path.join(MODEL_PATH, 'flat.tri')
        with open(cart3d_filename, 'w') as f:
            f.write(lines)

        log = get_logger(level='warning', encoding='utf-8')
        cart3d = read_cart3d(cart3d_filename, log=log, debug=False,
                            result_names=None)

        assert len(cart3d.points) == 5, 'npoints=%s' % len(cart3d.points)
        assert len(cart3d.elements) == 3, 'nelements=%s' % len(cart3d.elements)
        assert len(cart3d.regions) == 3, 'nregions=%s' % len(cart3d.regions)

        assert len(cart3d.loads) == 14, 'nloads=%s' % len(cart3d.loads)  # was 10
        assert len(cart3d.loads['Cp']) == 5, 'nCp=%s' % len(cart3d.loads['Cp'])

        outfile_name = os.path.join(MODEL_PATH, 'flat.bin.tri')
        cart3d.loads = None
        cart3d.write_cart3d(outfile_name, is_binary=True)
        cnormals = cart3d.get_normals()
        nnormals = cart3d.get_normals_at_nodes(cnormals)
        os.remove(cart3d_filename)
        os.remove(outfile_name)

    def test_cart3d_io_03(self):
        """read/write geometry in ascii/binary"""
        log = get_logger(level='warning', encoding='utf-8')
        infile_name = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')
        outfile_name = os.path.join(MODEL_PATH, 'threePlugs_out.tri')
        outfile_name_bin = os.path.join(MODEL_PATH, 'threePlugs_bin2.tri')
        outfile_name_bin_out = os.path.join(MODEL_PATH, 'threePlugs_bin_out.tri')

        cart3d = read_cart3d(infile_name, log=log, debug=False)
        cart3d.write_cart3d(outfile_name, is_binary=False)
        cart3d.write_cart3d(outfile_name_bin, is_binary=True)

        cart3d_ascii = read_cart3d(outfile_name, log=log, debug=False)
        check_array(cart3d.points, cart3d_ascii.points)
        check_array(cart3d.elements, cart3d_ascii.elements)

        cart3d_bin = read_cart3d(outfile_name_bin, log=log, debug=False)
        check_array(cart3d.points, cart3d_bin.points)
        check_array(cart3d.elements, cart3d_ascii.elements)
        #print(cart3d_bin.points)

        #print('---------------')
        #print(cart3d_bin.points)

        os.remove(outfile_name)
        os.remove(outfile_name_bin)

        cart3d.write_cart3d(outfile_name_bin_out, is_binary=False)
        os.remove(outfile_name_bin_out)

    def test_cart3d_to_stl(self):
        """convert to stl"""
        log = get_logger(level='warning', encoding='utf-8')
        cart3d_filename = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')
        stl_filename = os.path.join(MODEL_PATH, 'threePlugs.stl')

        argv = ['format_converter', 'cart3d', cart3d_filename,
                'stl', stl_filename, '-b', '--scale', '2.0']
        cmd_line_format_converter(argv=argv, quiet=True)
        #cart3d_to_stl_filename(cart3d_filename, stl_filename, log=log)
        os.remove(stl_filename)

    def test_cart3d_to_tecplot(self):
        """convert to tecplot"""
        log = get_logger(level='warning', encoding='utf-8')
        cart3d_filename = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')
        tecplot_filename = os.path.join(MODEL_PATH, 'threePlugs.plt')
        #cart3d_to_tecplot(cart3d_filename, tecplot_filename, log=log)

        argv = ['format_converter', 'cart3d', cart3d_filename,
                'tecplot', tecplot_filename, '--scale', '2.0']
        cmd_line_format_converter(argv=argv, quiet=True)
        os.remove(tecplot_filename)

    def test_cart3d_to_nastran_01(self):
        """convert to nastran small field"""
        log = get_logger(level='warning', encoding='utf-8')
        cart3d_filename = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')
        bdf_filename = os.path.join(MODEL_PATH, 'threePlugs.bdf')
        #cart3d_to_nastran_filename(cart3d_filename, bdf_filename, log=log)

        argv = ['format_converter', 'cart3d', cart3d_filename,
                'nastran', bdf_filename, '-b', '--scale', '2.0']
        cmd_line_format_converter(argv=argv, quiet=True)
        os.remove(bdf_filename)

    def test_cart3d_to_nastran_02(self):
        """convert to nastran large field"""
        log = get_logger(level='warning', encoding='utf-8')
        cart3d_filename = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')
        bdf_filename = os.path.join(MODEL_PATH, 'threePlugs2.bdf')
        model = cart3d_to_nastran_model(cart3d_filename, log=log)
        model.write_bdf(bdf_filename, size=16)
        self.assertAlmostEqual(model.nodes[1].xyz[0], 1.51971436,
                               msg='if this is 0.0, then the assign_type method had the float32 check removed')
        os.remove(bdf_filename)
        #model.write_bdf(out_filename=None, encoding=None, size=8,
                       #is_double=False,
                       #interspersed=True,
                       #enddata=None)

    #def test_cart3d_input_cntl(self):
        #"""tests the input.cntl reading"""
        #from pyNastran.converters.cart3d.input_cntl_reader import read_input_cntl
        #input_cntl_filename = os.path.join(MODEL_PATH, '')
        #read_input_cntl(input_cntl_filename, log=None, debug=False)

    def test_cart3d_input_c3d(self):
        """tests the input.c3d reading"""
        log = get_logger(level='warning', encoding='utf-8')
        input_c3d_filename = os.path.join(MODEL_PATH, 'input.c3d')
        read_input_c3d(input_c3d_filename, log=log, debug=False, stack=True)

    def test_bjet(self):
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'business_jet', 'bJet.a.tri')

        model = read_cart3d(geometry_filename, log=log)
        #print(model.points.min(axis=0))
        #print(model.points.max(axis=0))
        model.make_half_model(axis='z', remap_nodes=True)  # keep +z
        #print(model.points.min(axis=0))
        #print(model.points.max(axis=0))
        geometry_filename_half = os.path.join(MODEL_PATH, 'half_bjet.a.tri')
        model.write_cart3d(geometry_filename_half)
        read_cart3d(geometry_filename_half, log=log)

def check_array(points, points2):
    nnodes = points.shape[0]
    msg = ''
    nfailed = 0
    if not array_equal(points, points2):
        for nid in range(nnodes):
            p1 = points[nid]
            p2 = points2[nid]
            abs_sum_delta = sum(abs(p1-p2))
            if not allclose(abs_sum_delta, 0.0, atol=1e-6):
                msg += 'n=%s p1=%s p2=%s diff=%s\nsum(abs(p1-p2))=%s\n' % (
                    nid, str(p1), str(p2), str(p1-p2), abs_sum_delta)
                nfailed += 1
                if nfailed == 10:
                    break
    if msg:
        #print(msg)
        raise RuntimeError(msg)

if __name__ == '__main__':  # pragma: no cover
    import time
    time0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - time0))
