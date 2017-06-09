"""tests non-gui related Cart3d class/interface"""
from six.moves import range
import os
from numpy import array_equal, allclose
import unittest

import pyNastran
from pyNastran.converters.dev.avus.avus_grid import AvusGrid, read_avus
#from pyNastran.converters.avus.cart3d_to_nastran import cart3d_to_nastran_filename, cart3d_to_nastran_model
#from pyNastran.converters.avus.cart3d_to_stl import cart3d_to_stl_filename
#from pyNastran.converters.avus.cart3d_to_tecplot import cart3d_to_tecplot
#from pyNastran.converters.avus.input_c3d_reader import read_input_c3d
from pyNastran.utils.log import get_logger

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'converters', 'dev', 'avus')

class TestAvus(unittest.TestCase):

    def test_avus_io_01(self):
        """geometry"""
        zones = 1
        npoints = 7
        nfaces = 6
        ncells = 6
        mxppfs, mxfpcs = 1000, 1000

        line0 ="%s %s %s %s %s %s\n" % (zones, npoints, nfaces, ncells, mxppfs, mxfpcs)
        print(line0)
        lines = (
            line0 +
            "0.000000 0.000000 0.000000\n"
            "1.000000 0.000000 0.000000\n"
            "2.000000 0.000000 0.000000\n"
            "1.000000 1.000000 0.000000\n"
            "2.000000 1.000000 0.000000\n"
            "1.000000 -1.000000 0.000000\n"
            "2.000000 -1.000000 0.000000\n"
            "3  1 4 2\n"
            "3  2 4 5\n"
            "3  2 5 3\n"
            "3  2 6 1\n"
            "3  5 6 2\n"
            "3  5 5 2\n"
        )
        infile_name = os.path.join(test_path, 'flat_full.tri')
        with open(infile_name, 'w') as avus_file:
            avus_file.write(lines)

        log = get_logger(level='warning', encoding='utf-8')
        model = read_avus(infile_name, log=log, debug=False)
        assert len(model.nodes) == 7, 'nnodes=%s' % len(model.nodes)
        assert len(model.elements) == 6, 'nelements=%s' % len(model.elements)
        #assert len(model.regions) == 6, 'nregions=%s' % len(model.regions)
        #assert len(model.loads) == 0, 'nloads=%s' % len(model.loads)
        os.remove(infile_name)

    def test_avus_io_02(self):
        """geometry + results"""
        zones, npoints, nfaces, ncells, mxppfs, mxfpcs = -1, 5, 3, 3, -1, -1
        line0 ="%s %s %s %s %s %s\n" % (zones, npoints, nfaces, ncells, mxppfs, mxfpcs)
        lines = (
            line0 +
            "0. 0. 0.\n"
            "1. 0. 0.\n"
            "2. 0. 0.\n"
            "1. 1. 0.\n"
            "2. 1. 0.\n"
            "3  1 4 2\n"
            "3  2 4 5\n"
            "3  2 5 3\n"
        )
        avus_filename = os.path.join(test_path, 'flat.tri')
        with open(avus_filename, 'w') as avus_file:
            avus_file.write(lines)

        log = get_logger(level='warning', encoding='utf-8')
        model = read_avus(avus_filename, log=log, debug=False)

        assert len(model.nodes) == 5, 'nnodes=%s' % len(model.nodes)
        assert len(model.elements) == 3, 'nelements=%s' % len(model.elements)
        #assert len(model.regions) == 3, 'nregions=%s' % len(model.regions)

        #assert len(model.loads) == 14, 'nloads=%s' % len(model.loads)  # was 10
        #assert len(model.loads['Cp']) == 5, 'nCp=%s' % len(model.loads['Cp'])

        outfile_name = os.path.join(test_path, 'flat.bin.tri')
        model.loads = None
        model.write_avus(outfile_name, is_binary=True)
        os.remove(avus_filename)
        os.remove(outfile_name)

    #def test_avus_io_03(self):
        #"""read/write geometry in ascii/binary"""
        #log = get_logger(level='warning', encoding='utf-8')
        #infile_name = os.path.join(test_path, 'threePlugs.bin.tri')
        #outfile_name = os.path.join(test_path, 'threePlugs_out.tri')
        #outfile_name_bin = os.path.join(test_path, 'threePlugs_bin2.tri')
        #outfile_name_bin_out = os.path.join(test_path, 'threePlugs_bin_out.tri')

        #cart3d = read_cart3d(infile_name, log=log, debug=False)
        #cart3d.write_cart3d(outfile_name, is_binary=False)
        #cart3d.write_cart3d(outfile_name_bin, is_binary=True)

        #cart3d_ascii = read_cart3d(outfile_name, log=log, debug=False)
        #check_array(cart3d.points, cart3d_ascii.points)
        #check_array(cart3d.elements, cart3d_ascii.elements)

        #cart3d_bin = read_cart3d(outfile_name_bin, log=log, debug=False)
        #check_array(cart3d.points, cart3d_bin.points)
        #check_array(cart3d.elements, cart3d_ascii.elements)
        ##print(cart3d_bin.points)

        ##print('---------------')
        ##print(cart3d_bin.points)

        #os.remove(outfile_name)
        #os.remove(outfile_name_bin)

        #cart3d.write_cart3d(outfile_name_bin_out, is_binary=False)
        #os.remove(outfile_name_bin_out)

    #def test_avus_to_stl(self):
        #"""convert to stl"""
        #log = get_logger(level='warning', encoding='utf-8')
        #cart3d_filename = os.path.join(test_path, 'threePlugs.bin.tri')
        #stl_filename = os.path.join(test_path, 'threePlugs.stl')
        #cart3d_to_stl_filename(cart3d_filename, stl_filename, log=log)
        ##os.remove(stl_filename)

    #def test_avus_to_tecplot(self):
        #"""convert to tecplot"""
        #log = get_logger(level='warning', encoding='utf-8')
        #cart3d_filename = os.path.join(test_path, 'threePlugs.bin.tri')
        #tecplot_filename = os.path.join(test_path, 'threePlugs.plt')
        #cart3d_to_tecplot(cart3d_filename, tecplot_filename, log=log)
        ##os.remove(tecplot_filename)

    #def test_avus_to_nastran_01(self):
        #"""convert to nastran small field"""
        #log = get_logger(level='warning', encoding='utf-8')
        #cart3d_filename = os.path.join(test_path, 'threePlugs.bin.tri')
        #bdf_filename = os.path.join(test_path, 'threePlugs.bdf')
        #cart3d_to_nastran_filename(cart3d_filename, bdf_filename, log=log)
        #os.remove(bdf_filename)

    #def test_avus_to_nastran_02(self):
        #"""convert to nastran large field"""
        #log = get_logger(level='warning', encoding='utf-8')
        #cart3d_filename = os.path.join(test_path, 'threePlugs.bin.tri')
        #bdf_filename = os.path.join(test_path, 'threePlugs2.bdf')
        #model = cart3d_to_nastran_model(cart3d_filename, log=log)
        #model.write_bdf(bdf_filename, size=16)
        #self.assertAlmostEqual(model.nodes[1].xyz[0], 1.51971436,
                               #msg='if this is 0.0, then the assign_type method had the float32 check removed')
        #os.remove(bdf_filename)
        #model.write_bdf(out_filename=None, encoding=None, size=8,
                       #is_double=False,
                       #interspersed=True,
                       #enddata=None)

    #def test_cart3d_input_cntl(self):
        #"""tests the input.cntl reading"""
        #from pyNastran.converters.cart3d.input_cntl_reader import read_input_cntl
        #input_cntl_filename = os.path.join(test_path, '')
        #read_input_cntl(input_cntl_filename, log=None, debug=False)

    #def test_avus_input_c3d(self):
        #"""tests the input.c3d reading"""
        #log = get_logger(level='warning', encoding='utf-8')
        #input_c3d_filename = os.path.join(test_path, 'input.c3d')
        #read_input_c3d(input_c3d_filename, log=log, debug=False, stack=True)

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
    t0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - t0))
