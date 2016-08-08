from six.moves import range
import os
from numpy import array_equal, allclose
import unittest

import pyNastran
from pyNastran.converters.cart3d.cart3d import Cart3D
from pyNastran.converters.cart3d.cart3d_to_nastran import cart3d_to_nastran_filename, cart3d_to_nastran_model

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'converters', 'cart3d', 'models')

class TestCart3d(unittest.TestCase):

    def test_cart3d_io_01(self):
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
        infile_name = os.path.join(test_path, 'flat_full.tri')
        with open(infile_name, 'w') as f:
            f.write(lines)

        cart3d = Cart3D(log=None, debug=False)
        cart3d.read_cart3d(infile_name)
        assert len(cart3d.points) == 7, 'npoints=%s' % len(cart3d.points)
        assert len(cart3d.elements) == 6, 'nelements=%s' % len(cart3d.elements)
        assert len(cart3d.regions) == 6, 'nregions=%s' % len(cart3d.regions)
        assert len(cart3d.loads) == 0, 'nloads=%s' % len(cart3d.loads)

    def test_cart3d_io_02(self):
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
        infile_name = os.path.join(test_path, 'flat.tri')
        with open(infile_name, 'w') as f:
            f.write(lines)

        cart3d = Cart3D(log=None, debug=False)
        cart3d.read_cart3d(infile_name)

        assert len(cart3d.points) == 5, 'npoints=%s' % len(cart3d.points)
        assert len(cart3d.elements) == 3, 'nelements=%s' % len(cart3d.elements)
        assert len(cart3d.regions) == 3, 'nregions=%s' % len(cart3d.regions)

        assert len(cart3d.loads) == 14, 'nloads=%s' % len(cart3d.loads)  # was 10
        assert len(cart3d.loads['Cp']) == 5, 'nCp=%s' % len(cart3d.loads['Cp'])

        outfile_name = os.path.join(test_path, 'flat.bin.tri')
        cart3d.loads = None
        cart3d.write_cart3d(outfile_name, is_binary=True)
        cnormals = cart3d.get_normals()
        nnormals = cart3d.get_normals_at_nodes(cnormals)


    def test_cart3d_io_03(self):
        infile_name = os.path.join(test_path, 'threePlugs.bin.tri')
        outfile_name = os.path.join(test_path, 'threePlugs_out.tri')
        outfile_name_bin = os.path.join(test_path, 'threePlugs_bin2.tri')
        outfile_name_bin_out = os.path.join(test_path, 'threePlugs_bin_out.tri')
        cart3d = Cart3D(log=None, debug=False)

        cart3d.read_cart3d(infile_name)
        cart3d.write_cart3d(outfile_name, is_binary=False)
        cart3d.write_cart3d(outfile_name_bin, is_binary=True)

        cart3d_ascii = Cart3D(log=None, debug=False)
        cart3d_ascii.read_cart3d(outfile_name)
        check_array(cart3d.points, cart3d_ascii.points)

        cart3d_bin = Cart3D(log=None, debug=False)
        cart3d_bin.read_cart3d(outfile_name_bin)
        check_array(cart3d.points, cart3d_bin.points)

        os.remove(outfile_name)
        os.remove(outfile_name_bin)

        cart3d.write_cart3d(outfile_name_bin_out, is_binary=False)
        os.remove(outfile_name_bin_out)

    def test_cart3d_to_nastran_01(self):
        cart3d_filename = os.path.join(test_path, 'threePlugs.bin.tri')
        bdf_filename = os.path.join(test_path, 'threePlugs.bdf')
        cart3d_to_nastran_filename(cart3d_filename, bdf_filename)

    def test_cart3d_to_nastran_02(self):
        cart3d_filename = os.path.join(test_path, 'threePlugs.bin.tri')
        bdf_filename = os.path.join(test_path, 'threePlugs2.bdf')
        model = cart3d_to_nastran_model(cart3d_filename)
        model.write_bdf(bdf_filename, size=16)
        self.assertAlmostEqual(model.nodes[1].xyz[0], 1.51971436,
                               msg='if this is 0.0, then the assign_type method had the float32 check removed')
        #model.write_bdf(out_filename=None, encoding=None, size=8,
                       #is_double=False,
                       #interspersed=True,
                       #enddata=None)

def check_array(points, points2):
    nnodes = points.shape[0]
    if not array_equal(points, points2):
        for nid in range(nnodes):
            p1 = points[nid]
            p2 = points2[nid]
            abs_sum_delta = sum(abs(p1-p2))
            assert allclose(abs_sum_delta, 0.0, atol=1e-6), 'n=%s p1=%s p2=%s diff=%s\nsum(abs(p1-p2))=%s' % (nid, str(p1), str(p2), str(p1-p2), abs_sum_delta)


if __name__ == '__main__':  # pragma: no cover
    import time
    t0 = time.time()
    unittest.main()
    #test_1()
    #test_2()
    #test_3()
    print("dt = %s" % (time.time() - t0))
