from six.moves import range
import os
from numpy import array_equal, allclose
import unittest

import pyNastran
from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'converters', 'cart3d', 'models')

class TestCart3d(unittest.TestCase):

    def test_1(self):
        lines = """7 6
    0.000000 0.000000 0.000000
    1.000000 0.000000 0.000000
    2.000000 0.000000 0.000000
    1.000000 1.000000 0.000000
    2.000000 1.000000 0.000000
    1.000000 -1.000000 0.000000
    2.000000 -1.000000 0.000000
    1 4 2
    2 4 5
    2 5 3
    2 6 1
    5 6 2
    5 5 2
    1
    2
    3
    2
    4
    6
    """
        infileName = os.path.join(test_path, 'flat_full.tri')
        f = open(infileName, 'wb')
        f.write(lines)
        f.close()

        cart3d = Cart3DReader(log=None, debug=False)
        (points, elements, regions, loads) = cart3d.read_cart3d(infileName)
        assert len(points) == 7, 'npoints=%s' % len(points)
        assert len(elements) == 6, 'nelements=%s' % len(elements)
        assert len(regions) == 6, 'nregions=%s' % len(regions)
        assert len(loads) == 0, 'nloads=%s' % len(loads)

    def test_2(self):
        lines = """5 3 6
    0. 0. 0.
    1. 0. 0.
    2. 0. 0.
    1. 1. 0.
    2. 1. 0.
    1 4 2
    2 4 5
    2 5 3
    1
    2
    3
    1.
    1. 1. 1. 1. 1.
    2.
    2. 2. 2. 2. 2.
    3.
    3. 3. 3. 3. 3.
    4.
    4. 4. 4. 4. 4.
    5.
    5. 5. 5. 5. 5.

    """
        infileName = os.path.join(test_path, 'flat.tri')
        f = open(infileName, 'wb')
        f.write(lines)
        f.close()

        cart3d = Cart3DReader(log=None, debug=False)
        (points, elements, regions, loads) = cart3d.read_cart3d(infileName)

        assert len(points) == 5, 'npoints=%s' % len(points)
        assert len(elements) == 3, 'nelements=%s' % len(elements)
        assert len(regions) == 3, 'nregions=%s' % len(regions)

        assert len(loads) == 10, 'nloads=%s' % len(loads)
        assert len(loads['Cp']) == 5, 'nCp=%s' % len(loads['Cp'])

        outfileName = os.path.join(test_path, 'flat.bin.tri')
        cart3d.write_cart3d(outfileName, points, elements, regions, loads=None, is_binary=True)
        cnormals = cart3d.get_normals(points, elements)
        nnormals = cart3d.get_normals_at_nodes(points, elements, cnormals)


    def test_3(self):
        infileName = os.path.join(test_path, 'threePlugs.bin.tri')
        outfileName = os.path.join(test_path, 'threePlugs_out.tri')
        outfileName_bin = os.path.join(test_path, 'threePlugs_bin2.tri')
        outfileName_bin_out = os.path.join(test_path, 'threePlugs_bin_out.tri')
        cart3d = Cart3DReader(log=None, debug=False)

        (points, elements, regions, loads) = cart3d.read_cart3d(infileName)
        cart3d.write_cart3d(outfileName, points, elements, regions, loads=None, is_binary=False)
        cart3d.write_cart3d(outfileName_bin, points, elements, regions, loads=None, is_binary=True)

        (points2, elements2, regions2, loads2) = cart3d.read_cart3d(outfileName)
        check_array(points, points2)

        (points2, elements2, regions2, loads2) = cart3d.read_cart3d(outfileName_bin)
        check_array(points, points2)

        os.remove(outfileName)
        os.remove(outfileName_bin)

        cart3d.write_cart3d(outfileName_bin_out, points2, elements2, regions2, loads2, is_binary=False)
        os.remove(outfileName_bin_out)


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