#from pyNastran.converters.cart3d.cart3d_reader import generic_cart3d_reader
from six.moves import range
from cart3d_reader import Cart3DReader
from numpy import array_equal, allclose

def test_1():
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
    infileName = 'flat_full.tri'
    f = open(infileName, 'wb')
    f.write(lines)
    f.close()

    cart3d = Cart3DReader(log=None, debug=False)
    (points, elements, regions, loads) = cart3d.read_cart3d(infileName)
    assert len(points) == 7, 'npoints=%s' % len(points)
    assert len(elements) == 6, 'nelements=%s' % len(elements)
    assert len(regions) == 6, 'nregions=%s' % len(regions)
    assert len(loads) == 0, 'nloads=%s' % len(loads)

def test_2():
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
    infileName = 'flat.tri'
    f = open(infileName, 'wb')
    f.write(lines)
    f.close()

    cart3d = generic_cart3d_reader(infileName, log=None, debug=False)
    (points, elements, regions, loads) = cart3d.read_cart3d(infileName)

    assert len(points) == 5, 'npoints=%s' % len(points)
    assert len(elements) == 3, 'nelements=%s' % len(elements)
    assert len(regions) == 3, 'nregions=%s' % len(regions)

    assert len(loads) == 10, 'nloads=%s' % len(loads)
    assert len(loads['Cp']) == 5, 'nCp=%s' % len(loads['Cp'])

    outfileName = 'flat.bin.tri'
    cart3d.write_cart3d(outfileName, points, elements, regions, loads=None, is_binary=True)

def test_3():
    infileName = 'threePlugs.tri'
    outfileName = 'threePlugs_out.tri'
    outfileName_bin = 'threePlugs_bin.tri'
    outfileName_bin_out = 'threePlugs_bin_out.tri'
    cart3d = Cart3DReader(log=None, debug=False)

    (points, elements, regions, loads) = cart3d.read_cart3d(infileName)
    cart3d.write_cart3d(outfileName, points, elements, regions, loads=None, is_binary=False)
    cart3d.write_cart3d(outfileName_bin, points, elements, regions, loads=None, is_binary=True)

    (points2, elements2, regions2, loads2) = cart3d.read_cart3d(outfileName)
    check(points, points2)

    (points2, elements2, regions2, loads2) = cart3d.read_cart3d(outfileName_bin)
    check(points, points2)

    cart3d.write_cart3d(outfileName_bin_out, points2, elements2, regions2, loads2, is_binary=False)


def check(points, points2):
    nnodes, three = points.shape
    if not array_equal(points, points2):
        for nid in range(nnodes):
            p1 = points[nid]
            p2 = points2[nid]
            abs_sum_delta = sum(abs(p1-p2))
            assert allclose(abs_sum_delta, 0.0, atol=1e-6), 'n=%s p1=%s p2=%s diff=%s\nsum(abs(p1-p2))=%s' % (nid, str(p1), str(p2), str(p1-p2), abs_sum_delta)


if __name__ == '__main__':  # pragma: no cover
    import time
    t0 = time.time()
    #test_1()
    #test_2()
    test_3()
    print("dt = %s" % (time.time() - t0))