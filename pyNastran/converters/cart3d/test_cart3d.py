#from pyNastran.converters.cart3d.cart3d_reader import generic_cart3d_reader
from cart3d_reader import generic_cart3d_reader

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
    
    cart3d = generic_cart3d_reader(infileName, log=None, debug=False)
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


if __name__ == '__main__':
    test_1()
    test_2()