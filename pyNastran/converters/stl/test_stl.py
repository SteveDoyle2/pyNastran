from six.moves import range
import os
from numpy import array_equal, allclose
import unittest

import pyNastran
from pyNastran.converters.stl.stl import STL, read_stl
from pyNastran.converters.stl.stl_to_nastran import stl_to_nastran_filename

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'converters', 'stl')

class TestSTL(unittest.TestCase):

    def test_stl_io_01(self):
        lines = (
            'solid	testsphere\n'
            '    facet normal -0.13 -0.13 -0.98\n'
            '        outer loop\n'
            '            vertex 1.50000 1.50000 0.00000\n'
            '            vertex 1.50000 1.11177 0.05111\n'
            '            vertex 1.11177 1.50000 0.05111\n'
            '        endloop\n'
            '    endfacet\n'
            '    facet normal  0.13  0.13 -0.98\n'
            '        outer loop\n'
            '            vertex 1.50000 1.50000 0.00000\n'
            '            vertex 1.50000 1.88823 0.05111\n'
            '            vertex 1.88823 1.50000 0.05111\n'
            '        endloop\n'
            '    endfacet\n'
            'endsolid\n'
        )
        stl_filename = os.path.join(test_path, 'tris.stl')
        with open(stl_filename, 'w') as f:
            f.write(lines)

        stl = STL(log=None, debug=False)
        stl.read_stl(stl_filename)
        assert len(stl.nodes) == 6, 'nodes=%s' % len(stl.nodes)
        assert len(stl.elements) == 2, 'nelements=%s' % len(stl.elements)
        os.remove(stl_filename)

    def test_stl_io_02(self):
        lines = (
            'solid	testsphere\n'
            '    facet normal -0.13 -0.13 -0.98\n'
            '        outer loop\n'
            '            vertex 1.50000 1.50000 0.00000\n'
            '            vertex 1.50000 1.11177 0.05111\n'
            '            vertex 1.11177 1.50000 0.05111\n'
            '        endloop\n'
            '    endfacet\n'
            '    facet normal  0.13  0.13 -0.98\n'
            '        outer loop\n'
            '            vertex 1.50000 1.50000 0.00000\n'
            '            vertex 1.50000 1.88823 0.05111\n'
            '            vertex 1.88823 1.50000 0.05111\n'
            '        endloop\n'
            '    endfacet\n'
            'endsolid\n'
        )
        stl_filename = os.path.join(test_path, 'tris.stl')
        stl_out_filename = os.path.join(test_path, 'tris_out.stl')
        stl_bin_filename = os.path.join(test_path, 'tris_bin.stl')
        with open(stl_filename, 'w') as f:
            f.write(lines)

        stl = read_stl(stl_filename, log=None, debug=False)
        stl.write_stl(stl_out_filename, is_binary=False)
        stl_out = read_stl(stl_out_filename, log=None, debug=False)

        stl.write_stl(stl_bin_filename, is_binary=True)
        stl.write_stl(stl_bin_filename, is_binary=True, normalize_normal_vectors=True)
        stl_bin = read_stl(stl_bin_filename, log=None, debug=False)

        assert len(stl.nodes) == 6, 'nodes=%s' % len(stl.nodes)
        assert len(stl.elements) == 2, 'nelements=%s' % len(stl.elements)
        os.remove(stl_filename)
        os.remove(stl_out_filename)
        #os.remove(stl_bin_filename)

        #outfile_name = os.path.join(test_path, 'flat.bin.tri')
        #cnormals = stl.get_normals()
        #nnormals = stl.get_normals_at_nodes(cnormals)

    def test_stl_to_nastran_01(self):
        stl_filename = os.path.join(test_path, 'sphere.stl')
        bdf_filename = os.path.join(test_path, 'sphere.bdf')
        stl_to_nastran_filename(stl_filename, bdf_filename)
        os.remove(bdf_filename)


if __name__ == '__main__':  # pragma: no cover
    import time
    t0 = time.time()
    unittest.main()
    #test_1()
    #test_2()
    #test_3()
    print("dt = %s" % (time.time() - t0))
