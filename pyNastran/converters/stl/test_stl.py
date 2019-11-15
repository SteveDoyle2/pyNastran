import os
import unittest
import warnings
import shutil

import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.converters.stl.stl import read_stl
from pyNastran.converters.stl.stl_to_nastran import stl_to_nastran, stl_to_nastran_filename
from pyNastran.converters.stl.stl_to_cart3d import stl_to_cart3d
from pyNastran.converters.format_converter import cmd_line_format_converter

warnings.simplefilter('always')
np.seterr(all='raise')

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'converters', 'stl')


class TestSTL(unittest.TestCase):

    def test_stl_io_01(self):
        lines = (
            'solid  testsphere\n'
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
        log = get_logger(level='warning')
        #stl_filename = os.path.join(TEST_PATH, 'tris.stl')
        stl_filename = 'tris1.stl'
        with open(stl_filename, 'w') as stl_file:
            stl_file.write(lines)

        stl = read_stl(stl_filename, log=log, debug=False)
        stl.get_normals_at_nodes()
        scale = 1.0
        stl.scale_nodes(scale)
        stl.shift_nodes(0., 0., 0.)
        stl.log.info('end of shift')

        axes = 'xy'
        stl.flip_axes(axes, scale)
        stl.flip_axes(axes, scale)  # flip back
        stl.log.info('end of flip')

        #stl = STL(log=None, debug=False)
        #stl.read_stl(stl_filename)
        assert len(stl.nodes) == 6, 'nodes=%s' % len(stl.nodes)
        assert len(stl.elements) == 2, 'nelements=%s' % len(stl.elements)

        xyz = 'y'
        tol = 0.00001
        stl.log.info('mirror')
        stl.create_mirror_model(xyz, tol)
        assert len(stl.nodes) == 12, 'nodes=%s' % len(stl.nodes)
        assert len(stl.elements) == 4, 'nelements=%s' % len(stl.elements)
        stl.log.info('end of mirror')

        stl_filename2 = 'tris2.stl'
        stl_filename3 = 'tris3.stl'
        shutil.copyfile(stl_filename, stl_filename2)
        if os.path.exists(stl_filename3):
            os.remove(stl_filename3)
        argv = ['format_converter', 'stl', 'tris*.stl', 'stl', stl_filename3]
        cmd_line_format_converter(argv=argv, quiet=True)

        os.remove(stl_filename)
        os.remove(stl_filename2)
        os.remove(stl_filename3)

    def test_stl_io_02(self):
        lines = (
            'solid  testsphere\n'
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
        log = get_logger(level='warning')
        stl_filename = os.path.join(TEST_PATH, 'tris.stl')
        stl_out_filename = os.path.join(TEST_PATH, 'tris_out.stl')
        stl_bin_filename = os.path.join(TEST_PATH, 'tris_bin.stl')
        with open(stl_filename, 'w') as stl_file:
            stl_file.write(lines)

        stl = read_stl(stl_filename, log=log, debug=False)
        stl._get_normals_data(stl.elements)
        stl.write_stl(stl_out_filename, is_binary=False)
        stl_out = read_stl(stl_out_filename, log=log, debug=False)

        stl.write_stl(stl_bin_filename, is_binary=True)
        stl.write_stl(stl_bin_filename, is_binary=True, normalize_normal_vectors=True)
        stl_bin = read_stl(stl_bin_filename, log=log, debug=False)

        assert len(stl.nodes) == 6, 'nodes=%s' % len(stl.nodes)
        assert len(stl.elements) == 2, 'nelements=%s' % len(stl.elements)
        assert len(stl_out.elements) == 2, 'nelements=%s' % len(stl_out.elements)
        assert len(stl_bin.elements) == 2, 'nelements=%s' % len(stl_bin.elements)
        os.remove(stl_filename)
        os.remove(stl_out_filename)
        #os.remove(stl_bin_filename)

        #outfile_name = os.path.join(TEST_PATH, 'flat.bin.tri')
        #cnormals = stl.get_normals()
        #nnormals = stl.get_normals_at_nodes(cnormals)

    def test_stl_io_02_nan(self):
        lines = (
            'solid  testsphere\n'
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
            '    facet normal  0.13  0.13 -0.98\n'
            '        outer loop\n'
            '            vertex 0.00000 0.00000 0.00000\n'
            '            vertex 0.00000 0.00000 0.00000\n'
            '            vertex 0.00000 0.00000 0.00000\n'
            '        endloop\n'
            '    endfacet\n'
            'endsolid\n'
        )
        log = get_logger(level='warning')
        stl_filename = os.path.join(TEST_PATH, 'tris.stl')
        stl_out_filename = os.path.join(TEST_PATH, 'tris_out.stl')
        stl_bin_filename = os.path.join(TEST_PATH, 'tris_bin.stl')
        with open(stl_filename, 'w') as stl_file:
            stl_file.write(lines)

        stl = read_stl(stl_filename, log=log, debug=False)
        stl._get_normals_data(stl.elements)
        with self.assertRaises(RuntimeError):
            stl.write_stl(stl_out_filename, is_binary=False)
        stl.write_stl(stl_out_filename, is_binary=False, stop_on_failure=False)
        stl_out = read_stl(stl_out_filename, log=log, debug=False)

        stl.write_stl(stl_bin_filename, is_binary=True)
        stl.write_stl(stl_bin_filename, is_binary=True, normalize_normal_vectors=True)
        stl_bin = read_stl(stl_bin_filename, log=log, debug=False)

        assert len(stl.nodes) == 7, 'nodes=%s' % len(stl.nodes)
        assert len(stl.elements) == 3, 'nelements=%s' % len(stl.elements)
        assert len(stl_out.elements) == 3, 'nelements=%s' % len(stl_out.elements)
        assert len(stl_bin.elements) == 3, 'nelements=%s' % len(stl_bin.elements)
        os.remove(stl_filename)
        os.remove(stl_out_filename)
        os.remove(stl_bin_filename)

        #outfile_name = os.path.join(TEST_PATH, 'flat.bin.tri')
        #cnormals = stl.get_normals()
        #nnormals = stl.get_normals_at_nodes(cnormals)

    def test_stl_to_nastran_01(self):
        log = get_logger(level='warning')
        stl_filename = os.path.join(TEST_PATH, 'sphere.stl')
        bdf_filename_8 = os.path.join(TEST_PATH, 'sphere_8.bdf')
        bdf_filename_8b = os.path.join(TEST_PATH, 'sphere_8b.bdf')
        bdf_filename_16 = os.path.join(TEST_PATH, 'sphere_16.bdf')
        bdf_filename_double = os.path.join(TEST_PATH, 'sphere_double.bdf')
        stl_to_nastran_filename(stl_filename, bdf_filename_8, log=log)
        stl_to_nastran(stl_filename, bdf_filename_16, size=16, log=log)
        stl_to_nastran(stl_filename, bdf_filename_double, size=16, is_double=True, log=log)

        argv = ['format_converter', 'stl', stl_filename,
                'nastran', bdf_filename_8b]
        cmd_line_format_converter(argv=argv, quiet=True)

        os.remove(bdf_filename_8)
        os.remove(bdf_filename_8b)
        os.remove(bdf_filename_16)
        os.remove(bdf_filename_double)

    def test_stl_to_cart3d_01(self):
        log = get_logger(level='warning')
        stl_filename = os.path.join(TEST_PATH, 'sphere.stl')
        cart3d_filename = os.path.join(TEST_PATH, 'sphere.tri')
        unused_model = stl_to_cart3d(stl_filename, cart3d_filename=None, log=log)
        argv = ['format_converter', 'stl', stl_filename,
                'cart3d', cart3d_filename]
        cmd_line_format_converter(argv=argv, quiet=True)
        os.remove(cart3d_filename)

def main():  # pragma: no cover
    import time
    time0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - time0))

if __name__ == '__main__':  # pragma: no cover
    main()
