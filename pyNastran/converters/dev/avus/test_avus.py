"""tests non-gui related Avus class/interface"""
import os
import unittest
from cpylog import get_logger

import pyNastran
from pyNastran.converters.dev.avus.avus_grid import read_avus

PKG_PATH = pyNastran.__path__[0]
test_path = os.path.join(PKG_PATH, 'converters', 'dev', 'avus')

class TestAvus(unittest.TestCase):

    def test_avus_io_01(self):
        """geometry"""
        ndim = 3
        npatches = -1
        nzones = 1
        npoints = 7
        nfaces = 6
        ncells = 6
        mxppfs, mxfpcs = -1, -1

        line0 = '%s %s %s\n' % (ndim, nzones, npatches)
        line1 = "%s %s %s %s %s\n" % (npoints, nfaces, ncells, mxppfs, mxfpcs)
        lines = (
            line0 +
            line1 +
            "0.000000 0.000000 0.000000\n"
            "1.000000 0.000000 0.000000\n"
            "2.000000 0.000000 0.000000\n"
            "1.000000 1.000000 0.000000\n"
            "2.000000 1.000000 0.000000\n"
            "1.000000 -1.000000 0.000000\n"
            "2.000000 -1.000000 0.000000\n"
            "3  1 4 2 5 -42\n"
            "3  2 4 5 6 -42\n"
            "3  2 5 3 6 -42\n"
            "3  2 6 1 5 -42\n"
            "3  5 6 2 4 -42\n"
            "3  5 5 2 6 -42\n"
        )
        infile_name = os.path.join(test_path, 'flat_full.tri')
        with open(infile_name, 'w') as avus_file:
            avus_file.write(lines)

        log = get_logger(level='warning', encoding='utf-8')
        model = read_avus(infile_name, log=log, debug=False)
        assert len(model.nodes) == 7, 'nnodes=%s' % len(model.nodes)
        assert len(model.tri_elements) == 6, 'nelements=%s' % len(model.tri_elements)
        #assert len(model.regions) == 6, 'nregions=%s' % len(model.regions)
        #assert len(model.loads) == 0, 'nloads=%s' % len(model.loads)
        os.remove(infile_name)

    def test_avus_io_02(self):
        """geometry + results"""
        ndm = - 1
        nzones = 1
        npatches = -1
        npoints, nfaces, ncells, mxppfs, mxfpcs = 5, 3, 3, -1, -1
        line0 = '%s %s %s\n' % (ndm, nzones, npatches)
        line1 = "%s %s %s %s %s\n" % (npoints, nfaces, ncells, mxppfs, mxfpcs)
        lines = (
            line0 +
            line1 +
            # nodes
            "0. 0. 0.\n"
            "1. 0. 0.\n"
            "2. 0. 0.\n"
            "1. 1. 0.\n"
            "2. 1. 0.\n"

            # faces?
            "3  1 4 2 11 -12\n"
            "3  2 4 5 11 -12\n"
            "3  2 5 3 11 -12\n"
        )
        avus_filename = os.path.join(test_path, 'flat.tri')
        with open(avus_filename, 'w') as avus_file:
            avus_file.write(lines)

        log = get_logger(level='warning', encoding='utf-8')
        model = read_avus(avus_filename, log=log, debug=False)

        assert len(model.nodes) == 5, 'nnodes=%s' % len(model.nodes)
        assert len(model.tri_elements) == 3, 'nelements=%s' % len(model.tri_elements)
        #assert len(model.regions) == 3, 'nregions=%s' % len(model.regions)

        #assert len(model.loads) == 14, 'nloads=%s' % len(model.loads)  # was 10
        #assert len(model.loads['Cp']) == 5, 'nCp=%s' % len(model.loads['Cp'])

        outfile_name = os.path.join(test_path, 'flat.bin.tri')
        model.loads = None
        model.write_avus(outfile_name)
        os.remove(avus_filename)
        os.remove(outfile_name)


if __name__ == '__main__':  # pragma: no cover
    import time
    time0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - time0))
