import os
import unittest
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
#from pyNastran.bdf.bdf import BDF
from pyNastran.converters.dev.avus.avus_io import AvusIO

PKG_PATH = pyNastran.__path__[0]
test_path = os.path.join(PKG_PATH, 'converters', 'dev', 'avus')


class AvusGUI(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = AvusIO(self)


class TestAvusGUI(unittest.TestCase):

    def test_avus_io_01(self):
        """tests geometry"""
        ndim = 3
        zones = 1
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
        avus_filename = os.path.join(test_path, 'flat_full.tri')
        with open(avus_filename, 'w') as avus_file:
            avus_file.write(lines)

        test = AvusGUI()
        test.log = get_logger(level='warning', encoding='utf-8')
        test.model.load_avus_geometry(avus_filename)
        os.remove(avus_filename)

    def test_avus_geometry_01(self):
        """tests geometry"""
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

        test = AvusGUI()
        test.log = get_logger(level='warning', encoding='utf-8')
        test.model.load_avus_geometry(avus_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
