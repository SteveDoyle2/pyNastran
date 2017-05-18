import os
import unittest

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.fast.fast_io import FastIO
#from pyNastran.converters..cart3d import Cart3D
#from pyNastran.converters.nastran.nastran_to_cart3d import nastran_to_cart3d, nastran_to_cart3d_filename


pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'converters', 'fast')


class FastGUI(FastIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        FastIO.__init__(self)


class TestFastGUI(unittest.TestCase):

    def test_cart3d_geometry_01(self):
        geometry_filename = os.path.join(test_path, 'flow_demo1', 'om6inviscid.fgrid')
        dirname = None

        test = FastGUI()
        test.load_fast_geometry(geometry_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

