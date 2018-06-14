import os
import unittest

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.fast.fast_io import FastIO


PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'converters', 'fast')


class FastGUI(FastIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        FastIO.__init__(self, self)
        self.build_fmts(['fast'], stop_on_failure=True)


class TestFastGUI(unittest.TestCase):

    def test_fast_geometry_01(self):
        geometry_filename = os.path.join(TEST_PATH, 'flow_demo1', 'om6inviscid.fgrid')

        test = FastGUI()
        #test.load_fast_geometry(geometry_filename)
        test.on_load_geometry(geometry_filename, geometry_format='fast', raise_error=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
