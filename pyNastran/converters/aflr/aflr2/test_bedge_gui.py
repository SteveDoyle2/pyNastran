import os
import warnings
import unittest

import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.converters.aflr.aflr2.aflr2 import read_bedge
from pyNastran.converters.aflr.aflr2.bedge_io import BEdge_IO
from pyNastran.gui.testing_methods import FakeGUIMethods

warnings.simplefilter('always')
np.seterr(all='raise')

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'aflr', 'aflr2')


class BEdge_GUI(BEdge_IO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        BEdge_IO.__init__(self, self)
        self.build_fmts(['bedge'], stop_on_failure=True)


class TestBEdgeGUI(unittest.TestCase):
    def test_bedge_geometry(self):
        """tests the bedge gui"""
        bedge_filename = os.path.join(MODEL_PATH, 'm3.bedge')

        test = BEdge_GUI()
        test.log = get_logger(log=None, level='warning', encoding='utf-8')
        test.on_load_geometry(bedge_filename, geometry_format='bedge', raise_error=True)
        #test.load_bedge_geometry(bedge_filename)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
