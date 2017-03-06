from __future__ import print_function
import os
import unittest

import warnings
import numpy as np
warnings.simplefilter('always')
np.seterr(all='raise')

import pyNastran
from pyNastran.converters.aflr2.aflr2 import read_bedge
from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.aflr2.bedge_io import BEdge_IO
from pyNastran.utils.log import get_logger

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'aflr2')


class BEdge_GUI(BEdge_IO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        BEdge_IO.__init__(self)


class TestBEdgeGUI(unittest.TestCase):
    def test_bedge_geometry(self):
        """tests the bedge gui"""
        bedge_filename = os.path.join(model_path, 'm3.bedge')
        dirname = None

        test = BEdge_GUI()
        test.log = get_logger(log=None, level='warning', encoding='utf-8')
        test.load_bedge_geometry(bedge_filename, dirname)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
