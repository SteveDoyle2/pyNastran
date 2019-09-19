import os
import unittest
import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.avl.avl_io import AVL_IO


PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'avl', 'examples')


class AvlGUI(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = AVL_IO(self)
        self.build_fmts(['avl'], stop_on_failure=True)


class TestAvlGUI(unittest.TestCase):

    def test_avl_geometry_01(self):
        """tests the 737 model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'b737.avl')

        test = AvlGUI()
        test.log = log
        test.on_load_geometry(geometry_filename, geometry_format='avl', raise_error=True)

    def test_avl_geometry_02(self):
        """tests the bd model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'bd.avl')

        test = AvlGUI()
        test.log = log
        test.on_load_geometry(geometry_filename, geometry_format='avl', raise_error=True)

    def test_avl_geometry_03(self):
        """tests the greff model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'greff.avl')

        test = AvlGUI()
        test.log = log
        test.on_load_geometry(geometry_filename, geometry_format='avl', raise_error=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
