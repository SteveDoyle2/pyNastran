import os
from pathlib import Path
import unittest
import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.avl.avl import read_avl
from pyNastran.converters.avl.avl_io import AVL_IO


PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / 'converters' / 'avl' / 'examples'


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
        geometry_filename_out = os.path.join(MODEL_PATH, 'b737_out.avl')
        model = read_avl(geometry_filename)
        model.write_avl(geometry_filename_out)

        test = AvlGUI()
        test.log = log
        test.on_load_geometry(geometry_filename, geometry_format='avl', stop_on_failure=True)

    def test_avl_geometry_02(self):
        """tests the bd model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'bd.avl')
        geometry_filename_out = os.path.join(MODEL_PATH, 'bd_out.avl')
        model = read_avl(geometry_filename)
        model.write_avl(geometry_filename_out)

        test = AvlGUI()
        test.log = log
        test.on_load_geometry(geometry_filename, geometry_format='avl', stop_on_failure=True)

    def test_avl_geometry_03(self):
        """tests the greff model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'greff.avl')
        geometry_filename_out = os.path.join(MODEL_PATH, 'greff_out.avl')
        model = read_avl(geometry_filename)
        model.write_avl(geometry_filename_out)

        test = AvlGUI()
        test.log = log
        test.on_load_geometry(geometry_filename, geometry_format='avl', stop_on_failure=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
