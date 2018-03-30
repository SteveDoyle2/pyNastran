"""
tests Usm3d
"""
import os
import unittest

import pyNastran
from pyNastran.converters.usm3d.usm3d_reader import read_usm3d
from pyNastran.utils.log import get_logger


PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'usm3d', 'box')


class TestUsm3d(unittest.TestCase):
    """tests file read/write"""

    def test_usm3d_01(self):
        """tests the box.cogwg/box.flo turbulent model"""
        log = get_logger(level='error', encoding='utf-8')
        base_filename = os.path.join(MODEL_PATH, 'box')
        #flo_filename = os.path.join(MODEL_PATH, 'box.flo')

        model = read_usm3d(base_filename, log=log)
        basename = 'cat'
        model.write_usm3d(basename)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
