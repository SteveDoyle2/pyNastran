import os
import warnings
from pathlib import Path
import unittest

import numpy as np

from cpylog import get_logger
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.fld.fld_io import FLD_IO
import pyNastran

warnings.simplefilter('always')
np.seterr(all='raise')

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = Path(os.path.join(PKG_PATH, 'converters', 'fld'))


class FLD_GUI(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = FLD_IO(self)
        self.build_fmts(['fld'], stop_on_failure=True)


class TestFldGui(unittest.TestCase):
    def test_fld_geometry(self):
        log = get_logger(level='warning')
        fld_filename = TEST_PATH / 'test.fld'

        test = FLD_GUI()
        test.log = log
        test.on_load_geometry(fld_filename, geometry_format='fld', stop_on_failure=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
