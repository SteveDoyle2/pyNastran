import os
import unittest
import warnings
import shutil
from pathlib import Path

import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.converters.fld.fld import read_fld

warnings.simplefilter('always')
np.seterr(all='raise')

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = Path(os.path.join(PKG_PATH, 'converters', 'fld'))


class TestFLD(unittest.TestCase):

    def test_fld_io_01(self):
        log = get_logger(level='warning')
        fld_filename = TEST_PATH / 'test.fld'

        model = read_fld(fld_filename, log=log, debug=False)
        assert len(model.xyzp) == 3, model.xyzp.shape

def main():  # pragma: no cover
    import time
    time0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - time0))

if __name__ == '__main__':  # pragma: no cover
    main()
