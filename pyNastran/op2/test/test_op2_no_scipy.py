import unittest
from pathlib import Path

from cpylog import SimpleLogger
import pyNastran
from pyNastran.op2.op2 import read_op2

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').resolve()


class TestOp2NoScipy(unittest.TestCase):
    def test_op2_solid_bending_skip(self):
        log = SimpleLogger(level='warning')
        op2_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.op2'
        op2_filename_out = MODEL_PATH / 'solid_bending' / 'solid_bending_out.op2'
        model = read_op2(op2_filename, log=log)
        model.write_op2(op2_filename_out)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
