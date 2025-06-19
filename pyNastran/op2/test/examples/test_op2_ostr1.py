import unittest
from pathlib import Path
import pyNastran
from pyNastran.op2.op2 import OP2

PKG_PATH = Path(pyNastran.__path__[0])
TEST_PATH = PKG_PATH / 'op2' / 'test' / 'examples' / 'ostr1'
OP2_PATH = TEST_PATH / 'recovery.op2'

class TestReadOSTR1(unittest.TestCase):
    def test_ostr1(self):
        op2 = OP2()
        op2.read_op2(OP2_PATH)
