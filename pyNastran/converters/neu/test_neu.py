from pathlib import Path
import unittest
import pyNastran

from pyNastran.utils import print_bad_path
from pyNastran.converters.neu.neu import read_neu

PKG_PATH = Path(pyNastran.__path__[0])
assert PKG_PATH.exists(), print_bad_path(PKG_PATH)
DIRNAME = PKG_PATH / r'bdf' / 'cards' / 'aero' / 'examples' / 'flutter' / 'case5'
CURDIR = Path(__file__).parent

class TestNeu(unittest.TestCase):
    def test_all_8p2(self):
        neu_filename = CURDIR / 'all_elem_test_giv_eb.neu'
        model = read_neu(neu_filename, debug=None)
    def test_cp_anti(self):
        neu_filename = DIRNAME / 'cp2anti.neu'
        model = read_neu(neu_filename, debug=None)
    def test_flut_anti(self):
        neu_filename = DIRNAME / 'flut_anti.neu'
        model = read_neu(neu_filename)
    def test_f16_aero(self):
        neu_filename = DIRNAME / 'f16-aero.neu'
        model = read_neu(neu_filename)


if __name__ == '__main__':
    unittest.main()
