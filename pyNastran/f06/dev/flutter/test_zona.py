from pathlib import Path
import unittest
import pyNastran
from pyNastran.f06.dev.flutter.read_zona_out import read_zona_out
PKG_PATH = Path(pyNastran.__path__[0])
MODEL_DIR = PKG_PATH / 'bdf' / 'cards' / 'aero' / 'examples' / 'flutter'


class TestZona(unittest.TestCase):
    def test_zona_case1_out(self):
        zona_out_filename = MODEL_DIR / 'case1' / 'ha145e.out'
        responses, mass = read_zona_out(zona_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(mass) == 0, mass

    def test_zona_case2_out(self):
        zona_out_filename = MODEL_DIR / 'case2' / 'crop.out'
        responses, mass = read_zona_out(zona_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(mass) == 0, mass

    def test_zona_case3_out(self):
        zona_out_filename = MODEL_DIR / 'case3' / 'ha145fb.out'
        responses, mass = read_zona_out(zona_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(mass) == 0, mass

    def test_zona_case4_out(self):
        zona_out_filename = MODEL_DIR / 'case4' / 'ha145g.out'
        responses, mass = read_zona_out(zona_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(mass) == 0, mass

    def test_zona_case5_out(self):
        zona_out_filename = MODEL_DIR / 'case5' / 'f16ma41.out'
        responses, mass = read_zona_out(zona_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(mass) == 0, mass

    # def test_zona_case6_out_trim(self):
    #     zona_out_filename = MODEL_DIR / 'case6' / 'agard_trim.out'
    #     responses, mass = read_zona_out(zona_out_filename)
    #     assert len(responses) == 1, responses
    #     assert len(mass) == 0, mass

    def test_zona_case6_out_tran(self):
        zona_out_filename = MODEL_DIR / 'case6' / 'agardztran.out'
        responses, mass = read_zona_out(zona_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(mass) == 0, mass

    def test_zona_case7_out(self):
        zona_out_filename = MODEL_DIR / 'case7' / 'agardztaw.out'
        responses, mass = read_zona_out(zona_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(mass) == 0, mass


if __name__ == '__main__':
    unittest.main()
