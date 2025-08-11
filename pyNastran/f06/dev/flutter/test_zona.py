from pathlib import Path
import unittest
import pyNastran
from pyNastran.f06.dev.flutter.read_zona_out import read_zona_out
from pyNastran.f06.dev.flutter.read_zona_aic import read_zona_aic

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_DIR = PKG_PATH / 'bdf' / 'cards' / 'aero' / 'examples' / 'flutter'

try:
    import matplotlib
    IS_MATPLOTLIB = True
    matplotlib.use('Agg')
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

try:
    import pyNastran.f06.dev.flutter.nastran_utils
except ImportError:
    pass

try:
    import pyNastran.f06.dev.flutter.gui_flutter
except ImportError:
    pass

try:
    import pyNastran.f06.dev.flutter.utils
except ImportError:
    pass

try:
    import pyNastran.f06.dev.flutter.actions_builder
except ImportError:
    pass

try:
    import pyNastran.f06.dev.flutter.scalar_bar
except ImportError:
    pass


class TestZona(unittest.TestCase):
    def test_zona_case1_out(self):
        zona_out_filename = MODEL_DIR / 'case1' / 'ha145e.out'
        aic_filename = MODEL_DIR / 'case1' / 'HA145E_AIC.45'
        read_zona_aic(aic_filename)
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

    def test_zona_case4_f06(self):
        from pyNastran.f06.parse_geom import parse_f06_geom
        f06_filename = MODEL_DIR / 'case4' / 'ha145g.f06'
        system_lines, exec_lines, case_lines, bulk_lines = parse_f06_geom(f06_filename)

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
