import os
from pathlib import Path
import unittest
import numpy as np

import pyNastran
from pyNastran.f06.parse_geom import parse_f06_geom
from pyNastran.f06.dev.flutter.read_zaero_out import read_zaero_out
from pyNastran.f06.dev.flutter.read_zaero_aic import read_zaero_aic, write_zaero_aic
from pyNastran.f06.dev.flutter.op2_to_zaero import nastran_to_zaero_f06

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_DIR = PKG_PATH / "bdf" / "cards" / "aero" / "examples" / "flutter"
BDF_MODEL_DIR = PKG_PATH / '..' / 'models'

try:
    import matplotlib
    IS_MATPLOTLIB = True
    matplotlib.use("Agg")
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

try:
    import pyNastran.f06.dev.flutter.gui.gui_flutter
except ImportError:  # pragma: no cover
    pass

try:
    import pyNastran.f06.dev.flutter.gui.actions_builder
except ImportError:
    pass


class TestZaero(unittest.TestCase):
    def test_zaero_case1_aic(self):
        """Round-trip read/write of ZAERO AIC binary file.

        Checks: Mach, reduced frequencies (12 values), two symmetry
        conditions (symmetric + antisymmetric) each with shape
        (12, 24, 24) complex AIC matrices. Writer output re-reads
        with zero numerical error.
        """
        aic_filename = MODEL_DIR / "case1" / "HA145E_AIC.45"
        result = read_zaero_aic(aic_filename)

        assert result.mach == 0.45
        assert len(result.reduced_freqs) == 12
        np.testing.assert_allclose(result.reduced_freqs[0], 0.0)
        np.testing.assert_allclose(result.reduced_freqs[-1], 0.5)
        assert list(result.aic_matrices.keys()) == ["symmetric", "antisymmetric"]
        for key, aic in result.aic_matrices.items():
            assert aic.shape == (12, 24, 24), f"{key}: {aic.shape}"
            assert aic.dtype == np.complex128
        assert result.sym_flag == 2
        assert result.nrows == 24
        assert result.ncols == 24

        tmp_filename = MODEL_DIR / "case1" / "tmp_aic_roundtrip.45"
        try:
            write_zaero_aic(tmp_filename, result)
            result2 = read_zaero_aic(tmp_filename)

            assert result2.mach == result.mach
            np.testing.assert_array_equal(result2.reduced_freqs, result.reduced_freqs)
            for key in result.aic_matrices:
                np.testing.assert_array_equal(
                    result2.aic_matrices[key], result.aic_matrices[key]
                )
        finally:
            if os.path.exists(tmp_filename):
                os.remove(tmp_filename)

    def test_zaero_case1_out(self):
        zaero_out_filename = MODEL_DIR / "case1" / "ha145e.out"
        responses, data_dict = read_zaero_out(zaero_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(data_dict) == 2, data_dict

        matrices = data_dict["matrices"]
        assert "freq" in matrices, matrices
        assert "MHH" in matrices, matrices
        assert "BHH" in matrices, matrices
        assert "KHH" in matrices, matrices
        assert len(matrices["freq"]) == 4
        assert matrices["MHH"].shape == (4, 4)
        assert matrices["BHH"].shape == (4, 4)
        assert matrices["KHH"].shape == (4, 4)
        # print(matrices)

        bdf_model = data_dict["model"]
        assert bdf_model is not None
        assert bdf_model.is_zaero

    def test_zaero_case2_out(self):
        zaero_out_filename = MODEL_DIR / "case2" / "crop.out"
        responses, data_dict = read_zaero_out(zaero_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(data_dict) == 2, data_dict
        assert isinstance(data_dict["matrices"], dict), data_dict

    def test_zaero_case3_out(self):
        zaero_out_filename = MODEL_DIR / "case3" / "ha145fb.out"
        responses, data_dict = read_zaero_out(zaero_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(data_dict) == 2, data_dict
        assert isinstance(data_dict["matrices"], dict), data_dict

    def test_zaero_case4_f06(self):
        f06_filename = MODEL_DIR / "case4" / "ha145g.f06"
        system_lines, exec_lines, case_lines, bulk_lines = parse_f06_geom(f06_filename)

    def test_zaero_case4_out(self):
        zaero_out_filename = MODEL_DIR / "case4" / "ha145g.out"
        responses, data_dict = read_zaero_out(zaero_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(data_dict) == 2, data_dict
        assert isinstance(data_dict["matrices"], dict), data_dict

    def test_zaero_case5_out(self):
        zaero_out_filename = MODEL_DIR / "case5" / "f16ma41.out"
        responses, data_dict = read_zaero_out(zaero_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(data_dict) == 2, data_dict
        assert isinstance(data_dict["matrices"], dict), data_dict

    def test_zaero_case6_out_trim(self):
        zaero_out_filename = MODEL_DIR / 'case6' / 'agard_trim.out'
        responses, data_dict = read_zaero_out(zaero_out_filename)
        assert len(responses) == 0, responses
        assert len(data_dict) == 3, data_dict
        assert isinstance(data_dict["matrices"], dict), data_dict
        assert data_dict["model"] is not None

        trim = data_dict["trim"]
        assert trim["subcase"] == 1
        assert trim["mach"] == 0.954
        assert trim["dynamic_pressure"] == 1200.0
        assert "ALPHA" in trim["stability_derivatives"], trim["stability_derivatives"]
        assert "NZ" in trim["trim_results"], trim["trim_results"]
        assert "ALPHA" in trim["trim_variables"], trim["trim_variables"]
        assert len(trim["aero_forces"]) > 0, trim["aero_forces"]

    def test_zaero_case6_out_tran(self):
        zaero_out_filename = MODEL_DIR / "case6" / "agardztran.out"
        responses, data_dict = read_zaero_out(zaero_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(data_dict) == 2, data_dict
        assert isinstance(data_dict["matrices"], dict), data_dict

    def test_zaero_case7_out(self):
        zaero_out_filename = MODEL_DIR / "case7" / "agardztaw.out"
        responses, data_dict = read_zaero_out(zaero_out_filename, debug=None)
        assert len(responses) == 1, responses
        assert len(data_dict) == 2, data_dict
        assert isinstance(data_dict["matrices"], dict), data_dict

    def test_nastran_to_zaero(self):
        op2_filename = BDF_MODEL_DIR / "plate_py" / "plate_py.op2"
        bdf_filename = BDF_MODEL_DIR / "plate_py" / "plate_py.dat"
        nastran_to_zaero_f06(op2_filename, bdf_filename=bdf_filename)
        # nastran_to_zaero_f06(op2_filename, bdf_filename=None)

if __name__ == "__main__":  # pragma: no cover
    unittest.main()
