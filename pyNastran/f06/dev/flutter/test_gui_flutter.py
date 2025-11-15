from pathlib import Path
import unittest

import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.f06.dev.flutter.nastran_utils import (
    get_element_table, get_property_table, get_material_table)
from pyNastran.f06.dev.flutter.action import Action
from pyNastran.f06.dev.flutter.utils import (
    get_plot_flags,
    load_f06_op2, get_png_filename, get_plot_file,
    point_removal_str_to_point_removal,
    get_point_removal_str)

IS_DEV = pyNastran.DEV
if IS_DEV:
    from pyNastran.dev.bdf_vectorized3.bdf import read_bdf

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').absolute()
AERO_PATH = MODEL_PATH / 'aero'


class TestGuiFlutter(unittest.TestCase):
    def test_nastran_utils(self) -> None:
        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero.bdf'
        model = read_bdf(bdf_filename, debug='warning')
        get_element_table(model)
        get_property_table(model)
        get_material_table(model)

    def test_action(self) -> None:
        act = Action('cat', 'dog', show=True)
        str(act)

    def test_load_f06_op2(self) -> None:
        f06_filename = MODEL_PATH / 'aero' / '2_mode_flutter' / '0012_flutter.op2'
        log = SimpleLogger(level='warning')
        in_units = 'si'
        out_units = 'si'
        load_f06_op2(f06_filename, log, in_units, out_units, use_rhoref=False)

        op2_filename = MODEL_PATH / 'bwb' / 'bwb_saero.op2'
        load_f06_op2(op2_filename, log, in_units, out_units, use_rhoref=False)

    @unittest.skipIf(not IS_DEV, 'no flutter-dev')
    def test_plot_flutter_0012_dev(self):
        """tests load_f06_op2 and get_png_filename"""
        log = SimpleLogger(level='warning')
        dirname = AERO_PATH / '2_mode_flutter'
        f06_filename = dirname / '0012_flutter.f06'
        model, responses = load_f06_op2(f06_filename, log,
            in_units='si', out_units='si', use_rhoref=False)
        assert model is None, model
        assert len(responses) == 1, responses

        f06_filename2 = 'cat.f06'
        model, responses = load_f06_op2(
            f06_filename2, log,
            in_units='si', out_units='si', use_rhoref=False)
        assert model is None, model
        assert len(responses) == 0, responses

        f06_filename2 = 'cat.op2'
        model, responses = load_f06_op2(
            f06_filename2, log,
            in_units='si', out_units='si', use_rhoref=False)
        assert model is None, model
        assert len(responses) == 0, responses

    def test_plot_file(self):
        out = get_plot_file()

        base = 'base'
        x_plot_type = 'eas'
        plot_type = 'Vg'
        export_to_png = True
        png_filename0, png_filename = get_png_filename(
            base, x_plot_type, plot_type,
            export_to_png)
        assert png_filename0 == 'base_Vg.png', png_filename0
        assert png_filename == 'base_Vg.png', png_filename
        png_filename0, png_filename = get_png_filename(
            base, x_plot_type, 'x-'+plot_type,
            export_to_png)
        assert png_filename0 == 'base_eas-Vg.png', png_filename0
        assert png_filename == 'base_eas-Vg.png', png_filename

        png_filename0, png_filename = get_png_filename(
            base, 'x-'+x_plot_type, plot_type,
            export_to_png=False)
        assert png_filename0 == 'base_Vg.png', png_filename0
        assert png_filename is None, png_filename

    def test_get_plot_flags(self):
        # X_PLOT_TYPES = ['eas', 'tas', 'rho', 'q', 'mach', 'alt', 'kfreq', 'ikfreq', 'index']
        # PLOT_TYPES = ['x-damp-freq', 'x-damp-kfreq', 'root-locus', 'modal-participation', 'zimmerman']
        flags = get_plot_flags('x-damp-freq', 'eas')
        del flags

    def test_point_removal(self):
        log = SimpleLogger(level='debug')
        msg = '400:410,450:500'
        point_removal = [[400.0, 410.0], [450.0, 500.0]]
        point_removal_list = point_removal_str_to_point_removal(
            msg, log)
        assert np.allclose(point_removal, point_removal_list), point_removal_list

        out = get_point_removal_str(point_removal_list)
        assert out == msg, out

if __name__ == '__main__':
    unittest.main()
