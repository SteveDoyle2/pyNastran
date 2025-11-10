from pathlib import Path
import unittest

from cpylog import SimpleLogger

from pyNastran.dev.bdf_vectorized3.bdf import read_bdf
from pyNastran.f06.dev.flutter.nastran_utils import (
    get_element_table, get_property_table, get_material_table)
from pyNastran.f06.dev.flutter.action import Action
from pyNastran.f06.dev.flutter.utils import load_f06_op2

import pyNastran
PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').absolute()


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
        log = SimpleLogger(level='info')
        in_units = 'si'
        out_units = 'si'
        load_f06_op2(f06_filename, log, in_units, out_units, use_rhoref=False)

        op2_filename = MODEL_PATH / 'bwb' / 'bwb_saero.op2'
        load_f06_op2(op2_filename, log, in_units, out_units, use_rhoref=False)


if __name__ == '__main__':
    unittest.main()
