import os
import unittest
from pathlib import Path
import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
#from pyNastran.bdf.bdf import BDF
from pyNastran.converters.fluent.fluent_io import FluentIO
from pyNastran.converters.fluent.nastran_to_fluent import nastran_to_fluent


PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH /  '..' / 'models'
BWB_PATH = MODEL_PATH / 'bwb'
TEST_PATH = PKG_PATH / 'converters' / 'fluent'


class FluentGui(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = FluentIO(self)
        self.build_fmts(['fluent'], stop_on_failure=True)


class TestFluentGui(unittest.TestCase):

    def test_fluent_geometry_01(self):
        """tests the ascii three plugs model"""
        log = get_logger(level='warning', encoding='utf-8')
        #geometry_filename = MODEL_PATH / 'threePlugs.a.tri'

        nastran_filename = BWB_PATH / 'bwb_saero.bdf'
        vrt_filename2 = BWB_PATH / 'bwb_saero2.vrt'
        vrt_filename = BWB_PATH / 'bwb_saero.vrt'
        cel_filename = BWB_PATH / 'bwb_saero.cel'
        daten_filename = BWB_PATH / 'bwb_saero.daten'
        tecplot_filename = BWB_PATH / 'bwb_saero.plt'
        nastran_to_fluent(nastran_filename, vrt_filename, log=log)

        test = FluentGui()
        test.log = log
        test.on_load_geometry(vrt_filename, geometry_format='fluent', stop_on_failure=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
