import os
import unittest
from pathlib import Path
import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
#from pyNastran.bdf.bdf import BDF
from pyNastran.converters.fluent.fluent import Fluent
from pyNastran.converters.fluent.fluent_io import FluentIO
from pyNastran.converters.fluent.nastran_to_fluent import nastran_to_fluent
from pyNastran.converters.fluent.ugrid_to_fluent import ugrid_to_fluent_filename


PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH /  '..' / 'models'
BWB_PATH = MODEL_PATH / 'bwb'
TEST_PATH = PKG_PATH / 'converters' / 'fluent'
UGRID_PATH = PKG_PATH / 'converters' / 'aflr' / 'ugrid' / 'models'


class FluentGui(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = FluentIO(self)
        self.build_fmts(['fluent'], stop_on_failure=True)


class TestFluentGui(unittest.TestCase):

    def test_fluent_geometry_01(self):
        """tests the bwb model"""
        log = get_logger(level='warning', encoding='utf-8')
        #geometry_filename = MODEL_PATH / 'threePlugs.a.tri'

        nastran_filename = BWB_PATH / 'bwb_saero.bdf'
        #vrt_filename2 = BWB_PATH / 'bwb_saero2.vrt'
        vrt_filename = BWB_PATH / 'bwb_saero.vrt'
        #cel_filename = BWB_PATH / 'bwb_saero.cel'
        #daten_filename = BWB_PATH / 'bwb_saero.daten'
        tecplot_filename = BWB_PATH / 'bwb_saero.plt'
        nastran_to_fluent(nastran_filename, vrt_filename, log=log)

        log = get_logger(level='warning', encoding='utf-8')
        test = FluentGui()
        test.log = log
        test.on_load_geometry(
            vrt_filename, geometry_format='fluent',
            stop_on_failure=True)

    def test_fluent_gui_ugrid3d_gui_box(self):
        """simple UGRID3D box model"""
        ugrid_filename = UGRID_PATH / 'box.b8.ugrid'
        fluent_filename = UGRID_PATH / 'box.vrt'
        fluent_model = ugrid_to_fluent_filename(ugrid_filename, fluent_filename)

        log = get_logger(level='warning', encoding='utf-8')
        test = FluentGui()
        test.log = log
        test.on_load_geometry(
            fluent_filename, geometry_format='fluent', stop_on_failure=True)

    def test_fluent_gui_missing_nodes(self):
        model = Fluent()
        model.node_id = np.array([1, 2, 3, 4])
        model.xyz = np.array([
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 1., 0.],
            [0., 1., 0.],
        ])
        model.tris = np.array([
            [1, 10, 1, 2, 3],
        ])
        model.quads = np.array([
            [2, 12, 1, 2, 3, 4],
        ])
        model.element_id = np.array([1, 2])
        model.element_ids = np.array([1, 2])
        model.titles = ['ShellID', 'Pi']
        model.results = np.ones((len(model.element_id), 1)) * 3.14

        log = get_logger(level='warning', encoding='utf-8')
        test = FluentGui()
        test.log = log
        test.model.load_fluent_geometry(model)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
