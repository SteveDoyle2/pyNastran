import os
import unittest

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.openfoam.openfoam_io import OpenFoamIO
from pyNastran.utils.log import get_logger
from pyNastran.utils import print_bad_path

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'openfoam', 'models')


class OpenFoamGUI(OpenFoamIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        OpenFoamIO.__init__(self)


class TestOpenFoamGUI(unittest.TestCase):

    def test_openfoam_geometry_01(self):
        """tests the ascii three plugs model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'SnakeRiverCanyon', 'system', 'blockMeshDict')
        assert os.path.exists(geometry_filename), print_bad_path(geometry_filename)
        test = OpenFoamGUI()
        test.log = log
        test.load_openfoam_geometry_shell(geometry_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

