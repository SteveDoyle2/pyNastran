import os
import unittest

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.openfoam.openfoam_io import OpenFoamIO
from pyNastran.converters.cart3d.cart3d import Cart3D
from pyNastran.converters.nastran.nastran_to_cart3d import (
    nastran_to_cart3d, nastran_to_cart3d_filename)
from pyNastran.utils.log import get_logger
from pyNastran.utils import print_bad_path

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'openfoam', 'models')


class OpenFoamGUI(OpenFoamIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        OpenFoamIO.__init__(self)


class TestOpenFoamGUI(unittest.TestCase):

    def test_openfoam_geometry_01(self):
        """tests the ascii three plugs model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(model_path, 'SnakeRiverCanyon', 'system', 'blockMeshDict')
        assert os.path.exists(geometry_filename), print_bad_path(geometry_filename)
        test = OpenFoamGUI()
        test.log = log
        test.load_openfoam_geometry_shell(geometry_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

