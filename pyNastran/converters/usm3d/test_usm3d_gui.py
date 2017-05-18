import os
import unittest

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.usm3d.usm3d_io import Usm3dIO
#from pyNastran.converters.cart3d.cart3d import Cart3D
#from pyNastran.converters.nastran.nastran_to_cart3d import (
    #nastran_to_cart3d, nastran_to_cart3d_filename)
from pyNastran.utils.log import get_logger


pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'usm3d', 'box')


class Usm3dGUI(Usm3dIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        Usm3dIO.__init__(self)


class TestUsm3dGUI(unittest.TestCase):

    def test_cart3d_geometry_01(self):
        """tests the ascii three plugs model"""
        #log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(model_path, 'box.cogsg')
        flo_filename = os.path.join(model_path, 'box.flo')
        dirname = None

        test = Usm3dGUI()
        #test.log = log
        test.load_usm3d_geometry(geometry_filename, dirname)
        test.load_usm3d_results(flo_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

