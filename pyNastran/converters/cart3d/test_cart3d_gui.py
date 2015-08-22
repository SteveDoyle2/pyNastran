import os

from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.cart3d.cart3dIO import Cart3dIO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'cart3d', 'models')

import unittest


class Cart3dGUI(Cart3dIO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        Cart3dIO.__init__(self)


class TestCart3dGUI(unittest.TestCase):

    def test_cart3d_geometry(self):
        geometry_filename = os.path.join(model_path, 'threePlugs.a.tri')
        #out_filename = os.path.join(model_path, 'panair.out')
        dirname = None

        test = Cart3dGUI()
        #test.load_nastran_geometry(geometry_filename, None)
        test.load_cart3d_geometry(geometry_filename, dirname)

    def test_cart3d_geometry_02(self):
        geometry_filename = os.path.join(model_path, 'threePlugs.bin.tri')
        dirname = None

        test = Cart3dGUI()
        test.load_cart3d_geometry(geometry_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

