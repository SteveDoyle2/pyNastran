import os
import unittest

import pyNastran
from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
from pyNastran.converters.cart3d.cart3d import Cart3D
from pyNastran.converters.nastran.nastran_to_cart3d import nastran_to_cart3d, nastran_to_cart3d_filename


pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'cart3d', 'models')


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

    def test_nastran_to_cart3d_01(self):
        lines = (
            'SOL 101\n'
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1,,1.0,0.,0.\n'
            'GRID,2,,1.0,1.,0.\n'
            'GRID,3,,0.0,1.,0.\n'
            'CTRIA3,10,100,1,2,3\n'
            'PSHELL,100,1000,0.1\n'
            'MAT1,1000,3.0e7,,0.3\n'
            'ENDDATA\n'
        )
        bdf_filename = os.path.join(model_path, 'test01.bdf')
        cart3d_filename = os.path.join(model_path, 'test01.tri')
        cart3d_filename_out = os.path.join(model_path, 'test01_out.tri')
        with open(bdf_filename, 'w') as f:
            f.write(lines)

        nastran_to_cart3d_filename(bdf_filename, cart3d_filename)
        model = Cart3D()
        model.read_cart3d(cart3d_filename)
        model.write_cart3d(cart3d_filename_out)

        test = Cart3dGUI()
        test.load_cart3d_geometry(cart3d_filename, dirname=None)

    def test_nastran_to_cart3d_02(self):
        lines = (
            'SOL 101\n'
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1,,1.0,0.,0.\n'
            'GRID,2,,1.0,1.,0.\n'
            'GRID,19,,0.0,1.,0.\n'
            'CTRIA3,10,100,1,2,19\n'
            'PSHELL,100,1000,0.1\n'
            'MAT1,1000,3.0e7,,0.3\n'
            'ENDDATA\n'
        )
        bdf_filename = os.path.join(model_path, 'test02.bdf')
        cart3d_filename = os.path.join(model_path, 'test02.tri')
        cart3d_filename_out = os.path.join(model_path, 'test02_out.tri')
        with open(bdf_filename, 'w') as f:
            f.write(lines)

        nastran_to_cart3d_filename(bdf_filename, cart3d_filename)
        model = Cart3D()
        model.read_cart3d(cart3d_filename)
        model.write_cart3d(cart3d_filename_out)

        test = Cart3dGUI()
        test.load_cart3d_geometry(cart3d_filename, dirname=None)

    def test_nastran_to_cart3d_02(self):
        lines = (
            'SOL 101\n'
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,52,,1.0,0.,0.\n'
            'GRID,2,,1.0,1.,0.\n'
            'GRID,19,,0.0,1.,0.\n'
            'CTRIA3,10,100,52,2,19\n'
            'PSHELL,100,1000,0.1\n'
            'MAT1,1000,3.0e7,,0.3\n'
            'ENDDATA\n'
        )
        bdf_filename = os.path.join(model_path, 'test03.bdf')
        cart3d_filename = os.path.join(model_path, 'test03.tri')
        cart3d_filename_out = os.path.join(model_path, 'test03_out.tri')
        with open(bdf_filename, 'w') as f:
            f.write(lines)

        bdf = BDF(debug=False)
        bdf.read_bdf(bdf_filename)
        cart3d = nastran_to_cart3d(bdf)
        cart3d.write_cart3d(cart3d_filename)

        #model = Cart3D()
        #model.read_cart3d(cart3d_filename)
        #model.write_cart3d(cart3d_filename_out)

        #test = Cart3dGUI()
        #test.load_cart3d_geometry(cart3d_filename, dirname=None)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

