import os
import unittest
import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
from pyNastran.converters.cart3d.cart3d import Cart3D
from pyNastran.converters.nastran.nastran_to_cart3d import (
    nastran_to_cart3d, nastran_to_cart3d_filename)


PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'cart3d', 'models')


class Cart3dGUI(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = Cart3dIO(self)
        self.build_fmts(['cart3d'], stop_on_failure=True)


class TestCart3dGUI(unittest.TestCase):

    def test_cart3d_geometry_01(self):
        """tests the ascii three plugs model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'threePlugs.a.tri')

        test = Cart3dGUI()
        test.log = log
        test.on_load_geometry(geometry_filename, geometry_format='cart3d', raise_error=True)

    def test_cart3d_geometry_02(self):
        """tests the binary three plugs model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')

        test = Cart3dGUI()
        test.log = log
        test.model.load_cart3d_geometry(geometry_filename)

    def test_cart3d_geometry_03(self):
        """tests the business jet model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'business_jet', 'bJet.a.tri')

        test = Cart3dGUI()
        test.log = log
        test.model.load_cart3d_geometry(geometry_filename)

    def test_cart3d_bcs(self):
        """tests the power cube model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'power_cube', 'Components.i.tri')

        test = Cart3dGUI()
        test.log = log
        test.model.load_cart3d_geometry(geometry_filename)

    def test_cart3d_results(self):
        log = get_logger(level='warning', encoding='utf-8')
        lines = (
            "5 3 6\n"
            "0. 0. 0.\n"
            "1. 0. 0.\n"
            "2. 0. 0.\n"
            "1. 1. 0.\n"
            "2. 1. 0.\n"
            "1 4 2\n"
            "2 4 5\n"
            "2 5 3\n"
            "1\n"
            "2\n"
            "3\n"
            "1.\n"
            "1. 1. 1. 1. 1.\n"
            "2.\n"
            "2. 2. 2. 2. 2.\n"
            "3.\n"
            "3. 3. 3. 3. 3.\n"
            "4.\n"
            "4. 4. 4. 4. 4.\n"
            "5.\n"
            "5. 5. 5. 5. 5.\n"
        )
        test_path = os.path.join(PKG_PATH, 'converters', 'cart3d', 'models')
        geometry_filename = os.path.join(test_path, 'flat.tri')
        with open(geometry_filename, 'w') as cart3d_file:
            cart3d_file.write(lines)

        #geometry_filename = os.path.join(model_path, 'threePlugs.a.tri')
        #out_filename = os.path.join(model_path, 'panair.out')

        test = Cart3dGUI()
        test.log = log
        #test.load_nastran_geometry(geometry_filename)
        test.model.load_cart3d_geometry(geometry_filename)
        assert test.model.data_map is not None, test.model.data_map

        data = np.array([1., 2., 3.,])
        key = ('centroid', 'Node')
        out = test.model.data_map[key](data)
        #print(out)

        test.setup_fake_text_actors()
        icase = 0
        icase2 = icase + 1
        while icase2 < len(test.result_cases):
            #test.on_cycle_results(case=icase2, show_msg=True)
            unused_result_name = 'dummy'
            test._set_case(unused_result_name, icase2, explicit=False, cycle=False,
                           skip_click_check=False, min_value=None, max_value=None,
                           is_legend_shown=None, show_msg=True)
            icase2 += 1

        os.remove(geometry_filename)

    def test_nastran_to_cart3d_01(self):
        log = get_logger(level='warning', encoding='utf-8')
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
        bdf_filename = os.path.join(MODEL_PATH, 'test01.bdf')
        cart3d_filename = os.path.join(MODEL_PATH, 'test01.tri')
        cart3d_filename_out = os.path.join(MODEL_PATH, 'test01_out.tri')
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(lines)

        nastran_to_cart3d_filename(bdf_filename, cart3d_filename, log=log)
        model = Cart3D(log=log)
        model.read_cart3d(cart3d_filename)
        model.write_cart3d(cart3d_filename_out)

        test = Cart3dGUI()
        test.log = log
        test.model.load_cart3d_geometry(cart3d_filename)
        os.remove(bdf_filename)
        os.remove(cart3d_filename)
        os.remove(cart3d_filename_out)

    def test_nastran_to_cart3d_02(self):
        log = get_logger(level='warning', encoding='utf-8')
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
        bdf_filename = os.path.join(MODEL_PATH, 'test02.bdf')
        cart3d_filename = os.path.join(MODEL_PATH, 'test02.tri')
        cart3d_filename_out = os.path.join(MODEL_PATH, 'test02_out.tri')
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(lines)

        nastran_to_cart3d_filename(bdf_filename, cart3d_filename, log=log)
        model = Cart3D(log=log)
        model.read_cart3d(cart3d_filename)
        model.write_cart3d(cart3d_filename_out)

        test = Cart3dGUI()
        test.log = log
        test.model.load_cart3d_geometry(cart3d_filename)
        os.remove(bdf_filename)
        os.remove(cart3d_filename)
        os.remove(cart3d_filename_out)

    def test_nastran_to_cart3d_03(self):
        log = get_logger(level='warning', encoding='utf-8')
        lines = (
            'SOL 101\n'
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1,,1.0,0.,0.\n'
            'GRID,52,, 0., 0., 0.\n'
            'GRID,2 ,, 1., 0., 0.\n'
            'GRID,19,, 1., 1., 0.\n'
            'GRID,20,, 1., 0., 0.\n'
            'CTRIA3,10,100,1,2,19\n'
            'CQUAD4,11,100,52,2,19,20\n'
            'PSHELL,100,1000,0.1\n'
            'MAT1,1000,3.0e7,,0.3\n'
            'ENDDATA\n'
        )
        bdf_filename = os.path.join(MODEL_PATH, 'test03.bdf')
        cart3d_filename = os.path.join(MODEL_PATH, 'test03.tri')
        #cart3d_filename_out = os.path.join(MODEL_PATH, 'test03_out.tri')
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(lines)

        bdf = BDF(log=log, debug=False)
        bdf.read_bdf(bdf_filename)
        cart3d = nastran_to_cart3d(bdf, log=log)
        cart3d.write_cart3d(cart3d_filename)

        #model = Cart3D()
        #model.read_cart3d(cart3d_filename)
        #model.write_cart3d(cart3d_filename_out)

        #test = Cart3dGUI()
        #test.load_cart3d_geometry(cart3d_filename)
        os.remove(bdf_filename)
        os.remove(cart3d_filename)
        #os.remove(cart3d_filename_out)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
