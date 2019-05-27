import os
import unittest

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.abaqus.abaqus_io import AbaqusIO
#from pyNastran.converters.abaqus.abaqus import read_abaqus
from cpylog import get_logger
from pyNastran.converters.abaqus.test_unit_abaqus import make_model
PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'abaqus', 'models')


class AbaqusGui(AbaqusIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        AbaqusIO.__init__(self, self)
        self.build_fmts(['abaqus'], stop_on_failure=True)

class TestAbaqusGui(unittest.TestCase):
    def test_abaqus_1(self):
        """simple test"""
        lines = make_model()
        abaqus_filename = os.path.join(MODEL_PATH, 'abaqus.inp')
        with open(abaqus_filename, 'w') as abaqus_file:
            abaqus_file.write('\n'.join(lines))
        log = get_logger(level='error', encoding='utf-8')

        test = AbaqusGui()
        test.log = log
        #test.load_abaqus_geometry(abaqus_filename)
        test.on_load_geometry(abaqus_filename, geometry_format='abaqus', raise_error=True)
        os.remove(abaqus_filename)

    def test_abaqus_2(self):
        """two hex blocks with duplicate node ids"""
        abaqus_filename = os.path.join(MODEL_PATH, 'single_block.inp')
        log = get_logger(level='error', encoding='utf-8')

        test = AbaqusGui()
        test.log = log
        #test.load_abaqus_geometry(abaqus_filename)
        test.on_load_geometry(abaqus_filename, geometry_format='abaqus', raise_error=True)

if __name__ == '__main__':  #  pragma: no cover
    unittest.main()
