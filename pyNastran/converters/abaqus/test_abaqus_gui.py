from __future__ import print_function
import os
import unittest

from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.abaqus.abaqus_io import AbaqusIO
#from pyNastran.converters.abaqus.abaqus import read_abaqus
from pyNastran.utils.log import get_logger
from pyNastran.converters.abaqus.test_unit_abaqus import make_model


class AbaqusGui(AbaqusIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        AbaqusIO.__init__(self, self)
        self.build_fmts(['abaqus'], stop_on_failure=True)

class TestAbaqusGui(unittest.TestCase):
    def test_abaqus_1(self):
        """simple test"""
        lines = make_model()
        abaqus_filename = 'test.inp'
        with open(abaqus_filename, 'w') as abaqus_file:
            abaqus_file.write('\n'.join(lines))
        log = get_logger(level='warning', encoding='utf-8')

        test = AbaqusGui()
        test.log = log
        #test.load_abaqus_geometry(abaqus_filename)
        test.on_load_geometry(abaqus_filename, geometry_format='abaqus', raise_error=True)
        os.remove(abaqus_filename)

if __name__ == '__main__':  #  pragma: no cover
    unittest.main()
