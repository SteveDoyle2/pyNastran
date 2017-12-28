from __future__ import print_function
import os
import unittest

#import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.abaqus.abaqus_io import AbaqusIO
from pyNastran.converters.abaqus.abaqus import read_abaqus
from pyNastran.utils.log import get_logger


class AbaqusGui(AbaqusIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        AbaqusIO.__init__(self)

class TestAbaqusGui(unittest.TestCase):
    def test_abaqus_1(self):
        """simple test"""
        lines = [
            '*part, name=dummy',
            '*node',
            '1,0.,0.,0.',
            '2,1.,0.,0.',
            '3,1.,1.,0.',
            '4,0.,1.,0.',

            '*element, type=cpe3',
            '1,1,2,3',
            '*element, type=cpe4',
            '2,1,2,3,4',
            '*end part',
            #'*material, name=steel',
            #'42',
        ]
        abaqus_filename = 'test.inp'
        with open(abaqus_filename, 'w') as abaqus_file:
            abaqus_file.write('\n'.join(lines))
        log = get_logger(level='warning', encoding='utf-8')

        test = AbaqusGui()
        test.log = log
        test.load_abaqus_geometry(abaqus_filename)
        os.remove(abaqus_filename)

if __name__ == '__main__':  #  pragma: no cover
    unittest.main()
