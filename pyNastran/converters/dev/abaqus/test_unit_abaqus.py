from __future__ import print_function
import os
import unittest

from pyNastran.converters.dev.abaqus.abaqus import read_abaqus

class TestAbaqus(unittest.TestCase):
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
            '1,2,3',
            '*element, type=cpe4',
            '1,2,3,4',
            '*end part',
            #'*material, name=steel',
            #'42',
        ]
        abaqus_filename = 'test.inp'
        with open(abaqus_filename, 'w') as abaqus_file:
            abaqus_file.write('\n'.join(lines))
        read_abaqus(abaqus_filename, debug=False)
        os.remove(abaqus_filename)

if __name__ == '__main__':  #  pragma: no cover
    unittest.main()

