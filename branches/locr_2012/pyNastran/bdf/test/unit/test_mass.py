import unittest

import os
import pyNastran
from pyNastran.bdf.bdf import BDF

rootpath = pyNastran.__path__[0]
testpath = os.path.join(rootpath, 'bdf', 'test', 'unit')

class TestMass(unittest.TestCase):
    def test_mass1(self):  # passes
        mesh = BDF(debug=True, log=None)
        bdfname = os.path.join(testpath, 'test_mass.dat')
        mesh.readBDF(bdfname, includeDir=None, xref=True) 
        self.assertAlmostEqual(mesh.elements[1].Mass(),0.12)

if __name__ == '__main__':
    unittest.main()
