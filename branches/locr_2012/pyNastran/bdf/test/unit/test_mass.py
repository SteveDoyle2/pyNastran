from numpy import array, allclose
import unittest

from pyNastran.bdf.bdf import BDF


class TestMass(unittest.TestCase):
    def test_mass1(self):  # passes
        mesh = BDF(debug=True, log=None) 
        mesh.readBDF('test_mass.dat', includeDir=None, xref=True) 
        self.assertAlmostEqual(mesh.elements[1].Mass(),0.12)

if __name__ == '__main__':
    unittest.main()
