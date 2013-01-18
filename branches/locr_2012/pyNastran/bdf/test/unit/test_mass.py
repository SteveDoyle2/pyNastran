import unittest

import os
from numpy import array
import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.utils import object_methods

rootpath = pyNastran.__path__[0]
testpath = os.path.join(rootpath, 'bdf', 'test', 'unit')


class TestMass(unittest.TestCase):
        
    def verify_pcomp_element(self, element, mass, area, centroid, normal):
        #print object_methods(quad,'all')
        self.assertAlmostEqual(element.Mass(),mass)
        self.assertAlmostEqual(element.Area(),area)
        self.assertTrue(all(element.Centroid() == centroid))
        self.assertTrue(all(element.Normal()   == normal))
        with self.assertRaises(NotImplementedError): #  TODO: remove method (quad/tri)
            element.Volume()
        with self.assertRaises(NotImplementedError): #  TODO: remove method (quad/tri)
            element.Length()

    def verify_pshell_element(self, element, mass, area, centroid, normal, nsm):
        #print object_methods(quad,'all')
        self.assertAlmostEqual(element.Mass(),mass)
        self.assertAlmostEqual(element.Area(),area)
        self.assertAlmostEqual(element.Nsm(),nsm)
        self.assertTrue(all(element.Centroid() == centroid))
        self.assertTrue(all(element.Normal()   == normal))
        with self.assertRaises(NotImplementedError): #  TODO: remove method (quad/tri)
            element.Volume()
        with self.assertRaises(NotImplementedError): #  TODO: remove method (quad/tri)
            element.Length()

    def test_mass1(self):  # passes
        mesh = BDF(debug=True, log=None)
        bdfname = os.path.join(testpath, 'test_mass.dat')
        mesh.readBDF(bdfname, includeDir=None, xref=True)
        
        # quad - pcomp
        quad = mesh.elements[1]
        mass = 0.12
        area = 1.0
        centroid = array([.5, .5, 0.])
        normal = array([.0, .0, 1.])
        self.verify_pcomp_element(quad, mass, area, centroid, normal)
        
        # quad - pshell, nsm=0
        quad = mesh.elements[3]
        mass = 0.0125
        nsm = 0.
        self.verify_pshell_element(quad, mass, area, centroid, normal, nsm)

        # quad - pshell, nsm=1
        quad = mesh.elements[5]
        mass = 1.0125 # mass w/o nsm + 1.0 b/c area=1
        nsm = 1.
        self.verify_pshell_element(quad, mass, area, centroid, normal, nsm)

        # tri - pcomp
        tri = mesh.elements[2]
        mass = 0.06
        area = 0.5
        centroid = array([2/3., 1/3., 0.])
        normal = array([.0,  .0,    1.])
        self.verify_pcomp_element(tri, mass, area, centroid, normal)

        # tri - pshell, nsm=0
        tri = mesh.elements[4]
        mass = 0.00625
        nsm = 0.
        self.verify_pshell_element(tri, mass, area, centroid, normal, nsm)

        # tri - pshell, nsm=1
        tri = mesh.elements[6]
        mass = 0.50625 # mass w/o nsm + 0.5 b/c area=0.5
        nsm = 1.
        self.verify_pshell_element(tri, mass, area, centroid, normal, nsm)


if __name__ == '__main__':
    unittest.main()
