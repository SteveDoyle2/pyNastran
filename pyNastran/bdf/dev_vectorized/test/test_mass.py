import unittest

import os
from numpy import array
import pyNastran
from pyNastran.bdf.dev_vectorized.bdf import BDF
from pyNastran.utils import object_methods

rootpath = pyNastran.__path__[0]
testpath = os.path.join(rootpath, 'bdf', 'test', 'unit')


class TestMass(unittest.TestCase):

    def verify_pcomp_element(self, element, mass, area, centroid, normal):
        assert not isinstance(element, list), element
        #print object_methods(element, 'all')
        self.assertAlmostEqual(element.get_mass()[0], mass, msg='mass=%s expected=%s' % (element.get_mass()[0], mass))
        self.assertAlmostEqual(element.get_area()[0], area, msg='area=%s expected=%s' % (element.get_area()[0], area))
        self.assertTrue(all(element.get_centroid()[0] == centroid), msg='centroid=%s expected=%s' % (element.get_centroid()[0], centroid))
        self.assertTrue(all(element.get_normal()[0] == normal), msg='normal=%s expected=%s' % (element.get_normal()[0], normal))
        with self.assertRaises(AttributeError):
            element.Volume()
        with self.assertRaises(AttributeError):
            element.Length()

    def verify_pshell_element(self, element, mass, area, centroid, normal, nsm):
        assert not isinstance(element, list), element
        #print object_methods(element, 'all')
        self.assertAlmostEqual(element.get_mass()[0], mass, msg='mass=%s expected=%s' % (element.get_mass()[0], mass))
        self.assertAlmostEqual(element.get_area()[0], area, msg='area=%s expected=%s' % (element.get_area()[0], area))
        self.assertAlmostEqual(element.get_nsm()[0], nsm, msg='nsm=%s expected=%s' % (element.get_nsm()[0], nsm))
        self.assertTrue(all(element.get_centroid()[0] == centroid), msg='centroid=%s expected=%s' % (element.get_centroid()[0], centroid))
        self.assertTrue(all(element.get_normal()[0] == normal), msg='normal=%s expected=%s' % (element.get_normal()[0], normal))
        with self.assertRaises(AttributeError):
            element.Volume()
        with self.assertRaises(AttributeError):
            element.Length()

    def verify_psolid_element(self, element, mass, volume, centroid, rho, E=None, G=None, nu=None):
        assert not isinstance(element, list), element
        #print object_methods(element, 'all')
        self.assertAlmostEqual(element.get_volume(), volume, msg='volume=%s expected=%s' % (element.get_volume(), volume))
        self.assertTrue(all(element.get_centroid()[0] == centroid), msg='centroid=%s expected=%s' % (element.get_centroid(), centroid))
        if rho:
            #print "rho =", element.get_density(), rho
            self.assertAlmostEqual(element.get_density()[0], rho, msg='rho=%s expected=%s' % (element.get_density(), rho))
            #self.assertAlmostEqual(element.pid.Rho(), rho, msg='rho=%s expected=%s' % (element.pid.Rho(), rho))
            #self.assertEqual(element.pid.mid.type, 'MAT1', msg='mass=%s expected=%s' % (element.pid.mid.type, 'MAT1'))
            #self.assertAlmostEqual(element.pid.mid.Rho(), rho, msg='rho=%s expected=%s' % (element.pid.mid.Rho(), rho))

        self.assertAlmostEqual(element.get_mass(), mass, msg='mass=%s expected=%s' % (element.get_mass(), mass))

        #if E:
            #self.assertAlmostEqual(element.E(), E, msg='E=%s expected=%s' % (element.E(), E))
            #self.assertAlmostEqual(element.pid.E(), E, msg='E=%s expected=%s' % (element.pid.E(), E))
            #self.assertAlmostEqual(element.pid.mid.E(), E, msg='E=%s expected=%s' % (element.pid.mid.E(), E))
        #if G:
            #self.assertAlmostEqual(element.G(), G, msg='G=%s expected=%s' % (element.G(), G))
            #self.assertAlmostEqual(element.pid.G(), G, msg='G=%s expected=%s' % (element.pid.G(), G))
            #self.assertAlmostEqual(element.pid.mid.G(), G, msg='G=%s expected=%s' % (element.pid.mid.G(), G))
        #if nu:
            #self.assertAlmostEqual(element.Nu(), nu, msg='nu=%s expected=%s' % (element.Nu(), nu))
            #self.assertAlmostEqual(element.pid.Nu(), nu, msg='nu=%s expected=%s' % (element.pid.Nu(), nu))
            #self.assertAlmostEqual(element.pid.mid.Nu(), nu, msg='nu=%s expected=%s' % (element.pid.mid.Nu(), nu))

        with self.assertRaises(AttributeError):
            element.Area()
        with self.assertRaises(AttributeError):
            element.Length()

    def test_mass_shell_1(self):  # passes
        model = BDF(debug=False, log=None)
        bdfname = os.path.join(testpath, 'test_mass.dat')
        model.read_bdf(bdfname, include_dir=None, xref=True)

        # quad - pcomp
        quad = model.elements[1]
        mass = 0.12
        area = 1.0
        centroid = array([.5, .5, 0.])
        normal = array([.0, .0, 1.])
        self.verify_pcomp_element(quad, mass, area, centroid, normal)

        # quad - pshell, nsm=0
        quad = model.elements[3]
        mass = 0.0125
        nsm = 0.
        self.verify_pshell_element(quad, mass, area, centroid, normal, nsm)

        # quad - pshell, nsm=1
        quad = model.elements[5]
        mass = 1.0125 # mass w/o nsm + 1.0 b/c area=1
        nsm = 1.
        self.verify_pshell_element(quad, mass, area, centroid, normal, nsm)

        # tri - pcomp
        tri = model.elements[2]
        mass = 0.06
        area = 0.5
        centroid = array([2/3., 1/3., 0.])
        normal = array([.0,  .0,    1.])
        self.verify_pcomp_element(tri, mass, area, centroid, normal)

        # tri - pshell, nsm=0
        tri = model.elements[4]
        mass = 0.00625
        nsm = 0.
        self.verify_pshell_element(tri, mass, area, centroid, normal, nsm)

        # tri - pshell, nsm=1
        tri = model.elements[6]
        mass = 0.50625 # mass w/o nsm + 0.5 b/c area=0.5
        nsm = 1.
        self.verify_pshell_element(tri, mass, area, centroid, normal, nsm)

    def test_mass_solid_1(self):  # passes
        model = BDF(debug=False, log=None)
        bdfname = os.path.join(testpath, 'test_mass.dat')
        model.read_bdf(bdfname, include_dir=None, xref=True)

        # hexa - psolid - nsm = 0
        #print model.elements[7:8]
        #print model.elements[[7,8]]
        model.elements[7:9]
        model.elements[7:9:2]
        model.elements[1:100]
        print model.elements['cat']
        print model.get_elements('cat')
        #hexa = model.get_elements(7)
        #hexa = model.get_elements(7)
        #print hexa
        hexa = model.elements[7]
        mass = 0.2
        volume = 2. # l * w * h = 1 * 1 * 2
        rho = 0.1
        E = 1.0
        G = 2.0
        nu = 3.0
        centroid = array([0.5, 0.5, 1.0])
        self.verify_psolid_element(hexa, mass, volume, centroid, rho, E, G, nu)

        # tetra - psolid
        tetra = model.get_elements(8)
        mass = 1/30.
        volume = 1/3. # 1/3 * b * h = 1/3 * 0.5 * 2.0
        rho = 0.1
        E = 1.0
        G = 2.0
        nu = 3.0
        centroid = array([0.5, 0.25, 0.5])
        self.verify_psolid_element(tetra[0], mass, volume, centroid, rho, E, G, nu)

        # penta - psolid
        penta = model.get_elements(9)
        mass = 0.1
        volume = 1.0 # b * h = 0.5 * 2
        rho = 0.1
        E = 1.0
        G = 2.0
        nu = 3.0
        centroid = array([2/3., 1/3., 1.])
        self.verify_psolid_element(penta[0], mass, volume, centroid, rho, E, G, nu)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()