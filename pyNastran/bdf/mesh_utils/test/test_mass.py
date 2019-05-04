import unittest

import os
import numpy as np
import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.utils import object_methods

PKG_PATH = pyNastran.__path__[0]
mesh_utils_path = os.path.join(PKG_PATH, 'bdf', 'mesh_utils', 'test')


class TestMass(unittest.TestCase):

    def verify_pcomp_element(self, element, mass, area, centroid, normal):
        #print(object_methods(element, 'all'))
        self.assertAlmostEqual(element.Mass(), mass, msg='mass=%s expected=%s' % (element.Mass(), mass))
        self.assertAlmostEqual(element.Area(), area, msg='area=%s expected=%s' % (element.Area(), area))
        self.assertTrue(all(element.Centroid() == centroid), msg='centroid=%s expected=%s' % (element.Centroid(), centroid))
        self.assertTrue(all(element.Normal() == normal), msg='normal=%s expected=%s' % (element.Normal(), normal))
        with self.assertRaises(AttributeError):
            element.Volume()
        with self.assertRaises(AttributeError):
            element.Length()

    def verify_pshell_element(self, element, mass, area, centroid, normal, nsm):
        #print(object_methods(element, 'all'))
        self.assertAlmostEqual(element.Mass(), mass, msg='mass=%s expected=%s' % (element.Mass(), mass))
        self.assertAlmostEqual(element.Area(), area, msg='area=%s expected=%s' % (element.Area(), area))
        self.assertAlmostEqual(element.Nsm(), nsm, msg='nsm=%s expected=%s' % (element.Nsm(), nsm))
        self.assertTrue(all(element.Centroid() == centroid), msg='centroid=%s expected=%s' % (element.Centroid(), centroid))
        self.assertTrue(all(element.Normal() == normal), msg='normal=%s expected=%s' % (element.Normal(), normal))
        with self.assertRaises(AttributeError):
            element.Volume()
        with self.assertRaises(AttributeError):
            element.Length()

    def verify_psolid_element(self, element, mass, volume, centroid, rho, E=None, G=None, nu=None):
        #print(object_methods(element, 'all'))
        self.assertAlmostEqual(element.Volume(), volume, msg='volume=%s expected=%s' % (element.Volume(), volume))
        self.assertTrue(all(element.Centroid() == centroid), msg='centroid=%s expected=%s' % (element.Centroid(), centroid))
        if rho:
            self.assertAlmostEqual(element.Rho(), rho, msg='rho=%s expected=%s' % (element.Rho(), rho))
            self.assertAlmostEqual(element.pid_ref.Rho(), rho, msg='rho=%s expected=%s' % (element.pid_ref.Rho(), rho))
            self.assertEqual(element.pid_ref.mid_ref.type, 'MAT1', msg='mass=%s expected=%s' % (element.pid_ref.mid_ref.type, 'MAT1'))
            self.assertAlmostEqual(element.pid_ref.mid_ref.Rho(), rho, msg='rho=%s expected=%s' % (element.pid_ref.mid_ref.Rho(), rho))

        self.assertAlmostEqual(element.Mass(), mass, msg='mass=%s expected=%s' % (element.Mass(), mass))

        if E:
            self.assertAlmostEqual(element.E(), E, msg='E=%s expected=%s' % (element.E(), E))
            self.assertAlmostEqual(element.pid_ref.E(), E, msg='E=%s expected=%s' % (element.pid_ref.E(), E))
            self.assertAlmostEqual(element.pid_ref.mid_ref.E(), E, msg='E=%s expected=%s' % (element.pid_ref.mid_ref.E(), E))
        if G:
            self.assertAlmostEqual(element.G(), G, msg='G=%s expected=%s' % (element.G(), G))
            self.assertAlmostEqual(element.pid_ref.G(), G, msg='G=%s expected=%s' % (element.pid_ref.G(), G))
            self.assertAlmostEqual(element.pid_ref.mid_ref.G(), G, msg='G=%s expected=%s' % (element.pid_ref.mid_ref.G(), G))
        if nu:
            self.assertAlmostEqual(element.Nu(), nu, msg='nu=%s expected=%s' % (element.Nu(), nu))
            self.assertAlmostEqual(element.pid_ref.Nu(), nu, msg='nu=%s expected=%s' % (element.pid_ref.Nu(), nu))
            self.assertAlmostEqual(element.pid_ref.mid_ref.Nu(), nu, msg='nu=%s expected=%s' % (element.pid_ref.mid_ref.Nu(), nu))

        with self.assertRaises(AttributeError):
            element.Area()
        with self.assertRaises(AttributeError):
            element.Length()

    def test_mass_shell_1(self):  # passes
        model = BDF(debug=False, log=None)
        bdfname = os.path.join(mesh_utils_path, 'test_mass.dat')
        model.read_bdf(bdfname, xref=True)

        # quad - pcomp
        quad = model.elements[1]
        mass = 0.12
        area = 1.0
        centroid = np.array([.5, .5, 0.])
        normal = np.array([.0, .0, 1.])
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
        centroid = np.array([2/3., 1/3., 0.])
        normal = np.array([.0, .0, 1.])
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
        bdfname = os.path.join(mesh_utils_path, 'test_mass.dat')
        model.read_bdf(bdfname, xref=True)

        # hexa - psolid - nsm = 0
        hexa = model.elements[7]
        mass = 0.2
        volume = 2. # l * w * h = 1 * 1 * 2
        rho = 0.1
        E = 1.0
        G = 2.0
        nu = 3.0
        centroid = np.array([0.5, 0.5, 1.0])
        self.verify_psolid_element(hexa, mass, volume, centroid, rho, E, G, nu)

        # tetra - psolid
        tetra = model.elements[8]
        mass = 1/30.
        volume = 1/3. # 1/3 * b * h = 1/3 * 0.5 * 2.0
        rho = 0.1
        E = 1.0
        G = 2.0
        nu = 3.0
        centroid = np.array([0.5, 0.25, 0.5])
        self.verify_psolid_element(tetra, mass, volume, centroid, rho, E, G, nu)

        # penta - psolid
        penta = model.elements[9]
        mass = 0.1
        volume = 1.0 # b * h = 0.5 * 2
        rho = 0.1
        E = 1.0
        G = 2.0
        nu = 3.0
        centroid = np.array([2/3., 1/3., 1.])
        self.verify_psolid_element(penta, mass, volume, centroid, rho, E, G, nu)

        #wtmass = 0.00259
        mass, cg, I = model.mass_properties(element_ids=None, mass_ids=None,
                                           reference_point=None,
                                           sym_axis=None,
                                           scale=None)

        model.params['WTMASS'].values[0] = 1.0
        mass2, cg, I = model.mass_properties(element_ids=None, mass_ids=None,
                                             reference_point=None,
                                             sym_axis=None,
                                             scale=None)
        assert np.allclose(mass, 0.005311658333), 'mass=%s' % mass
        assert np.allclose(mass2, 2.050833333), 'mass2=%s' % mass2

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
