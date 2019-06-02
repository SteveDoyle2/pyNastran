import unittest

import os
from numpy import array, allclose
import pyNastran
from pyNastran.dev.bdf_vectorized.bdf import BDF
#from pyNastran.utils import object_methods

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'bdf', 'mesh_utils')


class TestMass(unittest.TestCase):

    def verify_pcomp_element(self, element, prop, nsm, thickness, mass, area, centroid, normal):
        assert not isinstance(element, list), element
        #print(object_methods(element, 'all'))
        self.assertAlmostEqual(element.get_area_by_element_id()[0], area, msg=
                               'PCOMP: area=%s expected=%s\n%s%s' % (element.get_area_by_element_id(), area, element, prop))
        self.assertAlmostEqual(element.get_thickness_by_element_id()[0], thickness,
                               msg='PCOMP: thickness=%s expected=%s\n%s%s' % (element.get_thickness_by_element_id(), thickness, element, prop))
        self.assertAlmostEqual(element.get_nonstructural_mass_by_element_id()[0], nsm,
                               msg='PCOMP: nsm=%s expected=%s\n%s%s' % (element.get_nonstructural_mass_by_element_id(), nsm, element, prop))
        self.assertTrue(all(element.get_centroid_by_element_id()[0] == centroid),
                        msg='PCOMP: centroid=%s expected=%s\n%s%s' % (element.get_centroid_by_element_id(), centroid, element, prop))
        self.assertTrue(all(element.get_normal_by_element_id()[0] == normal),
                        msg='PCOMP: normal=%s expected=%s\n%s%s' % (element.get_normal_by_element_id(), normal, element, prop))
        self.assertAlmostEqual(element.get_mass_by_element_id()[0], mass,
                               msg='PCOMP: mass=%s expected=%s\n%s%s' % (element.get_mass_by_element_id(), mass, element, prop))

    def verify_pshell_element(self, element, prop, mat, rho, mass, area, centroid, normal, nsm):
        assert not isinstance(element, list), element
        #print(object_methods(element, 'all'))
        #print("nsm = %s" % element.get_nonstructural_mass_by_element_id()[0])
        #print("density = %s" % element.get_density_by_element_id())
        if rho is not None:
            #print("%s.rho = %s %s" % (element.type, element.get_density_by_element_id(), rho))
            self.assertAlmostEqual(element.get_density_by_element_id()[0], rho, msg='PSHELL: rho=%s expected=%s\n%s%s%s' % (element.get_density_by_element_id(), rho, element, prop, mat))
        self.assertAlmostEqual(element.get_area_by_element_id()[0], area, msg='PSHELL: area=%s expected=%s\n%s%s%s' % (element.get_area_by_element_id(), area, element, prop, mat))
        self.assertAlmostEqual(element.get_nonstructural_mass_by_element_id()[0], nsm, msg='PSHELL: nsm=%s expected=%s\n%s%s%s' % (element.get_nonstructural_mass_by_element_id(), nsm, element, prop, mat))
        self.assertAlmostEqual(element.get_mass_by_element_id()[0], mass, msg='PSHELL: mass=%s expected=%s\n%s%s%s' % (element.get_mass_by_element_id()[0], mass, element, prop, mat))
        self.assertTrue(all(element.get_centroid_by_element_id()[0] == centroid), msg='PSHELL: centroid=%s expected=%s\n%s%s%s' % (element.get_centroid_by_element_id(), centroid, element, prop, mat))
        self.assertTrue(all(element.get_normal_by_element_id()[0] == normal), msg='PSHELL: normal=%s expected=%s\n%s%s%s' % (element.get_normal_by_element_id(), normal, element, prop, mat))

        assert mat is None

    def verify_psolid_element(self, element, mass, volume, centroid, rho, E=None, G=None, nu=None):
        #assert not isinstance(element, list), element
        #print(object_methods(element, 'all'))
        self.assertAlmostEqual(element.get_volume_by_element_id(), volume, msg='volume=%s expected=%s' % (element.get_volume_by_element_id(), volume))
        self.assertTrue(all(element.get_centroid_by_element_id()[0] == centroid), msg='centroid=%s expected=%s' % (element.get_centroid_by_element_id(), centroid))
        if rho:
            #print("%s.rho = %s %s" % (element.type, element.get_density_by_element_id(), rho))
            self.assertAlmostEqual(element.get_density_by_element_id()[0], rho, msg='rho=%s expected=%s' % (element.get_density_by_element_id(), rho))
            #self.assertAlmostEqual(element.pid.Rho(), rho, msg='rho=%s expected=%s' % (element.pid.Rho(), rho))
            #self.assertEqual(element.pid.mid.type, 'MAT1', msg='mass=%s expected=%s' % (element.pid.mid.type, 'MAT1'))
            #self.assertAlmostEqual(element.pid.mid.Rho(), rho, msg='rho=%s expected=%s' % (element.pid.mid.Rho(), rho))

        self.assertAlmostEqual(element.get_mass_by_element_id(), mass,
                               msg='mass=%s expected=%s' % (element.get_mass_by_element_id(), mass))

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
        bdfname = os.path.join(TEST_PATH, 'test_mass.dat')
        model.read_bdf(bdfname, xref=True)

        ###########
        # QUADS
        centroid = array([.5, .5, 0.])
        normal = array([.0, .0, 1.])
        ###########
        # quad - pcomp
        quad = model.elements[1]
        i = model.properties_shell.pcomp.get_property_index_by_property_id(quad.property_id)
        prop = model.properties_shell.pcomp[i]
        #mat = model.properties_shell.pshell[prop.material_id]
        mass = 0.12
        area = 1.0
        nsm = 0.
        thickness = 0.7
        rho1 = 0.1
        rho10 = 0.2
        t1 = 0.1
        t10 = 0.5
        # there are two layers of t1
        mpa = (2. * rho1 * t1 + rho10 * t10) + nsm
        mass2 = mpa * area
        assert allclose(mass, mass2), 'mass=%s mass2=%s diff=%s' % (mass, mass2, abs(mass - mass2))
        self.verify_pcomp_element(quad, prop, nsm, thickness, mass, area, centroid, normal)

        rho = None

        # quad - pshell, nsm=0
        eid = 3
        pid = 2
        quad = model.elements[eid]
        prop = model.properties_shell.pshell.slice_by_property_id(pid)
        #mat = model.materials[prop.material_id]  # good?
        mat = None
        rho = 0.1
        mass = 0.0125
        t = 0.125
        nsm = 0.
        area = 1.
        mass2 = area * (rho * t + nsm)
        assert allclose(mass, mass2), 'eid=%s pid=%s mass=%s mass2=%s diff=%s\n%s%s%s\nrho=%s A=%s t=%s nsm=%s' % (
            eid, pid, mass, mass2, abs(mass - mass2), quad, prop, mat, rho, area, t, nsm)
        centroid = array([.5, .5, 0.])
        normal = array([.0, .0, 1.])
        self.verify_pshell_element(quad, prop, mat, rho, mass, area, centroid, normal, nsm)

        # quad - pshell, nsm=1
        quad = model.elements[5]
        prop = model.properties_shell.pshell.slice_by_property_id(quad.property_id)
        #mat = model.properties_shell.pshell[prop.material_id]
        mat = None
        mass = 1.0125 # mass w/o nsm + 1.0 b/c area=1
        nsm = 1.
        self.verify_pshell_element(quad, prop, mat, rho, mass, area, centroid, normal, nsm)

        ###########
        # TRIS
        centroid = array([2., 1., 0.]) / 3.
        normal = array([.0, .0, 1.])
        ###########
        # tri - pcomp
        tri = model.elements[2]
        prop = model.properties_shell.pcomp.slice_by_property_id(tri.property_id)
        #mat = model.properties_shell.pshell[prop.material_id]
        mass = 0.06
        area = 0.5
        nsm = 0.
        thickness = 0.7
        self.verify_pcomp_element(tri, prop, nsm, thickness, mass, area, centroid, normal)

        # tri - pshell, nsm=0
        tri = model.elements[4]
        i = model.properties_shell.pshell.get_property_index_by_property_id(tri.property_id)
        prop = model.properties_shell.pshell[i]
        #mat = model.properties_shell.pshell[prop.material_id]
        mat = None
        mass = 0.00625
        nsm = 0.
        self.verify_pshell_element(tri, prop, mat, rho, mass, area, centroid, normal, nsm)

        # tri - pshell, nsm=1
        tri = model.elements[6]
        prop = model.properties_shell.pshell.slice_by_property_id(tri.property_id)
        #mat = model.properties_shell.pshell[prop.material_id]
        mat = None
        mass = 0.50625 # mass w/o nsm + 0.5 b/c area=0.5
        nsm = 1.
        self.verify_pshell_element(tri, prop, mat, rho, mass, area, centroid, normal, nsm)

    def test_bad_01(self):
        model = BDF(debug=False, log=None)
        bdfname = os.path.join(TEST_PATH, 'test_mass.dat')
        model.read_bdf(bdfname, xref=True)

        # this passes silently
        print(model.elements[['cat']])

        # this does not
        with self.assertRaises(TypeError):
            print(model.elements[None])

        #print('elements = %r' % model.elements['cat'])
        with self.assertRaises(KeyError):
            model.elements['cat']

    def test_combo_1(self):
        model = BDF(debug=False, log=None)
        bdfname = os.path.join(TEST_PATH, 'test_mass.dat')
        model.read_bdf(bdfname, xref=True)

        # these are valid
        eids, mass = model.elements.get_mass_by_element_id([8, 9])
        print('massA = %s' % mass)
        eids, mass = model.elements.get_mass_by_element_id(range(1, 10))
        print('massB = %s' % mass)

        # no analysis - out of range
        elements = model.elements[[100000, 100001]]
        eids, mass = model.elements.get_mass_by_element_id(range(100000, 100005))
        print('elementsC = %s' % eids)
        print('massC = %s' % mass)

        eids, mass = model.elements.get_mass_by_element_id(range(-10, -5))
        print('elementsD = %s' % eids)
        print('massD = %s' % mass)

        print('-------------------------')
        eids, mass = model.elements.get_mass_by_element_id(range(-3, 20))
        print('massE = %s' % mass)
        print('eidsE = %s' % eids)
        print('\neid   mass')
        print('----------')
        for eidi, massi in zip(eids, mass):
            print('%-5s %-5s' % (eidi, massi))
        #print(vstack([mass, eids]))

    def test_mass_solid_1(self):  # passes
        model = BDF(debug=False, log=None)
        bdfname = os.path.join(TEST_PATH, 'test_mass.dat')
        model.read_bdf(bdfname, xref=True)

        # hexa - psolid - nsm = 0
        #print(model.elements[7:8])
        #print(model.elements[[7,8]])
        model.elements[7:9]
        model.elements[7:9:2]
        model.elements[1:100]
        #hexa = model.get_elements(7)
        #hexa = model.get_elements(7)
        #print(hexa)
        hexa = model.elements[7]
        print('hexa =', hexa)
        mass = 0.2
        volume = 2. # l * w * h = 1 * 1 * 2
        rho = 0.1
        E = 1.0
        G = 2.0
        nu = 3.0
        centroid = array([0.5, 0.5, 1.0])
        self.verify_psolid_element(hexa, mass, volume, centroid, rho, E, G, nu)

        # tetra - psolid
        tetra = model.elements[8]
        mass = 1/30.
        volume = 1/3. # 1/3 * b * h = 1/3 * 0.5 * 2.0
        rho = 0.1
        E = 1.0
        G = 2.0
        nu = 3.0
        centroid = array([0.5, 0.25, 0.5])
        self.verify_psolid_element(tetra, mass, volume, centroid, rho, E, G, nu)

        # penta - psolid
        penta = model.elements[9]
        mass = 0.1
        volume = 1.0 # b * h = 0.5 * 2
        rho = 0.1
        E = 1.0
        G = 2.0
        nu = 3.0
        centroid = array([2/3., 1/3., 1.])
        self.verify_psolid_element(penta, mass, volume, centroid, rho, E, G, nu)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
