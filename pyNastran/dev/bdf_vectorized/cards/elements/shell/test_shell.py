from collections import defaultdict
from io import StringIO
import unittest

import numpy as np
from numpy import array, allclose

from pyNastran.dev.bdf_vectorized.bdf import BDF
#from pyNastran.dev.bdf_vectorized.cards.materials.mat1 import MAT1
#from pyNastran.dev.bdf_vectorized.cards.elements.shell.pcomp import PCOMP
from pyNastran.dev.bdf_vectorized.cards.elements.shell.pshell import PSHELL

class TestShellsV1(unittest.TestCase):
    def _make_cquad4(self, model, rho, nu, G, E, t, nsm):
        eid = 10
        pid = 20
        mid = 30
        n1 = 1
        n2 = 2
        n3 = 3
        n4 = 4
        A = 2.
        mid2 = mid3 = mid4 = twelveIt3 = tst = z1 = z2 = None

        mass = A * (t * rho + nsm)
        cards = [
            ['grid', n1, 0, 0., 0., 0.],
            ['grid', n2, 0, 2., 0., 0.],
            ['grid', n3, 0, 2., 1., 0.],
            ['grid', n4, 0, 0., 1., 0.],
            ['cquad4', eid, pid, n1, n2, n3, n4],
            ['pshell', pid, mid, t, mid2, twelveIt3, mid3, tst, nsm, z1, z2],
            ['mat1', mid, E, G, nu, rho],
        ]
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.build()
        cquad4 = model.elements[eid]

        # cquad4 / pshell
        self.assertEqual(cquad4.get_element_id_by_element_index(), eid)
        self.assertEqual(cquad4.get_property_id_by_element_index(), pid)
        #self.assertEqual(cquad4.Mid(), mid)
        #self.assertEqual(cquad4.Nsm(), nsm)
        self.assertEqual(cquad4.get_mass_by_element_id(), mass)
        self.assertAlmostEquals(cquad4.get_mass_per_area_by_element_id(), mass / A)
        self.assertEqual(cquad4.get_area_by_element_id(), A)
        self.assertEqual(cquad4.get_thickness_by_element_id(), t)
        #self.assertEqual(cquad4.Rho(), rho)  # removed because of PCOMP

    def _make_ctria3(self, model, rho, nu, G, E, t, nsm):
        eid = 10
        pid = 20
        mid = 30
        n1 = 1
        n2 = 2
        n3 = 3
        mid2 = mid3 = mid4 = twelveIt3 = tst = z1 = z2 = None
        z0 = sb = ft = Tref = ge = lam = None
        sout = None
        theta0 = 0.
        theta1 = 30.
        theta2 = 60.
        theta3 = 90.
        A = 2.
        cards = [
            ['grid', n1, 0, 0., 0., 0.],
            ['grid', n2, 0, 4., 0., 0.],
            ['grid', n3, 0, 4., 1., 0.],
            ['ctria3', eid, pid, n1, n2, n3],  # A = 1/2 * 4 * 1 = 2.
            ['pshell', pid, mid, t, mid2, twelveIt3, mid3, tst, nsm, z1, z2, mid4],

            ['ctria3', eid + 1, pid + 1, n1, n2, n3],  # A = 1/2 * 4 * 1 = 2.
            [
                'pcomp', pid + 1, z0, nsm, sb, ft, Tref, ge, lam,
                mid, t, theta0, sout,
                mid, 2 * t, theta1, sout,
                mid, 3 * t, theta2, sout,
                mid, 4 * t, theta3, sout,
                ],
            ['mat1', mid, E, G, nu, rho],
        ]
        card_count = defaultdict(int)
        for card in cards:
            card_count[card[0].upper()] += 1
        model.allocate(card_count)

        for card in cards:
            model.add_card(card, card[0], is_list=True)
        model.build()

        # ctria3 / pshell
        ctria3 = model.elements[eid]
        mass = A * (t * rho + nsm)
        self.assertEqual(ctria3.get_element_id_by_element_index(), eid)
        self.assertEqual(ctria3.get_property_id_by_element_index(), pid)
        #self.assertEqual(ctria3.Mid(), mid)
        #self.assertEqual(ctria3.Nsm(), nsm)
        self.assertEqual(ctria3.get_mass_by_element_id(), mass)
        self.assertAlmostEquals(ctria3.get_mass_per_area_by_element_id(), mass / A)
        self.assertEqual(ctria3.get_area_by_element_id(), A)
        self.assertEqual(ctria3.get_thickness_by_element_id(), t)
        #self.assertEqual(ctria3.MassPerArea(), mass / A)

        # removed because of PCOMP
        # also no E, G, J, Nu, for the same reason
        # what about Mid
        #self.assertEqual(ctria3.Rho(), rho)


        # pshell
        pshell = model.properties[pid]
        assert isinstance(pshell, PSHELL), type(pshell)
        self.assertEqual(pshell.get_property_id_by_property_index(), pid)
        #self.assertEqual(pshell.Mid(), mid)
        #self.assertEqual(pshell.Nsm(), nsm)
        #self.assertEqual(pshell.Thickness(), t)
        #self.assertEqual(pshell.Rho(), rho)
        self.assertEqual(pshell.z1[0], -t / 2.)
        self.assertEqual(pshell.z2[0], t / 2.)

        # ctria3 / pcomp
        ctria3 = model.elements[eid + 1]
        mass = A * (10 * t * rho + nsm)
        self.assertEqual(ctria3.get_element_id_by_element_index(), eid + 1)
        self.assertEqual(ctria3.get_property_id_by_element_id(), pid + 1)
        #self.assertEqual(ctria3.Mid(), mid)
        #self.assertEqual(ctria3.Nsm(), nsm)
        self.assertAlmostEquals(ctria3.get_mass_by_element_id(), mass)
        self.assertAlmostEquals(ctria3.get_mass_per_area_by_element_id(), mass / A)
        self.assertEqual(ctria3.get_area_by_element_id(), A)
        self.assertEqual(ctria3.get_thickness_by_element_id(), 10 * t)
        #self.assertEqual(ctria3.Rho(), rho)

        # pcomp
        pcomp_pid = pid + 1
        pcomp = model.properties.properties_shell.pcomp.slice_by_property_id(pcomp_pid)
        #print('pcomp =', type(pcomp))
        #self.assertEqual(pcomp.get_property_id()[0], pcomp_pid)
        self.assertEqual(pcomp.get_property_id_by_property_index(), pcomp_pid)
        self.assertEqual(pcomp.get_nplies_by_property_id(), 4)
        self.assertEqual(pcomp.get_nplies_by_property_index(), 4)

        self.assertEqual(pcomp.get_material_id_by_property_id_ply(pcomp_pid, 0), mid)
        self.assertEqual(pcomp.get_nonstructural_mass_by_property_id(), nsm)
        self.assertEqual(pcomp.get_nonstructural_mass_by_property_index(), nsm)

        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.get_material_id_by_property_id_ply(pcomp_pid, -1), mid)
        self.assertTrue(all(pcomp.get_material_ids_by_property_id(pcomp_pid)[0] == [mid] * 4))
        self.assertEqual(pcomp.get_material_id_by_property_id_ply(pcomp_pid, 0), mid)
        self.assertEqual(pcomp.get_material_id_by_property_id_ply(pcomp_pid, 1), mid)
        self.assertEqual(pcomp.get_material_id_by_property_id_ply(pcomp_pid, 2), mid)
        self.assertEqual(pcomp.get_material_id_by_property_id_ply(pcomp_pid, 3), mid)
        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.get_material_id_by_property_id_ply(pcomp_pid, 4), mid)

        #with self.assertRaises(IndexError):
            #self.assertEqual(pcomp.get_thickness_by_property_id_ply(pcomp_pid, -1), t)
        self.assertEqual(pcomp.get_thickness_by_property_id(), 10 * t)
        self.assertEqual(pcomp.get_thickness_by_property_id_ply(pcomp_pid, 0), t)
        self.assertEqual(pcomp.get_thickness_by_property_id_ply(pcomp_pid, 1), 2 * t)
        self.assertTrue(np.allclose(pcomp.get_thickness_by_property_id_ply(pcomp_pid, 2), 3 * t))
        self.assertEqual(pcomp.get_thickness_by_property_id_ply(pcomp_pid, 3), 4 * t)
        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.get_thickness_by_property_id_ply(pcomp_pid, 4), 5*t)

        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.get_density_by_property_id_ply(pcomp_pid, -1), rho)
        self.assertEqual(pcomp.get_density_by_property_id_ply(pcomp_pid, 0), rho)
        self.assertEqual(pcomp.get_density_by_property_id_ply(pcomp_pid, 1), rho)
        self.assertEqual(pcomp.get_density_by_property_id_ply(pcomp_pid, 2), rho)
        self.assertEqual(pcomp.get_density_by_property_id_ply(pcomp_pid, 3), rho)
        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.get_density_by_property_id_ply(pcomp_pid, 4), rho)

        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.get_theta_by_property_id_ply(pcomp_pid, -1), 0.)
        self.assertEqual(pcomp.get_theta_by_property_id_ply(pcomp_pid, 0), 0.)
        self.assertEqual(pcomp.get_theta_by_property_id_ply(pcomp_pid, 1), 30.)
        self.assertEqual(pcomp.get_theta_by_property_id_ply(pcomp_pid, 2), 60.)
        self.assertEqual(pcomp.get_theta_by_property_id_ply(pcomp_pid, 3), 90.)
        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.get_theta_by_property_id_ply(pcomp_pid, 4), rho)
        self.assertEqual(pcomp.z0, -10*t/2.)

    def test_pshell_01(self):
        """tests a CQUAD4 and a PSHELL"""

        rho = 0.1
        nu = 0.3
        G = None
        E = 1e7
        t = 0.3
        nsm = 0.0

        card_count = {
            'GRID': 4,
            'CQUAD4': 1,
            #'CTRIA3': 1,
            'PSHELL': 1,
            'PCOMP': 1,
            'MAT1': 1,
            'MAT8': 1,
        }
        #print('starting BDF1')
        model = BDF(debug=True)
        model.allocate(card_count)
        self._make_cquad4(model, rho, nu, G, E, t, nsm)

        card_count1 = {
            'GRID': 3,
            #'CQUAD4': 1,
            'CTRIA3': 2,
            'PSHELL': 1,
            'PCOMP': 1,
            'MAT1': 1,
            'MAT8': 1,
        }
        #print('starting BDF2')
        model = BDF(debug=True)
        model.allocate(card_count1)
        self._make_ctria3(model, rho, nu, G, E, t, nsm)

        card_count = {
            'GRID': 4,
            'CQUAD4': 1,
            #'CTRIA3': 2,
            'PSHELL': 1,
            'PCOMP': 1,
            'MAT1': 1,
            'MAT8': 1,
        }
        #print('starting BDF3')
        nsm = 1.0
        model = BDF(debug=False)
        model.allocate(card_count)
        self._make_cquad4(model, rho, nu, G, E, t, nsm)

        card_count = {
            'GRID': 3,
            #'CQUAD4': 1,
            'CTRIA3': 2,
            'PSHELL': 1,
            'PCOMP': 1,
            'MAT1': 1,
            'MAT8': 1,
        }
        #print('starting BDF4')
        model = BDF(debug=False)
        model.allocate(card_count)
        self._make_ctria3(model, rho, nu, G, E, t, nsm)


    def test_pcomp_01(self):
        """
        asymmetrical, nsm=0.0 and nsm=1.0
        """
        #self.pid = data[0]
        #self.z0 = data[1]
        #self.nsm = data[2]
        #self.sb = data[3]
        #self.ft = data[4]
        #self.tref = data[5]
        #self.ge = data[6]
        #self.lam = data[7]
        #Mid = data[8]
        #T = data[9]
        #Theta = data[10]
        #Sout = data[11]

        pid = 1
        z0 = 0.
        nsm = 0.
        sb = 0.
        ft = 'HILL'
        tref = 0.
        ge = 0.
        #lam = 'NO'  # is_symmetrical YES/NO
        lam = 'BEND'  # is_symmetrical YES/NO
        Mid = [1, 2, 3]
        Theta = [0., 10., 20.]
        T = [.1, .2, .3]
        Sout = ['YES', 'YES', 'NO']  # 0-NO, 1-YES
        card_lines = [
            'PCOMP', pid, z0, nsm, sb, ft, tref, ge, lam,
            Mid[0], T[0], Theta[0], Sout[0],
            Mid[1], T[1], Theta[1], Sout[1],
            Mid[2], T[2], Theta[2], Sout[2],
        ]
        model = BDF(debug=False)
        card_count = {
            'PCOMP' : 1,
            'MAT1' : 3,
        }
        model.allocate(card_count)

        model.add_card(card_lines, 'PCOMP', comment='', is_list=True)


        # material...
        mid = 1
        E = 1e7
        G = None
        nu = None
        rho = 1.0
        a = None
        St = None
        Sc = None
        Ss = None
        mcsid = None
        mat1_a = ['MAT1', mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
        mat1_b = ['MAT1', mid + 1, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
        mat1_c = ['MAT1', mid + 2, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
        model.add_card(mat1_a, 'MAT1', comment='', is_list=True)
        model.add_card(mat1_b, 'MAT1', comment='', is_list=True)
        model.add_card(mat1_c, 'MAT1', comment='', is_list=True)
        #card = BDFCard(mat1)
        #m = MAT1(model)
        #m.allocate(1)
        #m.add(card)
        #m.build()

        model.build()
        #card = BDFCard(data)
        #p = PCOMP(model)
        #p = model.properties.pcomp
        #p.add(card)
        #p.build()
        p = model.properties_shell.pcomp
        m = model.materials.mat1
        #self.assertFalse(p.is_symmetrical)
        self.assertEqual(p.get_nplies_by_property_id(), 3)

        self.assertAlmostEqual(p.get_thickness_by_property_id(pid), 0.6)
        self.assertAlmostEqual(p.get_thickness_by_property_id_ply(pid, 0), 0.1)
        self.assertAlmostEqual(p.get_thickness_by_property_id_ply(pid, 1), 0.2)
        self.assertAlmostEqual(p.get_thickness_by_property_id_ply(pid, 2), 0.3)
        with self.assertRaises(IndexError):
            p.get_thickness_by_property_id_ply(pid, 3)

        self.assertAlmostEqual(p.get_theta_by_property_id_ply(pid, 0), 0.)
        self.assertAlmostEqual(p.get_theta_by_property_id_ply(pid, 1), 10.)
        self.assertAlmostEqual(p.get_theta_by_property_id_ply(pid, 2), 20.)
        with self.assertRaises(IndexError):
            p.get_theta_by_property_id_ply(pid, 3)

        self.assertEqual(p.get_material_id_by_property_id_ply(pid, 0), 1)
        self.assertEqual(p.get_material_id_by_property_id_ply(pid, 1), 2)
        self.assertEqual(p.get_material_id_by_property_id_ply(pid, 2), 3)
        with self.assertRaises(IndexError):
            p.get_material_id_by_property_id_ply(pid, 3)

        #print('get_material_id_by_property_id = ', p.get_material_id_by_property_id(pid))
        self.assertEqual(p.get_material_ids_by_property_id(pid)[0], 1)
        self.assertEqual(p.get_material_ids_by_property_id(pid)[0], 1)
        self.assertEqual(p.get_material_ids_by_property_id(pid)[1], 2)
        self.assertEqual(p.get_material_ids_by_property_id(pid)[2], 3)

        self.assertEqual(p.get_sout_by_property_id_ply(pid, 0), 'YES')
        self.assertEqual(p.get_sout_by_property_id_ply(pid, 1), 'YES')
        self.assertEqual(p.get_sout_by_property_id_ply(pid, 2), 'NO')
        with self.assertRaises(IndexError):
            p.get_sout_by_property_id_ply(pid, 3)

        #for iply in range(len(p.plies)):
            #mid = p.plies[iply][0]
            #p.plies[iply][0] = m # MAT1
            ##p.mids = [m, m, m]

        f = StringIO()

        #print(m.write_card(f, size=8, material_id=None))
        p.write_card(f)
        m.write_card(f)
        print(f.getvalue())

        #Mid
        self.assertAlmostEqual(p.get_material_id_by_property_id_ply(pid, 0), 1)
        self.assertAlmostEqual(p.get_material_id_by_property_id_ply(pid, 1), 2)
        self.assertAlmostEqual(p.get_material_id_by_property_id_ply(pid, 2), 3)
        with self.assertRaises(IndexError):
            p.get_density_by_property_id_ply(pid, 3)

        #Rho
        self.assertAlmostEqual(p.get_density_by_property_id_ply(pid, 0), 1.0)
        self.assertAlmostEqual(p.get_density_by_property_id_ply(pid, 1), 1.0)
        self.assertAlmostEqual(p.get_density_by_property_id_ply(pid, 2), 1.0)
        with self.assertRaises(IndexError):
            p.get_density_by_property_id_ply(pid, 3)

        # MassPerArea
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id(), 0.6)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id(pid), 0.6)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id([pid]), 0.6)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 0), 0.1)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 1), 0.2)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 2), 0.3)
        #with self.assertRaises(IndexError):
            #p.MassPerArea(3)

        #----------------------
        # change the nsm to 1.0
        p.set_nonstructural_mass_by_property_id(pid, 1.0)

        self.assertEqual(p.get_nonstructural_mass_by_property_id(), 1.0)
        # MassPerArea
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id(), 1.6)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 0, method='nplies'), 0.1+1/3.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 1, method='nplies'), 0.2+1/3.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 2, method='nplies'), 0.3+1/3.)

        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 0, method='t'), 0.1+1/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 1, method='t'), 0.2+2/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 2, method='t'), 0.3+3/6.)

        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 0, method='rho*t'), 0.1+1/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 1, method='rho*t'), 0.2+2/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 2, method='rho*t'), 0.3+3/6.)

        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 0, method='t'), 0.1+1/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 1, method='t'), 0.2+2/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 2, method='t'), 0.3+3/6.)
        with self.assertRaises(IndexError):
            p.get_mass_per_area_by_property_id_ply(pid, 3, method='nplies')

        z0 = p.get_z0_by_property_id(pid)

        z = p.get_z_locations_by_property_id(pid)
        z_expected = array([0., T[0], T[0]+T[1], T[0]+T[1]+T[2]])
        for za, ze in zip(z, z_expected):
            self.assertAlmostEqual(za, ze)

        #z0  =
        p.z0[0] = 1.0
        z_expected = 1.0 + z_expected
        z = p.get_z_locations_by_property_id(pid)
        for za, ze in zip(z, z_expected):
            self.assertAlmostEqual(za, ze)

    def test_pcomp_02(self):
        """
        symmetrical, nsm=0.0 and nsm=1.0
        """
        model = BDF(debug=False)
        pid = 1
        z0 = 0.
        nsm = 0.
        sb = 0.
        ft = 'HOFF'
        tref = 0.
        ge = 0.
        lam = 'SYM'  # is_symmetrical SYM
        Mid = [1, 2, 3]
        Theta = [0., 10., 20.]
        T = [.1, .2, .3]
        Sout = ['YES', 'YES', 'NO']  # 0-NO, 1-YES
        pcomp = ['PCOMP', pid, z0, nsm, sb, ft, tref, ge, lam,
                Mid[0], T[0], Theta[0], Sout[0],
                Mid[1], T[1], Theta[1], Sout[1],
                Mid[2], T[2], Theta[2], Sout[2]]

        #----

        mid = 1
        E = 3.0e7
        G = None
        nu = None
        rho = 1.0
        a = None
        St = None
        Sc = None
        Ss = None
        mcsid = None
        mat1_a = ['MAT1', mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
        mat1_b = ['MAT1', mid + 1, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
        mat1_c = ['MAT1', mid + 2, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]

        card_count = {
            'PCOMP' : 1,
            'MAT1' : 3,
        }
        model.allocate(card_count)
        model.add_card(pcomp, pcomp[0], is_list=True)
        model.add_card(mat1_a, 'MAT1', is_list=True)
        model.add_card(mat1_b, 'MAT1', is_list=True)
        model.add_card(mat1_c, 'MAT1', is_list=True)
        model.build()

        p = model.properties.properties_shell.pcomp[0]
        self.assertTrue(p.is_symmetrical_by_property_id())
        self.assertTrue(p.is_symmetrical_by_property_index())
        self.assertEqual(p.get_nplies_by_property_id(), 6)

        self.assertAlmostEqual(p.get_thickness_by_property_id()[0], 1.2)
        self.assertAlmostEqual(p.get_thickness_by_property_id_ply(pid, 0), 0.1)
        self.assertAlmostEqual(p.get_thickness_by_property_id_ply(pid, 1), 0.2)
        self.assertAlmostEqual(p.get_thickness_by_property_id_ply(pid, 2), 0.3)
        self.assertAlmostEqual(p.get_thickness_by_property_id_ply(pid, 3), 0.1)
        self.assertAlmostEqual(p.get_thickness_by_property_id_ply(pid, 4), 0.2)
        self.assertAlmostEqual(p.get_thickness_by_property_id_ply(pid, 5), 0.3)
        with self.assertRaises(IndexError):
            p.get_thickness_by_property_id_ply(pid, 6)

        self.assertAlmostEqual(p.get_theta_by_property_id_ply(pid, 0), 0.)
        self.assertAlmostEqual(p.get_theta_by_property_id_ply(pid, 1), 10.)
        self.assertAlmostEqual(p.get_theta_by_property_id_ply(pid, 2), 20.)
        self.assertAlmostEqual(p.get_theta_by_property_id_ply(pid, 3), 0.)
        self.assertAlmostEqual(p.get_theta_by_property_id_ply(pid, 4), 10.)
        self.assertAlmostEqual(p.get_theta_by_property_id_ply(pid, 5), 20.)
        with self.assertRaises(IndexError):
            p.get_theta_by_property_id_ply(pid, 6)

        self.assertEqual(p.get_material_id_by_property_id_ply(pid, 0), 1)
        self.assertEqual(p.get_material_id_by_property_id_ply(pid, 1), 2)
        self.assertEqual(p.get_material_id_by_property_id_ply(pid, 2), 3)
        self.assertEqual(p.get_material_id_by_property_id_ply(pid, 3), 1)
        self.assertEqual(p.get_material_id_by_property_id_ply(pid, 4), 2)
        self.assertEqual(p.get_material_id_by_property_id_ply(pid, 5), 3)
        with self.assertRaises(IndexError):
            p.get_material_id_by_property_id_ply(pid, 6)

        self.assertTrue(allclose(
            p.get_material_ids_by_property_id(pid),
            [1, 2, 3, 1, 2, 3])
        )

        self.assertEqual(p.get_sout_by_property_id_ply(pid, 0), 'YES')
        self.assertEqual(p.get_sout_by_property_id_ply(pid, 1), 'YES')
        self.assertEqual(p.get_sout_by_property_id_ply(pid, 2), 'NO')
        self.assertEqual(p.get_sout_by_property_id_ply(pid, 3), 'YES')
        self.assertEqual(p.get_sout_by_property_id_ply(pid, 4), 'YES')
        self.assertEqual(p.get_sout_by_property_id_ply(pid, 5), 'NO')
        with self.assertRaises(IndexError):
            p.get_sout_by_property_id_ply(pid, 6)

        #---------------

        #Rho
        self.assertAlmostEqual(p.get_density_by_property_id_ply(pid, 0), 1.0)
        self.assertAlmostEqual(p.get_density_by_property_id_ply(pid, 1), 1.0)
        self.assertAlmostEqual(p.get_density_by_property_id_ply(pid, 2), 1.0)
        self.assertAlmostEqual(p.get_density_by_property_id_ply(pid, 3), 1.0)
        self.assertAlmostEqual(p.get_density_by_property_id_ply(pid, 4), 1.0)
        self.assertAlmostEqual(p.get_density_by_property_id_ply(pid, 5), 1.0)
        with self.assertRaises(IndexError):
            p.get_density_by_property_id_ply(pid, 6)

        # MassPerArea
        #self.assertAlmostEqual(p.get_mass_per_area_by_property_id(), 1.2)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id(pid), 1.2)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 0), 0.1)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 1), 0.2)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 2), 0.3)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 3), 0.1)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 4), 0.2)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 5), 0.3)
        with self.assertRaises(IndexError):
            p.get_mass_per_area_by_property_id_ply(pid, 6)

        self.assertEqual(p.get_nonstructural_mass_by_property_id(pid), 0.0)
        #----------------------
        # change the nsm to 1.0
        p.set_nonstructural_mass_by_property_id(pid, 1.0)

        self.assertEqual(p.get_nonstructural_mass_by_property_id(pid), 1.0)
        # MassPerArea
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id(pid), 2.2)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 0, method='nplies'), 0.1+1/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 1, method='nplies'), 0.2+1/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 2, method='nplies'), 0.3+1/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 3, method='nplies'), 0.1+1/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 4, method='nplies'), 0.2+1/6.)
        self.assertAlmostEqual(p.get_mass_per_area_by_property_id_ply(pid, 5, method='nplies'), 0.3+1/6.)
        with self.assertRaises(IndexError):
            p.get_mass_per_area_by_property_id_ply(pid, 6)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
