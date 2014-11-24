from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import range
import unittest
from numpy import array

from pyNastran.bdf.dev_vectorized.bdf import BDF
from pyNastran.bdf.dev_vectorized.cards.materials.mat1 import MAT1
from pyNastran.bdf.dev_vectorized.cards.elements.shell.pcomp import PCOMP
from pyNastran.bdf.dev_vectorized.cards.elements.shell.pshell import PSHELL

class TestShells(unittest.TestCase):
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
        model.cross_reference()
        cquad4 = model.elements[eid]

        # cquad4 / pshell
        self.assertEquals(cquad4.get_element_id_by_element_index(), eid)
        self.assertEquals(cquad4.get_property_id_by_element_index(), pid)
        #self.assertEquals(cquad4.Mid(), mid)
        #self.assertEquals(cquad4.Nsm(), nsm)
        self.assertEquals(cquad4.get_mass_by_element_id(), mass)
        self.assertAlmostEquals(cquad4.get_mass_per_area_by_element_id(), mass / A)
        self.assertEquals(cquad4.get_area_by_element_id(), A)
        self.assertEquals(cquad4.get_thickness_by_element_id(), t)
        #self.assertEquals(cquad4.Rho(), rho)  # removed because of PCOMP

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
            ['ctria3', eid, pid, n1, n2, n3],   # A = 1/2 * 4 * 1 = 2.
            ['pshell', pid, mid, t, mid2, twelveIt3, mid3, tst, nsm, z1, z2, mid4],

            ['ctria3', eid + 1, pid + 1, n1, n2, n3],   # A = 1/2 * 4 * 1 = 2.
            ['pcomp', pid + 1, z0, nsm, sb, ft, Tref, ge, lam,
                mid, t,     theta0, sout,
                mid, 2 * t, theta1, sout,
                mid, 3 * t, theta2, sout,
                mid, 4 * t, theta3, sout,
                ],
            ['mat1',mid, E, G, nu, rho],
        ]
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        # ctria3 / pshell
        ctria3 = model.elements[eid]
        mass = A * (t * rho + nsm)
        self.assertEquals(ctria3.get_element_id_by_element_index(), eid)
        self.assertEquals(ctria3.get_property_id_by_element_index(), pid)
        #self.assertEquals(ctria3.Mid(), mid)
        #self.assertEquals(ctria3.Nsm(), nsm)
        self.assertEquals(ctria3.get_mass_by_element_id(), mass)
        self.assertAlmostEquals(ctria3.get_mass_per_area_by_element_id(), mass / A)
        self.assertEquals(ctria3.get_area_by_element_id(), A)
        self.assertEquals(ctria3.get_thickness_by_element_id(), t)
        #self.assertEquals(ctria3.MassPerArea(), mass / A)

        # removed because of PCOMP
        # also no E, G, J, Nu, for the same reason
        # what about Mid
        #self.assertEquals(ctria3.Rho(), rho)


        # pshell
        pshell = model.properties[pid]
        assert isinstance(pshell, PSHELL), type(pshell)
        self.assertEquals(pshell.get_property_id_by_property_index(), pid)
        #self.assertEquals(pshell.Mid(), mid)
        #self.assertEquals(pshell.Nsm(), nsm)
        #self.assertEquals(pshell.Thickness(), t)
        #self.assertEquals(pshell.Rho(), rho)
        self.assertEquals(pshell.z1[0], -t / 2.)
        self.assertEquals(pshell.z2[0],  t / 2.)

        # ctria3 / pcomp
        ctria3 = model.elements[eid + 1]
        mass = A * (10 * t * rho + nsm)
        self.assertEquals(ctria3.get_element_id_by_element_index(), eid + 1)
        self.assertEquals(ctria3.get_property_id_by_element_id(), pid + 1)
        #self.assertEquals(ctria3.Mid(), mid)
        #self.assertEquals(ctria3.Nsm(), nsm)
        self.assertAlmostEquals(ctria3.get_mass_by_element_id(), mass)
        self.assertAlmostEquals(ctria3.get_mass_per_area_by_element_id(), mass / A)
        self.assertEquals(ctria3.get_area_by_element_id(), A)
        self.assertEquals(ctria3.get_thickness_by_element_id(), 10 * t)
        #self.assertEquals(ctria3.Rho(), rho)

        # pcomp
        pcomp = model.properties[pid + 1]
        pcomp = pcomp[pid + 1]
        print('pcomp =', pcomp)
        #self.assertEquals(pcomp.get_property_id()[0], pid + 1)
        self.assertEquals(pcomp.Pid(), pid + 1)
        self.assertEquals(pcomp.nPlies(), 4)

        self.assertEquals(pcomp.Mid(0), mid)
        self.assertEquals(pcomp.Nsm(), nsm)

        with self.assertRaises(IndexError):
            self.assertEquals(pcomp.Mid(-1), mid)
        self.assertEquals(pcomp.Mids(), [mid] * 4)
        self.assertEquals(pcomp.Mid(0), mid)
        self.assertEquals(pcomp.Mid(1), mid)
        self.assertEquals(pcomp.Mid(2), mid)
        self.assertEquals(pcomp.Mid(3), mid)
        with self.assertRaises(IndexError):
            self.assertEquals(pcomp.Mid(4), mid)

        with self.assertRaises(IndexError):
            self.assertEquals(pcomp.Thickness(-1), t)
        self.assertEquals(pcomp.Thickness(), 10 * t)
        self.assertEquals(pcomp.Thickness(0), t)
        self.assertEquals(pcomp.Thickness(1), 2 * t)
        self.assertEquals(pcomp.Thickness(2), 3 * t)
        self.assertEquals(pcomp.Thickness(3), 4 * t)
        with self.assertRaises(IndexError):
            self.assertEquals(pcomp.Thickness(4), 5*t)

        with self.assertRaises(IndexError):
            self.assertEquals(pcomp.Rho(-1), rho)
        self.assertEquals(pcomp.Rho(0), rho)
        self.assertEquals(pcomp.Rho(1), rho)
        self.assertEquals(pcomp.Rho(2), rho)
        self.assertEquals(pcomp.Rho(3), rho)
        with self.assertRaises(IndexError):
            self.assertEquals(pcomp.Rho(4), rho)

        with self.assertRaises(IndexError):
            self.assertEquals(pcomp.Theta(-1), 0.)
        self.assertEquals(pcomp.Theta(0), 0.)
        self.assertEquals(pcomp.Theta(1), 30.)
        self.assertEquals(pcomp.Theta(2), 60.)
        self.assertEquals(pcomp.Theta(3), 90.)
        with self.assertRaises(IndexError):
            self.assertEquals(pcomp.Theta(4), rho)
        self.assertEquals(pcomp.z0, -10*t/2.)

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
        print('starting BDF1')
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
        print('starting BDF2')
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
        print('starting BDF3')
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
        print('starting BDF4')
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
        #self.TRef = data[5]
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
        ft = 0.
        TRef = 0.
        ge = 0.
        lam = 'NO'  # isSymmetrical YES/NO
        Mid = [1,2,3]
        Theta = [0.,10.,20.]
        T = [.1,.2,.3]
        Sout = [1, 1, 0]  # 0-NO, 1-YES
        data = [pid, z0, nsm, sb, ft, TRef, ge, lam, Mid, T, Theta, Sout]

        p = PCOMP(data=data)
        self.assertFalse(p.isSymmetrical())
        self.assertEqual(p.nPlies(), 3)

        self.assertAlmostEqual(p.Thickness(), 0.6)
        self.assertAlmostEqual(p.Thickness(0), 0.1)
        self.assertAlmostEqual(p.Thickness(1), 0.2)
        self.assertAlmostEqual(p.Thickness(2), 0.3)
        with self.assertRaises(IndexError):
            p.Thickness(3)

        self.assertAlmostEqual(p.Theta(0), 0.)
        self.assertAlmostEqual(p.Theta(1), 10.)
        self.assertAlmostEqual(p.Theta(2), 20.)
        with self.assertRaises(IndexError):
            p.Theta(3)

        self.assertEqual(p.Mid(0), 1)
        self.assertEqual(p.Mid(1), 2)
        self.assertEqual(p.Mid(2), 3)
        with self.assertRaises(IndexError):
            p.Mid(3)

        self.assertEqual(p.Mids(), [1,2,3])

        self.assertEqual(p.sout(0), 'YES')
        self.assertEqual(p.sout(1), 'YES')
        self.assertEqual(p.sout(2), 'NO')
        with self.assertRaises(IndexError):
            p.sout(3)

        # material...
        #self.mid = data[0]
        #self.e = data[1]
        #self.g = data[2]
        #self.nu = data[3]
        #self.rho = data[4]
        #self.a = data[5]
        #self.TRef = data[6]
        #self.ge = data[7]
        #self.St = data[8]
        #self.Sc = data[9]
        #self.Ss = data[10]
        #self.Mcsid = data[11]
        mid = 1
        E = None
        G = None
        nu = None
        rho = 1.0
        a = None
        St = None
        Sc = None
        Ss = None
        Mcsid = None
        mat1 = [mid,E,G,nu,rho,a,TRef, ge, St, Sc, Ss, Mcsid]
        m = MAT1(data=mat1)
        for iply in range(len(p.plies)):
            mid = p.plies[iply][0]
            p.plies[iply][0] = m # MAT1
            #p.mids = [m, m, m]

        #Rho
        self.assertAlmostEqual(p.Rho(0), 1.0)
        self.assertAlmostEqual(p.Rho(1), 1.0)
        self.assertAlmostEqual(p.Rho(2), 1.0)
        with self.assertRaises(IndexError):
            p.Rho(3)

        # MassPerArea
        self.assertAlmostEqual(p.MassPerArea(), 0.6)
        self.assertAlmostEqual(p.MassPerArea(0), 0.1)
        self.assertAlmostEqual(p.MassPerArea(1), 0.2)
        self.assertAlmostEqual(p.MassPerArea(2), 0.3)
        with self.assertRaises(IndexError):
            p.MassPerArea(3)

        #----------------------
        # change the nsm to 1.0
        p.nsm = 1.0

        self.assertEqual(p.Nsm(), 1.0)
        # MassPerArea
        self.assertAlmostEqual(p.MassPerArea(), 1.6)
        self.assertAlmostEqual(p.MassPerArea(0, method='nplies'), 0.1+1/3.)
        self.assertAlmostEqual(p.MassPerArea(1, method='nplies'), 0.2+1/3.)
        self.assertAlmostEqual(p.MassPerArea(2, method='nplies'), 0.3+1/3.)

        self.assertAlmostEqual(p.MassPerArea(0, method='rho*t'), 0.1+1/6.)
        self.assertAlmostEqual(p.MassPerArea(1, method='rho*t'), 0.2+2/6.)
        self.assertAlmostEqual(p.MassPerArea(2, method='rho*t'), 0.3+3/6.)

        self.assertAlmostEqual(p.MassPerArea(0, method='t'), 0.1+1/6.)
        self.assertAlmostEqual(p.MassPerArea(1, method='t'), 0.2+2/6.)
        self.assertAlmostEqual(p.MassPerArea(2, method='t'), 0.3+3/6.)
        with self.assertRaises(IndexError):
            p.MassPerArea(3, method='nplies')

        z = p.get_z_locations()
        z_expected = array([0., T[0], T[0]+T[1], T[0]+T[1]+T[2]])
        for za, ze in zip(z, z_expected):
            self.assertAlmostEqual(za, ze)

        #z0  =
        p.z0 = 1.0
        z_expected = 1.0 + z_expected
        z = p.get_z_locations()
        for za, ze in zip(z, z_expected):
            self.assertAlmostEqual(za, ze)

    def test_pcomp_02(self):
        """
        symmetrical, nsm=0.0 and nsm=1.0
        """
        pid = 1
        z0 = 0.
        nsm = 0.
        sb = 0.
        ft = 0.
        TRef = 0.
        ge = 0.
        lam = 'SYM'  # isSymmetrical SYM
        Mid = [1,2,3]
        Theta = [0.,10.,20.]
        T = [.1,.2,.3]
        Sout = [1, 1, 0]  # 0-NO, 1-YES
        data = [pid, z0, nsm, sb, ft, TRef, ge, lam, Mid, T, Theta, Sout]
        p = PCOMP(data=data)
        self.assertTrue(p.isSymmetrical())
        self.assertEqual(p.nPlies(), 6)

        self.assertAlmostEqual(p.Thickness(), 1.2)
        self.assertAlmostEqual(p.Thickness(0), 0.1)
        self.assertAlmostEqual(p.Thickness(1), 0.2)
        self.assertAlmostEqual(p.Thickness(2), 0.3)
        self.assertAlmostEqual(p.Thickness(3), 0.1)
        self.assertAlmostEqual(p.Thickness(4), 0.2)
        self.assertAlmostEqual(p.Thickness(5), 0.3)
        with self.assertRaises(IndexError):
            p.Thickness(6)

        self.assertAlmostEqual(p.Theta(0), 0.)
        self.assertAlmostEqual(p.Theta(1), 10.)
        self.assertAlmostEqual(p.Theta(2), 20.)
        self.assertAlmostEqual(p.Theta(3), 0.)
        self.assertAlmostEqual(p.Theta(4), 10.)
        self.assertAlmostEqual(p.Theta(5), 20.)
        with self.assertRaises(IndexError):
            p.Theta(6)

        self.assertEqual(p.Mid(0), 1)
        self.assertEqual(p.Mid(1), 2)
        self.assertEqual(p.Mid(2), 3)
        self.assertEqual(p.Mid(3), 1)
        self.assertEqual(p.Mid(4), 2)
        self.assertEqual(p.Mid(5), 3)
        with self.assertRaises(IndexError):
            p.Mid(6)

        self.assertEqual(p.Mids(), [1,2,3,1,2,3])

        self.assertEqual(p.sout(0), 'YES')
        self.assertEqual(p.sout(1), 'YES')
        self.assertEqual(p.sout(2), 'NO')
        self.assertEqual(p.sout(3), 'YES')
        self.assertEqual(p.sout(4), 'YES')
        self.assertEqual(p.sout(5), 'NO')
        with self.assertRaises(IndexError):
            p.sout(6)


        mid = 1
        E = None
        G = None
        nu = None
        rho = 1.0
        a = None
        St = None
        Sc = None
        Ss = None
        Mcsid = None
        mat1 = [mid,E,G,nu,rho,a,TRef, ge, St, Sc, Ss, Mcsid]
        m = MAT1(data=mat1)
        for iply in range(len(p.plies)):
            mid = p.plies[iply][0]
            p.plies[iply][0] = m # MAT1

        #Rho
        self.assertAlmostEqual(p.Rho(0), 1.0)
        self.assertAlmostEqual(p.Rho(1), 1.0)
        self.assertAlmostEqual(p.Rho(2), 1.0)
        self.assertAlmostEqual(p.Rho(3), 1.0)
        self.assertAlmostEqual(p.Rho(4), 1.0)
        self.assertAlmostEqual(p.Rho(5), 1.0)
        with self.assertRaises(IndexError):
            p.Rho(6)

        # MassPerArea
        self.assertAlmostEqual(p.MassPerArea(), 1.2)
        self.assertAlmostEqual(p.MassPerArea(0), 0.1)
        self.assertAlmostEqual(p.MassPerArea(1), 0.2)
        self.assertAlmostEqual(p.MassPerArea(2), 0.3)
        self.assertAlmostEqual(p.MassPerArea(3), 0.1)
        self.assertAlmostEqual(p.MassPerArea(4), 0.2)
        self.assertAlmostEqual(p.MassPerArea(5), 0.3)
        with self.assertRaises(IndexError):
            p.MassPerArea(6)

        self.assertEqual(p.Nsm(), 0.0)
        #----------------------
        # change the nsm to 1.0
        p.nsm = 1.0

        self.assertEqual(p.Nsm(), 1.0)
        # MassPerArea
        self.assertAlmostEqual(p.MassPerArea(), 2.2)
        self.assertAlmostEqual(p.MassPerArea(0, method='nplies'), 0.1+1/6.)
        self.assertAlmostEqual(p.MassPerArea(1, method='nplies'), 0.2+1/6.)
        self.assertAlmostEqual(p.MassPerArea(2, method='nplies'), 0.3+1/6.)
        self.assertAlmostEqual(p.MassPerArea(3, method='nplies'), 0.1+1/6.)
        self.assertAlmostEqual(p.MassPerArea(4, method='nplies'), 0.2+1/6.)
        self.assertAlmostEqual(p.MassPerArea(5, method='nplies'), 0.3+1/6.)
        with self.assertRaises(IndexError):
            p.MassPerArea(6)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()