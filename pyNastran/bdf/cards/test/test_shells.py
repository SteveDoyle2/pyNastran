from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import range
import unittest
from numpy import array

from pyNastran.bdf.bdf import PCOMP, MAT1, BDF

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
        cquad4 = model.Element(eid)
        node_ids = cquad4.node_ids
        assert node_ids == [n1, n2, n3, n4], node_ids

        # cquad4 / pshell
        self.assertEqual(cquad4.Eid(), eid)
        self.assertEqual(cquad4.Pid(), pid)
        self.assertEqual(cquad4.Mid(), mid)
        self.assertEqual(cquad4.Nsm(), nsm)
        self.assertEqual(cquad4.Mass(), mass)
        self.assertAlmostEqual(cquad4.MassPerArea(), mass / A)
        self.assertEqual(cquad4.Area(), A)
        self.assertEqual(cquad4.Thickness(), t)
        #self.assertEqual(cquad4.Rho(), rho)  # removed because of PCOMP

    def _make_ctria3(self, model, rho, nu, G, E, t, nsm):
        eid = 10
        pid = 20
        mid = 30
        n1 = 1
        n2 = 2
        n3 = 3
        mid2 = mid3 = mid4 = twelveIt3 = tst = z1 = z2 = None
        z0 = sb = ft = tref = ge = lam = None
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
            [
                'pcomp', pid + 1, z0, nsm, sb, ft, tref, ge, lam,
                mid, t, theta0, sout,
                mid, 2 * t, theta1, sout,
                mid, 3 * t, theta2, sout,
                mid, 4 * t, theta3, sout,
            ],
            ['mat1', mid, E, G, nu, rho],
        ]
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        # ctria3 / pshell
        ctria3 = model.Element(eid)
        node_ids = ctria3.node_ids
        assert node_ids == [n1, n2, n3], node_ids
        mass = A * (t * rho + nsm)
        self.assertEqual(ctria3.Eid(), eid)
        self.assertEqual(ctria3.Pid(), pid)
        self.assertEqual(ctria3.Mid(), mid)
        self.assertEqual(ctria3.Nsm(), nsm)
        self.assertEqual(ctria3.Mass(), mass)
        self.assertAlmostEqual(ctria3.MassPerArea(), mass / A)
        self.assertEqual(ctria3.Area(), A)
        self.assertEqual(ctria3.Thickness(), t)
        self.assertEqual(ctria3.MassPerArea(), mass / A)

        # removed because of PCOMP
        # also no E, G, J, Nu, for the same reason
        # what about Mid
        #self.assertEqual(ctria3.Rho(), rho)


        # pshell
        pshell = model.Property(pid)
        self.assertEqual(pshell.Pid(), pid)
        self.assertEqual(pshell.Mid(), mid)
        self.assertEqual(pshell.Nsm(), nsm)
        self.assertEqual(pshell.Thickness(), t)
        self.assertEqual(pshell.Rho(), rho)
        self.assertEqual(pshell.z1, -t / 2.)
        self.assertEqual(pshell.z2, t / 2.)

        # ctria3 / pcomp
        ctria3 = model.Element(eid + 1)
        mass = A * (10 * t * rho + nsm)
        self.assertEqual(ctria3.Eid(), eid + 1)
        self.assertEqual(ctria3.Pid(), pid + 1)
        #self.assertEqual(ctria3.Mid(), mid)
        self.assertEqual(ctria3.Nsm(), nsm)
        self.assertAlmostEqual(ctria3.Mass(), mass)
        self.assertAlmostEqual(ctria3.MassPerArea(), mass / A)
        self.assertEqual(ctria3.Area(), A)
        self.assertEqual(ctria3.Thickness(), 10 * t)
        #self.assertEqual(ctria3.Rho(), rho)

        # pcomp
        pcomp = model.Property(pid + 1)
        self.assertEqual(pcomp.Pid(), pid + 1)
        self.assertEqual(pcomp.nPlies(), 4)

        self.assertEqual(pcomp.Mid(0), mid)
        self.assertEqual(pcomp.Nsm(), nsm)

        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.Mid(-1), mid)
        self.assertEqual(pcomp.Mids(), [mid] * 4)
        self.assertEqual(pcomp.Mid(0), mid)
        self.assertEqual(pcomp.Mid(1), mid)
        self.assertEqual(pcomp.Mid(2), mid)
        self.assertEqual(pcomp.Mid(3), mid)
        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.Mid(4), mid)

        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.Thickness(-1), t)
        self.assertEqual(pcomp.Thickness(), 10 * t)
        self.assertEqual(pcomp.Thickness(0), t)
        self.assertEqual(pcomp.Thickness(1), 2 * t)
        self.assertAlmostEqual(pcomp.Thickness(2), 3 * t, places=8) # 0.3
        self.assertEqual(pcomp.Thickness(3), 4 * t)
        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.Thickness(4), 5*t)

        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.Rho(-1), rho)
        self.assertEqual(pcomp.Rho(0), rho)
        self.assertEqual(pcomp.Rho(1), rho)
        self.assertEqual(pcomp.Rho(2), rho)
        self.assertEqual(pcomp.Rho(3), rho)
        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.Rho(4), rho)

        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.Theta(-1), 0.)
        self.assertEqual(pcomp.Theta(0), 0.)
        self.assertEqual(pcomp.Theta(1), 30.)
        self.assertEqual(pcomp.Theta(2), 60.)
        self.assertEqual(pcomp.Theta(3), 90.)
        with self.assertRaises(IndexError):
            self.assertEqual(pcomp.Theta(4), rho)
        self.assertEqual(pcomp.z0, -10*t/2.)

    def test_PSHELL_01(self):
        """tests a CQUAD4 and a PSHELL"""

        rho = 0.1
        nu = 0.3
        G = None
        E = 1e7
        t = 0.3
        nsm = 0.0

        model = BDF(debug=False)
        self._make_cquad4(model, rho, nu, G, E, t, nsm)

        model = BDF(debug=False)
        self._make_ctria3(model, rho, nu, G, E, t, nsm)

        nsm = 1.0
        model = BDF(debug=False)
        self._make_cquad4(model, rho, nu, G, E, t, nsm)

        model = BDF(debug=False)
        self._make_ctria3(model, rho, nu, G, E, t, nsm)

    def test_CQUAD4_01(self):
        model = BDF(debug=False)
        eid = 10
        pid = 20
        mid = 30
        n1 = 1
        n2 = 2
        n3 = 3
        n4 = 4
        n5 = 5
        n6 = 6
        A = 2.
        t = rho = nsm = E = G = nu = 0.1
        mid2 = mid3 = mid4 = twelveIt3 = tst = z1 = z2 = None

        mass = A * (t * rho + nsm)
        cards = [
            ['grid', n1, 0, 0., 0., 0.],
            ['grid', n2, 0, 2., 0., 0.],
            ['grid', n3, 0, 2., 1., 0.],
            ['grid', n4, 0, 0., 1., 0.],
            ['grid', n5, 0, 0., 0., 0.],
            ['grid', n6, 0, 2., 0., 0.],

            ['cquad4', eid, pid, n1, n2, n3, n4],
            ['cquad4', eid+1, pid, n5, n6, n3, n4],
            ['pshell', pid, mid, t, mid2, twelveIt3, mid3, tst, nsm, z1, z2],
            ['mat1', mid, E, G, nu, rho],
        ]
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)

        # get node IDs without cross referencing
        eids = [10]
        nids = model.get_node_ids_with_elements(eids)
        assert nids == set([1, 2, 3, 4]), nids

        eids = [11]
        nids = model.get_node_ids_with_elements(eids)
        assert nids == set([3, 4, 5, 6]), nids

        eids = [10, 11]
        nids = model.get_node_ids_with_elements(eids)
        assert nids == set([1, 2, 3, 4, 5, 6]), nids

        # get node IDs with cross referencing
        model.cross_reference()
        eids = [10]
        nids = model.get_node_ids_with_elements(eids)
        assert nids == set([1, 2, 3, 4]), nids

        eids = [11]
        nids = model.get_node_ids_with_elements(eids)
        assert nids == set([3, 4, 5, 6]), nids

        eids = [10, 11]
        nids = model.get_node_ids_with_elements(eids)
        assert nids == set([1, 2, 3, 4, 5, 6]), nids



    def test_PCOMP_01(self):
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
        tref = 0.
        ge = 0.
        lam = 'NO' # isSymmetrical YES/NO
        Mid = [1, 2, 3]
        theta = [0., 10., 20.]
        T = [.1, .2, .3]
        sout = [1, 1, 0]  # 0-NO, 1-YES
        data = [pid, z0, nsm, sb, ft, tref, ge, lam, Mid, T, theta, sout]

        p = PCOMP.add_op2_data(data)
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

        self.assertEqual(p.Mids(), [1, 2, 3])

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
        mat1 = [mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, Mcsid]
        m = MAT1.add_op2_data(mat1)
        for iply in range(len(p.plies)):
            mid = p.plies[iply][0]
            p.mids[iply] = m # MAT1
            #p.mids = [m, m, m]
        p.mids_ref = p.mids

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

    def test_PCOMP_02(self):
        """
        symmetrical, nsm=0.0 and nsm=1.0
        """
        pid = 1
        z0 = 0.
        nsm = 0.
        sb = 0.
        ft = 0.
        tref = 0.
        ge = 0.
        lam = 'SYM'  # isSymmetrical SYM
        Mid = [1, 2, 3]
        theta = [0., 10., 20.]
        T = [.1, .2, .3]
        sout = [1, 1, 0]  # 0-NO, 1-YES
        data = [pid, z0, nsm, sb, ft, tref, ge, lam, Mid, T, theta, sout]
        p = PCOMP.add_op2_data(data)
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

        self.assertEqual(p.Mids(), [1, 2, 3, 1, 2, 3])

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
        mat1 = [mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, Mcsid]
        m = MAT1.add_op2_data(mat1)
        for iply in range(len(p.plies)):
            mid = p.plies[iply][0]
            p.mids[iply] = m # MAT1
        p.mids_ref = p.mids

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
