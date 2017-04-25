from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import unittest
from six.moves import range
from numpy import array

from pyNastran.bdf.bdf import PCOMP, MAT1, BDF
from pyNastran.bdf.cards.test.utils import save_load_deck


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
        z0_elem = 0.1
        mid2 = mid3 = mid4 = theta_mcid = twelveIt3 = tst = z1 = z2 = None
        #z0_prop =  None

        mass = A * (t * rho + nsm)
        cards = [
            ['grid', n1, 0, 0., 0., 0.],
            ['grid', n2, 0, 2., 0., 0.],
            ['grid', n3, 0, 2., 1., 0.],
            ['grid', n4, 0, 0., 1., 0.],
            ['cquad4', eid, pid, n1, n2, n3, n4, theta_mcid, z0_elem],
            ['pshell', pid, mid, t, mid2, twelveIt3, mid3, tst, nsm, z1, z2],
            ['mat1', mid, E, G, nu, rho],
        ]
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.validate()
        model._verify_bdf(xref=False)
        model.mass_properties_no_xref()

        model.cross_reference()
        model.mass_properties()
        model._verify_bdf(xref=True)
        cquad4 = model.Element(eid)
        pshell = model.Property(pid)
        node_ids = cquad4.node_ids
        assert node_ids == [n1, n2, n3, n4], node_ids

        # cquad4 / pshell
        self.assertEqual(cquad4.eid, eid)
        self.assertEqual(cquad4.Pid(), pid)
        self.assertEqual(cquad4.Mid(), mid)
        self.assertEqual(cquad4.Nsm(), nsm)
        self.assertEqual(cquad4.Mass(), mass)
        self.assertAlmostEqual(cquad4.MassPerArea(), mass / A)
        self.assertEqual(cquad4.Area(), A)
        self.assertEqual(cquad4.Thickness(), t)
        self.assertEqual(cquad4.zoffset, z0_elem)
        self.assertEqual(pshell.z1, -t/2.)
        #self.assertEqual(cquad4.Rho(), rho)  # removed because of PCOMP

    def _make_ctria3(self, model, rho, nu, G, E, t, nsm):
        eid = 10
        pid = 20
        mid = 30
        n1 = 1
        n2 = 2
        n3 = 3
        mid2 = mid3 = mid4 = theta_mcid = twelveIt3 = tst = z1 = z2 = None
        z0_elem = 0.1

        z0_prop = sb = ft = tref = ge = lam = None
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
            ['ctria3', eid, pid, n1, n2, n3, theta_mcid, z0_elem],   # A = 1/2 * 4 * 1 = 2.
            ['pshell', pid, mid, t, mid2, twelveIt3, mid3, tst, nsm, z1, z2, mid4],

            ['ctria3', eid + 1, pid + 1, n1, n2, n3, theta_mcid, z0_elem],   # A = 1/2 * 4 * 1 = 2.
            [
                'pcomp', pid + 1, z0_prop, nsm, sb, ft, tref, ge, lam,
                mid, t, theta0, sout,
                mid, 2 * t, theta1, sout,
                mid, 3 * t, theta2, sout,
                mid, 4 * t, theta3, sout,
            ],
            ['mat1', mid, E, G, nu, rho],
        ]
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.validate()
        model._verify_bdf(xref=False)
        model.cross_reference()
        model._verify_bdf(xref=True)

        # ctria3 / pshell
        ctria3 = model.Element(eid)
        node_ids = ctria3.node_ids
        assert node_ids == [n1, n2, n3], node_ids
        mass = A * (t * rho + nsm)
        self.assertEqual(ctria3.eid, eid)
        self.assertEqual(ctria3.Pid(), pid)
        self.assertEqual(ctria3.Mid(), mid)
        self.assertEqual(ctria3.Nsm(), nsm)
        self.assertEqual(ctria3.Mass(), mass)
        self.assertAlmostEqual(ctria3.MassPerArea(), mass / A)
        self.assertEqual(ctria3.Area(), A)
        self.assertEqual(ctria3.Thickness(), t)
        self.assertEqual(ctria3.MassPerArea(), mass / A)
        self.assertEqual(ctria3.zoffset, z0_elem)
        ctria3.raw_fields()

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
        self.assertEqual(ctria3.eid, eid + 1)
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

    def test_pshell_01(self):
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

    def test_cquad4_01(self):
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
        save_load_deck(model)


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
        #self.tref = data[6]
        #self.ge = data[7]
        #self.St = data[8]
        #self.Sc = data[9]
        #self.Ss = data[10]
        #self.mcsid = data[11]
        mid = 1
        E = None
        G = None
        nu = None
        rho = 1.0
        a = None
        St = None
        Sc = None
        Ss = None
        mcsid = None
        mat1 = [mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
        with self.assertRaises(ValueError):
            m = MAT1.add_op2_data(mat1)

        G = 42.
        mat1 = [mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
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

    def test_pcomp_02(self):
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
        mcsid = None
        mat1 = [mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
        with self.assertRaises(ValueError):
            m = MAT1.add_op2_data(mat1)

        G = 42.
        mat1 = [mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
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

    def test_cshear(self):
        """tests a PSHEAR/CSHEAR"""
        model = BDF(debug=False)
        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(3, xyz=[1., 1., 0.])
        model.add_grid(4, xyz=[0., 1., 0.])

        eid = 10
        pid = 20
        mid = 30
        t = 0.1
        nids = [1, 2, 3, 4]

        cshear = model.add_cshear(eid, pid, nids, comment='cshear')
        pshear = model.add_pshear(pid, mid, t, nsm=0., f1=0., f2=0., comment='')

        E = 30.e7
        G = None
        nu = 0.3
        mat1 = model.add_mat1(mid, E, G, nu, rho=0.1, comment='mat1')
        model.pop_parse_errors()
        model.validate()

        cshear.raw_fields()
        cshear.write_card(size=8)

        pshear.raw_fields()
        pshear.write_card(size=8)
        pshear.write_card(size=16)
        pshear.write_card(size=16, is_double=True)

        model.validate()
        model._verify_bdf(xref=False)
        model.cross_reference()
        model._verify_bdf(xref=True)
        model.mass_properties()

        cshear.write_card(size=8)
        pshear.write_card(size=8)
        #save_load_deck(model)

    def test_shells(self):
        """tests a CTRIA3/CQUAD4/PSHELL and CTRIA6/CQUAD8/CQUAD/PCOMP"""
        model = BDF(debug=False)
        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(3, xyz=[1., 1., 0.])
        model.add_grid(4, xyz=[0., 1., 0.])

        model.add_grid(5, xyz=[.5, 0., 0.])
        model.add_grid(6, xyz=[1., 0.5, 0.])
        model.add_grid(7, xyz=[.5, 1., 0.])
        model.add_grid(8, xyz=[0., .5, 0.])

        model.add_grid(9, xyz=[.5, .5, 0.])

        E = 30.e7
        G = None
        nu = 0.3
        model.add_mat1(1, E, G, nu, rho=0.1)
        model.add_mat1(2, E, G, nu, rho=0.1)
        model.add_mat1(3, E, G, nu, rho=0.1)

        pid = 1
        nids = [1, 2, 3]
        model.add_ctria3(1, pid, nids)
        nids = [1, 2, 3, 4]
        model.add_cquad4(2, pid, nids)
        model.add_pshell(pid, mid1=2, t=0.1)

        pid = 2
        nids = [1, 2, 3, 5, 6, 9]
        ctria6 = model.add_ctria6(3, pid, nids, comment='ctria6')

        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        cquad8 = model.add_cquad8(4, pid, nids, comment='cquad8')

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        cquad = model.add_cquad(5, pid, nids, comment='cquad')

        mids = [1, 2, 3]
        thicknesses = [0.1, 0.2, 0.3]
        pcomp = model.add_pcomp(pid, mids, thicknesses)

        assert pcomp.Thickness() == sum(thicknesses), thicknesses

        pcomp.lam = 'SYM'
        assert pcomp.Thickness() == sum(thicknesses)*2, thicknesses

        model.validate()

        ctria6.raw_fields()
        ctria6.write_card(size=8)

        cquad8.raw_fields()
        cquad8.write_card(size=8)

        cquad.raw_fields()
        cquad.write_card(size=8)

        pcomp.raw_fields()
        pcomp.write_card(size=8)
        pcomp.write_card(size=16)
        pcomp.write_card(size=16, is_double=True)


        model._verify_bdf(xref=False)
        #--------------------------------
        model.cross_reference()
        model._verify_bdf(xref=True)

        ctria6.raw_fields()
        ctria6.write_card(size=8)

        cquad8.raw_fields()
        cquad8.write_card(size=8)

        cquad.raw_fields()
        cquad.write_card(size=8)

        pcomp.raw_fields()
        pcomp.write_card(size=8)
        pcomp.write_card(size=16)
        pcomp.write_card(size=16, is_double=True)

    def test_trax(self):
        """tests a CTRAX3/CTRAX6/???"""
        model = BDF(debug=False)
        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(3, xyz=[1., 1., 0.])
        model.add_grid(4, xyz=[0., 1., 0.])

        model.add_grid(5, xyz=[.5, 0., 0.])
        model.add_grid(6, xyz=[1., 0.5, 0.])
        model.add_grid(7, xyz=[.5, 1., 0.])
        model.add_grid(8, xyz=[0., .5, 0.])

        model.add_grid(9, xyz=[.5, .5, 0.])

        mid1 = 1
        E = 30.e7
        G = None
        nu = 0.3
        model.add_mat1(mid1, E, G, nu, rho=0.1)
        #model.add_mat1(2, E, G, nu, rho=0.1)
        #model.add_mat1(3, E, G, nu, rho=0.1)

        pid = 1
        nids = [1, 2, 3]
        ctrax3 = model.add_ctrax3(1, pid, nids, theta=0., comment='ctrax3')
        #model.add_pshell(pid, mid1=2, t=0.1)

        psolid = model.add_psolid(pid, mid1, cordm=0, integ=None, stress=None,
                                  isop=None, fctn='SMECH', comment='psolid')

        pid = 2
        nids = [1, 2, 3, 5, 6, 9]
        ctrax6 = model.add_ctrax6(2, pid, nids, theta=0., comment='ctrax6')

        plsolid = model.add_plsolid(pid, mid1, stress_strain='GRID', ge=0.,
                                    comment='plsolid')

        #assert pcomp.Thickness() == sum(thicknesses), thicknesses

        #pcomp.lam = 'SYM'
        #assert pcomp.Thickness() == sum(thicknesses)*2, thicknesses

        model.validate()

        ctrax6.raw_fields()
        ctrax6.write_card(size=8)


        psolid.raw_fields()
        psolid.write_card(size=8)
        #psolid.write_card(size=16)
        #psolid.write_card(size=16, is_double=True)

        plsolid.raw_fields()
        plsolid.write_card(size=8)
        #plsolid.write_card(size=16)
        #plsolid.write_card(size=16, is_double=True)

        model._verify_bdf(xref=False)

        #--------------------------------
        model.cross_reference()
        model._verify_bdf(xref=True)

        ctrax3.raw_fields()
        ctrax3.write_card(size=8)

        ctrax6.raw_fields()
        ctrax6.write_card(size=8)

        #pcomp.raw_fields()
        #pcomp.write_card(size=8)
        #pcomp.write_card(size=16)
        #pcomp.write_card(size=16, is_double=True)

    def test_ctriar_cquadr(self):
        """tests a CTRIAR/PSHELL/MAT8"""
        model = BDF(debug=False)
        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(3, xyz=[1., 1., 0.])
        model.add_grid(4, xyz=[0., 1., 0.])
        eid = 6
        pid = 13
        nids = [1, 2, 3]
        ctriar = model.add_ctriar(eid, pid, nids, comment='ctriar')
        ctriar.raw_fields()
        ctriar.write_card(size=8, is_double=False)
        ctriar.write_card(size=16, is_double=False)
        ctriar.flipNormal()

        eid = 8
        nids = [1, 2, 3, 4]
        cquadr = model.add_cquadr(eid, pid, nids, comment='cquadr')
        cquadr.raw_fields()
        cquadr.write_card(size=8, is_double=False)
        cquadr.write_card(size=16, is_double=False)
        cquadr.flipNormal()

        mid = 42
        pshell = model.add_pshell(pid, mid1=mid, t=0.2)
        e11 = 1e7
        e22 = 1e6
        nu12 = 0.3
        mat8 = model.add_mat8(mid, e11, e22, nu12)
        model.validate()
        model._verify_bdf(xref=False)
        model.cross_reference()
        model._verify_bdf(xref=True)
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(model)

    def test_cplstn34(self):
        """tests a CPLSTN3, CPLSTN4/PSHELL/MAT8"""
        model = BDF(debug=False)
        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(3, xyz=[1., 1., 0.])
        model.add_grid(4, xyz=[0., 1., 0.])
        pid = 4
        eid = 3
        nids = [1, 2, 3, 4]
        cplstn4 = model.add_cplstn4(eid, pid, nids, comment='cplstn4')

        eid = 5
        nids = [1, 2, 3]
        mid = 10
        cplstn3 = model.add_cplstn3(eid, pid, nids, comment='cplstn3')
        pplane = model.add_pplane(pid, mid, t=0.1, nsm=0.,
                                  formulation_option=0, comment='pplane')
        E = 1e7
        G = None
        nu = 0.3
        mat1 = model.add_mat1(mid, E, G, nu)

        cplstn3.raw_fields()
        cplstn4.raw_fields()
        pplane.raw_fields()

        model.validate()
        model._verify_bdf(xref=False)
        cplstn3.write_card(size=8)
        cplstn4.write_card(size=8)
        pplane.write_card(size=8)
        model.cross_reference()
        model.pop_xref_errors()
        #cplstn3.write_card(size=8)
        #cplstn4.write_card(size=8)

        model.uncross_reference()
        model.safe_cross_reference()

    def test_cplstn68(self):
        """tests a CPLSTN6, CPLSTN8/PSHELL/MAT8"""
        model = BDF(debug=False)
        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(5, xyz=[.5, 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(6, xyz=[1., .5, 0.])
        model.add_grid(3, xyz=[1., 1., 0.])
        model.add_grid(7, xyz=[.5, 1., 0.])
        model.add_grid(4, xyz=[0., 1., 0.])
        model.add_grid(8, xyz=[0., .5, 0.])
        pid = 4
        eid = 3
        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        cplstn8 = model.add_cplstn8(eid, pid, nids, comment='cplstn8')

        eid = 5
        nids = [1, 2, 3, 4, 5, 6]
        mid = 10
        cplstn6 = model.add_cplstn6(eid, pid, nids, comment='cplstn6')
        pplane = model.add_pplane(pid, mid, t=0.1, nsm=0.,
                                  formulation_option=0, comment='pplane')
        E = 1e7
        G = None
        nu = 0.3
        mat1 = model.add_mat1(mid, E, G, nu)

        cplstn6.raw_fields()
        cplstn8.raw_fields()
        pplane.raw_fields()

        model.validate()
        model._verify_bdf(xref=False)
        cplstn6.write_card(size=8)
        cplstn8.write_card(size=8)
        pplane.write_card(size=8)
        model.cross_reference()
        model.pop_xref_errors()
        #cplstn3.write_card(size=8)
        #cplstn4.write_card(size=8)

        model.uncross_reference()
        model.safe_cross_reference()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
