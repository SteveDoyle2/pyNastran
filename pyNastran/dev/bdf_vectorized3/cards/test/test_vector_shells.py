"""defines various shell element tests"""
import os
from io import StringIO
import unittest
import numpy as np
from numpy import array

from cpylog import get_logger
from pyNastran.dev.bdf_vectorized3.bdf import BDF # PCOMP, MAT1, BDF, CTRIA3, CQUAD4
#from pyNastran.dev.bdf_vectorized3.cards.materials import MAT1
#from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
#from pyNastran.bdf.cards.materials import get_mat_props_S
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck
#from pyNastran.bdf.mesh_utils.mass_properties import (
    #mass_properties, mass_properties_no_xref, mass_properties_nsm)
from pyNastran.dev.bdf_vectorized3.cards.elements.shell_utils import tri_volume, quad_volume

#try:
    #import matplotlib
    #IS_MATPLOTLIB = True
    #from pyNastran.bdf.cards.elements.plot import plot_equivalent_lamina_vs_theta
#except ImportError:
IS_MATPLOTLIB = False


class TestAcoustic(unittest.TestCase):
    def test_chacbr_pacbar(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        chacbr = model.chacbr
        pacbar = model.pacbar

        eid = 2
        pid = 10
        nodes = [1, 2, 3, 4, 5, 6, 7, 8]
        idi = model.add_chacbr(eid, pid, nodes, comment='caabsf') # solid

        mback = 1.0
        mseptm = 2.
        freson = 3.
        kreson = 4.
        model.add_pacbar(pid, mback, mseptm, freson, kreson, comment='')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        model.add_grid(5, [0., 0., 0.])
        model.add_grid(6, [1., 0., 0.])
        model.add_grid(7, [1., 1., 0.])
        model.add_grid(8, [0., 1., 0.])
        #model.add_paabsf(pid, table_reactance_real=None, table_reactance_imag=None,
                         #s=1.0, a=1.0, b=0.0, k=0.0, rhoc=1.0, comment='')
        model.setup()

        assert len(chacbr.write(size=8)) > 0
        assert len(chacbr.write(size=16)) > 0
        assert len(pacbar.write(size=8)) > 0
        assert len(pacbar.write(size=16)) > 0

        save_load_deck(model, run_mass_properties=False, run_equivalence=False)

    def test_caabsf_paabsf(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        caabsf = model.caabsf
        paabsf = model.paabsf

        eid = 2
        pid = 10
        nodes = [1, 2, 3, 4]
        idi = model.add_caabsf(eid, pid, nodes, comment='caabsf')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_paabsf(pid, table_reactance_real=None, table_reactance_imag=None,
                         s=1.0, a=1.0, b=0.0, k=0.0, rhoc=1.0, comment='')
        model.setup()

        assert len(caabsf.write(size=8)) > 0
        assert len(caabsf.write(size=16)) > 0
        assert len(paabsf.write(size=8)) > 0
        assert len(paabsf.write(size=16)) > 0

        save_load_deck(model)

class TestShells(unittest.TestCase):
    def test_pshell(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        pid = 10
        model.add_pshell(pid, mid1=1, mid2=2, mid3=3, mid4=4, tst=3.14)
        mid = 1
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.setup()
        assert ' 3.14' in model.pshell.rstrip(), model.pshell.rstrip()

    def _make_cquad4(self, model, rho, nu, G, E, t, nsm):
        """CQUAD4/PSHELL/PCOMP"""
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
        model.setup()
        model.validate()
        model._verify_bdf(xref=False)
        model.mass()
        #mass_properties_no_xref(model)

        #model.cross_reference()
        #mass_properties(model)
        model._verify_bdf(xref=True)
        cquad4 = model.cquad4.slice_card_by_element_id(eid)
        cquad4.center_of_mass()
        pshell = model.pshell.slice_card_by_property_id(pid)
        node_ids = cquad4.nodes
        assert np.array_equal(node_ids[0], [n1, n2, n3, n4]), node_ids

        # cquad4 / pshell
        self.assertEqual(cquad4.element_id, eid)
        self.assertEqual(cquad4.property_id, pid)
        self.assertEqual(cquad4.mass(), mass)
        self.assertAlmostEqual(cquad4.mass_per_area(), mass / A)
        self.assertEqual(cquad4.area(), A)
        self.assertEqual(cquad4.total_thickness(), t)
        self.assertEqual(cquad4.zoffset, z0_elem)
        self.assertEqual(pshell.z[0, 0], -t/2.)

        #self.assertEqual(cquad4.Mid(), mid)
        self.assertEqual(cquad4.nsm_per_area(), nsm)
        cquad4.get_edge_axes()
        #self.assertEqual(cquad4.Rho(), rho)  # removed because of PCOMP

    def _make_ctria3(self, model, rho, nu, G, E, t, nsm):
        """CTRIA3/PSHELL/PCOMP"""
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
        model.setup()
        model.validate()
        model._verify_bdf(xref=False)
        #model.cross_reference()
        model._verify_bdf(xref=True)

        # ctria3 / pshell
        ctria3 = model.ctria3.slice_card_by_element_id(eid)
        ctria3.center_of_mass()
        assert ctria3.n > 0

        node_ids = ctria3.nodes[0, :]
        assert np.array_equal(node_ids, [n1, n2, n3]), node_ids
        mass = A * (t * rho + nsm)
        self.assertEqual(ctria3.element_id, eid)
        self.assertEqual(ctria3.property_id, pid)
        self.assertEqual(ctria3.mass(), mass)
        self.assertAlmostEqual(ctria3.mass_per_area(), mass / A)
        self.assertEqual(ctria3.area(), A)
        thicknessi = ctria3.total_thickness()
        assert np.allclose(thicknessi, t)
        self.assertEqual(ctria3.mass_per_area(), mass / A)
        self.assertEqual(ctria3.zoffset, z0_elem)
        #self.assertEqual(ctria3.mass_material_id(), mid)
        #self.assertEqual(ctria3.Nsm(), nsm)
        #ctria3.get_edge_axes()
        #ctria3.raw_fields()

        # removed because of PCOMP
        # also no E, G, J, Nu, for the same reason
        # what about Mid
        #self.assertEqual(ctria3.Rho(), rho)


        # pshell
        pshell = model.Property(pid)[0]
        self.assertEqual(pshell.property_id, pid)
        self.assertEqual(pshell.mass_material_id(), mid)
        self.assertEqual(pshell.nsm, nsm)
        self.assertEqual(pshell.total_thickness(), t)
        self.assertEqual(pshell.density(), rho)
        self.assertEqual(pshell.z1, -t / 2.)
        self.assertEqual(pshell.z2, t / 2.)

        # ctria3 / pcomp
        ctria3 = model.Element(eid + 1)[0]
        mass = A * (10 * t * rho + nsm)
        self.assertEqual(ctria3.element_id, eid + 1)
        self.assertEqual(ctria3.property_id, pid + 1)
        #self.assertEqual(ctria3.Mid(), mid)
        #self.assertEqual(ctria3.Nsm(), nsm)
        self.assertAlmostEqual(ctria3.mass().squeeze(), mass)
        self.assertAlmostEqual(ctria3.mass_per_area().squeeze(), mass / A)
        self.assertEqual(ctria3.area().squeeze(), A)
        self.assertEqual(ctria3.total_thickness().squeeze(), 10 * t)
        #self.assertEqual(ctria3.Rho(), rho)

        # pcomp
        pcomp: PCOMP = model.Property(pid + 1)[0]
        self.assertEqual(pcomp.property_id, pid + 1)
        self.assertEqual(pcomp.nplies, 4)
        self.assertEqual(pcomp.nplies_total, 4)

        #self.assertEqual(pcomp.Mid(0), mid)
        #self.assertEqual(pcomp.Nsm(), nsm)

        #with self.assertRaises(IndexError):
            #self.assertEqual(pcomp.Mid(-1), mid)
        #self.assertEqual(pcomp.Mids(), [mid] * 4)
        #self.assertEqual(pcomp.Mid(0), mid)
        #self.assertEqual(pcomp.Mid(1), mid)
        #self.assertEqual(pcomp.Mid(2), mid)
        #self.assertEqual(pcomp.Mid(3), mid)
        #with self.assertRaises(IndexError):
            #self.assertEqual(pcomp.Mid(4), mid)

        #with self.assertRaises(IndexError):
            #self.assertEqual(pcomp.Thickness(-1), t)
        #self.assertEqual(pcomp.Thickness(), 10 * t)
        #self.assertEqual(pcomp.Thickness(0), t)
        #self.assertEqual(pcomp.Thickness(1), 2 * t)
        #self.assertAlmostEqual(pcomp.Thickness(2), 3 * t, places=8) # 0.3
        #self.assertEqual(pcomp.Thickness(3), 4 * t)
        #with self.assertRaises(IndexError):
            #self.assertEqual(pcomp.Thickness(4), 5*t)

        #with self.assertRaises(IndexError):
            #self.assertEqual(pcomp.Rho(-1), rho)
        #self.assertEqual(pcomp.Rho(0), rho)
        #self.assertEqual(pcomp.Rho(1), rho)
        #self.assertEqual(pcomp.Rho(2), rho)
        #self.assertEqual(pcomp.Rho(3), rho)
        #with self.assertRaises(IndexError):
            #self.assertEqual(pcomp.Rho(4), rho)

        #with self.assertRaises(IndexError):
            #self.assertEqual(pcomp.Theta(-1), 0.)
        #self.assertEqual(pcomp.Theta(0), 0.)
        #self.assertEqual(pcomp.Theta(1), 30.)
        #self.assertEqual(pcomp.Theta(2), 60.)
        #self.assertEqual(pcomp.Theta(3), 90.)
        #with self.assertRaises(IndexError):
            #self.assertEqual(pcomp.Theta(4), rho)
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

    def test_pshell_02_tflag0(self):
        """"""
        model = BDF(debug=True, log=None, mode='msc')
        card_lines = ["PSHELL  1       1       .050    1               1"]
        model.add_card_lines(card_lines, 'PSHELL', comment='', has_none=True)
        card_lines = [
            "CQUAD4  1       1       1       2       10      9                       +M00000",
            "+M00000                  0.0     0.0     .041    .041",
        ]
        model.add_card_lines(card_lines, 'CQUAD4', comment='', has_none=True)

        mid = 1
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, alpha=0.0, tref=0.0,
                       ge=0.0, St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(10, [1., 1., 0.])
        model.add_grid(9, [0., 1., 0.])

        model.setup(run_geom_check=True)
        cquad4 = model.cquad4
        pshell = model.pshell
        assert np.allclose(pshell.t, 0.05)

        #The continuation is optional. If it is not supplied, then
        # T1 through T4 will be set equal to the value of T on the
        # PSHELL entry.

        #print(cquad4.get_stats())
        total_thickness = cquad4.total_thickness()
        #print('total_thickness=', total_thickness)
        #tlow = 0.041 / 2
        tavg = (0.041 + 0.05) / 2.
        #T      : array([[0., 0., 0.041, 0.041]])
        #tflag  : array([0], dtype=int8)
        assert np.allclose(tavg, total_thickness), (tavg, total_thickness)

    def test_cquad4_01(self):
        """CQUAD4/PSHELL"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 10
        pid = 20
        mid = 30
        n1 = 1
        n2 = 2
        n3 = 3
        n4 = 4
        n5 = 5
        n6 = 6
        #A = 2.
        t = rho = nsm = E = G = nu = 0.1
        mid2 = mid3 = mid4 = twelveIt3 = tst = z1 = z2 = None

        #mass = A * (t * rho + nsm)
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
        # nids = model.get_node_ids_with_elements(eids)
        # assert nids == {1, 2, 3, 4}, nids

        eids = [11]
        # nids = model.get_node_ids_with_elements(eids)
        # assert nids == {3, 4, 5, 6}, nids

        eids = [10, 11]
        # nids = model.get_node_ids_with_elements(eids)
        # assert nids == {1, 2, 3, 4, 5, 6}, nids

        params = [
            ('T', 1.0),
            (6, 2.0), # 12I/T3
            (8, 3.0), # 'TST'
        ]
        make_dvprel_optimization(model, params, 'PSHELL', pid)

        # get node IDs with cross referencing
        model.cross_reference()
        #model.update_model_by_desvars(xref=True)

        if 0:
            eids = [10]
            nids = model.get_node_ids_with_elements(eids)
            assert nids == {1, 2, 3, 4}, nids

            eids = [11]
            nids = model.get_node_ids_with_elements(eids)
            assert nids == {3, 4, 5, 6}, nids

            eids = [10, 11]
            nids = model.get_node_ids_with_elements(eids)
            assert nids == {1, 2, 3, 4, 5, 6}, nids
        save_load_deck(model)

    def test_pcomp_01(self):
        """asymmetrical, nsm=0.0 and nsm=1.0"""
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
        lam = 'NO' # is_symmetrical YES/NO
        Mid = [1, 2, 3]
        theta = [0., 10., 20.]
        T = [.1, .2, .3]
        sout = [1, 1, 0]  # 0-NO, 1-YES
        data = [pid, z0, nsm, sb, ft, tref, ge, lam, Mid, T, theta, sout]

        model = BDF(debug=True, log=None, mode='msc')
        p = model.pcomp
        p.add_op2_data(data)
        p.parse_cards()
        self.assertFalse(p.is_symmetrical)
        self.assertEqual(p.nplies, 3)
        self.assertEqual(p.nplies_total, 3)

        #self.assertAlmostEqual(p.Thickness(), 0.6)
        #self.assertAlmostEqual(p.Thickness(0), 0.1)
        #self.assertAlmostEqual(p.Thickness(1), 0.2)
        #self.assertAlmostEqual(p.Thickness(2), 0.3)
        #with self.assertRaises(IndexError):
            #p.Thickness(3)

        #self.assertAlmostEqual(p.Theta(0), 0.)
        #self.assertAlmostEqual(p.Theta(1), 10.)
        #self.assertAlmostEqual(p.Theta(2), 20.)
        #with self.assertRaises(IndexError):
            #p.Theta(3)

        #self.assertEqual(p.Mid(0), 1)
        #self.assertEqual(p.Mid(1), 2)
        #self.assertEqual(p.Mid(2), 3)
        #with self.assertRaises(IndexError):
            #p.Mid(3)

        #self.assertEqual(p.Mids(), [1, 2, 3])

        #self.assertEqual(p.sout(0), 'YES')
        #self.assertEqual(p.sout(1), 'YES')
        #self.assertEqual(p.sout(2), 'NO')
        #with self.assertRaises(IndexError):
            #p.sout(3)

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
        mat1 = model.mat1
        mat1_fields = [mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
        with self.assertRaises(AssertionError):
            m = mat1.add_op2_data(mat1_fields)
            #mat1.parse_cards()

        G = 42.
        mcsid = 0
        mat1_fields = [mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
        m = mat1.add_op2_data(mat1_fields)
        #for iply in range(len(p.plies)):
            #mid = p.plies[iply][0]
            #p.mids[iply] = m # MAT1
            ##p.mids = [m, m, m]
        #p.mids_ref = p.mids

        ##Rho
        #self.assertAlmostEqual(p.Rho(0), 1.0)
        #self.assertAlmostEqual(p.Rho(1), 1.0)
        #self.assertAlmostEqual(p.Rho(2), 1.0)
        #with self.assertRaises(IndexError):
            #p.Rho(3)

        ## MassPerArea
        #self.assertAlmostEqual(p.MassPerArea(), 0.6)
        #self.assertAlmostEqual(p.MassPerArea(0), 0.1)
        #self.assertAlmostEqual(p.MassPerArea(1), 0.2)
        #self.assertAlmostEqual(p.MassPerArea(2), 0.3)
        #with self.assertRaises(IndexError):
            #p.MassPerArea(3)

        #----------------------
        # change the nsm to 1.0
        p.nsm = np.array([1.0])

        #self.assertEqual(p.Nsm(), 1.0)
        # MassPerArea
        #self.assertAlmostEqual(p.MassPerArea(), 1.6)
        #self.assertAlmostEqual(p.MassPerArea(0, method='nplies'), 0.1+1/3.)
        #self.assertAlmostEqual(p.MassPerArea(1, method='nplies'), 0.2+1/3.)
        #self.assertAlmostEqual(p.MassPerArea(2, method='nplies'), 0.3+1/3.)

        #self.assertAlmostEqual(p.MassPerArea(0, method='rho*t'), 0.1+1/6.)
        #self.assertAlmostEqual(p.MassPerArea(1, method='rho*t'), 0.2+2/6.)
        #self.assertAlmostEqual(p.MassPerArea(2, method='rho*t'), 0.3+3/6.)

        #self.assertAlmostEqual(p.MassPerArea(0, method='t'), 0.1+1/6.)
        #self.assertAlmostEqual(p.MassPerArea(1, method='t'), 0.2+2/6.)
        #self.assertAlmostEqual(p.MassPerArea(2, method='t'), 0.3+3/6.)
        #with self.assertRaises(IndexError):
            #p.MassPerArea(3, method='nplies')

        z = p.get_z_locations()
        z_expected = array([0., T[0], T[0]+T[1], T[0]+T[1]+T[2]])
        assert np.allclose(z, z_expected), f'z={z} z_expected={z_expected}'

        #z0  =
        p.z0 = np.array([1.0])
        z_expected = 1.0 + z_expected
        z = p.get_z_locations()
        assert np.allclose(z, z_expected), f'z={z} z_expected={z_expected}'

    def test_pcomp_update(self):
        model = BDF(debug=False)
        pid = 1
        mids = [1, 2, 3]
        thicknesses = 1.
        thetas = 0.
        model.add_pcomp(
            pid, mids, thicknesses, thetas=None, souts=None,
            nsm=0., sb=0., ft=None, tref=0., ge=0.,
            lam=None, z0=None, comment='')
        model.setup(run_geom_check=False)

        pids = [1]
        thickness = np.array([
            [2., 3., 4.],
        ])
        theta = np.array([
            [20., 30., 40.],
        ])
        model.pcomp.update_layers(pids, thickness=thickness, theta=theta)
        print(model.pcomp)

    def test_pcomp_02(self):
        """symmetrical, nsm=0.0 and nsm=1.0"""
        pid = 1
        z0 = 0.
        nsm = 0.
        sb = 0.
        ft = 0.
        tref = 0.
        ge = 0.
        lam = 'SYM'  # is_symmetrical SYM
        Mid = [1, 2, 3]
        theta = [0., 10., 20.]
        T = [.1, .2, .3]
        sout = [1, 1, 0]  # 0-NO, 1-YES
        data = [pid, z0, nsm, sb, ft, tref, ge, lam, Mid, T, theta, sout]
        model = BDF(debug=True, log=None, mode='msc')
        pcomp = model.pcomp
        pcomp.add_op2_data(data)
        pcomp.parse_cards()
        p = pcomp
        self.assertTrue(p.is_symmetrical)
        self.assertEqual(p.nplies, 3)
        self.assertEqual(p.nplies_total, 6)

        #self.assertAlmostEqual(p.Thickness(), 1.2)
        #self.assertAlmostEqual(p.Thickness(0), 0.1)
        #self.assertAlmostEqual(p.Thickness(1), 0.2)
        #self.assertAlmostEqual(p.Thickness(2), 0.3)
        #self.assertAlmostEqual(p.Thickness(3), 0.3)
        #self.assertAlmostEqual(p.Thickness(4), 0.2)
        #self.assertAlmostEqual(p.Thickness(5), 0.1)
        #with self.assertRaises(IndexError):
            #p.Thickness(6)

        #self.assertAlmostEqual(p.Theta(0), 0.)
        #self.assertAlmostEqual(p.Theta(1), 10.)
        #self.assertAlmostEqual(p.Theta(2), 20.)
        #self.assertAlmostEqual(p.Theta(3), 20.)
        #self.assertAlmostEqual(p.Theta(4), 10.)
        #self.assertAlmostEqual(p.Theta(5), 0.)
        #with self.assertRaises(IndexError):
            #p.Theta(6)

        #self.assertEqual(p.Mid(0), 1)
        #self.assertEqual(p.Mid(1), 2)
        #self.assertEqual(p.Mid(2), 3)
        #self.assertEqual(p.Mid(3), 3)
        #self.assertEqual(p.Mid(4), 2)
        #self.assertEqual(p.Mid(5), 1)
        #with self.assertRaises(IndexError):
            #p.Mid(6)

        #self.assertEqual(p.Mids(), [1, 2, 3, 3, 2, 1])

        #self.assertEqual(p.sout(0), 'YES')
        #self.assertEqual(p.sout(1), 'YES')
        #self.assertEqual(p.sout(2), 'NO')
        #self.assertEqual(p.sout(3), 'NO')
        #self.assertEqual(p.sout(4), 'YES')
        #self.assertEqual(p.sout(5), 'YES')
        #with self.assertRaises(IndexError):
            #p.sout(6)


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
        with self.assertRaises(AssertionError):
            m = model.mat1.add_op2_data(mat1)

        G = 42.
        mcsid = 0
        mat1 = [mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid]
        m = model.mat1.add_op2_data(mat1)

        #for iply in range(len(p.plies)):
            #mid = p.plies[iply][0]
            #p.mids[iply] = m # MAT1
        #p.mids_ref = p.mids

        #Rho
        #self.assertAlmostEqual(p.Rho(0), 1.0)
        #self.assertAlmostEqual(p.Rho(1), 1.0)
        #self.assertAlmostEqual(p.Rho(2), 1.0)
        #self.assertAlmostEqual(p.Rho(3), 1.0)
        #self.assertAlmostEqual(p.Rho(4), 1.0)
        #self.assertAlmostEqual(p.Rho(5), 1.0)
        #with self.assertRaises(IndexError):
            #p.Rho(6)

        ## MassPerArea
        #self.assertAlmostEqual(p.MassPerArea(), 1.2)
        #self.assertAlmostEqual(p.MassPerArea(0), 0.1)
        #self.assertAlmostEqual(p.MassPerArea(1), 0.2)
        #self.assertAlmostEqual(p.MassPerArea(2), 0.3)
        #self.assertAlmostEqual(p.MassPerArea(3), 0.3)
        #self.assertAlmostEqual(p.MassPerArea(4), 0.2)
        #self.assertAlmostEqual(p.MassPerArea(5), 0.1)
        #with self.assertRaises(IndexError):
            #p.MassPerArea(6)

        #self.assertEqual(p.Nsm(), 0.0)
        #----------------------
        # change the nsm to 1.0
        p.nsm = 1.0

        #self.assertEqual(p.Nsm(), 1.0)
        # MassPerArea
        #self.assertAlmostEqual(p.MassPerArea(), 2.2)
        #self.assertAlmostEqual(p.MassPerArea(0, method='nplies'), 0.1+1/6.)
        #self.assertAlmostEqual(p.MassPerArea(1, method='nplies'), 0.2+1/6.)
        #self.assertAlmostEqual(p.MassPerArea(2, method='nplies'), 0.3+1/6.)
        #self.assertAlmostEqual(p.MassPerArea(3, method='nplies'), 0.3+1/6.)
        #self.assertAlmostEqual(p.MassPerArea(4, method='nplies'), 0.2+1/6.)
        #self.assertAlmostEqual(p.MassPerArea(5, method='nplies'), 0.1+1/6.)
        #with self.assertRaises(IndexError):
            #p.MassPerArea(6)

    def test_cshear(self):
        """tests a PSHEAR/CSHEAR"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        eid = 10
        pid = 20
        mid = 30
        t = 0.1
        nids = [1, 2, 3, 4]

        cshear_id = model.add_cshear(eid, pid, nids, comment='cshear')
        pshear_id = model.add_pshear(pid, mid, t, nsm=0., f1=0., f2=0., comment='')

        dvids = [1]
        coeffs = 1.0
        model.add_dvprel1(1, 'PSHEAR', pid, 'T', dvids, coeffs,
                          p_min=None, p_max=1e20,
                          c0=0.0, validate=True,
                          comment='')
        model.add_desvar(1, 'T', 10.0)

        E = 30.e7
        G = None
        nu = 0.3
        mat1_id = model.add_mat1(mid, E, G, nu, rho=0.1, comment='mat1')
        model.setup()
        model.pop_parse_errors()
        model.validate()

        #cshear.raw_fields()
        model.cshear.write(size=8)

        #pshear.raw_fields()
        model.pshear.write(size=8)
        model.pshear.write(size=16)
        model.pshear.write(size=16, is_double=True)

        model.validate()
        model._verify_bdf(xref=False)
        #model.cross_reference()
        #model._verify_bdf(xref=True)
        #mass_properties(model)
        model.mass()

        model.cshear.write(size=8)
        model.pshear.write(size=8)
        if 0:
            model.update_model_by_desvars()
        save_load_deck(model)

    def test_shells(self):
        """tests a CTRIA3/CQUAD4/PSHELL and CTRIA6/CQUAD8/CQUAD/PCOMP"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        model.add_grid(5, [.5, 0., 0.])
        model.add_grid(6, [1., 0.5, 0.])
        model.add_grid(7, [.5, 1., 0.])
        model.add_grid(8, [0., .5, 0.])

        model.add_grid(9, [.5, .5, 0.])

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
        ctria6_id = model.add_ctria6(3, pid, nids, comment='ctria6')

        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        cquad8_id = model.add_cquad8(4, pid, nids, comment='cquad8')

        #nids = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        #cquad = model.add_cquad(5, pid, nids, comment='cquad')

        mids = [1, 2, 3]
        thicknesses = [0.1, 0.2, 0.3]
        model.add_pcomp(pid, mids, thicknesses)
        model.setup()

        pcomp = model.pcomp.slice_card_by_property_id(pid)
        assert np.allclose(pcomp.total_thickness(), sum(thicknesses)), f'total={pcomp.total_thickness()}; thicknesses={thicknesses} sum={sum(thicknesses)}'
        assert np.allclose(pcomp.thickness, [0.1, 0.2, 0.3]), pcomp.thickness
        assert np.allclose(pcomp.theta, [0., 0., 0.]), pcomp.theta

        pcomp.lam[0] = 'SYM'
        assert np.allclose(pcomp.total_thickness(), sum(thicknesses)*2), thicknesses

        #assert np.allclose(pcomp.get_thicknesses(), [0.1, 0.2, 0.3, 0.3, 0.2, 0.1]), pcomp.get_thicknesses()
        #assert np.allclose(pcomp.get_thetas(), [0., 0., 0., 0., 0., 0.]), pcomp.get_thetas()
        #---------------------------------------------------
        model.validate()

        #ctria6.raw_fields()
        model.ctria6.write(size=8)

        #cquad8.raw_fields()
        model.cquad8.write(size=8)

        #model.cquad.raw_fields()
        model.cquad.write(size=8)

        #model.pcomp.raw_fields()
        model.pcomp.write(size=8)
        model.pcomp.write(size=16)
        model.pcomp.write(size=16, is_double=True)

        model._verify_bdf(xref=False)

        params = [('T1', 1.0), ('THETA1', 2.0), ('Z0', 3.0), ('SB', 4.0),
                  ('TREF', 0.0), ('GE', 0.1)]
        make_dvprel_optimization(model, params, 'PCOMP', pid)
        #--------------------------------
        model.cross_reference()
        model._verify_bdf(xref=True)

        model.get_area_breakdown(property_ids=None, stop_if_no_area=True)
        model.get_mass_breakdown(stop_if_no_mass=True)
        #model.get_mass_breakdown_detailed(property_ids=None, stop_if_no_mass=True)
        model.get_volume_breakdown(property_ids=None, stop_if_no_volume=True)
        if 0:
            model.update_model_by_desvars(xref=True)

        #ctria6.raw_fields()
        model.ctria6.write(size=8)

        #cquad8.raw_fields()
        model.cquad8.write(size=8)

        #cquad.raw_fields()
        model.cquad.write(size=8)

        #pcomp.raw_fields()
        model.pcomp.write(size=8)
        model.pcomp.write(size=16)
        model.pcomp.write(size=16, is_double=True)
        model._verify_bdf(xref=False)

    def test_ctriar_cquadr(self):
        """tests a CTRIAR/PSHELL/MAT8"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        eid = 6
        pid = 13
        nids = [1, 2, 3]
        ctriar_id = model.add_ctriar(eid, pid, nids, comment='ctriar')
        eid += 1
        model.add_ctriar(eid, pid, nids, theta_mcid=0, comment='ctriar')
        eid += 1
        model.add_ctriar(eid, pid, nids, theta_mcid=90., comment='ctriar')

        eid = 10
        nids = [1, 2, 3, 4]
        cquadr_id = model.add_cquadr(eid, pid, nids, comment='cquadr')
        eid += 1
        cquadr_id = model.add_cquadr(eid, pid, nids, theta_mcid=0, comment='cquadr')
        eid += 1
        cquadr_id = model.add_cquadr(eid, pid, nids, theta_mcid=45., comment='cquadr')
        #cquadr.raw_fields()

        mid = 42
        model.add_pshell(pid, mid1=mid, t=0.2)
        e11 = 1e7
        e22 = 1e6
        nu12 = 0.3
        model.add_mat8(mid, e11, e22, nu12, rho=10.)
        model.setup()

        #ctriar.raw_fields()
        model.ctriar.write(size=8)
        model.ctriar.write(size=8, is_double=False)
        model.ctriar.write(size=16, is_double=False)
        model.ctriar.flip_normal()

        model.cquadr.write(size=8)
        #cquadr.write_card(size=16, is_double=False)
        model.cquadr.flip_normal()

        model.ctriar.element_coordinate_system()
        model.ctriar.material_coordinate_system()
        model.cquadr.element_coordinate_system()
        model.cquadr.material_coordinate_system()

        model.validate()
        model._verify_bdf(xref=False)
        #model.cross_reference()
        #model._verify_bdf(xref=True)
        #model.uncross_reference()
        #model.safe_cross_reference()

        mass1 = model.ctriar.mass_breakdown()
        mass2 = model.cquadr.mass_breakdown()

        model.get_area_breakdown(property_ids=None, stop_if_no_area=True)
        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True)
        #model.get_mass_breakdown_detailed(property_ids=None, stop_if_no_mass=True)
        model.get_volume_breakdown(property_ids=None, stop_if_no_volume=True)

        save_load_deck(model)

    def test_cplstn34(self):
        """tests a CPLSTN3, CPLSTN4/PSHELL/MAT8"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        pid = 4
        eid = 3
        nids = [1, 2, 3, 4]
        model.add_cplstn4(eid, pid, nids, comment='cplstn4')
        cplstn4 = model.cplstn4
        cplstn4.flip_normal()

        eid = 5
        nids = [1, 2, 3]
        mid = 10
        model.add_cplstn3(eid, pid, nids, comment='cplstn3')
        cplstn3 = model.cplstn3
        cplstn3.flip_normal()

        pplane = model.add_pplane(
            pid, mid, t=0.1, nsm=0.,
            formulation_option=0, comment='pplane')
        E = 1e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.setup()
        #cplstn3.repr_fields()
        #cplstn4.repr_fields()

        #cplstn3.raw_fields()
        #cplstn4.raw_fields()
        #pplane.raw_fields()

        model.validate()
        model._verify_bdf(xref=False)
        model.cplstn3.write(size=8)
        model.cplstn4.write(size=8)
        model.pplane.write(size=8)
        model.cross_reference()
        model.cplstn3.write(size=16)
        model.cplstn4.write(size=16)

        #model.uncross_reference()
        #model.safe_cross_reference()
        save_load_deck(model)

    def test_cplstn68(self):
        """tests a CPLSTN6, CPLSTN8/PSHELL/MAT8"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(5, [.5, 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(6, [1., .5, 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(7, [.5, 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(8, [0., .5, 0.])
        pid = 4
        eid = 3
        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        model.add_cplstn8(eid, pid, nids, comment='cplstn8')
        #cplstn8 = model.cplstn8

        eid = 5
        nids = [1, 2, 3, 4, 5, 6]
        mid = 10
        model.add_cplstn6(eid, pid, nids, comment='cplstn6')
        #cplstn6 = model.cplstn6

        model.add_pplane(
            pid, mid, t=0.1, nsm=0.,
            formulation_option=0, comment='pplane')
        #pplane = model.pplane
        E = 1e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        mat1 = model.mat1
        model.setup()

        #cplstn6.raw_fields()
        #cplstn8.raw_fields()
        #pplane.raw_fields()

        model.validate()
        model._verify_bdf(xref=False)
        model.cplstn6.write(size=8)
        model.cplstn8.write(size=8)
        model.pplane.write(size=8)
        model.cross_reference()
        model.cplstn3.write(size=8)
        model.cplstn4.write(size=8)

        #model.uncross_reference()
        #model.safe_cross_reference()
        save_load_deck(model)

    def test_cplstn34(self):
        """tests a CPLSTN3, CPLSTN4/PSHELL/MAT8"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        pid = 4
        eid = 3
        nids = [1, 2, 3, 4]
        model.add_cplstn4(eid, pid, nids, comment='cplstn4')
        cplstn4 = model.cplstn4
        cplstn4.flip_normal()

        eid = 5
        nids = [1, 2, 3]
        mid = 10
        model.add_cplstn3(eid, pid, nids, comment='cplstn3')
        cplstn3 = model.cplstn3
        cplstn3.flip_normal()

        pplane = model.add_pplane(
            pid, mid, t=0.1, nsm=0.,
            formulation_option=0, comment='pplane')
        E = 1e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.setup()
        #cplstn3.repr_fields()
        #cplstn4.repr_fields()

        #cplstn3.raw_fields()
        #cplstn4.raw_fields()
        #pplane.raw_fields()

        model.validate()
        model._verify_bdf(xref=False)
        model.cplsts3.write(size=8)
        model.cplsts4.write(size=8)
        model.pplane.write(size=8)
        model.cross_reference()
        model.cplsts3.write(size=16)
        model.cplsts4.write(size=16)

        #model.uncross_reference()
        #model.safe_cross_reference()
        save_load_deck(model)

    def test_cplsts68(self):
        """tests a CPLSTN6, CPLSTN8/PSHELL/MAT8"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(5, [.5, 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(6, [1., .5, 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(7, [.5, 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(8, [0., .5, 0.])
        pid = 4
        eid = 3
        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        model.add_cplsts8(eid, pid, nids, comment='cplstn8')

        eid = 5
        nids = [1, 2, 3, 4, 5, 6]
        mid = 10
        model.add_cplsts6(eid, pid, nids, comment='cplstn6')

        model.add_pplane(
            pid, mid, t=0.1, nsm=0.,
            formulation_option=0, comment='pplane')
        #pplane = model.pplane
        E = 1e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        mat1 = model.mat1
        model.setup()

        #cplstn6.raw_fields()
        #cplstn8.raw_fields()
        #pplane.raw_fields()

        model.validate()
        model._verify_bdf(xref=False)
        model.cplsts6.write(size=8)
        model.cplsts8.write(size=8)
        model.pplane.write(size=8)
        model.cross_reference()
        model.cplstn3.write(size=8)
        model.cplstn4.write(size=8)

        #model.uncross_reference()
        #model.safe_cross_reference()
        save_load_deck(model)

    def test_ctrishell68(self):
        """tests a CPLSTN6, CPLSTN8/PSHELL/MAT8"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(5, [.5, 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(6, [1., .5, 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(7, [.5, 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(8, [0., .5, 0.])
        pid = 4
        eid = 3
        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        cquad8_id = model.add_cquad8(eid, pid, nids, comment='cquad8')

        eid = 5
        nids = [1, 2, 3, 4, 5, 6]
        mid = 10
        ctria6_id = model.add_ctria6(eid, pid, nids, comment='ctria6')

        plplane_id = model.add_plplane(pid, mid, cid=0,
                                       stress_strain_output_location='GRID', comment='')

        pid = 5
        pplane_id = model.add_pplane(pid, mid, t=0.1, nsm=0.,
                                     formulation_option=0, comment='pplane')
        E = 1e7
        G = None
        nu = 0.3
        mat1_id = model.add_mat1(mid, E, G, nu, rho=1.)

        #ctria6.raw_fields()
        #cquad8.raw_fields()
        #pplane.raw_fields()
        model.setup()

        model.validate()
        model._verify_bdf(xref=False)
        model.ctria6.write(size=8)
        model.cquad8.write(size=8)
        model.pplane.write(size=8)
        model.plplane.write(size=8)

        mass1 = model.ctria6.mass_breakdown()
        mass2 = model.cquad8.mass_breakdown()

        #model.cross_reference()

        #model.uncross_reference()
        #model.safe_cross_reference()
        save_load_deck(model, run_mass_properties=False, run_test_bdf=False)
        #mass_properties(model)

    def test_shear(self):
        """tests a CSHEAR, PSHEAR"""
        pid = 10
        pid_pshell = 11

        mid = 100
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        nsm = 10.0
        t = 1.0
        rho = 1.0
        cshear_id = model.add_cshear(10, pid, [1, 2, 3, 4],
                                     comment='cshear')

        cquad4_id = model.add_cquad4(14, pid_pshell, [1, 2, 3, 4],
                                     comment='cquad4')
        model.add_pshear(pid, mid, t=t,
                         nsm=nsm, f1=0., f2=0., comment='pshear')
        model.add_pshell(pid_pshell, mid1=mid, t=t, mid2=None, twelveIt3=1.0,
                         mid3=None, tst=0.833333,
                         nsm=nsm, z1=None, z2=None,
                         mid4=None, comment='')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=rho)
        model.validate()
        model.setup()

        xyz_cid0_expected = np.array([
            [0., 0., 0.], # 1
            [1., 0., 0.], # 2
            [1., 1., 0.], # 3
            [0., 1., 0.], # 4
        ])
        xyz_cid0 = model.grid.xyz_cid0()
        #model.cross_reference()
        assert np.array_equal(model.grid.node_id, [1, 2, 3, 4])
        assert np.allclose(xyz_cid0, xyz_cid0_expected)

        assert model.cshear.n > 0, model.cshear

        area = 1.0
        mass_expected = area * (rho * t + nsm)
        #mass = mass_properties(model, )[0]
        model.inertia()
        model.mass()
        mass = model.mass_sum()
        assert np.allclose(mass, mass_expected*2), 'mass_properties all: mass=%s mass_expected=%s' % (mass, mass_expected*2)

        #mass = mass_properties(model, element_ids=10)[0]
        mass = model.mass_sum(element_id=10)
        assert np.allclose(mass, mass_expected), 'mass_properties reduced: mass=%s mass_expected=%s' % (mass, mass_expected)

        if 0:
            mass = mass_properties_nsm(model)[0]
            assert np.allclose(mass, mass_expected*2), 'mass_properties_nsm all: mass=%s mass_expected=%s' % (mass, mass_expected*2)

            mass = mass_properties_nsm(model, element_ids=10)[0]
            assert np.allclose(mass, mass_expected), 'mass_properties_nsm reduced: mass=%s mass_expected=%s' % (mass, mass_expected)

        bdf_file = StringIO()
        model.write_bdf(bdf_file)
        #model.uncross_reference()
        #model.cross_reference()

        model.get_area_breakdown(property_ids=None, stop_if_no_area=True)
        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True)
        #model.get_mass_breakdown_detailed(property_ids=None, stop_if_no_mass=True)
        model.get_volume_breakdown(property_ids=None, stop_if_no_volume=True)

        cshear = model.cshear
        assert np.allclose(cshear.mass(), mass_expected), cshear.Mass()

        #model.uncross_reference()
        #model.safe_cross_reference()
        #model.uncross_reference()

        #bdf_file = model.write_bdf(bdf_file)

        save_load_deck(model)

    def test_ctria6_cquad8_cquad9(self):
        """tests a CQUAD8 and CQUAD9"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 1
        pid = 10
        mid = 100
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(5, [.5, 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(6, [1., .5, 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(7, [.5, 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(8, [0., .5, 0.])
        model.add_grid(9, [.5, .5, 0.])

        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        cquad8_id = model.add_cquad8(eid, pid, nids, theta_mcid=0., comment='cquad8')
        model.cquad8.flip_normal()

        eid = 2
        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        cquad_id = model.add_cquad(eid, pid, nids, theta_mcid=0., comment='cquad')
        model.add_pshell(pid, mid1=mid, t=1.0)

        eid = 3
        nids = [1, 2, 3, 5, 6, 9]
        ctria6_id = model.add_ctria6(eid, pid, nids, theta_mcid=0., comment='ctria6')
        model.ctria6.flip_normal()

        eid = 4
        cquad4_id = model.add_cquad4(eid, pid, [1, 2, 3, 4])
        model.cquad4.flip_normal()
        str(model.cquad4)

        eid = 5
        # relative
        cquad4 = model.add_cquad4(eid, pid, [1, 2, 3, 4],
                                  tflag=1, T1=2., T2=2., T3=2., T4=2.)
        str(model.cquad4)

        eid = 6
        ctria3_id = model.add_ctria3(eid, pid, [1, 2, 3])
        model.ctria3.flip_normal()
        str(model.ctria3)

        eid = 7
        # relative
        ctria3_id = model.add_ctria3(eid, pid, [1, 2, 3],
                                     tflag=1, T1=2., T2=2., T3=2.)
        str(model.ctria3)

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1)

        #model.cross_reference()

        model.grid.sort()
        model.setup()
        model.ctria3.flip_normal()
        model.cquad4.flip_normal()
        model.ctria6.flip_normal()
        model.cquad8.flip_normal()

        assert len(model.ctria6.centroid()[0]) == 3, model.ctria6.centroid()
        assert len(model.ctria6.center_of_mass()[0]) == 3, model.ctria6.center_of_mass()

        assert np.allclose(model.cquad8.mass(), 0.1), model.cquad8.mass()
        assert np.allclose(model.cquad.mass(), 0.1), model.cquad.mass()
        assert np.allclose(model.ctria6.mass(), 0.05), model.ctria6.mass()

        model.get_area_breakdown(property_ids=None, stop_if_no_area=True)
        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True)
        #model.get_mass_breakdown_detailed(property_ids=None, stop_if_no_mass=True)
        model.get_volume_breakdown(property_ids=None, stop_if_no_volume=True)

        save_load_deck(model)

    def test_pplane_pshln1(self):
        """tests the PPLANE, PSHELN1"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        eid = 10
        pid = 100
        mid = 1000

        # only when t=0.0
        model.add_pplane(pid, mid, t=0.0, nsm=0.0, formulation_option=0, comment='')
        pplane = model.pplane
        pshln1 = model.pshln1
        #with self.assertRaises(AssertionError):
            #t = pplane.total_thickness()

        model.add_pshln1(
            pid, mid1=0, mid2=0, analysis='ISH',
            behx=None, integration=None,
            behxh=None, integration_h=None, comment='')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1)

        card_lines = [
            'PSHLN1,1,,,NO',
            '+       C3      DCTN    LDK                                             +',
            '+       C4      DCT     L                                               +',
            '+       C8      MB      Q'
        ]
        model.add_card_lines(card_lines, 'PSHLN1', comment='', has_none=True)
        model.setup()
        #t = pplane.total_thickness()
        #assert np.allclose(t, 1.0), t
        save_load_deck(model, run_read_write=False)

    def test_pplane_pshln2(self):
        """tests the PPLANE, PSHELN2"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        eid = 10
        pid = 100
        mid = 1000

        # only when t=0.0
        model.add_pplane(pid, mid, t=0.0, nsm=0.0, formulation_option=0, comment='')
        pplane = model.pplane
        pshln2 = model.pshln2
        #with self.assertRaises(AssertionError):
            #t = pplane.total_thickness()

        model.add_pshln2(
            pid, mid, direct=1, thickness=1.0, analysis='ISH',
            behx=None, integration=None, behxh=None, integration_h=None,
            comment='pshln2')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1)

        model.setup()
        t = pplane.total_thickness()
        assert np.allclose(t, 1.0), t
        save_load_deck(model)

    def test_shell_mcid(self):
        """tests that mcids=0 are correctly identified as not 0.0 and thus not dropped"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        eid = 10
        pid = 100
        mid = 1000
        model.add_ctria3(eid, pid, [1, 2, 3], zoffset=0., theta_mcid=0, tflag=0,
                         T1=None, T2=None, T3=None,
                         comment='')

        eid = 11
        model.add_cquad4(eid, pid, [1, 2,3, 4], theta_mcid=0, zoffset=None, tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')

        model.add_pshell(pid, mid1=mid, t=0.1, mid2=mid, twelveIt3=1.0,
                         mid3=None, tst=0.833333,
                         nsm=0.0, z1=None, z2=None,
                         mid4=None, comment='')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1)
        model.setup()
        #print(model.elements[11])
        assert model.ctria3.slice_card_by_element_id(10).rstrip()         == 'CTRIA3        10     100       1       2       3       0'
        assert model.ctria3.slice_card_by_element_id(10).write().rstrip() == 'CTRIA3        10     100       1       2       3       0'
        assert model.cquad4.slice_card_by_element_id(11).rstrip()         == 'CQUAD4        11     100       1       2       3       4       0'

        #model.cross_reference()
        unused_ABD = model.pshell.get_ABD_matrices(theta_offset=0.)
        assert model.ctria3.slice_card_by_element_id(10).rstrip() == 'CTRIA3        10     100       1       2       3       0'
        assert model.cquad4.slice_card_by_element_id(11).rstrip() == 'CQUAD4        11     100       1       2       3       4       0'

        #model.uncross_reference()
        assert model.ctria3.slice_card_by_element_id(10).rstrip() == 'CTRIA3        10     100       1       2       3       0'
        assert model.cquad4.slice_card_by_element_id(11).rstrip() == 'CQUAD4        11     100       1       2       3       4       0'
        #model.safe_cross_reference()
        #model.uncross_reference()
        assert model.ctria3.slice_card_by_element_id(10).rstrip() == 'CTRIA3        10     100       1       2       3       0'
        assert model.cquad4.slice_card_by_element_id(11).rstrip() == 'CQUAD4        11     100       1       2       3       4       0'

        model2 = save_load_deck(model)
        #model2.Elements(10)[0].comment = ''
        assert model2.ctria3.slice_card_by_element_id(10).rstrip() == 'CTRIA3        10     100       1       2       3       0'
        assert model2.cquad4.slice_card_by_element_id(11).rstrip() == 'CQUAD4        11     100       1       2       3       4       0'

    def test_abd_1(self):
        """tests some ABD matrix functionality for a PCOMP"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        nids = [1, 2, 3, 4]
        eid = 1
        pid = 10
        mid = 20
        model.add_cquad4(eid, pid, nids, theta_mcid=0.0, zoffset=0.,
                         tflag=0, T1=None, T2=None, T3=None, T4=None, comment='')

        thetas = [0., 45., 90.]
        thicknesses = [0.1] * 3
        mids = len(thicknesses) * [mid]
        pcomp_id = model.add_pcomp(pid, mids, thicknesses, thetas=None,
                                   souts=None, nsm=0., sb=0., ft=None, tref=0., ge=0.,
                                   lam=None, z0=0., comment='')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        #--------------------------
        #e1_e2 = 40.
        #g12_e2 = 0.5
        #nu12 = 0.25
        #e22 = 30e6

        e1_e2 = 3.
        g12_e2 = 0.5
        nu12 = 0.25
        e22 = 1.

        e11 = e1_e2 * e22
        g12 = g12_e2 * e22

        mid8 = 8
        mat8 = model.add_mat8(
            mid8, e11, e22, nu12, g12=g12, g1z=1e8, g2z=1e8, rho=0.,
            a1=0., a2=0., tref=0.,
            Xt=0., Xc=None, Yt=0., Yc=None, S=0., ge=0.,
            F12=0., strn=0., comment='')
        #S = get_mat_props_S(mat8)

        pid8 = 8
        pcomp8 = model.add_pcomp(pid8, [mid8], [1.], thetas=[0.],
                                 souts=None, nsm=0., sb=0., ft=None, tref=0., ge=0.,
                                 lam=None, z0=0., comment='')

        model.pop_parse_errors()
        model.cross_reference()
        ABD = model.pcomp.get_ABD_matrices()

        thetad = np.linspace(0., 90., num=91)
        if IS_MATPLOTLIB:
            plot_equivalent_lamina_vs_theta(
                pcomp8, mat8, thetad, plot=True, show=False, close=True,
                png_filename_base='lamina')
            os.remove('lamina_stiffness.png')
            os.remove('lamina_nu.png')

    def test_abd_2(self):
        """tests some ABD matrix functionality for a PCOMP"""
        log = get_logger(level='warning')
        model = BDF(log=log)

        #--------------------------
        mid = 2
        #MAT8*                  2        179000.0          8110.0           0.317*
        #*                 4140.0
        e11 = 179000.
        e22 = 8110.
        nu12 = 0.317
        g12 = 4140.
        model.add_mat8(mid, e11, e22, nu12, g12=g12, g1z=1e8, g2z=1e8,
                       rho=0., a1=0., a2=0., tref=0., Xt=0., Xc=None, Yt=0., Yc=None,
                       S=0., ge=0., F12=0., strn=0., comment='')

        #--------------------------
        #PCOMP*              1601                                                *
        #*                                                                    SYM
        #*                      2           0.185            45.0             YES
        #*                      2           0.185            90.0             YES
        #*                      2           0.185           -45.0             YES
        #*                      2           0.185             0.0             YES
        #*                      2           0.185           -45.0             YES
        #*                      2           0.185            90.0             YES
        #*                      2           0.185            45.0             YES
        #*                      2           0.185             0.0             YES
        thetas = [45., 90., -45., 0., -45., 90., 45.0, 0.]
        thicknesses = [0.185] * len(thetas)
        mids = [2] * len(thetas)
        pid = 2
        pcomp_sym = model.add_pcomp(
            pid, mids, thicknesses, thetas=thetas,
            souts=None, nsm=0., sb=0., ft=None, tref=0., ge=0.,
            lam='SYM', z0=None, comment='')

        model.setup()
        pcomp_sym = model.pcomp.slice_card_by_property_id(pid)

        thetas = [45., 90., -45., 0., -45., 90., 45.0, 0.]
        thetas += thetas[::-1]
        #print(f'*****thetas = {thetas}')
        thicknesses = [0.185] * len(thetas)
        mids = [2] * len(thetas)
        pid = 3
        pcomp_asym = model.add_pcomp(
            pid, mids, thicknesses, thetas=thetas,
            souts=None, nsm=0., sb=0., ft=None, tref=0., ge=0.,
            lam=None, z0=None, comment='')
        model.setup()

        #model.pop_parse_errors()
        #model.cross_reference()
        ABD = model.pcomp.get_ABD_matrices()

        #print(f'pcomp_sym:\n{pcomp_sym}')
        ABD2 = pcomp_sym.get_ABD_matrices()
        if 0:
            is_balanced, is_symmetric = pcomp_sym.is_balanced_symmetric(debug=False)
            assert is_symmetric
            A = ABD2[:3, :3]
            B = ABD2[-3:, :3]
            D = ABD2[-3:, -3:]
            #print(ABD2)
            #print(f'A:\n{A}\n')
            #print(f'B:\n{B}\n')
            #print(f'D:\n{D}\n')

            #print('-------------')
            #print(f'pcomp_asym:\n{pcomp_asym}')
            ABD3 = pcomp_asym.get_ABD_matrices()
            is_balanced, is_symmetric = pcomp_asym.is_balanced_symmetric(debug=False)
            assert is_symmetric
            A = ABD3[:3, :3]
            B = ABD3[-3:, :3]
            D = ABD3[-3:, -3:]
            #print(ABD3)
            #print(f'A:\n{A}\n')
            #print(f'B:\n{B}\n')
            #print(f'D:\n{D}\n')
            B = ABD2[:3, -3:]
            assert np.allclose(0., B.sum()), B
            assert np.allclose(ABD2, ABD3), ABD2

    #def test_lax(self):
        #"""tests adding cards in a lax way"""
        #card = ['CTRIA3', '301471', '301003', '301466', '301468', '30071', '0', '0']
        #bdf_card = BDFCard(card)
        #str(bdf_card)
        #card = CTRIA3.add_card_lax(bdf_card)
        #str(card)

    def test_tri_theta(self):
        model = BDF(debug=False, log=None, mode='msc')
        rho = 0.
        G = None
        E = 3.0e7
        t = 0.1
        nu = 0.3
        nsm = 0.
        self._make_ctria3(model, rho, nu, G, E, t, nsm)

        eid = 10
        elem = model.Element(eid)[0]

        #['grid', n1, 0, 0., 0., 0.],
        #['grid', n2, 0, 4., 0., 0.],
        #['grid', n3, 0, 4., 1., 0.],

        # element coordinate system and material coordinate system are coincident
        dxyz, centroid, imat, jmat, normal = elem.material_coordinate_system()
        centroid_expected = [2.66666667, 0.33333333, 0.]
        imat_expected = [1., 0., 0.]
        jmat_expected = [0., 1., 0.]
        normal_expected = [0., 0., 1.]

        dxyz, centroid, imat, jmat, normal = elem.element_coordinate_system()
        assert np.allclose(centroid, centroid_expected), centroid
        assert np.allclose(imat, imat_expected), imat
        assert np.allclose(jmat, jmat_expected), jmat
        assert np.allclose(normal, normal_expected), normal

        assert np.allclose(centroid, centroid_expected), centroid
        assert np.allclose(imat, imat_expected), imat
        assert np.allclose(jmat, jmat_expected), jmat
        assert np.allclose(normal, normal_expected), normal

        # ----------------------------------------------------------------
        elem.theta[0] = 45.
        # theta doesn't change the element coordinate system
        dxyz, centroid, imat, jmat, normal = elem.element_coordinate_system()
        assert np.allclose(centroid, centroid_expected), centroid
        assert np.allclose(imat, imat_expected), imat
        assert np.allclose(jmat, jmat_expected), jmat
        assert np.allclose(normal, normal_expected), normal

        dxyz, centroid, imat, jmat, normal = elem.material_coordinate_system()

        sqrt2_2 = 2 ** 0.5 / 2
        imat_expected = [sqrt2_2, sqrt2_2, 0.]
        jmat_expected = [-sqrt2_2, sqrt2_2, 0.]
        normal_expected = [0., 0., 1.]
        assert np.allclose(centroid, centroid_expected), centroid
        assert np.allclose(imat, imat_expected), imat
        assert np.allclose(jmat, jmat_expected), jmat
        assert np.allclose(normal, normal_expected), normal

        # ----------------------------------------------------------------
        elem.theta[0] = 90.
        dxyz, centroid, imat, jmat, normal = elem.material_coordinate_system()

        imat_expected = [0., 1., 0.]
        jmat_expected = [-1., 0., 0.]
        normal_expected = [0., 0., 1.]
        assert np.allclose(centroid, centroid_expected), centroid
        assert np.allclose(imat, imat_expected), imat
        assert np.allclose(jmat, jmat_expected), jmat
        assert np.allclose(normal, normal_expected), normal

    def test_quad_theta(self):
        model = BDF(debug=False, log=None, mode='msc')
        rho = 0.
        G = None
        E = 3.0e7
        t = 0.1
        nu = 0.3
        nsm = 0.
        self._make_cquad4(model, rho, nu, G, E, t, nsm)

        eid = 10
        elem: CQUAD4 = model.Element(eid)[0]

        #['grid', n1, 0, 0., 0., 0.],
        #['grid', n2, 0, 2., 0., 0.],
        #['grid', n3, 0, 2., 1., 0.],
        #['grid', n4, 0, 0., 1., 0.],

        # element coordinate system and material coordinate system are coincident
        dxyz, centroid, imat, jmat, normal = elem.material_coordinate_system()
        centroid_expected = [1., 0.5, 0.]
        imat_expected = [1., 0., 0.]
        jmat_expected = [0., 1., 0.]
        normal_expected = [0., 0., 1.]

        dxyz, centroid, imat, jmat, normal = elem.element_coordinate_system()
        assert np.allclose(centroid, centroid_expected), centroid
        assert np.allclose(imat, imat_expected), imat
        assert np.allclose(jmat, jmat_expected), jmat
        assert np.allclose(normal, normal_expected), normal

        assert np.allclose(centroid, centroid_expected), centroid
        assert np.allclose(imat, imat_expected), imat
        assert np.allclose(jmat, jmat_expected), jmat
        assert np.allclose(normal, normal_expected), normal

        # ----------------------------------------------------------------
        elem.theta[0] = 45.
        # theta doesn't change the element coordinate system
        dxyz, centroid, imat, jmat, normal = elem.element_coordinate_system()
        assert np.allclose(centroid, centroid_expected), centroid
        assert np.allclose(imat, imat_expected), imat
        assert np.allclose(jmat, jmat_expected), jmat
        assert np.allclose(normal, normal_expected), normal

        # theta changed the material coordinate system
        dxyz, centroid, imat, jmat, normal = elem.material_coordinate_system()
        sqrt2_2 = 2 ** 0.5 / 2
        centroid_expected = [1., 0.5, 0.]
        imat_expected = [sqrt2_2, sqrt2_2, 0.]
        jmat_expected = [-sqrt2_2, sqrt2_2, 0.]
        normal_expected = [0., 0., 1.]
        assert np.allclose(centroid, centroid_expected), centroid
        assert np.allclose(imat, imat_expected), imat
        assert np.allclose(jmat, jmat_expected), jmat
        assert np.allclose(normal, normal_expected), normal
        # ----------------------------------------------------------------
        elem.theta[0] = 90.
        dxyz, centroid, imat, jmat, normal = elem.material_coordinate_system()

        imat_expected = [0., 1., 0.]
        jmat_expected = [-1., 0., 0.]
        normal_expected = [0., 0., 1.]
        assert np.allclose(centroid, centroid_expected), centroid
        assert np.allclose(imat, imat_expected), imat
        assert np.allclose(jmat, jmat_expected), jmat
        assert np.allclose(normal, normal_expected), normal

    def test_pshln2(self):
        model = BDF(debug=False)
        pid = 1
        mid = 2
        model.add_pshln2(pid, mid, direct=1, thickness=1.0, analysis='ISH',
                         behx=None, integration=None, behxh=None, integration_h=None, comment='')
        model.setup()
        model.pshln2.write()
        save_load_deck(model)

    def test_cplsts3(self):
        model = BDF(debug=False)
        cplsts3 = model.cplsts3
        cplsts4 = model.cplsts4
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        expected_tri_area = 0.5
        expected_quad_area = 1.0

        eid = 10
        pid = 11
        mid = 12
        nids = [1, 2, 3]
        model.add_cplsts3(eid, pid, nids, theta=0.0, comment='')

        eid += 1
        nids = [1, 2, 3, 4]
        model.add_cplsts4(eid, pid, nids, theta=0.0, comment='')

        t = 0.31
        rho = 0.13
        model.add_pplane(pid, mid, t=t, nsm=0.0, formulation_option=0, comment='')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=rho)
        model.setup()
        assert len(cplsts3) == 1, cplsts3
        assert len(cplsts4) == 1, cplsts4

        expected_tri_mass = t * expected_tri_area * rho
        expected_quad_mass = t * expected_quad_area * rho

        areai = cplsts3.area()
        massi = cplsts3.mass()
        assert np.allclose(areai, expected_tri_area), (areai, expected_tri_area)
        assert np.allclose(massi, expected_tri_mass), (massi, expected_tri_mass)

        areai = cplsts4.area()
        massi = cplsts4.mass()
        assert np.allclose(areai, expected_quad_area), (areai, expected_quad_area)
        assert np.allclose(massi, expected_quad_mass), (massi, expected_quad_mass)

        save_load_deck(model, run_remove_unused=False, run_mass_properties=False)
        x = 1

    def test_tri_volume(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        rho = 0.0
        nu = 0.3
        E = 3.0e7
        G = None
        t = 1.0
        nsm = 0.0
        self._make_ctria3(model, rho, nu, G, E, t, nsm)
        model.setup()

        ctria3 = model.ctria3
        nodes = ctria3.nodes
        nelements = len(ctria3)
        area = ctria3.area()
        average_thickness1 = 3.1
        average_thickness2 = np.array([1., 2., 3.])
        expected_volume1 = area * average_thickness1
        expected_volume2 = area * 3 # average_thickness2[1]

        dthickness1 = np.ones((nelements, 3), dtype='float64') * average_thickness1
        vol = tri_volume(model.grid, nodes, dthickness1)
        assert np.allclose(vol, expected_volume1), (vol, expected_volume1)

        dthickness2 = np.ones((nelements, 3), dtype='float64') * average_thickness2
        vol = tri_volume(model.grid, nodes, dthickness2)
        assert np.allclose(vol, expected_volume2), (vol, expected_volume2)

    def test_quad_volume(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        rho = 0.0
        nu = 0.3
        E = 3.0e7
        G = None
        t = 1.0
        nsm = 0.0
        self._make_cquad4(model, rho, nu, G, E, t, nsm)
        model.setup()

        cquad4 = model.cquad4
        nodes = cquad4.nodes
        nelements = len(cquad4)
        area = cquad4.area()
        average_thickness1 = 3.1
        average_thickness2 = np.array([1., 2., 4., 5.])
        expected_volume1 = area * average_thickness1
        expected_volume2 = area * 17 / 6 # average_thickness2[1]

        dthickness1 = np.ones((nelements, 4), dtype='float64') * average_thickness1
        vol = quad_volume(model.grid, nodes, dthickness1)
        assert np.allclose(vol, expected_volume1), (vol, expected_volume1)

        dthickness2 = np.ones((nelements, 4), dtype='float64') * average_thickness2
        vol = quad_volume(model.grid, nodes, dthickness2)
        assert np.allclose(vol, expected_volume2), (vol, expected_volume2)


class TestAxisymmetricShells(unittest.TestCase):
    def test_cquadx4(self):
        """tests a CQUADX4"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        cquadx4 = model.cquadx4
        eid = 1
        pid = 2
        mid = 3
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        cquadx4_id = model.add_cquadx4(eid, pid, [1, 2, 3, 4], theta=0.,
                                       comment='cquadx4')
        psolid_id = model.add_psolid(pid, mid, cordm=0, integ=None, stress=None,
                                     isop=None, fctn='SMECH', comment='psolid')
        E = 3.0e7
        G = None
        nu = 0.3
        mat1_id = model.add_mat1(mid, E, G, nu)
        model.setup()
        #model.cross_reference()

        cquadx4.base_nodes
        cquadx4.midside_nodes
        cquadx4.centroid()
        mass = model.mass_sum()
        area = model.area()
        assert np.allclose(mass, 0.0), mass  # TODO: not sure

        #model.uncross_reference()
        #model.safe_cross_reference()
        #model.uncross_reference()
        #bdf_file = model.write_bdf(bdf_file)

        save_load_deck(model, run_remove_unused=False, run_mass_properties=False)

    def test_axi_trax(self):
        """tests a CTRAX3/CTRAX6/???"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        ctrax3 = model.ctrax3
        ctrax6 = model.ctrax6

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        model.add_grid(5, [.5, 0., 0.])
        model.add_grid(6, [1., 0.5, 0.])
        model.add_grid(7, [.5, 1., 0.])
        model.add_grid(8, [0., .5, 0.])

        model.add_grid(9, [.5, .5, 0.])

        mid1 = 1
        E = 30.e7
        G = None
        nu = 0.3
        model.add_mat1(mid1, E, G, nu, rho=0.1)
        #model.add_mat1(2, E, G, nu, rho=0.1)
        #model.add_mat1(3, E, G, nu, rho=0.1)

        pid = 1
        nids = [1, 2, 3]
        ctrax3_id = model.add_ctrax3(1, pid, nids, theta=0., comment='ctrax3')
        #model.add_pshell(pid, mid1=2, t=0.1)

        psolid_id = model.add_psolid(pid, mid1, cordm=0, integ=None, stress=None,
                                     isop=None, fctn='SMECH', comment='psolid')

        pid = 2
        nids = [1, 2, 3, 5, 6, 9]
        ctrax6_id = model.add_ctrax6(2, pid, nids, theta=0., comment='ctrax6')

        plsolid_id = model.add_plsolid(pid, mid1, stress_strain='GRID', ge=0.,
                                       comment='plsolid')
        mathp_id = model.add_mathp(mid1)
        #assert pcomp.Thickness() == sum(thicknesses), thicknesses

        #pcomp.lam = 'SYM'
        #assert pcomp.Thickness() == sum(thicknesses)*2, thicknesses

        model.setup()
        model.validate()
        ctrax3.centroid()
        ctrax6.centroid()
        area = model.area()

        #ctrax6.raw_fields()
        model.ctrax6.write(size=8)


        #psolid.raw_fields()
        model.psolid.write(size=8)
        model.psolid.write(size=16)
        model.psolid.write(size=16, is_double=True)

        #plsolid.raw_fields()
        model.plsolid.write(size=8)
        model.plsolid.write(size=16)
        model.plsolid.write(size=16, is_double=True)

        model._verify_bdf(xref=False)

        #--------------------------------
        #model.cross_reference()
        #model._verify_bdf(xref=True)

        #ctrax3.raw_fields()
        model.ctrax3.write(size=8)

        #ctrax6.raw_fields()
        model.ctrax6.write(size=8)

        #pcomp.raw_fields()
        model.pcomp.write(size=8)
        model.pcomp.write(size=16)
        model.pcomp.write(size=16, is_double=True)
        save_load_deck(model, run_remove_unused=False, run_mass_properties=False, run_convert=False)

    def test_axi_cquadx8(self):
        """tests a CQUADX8"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        cquadx8 = model.cquadx8
        eid = 1
        pid = 10
        mid = 100
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(5, [.5, 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(6, [1., .5, 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(7, [.5, 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(8, [0., .5, 0.])
        model.add_grid(9, [.5, .5, 0.])
        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        model.add_cquadx8(eid, pid, nids, theta=0., comment='cquadx8')
        model.add_psolid(pid, mid)
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        model.setup()
        cquadx8.base_nodes
        cquadx8.midside_nodes
        cquadx8.mass()
        cquadx8.centroid()
        model.area()
        save_load_deck(model, run_mass_properties=False, run_test_bdf=False)

    def test_axi_cquadx(self):
        """tests a CQUADX"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        cquadx = model.cquadx
        eid = 1
        pid = 10
        mid = 100
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(5, [.5, 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(6, [1., .5, 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(7, [.5, 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(8, [0., .5, 0.])
        model.add_grid(9, [.5, .5, 0.])
        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        model.add_cquadx(eid, pid, nids, theta_mcid=0., comment='cquadx')
        #model.add_psolid(pid, mid)
        model.add_plplane(pid, mid, cid=0, stress_strain_output_location='GRID', comment='')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        model.setup()
        cquadx.base_nodes
        cquadx.midside_nodes
        cquadx.mass()
        cquadx.centroid()
        model.area()
        save_load_deck(model, run_remove_unused=False,
                       run_mass_properties=False, run_test_bdf=False)

    def test_axi_ctriax(self):
        """tests a CTRIAX"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        ctriax = model.ctriax
        eid = 1
        pid = 10
        mid = 100
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        #model.add_grid(4, [0., 1., 0.])
        model.add_grid(5, [.5, 0., 0.])
        model.add_grid(6, [1., .5, 0.])
        model.add_grid(9, [.5, .5, 0.])
        #model.add_psolid(pid, mid)

        # 4---7---3
        # |     / |
        # 8   9   6
        # |/      |
        # 1---5---2
        nids = [1, 2, 3, 5, 6, 9]
        model.add_ctriax(eid, pid, nids, theta_mcid=0., comment='ctriax')
        #model.add_plsolid(pid, mid, stress_strain='GRID', ge=0., comment='')
        model.add_plplane(pid, mid, cid=0, stress_strain_output_location='GRID', comment='')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        model.setup()

        ctriax.base_nodes
        ctriax.midside_nodes
        ctriax.centroid()
        #ctriax.mass()
        model.area()
        save_load_deck(model, run_remove_unused=False, run_mass_properties=False, run_test_bdf=False)

    def test_axi_ctriax6(self):
        """tests a CTRIAX6"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        ctriax6 = model.ctriax6
        #cquadx8 = model.cquadx8
        eid = 1
        pid = 10
        mid = 100
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        #model.add_grid(4, [0., 1., 0.])
        model.add_grid(5, [.5, 0., 0.])
        model.add_grid(6, [1., .5, 0.])
        #model.add_grid(7, [.5, 1., 0.])
        #model.add_grid(8, [0., .5, 0.])
        model.add_grid(9, [.5, .5, 0.])
        #model.add_psolid(pid, mid)

        # 4---7---3
        # |     / |
        # 8   9   6
        # |/      |
        # 1---5---2
        nids = [1, 5, 2, 6, 3, 9]
        model.add_ctriax6(eid, mid, nids, theta=0., comment='ctriax6')

        #model.add_plsolid(pid, mid, stress_strain='GRID', ge=0., comment='')
        model.add_plplane(pid, mid, cid=0, stress_strain_output_location='GRID', comment='')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        #model.cross_reference()

        model.setup()
        ctriax6.base_nodes
        ctriax6.midside_nodes
        #ctriax6.mass()
        ctriax6.mass()
        model.area()
        save_load_deck(model, run_mass_properties=False, run_test_bdf=False)


def make_dvcrel_optimization(model: BDF, params, element_type: str, eid: int,
                             i: int=1) -> int:
    """makes a series of DVCREL1 and a DESVAR"""
    j = i
    for ii, (name, desvar_value) in enumerate(params):
        j = i + ii
        dvids = [j]
        coeffs = [1.0]
        model.add_dvcrel1(j, element_type, eid, name, dvids, coeffs,
                          cp_min=None, cp_max=1e20,
                          c0=0.0, validate=True,
                          comment='')
        model.add_desvar(j, 'v%s' % name, desvar_value)
    return j + 1

def make_dvprel_optimization(model: BDF,
                             params: list[tuple[str | int, float]],
                             prop_type: str, pid: int,
                             i: int=1) -> int:
    """makes a series of DVPREL1 and a DESVAR"""
    j = i
    assert pid > 0, pid
    for ii, (name, desvar_value) in enumerate(params):
        j = i + ii
        desvar_ids = [j]
        coeffs = [1.0]
        model.add_dvprel1(j, prop_type, pid, name, desvar_ids, coeffs,
                          p_min=None, p_max=1e20,
                          c0=0.0, validate=True,
                          comment='')
        model.add_desvar(j, 'v%s' % name, desvar_value)
    return j + 1

def make_dvmrel_optimization(model: BDF, params, material_type: str, mid: int,
                             i: int=1) -> int:
    """makes a series of DVMREL1 and a DESVAR"""
    j = i
    for ii, (name, desvar_value) in enumerate(params):
        j = i + ii
        dvids = [j]
        coeffs = [1.0]
        model.add_dvmrel1(j, material_type, mid, name, dvids, coeffs,
                          mp_min=None, mp_max=1e20,
                          c0=0.0, validate=True,
                          comment='')
        model.add_desvar(j, 'v%s' % name, desvar_value)
    return j + 1


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
