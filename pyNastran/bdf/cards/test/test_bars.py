import os
import unittest

import numpy as np

from pyNastran.bdf.bdf import BDF, BDFCard, CBAR, PBAR, PBARL, GRID, MAT1
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.mesh_utils.mass_properties import (
    mass_properties, mass_properties_nsm)  #mass_properties_breakdown
from pyNastran.bdf.cards.test.utils import save_load_deck


class TestBars(unittest.TestCase):
    """test CBAR/PBAR/PBARL classes"""
    def test_pbar_1(self):
        """tests the PBAR BDF add"""
        area = 0.0
        i11 = 4.9e-2
        i22 = 5.5e-2
        i12 = 6.6e-2
        j = 7.7e-2
        nsm = 1.0
        fields = [
            u'PBAR', 1510998, 1520998, area, i11,
            i22, j, nsm, None, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, None, None, i12,
        ]
        card = print_card_8(fields)
        #print(card)
        card = print_card_8(fields)
        lines = card.split('\n')
        model = BDF(debug=False)
        card = model._process_card(lines)
        cardi = BDFCard(card)
        pbar = PBAR.add_card(cardi)
        pbar.raw_fields()
        self.assertEqual(pbar.A, area)
        self.assertEqual(pbar.i1, i11)
        self.assertEqual(pbar.i2, i22)
        self.assertEqual(pbar.i12, i12)
        self.assertEqual(pbar.j, j)
        self.assertEqual(pbar.k1, None)
        self.assertEqual(pbar.k2, None)
        self.assertEqual(pbar.nsm, nsm)
        assert np.allclose(pbar.Area(), area)
        assert np.allclose(pbar.I11(), i11)
        assert np.allclose(pbar.I22(), i22)
        assert np.allclose(pbar.I12(), i12)
        assert np.allclose(pbar.J(), j)
        assert np.allclose(pbar.Nsm(), nsm)

    def test_pbar_2(self):
        """tests the PBAR BDF add"""
        pid = 1
        mid = 2
        A = None
        I1 = I2 = None
        J = None
        nsm = None
        c1 = c2 = d1 = d2 = e1 = e2 = f1 = f2 = None
        k1 = k2 = None
        i12 = 3.
        fields = [
            'PBAR', pid, mid, A, I1, I2, J, nsm, None,
            c1, c2, d1, d2, e1, e2, f1, f2,
            k1, k2, i12
        ]
        card = print_card_8(fields)
        lines = card.split('\n')
        model = BDF(debug=False)
        card = model._process_card(lines)
        cardi = BDFCard(card)

        pbar = PBAR.add_card(cardi)
        self.assertEqual(pbar.pid, 1)
        self.assertEqual(pbar.mid, 2)
        self.assertEqual(pbar.A, 0.0)
        self.assertEqual(pbar.i1, 0.0)
        self.assertEqual(pbar.i2, 0.0)
        self.assertEqual(pbar.j, 0.0)
        self.assertEqual(pbar.nsm, 0.0)
        self.assertEqual(pbar.i12, 3.0)
        self.assertEqual(pbar.c1, 0.0)
        self.assertEqual(pbar.c2, 0.0)
        self.assertEqual(pbar.d1, 0.0)
        self.assertEqual(pbar.d2, 0.0)
        self.assertEqual(pbar.e1, 0.0)
        self.assertEqual(pbar.e2, 0.0)
        self.assertEqual(pbar.k1, None)
        self.assertEqual(pbar.k2, None)

        #--------------------------------------------------------
        A = 6.
        I1 = 5.
        I2 = 4.
        J = 3.
        nsm = 2.
        c1 = c2 = d1 = d2 = e1 = e2 = f1 = f2 = None
        k1 = k2 = 1e2
        i12 = 0.
        fields = [
            'PBAR', pid, mid, A, I1, I2, J, nsm, None,
            c1, c2, d1, d2, e1, e2, f1, f2,
            k1, k2, i12]
        card = print_card_8(fields)
        lines = card.split('\n')
        model = BDF(debug=False)
        card = model._process_card(lines)

        cardi = BDFCard(card)
        pbar = PBAR.add_card(cardi)

        self.assertEqual(pbar.pid, 1)
        self.assertEqual(pbar.mid, 2)
        self.assertEqual(pbar.A, 6.0)
        self.assertEqual(pbar.i1, 5.0)
        self.assertEqual(pbar.i2, 4.0)
        self.assertEqual(pbar.j, 3.0)
        self.assertEqual(pbar.nsm, 2.0)
        self.assertEqual(pbar.i12, 0.0)
        self.assertEqual(pbar.c1, 0.0)
        self.assertEqual(pbar.c2, 0.0)
        self.assertEqual(pbar.d1, 0.0)
        self.assertEqual(pbar.d2, 0.0)
        self.assertEqual(pbar.e1, 0.0)
        self.assertEqual(pbar.e2, 0.0)
        self.assertEqual(pbar.k1, 1e2)
        self.assertEqual(pbar.k2, 1e2)

    def test_pbar_3(self):
        """tests the PBAR validate"""
        pid = 42
        mid = 10
        i1 = -1.
        i2 = -2.
        i12 = -3.
        j = -4.
        pbar = PBAR(pid, mid, A=0., i1=i1, i2=i2, i12=i12, j=j, nsm=0., c1=0., c2=0.,
                    d1=0., d2=0., e1=0., e2=0., f1=0., f2=0., k1=1.e8,
                    k2=1.e8, comment='pbar')
        with self.assertRaises(ValueError):
            pbar.validate()

        pbar.i1 = 1.
        with self.assertRaises(ValueError):
            pbar.validate()

        pbar.i2 = 2.
        with self.assertRaises(ValueError):
            pbar.validate()

        pbar.j = 4.
        pbar.validate()

        model = BDF(debug=False)
        pbar = model.add_pbar(pid, mid, A=0., i1=2., i2=2., i12=1., j=4., nsm=0., c1=0., c2=0.,
                              d1=0., d2=0., e1=0., e2=0., f1=0., f2=0., k1=1.e8,
                              k2=1.e8, comment='pbar')
        pbar.validate()
        nids = [100, 101]
        eid = 1000
        x = [0., 0., 1.]
        g0 = None
        model.add_cbar(eid, pid, nids, x, g0, comment='cbar')
        model.add_grid(100, [0., 0., 0.])
        model.add_grid(101, [1., 0., 0.])
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        save_load_deck(model)

    def test_cbar_g0(self):
        """modification of test_cbeam_01"""
        model = BDF(debug=False)
        pid = 200
        mid = 6
        model.add_pbar(pid, mid, A=0., i1=2., i2=2., i12=1., j=4., nsm=0., c1=0., c2=0.,
                       d1=0., d2=0., e1=0., e2=0., f1=0., f2=0., k1=1.e8,
                       k2=1.e8, comment='pbar')

        eid = 100
        nids = [10, 20]
        x = None
        g0 = 30
        cbar = model.add_cbar(eid, pid, nids, x, g0, comment='cbar')
        cbar.write_card_16(is_double=False)

        E = 1.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(20, [0., 1., 0.])
        model.add_grid(30, [0., 2., 0.])
        model.cross_reference()

        save_load_deck(model)

    def test_pbarl_1(self):
        """tests the PBARL"""
        model = BDF(log=None, debug=False)
        pid = 4
        mid = 40
        group = 'group'
        Type = 'bad_type'
        dim = 42
        nsm = 0.5
        pbarl = PBARL(pid, mid, Type, dim, group=group, nsm=nsm, comment='comment')
        with self.assertRaises(ValueError): # Type
            pbarl.validate()

        pbarl.Type = 'TUBE'
        with self.assertRaises(TypeError): # dim
            pbarl.validate()

        pbarl.dim = [20.]
        with self.assertRaises(RuntimeError):
            pbarl.validate()

        pbarl.dim = [2., 1.]
        #with self.assertRaises(ValueError):
            #pbarl.validate()
        #pbarl.group = 'MSCBML0'

        pbarl.validate()
        str(pbarl)
        pbarl.write_card(size=8, is_double=False)
        pbarl.write_card(size=16, is_double=False)
        pbarl.write_card(size=16, is_double=True)
        model.properties[pid] = pbarl

        nid1 = 52
        xyz1 = [0., 0., 0.]
        model.nodes[nid1] = GRID(nid1, cp=0, xyz=xyz1)

        nid2 = 53
        xyz2 = [1., 0., 0.]
        model.nodes[nid2] = GRID(nid2, cp=0, xyz=xyz2)

        E = 30.0e7
        G = None
        nu = 0.3
        mat = MAT1(mid, E, G, nu, rho=1.0)
        model.materials[mid] = mat

        eid = 42
        x = None
        g0 = None
        cbar = CBAR(eid, pid, [nid1, nid2], x, g0, offt='GGG',
                    pa=0, pb=0, wa=None, wb=None, comment='')
        with self.assertRaises(ValueError):
            cbar.validate()
        cbar.x = [0., 1., 2.]
        cbar.validate()
        model.elements[eid] = cbar
        pbarl._verify(xref=False)

        model.validate()
        model.cross_reference()
        pbarl._verify(xref=True)
        assert np.allclose(cbar.Mass(), 9.9247779608), cbar.Mass()

        mat.rho = 0.
        assert np.allclose(cbar.Mass(), 0.5), cbar.Mass()

        scale = 'FR'
        x = [0.2, 0.4, 0.6, 0.8]
        model.add_cbarao(eid, scale, x, comment='cbarao')
        model.add_card(['CBARAO', eid+1, 'RF', 6, 0.1, 0.2], 'CBARAO')
        save_load_deck(model, run_quality=False, run_test_bdf=False)

    def test_bar_mass_1(self):
        """tests CBAR/PBAR mass"""
        model = BDF(debug=False)
        #model.case_control_deck = CaseControlDeck(case_control_lines)
        spc = ['SPC1', 123456, 123456, 1]
        grid1 = ['GRID', 1, None, 0., 0., 0.]
        grid2 = ['GRID', 2, None, 1., 0., 0.]
        #grid3 = ['GRID', 3, None, 1., 0., 0.]
        force = ['FORCE', 100, 1, 0, 2., 3., 4.]
        pid = 11
        mid = 12
        cbar = [
            'CBAR', 10, pid, 1, 2, 0., 1., 0., None,
        ]
        k1 = k2 = None
        area = 2.0
        rho = 3.
        nu = 0.3
        i1 = 2.1
        i2 = 1.2
        i12 = 0.1
        j = None
        nsm = 0.1
        pbar = [
            'PBAR', pid, mid, area, i1, i2, j, nsm,
            None, None, None, None, None, None, None, None,
            k1, k2, i12
        ]

        mat1 = ['MAT1', mid, 3.0e7, None, nu, rho]
        model.add_card(grid1, 'GRID')
        model.add_card(grid2, 'GRID')
        model.add_card(cbar, 'CBAR')
        model.add_card(pbar, 'PBAR')
        model.add_card(mat1, 'MAT1')
        model.add_card(spc, 'SPC1')
        model.add_card(force, 'FORCE')
        model.validate()
        model.cross_reference()

        mass, unused_cg, unused_I = mass_properties(
            model,
            element_ids=None, mass_ids=None,
            reference_point=None,
            sym_axis=None,
            scale=None)
        #print('cg* =', cg)
        L = 1.0
        mass_per_length = area * rho + nsm
        mass = L * mass_per_length

        #xcg = (0.0 * mass_a + 1.0 * mass_b) / (mass_a + mass_b)
        #print(mass_a, mass_b, xcg, mass_a + mass_b)
        #print('mass =', mass)
        #cbar = CBEAM()
        cbar = model.elements[10]
        pbar = model.properties[11]
        assert pbar.Nu() == nu, 'pbar.Nu()=%s nu=%s' % (pbar.Nu(), nu)
        assert pbar.Rho() == rho, 'pbar.Rho()=%s rho=%s' % (pbar.Rho(), rho)
        assert np.allclose(cbar.Length(), 1.0), cbar.Length()
        #assert np.allclose(cbar.Mass(), 10.25), cbar.Mass()
        #assert np.allclose(cbar.MassPerLength(), 10.25), cbar.MassPerLength()
        #assert np.allclose(mass, 10.25), mass

        case_control_lines = (
            'SOL 101\n'
            'CEND\n'
            'SUBCASE 1\n'
            '    STRESS(PLOT,SORT1,REAL) = ALL\n'
            '    SPC = 123456\n'
            '    LOAD = 100\n'
            'BEGIN BULK\n'
            'PARAM,GRDPNT,0\n'
            'PARAM,POST,-1\n'
            'PARAM   POSTEXT YES\n'
        )
        with open('cbar.bdf', 'w') as bdf_file:
            bdf_file.write(case_control_lines)
            model.write_bdf(bdf_file, enddata=True)
        model2 = BDF(debug=False)
        model2.read_bdf('cbar.bdf')

        model2._verify_bdf(xref=True)
        if not os.path.exists('cbar.op2') and 0:
            os.system('nastran scr=yes bat=no old=no cbar.bdf')
        os.remove('cbar.bdf')

        if 0:  # pragma: no cover
            from pyNastran.op2.op2 import OP2
            op2 = OP2()
            op2.read_op2('cbar.op2')
            #os.remove('cbar.op2')
            gpw = op2.grid_point_weight
            op2_mass = gpw.mass.max()
            assert np.allclose(op2_mass, mass), 'op2_mass=%s mass=%s' % (op2_mass, mass)
            #print('op2_mass=%s mass=%s' % (op2_mass, mass))
            unused_op2_cg = gpw.cg

            unused_cg = np.array([0.5, 0., 0.], dtype='float32')
            #print('cg =', op2_cg)

    def test_bar_mass_2(self):
        """CBAR/PBARL"""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0., 1., 0.])

        mid = 1
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.)

        #---------------------------------------------------------------
        eid = 1
        pid = 101
        nids = [1, 2]
        x = [0., 0., 1.]
        g0 = None
        unused_cbar = model.add_cbar(eid, pid, nids, x, g0, offt='GGG',
                                     pa=0, pb=0, wa=None, wb=None,
                                     comment='CBAR')
        Type = 'BOX'
        dim = [1., 2., 0.1, 0.1]
        #pbeaml = model.add_pbeaml(pid, mid, Type, xxb, dims, nsm=None,
                                  #so=None, comment='PBEAML')
        unused_pbarl = model.add_pbarl(pid, mid, Type, dim, group='MSCBML0', nsm=0.,
                                       comment='PBARL')
        #---------------------------------------------------------------
        eid = 2
        pid = 102
        x = None
        g0 = 3
        unused_cbar = model.add_cbar(eid, pid, nids, x, g0, offt='GGG',
                                     pa=0, pb=0, wa=None, wb=None,
                                     comment='CBAR')
        Type = 'BOX'
        dim = [1., 2., 0.1, 0.1]
        unused_pbarl = model.add_pbarl(pid, mid, Type, dim, group='MSCBML0', nsm=0.,
                                       comment='PBARL')
        #---------------------------------------------------------------
        eid = 3
        pid = 103
        #cbar = model.add_cbar(eid, pid, nids, x, g0, offt='GGG',
                              #pa=42, pb=5, wa=None, wb=None,
                              #comment='CBAR')
        unused_pbar = model.add_pbar(pid, mid, A=1., i1=0., i2=0., i12=0., j=0., nsm=0.1,
                                     c1=0., c2=0.,
                                     d1=0., d2=0.,
                                     e1=0., e2=0.,
                                     f1=0., f2=0.,
                                     k1=1.e8, k2=1.e8,
                                     comment='pbar')

        #G = 3.0e7
        #E = None
        #nu = 0.3
        #model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0,
                       #St=0.0, Sc=0.0, Ss=0.0, mcsid=0,
                       #comment='')
        #---------------------------------------------------------------
        model.validate()
        model.pop_parse_errors()
        model._verify_bdf(xref=False)

        model.cross_reference()
        model.pop_xref_errors()

        model._verify_bdf(xref=True)
        model.uncross_reference()

    def test_pbar_nsm(self):
        model = BDF(debug=False)
        pid = 1
        mid = 1
        nsm = 1.
        area = 2.0
        pbar = model.add_pbar(pid, mid, A=area, i1=0., i2=0., i12=0., j=0.,
                              nsm=nsm,
                              c1=0., c2=0., d1=0., d2=0.,
                              e1=0., e2=0., f1=0., f2=0.,
                              k1=1.e8, k2=1.e8,
                              comment='')

        E = 1.0
        G = None
        nu = 0.3
        mat1 = model.add_mat1(mid, E, G, nu)

        #----------------
        card_lines = [
            'PBAR           2       1      2.                              1.',
        ]
        model.add_card(card_lines, 'PBAR', comment='', is_list=False,
                       has_none=True)
        pbar2 = model.properties[2]
        #------------------
        model.cross_reference()

        assert pbar.Nsm() == 1.0
        assert pbar.Area() == 2.0

        # mass/L = area*rho + nsm
        assert pbar.MassPerLength() == 1.0

        # area = 2.0
        mat1.rho = 10.0
        assert pbar.MassPerLength() == 21.0, pbar.MassPerLength()
        assert pbar2.MassPerLength() == 21.0, pbar2.MassPerLength()

    def test_pbarl_nsm(self):
        model = BDF(debug=False)
        pid = 1
        mid = 1
        bar_type = 'BAR'
        dim = [1., 2.]  # area = 2.0
        pbarl = model.add_pbarl(pid, mid, bar_type, dim, group='MSCBML0', nsm=1.,
                                comment='')

        E = 1.0
        G = None
        nu = 0.3
        mat1 = model.add_mat1(mid, E, G, nu)

        #----------------
        card_lines = [
            'PBARL   2       1               BAR',
            '        1.0     2.0      1.0',
        ]
        model.add_card(card_lines, 'PBARL', comment='', is_list=False,
                       has_none=True)
        pbarl2 = model.properties[2]
        #------------------
        model.cross_reference()

        assert pbarl.Nsm() == 1.0
        assert pbarl.Area() == 2.0

        # mass/L = area*rho + nsm
        assert pbarl.MassPerLength() == 1.0

        # area = 2.0
        mat1.rho = 10.0
        assert pbarl.MassPerLength() == 21.0, pbarl.MassPerLength()
        assert pbarl2.MassPerLength() == 21.0, pbarl2.MassPerLength()

        loadcase_id = 10
        eid = 11
        load_type = 'FZ'
        x1 = 0.
        x2 = None
        p1 = 10.
        scale = 'FR'
        model.add_pload1(loadcase_id, eid, load_type, scale, x1, p1,
                         x2=x2, p2=None, comment='pload1')

        scale = 'LE'
        model.add_pload1(loadcase_id, eid, load_type, scale, x1, p1,
                         x2=x2, p2=None, comment='')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0., 1., 0.])
        x = None
        g0 = 3
        model.add_cbar(eid, pid, [1, 2], x, g0)
        model.cross_reference()

        p0 = 1
        eids = None
        nids = None
        force1, moment1 = model.sum_forces_moments(p0, loadcase_id,
                                                   include_grav=False, xyz_cid0=None)
        force2, moment2 = model.sum_forces_moments_elements(p0, loadcase_id, eids, nids,
                                                            include_grav=False, xyz_cid0=None)
        #print(force1, force2)
        assert np.allclose(force1, force2), force1
        assert np.allclose(moment1, moment2), moment1
        save_load_deck(model, xref='standard', punch=True)

    def test_baror(self):
        """tests a BAROR"""
        model = BDF(debug=False)
        n1 = 10
        n2 = 20
        model.add_grid(n1, [0., 0., 0.])
        model.add_grid(n2, [1., 0., 0.])

        pid = 2
        mid = 1
        bar_type = 'BAR'
        dim = [1., 2.]  # area = 2.0
        unused_pbarl = model.add_pbarl(pid, mid, bar_type, dim, group='MSCBML0', nsm=1.,
                                       comment='')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.)

        card_lines = ['BAROR', None, pid, None, None, 0.6, 2.9, -5.87, 'GOG']
        model.add_card(card_lines, 'BAROR', comment='BAROR', is_list=True,
                       has_none=True)

        eid = 1
        card_lines = ['CBAR', eid, pid, n1, n2]
        model.add_card(card_lines, 'CBAR', comment='', is_list=True, has_none=True)
        model.pop_parse_errors()
        save_load_deck(model)

    def test_baror_2(self):
        model = BDF(debug=False)
        pid = 12
        is_g0 = True
        g0 = 42
        x = None
        baror = model.add_baror(pid, is_g0, g0, x, offt='GGG', comment='baror')
        baror.raw_fields()
        baror.write_card(size=8)
        baror.write_card(size=16)
        save_load_deck(model)

    def test_cbend(self):
        """tests a CBEND"""
        model = BDF(debug=False)

        eid = 7
        pid = 10
        nids = [2, 3]
        g0 = 5
        x = None
        geom = 1
        cbend = model.add_cbend(eid, pid, nids, g0, x, geom, comment='cbend')
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        model.add_grid(5, [0., 0., 0.])
        #pbend = model.add_pbend(pid, mid, beam_type, A, i1, i2, j,
                                #c1, c2, d1, d2, e1, e2, f1, f2,
                                #k1, k2, nsm, rc, zc, delta_n, fsi,
                                #rm, t, p, rb, theta_b, comment='')

        cbend.validate()
        cbend.raw_fields()
        cbend.write_card()
        cbend.write_card(size=16)

        model.validate()
        model._verify_bdf(xref=False)
        model.pop_parse_errors()
        #model.cross_reference()
        #model.pop_xref_errors()

        #model._verify_bdf(xref=True)
        #model.uncross_reference()

    def test_cbeam3(self):
        """tests a CBEAM3"""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        model.add_grid(4, [0., 0., 0.])
        eid = 1
        pid = 2
        nids = [1, 2, 3]
        x = None
        g0 = 4
        cbeam3 = model.add_cbeam3(eid, pid, nids, x, g0, wa=None, wb=None, wc=None, tw=None, s=None, comment='cbeam3')
        cbeam3.raw_fields()

        A = 1.
        iz = 2.
        iy = 3.
        mid = 4
        pbeam3 = model.add_pbeam3(pid, mid, A, iz, iy, iyz=0., j=None,
                                  nsm=0., cy=0., cz=0., dy=0., dz=0., ey=0., ez=0., fy=0., fz=0., comment='')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1)
        str(cbeam3)
        pbeam3s = str(pbeam3)
        #print(pbeam3s)
        str(pbeam3s)
        card_lines = pbeam3s.split('\n')

        cbeam3._verify(xref=False)
        model.cross_reference()
        model.uncross_reference()

        del model.properties[pid]
        model.cards_to_read.add('PBEAM3')
        model.add_card(card_lines, 'PBEAM3', comment='', ifile=None, is_list=False, has_none=True)
        model.pop_parse_errors()
        model.pop_xref_errors()
        assert pbeam3 == model.properties[pid]

    def test_bar_area(self):
        """tests the PBARL"""
        model = BDF(log=None, debug=False)
        mid = 40
        group = 'group'
        nsm = 0.0

        shape_dims_area = [
            # name, dims, area, i1
            ('ROD', [2.], 4. * np.pi, 0.),
            ('TUBE', [5., 1.], 24. * np.pi, 0.),
            ('BAR', [2., 3.], 6., 0.),
            ('BOX', [2., 3., 0.5, 0.5], 4., 0.),

            ('L', [2., 3., 1., 1.], 4., 0.),
            ('CHAN', [10., 10., 1., 1.], 28., None),
            ('CHAN1', [9., 0.1, 8., 10.], 19., None),
            ('CHAN2', [1, 1., 9., 10.], 26., None),

            # new
            ('I', [1., 1., 1., 0.1, 0.1, 0.1], 0.28, None),
            ('I1', [0.1, 1., 0.5, 1.], 1.05, None),
            ('H', [1.0, 0.1, 1.0, 0.1], 0.2, None),

            ('Z', [0.5, 0.5, 0.5, 1.], 0.75, None),
            ('Z', [0.8, 0.5, 0.5, 1.], 0.90, None),
            ('Z', [0.5, 0.8, 0.5, 1.], 1.05, None),
            ('Z', [0.5, 0.5, 0.8, 1.], 0.60, None),
            ('Z', [0.5, 0.5, 0.5, 2.], 1.75, None),

            ('CHAN', [1., 1., 0.1, 0.1], 0.28, None),
            ('CHAN1', [0.5, 0.5, 0.5, 1.], 0.75, None),
            ('CHAN2', [0.1, 0.1, 1., 1.], 0.28, None),
            ('CROSS', [0.1, 0.1, 1., 0.1], 0.11, None),
            ('HEXA', [0.1, 1., 1.], 0.90, None),
            ('HEXA', [0.2, 1., 1.], 0.80, None),
            ('HEXA', [0.1, 2., 1.], 1.90, None),
            ('HEXA', [0.1, 1., 2.], 1.80, None),

            ('HAT', [1., 0.1, 1., 0.1], 0.30, None),
            ('HAT', [2., 0.1, 1., 0.1], 0.50, None),
            ('HAT', [1., 0.2, 1., 0.1], 0.56, None),
            ('HAT', [1., 0.1, 2., 0.1], 0.40, None),
            ('HAT', [1., 0.1, 1., 0.2], 0.32, None),

            ('HAT1', [3., 1., 1., 0.1, 0.1], 0.76, None),
            ('HAT1', [3., 2., 1., 0.1, 0.1], 0.96, None),
            ('HAT1', [3., 1., 2., 0.1, 0.1], 0.76, None),
            ('HAT1', [3., 1., 1., 0.2, 0.1], 1.18, None),
            ('HAT1', [3., 1., 1., 0.1, 0.2], 1.04, None),

            ('T', [10., 10., 3., 0.5], 33.5, None),
            ('T2', [10., 5., 0.5, 2.0], 14., None), #  ball,hall,tflange,tweb

            ('T', [1., 1., 0.1, 0.1], 0.19, None),
            ('T1', [1., 1., 0.1, 0.1], 0.20, None),
            ('T2', [1., 1., 0.1, 0.1], 0.19, None),

            ('DBOX', [2., 1., 1., 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, ], 0.64, None),
            ('DBOX', [2., 2., 1., 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, ], 0.94, None),
            ('DBOX', [2., 1., 2., 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, ], 0.64, None),

            ('DBOX', [2., 1., 1., 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, ], 0.72, None),
            ('DBOX', [2., 1., 1., 0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, ], 0.72, None),
            ('DBOX', [2., 1., 1., 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1, ], 0.72, None),
            ('DBOX', [2., 1., 1., 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, ], 0.725, None),
            ('DBOX', [2., 1., 1., 0.1, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, ], 0.725, None),
            ('DBOX', [2., 1., 1., 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.1, ], 0.725, None),
            ('DBOX', [2., 1., 1., 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, ], 0.725, None),
        ]
        pid = 1
        for bar_type, dims, areai, i1 in shape_dims_area:
            pbarl = PBARL(pid, mid, bar_type, dims, group=group, nsm=nsm, comment='comment')
            pbarl.validate()
            area2 = pbarl.Area()
            if i1 is not None:
                pbarl.I1()
                pbarl.I2()
                pbarl.I12()
            assert np.allclose(areai, area2), 'bar_type=%r dims=%s area=%s area_expected=%s' % (bar_type, dims, area2, areai)
            pid += 1

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
