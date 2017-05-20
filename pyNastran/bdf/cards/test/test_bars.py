from __future__ import print_function
import os
import unittest

from numpy import allclose, array

from pyNastran.bdf.bdf import BDF, BDFCard, CBAR, PBAR, PBARL, GRID, MAT1
from pyNastran.bdf.field_writer_8 import print_card_8

class TestBars(unittest.TestCase):
    """test CBAR/PBAR/PBARL classes"""
    def test_pbar_1(self):
        """tests the PBAR BDF add"""
        fields = [
            u'PBAR', 1510998, 1520998, 0.0, 4.9000000000000006e-14,
            4.9000000000000006e-14, 0.0, 0.0, None, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, None, None, 0.0
        ]
        card = print_card_8(fields)
        #print(card)
        card = print_card_8(fields)
        lines = card.split('\n')
        model = BDF(debug=False)
        card = model.process_card(lines)
        cardi = BDFCard(card)
        pbar = PBAR.add_card(cardi)
        self.assertEqual(pbar.A, 0.), pbar.A
        self.assertEqual(pbar.i12, 0.), pbar.i12
        self.assertEqual(pbar.k1, None), pbar.k1
        self.assertEqual(pbar.k2, None), pbar.k2
        #with self.assertRaises(AssertionError):  # A=0, I12=0, K1=0
            #pbar = PBAR(card2)

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
        card = model.process_card(lines)
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
        card = model.process_card(lines)

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
                    k2=1.e8, comment='cat')
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
        #pbarl.group = 'MSCBMLO'

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
        cbar.validate()
        model.elements[eid] = cbar
        pbarl._verify(xref=False)

        model.validate()
        model.cross_reference()
        pbarl._verify(xref=True)
        assert allclose(cbar.Mass(), 9.9247779608), cbar.Mass()

        mat.rho = 0.
        assert allclose(cbar.Mass(), 0.5), cbar.Mass()


    def test_bar_mass_1(self):
        """tests CBAR/PBAR mass"""
        model = BDF(debug=False)
        #model.case_control_deck = CaseControlDeck(case_control_lines)
        spc = ['SPC1', 123456, 123456, 1]
        grid1 = ['GRID', 1, None, 0., 0., 0.]
        grid2 = ['GRID', 2, None, 1., 0., 0.]
        grid3 = ['GRID', 3, None, 1., 0., 0.]
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

        mass, cg, I = model.mass_properties(
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
        assert allclose(cbar.Length(), 1.0), cbar.Length()
        #assert allclose(cbar.Mass(), 10.25), cbar.Mass()
        #assert allclose(cbar.MassPerLength(), 10.25), cbar.MassPerLength()
        #assert allclose(mass, 10.25), mass

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

        if 0:
            from pyNastran.op2.op2 import OP2
            op2 = OP2()
            op2.read_op2('cbar.op2')
            #os.remove('cbar.op2')
            gpw = op2.grid_point_weight
            op2_mass = gpw.mass.max()
            assert allclose(op2_mass, mass), 'op2_mass=%s mass=%s' % (op2_mass, mass)
            #print('op2_mass=%s mass=%s' % (op2_mass, mass))
            op2_cg = gpw.cg

            cg = array([0.5, 0., 0.], dtype='float32')
            #print('cg =', op2_cg)

    def test_bar_mass_2(self):
        """CBAR/PBARL"""
        model = BDF(debug=False)
        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(3, xyz=[0., 1., 0.])

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
        cbar = model.add_cbar(eid, pid, nids, x, g0, offt='GGG',
                              pa=0, pb=0, wa=None, wb=None,
                              comment='CBAR')
        Type = 'BOX'
        dim = [1., 2., 0.1, 0.1]
        #pbeaml = model.add_pbeaml(pid, mid, Type, xxb, dims, nsm=None,
                                  #so=None, comment='PBEAML')
        pbarl = model.add_pbarl(pid, mid, Type, dim, group='MSCBMLO', nsm=0.,
                                comment='PBARL')
        #---------------------------------------------------------------
        eid = 2
        pid = 102
        x = None
        g0 = 3
        cbar = model.add_cbar(eid, pid, nids, x, g0, offt='GGG',
                              pa=0, pb=0, wa=None, wb=None,
                              comment='CBAR')
        Type = 'BOX'
        dim = [1., 2., 0.1, 0.1]
        pbarl = model.add_pbarl(pid, mid, Type, dim, group='MSCBMLO', nsm=0.,
                                comment='PBARL')
        #---------------------------------------------------------------
        model.validate()
        model._verify_bdf(xref=False)
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()

        model._verify_bdf(xref=True)
        model.uncross_reference()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
