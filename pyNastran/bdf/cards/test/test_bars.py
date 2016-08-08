from __future__ import print_function
from six.moves import zip, StringIO
import os
import unittest

from numpy import allclose, array

from pyNastran.bdf.bdf import BDF, BDFCard, PBAR #, GRID, MAT1
from pyNastran.bdf.field_writer_8 import print_card_8

bdf = BDF(debug=False)
class TestBars(unittest.TestCase):
    def test_pbar_01(self):
        fields = [
            u'PBAR', 1510998, 1520998, 0.0, 4.9000000000000006e-14,
            4.9000000000000006e-14, 0.0, 0.0, None, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, None, None, 0.0
        ]
        card = print_card_8(fields)
        #print(card)
        card = print_card_8(fields)
        lines = card.split('\n')
        card = bdf.process_card(lines)
        cardi = BDFCard(card)
        pbar = PBAR.add_card(cardi)
        self.assertEqual(pbar.A, 0.), pbar.A
        self.assertEqual(pbar.i12, 0.), pbar.i12
        self.assertEqual(pbar.k1, None), pbar.k1
        self.assertEqual(pbar.k2, None), pbar.k2
        #with self.assertRaises(AssertionError):  # A=0, I12=0, K1=0
            #pbar = PBAR(card2)

    def test_pbar_02(self):
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
        card = bdf.process_card(lines)
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
        card = bdf.process_card(lines)

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

    def test_bar_mass_01(self):
        model = BDF(debug=False)
        #model.case_control_deck = CaseControlDeck(case_control_lines)
        spc = ['SPC1', 123456, 123456, 1]
        grid1 = ['GRID', 1, None, 0., 0., 0.]
        grid2 = ['GRID', 2, None, 1., 0., 0.]
        grid3 = ['GRID', 3, None, 1., 0., 0.]
        force = ['FORCE', 100, 1, 0, 2., 3., 4.]
        cbar = [
            'CAR', 10, 11, 1, 2, 0., 1., 0., None,
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
            'PBAR', 11, 12, area, i1, i2, j, nsm,
            None, None, None, None, None, None, None, None,
            k1, k2, i12
        ]

        mat1 = ['MAT1', 12, 3.0e7, None, nu, rho]
        model.add_card(grid1, 'GRID')
        model.add_card(grid2, 'GRID')
        model.add_card(cbar, 'CBAR')
        model.add_card(pbar, 'PBAR')
        model.add_card(mat1, 'MAT1')
        model.add_card(spc, 'SPC1')
        model.add_card(force, 'FORCE')
        model.cross_reference()

        mass, cg, I = model.mass_properties(
            element_ids=None, mass_ids=None,
            reference_point=None,
            sym_axis=None,
            num_cpus=1,
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

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
