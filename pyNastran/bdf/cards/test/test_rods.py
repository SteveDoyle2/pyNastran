from io import StringIO
from math import pi, sqrt
import unittest

import numpy as np

from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.bdf import CROD, CONROD, PROD, CTUBE, PTUBE, GRID, MAT1
from pyNastran.bdf.cards.test.test_shells import make_dvprel_optimization
from pyNastran.bdf.cards.test.utils import save_load_deck

#from pyNastran.bdf.field_writer_8 import print_card_8

bdf = BDF(debug=False)
class TestRods(unittest.TestCase):
    def test_crod_01(self):
        lines = ['CROD          10     100      10       2']
        card = bdf._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = CROD.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()
        self.assertEqual(card.eid, 10)
        self.assertEqual(card.Pid(), 100)
        node_ids = card.node_ids
        assert node_ids == [10, 2], node_ids # probably wrong

    def test_prod_nsm(self):
        model = BDF(debug=False)
        pid = 1
        mid = 1
        nsm = 1.
        A = 2.0
        prod = model.add_prod(pid, mid, A, j=0., c=0., nsm=nsm, comment='')
        #print(prod)

        E = 1.0
        G = None
        nu = 0.3
        mat1 = model.add_mat1(mid, E, G, nu)
        #----------------
        card_lines = [
            'PROD           2       1      2.                      1.',
        ]
        model.add_card(card_lines, 'PROD', comment='', is_list=False,
                       has_none=True)
        model.pop_parse_errors()
        prod2 = model.properties[2]
        #------------------
        model.cross_reference()

        assert prod.Nsm() == 1.0
        assert prod.Area() == 2.0

        # mass/L = area*rho + nsm
        assert prod.MassPerLength() == 1.0

        # area = 2.0
        mat1.rho = 10.0
        assert prod.MassPerLength() == 21.0, prod.MassPerLength()
        assert prod2.MassPerLength() == 21.0, prod2.MassPerLength()

    def test_ptube_nsm(self):
        """tests a PTUBE for the NSM field"""
        model = BDF(debug=False)
        pid = 1
        mid = 1
        nsm = 1.
        OD1 = 2.0
        ptube = model.add_ptube(pid, mid, OD1, t=None, nsm=nsm, OD2=None,
                                comment='ptube')
        #print(ptube)

        E = 1.0
        G = None
        nu = 0.3
        rho = 0.0
        mat1 = model.add_mat1(mid, E, G, nu, rho=rho)
        #----------------
        card_lines = [
            'PTUBE          2       1      2.              1.',
        ]
        model.add_card(card_lines, 'PTUBE', comment='', is_list=False,
                       has_none=True)
        model.pop_parse_errors()
        ptube2 = model.properties[2]
        #print(ptube.get_stats())
        #------------------
        model.cross_reference()

        assert ptube.Nsm() == 1.0

        area = pi
        mpl = area * rho + nsm
        assert ptube.Area() == area, ptube.Area()

        # mass/L = area*rho + nsm
        assert ptube.MassPerLength() == 1.0, ptube.MassPerLength()

        rho = 10.0
        mat1.rho = rho
        mpl = area * rho + nsm
        assert ptube2.MassPerLength() == mpl, ptube2.MassPerLength()

        eid = 100
        nids = [10, 12]
        ctube = model.add_ctube(eid, pid, nids, comment='ctube')
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(12, [2., 0., 0.])

        model.safe_cross_reference()
        model.uncross_reference()
        model.cross_reference()
        mass, cg, inertia = model.mass_properties(
            element_ids=None, mass_ids=None,
            reference_point=None, sym_axis=None,
            scale=None, inertia_reference='cg')

        # L * (rho * a + nsm)
        A = pi * (OD1**2) / 4.
        mass_expected = 2.0 * (rho * A + nsm)
        assert np.allclose(mass, mass_expected), mass
        save_load_deck(model)

    def test_conrod_01(self):
        eid = 10
        nid1 = 2
        nid2 = 3
        mid = 5
        A = 27.0
        lines = ['conrod,%i, %i, %i, %i, %f' % (eid, nid1, nid2, mid, A)]
        card = bdf._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = CONROD.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()
        self.assertEqual(card.eid, eid)
        self.assertEqual(card.Mid(), mid)
        node_ids = card.node_ids
        assert node_ids == [nid1, nid2], node_ids

    def test_conrod_nsm(self):
        model = BDF(debug=False)
        mid = 1
        nsm = 1.
        A = 2.0
        eid = 1
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        conrod = model.add_conrod(eid, mid, nids, A=A, j=0.0, c=0.0, nsm=nsm,
                               comment='')

        E = 1.0
        G = None
        nu = 0.3
        mat1 = model.add_mat1(mid, E, G, nu)
        #----------------
        card_lines = [
            'CONROD         2       1       2       1      2.                      1.',
        ]
        model.add_card(card_lines, 'CONROD', comment='', is_list=False,
                       has_none=True)
        model.pop_parse_errors()
        conrod2 = model.elements[2]

        #------------------
        model.cross_reference()

        assert conrod.Nsm() == 1.0
        assert conrod.Area() == 2.0

        # mass/L = area*rho + nsm
        assert conrod.MassPerLength() == 1.0

        # area = 2.0
        mat1.rho = 10.0
        assert conrod.MassPerLength() == 21.0, conrod.MassPerLength()
        assert conrod2.MassPerLength() == 21.0, conrod2.MassPerLength()

    def test_rod_mass_01(self):
        """tests a CROD and a CONROD mass"""
        eid = 10
        pid = 67
        nid1 = 2
        nid2 = 3
        mid = 5
        A = 4.0
        J = 42.
        c = 3.333
        nsm = 0.
        size = 8

        E = 1e7
        G = 2e7
        nu = 0.325
        rho = 0.1
        L = 1.5
        xyz1 = [1., 0., 0.]
        xyz2 = [2.5, 0., 0.]
        self.get_mass(nid1, nid2, xyz1, xyz2, eid, pid, mid, A, J, c, nsm, E, G, nu, rho, L)
        nsm = 1.
        self.get_mass(nid1, nid2, xyz1, xyz2, eid, pid, mid, A, J, c, nsm, E, G, nu, rho, L)

    def get_mass(self, nid1, nid2, xyz1, xyz2, eid, pid, mid, A, J, c, nsm, E, G, nu, rho, L):
        """tests a CROD and a CONROD"""
        model = BDF(debug=False)
        lines = ['conrod,%i, %i, %i, %i, %f, %f, %f, %f' % (eid, nid1, nid2, mid, A, J, c, nsm)]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        conrod = CONROD.add_card(cardi)
        model._add_element_object(conrod)
        card = model.elements[eid]
        node_ids = card.node_ids
        assert node_ids == [nid1, nid2], node_ids # probably wrong

        lines = ['crod,%i, %i, %i, %i' % (eid+1, pid, nid1, nid2)]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        crod = CROD.add_card(cardi)
        model._add_element_object(crod)
        card = model.elements[eid+1]
        node_ids = card.node_ids
        assert node_ids == [nid1, nid2], node_ids # probably wrong

        lines = ['ctube,%i, %i, %i, %i' % (eid+2, pid+1, nid1, nid2)]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        ctube = CTUBE.add_card(cardi)
        model._add_element_object(ctube)
        card = model.elements[eid+2]
        node_ids = card.node_ids
        assert node_ids == [nid1, nid2], node_ids # probably wrong

        lines = ['prod,%i, %i, %f, %f, %f, %f' % (pid, mid, A, J, c, nsm)]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        prod = PROD.add_card(cardi)
        model._add_property_object(prod)

        #self.pid = integer(card, 1, 'pid')
        #self.mid = integer(card, 2, 'mid')
        #self.OD1 = double(card, 3, 'OD1')
        #self.t = double_or_blank(card, 4, 't', self.OD1 / 2.)
        #self.nsm = double_or_blank(card, 5, 'nsm', 0.0)
        #self.OD2 = double_or_blank(card, 6, 'OD2', self.OD1)
        OD1 = sqrt(4*A/pi)
        t = 0.
        OD2 = OD1
        lines = ['ptube,%i, %i, %f, %f, %f, %f' % (pid+1, mid, OD1, t, nsm, OD2)]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        ptube = PTUBE.add_card(cardi)
        model._add_property_object(ptube)

        lines = ['mat1,%i, %.2e, %.2e, %f, %f' % (mid, E, G, nu, rho)]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        card = MAT1.add_card(cardi)
        model._add_structural_material_object(card)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid1, 0, xyz1[0], xyz1[1], xyz1[2])]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        card = GRID.add_card(cardi)
        model._add_node_object(card)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid2, 0, xyz2[0], xyz2[1], xyz2[2])]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        card = GRID.add_card(cardi)
        model._add_node_object(card)

        model.cross_reference()
        mass = L * (rho * A + nsm)

        # conrod
        self.assertEqual(conrod.eid, eid)
        self.assertEqual(conrod.Pid(), -10)
        self.assertEqual(conrod.Mid(), mid)
        self.assertEqual(conrod.Length(), L)
        self.assertEqual(conrod.Nsm(), nsm)
        self.assertEqual(conrod.Mass(), mass)
        self.assertEqual(conrod.E(), E)
        self.assertEqual(conrod.G(), G)
        self.assertEqual(conrod.Area(), A)
        self.assertEqual(conrod.J(), J)
        self.assertEqual(conrod.C(), c)
        self.assertEqual(conrod.Rho(), rho)

        # crod
        self.assertEqual(crod.eid, eid+1)
        self.assertEqual(crod.Pid(), pid)
        self.assertEqual(crod.Mid(), mid)
        self.assertEqual(crod.Length(), L)
        self.assertEqual(crod.Nsm(), nsm)
        self.assertEqual(crod.Mass(), mass)
        self.assertEqual(crod.E(), E)
        self.assertEqual(crod.G(), G)
        self.assertEqual(crod.Area(), A)
        self.assertEqual(crod.J(), J)
        self.assertEqual(crod.C(), c)
        self.assertEqual(crod.Rho(), rho)
        #self.assertEqual(crod.Nu(), nu)

        # prod
        self.assertEqual(prod.Pid(), pid)
        self.assertEqual(prod.Mid(), mid)
        self.assertEqual(prod.Nsm(), nsm)
        self.assertEqual(prod.E(), E)
        self.assertEqual(prod.G(), G)
        self.assertEqual(prod.Area(), A)
        self.assertEqual(prod.J(), J)
        self.assertEqual(prod.C(), c)
        self.assertEqual(prod.Rho(), rho)

        # ctube
        self.assertEqual(ctube.eid, eid+2)
        self.assertEqual(ctube.Pid(), pid+1)
        self.assertEqual(ctube.Mid(), mid)
        self.assertEqual(ctube.Length(), L)
        self.assertEqual(ctube.Nsm(), nsm)
        self.assertAlmostEqual(ctube.Mass(), mass, 5)
        self.assertEqual(ctube.E(), E)
        self.assertEqual(ctube.G(), G)
        self.assertAlmostEqual(ctube.Area(), A, 5)
        ctube.J()
        self.assertEqual(ctube.Rho(), rho)

        # ptube
        self.assertEqual(ptube.Pid(), pid+1)
        self.assertEqual(ptube.Mid(), mid)
        self.assertEqual(ptube.Nsm(), nsm)
        self.assertEqual(ptube.E(), E)
        self.assertEqual(ptube.G(), G)
        self.assertAlmostEqual(ptube.Area(), A, 5)
        ptube.J()
        self.assertEqual(ptube.Rho(), rho)

    def test_rod_opt(self):
        model = BDF(debug=False,)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        eid_rod = 1
        pid_rod = 1
        mid = 1
        A = 1.0
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_crod(eid_rod, pid_rod, [1, 2])
        model.add_prod(pid_rod, mid, A)
        model.add_mat1(mid, E, G, nu)

        eid_tube = 2
        pid_tube = 2
        model.add_ctube(eid_tube, pid_tube, [1, 2])
        model.add_ptube(pid_tube, mid, 1.0)
        params = [
            ('A', 2.0),
            ('J', 2.0),
        ]
        i = make_dvprel_optimization(model, params, 'PROD', pid_rod, i=1)

        params = [
            ('T', 2.0),
        ]
        i = make_dvprel_optimization(model, params, 'PTUBE', pid_tube, i)

        model.cross_reference()
        model.update_model_by_desvars()

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
