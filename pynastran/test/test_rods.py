from six.moves import StringIO
from math import pi, sqrt

import unittest

from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.bdf import CROD, CONROD, PROD, CTUBE, PTUBE, GRID, MAT1

#from pyNastran.bdf.field_writer_8 import print_card_8

bdf = BDF(debug=False)
class TestRods(unittest.TestCase):
    def test_crod_01(self):
        lines = ['CROD          10     100      10       2']
        card = bdf.process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = CROD.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()
        self.assertEqual(card.Eid(), 10)
        self.assertEqual(card.Pid(), 100)
        node_ids = card.node_ids
        assert node_ids == [10, 2], node_ids # probably wrong

    def test_conrod_01(self):
        eid = 10
        nid1 = 2
        nid2 = 3
        mid = 5
        A = 27.0
        lines = ['conrod,%i, %i, %i, %i, %f' % (eid, nid1, nid2, mid, A)]
        card = bdf.process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = CONROD.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()
        self.assertEqual(card.Eid(), eid)
        self.assertEqual(card.Mid(), mid)
        node_ids = card.node_ids
        assert node_ids == [nid1, nid2], node_ids

    def test_mass_01(self):
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
        card = model.process_card(lines)
        cardi = BDFCard(card)
        conrod = CONROD.add_card(cardi)
        model.add_element(conrod)
        card = model.elements[eid]
        node_ids = card.node_ids
        assert node_ids == [nid1, nid2], node_ids # probably wrong

        lines = ['crod,%i, %i, %i, %i' % (eid+1, pid, nid1, nid2)]
        card = model.process_card(lines)
        cardi = BDFCard(card)
        crod = CROD.add_card(cardi)
        model.add_element(crod)
        card = model.elements[eid+1]
        node_ids = card.node_ids
        assert node_ids == [nid1, nid2], node_ids # probably wrong

        lines = ['ctube,%i, %i, %i, %i' % (eid+2, pid+1, nid1, nid2)]
        card = model.process_card(lines)
        cardi = BDFCard(card)
        ctube = CTUBE.add_card(cardi)
        model.add_element(ctube)
        card = model.elements[eid+2]
        node_ids = card.node_ids
        assert node_ids == [nid1, nid2], node_ids # probably wrong

        lines = ['prod,%i, %i, %f, %f, %f, %f' % (pid, mid, A, J, c, nsm)]
        card = model.process_card(lines)
        cardi = BDFCard(card)
        prod = PROD.add_card(cardi)
        model.add_property(prod)

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
        card = model.process_card(lines)
        cardi = BDFCard(card)
        ptube = PTUBE.add_card(cardi)
        model.add_property(ptube)

        lines = ['mat1,%i, %.2e, %.2e, %f, %f' % (mid, E, G, nu, rho)]
        card = model.process_card(lines)
        cardi = BDFCard(card)
        card = MAT1.add_card(cardi)
        model.add_structural_material(card)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid1, 0, xyz1[0], xyz1[1], xyz1[2])]
        card = model.process_card(lines)
        cardi = BDFCard(card)
        card = GRID.add_card(cardi)
        model.add_node(card)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid2, 0, xyz2[0], xyz2[1], xyz2[2])]
        card = model.process_card(lines)
        cardi = BDFCard(card)
        card = GRID.add_card(cardi)
        model.add_node(card)

        model.cross_reference()
        mass = L * (rho * A + nsm)

        # conrod
        self.assertEqual(conrod.Eid(), eid)
        self.assertEqual(conrod.Pid(), None)
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
        self.assertEqual(crod.Eid(), eid+1)
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
        self.assertEqual(ctube.Eid(), eid+2)
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

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
