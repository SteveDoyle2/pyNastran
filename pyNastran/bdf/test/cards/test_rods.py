import StringIO
from math import pi, sqrt
from itertools import count

import unittest

from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.bdf import CROD, CONROD, PROD, CTUBE, PTUBE, GRID, MAT1

from pyNastran.bdf.fieldWriter import print_card

bdf = BDF(debug=False)
class TestRods(unittest.TestCase):
    def test_crod_01(self):
        lines = ['CROD          10     100      10       2']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CROD(card)
        card.write_bdf(size, 'dummy')
        card.rawFields()
        self.assertEquals(card.Eid(), 10)
        self.assertEquals(card.Pid(), 100)

    def test_conrod_01(self):
        eid = 10
        nid1 = 2
        nid2 = 3
        mid = 5
        A = 27.0
        lines = ['conrod,%i, %i, %i, %i, %f' % (eid, nid1, nid2, mid, A)]
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CONROD(card)
        card.write_bdf(size, 'dummy')
        card.rawFields()
        self.assertEquals(card.Eid(), eid)
        self.assertEquals(card.Mid(), mid)

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

    def get_mass(self,nid1, nid2, xyz1, xyz2, eid, pid, mid, A, J, c, nsm, E, G, nu, rho, L):
        """tests a CROD and a CONROD"""
        model = BDF(debug=False)
        lines = ['conrod,%i, %i, %i, %i, %f, %f, %f, %f' % (eid, nid1, nid2, mid, A, J, c, nsm)]
        card = model.process_card(lines)
        card = BDFCard(card)
        conrod = CONROD(card)
        model.add_element(conrod)

        lines = ['crod,%i, %i, %i, %i' % (eid+1, pid, nid1, nid2)]
        card = model.process_card(lines)
        card = BDFCard(card)
        crod = CROD(card)
        model.add_element(crod)

        lines = ['ctube,%i, %i, %i, %i' % (eid+2, pid+1, nid1, nid2)]
        card = model.process_card(lines)
        card = BDFCard(card)
        ctube = CTUBE(card)
        model.add_element(ctube)

        lines = ['prod,%i, %i, %f, %f, %f, %f' % (pid, mid, A, J, c, nsm)]
        card = model.process_card(lines)
        card = BDFCard(card)
        prod = PROD(card)
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
        card = BDFCard(card)
        ptube = PTUBE(card)
        model.add_property(ptube)

        lines = ['mat1,%i, %.2e, %.2e, %f, %f' % (mid, E, G, nu, rho)]
        card = model.process_card(lines)
        card = BDFCard(card)
        card = MAT1(card)
        model.add_structural_material(card)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid1, 0, xyz1[0], xyz1[1], xyz1[2])]
        card = model.process_card(lines)
        card = BDFCard(card)
        card = GRID(card)
        model.add_node(card)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid2, 0, xyz2[0], xyz2[1], xyz2[2])]
        card = model.process_card(lines)
        card = BDFCard(card)
        card = GRID(card)
        model.add_node(card)

        model.cross_reference()
        mass = L * (rho * A + nsm)

        # conrod
        self.assertEquals(conrod.Eid(), eid)
        self.assertEquals(conrod.Pid(), None)
        self.assertEquals(conrod.Mid(), mid)
        self.assertEquals(conrod.Length(), L)
        self.assertEquals(conrod.Nsm(), nsm)
        self.assertEquals(conrod.Mass(), mass)
        self.assertEquals(conrod.E(), E)
        self.assertEquals(conrod.G(), G)
        self.assertEquals(conrod.Area(), A)
        self.assertEquals(conrod.J(), J)
        self.assertEquals(conrod.C(), c)
        self.assertEquals(conrod.Rho(), rho)

        # crod
        self.assertEquals(crod.Eid(), eid+1)
        self.assertEquals(crod.Pid(), pid)
        self.assertEquals(crod.Mid(), mid)
        self.assertEquals(crod.Length(), L)
        self.assertEquals(crod.Nsm(), nsm)
        self.assertEquals(crod.Mass(), mass)
        self.assertEquals(crod.E(), E)
        self.assertEquals(crod.G(), G)
        self.assertEquals(crod.Area(), A)
        self.assertEquals(crod.J(), J)
        self.assertEquals(crod.C(), c)
        self.assertEquals(crod.Rho(), rho)
        #self.assertEquals(crod.Nu(), nu)

        # prod
        self.assertEquals(prod.Pid(), pid)
        self.assertEquals(prod.Mid(), mid)
        self.assertEquals(prod.Nsm(), nsm)
        self.assertEquals(prod.E(), E)
        self.assertEquals(prod.G(), G)
        self.assertEquals(prod.Area(), A)
        self.assertEquals(prod.J(), J)
        self.assertEquals(prod.C(), c)
        self.assertEquals(prod.Rho(), rho)

        # ctube
        self.assertEquals(ctube.Eid(), eid+2)
        self.assertEquals(ctube.Pid(), pid+1)
        self.assertEquals(ctube.Mid(), mid)
        self.assertEquals(ctube.Length(), L)
        self.assertEquals(ctube.Nsm(), nsm)
        self.assertAlmostEquals(ctube.Mass(), mass, 5)
        self.assertEquals(ctube.E(), E)
        self.assertEquals(ctube.G(), G)
        self.assertAlmostEquals(ctube.Area(), A, 5)
        ctube.J()
        self.assertEquals(ctube.Rho(), rho)

        # ptube
        self.assertEquals(ptube.Pid(), pid+1)
        self.assertEquals(ptube.Mid(), mid)
        self.assertEquals(ptube.Nsm(), nsm)
        self.assertEquals(ptube.E(), E)
        self.assertEquals(ptube.G(), G)
        self.assertAlmostEquals(ptube.Area(), A, 5)
        ptube.J()
        self.assertEquals(ptube.Rho(), rho)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
