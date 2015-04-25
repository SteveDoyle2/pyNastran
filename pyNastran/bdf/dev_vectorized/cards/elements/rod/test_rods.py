from six.moves import StringIO
from math import pi, sqrt
from itertools import count

import unittest

from pyNastran.bdf.dev_vectorized.bdf import BDF, BDFCard, CROD, CONROD
#from pyNastran.bdf.dev_vectorized.cards.elements.rod.rods import (
#    CROD, CONROD, PROD)
#CTUBE, PTUBE, GRID, MAT1

from pyNastran.bdf.fieldWriter import print_card

debug = True
class TestRods(unittest.TestCase):
    def test_crod_01(self):
        model = BDF(debug=debug)
        model.allocate({'CROD': 1})

        lines = ['CROD          10     100      10       2']
        model.add_card(lines, 'CROD', is_list=False)

        size = 8
        f = StringIO()
        #card = CROD(card)
        card = model.crod[10]
        card.write_bdf(f, size)
        #card.rawFields()
        self.assertEquals(card.get_element_id_by_element_index(), 10)
        self.assertEquals(card.get_property_id_by_element_index(), 100)

    def test_conrod_01(self):
        model = BDF(debug=debug)
        model.allocate({'CONROD': 1})
        eid = 10
        nid1 = 2
        nid2 = 3
        mid = 5
        A = 27.0
        lines = ['conrod,%i, %i, %i, %i, %f' % (eid, nid1, nid2, mid, A)]
        model.add_card(lines, 'CONROD', is_list=False)

        size = 8
        card = model.conrod[eid]
        f = StringIO()
        card.write_bdf(f, size)
        #card.rawFields()
        self.assertEquals(card.get_element_id_by_element_index(), eid)
        self.assertEquals(card.get_material_id_by_element_index(), mid)

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
        card_count = {
            'CONROD' : 1,
            'CTUBE' : 1,
            'PTUBE' : 1,
            'CROD' : 1,
            'PROD' : 1,
            'GRID' : 2,
            'MAT1' : 1,
        }
        model = BDF(debug=debug)
        model.allocate(card_count)
        lines = ['conrod,%i, %i, %i, %i, %f, %f, %f, %f' % (eid, nid1, nid2, mid, A, J, c, nsm)]
        model.add_card(lines, 'conrod', is_list=False)

        lines = ['crod,%i, %i, %i, %i' % (eid+1, pid, nid1, nid2)]
        model.add_card(lines, 'crod', is_list=False)

        #lines = ['ctube,%i, %i, %i, %i' % (eid+2, pid+1, nid1, nid2)]
        #model.add_card(lines, 'ctube', is_list=False)

        lines = ['prod,%i, %i, %f, %f, %f, %f' % (pid, mid, A, J, c, nsm)]
        model.add_card(lines, 'prod', is_list=False)

        OD1 = sqrt(4*A/pi)
        t = 0.
        OD2 = OD1
        #lines = ['ptube,%i, %i, %f, %f, %f, %f' % (pid+1, mid, OD1, t, nsm, OD2)]
        #model.add_card(lines, 'ptube', is_list=False)

        lines = ['mat1,%i, %.2e, %.2e, %f, %f' % (mid, E, G, nu, rho)]
        model.add_card(lines, 'mat1', is_list=False)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid1, 0, xyz1[0], xyz1[1], xyz1[2])]
        model.add_card(lines, 'grid', is_list=False)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid2, 0, xyz2[0], xyz2[1], xyz2[2])]
        model.add_card(lines, 'grid', is_list=False)

        model.cross_reference()
        mass = L * (rho * A + nsm)

        # conrod
        conrod = model.conrod[eid]
        self.assertEquals(conrod.get_element_id_by_element_index(), eid)
        #self.assertEquals(conrod.get_property_id_by_element_index(), None)
        self.assertEquals(conrod.get_material_id_by_element_index(), mid)
        self.assertEquals(conrod.get_length_by_element_index(), L)
        #self.assertEquals(conrod.Nsm(), nsm)
        self.assertEquals(conrod.get_mass_by_element_index(), mass)
        #self.assertEquals(conrod.E(), E)
        #self.assertEquals(conrod.G(), G)
        #self.assertEquals(conrod.Area(), A)
        #self.assertEquals(conrod.J(), J)
        #self.assertEquals(conrod.C(), c)
        #self.assertEquals(conrod.Rho(), rho)

        # crod
        crod = model.crod[eid+1]
        self.assertEquals(crod.get_element_id_by_element_index(), eid+1)
        self.assertEquals(crod.get_property_id_by_element_index(), pid)
        self.assertEquals(crod.get_material_id_by_element_index(), mid)
        self.assertEquals(crod.get_length_by_element_index(), L)
        #self.assertEquals(crod.Nsm(), nsm)
        self.assertEquals(crod.get_mass_by_element_index(), mass)
        #self.assertEquals(crod.E(), E)
        #self.assertEquals(crod.G(), G)
        #self.assertEquals(crod.Area(), A)
        #self.assertEquals(crod.J(), J)
        #self.assertEquals(crod.C(), c)
        #self.assertEquals(crod.Rho(), rho)
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
        if 0:
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
