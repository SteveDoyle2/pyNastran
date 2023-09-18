from math import pi, sqrt

import unittest

from pyNastran.dev.bdf_vectorized.bdf import BDF #, BDFCard, CROD, CONROD
#from pyNastran.dev.bdf_vectorized.cards.elements.rod.rods import (
#    CROD, CONROD, PROD)
#CTUBE, PTUBE, GRID, MAT1

#from pyNastran.bdf.fieldWriter import print_card

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
        card = model.crod.slice_by_element_id([10])
        card.write_card(f, size)
        #card.raw_fields()
        self.assertEqual(card.get_element_id_by_element_index(), 10)
        self.assertEqual(card.get_property_id_by_element_index(), 100)

    def test_conrod_01(self):
        model = BDF(debug=debug)
        model.allocate({'CONROD': 1})
        eid = 10
        nid1 = 2
        nid2 = 3
        mid = 5
        A = 27.0
        lines = [f'conrod,{eid:d}, {nid1:d}, {nid2:d}, {mid:d}, {A:f}']
        model.add_card(lines, 'CONROD', is_list=False)

        size = 8
        card = model.conrod.slice_by_element_id([eid])
        f = StringIO()
        card.write_card(f, size)
        #card.raw_fields()
        self.assertEqual(card.get_element_id_by_element_index(), eid)
        self.assertEqual(card.get_material_id_by_element_index(), mid)

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

        lines = ['ctube,%i, %i, %i, %i' % (eid+2, pid+1, nid1, nid2)]
        model.add_card(lines, 'ctube', is_list=False)

        lines = ['prod,%i, %i, %f, %f, %f, %f' % (pid, mid, A, J, c, nsm)]
        model.add_card(lines, 'prod', is_list=False)

        OD1 = sqrt(4 * A / pi)
        t = 0.
        OD2 = OD1
        lines = ['ptube,%i, %i, %f, %f, %f, %f' % (pid+1, mid, OD1, t, nsm, OD2)]
        model.add_card(lines, 'ptube', is_list=False)

        lines = ['mat1,%i, %.2e, %.2e, %f, %f' % (mid, E, G, nu, rho)]
        model.add_card(lines, 'mat1', is_list=False)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid1, 0, xyz1[0], xyz1[1], xyz1[2])]
        model.add_card(lines, 'grid', is_list=False)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid2, 0, xyz2[0], xyz2[1], xyz2[2])]
        model.add_card(lines, 'grid', is_list=False)

        model.build()
        mass = L * (rho * A + nsm)

        f = StringIO()
        model.write_bdf(out_filename=f, interspersed=True, size=8,
                       precision='single',
                       enddata=None)
        print(f.getvalue())
        #positions = model.get_positions()
        grid_cid0 = None

        # conrod
        conrod = model.conrod.slice_by_element_id(eid)
        self.assertEqual(conrod.get_element_id_by_element_index(), eid)
        #self.assertEqual(conrod.get_property_id_by_element_id(), None)
        self.assertEqual(conrod.get_material_id_by_element_id(eid), mid)
        self.assertEqual(conrod.get_length_by_element_index(i=None, grid_cid0=grid_cid0), L)
        #self.assertEqual(conrod.Nsm(), nsm)


        rhoi = conrod.get_density_by_element_id(eid)
        Ai = conrod.get_area_by_element_id(eid)
        Li = conrod.get_length_by_element_id(eid, grid_cid0=grid_cid0)
        nsmi = conrod.get_non_structural_mass_by_element_id(eid)
        massa = conrod.get_mass_by_element_index()
        mass_msg_conrod = 'mass = L * (rho * A + nsm)\n'
        mass_msg_conrod += 'L=%s expected=%s\n' % (Li, L)
        mass_msg_conrod += 'rho=%s expected=%s\n' % (rhoi, rho)
        mass_msg_conrod += 'A=%s expected=%s\n' % (Ai, A)
        mass_msg_conrod += 'nsm=%s expected=%s\n' % (nsmi, nsm)
        mass_msg_conrod += 'mass=%s actual=%s expected=%s\n' % (Li * (rhoi*Ai + nsmi), massa, mass)
        #mass_msg_conrod += 'mass=%s expected=%s\n' % (Li * (rhoi*Ai + nsmi), mass)

        self.assertEqual(massa, mass, mass_msg_conrod)
        #self.assertEqual(conrod.E(), E)
        #self.assertEqual(conrod.G(), G)
        #self.assertEqual(conrod.area(), A)
        #self.assertEqual(conrod.J(), J)
        #self.assertEqual(conrod.C(), c)
        #self.assertEqual(conrod.Rho(), rho)

        # crod
        crod_eid = eid + 1
        crod = model.crod.slice_by_element_id([crod_eid])
        self.assertEqual(crod.get_element_id_by_element_index(), crod_eid)
        self.assertEqual(crod.get_property_id_by_element_id(crod_eid), pid)
        self.assertEqual(crod.get_material_id_by_element_id(crod_eid), mid)
        rhoi = crod.get_density_by_element_id(crod_eid)
        Ai = crod.get_area_by_element_id(crod_eid)
        Li = crod.get_length_by_element_id(crod_eid, grid_cid0=grid_cid0)
        nsmi = crod.get_non_structural_mass_by_element_id(crod_eid)
        self.assertEqual(Li, L)
        #self.assertEqual(crod.Nsm(), nsm)


        massa = crod.get_mass_by_element_id(crod_eid)
        mass_msg_crod = 'mass = L * (rho * A + nsm)\n'
        mass_msg_crod += 'L=%s expected=%s\n' % (Li, L)
        mass_msg_crod += 'rho=%s expected=%s\n' % (rhoi, rho)
        mass_msg_crod += 'A=%s expected=%s\n' % (Ai, A)
        mass_msg_crod += 'nsm=%s expected=%s\n' % (nsmi, nsm)
        mass_msg_crod += 'mass=%s actual=%s expected=%s\n' % (Li * (rhoi*Ai + nsmi), massa, mass)
        self.assertEqual(massa, mass, mass_msg_crod)
        #self.assertEqual(crod.E(), E)
        #self.assertEqual(crod.G(), G)
        #self.assertEqual(crod.area(), A)
        #self.assertEqual(crod.J(), J)
        #self.assertEqual(crod.C(), c)
        #self.assertEqual(crod.Rho(), rho)
        #self.assertEqual(crod.Nu(), nu)

        # prod
        prod = model.prod.slice_by_property_id([pid])
        self.assertEqual(prod.property_id[0], pid)
        self.assertEqual(prod.get_material_id_by_property_id(pid), mid)
        self.assertEqual(prod.get_non_structural_mass_by_property_id(pid), nsm)
        self.assertEqual(prod.get_E_by_property_id(pid), E)
        self.assertEqual(prod.get_G_by_property_id(pid), G)
        self.assertEqual(prod.get_area_by_property_id(pid), A)
        self.assertEqual(prod.get_J_by_property_id(pid), J)
        self.assertEqual(prod.get_c_by_property_id(pid), c)
        self.assertEqual(prod.get_density_by_property_id(pid), rho)

        # ctube
        if 1:
            ctube_eid = eid + 2
            ptube_pid = pid + 1
            assert ctube_eid == 12, ctube_eid
            assert ptube_pid == 68, ptube_pid
            ctube = model.ctube.slice_by_element_id(ctube_eid)
            self.assertEqual(ctube.get_element_id_by_element_index(), ctube_eid)
            self.assertEqual(ctube.get_property_id_by_element_id(ctube_eid), ptube_pid)
            self.assertEqual(ctube.get_material_id_by_element_id(ctube_eid), mid)
            self.assertEqual(ctube.get_length_by_element_id(ctube_eid, grid_cid0), L)
            self.assertEqual(ctube.get_non_structural_mass_by_element_id(ctube_eid), nsm)
            self.assertAlmostEquals(ctube.get_mass_by_element_id(ctube_eid), mass, 5)
            self.assertEqual(ctube.get_E_by_element_id(ctube_eid), E)
            self.assertEqual(ctube.get_G_by_element_id(ctube_eid), G)
            self.assertAlmostEquals(ctube.get_area_by_element_id(ctube_eid), A, 5)
            ctube.get_J_by_element_id(ctube_eid)
            self.assertEqual(ctube.get_density_by_element_id(), rho)

            # ptube
            ptube = model.ptube.slice_by_property_id(pid+1)
            self.assertEqual(ptube.get_property_id_by_property_index(), pid+1)
            self.assertEqual(ptube.get_material_id_by_property_id(), mid)
            self.assertEqual(ptube.get_non_structural_mass_by_property_id(), nsm)
            self.assertEqual(ptube.get_E_by_property_id(), E)
            self.assertEqual(ptube.get_G_by_property_id(), G)
            self.assertAlmostEquals(ptube.get_area_by_property_id(), A, 5)
            ptube.get_J_by_property_id()
            self.assertEqual(ptube.get_density_by_property_id(), rho)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
