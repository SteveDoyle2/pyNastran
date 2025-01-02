#from io import StringIO
from math import pi, sqrt
import unittest

import numpy as np

from pyNastran.dev.bdf_vectorized3.bdf import BDF, BDFCard
from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_shells import make_dvprel_optimization
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck

#from pyNastran.bdf.field_writer_8 import print_card_8

class TestRods(unittest.TestCase):
    def test_crod_01(self):
        model = BDF(debug=False)
        lines = ['CROD          10     100      10       2']
        card = model._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        crod_id = model.crod.add_card(cardi, ifile=0)
        A = 1.0
        pid = 100
        mid = 1
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_prod(pid, mid, A)
        model.add_mat1(mid, E, G, nu)

        crod = model.crod
        crod.write(size, 'dummy')
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.setup()
        #crod.raw_fields()
        self.assertEqual(crod.element_id, 10)
        self.assertEqual(crod.property_id, 100)
        node_ids = crod.nodes[0]
        assert np.array_equal(node_ids, [10, 2]), node_ids # probably wrong
        assert np.array_equal(crod.get_edge_ids(), [(2, 10)])
        save_load_deck(model)

    def test_prod_nsm(self):
        model = BDF(debug=False)
        pid = 1
        mid = 1
        nsm = 1.
        A = 2.0
        prod_id = model.add_prod(pid, mid, A, j=0., c=0., nsm=nsm, comment='')
        #print(prod)

        E = 1.0
        G = None
        nu = 0.3
        mat1_id = model.add_mat1(mid, E, G, nu)
        mat1 = model.mat1
        #----------------
        card_lines = [
            'PROD           2       1      2.                      1.',
        ]
        model.add_card(card_lines, 'PROD', comment='', is_list=False,
                       has_none=True)
        model.pop_parse_errors()
        model.cross_reference()
        prod2 = model.Property(2)[0]
        #------------------
        prod = model.prod.slice_card_by_id(1)
        #assert prod.nsm() == 1.0
        assert prod.area() == 2.0

        # mass/L = area*rho + nsm
        assert prod.mass_per_length() == 1.0

        # area = 2.0
        mat1.rho = np.array([10.0])
        assert prod.mass_per_length() == 21.0, prod.mass_per_length()
        assert prod2.mass_per_length() == 21.0, prod2.mass_per_length()
        save_load_deck(model)

    def test_ptube_nsm(self):
        """tests a PTUBE for the NSM field"""
        model = BDF(debug=False)
        pid = 1
        mid = 1
        nsm = 1.
        OD1 = 2.0
        ptube_id = model.add_ptube(pid, mid, OD1, t=None, nsm=nsm, OD2=None,
                                   comment='ptube')
        #ptube = model.ptube
        #print(ptube)

        E = 1.0
        G = None
        nu = 0.3
        rho = 0.0
        mat1_id = model.add_mat1(mid, E, G, nu, rho=rho)
        mat1 = model.mat1
        #----------------
        card_lines = [
            'PTUBE          2       1      2.              1.',
        ]
        model.add_card(card_lines, 'PTUBE', comment='', is_list=False,
                       has_none=True)
        model.pop_parse_errors()
        model.setup(run_geom_check=True)
        ptube = model.Property(1)[0]
        ptube2 = model.Property(2)[0]
        #print(ptube.get_stats())
        #------------------
        model.cross_reference()

        #assert ptube.Nsm() == 1.0

        area = pi
        mpl = area * rho + nsm
        assert ptube.area() == area, ptube.area()

        # mass/L = area*rho + nsm
        assert ptube.mass_per_length() == 1.0, ptube.mass_per_length()

        rho = 10.0
        mat1.rho = np.array([rho])
        mpl = area * rho + nsm
        assert ptube2.mass_per_length() == mpl, ptube2.MassPerLength()

        eid = 100
        nids = [10, 12]
        ctube_id = model.add_ctube(eid, pid, nids, comment='ctube')
        ctube = model.ctube
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(12, [2., 0., 0.])

        #model.safe_cross_reference()
        #model.uncross_reference()
        model.cross_reference()
        assert np.array_equal(ctube.get_edge_ids(), [(10, 12)])
        #mass, cg, inertia = mass_properties(
            #model, element_ids=None, mass_ids=None,
            #reference_point=None, sym_axis=None,
            #scale=None, inertia_reference='cg')
        mass, cg, inertia = model.inertia_sum(element_id=None)

        # L * (rho * a + nsm)
        A = pi * (OD1**2) / 4.
        mass_expected = 2.0 * (rho * A + nsm)
        assert np.allclose(mass, mass_expected), mass
        model.setup()
        save_load_deck(model)

    def test_conrod_01(self):
        model = BDF(debug=False)
        eid = 10
        nid1 = 2
        nid2 = 3
        mid = 5
        A = 27.0
        lines = ['conrod,%i, %i, %i, %i, %f' % (eid, nid1, nid2, mid, A)]
        card = model._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        conrod_id = model.conrod.add_card(cardi, ifile=0)
        model.setup(run_geom_check=False)

        conrod = model.conrod
        conrod.write(size, 'dummy')
        #conrod.raw_fields()
        self.assertEqual(conrod.element_id, eid)
        self.assertEqual(conrod.material_id, mid)
        model.setup(run_geom_check=False)
        node_ids = conrod.nodes[0]
        assert np.array_equal(node_ids, [nid1, nid2]), node_ids
        assert np.array_equal(conrod.get_edge_ids(), [(2, 3)])

    def test_conrod_nsm(self):
        model = BDF(debug=False)
        mid = 1
        nsm = 1.
        A = 2.0
        eid = 1
        model.add_grid(1, [0., 0., 0.], comment='node1')
        model.add_grid(2, [1., 0., 0.], comment='node2')
        nids = [1, 2]
        conrod_id = model.add_conrod(eid, mid, nids, A=A, j=0.0, c=0.0, nsm=nsm,
                                     comment='')
        model.setup(run_geom_check=False)

        grid = model.grid.slice_card_by_id(2)
        conrod = model.conrod.slice_card_by_element_id(eid)
        assert np.array_equal(conrod.get_edge_ids(), [(1, 2)])

        E = 1.0
        G = None
        nu = 0.3
        mat1_id = model.add_mat1(mid, E, G, nu)
        mat1 = model.mat1
        #----------------
        card_lines = [
            'CONROD         2       1       2       1      2.                      1.',
        ]
        model.add_card(card_lines, 'CONROD', comment='', is_list=False,
                       has_none=True)
        model.pop_parse_errors()
        model.setup()
        conrod2 = model.Element(2)[0]
        assert np.array_equal(conrod2.get_edge_ids(), [(1, 2)])

        #------------------
        model.cross_reference()

        assert conrod.nonstructural_mass() == 1.0
        assert conrod.area() == 2.0

        # mass/L = area*rho + nsm
        assert conrod.mass_per_length() == 1.0

        # area = 2.0
        mat1.rho = np.array([10.0])
        assert conrod.mass_per_length() == 21.0, conrod.mass_per_length()
        assert conrod2.mass_per_length() == 21.0, conrod2.mass_per_length()

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
        conrod_id = model.conrod.add_card(cardi, ifile=0)
        conrod = model.conrod
        #add_methods = model._add_methods
        #add_methods._add_element_object(conrod)
        model.setup(run_geom_check=False)
        card = model.Element(eid)[0]
        node_ids = card.nodes[0]
        assert np.array_equal(node_ids, [nid1, nid2]), node_ids # probably wrong
        assert np.array_equal(conrod.get_edge_ids(), [(nid1, nid2)])

        lines = ['crod,%i, %i, %i, %i' % (eid+1, pid, nid1, nid2)]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        crod_id = model.crod.add_card(cardi, ifile=0)

        model.setup(run_geom_check=False)
        crod = model.crod
        #add_methods._add_element_object(crod)
        card = model.Element(eid+1)[0]
        node_ids = card.nodes[0]
        assert np.array_equal(node_ids, [nid1, nid2]), node_ids # probably wrong
        assert np.array_equal(crod.get_edge_ids(), [(nid1, nid2)])

        lines = ['ctube,%i, %i, %i, %i' % (eid+2, pid+1, nid1, nid2)]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        ctube_id = model.ctube.add_card(cardi, ifile=0)
        ctube = model.ctube

        model.setup(run_geom_check=False)
        card = model.Element(eid+2)[0]
        node_ids = card.nodes[0]
        assert np.array_equal(node_ids, [nid1, nid2]), node_ids # probably wrong
        assert np.array_equal(ctube.get_edge_ids(), [(nid1, nid2)])

        lines = ['prod,%i, %i, %f, %f, %f, %f' % (pid, mid, A, J, c, nsm)]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        prod_i = model.prod.add_card(cardi, ifile=0)
        prod = model.prod
        ptube = model.ptube

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
        ptube_id = ptube.add_card(cardi, ifile=0)
        #add_methods._add_property_object(ptube)

        lines = ['mat1,%i, %.2e, %.2e, %f, %f' % (mid, E, G, nu, rho)]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        card = model.mat1.add_card(cardi, ifile=0)
        mat1 = model.mat1
        #add_methods._add_structural_material_object(card)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid1, 0, xyz1[0], xyz1[1], xyz1[2])]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        card = model.grid.add_card(cardi, ifile=0)
        #add_methods._add_node_object(card)

        lines = ['grid,%i, %i, %f, %f, %f' % (nid2, 0, xyz2[0], xyz2[1], xyz2[2])]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        card = model.grid.add_card(cardi, ifile=0)
        #add_methods._add_node_object(card)

        model.cross_reference()
        mass = L * (rho * A + nsm)

        # conrod
        self.assertEqual(conrod.element_id, eid)
        self.assertEqual(conrod.material_id, mid)
        #self.assertEqual(conrod.Pid(), -10)
        self.assertEqual(conrod.length(), L)
        self.assertEqual(conrod.nonstructural_mass(), nsm)
        self.assertEqual(conrod.mass(), mass)
        #self.assertEqual(conrod.E(), E)
        #self.assertEqual(conrod.G(), G)
        self.assertEqual(conrod.area(), A)
        self.assertEqual(conrod.J, J)
        self.assertEqual(conrod.c, c)
        #self.assertEqual(conrod.Rho(), rho)

        # crod
        self.assertEqual(crod.element_id, eid+1)
        self.assertEqual(crod.property_id, pid)
        #self.assertEqual(crod.material_id(), mid)
        self.assertEqual(crod.length(), L)
        self.assertEqual(crod.nonstructural_mass(), nsm)
        self.assertEqual(crod.mass(), mass)
        #self.assertEqual(crod.E(), E)
        #self.assertEqual(crod.G(), G)
        self.assertEqual(crod.area(), A)
        #self.assertEqual(crod.J(), J)
        #self.assertEqual(crod.C(), c)
        #self.assertEqual(crod.Rho(), rho)
        #self.assertEqual(crod.Nu(), nu)

        # prod
        self.assertEqual(prod.property_id, pid)
        self.assertEqual(prod.material_id, mid)
        self.assertEqual(prod.nsm, nsm)
        #self.assertEqual(prod.E(), E)
        #self.assertEqual(prod.G(), G)
        self.assertEqual(prod.area(), A)
        #self.assertEqual(prod.J(), J)
        self.assertEqual(prod.c, c)
        #self.assertEqual(prod.Rho(), rho)

        # ctube
        self.assertEqual(ctube.element_id, eid+2)
        self.assertEqual(ctube.property_id, pid+1)
        #self.assertEqual(ctube.material_id, mid)
        self.assertEqual(ctube.length(), L)
        self.assertEqual(ctube.nonstructural_mass(), nsm)
        assert np.allclose(ctube.mass(), mass)
        #self.assertEqual(ctube.E(), E)
        #self.assertEqual(ctube.G(), G)
        assert np.allclose(ctube.area(), A)
        #ctube.J()
        #self.assertEqual(ctube.Rho(), rho)

        # ptube
        self.assertEqual(ptube.property_id, pid+1)
        self.assertEqual(ptube.material_id, mid)
        self.assertEqual(ptube.nsm, nsm)
        #self.assertEqual(ptube.E(), E)
        #self.assertEqual(ptube.G(), G)
        assert np.allclose(ptube.area(), A)
        #ptube.J()
        #self.assertEqual(ptube.Rho(), rho)

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
        if 0:
            model.update_model_by_desvars()

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
