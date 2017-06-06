from __future__ import print_function
import os
import copy
import unittest

from pyNastran.bdf.bdf import read_bdf, BDF, BDFCard
from pyNastran.bdf.cards.elements.solid import (
    #CTETRA4, CHEXA8, CPENTA6,
    #CTETRA10, CHEXA20,
    CPENTA15
)


class TestSolids(unittest.TestCase):

    def test_cpenta_01(self):
        """tests a cpenta15"""
        lines = [
            'CPENTA,85,22,201,202,203,205,206,207,+PN2',
            '+PN2,209,210,217,  ,  ,  ,213,214,',
            ',218'
        ]
        bdf = BDF(debug=False)
        card = bdf.process_card(lines)
        card = BDFCard(card)

        solid = CPENTA15.add_card(card, comment='cpenta15')
        solid.write_card(size=8, is_double=False)
        solid.write_card(size=16, is_double=False)
        solid.raw_fields()

        node_ids = solid.node_ids
        assert node_ids == [201, 202, 203, 205, 206, 207,
                            209, 210, 217, None, None, None, 213, 214, 218], node_ids
        nids = [201, 202, 203, 205, 206, 207,
                209, 210, 217, None, None, None, 213, 214, 218]
        CPENTA15.add_card(card, comment='spike')
        eid = 85
        pid = 22
        bdf.add_cpenta(eid, pid, nids, comment='spike')

    def test_cpenta_01b(self):
        pass
          # ..todo:: this is wrong...
        #lines = [  # this will fail!
        #    'CPENTA,85,22,201,202,203,205,206,207,+PN2',
        #    '+PN2,209,210,217,  ,  ,  ,213,214,218'
        #]

    #def test_solid_02(self):
        #"""checks nonlinear static solid material"""
        #pass

    def test_solid_01(self):
        """checks linear static solid material"""
        mid = 2
        pid = 4
        rho = 0.1
        cards = [
            #$ Solid Nodes
            ['GRID', 11, 0, 0., 0., 0., 0],
            ['GRID', 12, 0, 1., 0., 0., 0],
            ['GRID', 13, 0, 1., 1., 0., 0],
            ['GRID', 14, 0, 0., 1., 0., 0],

            ['GRID', 15, 0, 0., 0., 2., 0],
            ['GRID', 16, 0, 1., 0., 2., 0],
            ['GRID', 17, 0, 1., 1., 2., 0],
            ['GRID', 18, 0, 0., 1., 2., 0],

            # Solids
            ['CHEXA', 7, pid, 11, 12, 13, 14, 15, 16, 17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID', 21, 0, 0., 0., 0., 0,],
            ['GRID', 22, 0, 1., 0., 0., 0,],
            ['GRID', 23, 0, 1., 1., 0., 0,],
            ['GRID', 24, 0, 0., 0., 2., 0,],
            ['GRID', 25, 0, 1., 0., 2., 0,],
            ['GRID', 26, 0, 1., 1., 2., 0,],
            ['CPENTA', 9, pid, 21, 22, 23, 24, 25, 26],

            # static
            ['PSOLID', pid, mid, 0],
            ['MAT1', mid, 1.0, 2.0, 3.0, rho]
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        # CTETRA
        eid = 8
        mid = 2
        pid = 4
        nsm = 0.
        V = 1. / 3.
        rho = 0.1
        self.check_solid(model, eid, 'CTETRA', pid, 'PSOLID', mid, 'MAT1', nsm, rho, V)

        eid = 9
        V = 1.0
        self.check_solid(model, eid, 'CPENTA', pid, 'PSOLID', mid, 'MAT1', nsm, rho, V)

        eid = 7
        V = 2.0
        self.check_solid(model, eid, 'CHEXA', pid, 'PSOLID', mid, 'MAT1', nsm, rho, V)

    def test_solid_02(self):
        mid = 2
        pid = 4
        rho = 0.1
        cards = [
            #$ Solid Nodes
            ['GRID', 11, 0, 0., 0., 0., 0],
            ['GRID', 12, 0, 1., 0., 0., 0],
            ['GRID', 13, 0, 1., 1., 0., 0],
            ['GRID', 14, 0, 0., 1., 0., 0],

            ['GRID', 15, 0, 0., 0., 2., 0],
            ['GRID', 16, 0, 1., 0., 2., 0],
            ['GRID', 17, 0, 1., 1., 2., 0],
            ['GRID', 18, 0, 0., 1., 2., 0],

            # Solids
            ['CHEXA', 7, pid, 11, 12, 13, 14, 15, 16, 17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID', 21, 0, 0., 0., 0., 0,],
            ['GRID', 22, 0, 1., 0., 0., 0,],
            ['GRID', 23, 0, 1., 1., 0., 0,],
            ['GRID', 24, 0, 0., 0., 2., 0,],
            ['GRID', 25, 0, 1., 0., 2., 0,],
            ['GRID', 26, 0, 1., 1., 2., 0,],
            ['CPENTA', 9, pid, 21, 22, 23, 24, 25, 26],

            # hyperelastic
            ['PLSOLID', pid, mid, 'GRID'],
            ['MATHP', mid, None, None, None, rho],
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        # CTETRA
        eid = 8
        nsm = 0.
        V = 1. / 3.
        self.check_solid(model, eid, 'CTETRA', pid, 'PLSOLID', mid, 'MATHP', nsm, rho, V)

        eid = 9
        V = 1.0
        self.check_solid(model, eid, 'CPENTA', pid, 'PLSOLID', mid, 'MATHP', nsm, rho, V)

        eid = 7
        V = 2.0
        self.check_solid(model, eid, 'CHEXA', pid, 'PLSOLID', mid, 'MATHP', nsm, rho, V)

    def test_solid_03(self):
        """checks linear static solid material"""
        mid = 2
        pid = 4
        rho = 0.1
        cards = [
            #$ Solid Nodes
            ['GRID', 11, 0, 0., 0., 0., 0],
            ['GRID', 12, 0, 1., 0., 0., 0],
            ['GRID', 13, 0, 1., 1., 0., 0],
            ['GRID', 14, 0, 0., 1., 0., 0],

            ['GRID', 15, 0, 0., 0., 2., 0],
            ['GRID', 16, 0, 1., 0., 2., 0],
            ['GRID', 17, 0, 1., 1., 2., 0],
            ['GRID', 18, 0, 0., 1., 2., 0],

            # Solids
            ['CHEXA', 7, pid, 11, 12, 13, 14, 15, 16, 17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID', 21, 0, 0., 0., 0., 0,],
            ['GRID', 22, 0, 1., 0., 0., 0,],
            ['GRID', 23, 0, 1., 1., 0., 0,],
            ['GRID', 24, 0, 0., 0., 2., 0,],
            ['GRID', 25, 0, 1., 0., 2., 0,],
            ['GRID', 26, 0, 1., 1., 2., 0,],
            ['CPENTA', 9, pid, 21, 22, 23, 24, 25, 26],

            # static
            ['PSOLID', pid, mid, 0],
            ['MAT1', mid, 1.0, 2.0, 3.0, rho],
            ['MATS1', mid, None, 'PLASTIC', 0.0, 1, 1, 100000., ],
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

    def test_solid_04(self):
        """checks linear static solid material"""
        mid = 2
        pid = 4
        rho = 0.1
        table_id = 42
        cards = [
            #$ Solid Nodes
            ['GRID', 11, 0, 0., 0., 0., 0],
            ['GRID', 12, 0, 1., 0., 0., 0],
            ['GRID', 13, 0, 1., 1., 0., 0],
            ['GRID', 14, 0, 0., 1., 0., 0],

            ['GRID', 15, 0, 0., 0., 2., 0],
            ['GRID', 16, 0, 1., 0., 2., 0],
            ['GRID', 17, 0, 1., 1., 2., 0],
            ['GRID', 18, 0, 0., 1., 2., 0],

            # Solids
            ['CHEXA', 7, pid, 11, 12, 13, 14, 15, 16, 17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID', 21, 0, 0., 0., 0., 0,],
            ['GRID', 22, 0, 1., 0., 0., 0,],
            ['GRID', 23, 0, 1., 1., 0., 0,],
            ['GRID', 24, 0, 0., 0., 2., 0,],
            ['GRID', 25, 0, 1., 0., 2., 0,],
            ['GRID', 26, 0, 1., 1., 2., 0,],
            ['CPENTA', 9, pid, 21, 22, 23, 24, 25, 26],

            # static
            ['PSOLID', pid, mid, 0],
            ['MAT1', mid, 1.0, 2.0, 3.0, rho],
            ['MATS1', mid, table_id, 'PLASTIC', 0.0, 1, 1, 100000., ],
            #['TABLEST'],
            ['TABLES1', table_id, 1, None, None, None, None, None, None,
             1.0, 10.0, 2.0, 10.0, 'ENDT'],
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        mat = model.Material(mid)
        mat.E()

    def test_solids_ctetra4(self):
        """tests a CTETRA4"""
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model = BDF(debug=False)
        model.add_grid(11, xyz=[0., 0., 0.])
        model.add_grid(12, xyz=[1., 0., 0.])
        model.add_grid(13, xyz=[1., 1., 0.])
        model.add_grid(15, xyz=[0., 0., 2.])
        model.add_psolid(pid, mid)
        mat1 = model.add_mat1(mid, E, G, nu, rho=0.1)
        nids = [11, 12, 13, 15]
        ctetra = model.add_ctetra(eid, pid, nids, comment='ctetra')
        end_checks(model)

    def test_solids_ctetra4_mat9(self):
        """tests a CTETRA4"""
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model = BDF(debug=False)
        model.add_grid(11, xyz=[0., 0., 0.])
        model.add_grid(12, xyz=[1., 0., 0.])
        model.add_grid(13, xyz=[1., 1., 0.])
        model.add_grid(15, xyz=[0., 0., 2.])
        model.add_psolid(pid, mid)
        mat9 = model.add_mat9(mid, G11=1000., rho=0.2)
        mat9.raw_fields()

        nids = [11, 12, 13, 15]
        model.add_ctetra(eid, pid, nids, comment='ctetra')
        end_checks(model)

    def test_solids_ctetra10(self):
        """tests a CTETRA10"""
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model = BDF(debug=False)
        g110 = model.add_grid(110, xyz=[0., 0., 0.])
        g120 = model.add_grid(120, xyz=[1., 0., 0.])
        g130 = model.add_grid(130, xyz=[1., 1., 0.])
        g140 = model.add_grid(140, xyz=[0., 0., 2.])

        model.add_grid(111, xyz=g110.xyz+g120.xyz)
        model.add_grid(112, xyz=g120.xyz+g130.xyz)
        model.add_grid(113, xyz=g130.xyz+g110.xyz)

        model.add_grid(121, xyz=g110.xyz+g140.xyz)
        model.add_grid(122, xyz=g120.xyz+g140.xyz)
        model.add_grid(123, xyz=g130.xyz+g140.xyz)

        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [
            110, 120, 130, 140,
            111, 112, 113,
            121, 122, 123
        ]
        model.add_ctetra(eid, pid, nids, comment='ctetra10')
        end_checks(model)

    def test_solids_cpyram5(self):
        """tests a CPYRAM5"""
        model = BDF(debug=False)
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model.add_grid(10, xyz=[0., 0., 0.])
        model.add_grid(20, xyz=[1., 0., 0.])
        model.add_grid(30, xyz=[1., 1., 0.])
        model.add_grid(40, xyz=[0., 0., 2.])
        model.add_grid(50, xyz=[1., 1., 2.])
        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [10, 20, 30, 40, 50]
        model.add_cpyram(eid, pid, nids, comment='cpyram')

        nids2 = copy.deepcopy(nids)
        nids2.append(None)
        eid += 1
        elem2 = model.add_cpyram(eid, pid, nids2, comment='cpyram13')
        elem2.write_card(size=8)
        elem2.write_card(size=16)
        elem2.write_card_16(is_double=False)

        end_checks(model)

    def test_solids_cpenta(self):
        """tests a CPENTA6"""
        model = BDF(debug=False)
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model.add_grid(21, xyz=[0., 0., 0.])
        model.add_grid(22, xyz=[1., 0., 0.])
        model.add_grid(23, xyz=[1., 1., 0.])
        model.add_grid(24, xyz=[0., 0., 2.])
        model.add_grid(25, xyz=[1., 0., 2.])
        model.add_grid(26, xyz=[1., 1., 2.])
        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [21, 22, 23, 24, 25, 26]
        model.add_cpenta(eid, pid, nids, comment='cpenta')

        nids2 = copy.deepcopy(nids)
        nids2.append(None)
        eid += 1
        elem2 = model.add_cpenta(eid, pid, nids2, comment='cpenta15')
        elem2.write_card(size=8)
        elem2.write_card(size=16)
        elem2.write_card_16(is_double=False)
        end_checks(model)

    def test_solids_chexa(self):
        """tests a CHEXA8"""
        model = BDF(debug=False)
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model.add_grid(11, xyz=[0., 0., 0.])
        model.add_grid(12, xyz=[1., 0., 0.])
        model.add_grid(13, xyz=[1., 1., 0.])
        model.add_grid(14, xyz=[0., 1., 0.])

        model.add_grid(15, xyz=[0., 0., 2.])
        model.add_grid(16, xyz=[1., 0., 2.])
        model.add_grid(17, xyz=[1., 1., 2.])
        model.add_grid(18, xyz=[0., 1., 2.])
        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [11, 12, 13, 14, 15, 16, 17, 18]
        elem = model.add_chexa(eid, pid, nids, comment='chexa')
        elem.write_card(size=8)
        elem.write_card(size=16)
        elem.write_card_16(is_double=False)

        nids2 = copy.deepcopy(nids)
        nids2.append(None)
        eid += 1
        elem2 = model.add_chexa(eid, pid, nids2, comment='chexa16')
        elem2.write_card(size=8)
        elem2.write_card(size=16)
        elem2.write_card_16(is_double=False)

        end_checks(model)
        elem = model.elements[eid]
        assert elem.Mass() > 0, elem.Mass()

    def check_solid(self, model, eid, etype, pid, ptype, mid, mtype, nsm, rho, V):
        """checks that various solid methods work"""
        mass = rho * V
        element = model.elements[eid]
        element.node_ids
        assert pid in model.properties, 'pid is missing for\n%s' % str(element)
        self.assertEqual(element.type, etype)
        self.assertEqual(element.eid, eid)
        self.assertEqual(element.pid.type, ptype)
        self.assertEqual(element.Pid(), pid)
        self.assertEqual(element.pid.mid.type, mtype)
        self.assertEqual(element.Mid(), mid)
        self.assertEqual(element.Volume(), V)
        self.assertEqual(element.Mass(), mass)


def end_checks(model):
    """various checks"""
    model.validate()
    model._verify_bdf(xref=False)
    model.cross_reference()
    model._verify_bdf(xref=True)
    model.uncross_reference()
    model.cross_reference()
    model.pop_xref_errors()
    mass, cg, inertia = model.mass_properties()
    assert mass > 0, 'mass=%s, cg=%s, inertia=%s' % (mass, cg, inertia)

    bdf_filename = 'solid_test.bdf'
    model.write_bdf(bdf_filename)
    model2 = read_bdf(bdf_filename, debug=False)
    os.remove(bdf_filename)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
