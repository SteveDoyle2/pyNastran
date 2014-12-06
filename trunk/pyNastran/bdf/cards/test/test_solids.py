import unittest

from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.bdf import CPENTA15

bdf = BDF(debug=False)

class TestSolids(unittest.TestCase):

    def test_cpenta_01(self):
        lines = ['CPENTA,85,22,201,202,203,205,206,207,+PN2',
                 '+PN2,209,210,217,  ,  ,  ,213,214,218']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CPENTA15(card)
        card.write_bdf(size, 'dummy')
        card.raw_fields()

    def test_solid_02(self):
        """checks nonlinear static solid material"""
        pass

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
            ['CHEXA', 7, pid, 11, 12, 13, 14, 15, 16,  17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID', 21, 0, 0., 0., 0.,  0,],
            ['GRID', 22, 0, 1., 0., 0.,  0,],
            ['GRID', 23, 0, 1., 1., 0.,  0,],
            ['GRID', 24, 0, 0., 0., 2.,  0,],
            ['GRID', 25, 0, 1., 0., 2.,  0,],
            ['GRID', 26, 0, 1., 1., 2.,  0,],
            ['CPENTA', 9, pid, 21, 22, 23, 24, 25, 26],

            # static
            ['PSOLID', pid, mid, 0],
            ['MAT1',  mid, 1.0, 2.0, 3.0, rho]
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

    def check_solid(self, model, eid, etype, pid, ptype, mid, mtype, nsm, rho, V):
        mass = rho * V
        element = model.elements[eid]
        self.assertEquals(element.type, etype)
        self.assertEquals(element.Eid(), eid)
        self.assertEquals(element.pid.type, ptype)
        self.assertEquals(element.Pid(), pid)
        self.assertEquals(element.pid.mid.type, mtype)
        self.assertEquals(element.Mid(), mid)
        self.assertEquals(element.Volume(), V)
        self.assertEquals(element.Mass(), mass)

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
            ['CHEXA', 7, pid, 11, 12, 13, 14, 15, 16,  17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID', 21, 0, 0., 0., 0.,  0,],
            ['GRID', 22, 0, 1., 0., 0.,  0,],
            ['GRID', 23, 0, 1., 1., 0.,  0,],
            ['GRID', 24, 0, 0., 0., 2.,  0,],
            ['GRID', 25, 0, 1., 0., 2.,  0,],
            ['GRID', 26, 0, 1., 1., 2.,  0,],
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
            ['CHEXA',  7, pid, 11, 12, 13, 14, 15, 16,  17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID',  21, 0, 0., 0., 0.,  0,],
            ['GRID',  22, 0, 1., 0., 0.,  0,],
            ['GRID',  23, 0, 1., 1., 0.,  0,],
            ['GRID',  24, 0, 0., 0., 2.,  0,],
            ['GRID',  25, 0, 1., 0., 2.,  0,],
            ['GRID',  26, 0, 1., 1., 2.,  0,],
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

    def test_solid_03(self):
        """checks linear static solid material"""
        mid = 2
        pid = 4
        rho = 0.1
        tableID = 42
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
            ['CHEXA',  7, pid, 11, 12, 13, 14, 15, 16,  17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID',  21, 0, 0., 0., 0.,  0,],
            ['GRID',  22, 0, 1., 0., 0.,  0,],
            ['GRID',  23, 0, 1., 1., 0.,  0,],
            ['GRID',  24, 0, 0., 0., 2.,  0,],
            ['GRID',  25, 0, 1., 0., 2.,  0,],
            ['GRID',  26, 0, 1., 1., 2.,  0,],
            ['CPENTA', 9, pid, 21, 22, 23, 24, 25, 26],

            # static
            ['PSOLID', pid, mid, 0],
            ['MAT1', mid, 1.0, 2.0, 3.0, rho],
            ['MATS1', mid, tableID, 'PLASTIC', 0.0, 1, 1, 100000., ],
            #['TABLEST'],
            ['TABLES1', tableID, 1, None, None, None, None, None, None,
            1.0, 10.0, 2.0, 10.0, 'ENDT'],
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        mat = model.Material(mid)
        mat.E()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
