from six.moves import StringIO
import unittest

from pyNastran.bdf.dev_vectorized.bdf import BDF, BDFCard
from pyNastran.bdf.dev_vectorized.cards.elements.solid.elements_solid import (
    ElementsSolid, CTETRA4, CPENTA6, CHEXA8, CTETRA10, CPENTA15, CHEXA20)

class TestSolids(unittest.TestCase):

    def test_cpenta_01(self):
        model = BDF()
        lines = ['CPENTA,85,22,201,202,203,205,206,207']
        card = model.process_card(lines)
        card = BDFCard(card)

        size = 8
        f = StringIO()
        penta = CPENTA6(model)
        penta.allocate(1)
        penta.add(card)
        penta.write_bdf(f, size)
        #card.rawFields()
        print(f.getvalue())

    def test_cpenta_02(self):
        model = BDF()
        lines = ['CPENTA,85,22,201,202,203,205,206,207,+PN2',
                 '+PN2,209,210,217,  ,  ,  ,213,214,218']
        card = model.process_card(lines)
        card = BDFCard(card)

        size = 8
        f = StringIO()
        penta = CPENTA15(model)
        penta.allocate(1)
        penta.add(card)
        penta.write_bdf(f, size)
        #card.rawFields()
        print(f.getvalue())

    def test_chexa_01(self):
        model = BDF()
        lines = ['CHEXA,85,22,201,202,203,205,206,207,+PN2',
                 '+PN2,209,210']
        card = model.process_card(lines)
        card = BDFCard(card)

        size = 8
        f = StringIO()
        hexa = CHEXA8(model)
        hexa.allocate(1)
        hexa.add(card)
        hexa.write_bdf(f, size)
        #card.rawFields()
        print(f.getvalue())

    def test_chexa_02(self):
        model = BDF()
        lines = ['CHEXA,85,22,201,202,203,205,206,207,+PN2',
                 '+PN2,209,210,217,  ,  ,  ,213,214,218']
        card = model.process_card(lines)
        card = BDFCard(card)

        size = 8
        f = StringIO()
        hexa = CHEXA20(model)
        hexa.allocate(1)
        hexa.add(card)
        hexa.write_bdf(f, size)
        #card.rawFields()
        print(f.getvalue())

    def test_ctetra_01(self):
        model = BDF()
        lines = ['CTETRA,85,22,201,202,203,205']
        card = model.process_card(lines)
        card = BDFCard(card)

        size = 8
        f = StringIO()
        hexa = CTETRA4(model)
        hexa.allocate(1)
        hexa.add(card)
        hexa.write_bdf(f, size)
        #card.rawFields()
        print(f.getvalue())

    def test_ctetra_02(self):
        model = BDF()
        lines = ['CTETRA,85,22,201,202,203,205,206,207,+PN2',
                 '+PN2,209,210,217']
        card = model.process_card(lines)
        card = BDFCard(card)

        size = 8
        f = StringIO()
        hexa = CTETRA10(model)
        hexa.allocate(1)
        hexa.add(card)
        hexa.write_bdf(f, size)
        #card.rawFields()
        print(f.getvalue())

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
        card_count = {
            'GRID' : 14,
            'CTETRA4': 1,
            'CPENTA6': 1,
            'CHEXA8' : 1,
            'PSOLID': 1,
            'MAT1': 1,
        }
        model = BDF(debug=True)
        model.allocate(card_count)
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
        self.check_solid(model, eid, 'CTETRA4', pid, 'PSOLID', mid, 'MAT1', nsm, rho, V)

        eid = 9
        V = 1.0
        self.check_solid(model, eid, 'CPENTA6', pid, 'PSOLID', mid, 'MAT1', nsm, rho, V)

        eid = 7
        V = 2.0
        self.check_solid(model, eid, 'CHEXA8', pid, 'PSOLID', mid, 'MAT1', nsm, rho, V)

    def check_solid(self, model, eid, etype, pid, ptype, mid, mtype, nsm, rho, V):
        mass = rho * V
        element = model.elements[eid]
        #assert isinstance(element, ElementsSolid), type(element)
        assert isinstance(element, (CTETRA4, CPENTA6, CHEXA8,
                                    CTETRA10, CPENTA15, CHEXA20)), type(element)
        self.assertEquals(element.type, etype)
        self.assertEquals(element.n, 1)
        self.assertEquals(element.i, 1)
        self.assertEquals(element.get_element_id_by_element_index(), eid)
        #self.assertEquals(element.pid.type, ptype)
        self.assertEquals(element.get_property_id_by_element_index(), pid)
        self.assertEquals(element.get_property_id_by_element_id([eid]), pid)
        #self.assertEquals(element.pid.mid.type, mtype)
        self.assertEquals(element.get_material_id_by_element_id([eid]), mid)
        #self.assertEquals(element.get_material_id_by_element_index(), mid)
        self.assertEquals(element.get_volume_by_element_id([eid]), V)
        #self.assertEquals(element.get_volume_by_element_index(), V)
        self.assertEquals(element.get_mass_by_element_id(), mass)
        #self.assertEquals(element.get_mass_by_element_index(), mass)

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
        card_count = {
            'GRID' : 14,
            'CTETRA4': 1,
            'CPENTA6': 1,
            'CHEXA8' : 1,
            'PLSOLID': 1,
            'MATHP': 1,
        }
        model = BDF(debug=True)
        model.allocate(card_count)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        # CTETRA
        eid = 8
        nsm = 0.
        V = 1. / 3.
        self.check_solid(model, eid, 'CTETRA4', pid, 'PLSOLID', mid, 'MATHP', nsm, rho, V)

        eid = 9
        V = 1.0
        self.check_solid(model, eid, 'CPENTA6', pid, 'PLSOLID', mid, 'MATHP', nsm, rho, V)

        eid = 7
        V = 2.0
        self.check_solid(model, eid, 'CHEXA8', pid, 'PLSOLID', mid, 'MATHP', nsm, rho, V)

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
        card_count = {
            'GRID' : 14,
            'CTETRA4': 1,
            'CPENTA6': 1,
            'CHEXA8' : 1,
            'PSOLID': 1,
            'MAT1': 1,
            'MATS1': 1,
        }
        model = BDF(debug=True)
        model.allocate(card_count)
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
        card_count = {
            'GRID' : 14,
            'CPENTA6': 1,
            'CTETRA4': 1,
            'PSOLID': 1,
            'CHEXA8' : 1,
            'MAT1': 1,
            'MATS1': 1,
            'TABLES1': 1,
        }
        model = BDF(debug=False)
        model.allocate(card_count)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        mat = model.materials[mid]
        print('----MAT----', type(mat))
        print(mat)
        print('E = %s' % mat.get_E_by_material_index())
        print('E = %s' % mat.get_E_by_material_id())



if __name__ == '__main__':  # pragma: no cover
    unittest.main()
