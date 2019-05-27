import os
import unittest

import pyNastran
from pyNastran.bdf.bdf import BDF, GRID

pkg_path = pyNastran.__path__[0]
mesh_utils_path = os.path.join(pkg_path, 'bdf', 'mesh_utils', 'test')



class TestOpenMDAO(unittest.TestCase):

    def test_openmdao_good_1(self):
        """we replace WTMASS, YES_A, YES_B, YES_C with values"""
        updates = [
            #['MAT1', 3, 10.0],  # 3 is E -> set to 10.0
            #['MAT1', 4, 10.0],  # 3 is G -> set to 10.0
            ['GRID', 1, 3, 10.0],  # 3 is x1 -> set to 10.0
            ['GRID', 1, 4, 20.0],  # 4 is x2 -> set to 20.0
            ['CPENTA', 9, 2, 10],  # 2 is property_id -> set to 10
            ['CPENTA', 9, 3, 20],  # 3 is node1 -> set to 20
            ['PSOLID', 4, 1, 2],   # 1 is material_id
            ['PARAM', 'WTMASS', 1, 'WTMASs'],  # key
            ['PARAM', 'WTMASS', 2, 0.0025],  # value1
            ['PCOMP', 1, 2, 1.],
            ['PCOMP', 1, 3, 2.],
            ['CTETRA', 8, 3, 1], # nid[0]
            ['CTETRA', 8, 4, 2], # nid[1]
            ['CTETRA', 8, 5, 3], # nid[2]
            ['CTETRA', 8, 6, 4], # nid[3]
        ]
        #GRID           1       0      0.      0.      0.       0
        #GRID           2       0      1.      0.      0.       0
        #GRID           3       0      1.      1.      0.       0
        #GRID           4       0      0.      1.      0.       0
        #CPENTA         9       4      21      22      23      24      25      26
        #PSOLID   4       1       0
        #CTETRA         8       4      11      12      13      15

        bdf_filename = os.path.join(mesh_utils_path, 'test_mass.dat')

        model = BDF(debug=False)
        model.read_bdf(bdf_filename)
        pcomp_updates = [
            ['PCOMP', 1, 15, 'YES_A', 'souts_0'],
            ['PCOMP', 1, 19, 'YES_B', 'souts_1'],

            ['PCOMP', 1, 25, 'YES_C', 'souts_2'],
            #['PCOMP', 1, 29, 'YES_D', 'souts_3'],
        ]
        for iupdate in updates:
            card_type, itype, ifield, value = iupdate
            card = model.update_card(card_type, itype, ifield, value)

        for iupdate in pcomp_updates:
            card_type, itype, ifield, value, field_name = iupdate
            card = model.update_card(card_type, itype, ifield, value)
            if '_' in field_name:
                field_name2, index = field_name.split('_')
                index = int(index)
                actual = getattr(card, field_name2)[index]
                assert actual == value, 'field_name=%r ifield=%s value=%s actual=%s\n%s' % (
                    field_name, ifield, value, actual, card.print_raw_card())
            #if card_type == 'PCOMP':
                #print(card)

    def test_openmdao_bad_1(self):
        updates = [  # KeyError
            ['GRID', 1, 0, 20.0],
            ['GRID', 1, 10, 20.0],
            ['CHEXA', 100, 3, 20],
            ['FAKECARD', 100, 3, 20],
        ]
        bdf_filename = os.path.join(mesh_utils_path, 'test_mass.dat')

        model = BDF(debug=False)
        model.read_bdf(bdf_filename)

        for iupdate in updates:
            card_type, itype, ifield, value = iupdate
            self.assertRaises(KeyError, model.update_card, card_type, itype, ifield, value)
            #print("tried to apply %r=%s" % (card_type, itype))

    def test_openmdao_bad_2(self):
        updates = [  # IndexError
            ['PARAM', 'WTMASS', 3, 0.005],  # value2; invalid b/c WTMASS
            ['PCOMP', 1, 26, 'MID_D'],  # too many plies
            ['PCOMP', 1, 27, 'T_D'],  # too many plies
            ['PCOMP', 1, 28, 'THETA_D'],  # too many plies
            ['PCOMP', 1, 29, 'YES_D'],  # too many plies
        ]
        bdf_filename = os.path.join(mesh_utils_path, 'test_mass.dat')

        model = BDF(debug=False)
        model.read_bdf(bdf_filename)

        for iupdate in updates:
            card_type, itype, ifield, value = iupdate
            self.assertRaises(IndexError, model.update_card, card_type, itype, ifield, value)
            #print("tried to apply %r=%s" % (card_type, itype))

    def test_openmaod_bad_3(self):
        params_bad = {1 : '10'}
        params_good = {'cat' : '10'}
        model = BDF(debug=False)

        with self.assertRaises(TypeError):
            model.set_dynamic_syntax(params_bad)

        model.set_dynamic_syntax(params_good)
        val = model._parse_dynamic_syntax('%cat')
        self.assertEqual('10', val)
        #with self.assertRaises(KeyError):
            #self._parse_dynamic_syntax(key)

    def test_openmdao_field_1(self):
        model = BDF(debug=False)
        node = GRID(57, cp=10, xyz=[1., 2., 3.], cd=7, ps='', seid=0, comment='')
        model.nodes[57] = node
        #print()
        str(node)
        model.update_card('GRID', 57, 2, 11) # cp
        str(node)
        model.update_card('GRID', 57, 2, 12) # cp
        str(node)
        model.update_card('GRID', 57, 3, 13.) # x
        str(node)
        model.update_card('GRID', 57, 4, 14.) # y
        str(node)
        model.update_card('GRID', 57, 4, 14.) # y
        str(node)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
