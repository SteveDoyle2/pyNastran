"""no special requirements"""
import unittest

from pyNastran.gui.menus.groups_modify.groups import Group, _get_collapsed_text
from pyNastran.gui.menus.groups_modify.utils import get_groups_sorted_by_name

class TestGroupsModify(unittest.TestCase):
    def test_groups_1(self):
        nodes_dict = {'node_str': '', 'nodes_pound': '100',}
        data = {
            'font_size' : 8,
            0 : Group(name='main', element_str='1:#', elements_pound='103', **nodes_dict, editable=False),
            1 : Group(name='wing', element_str='1:10', elements_pound='103', **nodes_dict, editable=True),
            2 : Group(name='fuselage1', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            3 : Group(name='fuselage2', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            4 : Group(name='fuselage3', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            5 : Group(name='fuselage4', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            6 : Group(name='fuselage5', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            7 : Group(name='fuselage6', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            8 : Group(name='fuselage7', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            9 : Group(name='fuselage8', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            10 : Group(name='fuselage9', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            11 : Group(name='fuselage10', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            12 : Group(name='fuselage11', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            13 : Group(name='fuselage12', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            14 : Group(name='fuselage13', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            15 : Group(name='fuselage14', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
        }
        groups = get_groups_sorted_by_name(data)
        assert groups == ['main', 'wing',
                          'fuselage1', 'fuselage2', 'fuselage3', 'fuselage4', 'fuselage5',
                          'fuselage6', 'fuselage7', 'fuselage8', 'fuselage9', 'fuselage10',
                          'fuselage11', 'fuselage12', 'fuselage13', 'fuselage14',], groups

if __name__ == '__main__':
    unittest.main()
