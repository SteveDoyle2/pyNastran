"""no special requirements"""
import unittest
import numpy as np

from pyNastran.gui.menus.groups_modify.groups import Group, NodeGroup, _get_collapsed_text
from pyNastran.gui.menus.groups_modify.utils import get_groups_sorted_by_name

class TestGroupsModify(unittest.TestCase):
    def test_groups_1(self):
        nodes_dict = {'node_str': '', 'nodes_pound': '100',}
        data = {
            'font_size' : 8,
            0 : Group(name='main', element_str='1:#', elements_pound='103', **nodes_dict, editable=False),
            1 : Group(name='wing', element_str='1:10', elements_pound='103', **nodes_dict, editable=True),
            2 : Group(name='fuselage1', element_str=[3, 33], elements_pound='103', **nodes_dict, editable=True),
            3 : Group(name='fuselage2', element_str='50:60', elements_pound='103', **nodes_dict, editable=True),
            4 : Group(name='fuselage3', element_str='50:60', elements_pound='103', node_str=[1,2,3], nodes_pound='103', editable=True),
        }
        assert np.array_equal(data[1].element_ids, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        assert np.array_equal(data[4].node_ids, [1, 2, 3])
        str(data[4])

        data[3].element_ids = np.array([24, 29, 28], dtype='int32')
        assert np.array_equal(data[3].element_ids, [24, 28, 29])
        data[3].node_ids = np.array([4, 9, 8], dtype='int32')
        assert np.array_equal(data[3].node_ids, [4, 8, 9])

        groups = get_groups_sorted_by_name(data)
        assert groups == ['main', 'wing',
                          'fuselage1', 'fuselage2', 'fuselage3',], groups

    def test_groups_2(self):
        data = {
            'font_size': 8,
            0: NodeGroup(name='main', node_str='1:#', nodes_pound='103', editable=False),
            1: NodeGroup(name='wing', node_str='1:10', nodes_pound='103', editable=True),
            2: NodeGroup(name='fuselage1', node_str='50:60', nodes_pound='103', editable=True),
            3: NodeGroup(name='fuselage2', node_str='50:60', nodes_pound='103', editable=True),
            4: NodeGroup(name='fuselage3', node_str='50:60', nodes_pound='103', editable=True),
            5: NodeGroup(name='fuselage13', node_str=[1,2,3], nodes_pound='103', editable=True),
        }
        assert np.array_equal(data[1].node_ids, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        assert np.array_equal(data[5].node_ids, [1, 2, 3])
        data[5].node_ids = np.array([4, 9, 8], dtype='int32')
        assert np.array_equal(data[5].node_ids, [4, 8, 9])
        print(data[5])
        str(data[4])

        groups = get_groups_sorted_by_name(data)
        assert groups == ['main', 'wing',
                          'fuselage1', 'fuselage2', 'fuselage3', 'fuselage13'], groups


if __name__ == '__main__':
    unittest.main()
