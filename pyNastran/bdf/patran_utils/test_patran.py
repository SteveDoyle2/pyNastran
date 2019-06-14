import os
import unittest
import numpy as np

from pyNastran.bdf.patran_utils.colon_syntax import (
    parse_patran_syntax, parse_patran_syntax_dict, parse_patran_syntax_dict_map,
    write_patran_syntax_dict)
from pyNastran.bdf.patran_utils.read_patran_custom_results import read_patran
import pyNastran

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models', 'patran_fmt')


class TestPatranSyntax(unittest.TestCase):
    def test_bdf_utils_01(self):
        """tests parse_patran_syntax"""
        msg = '1:10  14:20:2  50:40:-1'
        output = parse_patran_syntax(msg, pound=None)
        expected = np.array(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
             14, 16, 18, 20,
             40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]
        )
        error_msg = 'expected equal; A-B=%s; B-A=%s' % (
            np.setdiff1d(output, expected), np.setdiff1d(expected, output))
        assert np.array_equal(output, expected), error_msg

        msg = '1:#'
        output = parse_patran_syntax(msg, pound=5)
        assert np.array_equal(output, [1, 2, 3, 4, 5])

        msg = '#:1'
        with self.assertRaises(ValueError):
            output = parse_patran_syntax(msg, pound=None)
        #assert array_equal(output, [1, 2, 3, 4, 5])

        msg = '1:#'
        output = parse_patran_syntax(msg, pound='5')
        assert np.array_equal(output, [1, 2, 3, 4, 5])

        # should this raise an error?
        msg = '#:1'
        #with self.assertRaises(ValueError):
        output = parse_patran_syntax(msg, pound='5')

    def test_bdf_utils_02(self):
        """tests parse_patran_syntax_dict"""
        msg = 'n 1:10  14:20:2  50:40:-1 e 10 20'
        expected_nodes = np.array(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
             14, 16, 18, 20,
             40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]
        )
        expected_elements = np.array([10, 20])

        output_dict = parse_patran_syntax_dict(msg, pound_dict=None)
        assert np.array_equal(output_dict['n'], expected_nodes)
        assert np.array_equal(output_dict['e'], expected_elements)


        # messing with the order a bit
        msg = 'n 1:#  e 2:#:4 junk 1 7 12:#:2 3'
        expected_nodes = np.array([1, 2, 3, 4, 5])
        expected_elements = np.array([2, 6, 10, 14])
        expected_junk = np.array([1, 3, 7, 12, 14, 16, 18, 20])
        pound_dict = {
            'n' : 5,
            'e' : 14,
            'junk' : 20.,
        }

        output_dict = parse_patran_syntax_dict(msg, pound_dict=pound_dict)
        assert np.array_equal(output_dict['n'], expected_nodes)
        assert np.array_equal(output_dict['e'], expected_elements)

        error_msg = 'expected equal; A-B=%s; B-A=%s' % (
            np.setdiff1d(expected_junk, output_dict['junk']),
            np.setdiff1d(output_dict['junk'], expected_junk))
        assert np.array_equal(output_dict['junk'], expected_junk), error_msg

        msg = write_patran_syntax_dict({'e' : [2, 6, 10, 14]})
        assert msg == 'e 2:14:4', 'msg=%r' % msg

        msg = write_patran_syntax_dict({'e' : [1, 2, 6, 10, 14]})
        assert msg == 'e 1 2 6:14:4', 'msg=%r' % msg

        msg = write_patran_syntax_dict(
            {
                'n' : [1, 2, 6, 10, 14],
                'e' : [1, 2, 6, 10, 14],
            },
        )
        assert msg == 'e 1 2 6:14:4 n 1 2 6:14:4', 'msg=%r' % msg

        out = parse_patran_syntax_dict('')
        assert len(out) == 0, 'out=%s' % out

    def test_bdf_utils_03(self):
        """tests parse_patran_syntax_dict"""
        node_sets = "e 1:3 n 2:6:2 Node 10:13 N 15 coord 1:10"
        type_map = {
            'n' : 'Node',
            'Node' : 'Node',
            'e' : 'Element',
            'Elm' : 'Element',
            'Element' : 'Element',
        }

        data = parse_patran_syntax_dict(node_sets, type_map)
        data_expected = {
            'Element' : np.array([1, 2, 3]),
            'Node' : np.array([2, 4, 6, 10, 11, 12, 13, 15]),
        }

        data = parse_patran_syntax_dict_map(node_sets, type_map, msg='')
        assert len(data.keys()) == len(data_expected.keys()), 'data.keys=%s data_expected.keys=%s' % (data.keys(), data_expected.keys())
        for key, value in sorted(data.items()):
            assert key in data_expected, 'cant find key=%r' % key
            value_expected = data_expected[key]
            assert np.array_equal(value, value_expected), 'key=%r\nvalue=%r\nexpected=%r' % (key, value, value_expected)


class TestPatran(unittest.TestCase):
    def test_read_patran(self):
        """tests the nodal patran read function"""
        #bdf_filename = os.path.join(MODEL_PATH, '0012_20.bdf')
        node_filename = os.path.join(MODEL_PATH, 'normals.nod')
        data = read_patran(node_filename, fdtype='float64', idtype='int32')
        nids = data['nids']
        assert nids.min() == 1, nids
        assert nids.max() == 2214, nids


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
