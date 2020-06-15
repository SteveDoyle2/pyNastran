import unittest
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.subcase_utils import expand_thru_case_control
from pyNastran.bdf.cards.expand_card import expand_thru, expand_thru_by


class TestBaseCard(unittest.TestCase):
    """Tests methods used by ``BaseCard``"""
    def test_expand_floats(self):
        """tests expand"""
        values1 = expand_thru_case_control(['500. 2000.'])
        values2 = expand_thru_case_control([500., 2000.])
        assert values1 == [500., 2000.], values1
        assert values2 == [500., 2000.], values2

    def test_expand_thru(self):
        """tests expand_thru"""
        values1 = expand_thru_case_control(['include 14,15, include 18 thru 24'])
        assert values1 == [14, 15, 18, 19, 20, 21, 22, 23, 24], values1

        values1 = expand_thru_case_control(['include 10 thru 17 except 14 thru 15, include 43 thru 45'])
        assert values1 == [10, 11, 12, 13, 16, 17, 43, 44, 45], values1

        values1 = expand_thru_case_control(['1', '2', '21', 'THRU', '25', '31', 'THRU', '37'])
        values2 = expand_thru_case_control(['1', '2', '21', 'THRU', '25', '31', 'THRU', '37 '])
        values3 = expand_thru_case_control([1, 2, 21, 'THRU', 25, 31, 'THRU', 37])
        assert values1 == [1, 2, 21, 22, 23, 24, 25, 31, 32, 33, 34, 35, 36, 37], values1
        assert values2 == [1, 2, 21, 22, 23, 24, 25, 31, 32, 33, 34, 35, 36, 37], values2
        assert values3 == [1, 2, 21, 22, 23, 24, 25, 31, 32, 33, 34, 35, 36, 37], values3

        values1 = expand_thru_case_control(['1 , 9 THRU 12 , 1040'])
        values2 = expand_thru_case_control(['1 , 9 THRU 12 , 1040'])
        assert values1 == [1, 9, 10, 11, 12, 1040], values1
        assert values2 == [1, 9, 10, 11, 12, 1040], values2

        values1 = expand_thru(['1', 'THRU', '10'])
        values2 = expand_thru_case_control(['1', 'THRU', '10'])
        values3 = expand_thru_case_control(['1 THRU 10'])
        values4 = expand_thru_case_control(['1 THRU 10 '])
        assert values1 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values1
        assert values2 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values2
        assert values3 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values3
        assert values4 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values4


        values1 = expand_thru(['1', 'thru', '10'])
        values2 = expand_thru_case_control(['1', 'thru', '10'])
        values3 = expand_thru_case_control(['1 thru 10'])
        values4 = expand_thru_case_control(['1, thru, 10'])
        assert values1 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values1
        assert values2 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values2
        assert values3 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values3
        assert values4 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values4

        values1 = expand_thru([1, 2, 3, 4, 5])
        values2 = expand_thru_case_control([1, 2, 3, 4, 5])
        values3 = expand_thru_case_control(['1, 2, 3, 4, 5'])
        values4 = expand_thru_case_control(['1 2 3 4 5'])
        values5 = expand_thru_case_control(['1 2 3 4 5 '])
        assert values1 == [1, 2, 3, 4, 5], values1
        assert values2 == [1, 2, 3, 4, 5], values2
        assert values3 == [1, 2, 3, 4, 5], values3
        assert values4 == [1, 2, 3, 4, 5], values4
        assert values5 == [1, 2, 3, 4, 5], values5

        values1 = expand_thru_by(['1', 'THRU', '10', 'BY', 2])
        values2 = expand_thru_case_control(['1', 'THRU', '10', 'BY', 2])
        values3 = expand_thru_case_control(['1 THRU 10 BY 2'])
        values4 = expand_thru_case_control(['1 THRU 10 BY 2 '])
        assert values1 == [1, 3, 5, 7, 9, 10], values1
        assert values2 == [1, 3, 5, 7, 9, 10], values2
        assert values3 == [1, 3, 5, 7, 9, 10], values3
        assert values4 == [1, 3, 5, 7, 9, 10], values4

        values1 = expand_thru_by(['1', 'thru', '10', 'by', 2])
        values2 = expand_thru_case_control(['1', 'thru', '10', 'by', 2])
        values3 = expand_thru_case_control(['1 thru 10 by 2'])
        values4 = expand_thru_case_control(['1 thru 10, by, 2,'])
        assert values1 == [1, 3, 5, 7, 9, 10], values1
        assert values2 == [1, 3, 5, 7, 9, 10], values2
        assert values3 == [1, 3, 5, 7, 9, 10], values3
        assert values4 == [1, 3, 5, 7, 9, 10], values4

        values1 = expand_thru_by([1, 2, 3, 4, 5])
        values2 = expand_thru_case_control([1, 2, 3, 4, 5])
        assert values1 == [1, 2, 3, 4, 5], values1
        assert values2 == [1, 2, 3, 4, 5], values2

        values1 = expand_thru_case_control(['1,5,THRU,11,EXCEPT,7,8,13'])
        values2 = expand_thru_case_control([1, 5, 'THRU', 11, 'EXCEPT', 7, 8, 13])
        assert values1 == [1, 5, 6, 9, 10, 11, 13], values1
        assert values2 == [1, 5, 6, 9, 10, 11, 13], values2

    def test_collapse_thru(self):
        """
        tests collapse_thru method used by SETx cards
        """
        data = [1, 2, 3, 4, 5, 10]
        expected = [1, u'THRU', 5, 10]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 17]
        expected = [1, 3, 4, 5, 6, 17]
        msg = 'expected=%s actual=%s' % (expected, collapse_thru_by(data))
        self.assertEqual(collapse_thru_by(data), expected, msg)

        data = [1, 3, 4, 5, 6, 7, 17]
        expected = [1, 3, 4, 'THRU', 7, 17]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 6, 8, 10, 12, 14, 17]
        expected = [1, 3, 4, 'THRU', 14, 'BY', 2, 17]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 101]
        expected = [1, 3, 4, 5, 6, 8, 'THRU', 22, 'BY', 2, 101]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5]
        expected = [1, 'THRU', 5]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [5]
        expected = [5]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5, 7, 9, 11, 12, 14, 16]
        expected = [1, 'THRU', 5,
                    7, 9, 11,
                    12, 14, 16]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2]
        expected = [1, 2]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 5, 7, 9, 11]
        expected = [1, 'THRU', 11, 'BY', 2]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4]
        expected = [1, 'THRU', 4]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3]
        expected = [1, 2, 3]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5, 6, 7, 8]
        expected = [1, 'THRU', 8]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))


if __name__ == '__main__':   # pragma: no cover
    unittest.main()
