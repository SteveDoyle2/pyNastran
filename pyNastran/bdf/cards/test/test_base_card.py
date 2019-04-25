from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import unittest
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.cards.base_card import expand_thru, expand_thru_by #, expand_thru_exclude

class TestBaseCard(unittest.TestCase):
    """Tests methods used by ``BaseCard``"""
    def test_expand_thru(self):
        """tests expand_thru"""
        values = expand_thru(['1', 'THRU', '10'])
        assert values == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values

        values = expand_thru(['1', 'thru', '10'])
        assert values == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values

        values = expand_thru([1, 2, 3, 4, 5])
        assert values == [1, 2, 3, 4, 5], values

        values = expand_thru_by(['1', 'THRU', '10', 'BY', 2])
        assert values == [1, 3, 5, 7, 9, 10], values

        values = expand_thru_by(['1', 'thru', '10', 'by', 2])
        assert values == [1, 3, 5, 7, 9, 10], values

        values = expand_thru_by([1, 2, 3, 4, 5])
        assert values == [1, 2, 3, 4, 5], values

        #values = expand_thru_exclude(['1', 'thru', '5', 'exclude', 2])
        #assert values == [1, 3, 4, 5], values

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


if __name__ == '__main__':
    unittest.main()
