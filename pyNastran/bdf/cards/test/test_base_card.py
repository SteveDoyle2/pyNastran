import unittest
from pyNastran.bdf.cards.base_card import expand_thru, expand_thru_by #, expand_thru_exclude

class TestBaseCard(unittest.TestCase):
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

if __name__ == '__main__':
    unittest.main()
