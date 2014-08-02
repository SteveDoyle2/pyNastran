import unittest

from pyNastran.bdf.fieldWriter import print_int_card_blocks


class TestSets(unittest.TestCase):
    def test_set3_01(self):
        fields_blocks = [
            'SET1',
            [['a', 1.0, 3], False], # these are not all integers
            [[1, 2, 3], True], # these are all integers
        ]
        msg = print_int_card_blocks(fields_blocks)
        self.assertTrue('SET1           a      1.       3       1       2       3\n', msg)

if __name__ == '__main__':
    unittest.main()
