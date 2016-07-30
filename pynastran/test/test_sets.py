import unittest

from pyNastran.bdf.field_writer_8 import print_int_card_blocks


class TestSets(unittest.TestCase):
    def test_set3_01(self):
        fields_blocks = [
            'SET1',
            [['a', 1.0, 3], False], # these are not all integers
            [[1, 2, 3], True], # these are all integers
        ]
        msg = print_int_card_blocks(fields_blocks)
        self.assertEqual('SET1           a      1.       3       1       2       3\n', msg)

        fields_blocks = [
            'SET1',
            [['a', 1.0, 3], False], # these are not all integers
            [[1, 2, 3, 5, 4], True], # these are all integers
        ]
        msg2 = print_int_card_blocks(fields_blocks)
        #print('%r' % msg2)
        self.assertEqual('SET1           a      1.       3       1       2       3       5       4\n', msg2)

        fields_blocks = [
            'SET1',
            [['a', 1.0, 3], False], # these are not all integers
            [[1, 2, 3, 5, 4, 6], True], # these are all integers
        ]
        msg3 = print_int_card_blocks(fields_blocks)
        #print('%r' % msg3)
        self.assertEqual('SET1           a      1.       3       1       2       3       5       4\n'
                         '               6\n', msg3)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
