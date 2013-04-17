import unittest
from pyNastran.bdf.bdf import BDFCard, SUPORT1

class TestConstraints(unittest.TestCase):
    def test_support1_01(self):
        card = ['SUPORT1', '1']
        card = BDFCard(card)
        with self.assertRaises(AssertionError):  # too short
            SUPORT1(card)
        
        card = ['SUPORT1', '1', '2']
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '0')

        card = ['SUPORT1', '1', '2', '432']
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '234')

        card = ['SUPORT1', '1', '2', '432', 3]
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[1], 3)

        card = ['SUPORT1', '1', '2', None, 3]
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(con.Cs[1], '0')

        card = ['SUPORT1', '1', '2', '432', 3, '1325']
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[1], '1235')

if __name__ == '__main__':
    unittest.main()