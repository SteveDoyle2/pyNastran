import StringIO
import unittest

from itertools import izip, count

from pyNastran.bdf.bdf import BDF, BDFCard, PBAR #, GRID, MAT1

from pyNastran.bdf.fieldWriter import print_card

bdf = BDF(debug=False)
class TestBars(unittest.TestCase):
    def test_pbar_01(self):
        fields = [u'PBAR', 1510998, 1520998, 0.0, 4.9000000000000006e-14, 4.9000000000000006e-14, 0.0, 0.0, None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None, None, 0.0]
        card = print_card(fields)
        #print(card)
        card = print_card(fields)
        lines = card.split('\n')
        card = bdf.process_card(lines)
        card2 = BDFCard(card)
        with self.assertRaises(AssertionError):  # A=0, I12=0, K1=0
            pbar = PBAR(card2)

    def test_pbar_02(self):
        pid = 1
        mid = 2
        A = None
        I1 = I2 = None
        J = None
        nsm = None
        c1=c2=d1=d2=e1=e2=f1=f2=None
        k1=k2=None
        i12 = 3.
        fields = ['PBAR', pid, mid, A, I1, I2, J, nsm, None,
                          c1, c2, d1, d2, e1, e2, f1, f2,
                          k1, k2, i12]
        card = print_card(fields)
        lines = card.split('\n')
        card = bdf.process_card(lines)
        #print(card)
        card2 = BDFCard(card)
        #print(card2)
        pbar = PBAR(card2)
        self.assertEqual(pbar.pid, 1)
        self.assertEqual(pbar.mid, 2)
        self.assertEqual(pbar.A, 0.0)
        self.assertEqual(pbar.i1, 0.0)
        self.assertEqual(pbar.i2, 0.0)
        self.assertEqual(pbar.j, 0.0)
        self.assertEqual(pbar.nsm, 0.0)
        self.assertEqual(pbar.i12, 3.0)
        self.assertEqual(pbar.C1, 0.0)
        self.assertEqual(pbar.C2, 0.0)
        self.assertEqual(pbar.D1, 0.0)
        self.assertEqual(pbar.D2, 0.0)
        self.assertEqual(pbar.E1, 0.0)
        self.assertEqual(pbar.E2, 0.0)
        self.assertEqual(pbar.K1, 1e8)
        self.assertEqual(pbar.K2, 1e8)


if __name__ == '__main__':
    unittest.main()
