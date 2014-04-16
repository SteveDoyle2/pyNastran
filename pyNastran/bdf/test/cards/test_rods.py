import StringIO
import unittest

from itertools import izip, count

from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.bdf import CROD, CONROD

from pyNastran.bdf.fieldWriter import print_card

bdf = BDF()
class TestRods(unittest.TestCase):
    def test_crod_01(self):
        lines = ['CROD          10     100      10       2']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CROD(card)
        card.write_bdf(size, 'dummy')
        card.rawFields()
        self.assertEquals(card.Eid(), 10)
        self.assertEquals(card.Pid(), 100)

    def test_conrod_01(self):
        eid = 10
        nid1 = 2
        nid2 = 3
        mid = 5
        A = 27.0
        lines = ['conrod,%i, %i, %i, %i, %f' % (eid, nid1, nid2, mid, A)]
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CONROD(card)
        card.write_bdf(size, 'dummy')
        card.rawFields()
        self.assertEquals(card.Eid(), eid)
        self.assertEquals(card.Mid(), mid)


if __name__ == '__main__':
    unittest.main()
