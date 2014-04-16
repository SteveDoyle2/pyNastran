import unittest

from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.bdf import CGAP, PGAP, CDAMP1, CBUSH

bdf = BDF()

class TestElements(unittest.TestCase):


    def test_cbush_01(self):
        lines = ['cbush,101,102,1,,,,,0']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CBUSH(card)
        self.assertEquals(card.Eid(), 101)
        self.assertEquals(card.Pid(), 102)
        card.write_bdf(size, 'dummy')
        card.rawFields()

    def test_cdamp1_01(self):
        lines = ['CDAMP1, 2001, 20, 1001, 1']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CDAMP1(card)
        self.assertEquals(card.Eid(), 2001)
        self.assertEquals(card.Pid(), 20)
        card.write_bdf(size, 'dummy')
        card.rawFields()

    def test_cgap_01(self):
        lines = ['CGAP    899     90      21      99      0.      1.      0.      0']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CGAP(card)
        self.assertEquals(card.Eid(), 899)
        self.assertEquals(card.Pid(), 90)
        card.write_bdf(size, 'dummy')
        card.rawFields()

    def test_pgap_01(self):
        lines = ['PGAP    90                      1.E+5']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = PGAP(card)
        card.write_bdf(size, 'dummy')
        self.assertEquals(card.Pid(), 90)
        card.rawFields()


if __name__ == '__main__':
    unittest.main()
