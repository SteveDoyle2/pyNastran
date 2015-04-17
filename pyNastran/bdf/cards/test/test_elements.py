import unittest

from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.bdf import CGAP, PGAP, CDAMP1, CBUSH

bdf = BDF(debug=False)

class TestElements(unittest.TestCase):

    def test_cbush_01(self):
        lines = ['cbush,101,102,1,,,,,0']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CBUSH(card)
        self.assertEqual(card.Eid(), 101)
        self.assertEqual(card.Pid(), 102)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_cdamp1_01(self):
        lines = ['CDAMP1, 2001, 20, 1001, 1']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CDAMP1(card)
        self.assertEqual(card.Eid(), 2001)
        self.assertEqual(card.Pid(), 20)
        node_ids = card.node_ids
        assert node_ids == [1001, 0], node_ids
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_cgap_01(self):
        lines = ['CGAP    899     90      21      99      0.      1.      0.      0']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CGAP(card)
        node_ids = card.node_ids
        assert node_ids == [21, 99], node_ids
        self.assertEqual(card.Eid(), 899)
        self.assertEqual(card.Pid(), 90)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_pgap_01(self):
        lines = ['PGAP    90                      1.E+5']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = PGAP(card)
        card.write_card(size, 'dummy')
        self.assertEqual(card.Pid(), 90)
        card.raw_fields()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
