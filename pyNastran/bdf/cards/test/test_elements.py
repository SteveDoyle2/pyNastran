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
        elem = CBUSH.add_card(card)
        self.assertEqual(elem.Eid(), 101)
        self.assertEqual(elem.Pid(), 102)
        elem.write_card(size, 'dummy')
        elem.raw_fields()

    def test_cdamp1_01(self):
        lines = ['CDAMP1, 2001, 20, 1001, 1']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        elem = CDAMP1.add_card(card)
        self.assertEqual(elem.Eid(), 2001)
        self.assertEqual(elem.Pid(), 20)
        node_ids = elem.node_ids
        assert node_ids == [1001, None], node_ids
        elem.write_card(size, 'dummy')
        elem.raw_fields()

    def test_cgap_01(self):
        lines = ['CGAP    899     90      21      99      0.      1.      0.      0']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        elem = CGAP.add_card(card)
        node_ids = elem.node_ids
        assert node_ids == [21, 99], node_ids
        self.assertEqual(elem.Eid(), 899)
        self.assertEqual(elem.Pid(), 90)
        elem.write_card(size, 'dummy')
        elem.raw_fields()

    def test_pgap_01(self):
        lines = ['PGAP    90                      1.E+5']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        elem = PGAP.add_card(card)
        elem.write_card(size, 'dummy')
        self.assertEqual(elem.Pid(), 90)
        elem.raw_fields()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
