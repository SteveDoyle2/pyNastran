import unittest

from pyNastran.bdf.bdf import BDF, BDFCard, PELAS

bdf = BDF(debug=False)
class TestSprings(unittest.TestCase):
    def test_pelas_01(self):
        lines = ['pelas, 201, 1.e+5']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = PELAS(card)
        card.write_card(size, 'dummy')
        card.raw_fields()
        self.assertEqual(card.Pid(), 201)
        self.assertEqual(card.K(), 1e5)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
