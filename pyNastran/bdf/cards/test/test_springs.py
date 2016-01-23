import unittest

from pyNastran.bdf.bdf import BDF, BDFCard, PELAS

bdf = BDF(debug=False)
class TestSprings(unittest.TestCase):
    def test_pelas_01(self):
        lines = ['pelas, 201, 1.e+5']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        elem = PELAS()
        elem.add_card(card)
        elem.write_card(size, 'dummy')
        elem.raw_fields()
        self.assertEqual(elem.Pid(), 201)
        self.assertEqual(elem.K(), 1e5)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
