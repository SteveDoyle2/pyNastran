import unittest

from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.bdf import CPENTA15

bdf = BDF()

class TestSolids(unittest.TestCase):

    def test_cpenta_01(self):
        lines = ['CPENTA,85,22,201,202,203,205,206,207,+PN2',
                 '+PN2,209,210,217,  ,  ,  ,213,214,218']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CPENTA15(card)
        card.write_bdf(size, 'dummy')
        card.rawFields()


if __name__ == '__main__':
    unittest.main()
