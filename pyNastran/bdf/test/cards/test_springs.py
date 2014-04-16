import unittest

from pyNastran.bdf.bdf import BDF, BDFCard, PELAS
#from pyNastran.bdf.fieldWriter import print_card

bdf = BDF()
class TestSprings(unittest.TestCase):
    def test_pelas_01(self):
        lines = ['pelas, 201, 1.e+5']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = PELAS(card)
        card.write_bdf(size, 'dummy')
        card.rawFields()
        self.assertEquals(card.Pid(), 201)
        self.assertEquals(card.K(), 1e5)

if __name__ == '__main__':
    unittest.main()
