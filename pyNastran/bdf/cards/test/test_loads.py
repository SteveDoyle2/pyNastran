import unittest

from pyNastran.bdf.bdf import BDF, BDFCard, DAREA, PLOAD4
from pyNastran.bdf.bdf import SET1, AESTAT, DMI, DMIG

bdf = BDF(debug=False)

class TestLoads(unittest.TestCase):
    def test_darea_01(self):
        #
        #DAREA SID P1 C1 A1  P2 C2 A2
        #DAREA 3   6   2 8.2 15 1  10.1
        lines = ['DAREA,3,6,2,8.2,15,1,10.1']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = DAREA(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_pload4_01(self):
        lines = ['PLOAD4  1000    1       -60.    -60.    60.             1']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = PLOAD4(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_pload4_02(self):
        lines = ['PLOAD4  1       101     1.                              10000   10011']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = PLOAD4(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_aestat_01(self):
        lines = ['AESTAT  502     PITCH']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = AESTAT(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_dmi_01(self):
        lines = ['DMI,Q,0,6,1,0,,4,4']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = DMI(card)
        card.write_card(size, 'dummy')
        #card.rawFields()

    def test_set1_01(self):
        lines = ['SET1,    1100,    100,     101']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = SET1(card)
        card.write_card(size, 'dummy')
        card.raw_fields()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

