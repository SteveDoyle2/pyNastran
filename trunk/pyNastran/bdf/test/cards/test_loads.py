import unittest

from pyNastran.bdf.bdf import BDF, BDFCard, DAREA, SUPORT, SUPORT1, CGAP, PGAP, CDAMP1, SET1, AESTAT, DMI, DMIG
from pyNastran.bdf.fieldWriter import print_card

bdf = BDF()

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
        card.write_bdf(size, 'dummy')
        card.rawFields()

    def test_support_01(self):
        lines = ['SUPORT, 1,      126']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = SUPORT(card)
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

    def test_aestat_01(self):
        lines = ['AESTAT  502     PITCH']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = AESTAT(card)
        card.write_bdf(size, 'dummy')
        card.rawFields()

    def test_dmi_01(self):
        lines = ['DMI,Q,0,6,1,0,,4,4']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = DMI(card)
        card.write_bdf(size, 'dummy')
        #card.rawFields()

    def test_dmig_01(self):
        lines = ['DMIG    ENFORCE 0       1       1       0']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = DMIG(card)
        card.write_bdf(size, 'dummy')
        #card.rawFields()

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

    def test_set1_01(self):
        lines = ['SET1,    1100,    100,     101']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = SET1(card)
        card.write_bdf(size, 'dummy')
        card.rawFields()


if __name__ == '__main__':
    unittest.main()

