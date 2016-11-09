import unittest
from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.cards.bdf_tables import (
    TABLED1, TABLED2, TABLED3, TABLED4,
    TABLEM1, TABLEM2, TABLEM3, TABLEM4,
    TABDMP1, #TABLES1, TABLEST, TABRND1, TABRNDG,
)
from pyNastran.bdf.field_writer_8 import print_card_8

model = BDF(debug=False)
class TestTables(unittest.TestCase):

    def test_tabdmp1_01(self):
        lines = ['TABDMP1,100,,,,,,,,+',
                 '+,1e-3,.02,200.,.02,ENDT',]
        card = model.process_card(lines)
        card = BDFCard(card)
        #print(card)
        card2 = TABDMP1.add_card(card, comment='table')
        assert card2._comment == '$table\n', '%r' % card2._comment
        assert card2.comment == '$table\n', '%r' % card2.comment
        fields = card2.raw_fields()
        msg = card2.write_card(size=8).rstrip()
        #print(msg)
        lines_expected = [
            '$table',
            'TABDMP1      100       G',
            '            .001     .02    200.     .02    ENDT']
            #'            1E-3    0.02   200.0    0.02    ENDT']
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

        msg = card2.write_card(size=16).rstrip()
        #print(msg)
        lines_expected = [
            '$table',
            'TABDMP1*             100               G',
            '*',
            '*                   .001             .02            200.             .02',
            '*                   ENDT'
        ]
        lines_actual = [line.rstrip() for line in msg.rstrip().split('\n')]
        msg = '\n%r\n\n%r\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

    def test_tabled1(self):
        lines = [
            'TABLED1, 32',
            ',-3.0, 6.9, 2.0, 5.6, 3.0, 5.6, ENDT',
        ]
        card = model.process_card(lines)
        card = BDFCard(card)
        card2 = TABLED1.add_card(card)
        #print(card2)

    def test_tabled2(self):
        lines = [
            'TABLED2, 15, -10.5',
            ',1.0, -4.5, 2.0, -4.2, 2.0, 2.8, 7.0, 6.5',
            ',SKIP, SKIP, 9.0, 6.5, ENDT',
        ]
        card = model.process_card(lines)
        card = BDFCard(card)
        card2 = TABLED2.add_card(card)
        #print(card2)

    def test_tabled3(self):
        lines = [
            'TABLED3, 62, 126.9, 30.0',
            ',2.9, 2.9, 3.6, 4.7, 5.2, 5.7, ENDT',
        ]
        card = model.process_card(lines)
        card = BDFCard(card)
        card2 = TABLED3.add_card(card)
        #print(card2)

    def test_tabled4(self):
        lines = [
            'TABLED4, 28, 0.0, 1.0, 0.0, 100.',
            ',2.91, -0.0329, 6.51-5, 0.0, -3.4-7, ENDT',
        ]
        card = model.process_card(lines)
        card = BDFCard(card)
        card2 = TABLED4.add_card(card)
        #print(card2)

    #def test_tableht(self):
        #lines = [
            #'TABLEHT, 85',
            #'10.0, 101, 25.0, 102, 40.0, 110, ENDT',
        #]
        #card = model.process_card(lines)
        #card = BDFCard(card)
        #card2 = TABLEHT.add_card(card)

    def test_tablem1(self):
        lines = [
            'TABLEM1, 32',
            '-3.0, 6.9, 2.0, 5.6, 3.0, 5.6, ENDT',
        ]
        card = model.process_card(lines)
        card = BDFCard(card)
        card2 = TABLEM1.add_card(card)
        #print(card2)

    def test_tablem2(self):
        lines = [
            'TABLEM2, 15, -10.5',
            ',1.0, -4.5, 2.0, -4.5, 2.0, 2.8, 7.0, 6.5',
            ',SKIP, SKIP, 9.0, 6.5, ENDT',
        ]
        card = model.process_card(lines)
        card = BDFCard(card)
        card2 = TABLEM2.add_card(card)
        #print(card2)

    def test_tablem3(self):
        lines = [
            'TABLEM3, 62, 126.9, 30.0',
            ',2.9, 2.9, 3.6, 4.7, 5.2, 5.7, ENDT',
        ]
        card = model.process_card(lines)
        card = BDFCard(card)
        card2 = TABLEM3.add_card(card)
        #print(card2)

    def test_tablem4(self):
        lines = [
            'TABLEM4, 28, 0.0, 1.0, 0.0, 100.',
            ',2.91, -0.0329, 6.51-5, 0.0, -3.4-7, ENDT',
        ]
        card = model.process_card(lines)
        card = BDFCard(card)
        card2 = TABLEM4.add_card(card)
        #print(card2)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
