import unittest
from pyNastran.bdf.bdf import BDF, BDFCard, TABDMP1
from pyNastran.bdf.fieldWriter import print_card

bdf = BDF()
class TestTables(unittest.TestCase):


    def test_tabdmp1_01(self):
        lines = ['TABDMP1,100,,,,,,,,+',
                 '+,1e-3,.02,200.,.02,ENDT',]
        card = bdf.process_card(lines)
        card = BDFCard(card)
        #print(card)
        card2 = TABDMP1(card)
        fields = card2.rawFields()
        msg = print_card(fields).rstrip()
        #print(msg)
        lines_expected = ['TABDMP1      100       G',
                           #2345678#2345678#2345678#2345678#2345678
                          '            .001     .02    200.     .02    ENDT']
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)


if __name__ == '__main__':
    unittest.main()