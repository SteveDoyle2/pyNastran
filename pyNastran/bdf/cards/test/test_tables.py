import unittest
from pyNastran.bdf.bdf import BDF, BDFCard, TABDMP1
from pyNastran.bdf.field_writer_8 import print_card_8

bdf = BDF(debug=False)
class TestTables(unittest.TestCase):

    def test_tabdmp1_01(self):
        lines = ['TABDMP1,100,,,,,,,,+',
                 '+,1e-3,.02,200.,.02,ENDT',]
        card = bdf.process_card(lines)
        card = BDFCard(card)
        #print(card)
        card2 = TABDMP1.add_card(card)
        fields = card2.raw_fields()
        msg = print_card_8(fields).rstrip()
        #print(msg)
        lines_expected = [
            'TABDMP1      100       G',
            #'            .001     .02    200.     .02    ENDT']
            '            1E-3    0.02   200.0    0.02    ENDT']
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
