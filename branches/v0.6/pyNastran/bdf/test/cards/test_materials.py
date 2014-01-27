import unittest

from pyNastran.bdf.bdf import BDF, BDFCard, MAT8
from pyNastran.bdf.fieldWriter import print_card

bdf = BDF()

class TestMaterials(unittest.TestCase):

    def test_mat8_01(self):  # should fail...
        #lines = [  # fails???
        #    'MAT8*    4700007        1675.47         1675.47         .33             *   LHIG',
        #    '*   LHIG28.2            210000.         78000.                          *   LHIH',
        #    '*   LHIH1.32-5          1.32-5          75.             1.943           *   LHII',
        #    '*   LHII1.943           1.943           1.943           3.35',
        #]
        lines = [  # fails
            'MAT8*    4700010        2.83+6          1.14+6          .55             *   LHIJ',
            '*   LHIJ717000.         285194.         285194.                         *   LHIK',
            '*   LHIK9.17-6          2.606-5         70.                             *   LHIL',
            '*   LHIL',
        ]
        lines_expected = [
            'MAT8*    4700010        2.83+6          1.14+6          .55             *',
            '*       717000.         285194.         285194.                         *',
            '*       9.17-6          2.606-5         70.                             *',
            '*       ',
        ]

        card = bdf.process_card(lines)
        #print(print_card(card))
        card = BDFCard(card)
        #print("card =", card)
        #with self.assertRaises(RuntimeError):  # temporary RuntimeError
        card2 = MAT8(card)

        fields = card2.rawFields()
        msg = print_card(fields)
        #f = StringIO.StringIO()
        size = 16
        msg = card2.write_bdf(size, 'dummy')
        #msg = f.getvalue()
        print(msg)

        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            msg =  'actual   = %r\n' % actual
            msg += 'expected = %r' % expected
            self.assertEqual(actual, expected, msg)

if __name__ == '__main__':
    unittest.main()
