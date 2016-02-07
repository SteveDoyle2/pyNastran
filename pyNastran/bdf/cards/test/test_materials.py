import unittest

from pyNastran.bdf.bdf import BDF, BDFCard, MAT1, MAT8, MAT11
from pyNastran.bdf.field_writer_8 import print_card_8

bdf = BDF(debug=False)

class TestMaterials(unittest.TestCase):
    def test_mat1_01(self):
        #
        #MAT5           1    700.    300.    900.    400.    200.    600.     90.+
        #+             .1
        mid = 1
        E = 2e7
        G = 3e7
        nu = 0.4
        rho = 1e-50
        fields = ['MAT1', mid, E, G, nu, rho]

        card = BDFCard(fields)

        mat1 = MAT1.add_card(card)
        self.assertEqual(mid, mat1.Mid())
        self.assertEqual(E, mat1.E())
        self.assertEqual(G, mat1.G())
        self.assertEqual(nu, mat1.Nu())
        self.assertEqual(rho, mat1.Rho())

        size = 8
        msg = mat1.write_card(size, 'dummy')
        self.assertEqual(msg,
                         'MAT1           1    2.+7    3.+7      .4   1.-50\n')

        size = 16
        expected = 'MAT1*                  1       20000000.       30000000.              .4\n*                  1.-50\n'
        actual = mat1.write_card(size, is_double=False)
        msg = 'actual=\n%s\nexpected=\n%s' % (actual, expected)
        self.assertEqual(actual, expected, msg)

    def test_mat5_01(self):
        #
        #MAT5           1    700.    300.    900.    400.    200.    600.     90.+
        #+             .1
        pass

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
            'MAT8*            4700010        2830000.        1140000.             .55',
            '*                717000.         285194.         285194.',
            '*              .00000917       .00002606             70.',
            '*',
        ]

        card = bdf.process_card(lines)
        #print(print_card_8(card))
        cardi = BDFCard(card)
        #print("card =", card)
        #with self.assertRaises(RuntimeError):  # temporary RuntimeError
        card2 = MAT8.add_card(cardi)

        fields = card2.raw_fields()
        msg = print_card_8(fields)
        #f = StringIO.StringIO()
        size = 16
        msg = card2.write_card(size, 'dummy')
        #msg = f.getvalue()
        #print(msg)

        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        #print(msg)
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            msg = 'actual   = %r\n' % actual
            msg += 'expected = %r' % expected
            self.assertEqual(actual, expected, msg)

    def test_mat11_01(self):
        lines = [
            'MAT11          1    1.+75000000. 700000.      .1     .13     .267000000.+',
            '+       9000000.3000000.      .1    1.-5    7.-6    8.-6     50.',
        ]
        lines_expected = [
            'MAT11          1    1.+75000000. 700000.      .1     .13     .267000000.',
            '        9000000.3000000.      .1  .00001 .000007 .000008     50.'
        ]
        card = bdf.process_card(lines)
        cardi = BDFCard(card)
        card2 = MAT11.add_card(cardi)

        fields = card2.raw_fields()
        msg = print_card_8(fields)
        #f = StringIO.StringIO()
        size = 8
        msg = card2.write_card(size, 'dummy')
        #msg = f.getvalue()
        #print(msg)

        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        #print(msg)
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            msg = '\nactual   = %r\n' % actual
            msg += 'expected =  %r' % expected
            self.assertEqual(actual, expected, msg)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
