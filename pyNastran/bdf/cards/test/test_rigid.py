import unittest
from pyNastran.bdf.bdf import BDF, BDFCard, RBE1, RBE2, RBE3
from pyNastran.bdf.field_writer_8 import print_card_8

bdf = BDF(debug=False)
class TestRigid(unittest.TestCase):

    def test_rbe3_01(self):
        lines = [
            'rbe3,6, ,3,123456,1.0,123456,41,4,+rbe3',
            '+rbe3,alpha,2.0e-4',
        ]
        card = bdf.process_card(lines)
        card = BDFCard(card)
        rbe = RBE3.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()
        lines_expected = [
            'RBE3           6               3  123456      1.  123456      41       4',
            '           ALPHA   .0002'
        ]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

    #-------------------------------------------------------------------------
    def test_rbe2_01(self):
        lines = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898 117748',
            '+         117765  117764  117763  109821  117743  117744  117750 117751',
            '+         117745  117746  101902    1.-6',
        ]
        card = bdf.process_card(lines)
        card = BDFCard(card)
        rbe = RBE2.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()
        #print(msg)
        lines_expected = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898  117748',
            '          117765  117764  117763  109821  117743  117744  117750  117751',
            '          117745  117746  101902 .000001'
        ]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

    def test_rbe2_02(self):
        lines = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898 117748',
            '+         117765  117764  117763  109821  117743  117744  117750 117751',
            '+         117745  117746  101902    ',
        ]
        card = bdf.process_card(lines)
        card = BDFCard(card)
        rbe = RBE2.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()
        lines_expected = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898  117748',
            '          117765  117764  117763  109821  117743  117744  117750  117751',
            '          117745  117746  101902      0.'
        ]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)
    #-------------------------------------------------------------------------
    def test_rbe1_01(self):
        lines = [
            'RBE1    10201   10201   123     10202   456',
            '           UM   10201   456     10202   123',
        ]

        card = bdf.process_card(lines)
        #print(print_card_8(card))
        card = BDFCard(card)
        rbe = RBE1.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()

        lines_expected = [
            'RBE1       10201   10201     123   10202     456',
            '              UM   10201     456   10202     123'
        ]

        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

    def test_rbe1_02(self):
        lines = [
            'RBE1        1001    1000  123456',
            '              UM    1002     123    1003     123    1004     123',
            '                    1005     123    1006     123    1008     123',
            '                    1009     123    1010     123    1011     123',
            '                    1012     123',
        ]
        card = bdf.process_card(lines)
        #print(print_card_8(card))
        card = BDFCard(card)
        rbe = RBE1.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()

        lines_expected = [
            'RBE1        1001    1000  123456',
            '              UM    1002     123    1003     123    1004     123',
            '                    1005     123    1006     123    1008     123',
            '                    1009     123    1010     123    1011     123',
            '                    1012     123',
        ]


        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

    def test_rbe1_03(self):
        lines = [
            'rbe1,46,3,123456, , , , , ,+rbe46',
            '+rbe46,UM,4,123456,5,123456,2.0-6'
        ]
        card = bdf.process_card(lines)
        card = BDFCard(card)
        rbe = RBE1.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()

        lines_expected = [
            'RBE1          46       3  123456',
            '              UM       4  123456       5  123456 .000002'
        ]

        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
