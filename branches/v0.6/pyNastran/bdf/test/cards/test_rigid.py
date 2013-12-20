import unittest
from pyNastran.bdf.bdf import BDF, BDFCard, RBE1, RBE2
from pyNastran.bdf.fieldWriter import print_card

bdf = BDF()
class TestRigid(unittest.TestCase):

    def test_rbe2_01(self):
        lines = ['RBE2      100045  166007  123456  117752  101899  117766  101898 117748',
                 '+         117765  117764  117763  109821  117743  117744  117750 117751',
                 '+         117745  117746  101902    1.-6',]
        card = bdf.process_card(lines)
        card = BDFCard(card)
        #print(card)
        card2 = RBE2(card)
        fields = card2.rawFields()
        msg = print_card(fields).rstrip()
        #print(msg)
        lines_expected = ['RBE2      100045  166007  123456  117752  101899  117766  101898  117748',
                          '          117765  117764  117763  109821  117743  117744  117750  117751',
                          '          117745  117746  101902 .000001']
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected)
    def test_rbe2_02(self):
        lines = ['RBE2      100045  166007  123456  117752  101899  117766  101898 117748',
                 '+         117765  117764  117763  109821  117743  117744  117750 117751',
                 '+         117745  117746  101902    ',]
        card = bdf.process_card(lines)
        card = BDFCard(card)
        #print(card)
        card2 = RBE2(card)
        fields = card2.rawFields()
        msg = print_card(fields).rstrip()
        #print(msg)
        lines_expected = ['RBE2      100045  166007  123456  117752  101899  117766  101898  117748',
                          '          117765  117764  117763  109821  117743  117744  117750  117751',
                          '          117745  117746  101902      0.']
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected)
    #-------------------------------------------------------------------------
    def test_rbe1_01(self):
        lines = ['RBE1    10201   10201   123     10202   456',
                 '           UM   10201   456     10202   123',]

        card = bdf.process_card(lines)
        #print(print_card(card))
        card = BDFCard(card)
        #print(card)
        card2 = RBE1(card)
        fields = card2.rawFields()
        msg = print_card(fields).rstrip()
        #print(msg)

        lines_expected = ['RBE1       10201   10201     123   10202     456',
                          '              UM   10201     456   10202     123']

        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected)

    def test_rbe1_02(self):
        lines = ['RBE1        1001    1000  123456',
                 '              UM    1002     123    1003     123    1004     123',
                 '                    1005     123    1006     123    1008     123',
                 '                    1009     123    1010     123    1011     123',
                 '                    1012     123',]
        card = bdf.process_card(lines)
        #print(print_card(card))
        card = BDFCard(card)
        #print(card)
        card2 = RBE1(card)
        fields = card2.rawFields()
        msg = print_card(fields).rstrip()
        #print(msg)

        lines_expected = ['RBE1        1001    1000  123456',
                          '              UM    1002     123    1003     123    1004     123',
                          '                    1005     123    1006     123    1008     123',
                          '                    1009     123    1010     123    1011     123',
                          '                    1012     123',]


        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected)

if __name__ == '__main__':
    unittest.main()