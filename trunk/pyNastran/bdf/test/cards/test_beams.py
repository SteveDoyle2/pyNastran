import unittest
from pyNastran.bdf.bdf import BDF, BDFCard, PBEAM
from pyNastran.bdf.fieldWriter import print_card

bdf = BDF()
class TestBeams(unittest.TestCase):

    def test_pbeam_01(self):
        lines =['PBEAM,39,6,2.9,3.5,5.97',
                '     ,  , ,2.0,-4.0',
                '     ,YES,1.0,5.3,56.2,78.6',
                '     ,   ,   ,2.5,-5.0',
                '     ,   ,   ,1.1,    ,2.1,,0.21',
                '     ,   ,   ,   ,    ,0.5,,0.0',]
        card = bdf.process_card(lines)
        #print print_card(card)
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.rawFields()
        msg = print_card(fields).rstrip()

        lines_expected = ['PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
                          '              0.      0.      2.     -4.      0.      0.      0.      0.',
                          '             YES      1.     5.3    56.2    78.6    56.2      0.      0.',
                          '              0.      0.      0.     2.5     -5.      0.      0.      0.',
                          '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
                          '              0.      0.      0.      0.      .5      .5      0.      0.']
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEquals(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEquals(actual, expected)

    def test_pbeam_02(self):
        lines =['PBEAM,39,6,2.9,3.5,5.97',
                '     ,  , ,2.0,-4.0',
                '     ,YES,1.0,5.3,56.2,78.6',
                '     ,   ,   ,2.5,-5.0',
                '     ,YES,1.0,5.3,56.2,78.6',
                '     ,   ,   ,2.5,-5.0',
                '     ,YES,1.0,5.3,56.2,78.6',
                '     ,   ,   ,2.5,-5.0',
                '     ,   ,   ,1.1,    ,2.1,,0.21',
                '     ,   ,   ,   ,    ,0.5,,0.0',]

        card = bdf.process_card(lines)
        #print print_card(card)
        card = BDFCard(card)
        #print "card =", card
        card2 = PBEAM(card)
        fields = card2.rawFields()
        msg = print_card(fields)
        #print msg

        lines_expected = ['PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
                          '              0.      0.      2.     -4.      0.      0.      0.      0.',
                          '             YES      1.     5.3    56.2    78.6    56.2      0.      0.',
                          '              0.      0.      0.     2.5     -5.      0.      0.      0.',
                          '             YES      1.     5.3    56.2    78.6    56.2      0.      0.',
                          '              0.      0.      0.     2.5     -5.      0.      0.      0.',
                          '             YES      1.     5.3    56.2    78.6    56.2      0.      0.',
                          '              0.      0.      0.     2.5     -5.      0.      0.      0.',
                          '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
                          '              0.      0.      0.      0.      .5      .5      0.      0.']
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEquals(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEquals(actual, expected)

    def test_pbeam_03(self):
        lines =['PBEAM,39,6,2.9,3.5,5.97',
                '     ,  , ,2.0,-4.0',
                '     ,YES,1.0,5.3,56.2,78.6',
                '     ,   ,   ,2.5,-5.0',
                '     ,   ,   ,1.1,    ,2.1,,0.21',
                '     ,   ,   ,   ,    ,0.5,,0.0',]

        card = bdf.process_card(lines)
        #print print_card(card)
        card = BDFCard(card)
        #print "card =", card
        card2 = PBEAM(card)
        fields = card2.rawFields()
        msg = print_card(fields)
        #print msg

        lines_expected = ['PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
                          '              0.      0.      2.     -4.      0.      0.      0.      0.',
                          '             YES      1.     5.3    56.2    78.6    56.2      0.      0.',
                          '              0.      0.      0.     2.5     -5.      0.      0.      0.',
                          '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
                          '              0.      0.      0.      0.      .5      .5      0.      0.',]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEquals(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEquals(actual, expected)

    def test_pbeam_04(self):
        lines =['PBEAM,39,6,2.9,3.5,5.97',
                '     ,  , ,2.0,-4.0',
                '     ,   ,   ,1.1,    ,2.1,,0.21',
                '     ,   ,   ,   ,    ,0.5,,0.0',]

        card = bdf.process_card(lines)
        #print print_card(card)
        card = BDFCard(card)
        #print "card =", card
        card2 = PBEAM(card)
        fields = card2.rawFields()
        msg = print_card(fields)
        #print msg

        lines_expected = ['PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
                          '              0.      0.      2.     -4.      0.      0.      0.      0.',
                          '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
                          '              0.      0.      0.      0.      .5      .5      0.      0.',]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEquals(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEquals(actual, expected)

    def test_pbeam_05(self):
        lines =['PBEAM,39,6,2.9,3.5,5.97',
                '     ,  , ,2.0,-4.0',]

        card = bdf.process_card(lines)
        print print_card(card)
        card = BDFCard(card)
        #print "card =", card
        card2 = PBEAM(card)
        fields = card2.rawFields()
        msg = print_card(fields)
        print msg

        lines_expected = ['PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
                          '              0.      0.      2.     -4.      0.      0.      0.      0.',
                          '              1.      1.      0.      0.      0.      0.      0.      0.',
                          '              0.      0.      0.      0.      0.      0.      0.      0.',]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEquals(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            msg =  'actual   = %r\n' % actual
            msg += 'expected = %r' % expected
            self.assertEquals(actual, expected, msg)

    def test_pbeam_06(self):
        lines =['PBEAM   1       1       1.      60.     1.                              PBEAM1',
                '+BEAM1  5.              -5.                                             PBEAM2',
                '+BEAM2  YES     1.      2.      240.                                    PBEAM3',
                '+BEAM3  10.             -10.                                            PBEAM4',
                '+BEAM4                  -.666667',]

        card = bdf.process_card(lines)
        print print_card(card)
        card = BDFCard(card)
        #print "card =", card
        card2 = PBEAM(card)
        fields = card2.rawFields()
        msg = print_card(fields)
        #print msg

        lines_expected = ['PBEAM          1       1      1.     60.      1.      0.      0.      0.',
                          '              5.      0.     -5.      0.      0.      0.      0.      0.',
                          '             YES      1.      2.    240.      0.    240.      0.      0.',
                          '              0.     10.      0.    -10.      0.      0.      0.      0.',
                          '              1.      1.-.666667      0.      0.      0.      0.      0.',
                          '              0.      0.      0.      0.      0.      0.      0.      0.',
        ]
        print '\n'.join(lines_expected)
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEquals(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            msg =  'actual   = %r\n' % actual
            msg += 'expected = %r' % expected
            self.assertEquals(actual, expected, msg)

    def _test_pbeam_07(self):

        lines = ['PBEAM   100     100     1.00    10.     1.0                             +Z1',
                 '+Z1     NO      1.0                                                     +Z4',
                 '+Z4     0.0     0.0',]
        card = bdf.process_card(lines)
        print print_card(card)
        card = BDFCard(card)
        #print "card =", card
        #with self.assertRaises(RuntimeError):  # temporary RuntimeError
        card2 = PBEAM(card)
        
        if 0:
            fields = card2.rawFields()
            msg = print_card(fields)
            print msg

            lines_actual = msg.rstrip().split('\n')
            msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
            msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
            self.assertEquals(len(lines_actual), len(lines_expected), msg)
            for actual, expected in zip(lines_actual, lines_expected):
                msg =  'actual   = %r\n' % actual
                msg += 'expected = %r' % expected
                self.assertEquals(actual, expected, msg)

if __name__ == '__main__':
    unittest.main()