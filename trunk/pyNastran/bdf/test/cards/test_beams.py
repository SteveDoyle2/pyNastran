import StringIO
import unittest

from itertools import izip, count

from pyNastran.bdf.bdf import BDF, BDFCard, PBEAM, PBAR
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

        lines_expected = ['PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
                          '              0.      0.      2.     -4.      0.      0.      0.      0.',
                          '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
                          '              0.      0.     2.5     -5.      0.      0.      0.      0.',
                          '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
                          '              0.      0.      0.      0.      .5      .5      0.      0.']
        self._compare(fields, lines_expected)

    def _compare(self, fields, lines_expected):
        msg = print_card(fields).rstrip()
        lines_actual = msg.rstrip().split('\n')
        msgA = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msgA += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msgA)
        for i, actual, expected in izip(count(), lines_actual, lines_expected):
            actual = str(actual)
            expected = str(expected)
            msg = msgA + '\ni=%s' % i + '\nactual   = %r\n' % actual
            msg += 'expected = %r' % expected
            self.assertEqual(actual, expected, msg)

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
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.rawFields()

        lines_expected = ['PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
                          '              0.      0.      2.     -4.      0.      0.      0.      0.',
                          '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
                          '              0.      0.     2.5     -5.      0.      0.      0.      0.',
                          '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
                          '              0.      0.     2.5     -5.      0.      0.      0.      0.',
                          '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
                          '              0.      0.     2.5     -5.      0.      0.      0.      0.',
                          '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
                          '              0.      0.      0.      0.      .5      .5      0.      0.']
        self._compare(fields, lines_expected)

    def test_pbeam_03(self):
        lines =['PBEAM,39,6,2.9,3.5,5.97',
                '     ,  , ,2.0,-4.0',
                '     ,YES,1.0,5.3,56.2,78.6',
                '     ,   ,   ,2.5,-5.0',
                '     ,   ,   ,1.1,    ,2.1,,0.21',
                '     ,   ,   ,   ,    ,0.5,,0.0',]

        card = bdf.process_card(lines)
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.rawFields()

        lines_expected = ['PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
                          '              0.      0.      2.     -4.      0.      0.      0.      0.',
                          '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
                          '              0.      0.     2.5     -5.      0.      0.      0.      0.',
                          '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
                          '              0.      0.      0.      0.      .5      .5      0.      0.',]
        self._compare(fields, lines_expected)

    def test_pbeam_04(self):
        lines =['PBEAM,39,6,2.9,3.5,5.97',
                '     ,  , ,2.0,-4.0',
                '     ,   ,   ,1.1,    ,2.1,,0.21',
                '     ,   ,   ,   ,    ,0.5,,0.0',]

        card = bdf.process_card(lines)
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.rawFields()
        lines_expected = ['PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
                          '              0.      0.      2.     -4.      0.      0.      0.      0.',
                          '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
                          '              0.      0.      0.      0.      .5      .5      0.      0.',]
        self._compare(fields, lines_expected)

    def test_pbeam_05(self):
        lines =['PBEAM,39,6,2.9,3.5,5.97',
                '     ,  , ,2.0,-4.0',]

        card = bdf.process_card(lines)
        #print(print_card(card))
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.rawFields()
        msg = print_card(fields)

        lines_expected = ['PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
                          '              0.      0.      2.     -4.      0.      0.      0.      0.',
                          '              1.      1.      0.      0.      0.      0.      0.      0.',
                          '              0.      0.      0.      0.      0.      0.      0.      0.',]
        self._compare(fields, lines_expected)

    def test_pbeam_06(self):
        lines =['PBEAM   1       1       1.      60.     1.                              PBEAM1',
                '+BEAM1  5.              -5.                                             PBEAM2',
                '+BEAM2  YES     1.      2.      240.                                    PBEAM3',
                '+BEAM3  10.             -10.                                            PBEAM4',
                '+BEAM4                  -.666667',]

        card = bdf.process_card(lines)
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.rawFields()

        lines_expected = ['PBEAM          1       1      1.     60.      1.      0.      0.      0.',
                          '              5.      0.     -5.      0.      0.      0.      0.      0.',
                          '             YES      1.      2.    240.      0.      0.      0.      0.',
                          '             10.      0.    -10.      0.      0.      0.      0.      0.',
                          '              1.      1.-.666667      0.      0.      0.      0.      0.',
                          '              0.      0.      0.      0.      0.      0.      0.      0.',
        ]
        self._compare(fields, lines_expected)

    def _test_pbeam_07(self):
        lines = ['PBEAM   100     100     1.00    10.     1.0                             +Z1',
                 '+Z1     NO      1.0                                                     +Z4',
                 '+Z4     0.0     0.0',]
        card = bdf.process_card(lines)
        card = BDFCard(card)
        #with self.assertRaises(RuntimeError):  # temporary RuntimeError
        card2 = PBEAM(card)

        if 0:
            fields = card2.rawFields()
            msg = print_card(fields)

            lines_actual = msg.rstrip().split('\n')
            msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
            msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
            self.assertEqual(len(lines_actual), len(lines_expected), msg)
            for actual, expected in zip(lines_actual, lines_expected):
                msg =  '\nactual   = %r\n' % actual
                msg += 'expected = %r' % expected
                self.assertEqual(actual, expected, msg)

    def test_pbeam_08(self):  # should fail...
        lines = [
        'PBEAM*   4570049         4570010        .12             2.56-4          *    HRY',
        '*    HRY.005625                         8.889-4         6.4444-7        *    HRZ',
        '*    HRZ-.04            -.75            .04             -.75            *    HSA',
        '*    HSA.04             .75             -.04            .75             *    HSB',
        '*    HSB YES            1.              .12             2.56-4          *    HSC',
        '*    HSC.005625                         8.889-4         6.4444-7        *    HSD',
        '*    HSD-.04            -.75            .04             -.75            *    HSE',
        '*    HSE.04             .75             -.04            .75             *    HSF',
        '*    HSF.853433         .849842                                         *    HSG',
        '*    HSG',
        ]
        lines_expected = [
        'PBEAM*           4570049         4570010             .12         .000256',
        '*                .005625                        .0008889    .00000064444',
        '*                   -.04            -.75             .04            -.75',
        '*                    .04             .75            -.04             .75',
        '*                    YES              1.             .12         .000256',
        '*                .005625                        .0008889    .00000064444',
        '*                   -.04            -.75             .04            -.75',
        '*                    .04             .75            -.04             .75',
        '*                .853433         .849842',
        '*',
        ]

        card = bdf.process_card(lines)
        card = BDFCard(card)
        #with self.assertRaises(RuntimeError):  # temporary RuntimeError
        card2 = PBEAM(card)

        if 1:
            fields = card2.rawFields()
            msg = print_card(fields)
            size = 16
            msg = card2.write_bdf(size, 'dummy')

            lines_actual = msg.rstrip().split('\n')
            msgA = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
            msgA += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
            self.assertEqual(len(lines_actual), len(lines_expected), msg)
            for actual, expected in zip(lines_actual, lines_expected):
                actual = str(actual)
                expected = str(expected)
                msg =  msgA + '\nactual   = %r\n' % actual
                msg += 'expected = %r' % expected
                self.assertEqual(actual, expected, msg)

    def test_pbeam_09(self):
        #pid=
        fields = [u'PBEAM', 4570049, 4570010, 0.12, 0.000256, 0.005625, None,
            0.0008889, 6.4444e-07, -0.04, -0.75, 0.04, -0.75, 0.04, 0.75,
            -0.04, 0.75, 'YES', 1.0, 0.12, 0.000256, 0.005625, 0.000256,
            None, 0.0008889, 6.4444e-07, -0.04, -0.75, 0.04, -0.75, 0.04,
            0.75, -0.04, 0.853433, 0.849842]
        #fields = [u'PBAR', 1510998, 1520998, 0.0, 4.9000000000000006e-14, 4.9000000000000006e-14, 0.0, 0.0, None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None, None, 0.0]
        card = print_card(fields)
        #print card
        card = print_card(fields)
        lines = card.split('\n')
        card = bdf.process_card(lines)
        card2 = BDFCard(card)
        #with self.assertRaises(AssertionError):  # A=0, I12=0, K1=0
        pbeam = PBEAM(card2)
        fields2 = pbeam.reprFields()
        assert fields == fields

class TestBars(unittest.TestCase):
    def test_pbar_01(self):
        #pid=
        fields = [u'PBAR', 1510998, 1520998, 0.0, 4.9000000000000006e-14, 4.9000000000000006e-14, 0.0, 0.0, None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None, None, 0.0]
        card = print_card(fields)
        #print card
        card = print_card(fields)
        lines = card.split('\n')
        card = bdf.process_card(lines)
        card2 = BDFCard(card)
        with self.assertRaises(AssertionError):  # A=0, I12=0, K1=0
            pbar = PBAR(card2)

    def test_pbar_02(self):
        pid = 1
        mid = 2
        A = None
        I1 = I2 = None
        J = None
        nsm = None
        c1=c2=d1=d2=e1=e2=f1=f2=None
        k1=k2=None
        i12 = 3.
        fields = ['PBAR', pid, mid, A, I1, I2, J, nsm, None,
                          c1, c2, d1, d2, e1, e2, f1, f2,
                          k1, k2, i12]
        card = print_card(fields)
        lines = card.split('\n')
        card = bdf.process_card(lines)
        #print(card)
        card2 = BDFCard(card)
        #print(card2)
        pbar = PBAR(card2)
        self.assertEqual(pbar.pid, 1)
        self.assertEqual(pbar.mid, 2)
        self.assertEqual(pbar.A, 0.0)
        self.assertEqual(pbar.i1, 0.0)
        self.assertEqual(pbar.i2, 0.0)
        self.assertEqual(pbar.j, 0.0)
        self.assertEqual(pbar.nsm, 0.0)
        self.assertEqual(pbar.i12, 3.0)
        self.assertEqual(pbar.C1, 0.0)
        self.assertEqual(pbar.C2, 0.0)
        self.assertEqual(pbar.D1, 0.0)
        self.assertEqual(pbar.D2, 0.0)
        self.assertEqual(pbar.E1, 0.0)
        self.assertEqual(pbar.E2, 0.0)
        self.assertEqual(pbar.K1, 1e8)
        self.assertEqual(pbar.K2, 1e8)


if __name__ == '__main__':
    unittest.main()