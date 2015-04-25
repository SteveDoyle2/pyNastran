from six.moves import zip, StringIO
import unittest

from itertools import count

from pyNastran.bdf.bdf import BDF, BDFCard, PBEAM, PBAR, CBEAM, GRID, MAT1
from pyNastran.bdf.bdf import CROD, CONROD
from pyNastran.bdf.bdf import PELAS

from pyNastran.bdf.field_writer_8 import print_card_8

bdf = BDF(debug=False)
class TestBeams(unittest.TestCase):
    def test_pbeam_01(self):
        lines = [
            'PBEAM,39,6,2.9,3.5,5.97',
            '     ,  , ,2.0,-4.0',
            '     ,YES,1.0,5.3,56.2,78.6',
            '     ,   ,   ,2.5,-5.0',
            '     ,   ,   ,1.1,    ,2.1,,0.21',
            '     ,   ,   ,   ,    ,0.5,,0.0',
        ]
        card = bdf.process_card(lines)
        #print(print_card_8(card))
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.raw_fields()

        lines_expected = [
            'PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
            '              0.      0.      2.     -4.      0.      0.      0.      0.',
            '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
            '              0.      0.     2.5     -5.      0.      0.      0.      0.',
            '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
            '              0.      0.      0.      0.      .5      .5      0.      0.'
        ]
        self._compare(fields, lines_expected)

    def _compare(self, fields, lines_expected):
        msg = print_card_8(fields).rstrip()
        lines_actual = msg.rstrip().split('\n')
        msgA = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msgA += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msgA)
        for i, actual, expected in zip(count(), lines_actual, lines_expected):
            actual = str(actual)
            expected = str(expected)
            msg = msgA + '\ni=%s' % i + '\nactual   = %r\n' % actual
            msg += 'expected = %r' % expected
            self.assertEqual(actual, expected, msg)

    def test_pbeam_02(self):
        lines = [
            'PBEAM,39,6,2.9,3.5,5.97',
            '     ,  , ,2.0,-4.0',
            '     ,YES,1.0,5.3,56.2,78.6',
            '     ,   ,   ,2.5,-5.0',
            '     ,YES,1.0,5.3,56.2,78.6',
            '     ,   ,   ,2.5,-5.0',
            '     ,YES,1.0,5.3,56.2,78.6',
            '     ,   ,   ,2.5,-5.0',
            '     ,   ,   ,1.1,    ,2.1,,0.21',
            '     ,   ,   ,   ,    ,0.5,,0.0',
        ]

        card = bdf.process_card(lines)
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.raw_fields()

        lines_expected = [
            'PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
            '              0.      0.      2.     -4.      0.      0.      0.      0.',
            '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
            '              0.      0.     2.5     -5.      0.      0.      0.      0.',
            '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
            '              0.      0.     2.5     -5.      0.      0.      0.      0.',
            '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
            '              0.      0.     2.5     -5.      0.      0.      0.      0.',
            '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
            '              0.      0.      0.      0.      .5      .5      0.      0.'
        ]
        self._compare(fields, lines_expected)

    def test_pbeam_03(self):
        lines = [
            'PBEAM,39,6,2.9,3.5,5.97',
            '     ,  , ,2.0,-4.0',
            '     ,YES,1.0,5.3,56.2,78.6',
            '     ,   ,   ,2.5,-5.0',
            '     ,   ,   ,1.1,    ,2.1,,0.21',
            '     ,   ,   ,   ,    ,0.5,,0.0',
        ]

        card = bdf.process_card(lines)
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.raw_fields()

        lines_expected = [
            'PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
            '              0.      0.      2.     -4.      0.      0.      0.      0.',
            '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
            '              0.      0.     2.5     -5.      0.      0.      0.      0.',
            '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
            '              0.      0.      0.      0.      .5      .5      0.      0.',
        ]
        self._compare(fields, lines_expected)

    def test_pbeam_04(self):
        lines = [
            'PBEAM,39,6,2.9,3.5,5.97',
            '     ,  , ,2.0,-4.0',
            '     ,   ,   ,1.1,    ,2.1,,0.21',
            '     ,   ,   ,   ,    ,0.5,,0.0',
        ]

        card = bdf.process_card(lines)
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.raw_fields()
        lines_expected = [
            'PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
            '              0.      0.      2.     -4.      0.      0.      0.      0.',
            '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
            '              0.      0.      0.      0.      .5      .5      0.      0.',
        ]
        self._compare(fields, lines_expected)

    def test_pbeam_05(self):
        lines = [
            'PBEAM,39,6,2.9,3.5,5.97',
            '     ,  , ,2.0,-4.0',
        ]

        card = bdf.process_card(lines)
        #print(print_card_8(card))
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.raw_fields()
        msg = print_card_8(fields)

        lines_expected = [
            'PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
            '              0.      0.      2.     -4.      0.      0.      0.      0.',
            '              1.      1.      0.      0.      0.      0.      0.      0.',
            '              0.      0.      0.      0.      0.      0.      0.      0.',
        ]
        self._compare(fields, lines_expected)

    def test_cbeam_05(self):
        # modification of test_pbeam_05
        model = BDF(debug=False)
        lines = ['PBEAM,3,6,2.9,3.5,5.97,0.4,3.14',
                 '     , , ,2.0,-4.0',]
        card = model.add_card(lines, 'PBEAM', is_list=False)

        lines = ['CBEAM         10       3      1       2      0.01.000000     0.0']
        model.add_card(lines, 'CBEAM', is_list=False)

        lines = ['MAT1, 6, 1.0e7,,0.3']
        model.add_card(lines, 'MAT1', is_list=False)

        lines = ['GRID,1,,0.,0.,0.']
        model.add_card(lines, 'GRID', is_list=False)

        lines = ['GRID,2,,0.,0.,0.']
        model.add_card(lines, 'GRID', is_list=False)
        model.cross_reference()

        cbeam = model.Element(10)
        #cbeam = model.elements[10]
        #print("Area = ", cbeam.Area())
        #print("I11 = ", cbeam.I11())
        #print("I22 = ", cbeam.I22())
        #print("I12 = ", cbeam.I12())
        #print("J = ", cbeam.J())

        #print("Area = ", cbeam.Area())
        #print("I11 = ", cbeam.I1())
        #print("I22 = ", cbeam.I2())
        #print("I12 = ", cbeam.I12())
        #print("J = ", cbeam.J())
        node_ids = cbeam.node_ids
        assert node_ids == [1, 2], node_ids
        self.assertEqual(cbeam.Area(), 2.9)
        self.assertEqual(cbeam.I11(), 3.5)
        self.assertEqual(cbeam.I22(), 5.97)
        self.assertEqual(cbeam.I12(), 0.4)
        self.assertEqual(cbeam.J(), 3.14)

    def test_pbeam_06(self):
        lines = [
            'PBEAM   1       1       1.      60.     1.                              PBEAM1',
            '+BEAM1  5.              -5.                                             PBEAM2',
            '+BEAM2  YES     1.      2.      240.                                    PBEAM3',
            '+BEAM3  10.             -10.                                            PBEAM4',
            '+BEAM4                  -.666667',
        ]

        card = bdf.process_card(lines)
        card = BDFCard(card)
        card2 = PBEAM(card)
        fields = card2.raw_fields()

        lines_expected = [
            'PBEAM          1       1      1.     60.      1.      0.      0.      0.',
            '              5.      0.     -5.      0.      0.      0.      0.      0.',
            '             YES      1.      2.    240.      0.      0.      0.      0.',
            '             10.      0.    -10.      0.      0.      0.      0.      0.',
            '              1.      1.-.666667      0.      0.      0.      0.      0.',
            '              0.      0.      0.      0.      0.      0.      0.      0.',
        ]
        self._compare(fields, lines_expected)

    def test_pbeam_07(self):
        lines = [
            'PBEAM   100     100     1.00    10.     1.0                             +Z1',
            '+Z1     NO      1.0                                                     +Z4',
            '+Z4     0.0     0.0',
        ]
        card = bdf.process_card(lines)
        card = BDFCard(card)
        with self.assertRaises(SyntaxError):  # .. todo:: is this correct?
            card2 = PBEAM(card)

        #if 0:
            #fields = card2.raw_fields()
            #msg = print_card_8(fields)

            #lines_actual = msg.rstrip().split('\n')
            #msg = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
            #msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
            #self.assertEqual(len(lines_actual), len(lines_expected), msg)
            #for actual, expected in zip(lines_actual, lines_expected):
                #msg =  '\nactual   = %r\n' % actual
                #msg += 'expected = %r' % expected
                #self.assertEqual(actual, expected, msg)

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
            fields = card2.raw_fields()
            msg = print_card_8(fields)
            size = 16
            msg = card2.write_card(size, 'dummy')

            lines_actual = msg.rstrip().split('\n')
            msgA = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
            msgA += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
            self.assertEqual(len(lines_actual), len(lines_expected), msg)
            for actual, expected in zip(lines_actual, lines_expected):
                actual = str(actual)
                expected = str(expected)
                msg = msgA + '\nactual   = %r\n' % actual
                msg += 'expected = %r' % expected
                self.assertEqual(actual, expected, msg)

    def test_pbeam_09(self):
        fields = [
            u'PBEAM', 4570049, 4570010, 0.12, 0.000256, 0.005625, None,
            0.0008889, 6.4444e-07, -0.04, -0.75, 0.04, -0.75, 0.04, 0.75,
            -0.04, 0.75, 'YES', 1.0, 0.12, 0.000256, 0.005625, 0.000256,
            None, 0.0008889, 6.4444e-07, -0.04, -0.75, 0.04, -0.75, 0.04,
            0.75, -0.04, 0.853433, 0.849842
        ]
        #fields = [u'PBAR', 1510998, 1520998, 0.0, 4.9000000000000006e-14,
        #4.9000000000000006e-14, 0.0, 0.0, None, 0.0, 0.0, 0.0, 0.0,
        #0.0, 0.0, 0.0, 0.0, None, None, 0.0]
        card = print_card_8(fields)
        #print(card)
        card = print_card_8(fields)
        lines = card.split('\n')
        card = bdf.process_card(lines)
        card2 = BDFCard(card)
        #with self.assertRaises(AssertionError):  # A=0, I12=0, K1=0
        pbeam = PBEAM(card2)
        fields2 = pbeam.repr_fields()
        assert fields == fields


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
