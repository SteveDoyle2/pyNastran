"""
tests:
  - CBEAM, PBEAM, PBEAML, PBCOMP, PBMSECT
  - CBEAM3, PBEAM3

"""
import os
import unittest
from itertools import count

import numpy as np

try:
    import matplotlib
    IS_MATPLOTLIB = True
except ImportError:  # pragma: no cover
    IS_MATPLOTLIB = False

if IS_MATPLOTLIB:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

from pyNastran.bdf.bdf import BDF, BDFCard, PBEAM, PBEND, PBMSECT, PBRSECT
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.cards.properties.bars import get_beam_sections

class TestBeams(unittest.TestCase):
    """
    tests:
      - CBEAM, PBEAM, PBEAML, PBCOMP, PBMSECT
      - CBEAM3, PBEAM3

    """
    def test_pbeam_01(self):
        """tests a nasty PBEAM"""
        lines = [
            'PBEAM,39,6,2.9,3.5,5.97',
            '     ,  , ,2.0,-4.0',
            '     ,YES,1.0,5.3,56.2,78.6',
            '     ,   ,   ,2.5,-5.0',
            '     ,   ,   ,1.1,    ,2.1,,0.21',
            '     ,   ,   ,   ,    ,0.5,,0.0',
        ]
        bdf = BDF(debug=False)
        card = bdf._process_card(lines)
        #print(print_card_8(card))
        cardi = BDFCard(card)
        card2 = PBEAM.add_card(cardi)
        fields = card2.raw_fields()

        lines_expected = [
            'PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
            '              0.      0.      2.     -4.      0.      0.      0.      0.',
            '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
            '              0.      0.     2.5     -5.      0.      0.      0.      0.',
            '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
            '              0.      0.      0.      0.      .5      0.      0.      0.'
        ]
        self._compare(fields, lines_expected)

    def _compare(self, fields, lines_expected):
        msg = print_card_8(fields).rstrip()
        lines_actual = msg.rstrip().split('\n')
        msg_a = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg_a += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg_a)
        for i, actual, expected in zip(count(), lines_actual, lines_expected):
            actual = str(actual)
            expected = str(expected)
            msg = msg_a + '\ni=%s' % i + '\nactual   = %r\n' % actual
            msg += 'expected = %r' % expected
            self.assertEqual(actual, expected, msg)

    def test_pbeam_02(self):
        """tests a nasty PBEAM"""
        bdf = BDF(debug=False)
        lines = [
            'PBEAM,39,6,2.9,3.5,5.97',
            '     ,  , ,2.0,-4.0',
            '     ,YES,1.0,5.3,56.2,78.6',
            '     ,   ,   ,2.5,-5.0',
            #'     ,YES,1.0,5.3,56.2,78.6',
            #'     ,   ,   ,2.5,-5.0',
            #'     ,YES,1.0,5.3,56.2,78.6',
            #'     ,   ,   ,2.5,-5.0',
            '     ,   ,   ,1.1,    ,2.1,,0.21',
            '     ,   ,   ,   ,    ,0.5,,0.0',
        ]

        card = bdf._process_card(lines)
        cardi = BDFCard(card)
        card2 = PBEAM.add_card(cardi)
        fields = card2.raw_fields()

        lines_expected = [
            'PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
            '              0.      0.      2.     -4.      0.      0.      0.      0.',
            '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
            '              0.      0.     2.5     -5.      0.      0.      0.      0.',
            #'             YES      1.     5.3    56.2    78.6      0.      0.      0.',
            #'              0.      0.     2.5     -5.      0.      0.      0.      0.',
            #'             YES      1.     5.3    56.2    78.6      0.      0.      0.',
            #'              0.      0.     2.5     -5.      0.      0.      0.      0.',
            '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
            '              0.      0.      0.      0.      .5      0.      0.      0.'
        ]
        self._compare(fields, lines_expected)

    def test_pbeam_03(self):
        """tests a nasty PBEAM"""
        bdf = BDF(debug=False)
        lines = [
            'PBEAM,39,6,2.9,3.5,5.97',
            '     ,  , ,2.0,-4.0',
            '     ,YES,1.0,5.3,56.2,78.6',
            '     ,   ,   ,2.5,-5.0',
            '     ,   ,   ,1.1,    ,2.1,,0.21',
            '     ,   ,   ,   ,    ,0.5,,0.0',
        ]

        card = bdf._process_card(lines)
        cardi = BDFCard(card)
        card2 = PBEAM.add_card(cardi)
        fields = card2.raw_fields()

        lines_expected = [
            'PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
            '              0.      0.      2.     -4.      0.      0.      0.      0.',
            '             YES      1.     5.3    56.2    78.6      0.      0.      0.',
            '              0.      0.     2.5     -5.      0.      0.      0.      0.',
            '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
            '              0.      0.      0.      0.      .5      0.      0.      0.',
        ]
        self._compare(fields, lines_expected)

    def test_pbeam_04(self):
        """tests a nasty PBEAM"""
        bdf = BDF(debug=False)
        lines = [
            'PBEAM,39,6,2.9,3.5,5.97',
            '     ,  , ,2.0,-4.0',
            '     ,   ,   ,1.1,    ,2.1,,0.21',
            '     ,   ,   ,   ,    ,0.5,,0.0',
        ]

        card = bdf._process_card(lines)
        cardi = BDFCard(card)
        card2 = PBEAM.add_card(cardi)
        fields = card2.raw_fields()
        lines_expected = [
            'PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
            '              0.      0.      2.     -4.      0.      0.      0.      0.',
            '              1.      1.     1.1      0.     2.1     2.1     .21     .21',
            '              0.      0.      0.      0.      .5      0.      0.      0.',
        ]
        self._compare(fields, lines_expected)

    def test_pbeam_05(self):
        """tests a nasty PBEAM"""
        bdf = BDF(debug=False)
        lines = [
            'PBEAM,39,6,2.9,3.5,5.97',
            '     ,  , ,2.0,-4.0',
        ]

        card = bdf._process_card(lines)
        #print(print_card_8(card))
        cardi = BDFCard(card)
        card2 = PBEAM.add_card(cardi)
        fields = card2.raw_fields()
        #msg = print_card_8(fields)

        lines_expected = [
            'PBEAM         39       6     2.9     3.5    5.97      0.      0.      0.',
            '              0.      0.      2.     -4.      0.      0.      0.      0.',
            '              1.      1.      0.      0.      0.      0.      0.      0.',
            '              0.      0.      0.      0.      0.      0.      0.      0.',
        ]
        self._compare(fields, lines_expected)

    def test_pbeam_06(self):
        """tests a nasty PBEAM"""
        bdf = BDF(debug=False)
        lines = [
            'PBEAM   1       1       1.      60.     1.                              PBEAM1',
            '+BEAM1  5.              -5.                                             PBEAM2',
            '+BEAM2  YES     1.      2.      240.                                    PBEAM3',
            '+BEAM3  10.             -10.                                            PBEAM4',
            '+BEAM4                  -.666667',
        ]

        card = bdf._process_card(lines)
        cardi = BDFCard(card)
        card2 = PBEAM.add_card(cardi)
        fields = card2.raw_fields()

        lines_expected = [
            'PBEAM          1       1      1.     60.      1.      0.      0.      0.',
            '              5.      0.     -5.      0.      0.      0.      0.      0.',
            '             YES      1.      2.    240.      1.      0.      0.      0.',
            '             10.      0.    -10.      0.      0.      0.      0.      0.',
            '              1.      1.-.666667      0.      0.      0.      0.      0.',
            '              0.      0.      0.      0.      0.      0.      0.      0.',
        ]
        self._compare(fields, lines_expected)

    def test_pbeam_07(self):
        """tests a nasty PBEAM"""
        bdf = BDF(debug=False)
        lines = [
            'PBEAM   100     100     1.00    10.     1.0                             +Z1',
            '+Z1     NO      1.0                                                     +Z4',
            '+Z4     0.0     0.0',
        ]
        card = bdf._process_card(lines)
        cardi = BDFCard(card)
        PBEAM.add_card(cardi)

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

    def test_pbeam_08(self):
        """tests a nasty PBEAM"""
        bdf = BDF(debug=False)
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

        card = bdf._process_card(lines)
        cardi = BDFCard(card)
        card2 = PBEAM.add_card(cardi)

        fields = card2.raw_fields()
        msg = print_card_8(fields)
        size = 16
        msg = card2.write_card(size, 'dummy')

        lines_actual = msg.rstrip().split('\n')
        msg_a = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg_a += 'nlines_actual=%i nlines_expected=%i' % (
            len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            actual = str(actual)
            expected = str(expected)
            msg = msg_a + '\nactual   = %r\n' % actual
            msg += 'expected = %r' % expected
            self.assertEqual(actual, expected, msg)

    def test_pbeam_09(self):
        """tests a nasty PBEAM"""
        bdf = BDF(debug=False)
        fields = [
            u'PBEAM', 4570049, 4570010, 0.12, 0.000256, 0.005625, None, 0.0008889, 6.4444e-07,
            -0.04, -0.75, 0.04, -0.75, 0.04, 0.75, -0.04, 0.75,
            'YES', 1.0, 0.12, 0.000256, 0.005625, 0.000256, None, 0.0008889,
            6.4444e-07, -0.04, -0.75, 0.04, -0.75, 0.04, 0.75, -0.04,
            0.853433, 0.849842
        ]
        #fields = [u'PBAR', 1510998, 1520998, 0.0, 4.9000000000000006e-14,
        #4.9000000000000006e-14, 0.0, 0.0, None, 0.0, 0.0, 0.0, 0.0,
        #0.0, 0.0, 0.0, 0.0, None, None, 0.0]
        card = print_card_8(fields)
        #print(card)
        card = print_card_8(fields)
        lines = card.split('\n')
        card = bdf._process_card(lines)
        cardi = BDFCard(card)
        #with self.assertRaises(AssertionError):  # A=0, I12=0, K1=0
        pbeam = PBEAM.add_card(cardi)
        pbeam.repr_fields()
        assert fields == fields

    def test_pbeam_10(self):
        """
        The number of continuation lines on a PBEAM for a NO/YESA is 0.
        The number of continuation lines on a PBEAM for a YES is 1.
        """
        model = BDF(debug=False)

        # correct - standard NO with full list of values
        lines = [
            'PBEAM          1       1 5.094+7 289940.1.6043+7         271610. 3.73058',
            '            0.',
            '              NO     1.0 .7489+7 238250.1.3182+7   1.-12 223170.3.458069',
            '              0.      0.           .872    .718',
            '              0.  .33936      0. .31983',
        ]
        #lines_expected = lines
        model.add_card(lines, 'PBEAM', is_list=False)

        # correct - 2nd continuation line
        lines = [
            'PBEAM          2       1 5.094+7 289940.1.6043+7         271610. 3.73058',
            '+',
            '              NO     1.0 .7489+7 238250.1.3182+7   1.-12 223170.3.458069',
            '              0.      0.           .872    .718',
            '              0.  .33936      0. .31983',
        ]
        #lines_expected = lines
        model.add_card(lines, 'PBEAM', is_list=False)

        # BAD x/xb value of 1.4 (should be 0 to 1.0)
        lines = [
            'PBEAM          2       1 5.094+7 289940.1.6043+7         271610. 3.73058',
            '+',
            '              NO     1.4 .7489+7 238250.1.3182+7   1.-12 223170.3.458069',
            '              0.      0.           .872    .718',
            '              0.  .33936      0. .31983',
        ]
        #lines_expected = lines
        with self.assertRaises(AssertionError):
            model.add_card(lines, 'PBEAM', is_list=False)

        # error - 3 lines after NO
        lines = [
            'PBEAM          3       1 5.094+7 289940.1.6043+7         271610. 3.73058',
            '            0.',
            '              NO     1.0 .7489+7 238250.1.3182+7   1.-12 223170.3.458069',
            '            0.',
            '              0.      0.           .872    .718',
            '              0.  .33936      0. .31983',
        ]
        with self.assertRaises(RuntimeError):
            model.add_card(lines, 'PBEAM', is_list=False)

        # correct - skipped 2nd line
        lines = [
            'PBEAM          4       1 5.094+7 289940.1.6043+7         271610. 3.73058',
            '              NO     1.0 .7489+7 238250.1.3182+7   1.-12 223170.3.458069',
            '              0.      0.           .872    .718',
            '              0.  .33936      0. .31983',
        ]
        #lines_expected = lines
        model.add_card(lines, 'PBEAM', is_list=False)

        # correct - skipped 2nd line and last line
        lines = [
            'PBEAM          5       1 5.094+7 289940.1.6043+7         271610. 3.73058',
            '              NO     1.0 .7489+7 238250.1.3182+7   1.-12 223170.3.458069',
            '              0.      0.           .872    .718',
        ]
        #lines_expected = lines
        model.add_card(lines, 'PBEAM', is_list=False)

        # correct - skipped 2nd line and last 2 lines
        lines = [
            'PBEAM          6       1 5.094+7 289940.1.6043+7         271610. 3.73058',
            '              NO     1.0 .7489+7 238250.1.3182+7   1.-12 223170.3.458069',
        ]
        #lines_expected = lines
        model.add_card(lines, 'PBEAM', is_list=False)

        # correct - single line
        lines = [
            'PBEAM          7       1 5.094+7 289940.1.6043+7         271610. 3.73058',
        ]
        #lines_expected = lines
        model.add_card(lines, 'PBEAM', is_list=False)

    def test_pbeam_11(self):
        """tests a nasty PBEAM"""
        model = BDF(debug=False)
        lines = [
            #'PBEAM         1       1    1.e8    1.e8    1.e8              10.      1.',
            #'             NO      1.    3.e5    3.e8    3.e8              10.      1.',
            'PBEAM          1       1    1.+8    1.+8    1.+8             10.      1.',
            '              NO      1. 300000.    3.+8    3.+8             10.      1.',
        ]
        model.add_card(lines, 'PBEAM', is_list=False)
        lines_expected = [
            'PBEAM          1       1    1.+8    1.+8    1.+8             10.      1.',
            '+',
            '              NO      1. 300000.    3.+8    3.+8             10.      1.',
        ]
        prop = model.properties[1]
        #print(prop.raw_fields())
        lines_actual = prop.write_card().split('\n')
        msg_a = ''
        for line_expected, line_actual in zip(lines_expected, lines_actual):
            #assert line_expected == line_actual, line_actual
            actual = str(line_actual)
            expected = str(line_expected)
            msg = msg_a + '\nactual   = %r\n' % actual
            msg += 'expected = %r' % expected
            self.assertEqual(actual, expected, msg)

    def test_pbeaml_01(self):
        model = BDF(debug=False)
        model.validate()

        fields = [
            'PBEAML', '622', '623', None, 'BAR', None, None, None, None,
            '.238', '.238', None,
            'NO', '.166', '.238', '.238', None,
            'NO', '.510', '.126', '.126', None,
            'NO', '.565', '.117', '.117', None,
            'NO', '.619', '.126', '.126', None,
            'NO', '.963', '.238', '.238', None,
            'NO', '1.0', '.238', '.238', None,
            'NO']
        model.add_card(fields, fields[0])

        fields = [
            'PBEAML', '623', '623', None, 'BAR', None, None, None, None,
            '.238', '.238', None,
            'NO', '.166', '.238', '.238', None,
            'NO', '.510', '.126', '.126', None,
            'NO', '.565', '.117', '.117', None,
            'NO', '.619', '.126', '.126', None,
            'NO', '.963', '.238', '.238', None,
            'NO', '1.0', '.238', '.238', None,
        ]
        model.add_card(fields, fields[0])

        #fields = [
            #'PBEAML', 625, 624,   None, 'BAR',  None,  None,  None,  None,
                    #0.238, 0.238, None, 'NO',   0.166, 0.238, 0.238, None,
                    #'NO',  0.51,  0.126, 0.126, None,  'NO',  0.565, 0.177,
                    #0.117, None, 'NO',   0.619, 0.126, 0.126, None,  'NO',
                    #0.963, 0.238, 0.238, None,  'NO',  1.0,   0.238, 0.238,
                    #None, 'NO',
        #]
        #model.add_card(fields, fields[0])

        model.pop_parse_errors()
        str(model.properties[622])
        str(model.properties[623])

    def test_cbeam_01(self):
        """modification of test_pbeam_05"""
        model = BDF(debug=False)
        lines = ['PBEAM,3,6,2.9,3.5,5.97,0.4,3.14',
                 '     , , ,2.0,-4.0',]
        model.add_card(lines, 'PBEAM', is_list=False)

        lines = ['CBEAM         10       3      1       2      0.01.000000     0.0']
        model.add_card(lines, 'CBEAM', is_list=False)

        mid = 6
        E = 1.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
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

    def test_cbeam_g0(self):
        """modification of test_cbeam_01"""
        model = BDF(debug=False)
        lines = ['PBEAM,200,6,2.9,3.5,5.97,0.4,3.14',
                 '     , , ,2.0,-4.0',]
        model.add_card(lines, 'PBEAM', is_list=False)

        eid = 100
        pid = 200
        nids = [10, 20]
        x = None
        g0 = 30
        model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None, pa=0,
                        pb=0, wa=None, wb=None, sa=0,
                        sb=0, comment='')

        mid = 6
        E = 1.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(20, [1., 0., 0.])
        model.add_grid(30, [0., 1., 0.])
        model.cross_reference()

        save_load_deck(model)

    def test_cbeam_pbeaml(self):
        """CBEAM/PBEAML"""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0., 1., 0.])

        mid = 1
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.)

        #---------------------------------------------------------------
        eid = 1
        pid = 101
        nids = [1, 2]
        x = [0., 0., 1.]
        g0 = None
        unused_cbeam1 = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                        pa=0, pb=0, wa=None, wb=None,
                                        sa=0, sb=0, comment='CBEAM')
        beam_type = 'BOX'
        xxb = [0.]
        dims = [[1., 2., 0.1, 0.1], [1., 2., 0.1, 0.1]]
        pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, nsm=None,
                                  so=None, comment='PBEAML')
        #---------------------------------------------------------------
        eid = 2
        pid = 102
        x = None
        g0 = 3
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                pa=0, pb=0, wa=None, wb=None,
                                sa=0, sb=0, comment='CBEAM')
        beam_type = 'BOX'
        xxb = [0.]
        dims = [[1., 2., 0.1, 0.1]]
        pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, nsm=None,
                                  so=None, comment='PBEAML')
        #---------------------------------------------------------------
        eid = 3
        pid = 103
        x = [0., 0., 1.]
        g0 = None
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                pa=0, pb=0, wa=None, wb=None,
                                sa=0, sb=0, comment='CBEAM')
        beam_type = 'BAR'
        xxb = [0.]
        dims = [[2., 3.]]
        pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, nsm=None,
                                  so=None, comment='PBEAML')
        #---------------------------------------------------------------
        eid = 4
        pid = 104
        x = [0., 0., 1.]
        g0 = None
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                pa=0, pb=0, wa=None, wb=None,
                                sa=0, sb=0, comment='CBEAM')
        beam_type = 'ROD'
        xxb = [0.]
        dims = [[2.]]
        pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, nsm=None,
                                  so=None, comment='PBEAML')
        #---------------------------------------------------------------
        eid = 5
        pid = 105
        x = [0., 0., 1.]
        g0 = None
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                pa=0, pb=0, wa=None, wb=None,
                                sa=0, sb=0, comment='CBEAM')
        beam_type = 'TUBE'
        xxb = [0.]
        dims = [[2., 1.]]
        pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, nsm=None,
                                  so=None, comment='PBEAML')
        #---------------------------------------------------------------
        eid = 6
        pid = 106
        x = [0., 0., 1.]
        g0 = None
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                pa=0, pb=0, wa=None, wb=None,
                                sa=0, sb=0, comment='CBEAM')
        beam_type = 'L'
        xxb = [0.]
        dims = [[2., 2., 0.1, 0.1]]
        pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, nsm=None,
                                  so=None, comment='PBEAML')
        #---------------------------------------------------------------
        eid = 7
        pid = 107
        x = [0., 0., 1.]
        g0 = None
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                pa=0, pb=0, wa=None, wb=None,
                                sa=0, sb=0, comment='CBEAM')
        beam_type = 'T'
        xxb = [0.]
        dims = [[1., 2., 0.1, 0.2]]
        pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, nsm=None,
                                  so=None, comment='PBEAML')
        #---------------------------------------------------------------
        eid = 8
        pid = 108
        x = [0., 0., 1.]
        g0 = None
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                pa=0, pb=0, wa=None, wb=None,
                                sa=0, sb=0, comment='CBEAM')
        beam_type = 'T1'
        xxb = [0.]
        dims = [[1., 2., 0.2, 0.1]]
        pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, nsm=None,
                                  so=None, comment='PBEAML')
        #---------------------------------------------------------------
        eid = 9
        pid = 109
        x = [0., 0., 1.]
        g0 = None
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                pa=0, pb=0, wa=None, wb=None,
                                sa=0, sb=0, comment='CBEAM')
        beam_type = 'T2'
        xxb = [0.]
        dims = [[1., 2., 0.2, 0.1]]
        pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, nsm=None,
                                  so=None, comment='PBEAML')
        #---------------------------------------------------------------
        eid = 10
        pid = 110
        x = [0., 0., 1.]
        g0 = None
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                pa=42, pb=5, wa=None, wb=None,
                                sa=0, sb=0, comment='CBEAM')
        assert cbeam.is_offt is True, cbeam.is_offt
        assert cbeam.is_bit is False, cbeam.is_bit

        #pid : int
            #property id
        #mid : int
            #material id
        #xxb : List[float]
            #The percentage locations along the beam [0., ..., 1.]
        #so : List[str]
            #YES, YESA, NO
        #area : List[float]
            #area
        #i1, i2, i12, j : List[float]
            #moments of inertia
        #nsm : List[float]
            #nonstructural mass per unit length
        #c1/c2, d1/d2, e1/e2, f1/f2 : List[float]
           #the y/z locations of the stress recovery points
           #c1 - point C.y
           #c2 - point C.z
        xxb = [0.]
        so = [2.]
        area = [3.]
        i1 = [4.]
        i2 = [5.]
        i12 = [0.6]
        j = [7.]
        nsm = [0.08]

        c1 = [0.]
        c2 = [0.]
        d1 = [0.]
        d2 = [0.]
        e1 = [0.]
        e2 = [0.]
        f1 = [0.]
        f2 = [0.]
        unused_pbeam = model.add_pbeam(
            pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
            c1, c2, d1, d2,
            e1, e2, f1, f2,
            k1=1., k2=1., s1=0., s2=0.,
            nsia=0., nsib=None, cwa=0., cwb=None,
            m1a=0., m2a=0., m1b=None, m2b=None,
            n1a=0., n2a=0., n1b=None, n2b=None,
            #m1a=0., m2a=None, m1b=0., m2b=None,
            #n1a=0., n2a=None, n1b=0., n2b=None,
            comment='pbeam')
        #print(pbeam)
        #---------------------------------------------------------------

        model.validate()
        model.pop_parse_errors()
        model.pop_xref_errors()
        pbeaml.raw_fields()
        pbeaml.write_card_16()
        model._verify_bdf(xref=False)

        #----------------
        model.cross_reference()

        model._verify_bdf(xref=True)
        model.mass_properties()
        pids_to_area = model.get_area_breakdown(property_ids=None, sum_bar_area=False)
        for pid, area in sorted(pids_to_area.items()):
            assert area > 0., 'pid=%s area=%s' % (pid, area)

        pids_to_mass, _mass_type_to_mass = model.get_mass_breakdown(property_ids=None)
        for pid, mass in sorted(pids_to_mass.items()):
            assert mass > 0., 'pid=%s mass=%s' % (pid, mass)

        #--------------------
        model.uncross_reference()
        model._verify_bdf(xref=False)
        model.safe_cross_reference()
        model.uncross_reference()

    def test_cbeam_bit(self):
        """tests an BIT field on the CBEAM"""
        model = BDF(debug=False)

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])

        eid = 1
        pid = 1
        bit = 42
        nids = [1, 2]
        x = [0., 0., 1.]
        g0 = None
        bit = None
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt=None, bit=bit, pa=0,
                                pb=0, wa=None, wb=None, sa=0,
                                sb=0, comment='')
        with self.assertRaises(RuntimeError):
            cbeam.validate()
        del model.elements[eid]

        bit = 12
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt=None, bit=bit, pa=0,
                                pb=0, wa=None, wb=None, sa=0,
                                sb=0, comment='')
        with self.assertRaises(AssertionError):
            cbeam.raw_fields()
        del model.elements[eid]

        bit = 42.
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt=None, bit=bit, pa=0,
                                pb=0, wa=None, wb=None, sa=0,
                                sb=0, comment='')
        assert cbeam.is_offt is False, cbeam.is_offt
        assert cbeam.is_bit is True, cbeam.is_bit
        mid = 1
        beam_type = 'BAR'
        dim = [1., 2.]  # area = 2.0
        nsm = 1.
        xxb = [0., 1.]
        dims = [dim, dim]
        model.add_pbeaml(pid, mid, beam_type, xxb, dims, so=None,
                         nsm=[1.0],
                         group='MSCBML0',
                         comment='')

        pid += 1

        so = ['YES', 'YESA', 'NO', 'YES']
        xxb = [0., 0.5, 0.75, 1.]
        area = [1., 2., 3., 4.]
        i1 = [2., 3., 4., 5.]
        i2 = [3., 4., .5, .6]
        i12 = [0.1, 0.2, 0.3, 0.4]
        j = [5., 6., 6.5, 7.0]
        nsm = None
        pbeam = model.add_pbeam(pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                                c1=None, c2=None, d1=None, d2=None,
                                e1=None, e2=None, f1=None, f2=None,
                                k1=1., k2=1., s1=0., s2=0.,
                                nsia=0., nsib=None, cwa=0., cwb=None,
                                m1a=0., m2a=0., m1b=None, m2b=None,
                                n1a=0., n2a=0., n1b=None, n2b=None,
                                comment='')
        str(pbeam)

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.)
        save_load_deck(model, run_mass_properties=False, run_test_bdf=False)


    def test_beam_mass_01(self):
        """tests a CBEAM/PBEAM and gets the mass_properties"""
        model = BDF(debug=False)
        #model.case_control_deck = CaseControlDeck(case_control_lines)
        spc = ['SPC1', 123456, 123456, 1]
        grid1 = ['GRID', 1, None, 0., 0., 0.]
        grid2 = ['GRID', 2, None, 1., 0., 0.]
        #grid3 = ['GRID', 3, None, 1., 0., 0.]
        force = ['FORCE', 100, 1, 0, 2., 3., 4.]

        pid = 11
        cbeam = [
            'CBEAM', 10, pid, 1, 2, 0., 1., 0., None,
        ]
        mid = 12
        nsm_offset_a = [0., 0., 0.]
        nsm_offset_b = [0., 50., 0.]
        k1 = k2 = None
        area1 = 2.0
        area2 = 1.0
        rho = 3.
        nsm_a = 0.
        nsm_b = 0.
        nu = 0.3
        pbeam = ([
            'PBEAM', pid, mid,
            area1, 2.1, 1.2, 1.3, None, nsm_a,
            None, None, None, None, None, None, None, None,

            'YES', 1.0, area2, 2.1, 1.2, 1.3, None, nsm_b,
            None, None, None, None, None, None, None, None,

            # 100s are the NSIa, NSIb (moment of inertia per unit length)
            k1, k2, None, None, 100., 100., None, None] +

                 # Nones are neutral axis offset
                 nsm_offset_a + nsm_offset_b + [None, None, None, None])
        #print('\n' + print_card_8(pbeam))

        mat1 = ['MAT1', 12, 3.0e7, None, nu, rho]
        model.add_card(grid1, 'GRID')
        model.add_card(grid2, 'GRID')
        #model.add_card(grid3, 'GRID')
        model.add_card(cbeam, 'CBEAM')
        model.add_card(pbeam, 'PBEAM')
        model.add_card(mat1, 'MAT1')
        model.add_card(spc, 'SPC1')
        model.add_card(force, 'FORCE')
        model.cross_reference()
        #print(model.properties[11])

        unused_mass, unused_cg, unused_I = model.mass_properties(
            element_ids=None, mass_ids=None,
            reference_point=None,
            sym_axis=None,
            scale=None)
        #print('cg* =', cg)
        L = 1.0
        area = (area1 + area2) / 2.
        nsm = (nsm_a + nsm_b) / 2.
        mass_per_length = area * rho + nsm
        unused_mass = L * mass_per_length

        mass_a = L / 2. * (area1 * rho + nsm_a)
        mass_b = L / 2. * (area2 * rho + nsm_b)
        #xcg = (0.0 * mass_a + 1.0 * mass_b) / (mass_a + mass_b)

        #print(mass_a, mass_b, xcg, mass_a + mass_b)
        #print('mass =', mass)
        #cbeam = CBEAM()
        cbeam = model.elements[10]
        pbeam = model.properties[11]
        cg = cbeam.center_of_mass()
        centroid = cbeam.Centroid()
        assert pbeam.Nu() == nu, 'pbeam.Nu()=%s nu=%s' % (pbeam.Nu(), nu)
        assert pbeam.Rho() == rho, 'pbeam.Rho()=%s rho=%s' % (pbeam.Rho(), rho)
        assert np.allclose(cbeam.Length(), 1.0), cbeam.Length()
        assert np.allclose(cg[0], 0.5), 'cg=%s' % cg
        assert np.allclose(centroid[0], 0.5), 'centroid=%s' % centroid
        #assert allclose(cbeam.Mass(), 10.25), cbeam.Mass()
        #assert allclose(cbeam.MassPerLength(), 10.25), cbeam.MassPerLength()
        #assert allclose(mass, 10.25), mass

        with open('pbeam12.bdf', 'w') as bdf_file:
            case_control_lines = (
                'SOL 101\n'
                'CEND\n'
                'SUBCASE 1\n'
                '    STRESS(PLOT,SORT1,REAL) = ALL\n'
                '    SPC = 123456\n'
                '    LOAD = 100\n'
                'BEGIN BULK\n'
                'PARAM,GRDPNT,0\n'
                'PARAM,POST,-1\n'
                'PARAM   POSTEXT YES\n'
            )
            bdf_file.write(case_control_lines)
            model.write_bdf(bdf_file, enddata=True)

        model2 = BDF(debug=False)
        model2.read_bdf('pbeam12.bdf')

        #if not os.path.exists('pbeam12.op2') and 0:  # pragma: no cover
            #os.system('nastran scr=yes bat=no old=no pbeam12.bdf')
        os.remove('pbeam12.bdf')

        #if 0:  # pragma: no cover
            #from pyNastran.op2.op2 import OP2
            #op2 = OP2()
            #op2.read_op2('pbeam12.op2')
            ##os.remove('pbeam12.op2')
            #gpw = op2.grid_point_weight
            #op2_mass = gpw.mass.max()
            #assert op2_mass == mass, 'op2_mass=%s mass=%s' % (op2_mass, mass)
            #print('op2_mass=%s mass=%s' % (op2_mass, mass))
            #op2_cg = gpw.cg
            #cg = array([0.5, 0., 0.], dtype='float32')
            #print('cg =', op2_cg)

    def test_pbeam_nsm(self):
        """tests a PBEAM with nonstructural mass"""
        model = BDF(debug=False)
        pid = 1
        mid = 1
        nsm = [1., 1.]
        xxb = [0., 1.]
        so = ['YES', 'YES']

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        x = [0., 1., 0.]
        g0 = None

        area = [2.0, 2.0]
        i1 = i2 = j = [0.1, 0.1]
        i12 = [0.01, 0.01]
        c1 = c2 = d1 = d2 = e1 = e2 = f1 = f2 = [0., 0.]
        pbeam_a1 = model.add_pbeam(pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                                   c1, c2, d1, d2,
                                   e1, e2, f1, f2,
                                   k1=1., k2=1.,
                                   s1=0., s2=0.,
                                   nsia=0., nsib=None,
                                   cwa=0., cwb=None,
                                   m1a=0., m2a=0., m1b=None, m2b=None,
                                   n1a=0., n2a=0., n1b=None, n2b=None,
                                   comment='')

        nsm = [1., 1.]
        pid_pbeam_nsm = 30
        pbeam_b1 = model.add_pbeam(pid_pbeam_nsm, mid, xxb, so, area, i1, i2, i12, j, nsm,
                                   c1, c2, d1, d2,
                                   e1, e2, f1, f2,
                                   k1=1., k2=1.,
                                   s1=0., s2=0.,
                                   nsia=10., nsib=10.,
                                   cwa=0., cwb=None,
                                   # cg location at A/B (1.,1.)
                                   m1a=1., m2a=1.,
                                   m1b=None, m2b=None,
                                   # neutral axis at A/B (0., 0.)
                                   n1a=0., n2a=0.,
                                   n1b=None, n2b=None,
                                   comment='')
        eid = 42
        model.add_cbeam(eid, pid_pbeam_nsm, [1, 2], x, g0, offt='GGG', bit=None,
                        pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment='')

        pbeam_a1.validate()
        pbeam_b1.validate()

        E = 1.0
        G = None
        nu = 0.3
        mat1 = model.add_mat1(mid, E, G, nu)
        #----------------
        card_lines = [
            'PBEAM          2       1             BAR',
            '              1.      2.      1.',
        ]
        model.add_card(card_lines, 'PBEAML', comment='', is_list=False,
                       has_none=True)
        model.pop_parse_errors()
        pbeam_a2 = model.properties[2]
        #------------------
        model.cross_reference()
        model.pop_xref_errors()

        assert pbeam_a1.Nsm() == 1.0
        assert pbeam_a1.Area() == 2.0

        # mass/L = area*rho + nsm
        assert pbeam_a1.MassPerLength() == 1.0
        assert pbeam_b1.MassPerLength() == 1.0, pbeam_b1.MassPerLength() # should be 10

        mass, cg1, unused_inertia = model.mass_properties(
            element_ids=eid,
            mass_ids=None,
            reference_point=None,
            sym_axis=None,
            scale=None)
        assert mass == 1.0, mass
        #print('cg1=%s' % cg)
        assert np.allclose(cg1, [0.5, 1., 1.]), cg1

        mass, cg2, unused_inertia = model.mass_properties_nsm(
            element_ids=eid, mass_ids=None,
            nsm_id=None,
            reference_point=None,
            sym_axis=None,
            scale=None,
            xyz_cid0_dict=None,
            debug=False)
        #print('mass =', mass)
        #print('cg2=%s' % cg2)
        assert mass == 1.0, mass
        assert np.allclose(cg2, [0.5, 1., 1.]), cg2

        # area = 2.0
        mat1.rho = 10.0
        assert pbeam_a1.MassPerLength() == 21.0, pbeam_a1.MassPerLength()
        assert pbeam_a2.MassPerLength() == 21.0, pbeam_a2.MassPerLength()

        save_load_deck(model)

    def test_pbeaml_nsm(self):
        """tests the PBEAML with non structural mass"""
        model = BDF(debug=False)
        pid = 1
        mid = 1
        beam_type = 'BAR'
        dim = [1., 2.]  # area = 2.0
        unused_nsm = 1.
        xxb = [0., 1.]
        dims = [dim, dim]
        pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, so=None,
                                  nsm=[1.0],
                                  group='MSCBML0',
                                  comment='')
        E = 1.0
        G = None
        nu = 0.3
        mat1 = model.add_mat1(mid, E, G, nu)
        #----------------
        card_lines = [
            'PBEAML         2       1             BAR',
            '              1.      2.      1.',
        ]
        model.add_card(card_lines, 'PBEAML', comment='', is_list=False,
                       has_none=True)
        model.pop_parse_errors()
        pbeaml2 = model.properties[2]
        #------------------
        model.cross_reference()

        assert pbeaml.Nsm() == 1.0
        assert pbeaml.Area() == 2.0

        # mass/L = area*rho + nsm
        assert pbeaml.MassPerLength() == 1.0

        # area = 2.0
        mat1.rho = 10.0
        assert pbeaml.MassPerLength() == 21.0, pbeaml.MassPerLength()
        assert pbeaml2.MassPerLength() == 21.0, pbeaml2.MassPerLength()

        eid = 1
        nids = [1, 2]
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0., 0., 1.])
        x = None
        g0 = 3
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG',
                                bit=None, pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment='')
        model.uncross_reference()
        model.cross_reference()

        cg = cbeam.center_of_mass()
        centroid = cbeam.Centroid()
        assert np.allclose(cg, [0.5, 0., 0.]), cg
        assert np.allclose(centroid, [0.5, 0., 0.]), centroid
        save_load_deck(model)

    def test_pbeam_opt(self):
        """tests a PBEAM with DVPREL1"""
        model = BDF(debug=False)
        pid = 1
        mid = 1
        nsm = [1., 1.]
        xxb = [0., 1.]
        so = ['YES', 'YES']

        area = [2.0, 2.0]
        i1 = i2 = j = [0.1, 0.1]
        i12 = [0.01, 0.01]
        c1 = c2 = d1 = d2 = e1 = e2 = f1 = f2 = [0., 0.]
        pbeam = model.add_pbeam(pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                                c1, c2, d1, d2,
                                e1, e2, f1, f2,
                                k1=1., k2=1.,
                                s1=0., s2=0.,
                                nsia=0., nsib=None,
                                cwa=0., cwb=None,
                                m1a=0., m2a=0., m1b=None, m2b=None,
                                n1a=0., n2a=0., n1b=None, n2b=None,
                                comment='')
        pbeam.validate()

        E = 1.0
        G = None
        nu = 0.3
        unused_mat1 = model.add_mat1(mid, E, G, nu)
        prop_type = 'PBEAM'

        desvar_id = 1
        xinit = 1.0
        label = 'VAR1'
        model.add_desvar(desvar_id, label, xinit, xlb=-1e20, xub=1e20,
                         delx=None, ddval=None,
                         comment='')
        desvar_id = 2
        label = 'VAR2'
        model.add_desvar(desvar_id, label, xinit, xlb=-1e20, xub=1e20,
                         delx=None, ddval=None,
                         comment='')
        for i in range(15):
            oid = i + 1
            pname_fid = -8 - i
            if pname_fid in [-11, -12, -13, -22, -23]:
                continue
            dvids = desvar_id
            coeffs = 1.0
            model.add_dvprel1(oid, prop_type, pid, pname_fid, dvids, coeffs,
                              p_min=None, p_max=1e20,
                              c0=0.0, validate=True,
                              comment='')
        oid = 20
        dependent_desvar = 1
        independent_desvars = 2
        coeffs = 1.0
        model.add_dlink(oid, dependent_desvar, independent_desvars, coeffs,
                        c0=0., cmult=1., comment='')
        model.cross_reference()
        model.update_model_by_desvars()
        save_load_deck(model, run_save_load_hdf5=True)

    def test_pbcomp(self):
        """tests a PBCOMP"""
        model = BDF(debug=False)
        pid = 100
        mid = 101
        y = [0., 1.]
        z = [0., 1.]
        c = [0., 1.]
        mids = [mid, mid]
        pbcomp = model.add_pbcomp(pid, mid, y, z, c, mids, area=0.0,
                                  i1=0.0, i2=0.0, i12=0.0, j=0.0, nsm=0.0,
                                  k1=1.0, k2=1.0, m1=0.0, m2=0.0, n1=0.0, n2=0.0,
                                  symopt=0, comment='pbcomp')
        pbcomp.validate()
        pbcomp.raw_fields()
        pbcomp.write_card()
        pbcomp.write_card(size=16)

        eid = 10
        nids = [1, 5]
        x = [0., 0., 1.]
        g0 = None
        cbeam = model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                                pa=0, pb=0, wa=None, wb=None, sa=0, sb=0,
                                comment='')
        cbeam.raw_fields()
        cbeam.write_card_16()

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(5, [1., 0., 0.])
        model.validate()
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()
        save_load_deck(model, run_convert=False, run_renumber=False)

    def test_beamor(self):
        """tests a BEAMOR"""
        model = BDF(debug=False)
        n1 = 10
        n2 = 20
        model.add_grid(n1, [0., 0., 0.])
        model.add_grid(n2, [1., 0., 0.])

        pid = 2
        mid = 1
        beam_type = 'BAR'
        dim = [1., 2.]  # area = 2.0
        unused_nsm = 1.
        xxb = [0., 1.]
        dims = [dim, dim]
        unused_pbeaml = model.add_pbeaml(pid, mid, beam_type, xxb, dims, so=None,
                                         nsm=[1.0],
                                         group='MSCBML0',
                                         comment='')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.)

        card_lines = ['BEAMOR', None, pid, None, None, 0.6, 2.9, -5.87, 'GOG']
        model.add_card(card_lines, 'BEAMOR', comment='BEAMOR', is_list=True,
                       has_none=True)

        eid = 1
        card_lines = ['CBEAM', eid, pid, n1, n2]
        model.add_card(card_lines, 'CBEAM', comment='', is_list=True,
                       has_none=True)
        model.pop_parse_errors()

    def test_pbeam3(self):
        """tests a PBEAM3"""
        model = BDF(debug=False)
        pid = 10
        mid = 11
        A = 1.0
        iz = 2.0
        iy = 3.0
        iyz = 0.4
        j = 5.0
        #form = 'GS'
        #options = {'OUTP' : 2}
        #pbrsect = model.add_pbrsect(pid, mid, form, options, comment='pbrsect')
        pbeam3 = model.add_pbeam3(pid, mid, A, iz, iy, iyz, j, nsm=0.,
                                  cy=0., cz=0., dy=0., dz=0., ey=0., ez=0., fy=0., fz=0.,
                                  comment='pbeam3')
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.2)

        pbeam3.validate()
        #pbeam3.raw_fields()
        pbeam3.write_card()
        pbeam3.write_card(size=16)
        model.cross_reference()
        model.pop_xref_errors()

        pbeam3.write_card()
        pbeam3.write_card(size=16)

        model.uncross_reference()
        pbeam3.write_card()
        pbeam3.write_card(size=16)

    def test_cbend(self):
        model = BDF(debug=False)
        model.add_grid(10, [1., 0., 0.])
        model.add_grid(11, [.707, .707, 0.])
        model.add_grid(12, [0., 1., 0.])
        eid = 2
        pid = 3
        nids = [10, 11]
        g0 = 12
        geom = 1

        mid = 100

        # The center of curvature lies on the line AO (or its
        # extension) or vector v.
        #
        #      0
        #     /
        #    +
        #   / \
        #  A   R   B
        #    ------
        #      --   <---- curve A-B
        #
        x = None
        cbend = model.add_cbend(eid, pid, nids, g0, x, geom, comment='cbend')

        A = 1.0
        i1 = 2.0
        i2 = 3.0
        j = 4.0
        pbend1 = PBEND.add_beam_type_1(
            pid, mid, A, i1, i2, j,
            rb=None, theta_b=None,
            c1=0., c2=0., d1=0., d2=0., e1=0., e2=0., f1=0., f2=0.,
            k1=None, k2=None, nsm=0.,
            rc=0., zc=0., delta_n=0., comment='pbend1')

        fsi = 1
        rm = 0.2
        t = 0.4
        #p = 0.5
        #rb = 0.6
        #theta_b = 0.7
        pbend2 = PBEND.add_beam_type_2(
            pid, mid, fsi, rm, t,
            p=None, rb=None, theta_b=None,
            nsm=0., rc=0., zc=0., comment='')
        model.properties[pid] = pbend1

        #model.add_pbend(pid, mid, beam_type, A, i1, i2, j,
                        #c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                        #nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b, comment='')
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.2)
        model.validate()
        model.cross_reference()
        model.pop_xref_errors()

        model.uncross_reference()

        cbend.write_card()
        cbend.write_card(size=16)

        pbend1.write_card()
        pbend1.write_card(size=16)

        pbend2.write_card()
        pbend2.write_card(size=16)

        save_load_deck(model, punch=True, run_convert=False)

    def test_pbrsect(self):
        """tests a PBRSECT"""
        model = BDF(debug=False)
        pid = 2
        mid = 3
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.2)

        outp_id = 10
        ids = [1, 2]
        model.add_point(1, [0., 0., 0.], cp=0, comment='point')
        model.add_point(2, [1., 0., 0.], cp=0, comment='point')
        model.add_set1(outp_id, ids, is_skin=False, comment='')

        form = 'GS'
        options = []
        pbrsect = model.add_pbrsect(pid, mid, form, options, comment='pbrsect')
        with self.assertRaises(AssertionError):
            pbrsect.validate()

        #------------------------------------
        form = 'GS'
        options = [
            ['OUTP', 10],
            ['INP', 20],
        ]
        pbrsect = model.add_pbrsect(pid, mid, form, options, comment='pbrsect')
        pbrsect.validate()
        pbrsect.write_card()

        #------------------------------------
        card = [
            'PBRSECT 4       3       GS',
            '        OUTP=10,INP=20',
        ]
        PBRSECT.add_card(card, comment='')

        #------------------------------------
        card = [
            'PBRSECT 4\t3\tGS',
            '        OUTP=10,INP=20',
        ]
        PBRSECT.add_card(card, comment='')

        #------------------------------------
        card = [
            'PBRSECT 4       3       GS',
            '        OUTP=10,BRP=20,T=1.0,T(11)=1.2, NSM=0.01',
        ]
        pbrsect = PBRSECT.add_card(card, comment='')
        #pbrsect.cross_reference(model)

        #------------------------------------
        card = [
            'PBRSECT 5       3       CP',
            '        OUTP=10,BRP=20,T=1.0,T(11)=[1.2,PT=(123,204)], NSM=0.01',
        ]
        pbrsect = PBRSECT.add_card(card, comment='')
        pbrsect.write_card()

    def test_pbrsect_1(self):
        model = BDF(debug=False, log=None, mode='msc')
        mid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.2)

        #POINT   201             0.5     0.5
        #POINT   202             10.5    0.5
        #POINT   203             19.5    0.5
        #POINT   204             19.5    18.9
        #POINT   205             10.5    18.9
        #POINT   206             0.5     18.9
        #POINT   224             10.5    10.0
        model.add_point(201, [0.5, 0.5, 0.])
        model.add_point(202, [10.5, 0.5, 0.])
        model.add_point(203, [19.5, 0.5, 0.])
        model.add_point(204, [19.5, 18.9, 0.])
        model.add_point(205, [10.5, 18.9, 0.])
        model.add_point(206, [0.5, 18.9, 0.])
        model.add_point(224, [10.5, 10.0, 0.])

        #SET1    201     201     THRU    206
        #SET1    202     202     224     205
        model.add_set1(201, [201, 202, 203, 204, 205, 206])
        set_card = model.add_set1(202, [202, 224, 205])
        assert set_card.ids == [202, 224, 205]

        pid = 32
        form = 'CP'
        options = [
            ['OUTP', 201],
            ['T', 1.0],
            ['BRP', 202],
            ['T(11)', '[2.0, PT=(202,224)]'],
            ['T(12)', '[4.0, PT=(224,205)]'],
        ]
        model.add_pbrsect(pid, mid, form, options)
        model.cross_reference()
        #PBRSECT 32      10      CP
            #OUTP=201,T=1.0,BRP=202,
            #T(11)=[2.0,PT=(202,224)],T(12)=[1.2,PT=(224,205)]

        #model2 = BDF(debug=True, log=None, mode='msc')
        #model2.read_bdf(r'C:\MSC.Software\MSC.Nastran\msc20051\nast\tpl\zbr3.dat')
        if IS_MATPLOTLIB:
            for pid, prop in model.properties.items():
                prop.plot(model, show=False)
            plt.close()

    def test_pbrsect_2(self):
        model = BDF(debug=False)
        pid = 2
        mid = 3
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.2)
        model.add_point(1001, [0., 0., 0.])
        model.add_point(1002, [22.5, 0., 0.])
        model.add_point(1003, [22.5, 45., 0.])
        model.add_point(1004, [0., 45., 0.])
        model.add_point(1005, [45., 0., 0.])
        model.add_point(1006, [45., 45., 0.])
        #point, 1001,,   0.,   0.
        #point, 1002,,  22.5,  0.
        #point, 1003,,  22.5, 45.
        #point, 1004,,   0.,  45.
        #point, 1005,,  45.,   0.
        #point, 1006,,  45.,  45.
        model.add_set1(1000, [1001, 1002, 1003, 1004])
        model.add_set1(2000, [1002, 1005])
        model.add_set1(3000, [1003, 1006])
        set3 = model.add_set3(3001, 'POINT', [1006, 1001])
        assert set3.ids == [1006, 1001]
        #set3, 1000, point, 1001, 1002, 1003, 1004
        #set3, 2000, point, 1002, 1005
        #set3, 3000, point, 1003, 1006
        options = [
            ['OUTP', 1000],
            ['BRP(1)', 2000],
            ['BRP(2)', 3000],
            ['T(1)', 5.],
        ]
        prop = model.add_pbmsect(pid, mid, 'OP', options)

        #pbmsect, 100, 200, op
            #outp=1000, brp(1)=2000, brp(2)=3000
            #t(1)=5.
        model.cross_reference()
        if IS_MATPLOTLIB:
            prop.plot(model, show=False)
            plt.close()

    def test_pbmsect_2(self):
        """
        PBMSECT bulk data entry is utilized to describe the shape of I section
        and PARAM,ARBMSTYP is used to control the selection of formulation.
        Note that default value for PARAM,ARBMSTYP select VKI formulation
        to compute sectional properties of arbitrary cross section with
        isotropic material. However, PARAM,ARBMSTYP,TIMISHEN must be present
        in the bulk data section if PBMSECT entry with Core and/or Layer
        keywords exists in the deck.

        ftp://ftp.uni-siegen.de/lokales/div/MSC/mdnastran_r3_doc_application.pdf
        """
        model = BDF(debug=False)
        pid = 2
        mid = 3
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.2)
        # to select VAM
        #PARAM,ARBMSTYPE,TIMOSHEN

        #$ Composite case
        #PBMSECT 32 1 OP 0.015
        #OUTP=101,C=101,brp=103,c(1)=[201,pt=(15,34)]
        #---------------------
        #$ Section profile
        model.add_point(1, [-0.5, .23, 0.])
        model.add_point(2, [0.00, .23, 0.])
        model.add_point(3, [0.50, .23, 0.])
        model.add_point(4, [-0.5, -.23, 0.])
        model.add_point(5, [0.00, -.23, 0.])
        model.add_point(6, [0.50, -.23, 0.])
        #point 1 -0.50 0.23
        #point 2 0.00 0.23
        #point 3 0.50 0.23
        #point 4 -0.50 -0.23
        #point 5 0.00 -0.23
        #point 6 0.50 -0.23

        #.......2.......3.......4.......5.......6.......7.......8.......9.......10.....
        model.add_set1(101, [1, 2, 5, 6])
        model.add_set1(201, [2, 3])
        model.add_set1(102, [5, 4])
        #SET1 101 1 2 5 6
        #SET1 201 2 3
        #SET1 102 5 4

        # Ply properties
        #.......2.......3.......4.......5.......6.......7.......8.......9.......10.....
        #$MAT8 501 20.59e6 1.42e6 0.42 0.89e6 0.89e6 0.89e6
        #$MAT1 501 1.+7 .3

        # isotropic case using T keyword
        options = [
            ['OUTP', 101],
            ['T', 0.04],
            ['BRP(1)', 201],
            ['BRP(3)', 102],
        ]
        prop1 = model.add_pbmsect(pid, mid, 'OP', options)
        #PBMSECT 31 1 OP
        #OUTP=101,t=0.04,BRP(1)=201,BRP(3)=102

        # isotropic case using C and MAT1
        options = [
            ['OUTP', 101],
            ['CORE', 301],
            ['CORE(1)', '[101,PT=(1,2)]'],
            ['CORE(2)', '[202,PT=(5,6)]'],
            ['BRP(1)', 201],
            ['CORE(3)', '[201,PT=(2,3)]'],
            ['BRP(3)', 102],
            ['C(3)', '[102,PT=(5,4)]'],
        ]
        prop2 = model.add_pbmsect(pid+1, mid, 'OP', options)
        str(prop1)
        str(prop2)
        model.cross_reference()
        #PBMSECT 32 OP +
        #OUTP=101,CORE=301,CORE(1)=[101,PT=(1,2)],CORE(2)=[202,PT=(5,6)],+
        #BRP(1)=201,CORE(3)=[201,PT=(2,3)], +
        #BRP(3)=102,CORE(3)=[102,PT=(5,4)]

        if IS_MATPLOTLIB:
            prop1.plot(model, figure_id=pid, show=False)
            #plt.close()
            prop2.plot(model, figure_id=pid+1, show=False)
            #plt.close()
        #save_load_deck(model, run_convert=False, run_remove_unused=False)

    def test_pbmsect_1(self):
        model = BDF(debug=False, log=None, mode='msc')
        mid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.2)

        model.add_point(111, [0.0, 0.00, 0.,])
        model.add_point(112, [20., 0.00, 0.,])
        model.add_point(113, [20., 19.4, 0.,])
        model.add_point(114, [0.0, 19.4, 0.,])
        model.add_point(125, [1.0, 1.00, 0.,])
        model.add_point(126, [10., 1.00, 0.,])
        model.add_point(127, [10., 18.4, 0.,])
        model.add_point(128, [1.0, 18.4, 0.,])
        model.add_point(132, [11., 1.00, 0.,])
        model.add_point(133, [19., 1.00, 0.,])
        model.add_point(134, [19., 18.4, 0.,])
        model.add_point(135, [11., 18.4, 0.,])
        #POINT   111             0.0     0.0
        #POINT   112             20.0    0.0
        #POINT   113             20.0    19.4
        #POINT   114             0.0     19.4
        #POINT   125             1.0     1.0
        #POINT   126             10.0    1.0
        #POINT   127             10.0    18.4
        #POINT   128             1.0     18.4
        #POINT   132             11.0    1.0
        #POINT   133             19.0    1.0
        #POINT   134             19.0    18.4
        #POINT   135             11.0    18.4

        model.add_set1(111, [111, 'THRU', 114])
        model.add_set1(122, [128, 127, 126, 125, 128])
        model.add_set1(133, [132, 133, 134, 135])
        #SET1    111     111     THRU    114
        #SET1    122     128     127     126     125     128
        #SET1    133     132     133     134     135

        options = [
            ['OUTP', 111],
            ['INP', 122],
            ['INP', 133],
        ]
        pid = 31
        form = 'GS'
        prop = model.add_pbmsect(pid, mid, form, options)
        model.cross_reference()
        #print(prop)
        if IS_MATPLOTLIB:
            prop.plot(model, figure_id=2, show=False)
            plt.close()
        #pbmsect  31     10      GS
            #OUTP=111,inp=122,inp=133
        #DESVAR  31      DBOXGS	15.0	0.01	60.0
        #DVPREL1	112	pbmsect	31      W
            #31      1.0

    def test_pbmsect(self):
        model = BDF(debug=False)
        pid = 2
        mid = 3
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.2)

        #PBMSECT 32      10      OP
            #OUTP=101,T=0.1,brp=102,brp=103,brp=104,nsm=0.015

        outp_id = 10
        ids = [1, 2]
        model.add_point(1, [0., 0., 0.], cp=0, comment='point')
        model.add_point(2, [1., 0., 0.], cp=0, comment='point')
        model.add_set1(outp_id, ids, is_skin=False, comment='')

        form = 'GS'
        options = [
            ['OUTP', outp_id],
        ]
        pbmsect = model.add_pbmsect(pid, mid, form, options, comment='pbmsect')
        model.cross_reference()
        if IS_MATPLOTLIB:
            #pbmsect.cross_reference(model)
            pbmsect.plot(model, show=False)
            plt.close()


        pbmsect.validate()
        pbmsect.write_card()
        card = [
            'PBMSECT 4       3       GS',
            '        OUTP=10,INP=20',
        ]
        PBMSECT.add_card(card, comment='')

        card = [
            'PBMSECT 4\t3\tGS',
            '        OUTP=10,INP=20',
        ]
        PBMSECT.add_card(card, comment='')

        card = [
            'PBMSECT 4       3       GS',
            '        OUTP=10,BRP=20,T=1.0,T(11)=1.2, NSM=0.01',
        ]
        pbmsect = PBMSECT.add_card(card, comment='')
        # doesn't work because missing SET=20 (from BRP=20)
        #pbmsect.cross_reference(model)

        line = 'OUTP=10,BRP=20,T=1.0,T(11)=[1.2,PT=(123,204)], NSM=0.01'
        sections = get_beam_sections(line)
        msg = '\nsections: %s\n' % sections
        msg += "expected: [u'OUTP=10', u'BRP=20', u'T=1.0', u'T(11)=[1.2,PT=(123,204)', 'NSM=0.01']"
        assert sections == ['OUTP=10', 'BRP=20', 'T=1.0', 'T(11)=[1.2,PT=(123,204)]', 'NSM=0.01'], msg

        card = [
            'PBMSECT 5       3       CP',
            '        OUTP=10,BRP=20,T=1.0,T(11)=[1.2,PT=(123,204)], NSM=0.01',
        ]
        pbmsect = PBMSECT.add_card(card, comment='')
        pbmsect.write_card()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
