import random
import unittest

import numpy as np
from pyNastran.bdf.field_writer import print_card
from pyNastran.bdf.field_writer_8 import (print_field_8, print_float_8,
                                          set_default_if_blank, print_int_card,
                                          set_blank_if_default, is_same, print_card_8,
                                          print_scientific_8,
                                          print_int_default_8, print_float_default_8,
                                          set_string8_blank_if_default)
from pyNastran.bdf.field_writer_16 import (print_field_16, print_card_16, print_float_16,
                                           print_scientific_16, set_string16_blank_if_default)
from pyNastran.bdf.field_writer_double import print_card_double, print_field_double, print_scientific_double


from pyNastran.bdf.bdf_interface.assign_type import interpret_value


class Testfield_writer_8(unittest.TestCase):

    def test_print_int_card(self):
        card = ['SET1', 1, 2, 3, 4, 5, 6, 7]
        outa = print_int_card(card)
        outb = print_card_8(card)
        # print(f'outa = {outa!r}')
        # print(f'outb = {outb!r}')
        expected = 'SET1           1       2       3       4       5       6       7\n'
        assert outa == expected
        assert outb == expected

        card = ['SET1', 1, 2, 3, 4, 5, 6, 7, 8]
        outc = print_int_card(card)
        outd = print_card_8(card)
        # print(f'outc = {outd!r}')
        # print(f'outd = {outc!r}')
        expected = 'SET1           1       2       3       4       5       6       7       8\n'
        assert outc == expected
        assert outd == expected

        card = ['SET1', 1, 2, 3, 4, 5, 6, 7, 8, 9]
        oute = print_int_card(card)
        outf = print_card_8(card)
        # print(f'oute = {oute!r}')
        # print(f'outf = {outf!r}')
        expected = 'SET1           1       2       3       4       5       6       7       8\n               9\n'
        assert oute == expected
        assert outf == expected

        card = ['SET1', 1, 2, 3, 4, 5, 6, 7, 8,
                9, 10, 11, 12, 13, 14, 15]
        assert len(card) % 8 == 0
        outg = print_int_card(card)
        outh = print_card_8(card)
        # print(f'outg = {outc!r}')
        # print(f'outh = {outh!r}')
        expected = (
            'SET1           1       2       3       4       5       6       7       8\n'
            '               9      10      11      12      13      14      15\n')
        assert outg == expected
        assert outh == expected

    def test_field_vals_8_edge_cases(self):
        self.assertEqual(print_field_8(-999999.56), '  -10.+5',
                         print_field_8(-999999.56))

    def test_field_vals_8(self):
        self.assertEqual(print_field_8(1e20), '   1.+20',
                         print_field_8(1e20))
        self.assertEqual(print_field_8(-.723476), '-.723476',
                         print_field_8(-.723476))
        self.assertEqual(print_field_8(125000.), ' 125000.',
                         print_field_8(125000.))
        self.assertEqual(print_field_8(12500000.), '  1.25+7',
                         print_field_8(12500000.))
        self.assertEqual(print_field_8(47.77267), '47.77267',
                         print_field_8(47.77267))
        self.assertEqual(print_field_8(.001), '    .001',
                         print_field_8(.001))
        self.assertEqual(print_field_8(.0000001), '.0000001',
                         print_field_8(.0000001))
        self.assertEqual(print_field_8(-5.007e-3), '-5.007-3',
                         print_field_8(-5.007e-3))
        self.assertEqual(print_field_8(-0.0748662), '-.074866',
                         print_field_8(-0.0748662))

    def test_field_random(self):
        for i in range(100):
            a = random.uniform(-20, 20)
            a2 = 10 ** a
            compare(a2)
            compare(-a2)

    def test_strings_8(self):
        self.assertEqual(print_field_8(None), '        ',
                         print_field_8(None))
        self.assertEqual(print_field_8('asdf'), '    asdf',
                         print_field_8('asdf'))
        self.assertEqual(print_field_8('  asdf  '), '  asdf  ',
                         print_field_8('  asdf  '))
        self.assertRaises(RuntimeError, print_field_8, '  asdf   ')

    def test_strings_16(self):
        self.assertEqual(print_field_16(None), '                ',
                         print_field_16(None))
        self.assertEqual(print_field_16('asdf'), '            asdf',
                         print_field_16('asdf'))
        self.assertEqual(print_field_16('  asdf  '), '          asdf  ',
                         print_field_16('  asdf  '))
        self.assertRaises(RuntimeError, print_field_16, '          asdf   ')

    def test_field_defaults(self):
        self.assertEqual(set_blank_if_default(0.0, 0.0), None,
                         set_blank_if_default(0.0, 0.0))

        self.assertEqual(set_blank_if_default(1.0, 0.0), 1.0,
                         set_blank_if_default(1.0, 0.0))

        # None
        self.assertEqual(set_default_if_blank(None, None), None,
                         set_default_if_blank(None, None))

        # floats
        self.assertEqual(set_default_if_blank(4.0, None), 4.0,
                         set_default_if_blank(4.0, None))
        self.assertEqual(set_default_if_blank(None, 4.0), 4.0,
                         set_default_if_blank(None, 4.0))

        # ints
        self.assertEqual(set_default_if_blank(4, None), 4,
                         set_default_if_blank(4, None))
        self.assertEqual(set_default_if_blank(None, 4), 4,
                         set_default_if_blank(None, 4))

        # strings
        self.assertEqual(set_default_if_blank('dummy', 'GGG'), 'dummy',
                         set_default_if_blank('dummy', 'GGG'))
        self.assertEqual(set_default_if_blank(None, 'GGG'), 'GGG',
                         set_default_if_blank(None, 'GGG'))

        #set_default_if_blank

    def test_ints_8(self):
        self.assertEqual(print_field_8(1), '       1', 'a')
        self.assertEqual(print_field_8(12345678), '12345678', 'b')
        self.assertRaises(RuntimeError, print_field_8, 123456789)
        self.assertEqual(print_field_8('12345678'), '12345678', 'c')
        self.assertEqual(print_field_8('1       '), '1       ',
                         '|%s|' %(print_field_8('1       ')))

    def test_ints_16(self):
        self.assertEqual(print_field_16(1), '               1', 'a')
        self.assertEqual(print_field_16(12345678), '        12345678', 'b')
        self.assertEqual(print_field_16(1234567890123456), '1234567890123456', 'c')
        self.assertRaises(RuntimeError, print_field_16, 12345678901234567)

        #msg = print_field_16('12345678        ')
        #msg = '|%s| len(msg)=%s' %(msg, len(msg))
        #self.assertEqual(print_field_16('12345678'), '12345678        ',msg)
        self.assertEqual(print_field_16('1               '), '1               ',
                         '|%s|' % (print_field_8('1       ')))

    def test_floats_positive_8(self):
        # ideal
        #self.assertEqual(print_float_8(-.003607), '-.003607',
        #                 print_float_8(-.003607))

        # actual
        self.assertEqual(print_float_8(-.003607), '-3.607-3',
                         print_float_8(-.003607))

        # good
        self.assertEqual(print_float_8(1.2), '     1.2',
                         print_float_8(1.2))
        self.assertEqual(print_float_8(0.5), '      .5',
                         print_float_8(0.5))
        self.assertEqual(print_float_8(-0.5), '     -.5',
                         print_float_8(-0.5))

        self.assertEqual(print_field_8(1.2), '     1.2', 'a')
        self.assertEqual(print_field_8(1.23456789), '1.234568', 'b')
        self.assertEqual(print_field_8(12.234568), '12.23457', 'c')
        self.assertEqual(print_field_8(123.23457), '123.2346', 'd')
        self.assertEqual(print_field_8(1234.23468), '1234.235', 'e')
        self.assertEqual(print_field_8(12345.238), '12345.24', 'f')
        self.assertEqual(print_field_8(123456.28), '123456.3', 'g')
        self.assertEqual(print_field_8(1234567.25), '1234567.',
                         print_field_8(1234567.25))  # 1.2346+6
        self.assertEqual(print_field_8(12345678.), '1.2346+7',
                         '%r' % print_field_8(12345678.))
        self.assertEqual(print_field_8(
            123456789.), '1.2346+8', '%r' % print_field_8(12345678.))

        self.assertEqual(print_field_8(0.1), '      .1',
                         'A|%s|' % print_field_8(0.1))
        self.assertEqual(print_field_8(0.0001), '   .0001',
                         'B|%s|' % print_field_8(0.0001))
        self.assertEqual(print_field_8(0.00001), '  .00001',
                         'C|%s|' % print_field_8(0.00001))
        self.assertEqual(print_field_8(0.000001), ' .000001',
                         'D|%s|' % print_field_8(0.000001))
        self.assertEqual(print_field_8(0.0000001), '.0000001',
                         'E|%s|' % print_field_8(0.0000001))
        self.assertEqual(print_field_8(0.00000012), '   1.2-7',
                         'F|%s|' % print_field_8(0.00000012))
        self.assertEqual(print_field_8(0.000748519), '7.4852-4',
                         'G|%s|' % print_field_8(0.000748519))
        self.assertEqual(print_field_8(0.12345678), '.1234568',
                         'H|%s|' % print_field_8(0.12345678))
        self.assertEqual(print_field_8(0.00012349), '1.2349-4',
                         'I|%s|' % print_field_8(0.00012349))
        self.assertEqual(print_field_8(0.000012349), '1.2349-5',
                         'J|%s|' % print_field_8(0.000012349))
        self.assertEqual(print_field_8(5e-08), '    5.-8',
                         'K|%s|' % print_field_8(5e-08))

        self.assertEqual(print_field_8(1e-20), '   1.-20',
                         'L|%s|' % print_field_8(1e-20))

        self.assertEqual(print_field_8(0.0), '      0.',
                         print_field_8(0.0))
        self.assertEqual(print_field_8(1.0), '      1.',
                         print_field_8(1.0))

    def test_floats_negative_8(self):
        self.assertEqual(print_field_8(-1.2), '    -1.2',
                         print_field_8(-1.2))
        self.assertEqual(print_field_8(-1.23456789), '-1.23457',
                         print_field_8(-1.23456789))
        self.assertEqual(print_field_8(-12.234568), '-12.2346',
                         print_field_8(-12.234568))
        self.assertEqual(print_field_8(-123.23457), '-123.235',
                         print_field_8(-123.23457))
        self.assertEqual(print_field_8(-1234.23468), '-1234.23',
                         print_field_8(-1234.23468))
        self.assertEqual(print_field_8(-12345.238), '-12345.2',
                         print_field_8(-12345.238))
        self.assertEqual(print_field_8(-123456.28), '-123456.',
                         print_field_8(-123456.28))
        self.assertEqual(print_field_8(-1234567.25), '-1.235+6',  # is this right?
                         print_field_8(-1234567.25))
        self.assertEqual(print_field_8(-12345678.), '-1.235+7',  # is this right?
                         '|%s|' % print_field_8(-12345678.))
        self.assertEqual(print_field_8(-123456789.), '-1.235+8',  # is this right?
                         '|%s|' % print_field_8(-123456789.))

        self.assertEqual(print_field_8(-0.1), '     -.1',
                         'A|%s|' % print_field_8(-0.1))
        self.assertEqual(print_field_8(-0.0001), '  -.0001',
                         'B|%s|' % print_field_8(-0.0001))
        self.assertEqual(print_field_8(-0.00001), ' -.00001',
                         'C|%s|' % print_field_8(-0.00001))
        self.assertEqual(print_field_8(-0.000001), '   -1.-6',
                         'D|%s|' % print_field_8(-0.000001))
        self.assertEqual(print_field_8(-0.0000001), '   -1.-7',
                         'E|%s|' % print_field_8(-0.0000001))
        self.assertEqual(print_field_8(-0.00000012), '  -1.2-7',
                         'F|%s|' % print_field_8(-0.00000012))
        self.assertEqual(print_field_8(-0.000748519), '-7.485-4',
                         'G|%s|' % print_field_8(-0.000748519))
        self.assertEqual(print_field_8(-0.12345678), '-.123457',
                         'H|%s|' % print_field_8(-0.12345678))
        self.assertEqual(print_field_8(-0.00012349), '-1.235-4',
                         'I|%s|' % print_field_8(-0.00012349))
        self.assertEqual(print_field_8(-0.000012349), '-1.235-5',
                         'J|%s|' % print_field_8(-0.000012349))
        self.assertEqual(print_field_8(-1e-5), ' -.00001',
                         'K|%s|' % print_field_8(-1e-5))
        self.assertEqual(print_field_8(-1e-6), '   -1.-6',
                         'L|%s|' % print_field_8(-1e-6))
        self.assertEqual(print_field_8(-1e-7), '   -1.-7',
                         'M|%s|' % print_field_8(-1e-7))
        self.assertEqual(print_field_8(-1e-20), '  -1.-20',
                         'N|%s|' % print_field_8(-1e-20))

    def test_float_8_edge_cases(self):
        """Test print_float_8 at notation boundaries and precision limits."""
        expected_value = [
            # Boundary between fixed and scientific (positive)
            ('999999.9', 999999.9),
            ('1000000.', 1000000.0),
            ('9999999.', 9999999.0),
            # Boundary between fixed and scientific (negative)
            ('  -10.+5', -999999.9),
            ('   -1.+6', -1000000.0),
            ('  -10.+6', -9999999.0),
            # Near-zero threshold (5e-8 boundary)
            ('    5.-8', 5e-8),
            ('   4.9-8', 4.9e-8),
            ('   -5.-8', -5e-8),
            ('   -5.-7', -5e-7),
            # NaN and zero
            ('        ', float('nan')),
            ('      0.', 0.0),
            ('      0.', -0.0),
            # Powers of 10 (notation transitions)
            ('    .001', 0.001),
            ('   .0001', 0.0001),
            ('     .01', 0.01),
            ('      .1', 0.1),
            ('      1.', 1.0),
            ('     10.', 10.0),
            ('    100.', 100.0),
            ('   1000.', 1000.0),
            ('  10000.', 10000.0),
            (' 100000.', 100000.0),
            ('   -.001', -0.001),
            ('  -.0001', -0.0001),
            ('    -.01', -0.01),
            ('     -.1', -0.1),
            ('     -1.', -1.0),
            ('    -10.', -10.0),
            ('   -100.', -100.0),
            ('  -1000.', -1000.0),
            (' -10000.', -10000.0),
            ('-100000.', -100000.0),
            # Max precision at each magnitude
            ('9.999999', 9.999999),
            ('99.99999', 99.99999),
            # Rounding at field boundary
            ('1234568.', 1234567.5),
            ('-1.235+6', -1234567.5),
            ('  -10.+6', -9999999.5),
            # Rounding overflow: must not lose decimal point
            ('   10.+6', 9999999.5),
            ('   10.+6', 9999999.9),
        ]
        for expected, value in expected_value:
            actual = print_float_8(value)
            self.assertEqual(len(actual), 8,
                             f'len={len(actual)} for value={value}: |{actual}|')
            self.assertEqual(actual, expected,
                             f'value={value}: actual=|{actual}| expected=|{expected}|')

    def test_print_card_8(self):
        self.assertEqual(print_card_8(['GRID', 1]), 'GRID           1\n')
        self.assertEqual(print_card_8(['GRID', 1, None, None, None, None, None,
                                       None, None, 1]),
                         'GRID           1\n               1\n',
                         print_card_8(['GRID', 1, None, None, None, None, None,
                                       None, None, 1]))

        self.assertEqual(print_card_8(['GRID', 1, None, None, None, None, None,
                                       None, None, 1, None]),
                         'GRID           1\n               1\n',
                         print_card_8(['GRID', 1, None, None, None, None, None,
                                       None, None, 1, None]))

        self.assertEqual(print_card_8(['GRID', 1, None, None, None, None, None,
                                       None, None, None, None, None, None, None,
                                       None, None, None, 1, None]),
                         'GRID           1\n+\n               1\n',
                         print_card_8(['GRID', 1, None, None, None, None, None,
                                       None, None, None, None, None, None, None,
                                       None, None, None, 1, None]))

    def test_same(self):
        self.assertTrue(is_same(1.0, 1.000))
        self.assertFalse(is_same(1.0, 1e-15 + 1))
        self.assertTrue(is_same('MAT', 'MAT'))
        self.assertFalse(is_same('MAT', 'MAT1'))

    def test_card_double(self):
        fields = ['GRID', 1, None, 120.322, -4.82872, 1.13362]
        card_expected = (
            'GRID*                  1                1.2032200000D+02-4.828720000D+00\n'
            '*       1.1336200000D+00\n'
        )
        card1 = print_card(fields, size=16, is_double=True)
        card2 = print_card_double(fields)
        self.assertEqual(card1, card_expected)
        self.assertEqual(card2, card_expected)

        #card = print_card_double(['CHEXA', 1, 2, 2, 3, 4, 1, 8, 5,
                                           #6, 7])
        #card_expected = ''
        #card_expected += 'CHEXA*                 1               2               2               3\n'
        #card_expected +='*                      4               1               8               5\n'
        #card_expected +='*                      6               7\n'
        #card_expected += '*       1.1336200000D+00\n'

        #card = print_card_double(['CTETRA',6437,1,14533,5598,1577,9976,42364,5599,42365,42363,30022,12904])
        #card_expected = ''

        #card_expected += 'CTETRA*             6437               1           14533            5598\n'
        #card_expected += '*                   1577            9976           42364            5599\n'
        #card_expected += '*                  42365           42363           30022           12904\n'
        #card_expected += '*       \n'

        #card_expected +='CTETRA*  6437            1               14533           5598\n'
        #card_expected +='*        1577            9976            42364           5599\n'
        #card_expected +='*        42365           42363           30022           12904\n'
        #card_expected +='*       \n'
        #self.assertEqual(card, card_expected)

        #=============================
        #           mid   E1      E2      nu12 G12      G1z      G2z      rho
        fields = [
            'MAT8', 6, 1.7e+7, 1.7e+7, .98, 340000., 180000., 180000., 0.0001712,
            # a1    a2    tref
            None, 71.33]
        card_expected = (
            'MAT8           6   1.7+7   1.7+7     .98 340000. 180000. 180000..0001712\n'
            '                   71.33\n'
        )
        card1 = print_card_8(fields)
        card2 = print_card(fields)
        self.assertEqual(card1, card_expected)
        self.assertEqual(card2, card_expected)

        fields = ['MAT8', 6, 1.7e+7, 1.7e+7, .98, 340000., 180000., 180000., 0.0001712,
                  # a1    a2    tref
                  None, 71.33]
        card1 = print_card_16(fields)
        card2 = print_card(fields, size=16)

        card_expected = (
            'MAT8*                  6       17000000.       17000000.             .98\n'
            '*                340000.         180000.         180000.        .0001712\n'
            '*                                  71.33\n'
            '*\n'
        )
        # bad
        self.assertEqual(card1, card_expected)
        self.assertEqual(card2, card_expected)
        #=============================
        #layer = mid, t, theta sout
        #                      mid E1              E2                  nu12
        #'MAT8*                  6 1.7000000000+07 1.7000000000+07     0.980000000+ '
        #                      G12 G1z             G2z                 rho
        #'*        340000.00000000 180000.00000000 180000.00000000     0.000171204+ '
        #                       a1              a2 Tref                Xt
        #'*                    0.0             0.0    71.330000000                + '
        #                      Xc               Yt Yc                  S
        #'*                                                                       + '
        #                       GE             F12 STRN
        #'*                    0.0             0.0                                + '
        #'*

        #====
        #'MAT8*                  6 1.7000000000+07 1.7000000000+07     0.980000000+ ' - 1
        #'*        340000.00000000 180000.00000000 180000.00000000     0.000171204+ ' - 2
        #'*                    0.0             0.0    71.330000000                + ' - 3
        #'*                                                                       + ' - 4
        #'*                    0.0             0.0                                + ' - 5
        #'*                                                                           - 6
        #print('1')

    def test_float_8(self):
        small_exponent = -17
        large_exponent = 17

        nums = (
            [0., 0.000034, -0.000034] +
            [9./11 * 10**x for x in range(small_exponent, large_exponent+1)] +
            [-9./11 * 10**x for x in range(small_exponent, large_exponent+1)])

        expected = [
            '      0.', ' .000034', '  -3.4-5',

            '8.182-18', '8.182-17', '8.182-16', '8.182-15', '8.182-14', '8.182-13',
            '8.182-12', '8.182-11', '8.182-10', '8.1818-9', '8.1818-8', '8.1818-7',
            '8.1818-6', '8.1818-5', '8.1818-4', '.0081818', '.0818182', '.8181818',
            '8.181818', '81.81818', '818.1818', '8181.818', '81818.18', '818181.8',
            '8181818.', '8.1818+7', '8.1818+8', '8.1818+9', '8.182+10', '8.182+11',
            '8.182+12', '8.182+13', '8.182+14', '8.182+15', '8.182+16',

            '-8.18-18', '-8.18-17', '-8.18-16', '-8.18-15', '-8.18-14', '-8.18-13',
            '-8.18-12', '-8.18-11', '-8.18-10', '-8.182-9', '-8.182-8', '-8.182-7',
            '-8.182-6', '-8.182-5', '-8.182-4', '-8.182-3', '-.081818', '-.818182',
            '-8.18182', '-81.8182', '-818.182', '-8181.82', '-81818.2', '-818182.',
            '-8.182+6', '-8.182+7', '-8.182+8', '-8.182+9', '-8.18+10', '-8.18+11',
            '-8.18+12', '-8.18+13', '-8.18+14', '-8.18+15', '-8.18+16',
        ]
        for num, expectedi in zip(nums, expected):
            output = print_float_8(num)
            self.assertEqual(len(output), 8, msg='output=%r len(output)=%i' % (output, len(output)))
            self.assertEqual(output, expectedi, msg='num=%s output=%r expected=%r' % (num, output, expectedi))

    def test_float_8_many(self):
        for istart in np.arange(-13, 13):
            #print(istart)
            nums = np.logspace(istart, istart+1, num=1000, endpoint=True, base=10.0, dtype=None)
            for num in nums:
                output = print_float_8(num)
                self.assertEqual(len(output), 8, msg='output=%r len(output)=%i' % (output, len(output)))
                #self.assertEqual(output, expectedi, msg='num=%s output=%r expected=%r' % (x, output, expectedi))
                output = print_scientific_8(num)
                self.assertEqual(len(output), 8, msg='output=%r len(output)=%i' % (output, len(output)))

    def test_scientific_8(self):
        expected_num = [
            ('      0.', 0.),
            ('    1.-1', 0.1),
            ('   1.1-1', 0.11),
            ('1.2345-1', 0.123451),
            ('1.2346-1', 0.123459),
            # positive values with positive exponent
            ('   1.5+3', 1.5e3),
            ('  1.25+7', 1.25e7),
            ('   1.+20', 1e20),
            ('    1.+0', 1.0),
            ('    1.+1', 10.0),
            # negative values
            ('  -1.5+3', -1.5e3),
            (' -9.31-4', -9.31e-4),
            ('-5.007-3', -5.007e-3),
            ('  -10.+5', -999999.56),
            ('  -1.+20', -1e20),
            # small positive values
            ('    1.-9', 1e-9),
            ('    5.-8', 5e-8),
            (' 5.007-3', 5.007e-3),
            # large exponents (2+ digits)
            ('1.235-10', 1.23456e-10),
            (' 1.23+10', 1.23e10),
            ('-1.23+10', -1.23e10),
            ('-1.23-10', -1.23e-10),
            # exact powers of 10
            ('    1.+2', 100.0),
            ('    1.-2', 0.01),
            ('    1.-5', 1e-5),
            ('    1.+9', 1e9),
        ]
        for expected, num in expected_num:
            actual = print_scientific_8(num)
            assert len(actual) == 8, f'len={len(actual)} actual={actual!r} num={num}'
            assert actual == expected, f'actual={actual!r} expected={expected!r} num={num}'

    def test_scientific_16_edge_cases(self):
        expected_num = [
            ('              0.', 0.),
            ('           1.5+3', 1.5e3),
            ('          -1.5+3', -1.5e3),
            ('           1.+20', 1e20),
            ('          -1.+20', -1e20),
            ('            1.+0', 1.0),
            ('            1.-1', 0.1),
            ('            1.-9', 1e-9),
            ('       1.2345-10', 1.2345e-10),
            ('      -1.2345-10', -1.2345e-10),
            ('         1.23+10', 1.23e10),
            ('        -1.23+10', -1.23e10),
            ('            1.+2', 100.0),
            ('            1.-2', 0.01),
        ]
        for expected, num in expected_num:
            actual = print_scientific_16(num)
            assert len(actual) == 16, f'len={len(actual)} actual={actual!r} num={num}'
            assert actual == expected, f'actual={actual!r} expected={expected!r} num={num}'

    def test_scientific_16(self):
        small_exponent = -17
        large_exponent = 17

        nums = (
            [0.,
             0.000034, -0.000034] +
            [9./11 * 10**x for x in range(small_exponent, large_exponent+1)] +
            [-9./11 * 10**x for x in range(small_exponent, large_exponent+1)])

        expected = [
            '              0.',
            '           3.4-5', '          -3.4-5',

            '8.18181818182-18', '8.18181818182-17', '8.18181818182-16', '8.18181818182-15',
            '8.18181818182-14', '8.18181818182-13', '8.18181818182-12', '8.18181818182-11',
            '8.18181818182-10', '8.181818181818-9', '8.181818181818-8', '8.181818181818-7',
            '8.181818181818-6', '8.181818181818-5', '8.181818181818-4', '8.181818181818-3',
            '8.181818181818-2', '8.181818181818-1', '8.181818181818+0',

            '8.181818181818+1', '8.181818181818+2', '8.181818181818+3', '8.181818181818+4',
            '8.181818181818+5', '8.181818181818+6', '8.181818181818+7', '8.181818181818+8',
            '8.181818181818+9', '8.18181818182+10', '8.18181818182+11', '8.18181818182+12',
            '8.18181818182+13', '8.18181818182+14', '8.18181818182+15', '8.18181818182+16',

            '-8.1818181818-18', '-8.1818181818-17', '-8.1818181818-16', '-8.1818181818-15',
            '-8.1818181818-14', '-8.1818181818-13', '-8.1818181818-12', '-8.1818181818-11',
            '-8.1818181818-10', '-8.18181818182-9', '-8.18181818182-8', '-8.18181818182-7',
            '-8.18181818182-6', '-8.18181818182-5', '-8.18181818182-4', '-8.18181818182-3',
            '-8.18181818182-2', '-8.18181818182-1', '-8.18181818182+0',

            '-8.18181818182+1', '-8.18181818182+2', '-8.18181818182+3', '-8.18181818182+4',
            '-8.18181818182+5', '-8.18181818182+6', '-8.18181818182+7', '-8.18181818182+8',
            '-8.18181818182+9', '-8.1818181818+10', '-8.1818181818+11', '-8.1818181818+12',
            '-8.1818181818+13', '-8.1818181818+14', '-8.1818181818+15', '-8.1818181818+16', ]
        for num, expectedi in zip(nums, expected):
            output = print_scientific_16(num)
            self.assertEqual(len(output), 16, msg='output=%r len(output)=%i' % (output, len(output)))
            self.assertEqual(output, expectedi, msg='num=%s output=%r expected=%r' % (num, output, expectedi))
            #print('%16s %r' % (x, output))

    def test_float_16(self):
        small_exponent = -17
        large_exponent = 17

        nums = (
            [0., 0.000034, -0.000034,
             0.000000000000000000000001, -0.000000000000000000000001] +
            [9./11 * 10**x for x in range(small_exponent, large_exponent+1)] +
            [-9./11 * 10**x for x in range(small_exponent, large_exponent+1)])

        expected = [
            '              0.', '         .000034', '        -.000034',
            '           1.-24', '          -1.-24',

            '8.18181818182-18', '8.18181818182-17', '8.18181818182-16', '8.18181818182-15',
            '8.18181818182-14', '8.18181818182-13', '8.18181818182-12', '8.18181818182-11',
            '8.18181818182-10', '8.181818181818-9', '8.181818181818-8', '8.181818181818-7',
            '8.181818181818-6', '8.181818181818-5', '8.181818181818-4', '.008181818181818',
            '.081818181818182', '.818181818181818', '8.18181818181818', '81.8181818181818',
            '818.181818181818', '8181.81818181818', '81818.1818181818', '818181.818181818',
            '8181818.18181818', '81818181.8181818', '818181818.181818', '8181818181.81818',
            '81818181818.1818', '818181818181.818', '8181818181818.18', '81818181818181.8',
            '818181818181818.', '8.18181818182+15', '8.18181818182+16',

            '-8.1818181818-18', '-8.1818181818-17', '-8.1818181818-16', '-8.1818181818-15',
            '-8.1818181818-14', '-8.1818181818-13', '-8.1818181818-12', '-8.1818181818-11',
            '-8.1818181818-10', '-8.18181818182-9', '-8.18181818182-8', '-8.18181818182-7',
            '-8.18181818182-6', '-8.18181818182-5', '-8.18181818182-4', '-8.18181818182-3',
            '-.08181818181818', '-.81818181818182', '-8.1818181818182', '-81.818181818182',
            '-818.18181818182', '-8181.8181818182', '-81818.181818182', '-818181.81818182',
            '-8181818.1818182', '-81818181.818182', '-818181818.18182', '-8181818181.8182',
            '-81818181818.182', '-818181818181.82', '-8181818181818.2', '-81818181818182.',
            '-8.1818181818+14', '-8.1818181818+15', '-8.1818181818+16',
        ]
        for num, expectedi in zip(nums, expected):
            output = print_float_16(num)
            self.assertEqual(len(output), 16, msg='output=%r len(output)=%i' % (output, len(output)))
            self.assertEqual(output, expectedi, msg='num=%s output=%r expected=%r' % (num, output, expectedi))

        nums = [0.99999999999999 * 10**x for x in range(small_exponent, large_exponent+1)]
        unused_positive_output = [print_float_16(x) for x in nums]
        unused_negative_output = [print_float_16(-x) for x in nums]

    def test_float_16_edge_cases(self):
        """Test print_float_16 at notation boundaries and precision limits."""
        expected_value = [
            # Zero / NaN
            ('              0.', 0.0),
            ('              0.', -0.0),
            ('                ', float('nan')),
            # Small values near scientific threshold (5e-16)
            ('           1.-16', 1e-16),
            ('          -1.-16', -1e-16),
            ('           5.-16', 5e-16),
            ('          -5.-16', -5e-16),
            # Fixed-point small values
            ('          .00001', 1e-5),
            ('         -.00001', -1e-5),
            ('            .001', 0.001),
            ('           -.001', -0.001),
            ('         .000034', 0.000034),
            ('        -.000034', -0.000034),
            # Fixed-point values at various magnitudes
            ('              1.', 1.0),
            ('             -1.', -1.0),
            ('         100000.', 1e5),
            ('        -100000.', -1e5),
            # Max fixed-point before scientific
            ('999999999999999.', 999999999999999.),
            # Scientific notation (large)
            ('           1.+15', 1e15),
            ('          -1.+15', -1e15),
            ('           1.+16', 1e16),
            ('          -1.+16', -1e16),
            # Full precision
            ('1.23456789012345', 1.23456789012345),
            ('-1.2345678901235', -1.23456789012345),
        ]
        for expected, value in expected_value:
            actual = print_float_16(value)
            self.assertEqual(len(actual), 16,
                             f'len={len(actual)} for value={value}: |{actual}|')
            self.assertEqual(actual, expected,
                             f'value={value}: actual=|{actual}| expected=|{expected}|')


    def test_print_int_default_8(self):
        """Test print_int_default_8: returns blank when value==default, else 8-char int."""
        # value == default -> blank
        self.assertEqual(print_int_default_8(0, 0), '        ')
        self.assertEqual(print_int_default_8(5, 5), '        ')
        self.assertEqual(print_int_default_8(-1, -1), '        ')
        # value != default -> formatted int
        self.assertEqual(print_int_default_8(1, 0), '       1')
        self.assertEqual(print_int_default_8(0, 1), '       0')
        self.assertEqual(print_int_default_8(-42, 0), '     -42')
        self.assertEqual(print_int_default_8(12345678, 0), '12345678')

    def test_print_float_default_8(self):
        """Test print_float_default_8: returns blank when value==default, else print_float_8."""
        # value == default -> blank
        self.assertEqual(print_float_default_8(0.0, 0.0), '        ')
        self.assertEqual(print_float_default_8(1.5, 1.5), '        ')
        self.assertEqual(print_float_default_8(-3.0, -3.0), '        ')
        # value != default -> formatted float
        self.assertEqual(print_float_default_8(1.5, 0.0), '     1.5')
        self.assertEqual(print_float_default_8(0.0, 1.0), '      0.')
        self.assertEqual(print_float_default_8(-0.5, 0.0), '     -.5')

    def test_set_string8_blank_if_default(self):
        """Test set_string8_blank_if_default: returns 8-char blank or formatted value."""
        # value == default -> 8 spaces
        self.assertEqual(set_string8_blank_if_default(0.0, 0.0), '        ')
        self.assertEqual(set_string8_blank_if_default(1, 1), '        ')
        # value != default -> 8-char formatted
        self.assertEqual(set_string8_blank_if_default(5.0, 0.0), '     5.0')
        self.assertEqual(set_string8_blank_if_default(42, 0), '      42')
        self.assertEqual(set_string8_blank_if_default('ABC', 'DEF'), '     ABC')

    def test_set_string16_blank_if_default(self):
        """Test set_string16_blank_if_default: returns 16-char blank or formatted value."""
        # value == default -> 16 spaces
        self.assertEqual(set_string16_blank_if_default(0.0, 0.0), '                ')
        self.assertEqual(set_string16_blank_if_default(1, 1), '                ')
        # value != default -> 16-char formatted
        self.assertEqual(set_string16_blank_if_default(5.0, 0.0), '             5.0')
        self.assertEqual(set_string16_blank_if_default(42, 0), '              42')
        self.assertEqual(set_string16_blank_if_default('ABC', 'DEF'), '             ABC')

    def test_print_field_8_numpy_types(self):
        """Test print_field_8 with numpy int32, int64, float32, float64."""
        # numpy integers
        self.assertEqual(print_field_8(np.int32(7)), '       7')
        self.assertEqual(print_field_8(np.int64(-100)), '    -100')
        self.assertEqual(print_field_8(np.int32(12345678)), '12345678')
        # numpy floats
        self.assertEqual(len(print_field_8(np.float32(1.5))), 8)
        self.assertEqual(print_field_8(np.float64(0.0)), '      0.')
        self.assertEqual(print_field_8(np.float64(125000.)), ' 125000.')
        # numpy float32 precision is lower, just verify length
        field = print_field_8(np.float32(3.14))
        self.assertEqual(len(field), 8)

    def test_print_field_16_numpy_types(self):
        """Test print_field_16 with numpy int32, int64, float32."""
        # numpy integers
        self.assertEqual(print_field_16(np.int32(7)), '               7')
        self.assertEqual(print_field_16(np.int64(-100)), '            -100')
        self.assertEqual(print_field_16(np.int64(1234567890123456)), '1234567890123456')
        # numpy floats
        self.assertEqual(print_field_16(np.float32(0.0)), '              0.')
        field = print_field_16(np.float32(3.14))
        self.assertEqual(len(field), 16)

    def test_print_card_8_empty_string_field(self):
        """Test that '' fields are treated as blanks in print_card_8."""
        # '' should produce same result as None (blank field)
        card_none = print_card_8(['GRID', 1, None, 0.0])
        card_empty = print_card_8(['GRID', 1, '', 0.0])
        self.assertEqual(card_none, card_empty)

    def test_print_card_8_mixed_types(self):
        """Test print_card_8 with floats, strings, ints, and None together."""
        # GRID card: name, nid, cp, x1, x2, x3, cd, ps, seid
        card = ['GRID', 1, None, 1.5, -2.0, 3.0, None, None, None]
        result = print_card_8(card)
        self.assertIn('GRID', result)
        self.assertIn('1.5', result)
        self.assertIn('-2.', result)
        self.assertIn('3.', result)
        # Verify field widths by checking total line length (80 chars with newline)
        lines = result.strip('\n').split('\n')
        for line in lines:
            self.assertLessEqual(len(line), 80)

        # MAT1 card: mixed ints and floats
        card = ['MAT1', 1, 3.0e7, None, 0.3, 0.1, None, None, None, 100.0]
        result = print_card_8(card)
        lines = result.strip('\n').split('\n')
        for line in lines:
            self.assertLessEqual(len(line), 80)

    def test_print_card_16_dedicated(self):
        """Dedicated tests for print_card_16 with is_double=False."""
        # Simple card: 1 continuation line
        card = ['GRID', 1, None, 1.0, 2.0, 3.0]
        result = print_card_16(card)
        self.assertTrue(result.startswith('GRID*'))
        lines = result.strip('\n').split('\n')
        for line in lines:
            self.assertLessEqual(len(line), 80)
        # Must contain continuation markers
        self.assertIn('*', result)

        # Card with multiple continuation lines
        card = ['CTETRA', 6437, 1, 14533, 5598, 1577, 9976, 42364, 5599, 42365, 42363, 30022, 12904]
        result = print_card_16(card)
        self.assertTrue(result.startswith('CTETRA*'))
        lines = result.strip('\n').split('\n')
        for line in lines:
            self.assertLessEqual(len(line), 80)

    def test_print_card_16_wipe_fields_false(self):
        """Test print_card_16 with wipe_fields=False preserves trailing blanks."""
        # With wipe_fields=True (default), trailing Nones are removed
        card = ['GRID', 1, None, 0.0, 0.0, 0.0, None, None, None]
        result_wipe = print_card_16(card, wipe_fields=True)
        result_no_wipe = print_card_16(card, wipe_fields=False)
        # wipe_fields=False should produce output at least as long
        # (trailing Nones remain as blank fields)
        self.assertGreaterEqual(len(result_no_wipe), len(result_wipe))

    def test_float_16_many(self):
        """Brute-force sweep for print_float_16 verifying 16-char width (analogous to test_float_8_many)."""
        for istart in np.arange(-13, 13):
            nums = np.logspace(istart, istart + 1, num=1000, endpoint=True, base=10.0)
            for num in nums:
                output = print_float_16(num)
                self.assertEqual(len(output), 16,
                                 msg='output=%r len(output)=%i num=%s' % (output, len(output), num))
                output_neg = print_float_16(-num)
                self.assertEqual(len(output_neg), 16,
                                 msg='output=%r len(output)=%i num=%s' % (output_neg, len(output_neg), -num))
                output_sci = print_scientific_16(num)
                self.assertEqual(len(output_sci), 16,
                                 msg='output=%r len(output)=%i num=%s' % (output_sci, len(output_sci), num))

    def test_print_scientific_double(self):
        """Test print_scientific_double produces 16-char D-notation fields."""
        expected_value = [
            # zero
            ('0.0000000000D+00', 0.0),
            # negative zero normalizes to positive
            ('0.0000000000D+00', -0.0),
            # positive values
            ('1.0000000000D+00', 1.0),
            ('1.2500000000D+02', 125.0),
            ('5.0000000000D-08', 5e-8),
            ('1.0000000000D+20', 1e20),
            ('9.9999999000D+00', 9.9999999),
            # negative values (1 fewer mantissa digit due to sign)
            ('-1.000000000D+00', -1.0),
            ('-1.250000000D+02', -125.0),
            ('-5.000000000D-08', -5e-8),
            ('-1.000000000D+20', -1e20),
            ('-9.999999900D+00', -9.9999999),
        ]
        for expected, value in expected_value:
            actual = print_scientific_double(value)
            self.assertEqual(len(actual), 16,
                             f'len={len(actual)} for value={value}: |{actual}|')
            self.assertEqual(actual, expected,
                             f'value={value}: actual=|{actual}| expected=|{expected}|')

    def test_print_scientific_double_many(self):
        """Brute-force sweep for print_scientific_double verifying 16-char width and D-notation."""
        for istart in np.arange(-13, 13):
            nums = np.logspace(istart, istart + 1, num=200, endpoint=True, base=10.0)
            for num in nums:
                output = print_scientific_double(num)
                self.assertEqual(len(output), 16,
                                 msg='output=%r len=%i num=%s' % (output, len(output), num))
                self.assertIn('D', output,
                              msg='output=%r num=%s missing D-notation' % (output, num))
                output_neg = print_scientific_double(-num)
                self.assertEqual(len(output_neg), 16,
                                 msg='output=%r len=%i num=%s' % (output_neg, len(output_neg), -num))
                self.assertIn('D', output_neg)

    def test_print_field_double(self):
        """Test print_field_double for int, float, None, and string inputs."""
        # integers -> right-justified 16-char
        self.assertEqual(print_field_double(1), '               1')
        self.assertEqual(print_field_double(12345678), '        12345678')
        # None -> 16 spaces
        self.assertEqual(print_field_double(None), '                ')
        # string -> right-justified 16-char
        self.assertEqual(print_field_double('THRU'), '            THRU')
        # floats -> scientific D-notation
        self.assertEqual(print_field_double(0.0), '0.0000000000D+00')
        self.assertEqual(print_field_double(1.0), '1.0000000000D+00')
        self.assertEqual(print_field_double(-1.0), '-1.000000000D+00')
        # all results must be 16 chars
        for val in [0.0, 1.0, -1.0, 1e-10, -1e-10, 1e20, -1e20, 3.14159, -0.001]:
            field = print_field_double(val)
            self.assertEqual(len(field), 16, f'value={val} field=|{field}|')

    def test_print_card_double_multi_continuation(self):
        """Test print_card_double with many fields producing multiple continuation lines."""
        # CHEXA: 2 continuation lines of data
        card = ['CHEXA', 1, 2, 10, 20, 30, 40, 50, 60, 70, 80]
        result = print_card_double(card)
        self.assertTrue(result.startswith('CHEXA*'))
        lines = result.strip('\n').split('\n')
        for line in lines:
            self.assertLessEqual(len(line), 80)
        # should have continuation markers
        cont_lines = [l for l in lines if l.startswith('*')]
        self.assertGreaterEqual(len(cont_lines), 1)

    def test_print_card_double_wipe_fields_false(self):
        """Test print_card_double with wipe_fields=False."""
        card = ['GRID', 1, None, 0.0, 0.0, 0.0, None, None, None]
        result_wipe = print_card_double(card, wipe_fields=True)
        result_no_wipe = print_card_double(card, wipe_fields=False)
        self.assertGreaterEqual(len(result_no_wipe), len(result_wipe))


def compare(value_in):
    field = print_field_8(value_in)
    val = interpret_value(field)
    if val != 0:
        p = (val - value_in) / val
        if p > 0.01:  # pragma: no cover
            raise ValueError('val=%s value_in=%s' % (val, value_in))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
