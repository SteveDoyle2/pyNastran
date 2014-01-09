from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import random
import unittest

from pyNastran.bdf.fieldWriter import (print_field, print_float_8,
                                       set_default_if_blank,
                                       set_blank_if_default, is_same, print_card_8)
from pyNastran.bdf.fieldWriter16 import print_field_16
from pyNastran.bdf.field_writer_double import print_card_double


from pyNastran.bdf.bdfInterface.assign_type import interpret_value


class TestFieldWriter(unittest.TestCase):

    def test_field_vals_8(self):
        self.assertEqual(print_field(1e20),     '   1.+20',
                         print_field(1e20))
        self.assertEqual(print_field(-.723476), '-.723476',
                         print_field(-.723476))
        self.assertEqual(print_field(125000.),  ' 125000.',
                         print_field(125000.))
        self.assertEqual(print_field(12500000.),  '  1.25+7',
                         print_field(12500000.))
        self.assertEqual(print_field(47.77267),  '47.77267',
                         print_field(47.77267))
        self.assertEqual(print_field(.001),  '    .001',
                         print_field(.001))
        self.assertEqual(print_field(.0000001),  '.0000001',
                         print_field(.0000001))
        self.assertEqual(print_field(-5.007e-3),  '-5.007-3',
                         print_field(-5.007e-3))

        self.assertEqual(print_field(-0.0748662),  '-.074866',
                         print_field(-0.0748662))


    #@unittest.skip('skipping test_2')
    def test_2(self):
        for i in xrange(100):
            a = random.uniform(-20, 20)
            a2 = 10 ** a
            compare(a2)
            compare(-a2)

    def test_strings_8(self):
        self.assertEqual(print_field(None), '        ',
                         print_field(None))
        self.assertEqual(print_field('asdf'), '    asdf',
                         print_field('asdf'))
        self.assertEqual(print_field('  asdf  '), '  asdf  ',
                         print_field('  asdf  '))
        self.assertRaises(RuntimeError, print_field, '  asdf   ')

    def test_strings_16(self):
        self.assertEqual(print_field_16(None),         '                ',
                         print_field_16(None))
        self.assertEqual(print_field_16('asdf'),       '            asdf',
                         print_field_16('asdf'))
        self.assertEqual(print_field_16('  asdf  '),   '          asdf  ',
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
        self.assertEqual(print_field(1), '       1', 'a')
        self.assertEqual(print_field(12345678), '12345678', 'b')
        self.assertRaises(RuntimeError, print_field, 123456789)
        self.assertEqual(print_field('12345678'), '12345678', 'c')
        self.assertEqual(print_field('1       '),'1       ',
                         '|%s|' %(print_field('1       ')))

    def test_ints_16(self):
        self.assertEqual(print_field_16(1), '               1', 'a')
        self.assertEqual(print_field_16(12345678), '        12345678', 'b')
        self.assertEqual(print_field_16(1234567890123456), '1234567890123456', 'c')
        self.assertRaises(RuntimeError, print_field_16, 12345678901234567)

        #msg = print_field_16('12345678        ')
        #msg = '|%s| len(msg)=%s' %(msg, len(msg))
        #self.assertEqual(print_field_16('12345678'), '12345678        ',msg)
        self.assertEqual(print_field_16('1               '),'1               ',
                          '|%s|' %(print_field('1       ')))

    def test_floats_positive_8(self):
        tol = 1.0
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

        self.assertEqual(print_field(1.2), '     1.2', 'a')
        self.assertEqual(print_field(1.23456789), '1.234568', 'b')
        self.assertEqual(print_field(12.234568), '12.23457', 'c')
        self.assertEqual(print_field(123.23457), '123.2346', 'd')
        self.assertEqual(print_field(1234.23468), '1234.235', 'e')
        self.assertEqual(print_field(12345.238), '12345.24', 'f')
        self.assertEqual(print_field(123456.28), '123456.3', 'g')
        self.assertEqual(print_field(1234567.25), '1234567.',
                         print_field(1234567.25))  # 1.2346+6
        self.assertEqual(print_field(12345678.), '1.2346+7',
                         '%r' % print_field(12345678.))
        self.assertEqual(print_field(
            123456789.), '1.2346+8', '%r' % print_field(12345678.))

        self.assertEqual(print_field(0.1), '      .1',
                         'A|%s|' % print_field(0.1))
        self.assertEqual(print_field(0.0001), '   .0001',
                         'B|%s|' % print_field(0.0001))
        self.assertEqual(print_field(0.00001), '  .00001',
                         'C|%s|' % print_field(0.00001))
        self.assertEqual(print_field(0.000001), ' .000001',
                         'D|%s|' % print_field(0.000001))
        self.assertEqual(print_field(0.0000001), '.0000001',
                         'E|%s|' % print_field(0.0000001))
        self.assertEqual(print_field(0.00000012), '   1.2-7',
                         'F|%s|' % print_field(0.00000012))
        self.assertEqual(print_field(0.000748519), '7.4852-4',
                         'G|%s|' % print_field(0.000748519))
        self.assertEqual(print_field(0.12345678), '.1234568',
                         'H|%s|' % print_field(0.12345678))
        self.assertEqual(print_field(0.00012349), '1.2349-4',
                         'I|%s|' % print_field(0.00012349))
        self.assertEqual(print_field(0.000012349), '1.2349-5',
                         'J|%s|' % print_field(0.000012349))
        self.assertEqual(print_field(5e-08),       '    5.-8',
                         'K|%s|' % print_field(5e-08))

        self.assertEqual(print_field(1e-20), '   1.-20',
                         'L|%s|' % print_field(1e-20))

        self.assertEqual(print_field(0.0),  '      0.',
                         print_field(0.0))
        self.assertEqual(print_field(1.0),  '      1.',
                         print_field(1.0))

    def test_floats_negative_8(self):
        self.assertEqual(print_field(-1.2), '    -1.2',
                         print_field(-1.2))
        self.assertEqual(print_field(-1.23456789), '-1.23457',
                         print_field(-1.23456789))
        self.assertEqual(print_field(-12.234568),  '-12.2346',
                         print_field(-12.234568))
        self.assertEqual(print_field(-123.23457),  '-123.235',
                         print_field(-123.23457))
        self.assertEqual(print_field(-1234.23468), '-1234.23',
                         print_field(-1234.23468))
        self.assertEqual(print_field(-12345.238),  '-12345.2',
                         print_field(-12345.238))
        self.assertEqual(print_field(-123456.28),  '-123456.',
                         print_field(-123456.28))
        self.assertEqual(print_field(-1234567.25), '-1.235+6', # is this right?
                         print_field(-1234567.25))
        self.assertEqual(print_field(-12345678.),  '-1.235+7', # is this right?
                         '|%s|' % print_field(-12345678.))
        self.assertEqual(print_field(-123456789.), '-1.235+8', # is this right?
                         '|%s|' % print_field(-123456789.))

        self.assertEqual(print_field(-0.1), '     -.1',
                         'A|%s|' % print_field(-0.1))
        self.assertEqual(print_field(-0.0001), '  -.0001',
                         'B|%s|' % print_field(-0.0001))
        self.assertEqual(print_field(-0.00001), ' -.00001',
                         'C|%s|' % print_field(-0.00001))
        self.assertEqual(print_field(-0.000001), '   -1.-6',
                         'D|%s|' % print_field(-0.000001))
        self.assertEqual(print_field(-0.0000001), '   -1.-7',
                         'E|%s|' % print_field(-0.0000001))
        self.assertEqual(print_field(-0.00000012), '  -1.2-7',
                         'F|%s|' % print_field(-0.00000012))
        self.assertEqual(print_field(-0.000748519), '-7.485-4',
                         'G|%s|' % print_field(-0.000748519))
        self.assertEqual(print_field(-0.12345678), '-.123457',
                         'H|%s|' % print_field(-0.12345678))
        self.assertEqual(print_field(-0.00012349), '-1.235-4',
                         'I|%s|' % print_field(-0.00012349))
        self.assertEqual(print_field(-0.000012349), '-1.235-5',
                         'J|%s|' % print_field(-0.000012349))
        self.assertEqual(print_field(-1e-5),  ' -.00001',
                         'K|%s|' % print_field(-1e-5))
        self.assertEqual(print_field(-1e-6),  '   -1.-6',
                         'L|%s|' % print_field(-1e-6))
        self.assertEqual(print_field(-1e-7),  '   -1.-7',
                         'M|%s|' % print_field(-1e-7))
        self.assertEqual(print_field(-1e-20), '  -1.-20',
                         'N|%s|' % print_field(-1e-20))

    def test_print_card_8(self):
        self.assertEqual(print_card_8(['GRID',1]),'GRID           1\n')
        self.assertEqual(print_card_8(['GRID', 1, None, None, None, None, None,
                                     None, None, 1]),
                         'GRID           1\n               1\n',
                         print_card_8(['GRID', 1, None, None, None, None, None,
                                     None, None, 1]))

        self.assertEqual(print_card_8(['GRID', 1, None, None, None, None, None,
                                     None, None,
                                     1, None]),
                         'GRID           1\n               1\n',
                         print_card_8(['GRID', 1, None, None, None, None, None,
                                     None, None,
                                     1, None]))

        self.assertEqual(print_card_8(['GRID', 1, None, None, None, None, None,
                                     None, None,
                                     None,None,None,None,None,None,None,None,
                                     1, None]),
                         'GRID           1\n+\n               1\n',
                         print_card_8(['GRID', 1, None, None, None, None, None,
                                     None, None,
                                     None,None,None,None,None,None,None,None,
                                     1, None]))

    def test_same(self):
        self.assertTrue(is_same(1.0, 1.000))
        self.assertFalse(is_same(1.0, 1e-15 + 1))
        self.assertTrue(is_same('MAT', 'MAT'))
        self.assertFalse(is_same('MAT', 'MAT1'))

    def test_card_double(self):
        card = print_card_double(['GRID', 1, None, 120.322,-4.82872,1.13362,])
        card_expected = 'GRID*                  1                1.2032200000D+02-4.828720000D+00\n'
        card_expected += '*       1.1336200000D+00\n'
        self.assertEqual(card, card_expected)
        print('1')

def compare(valueIn):
    #print "a = |%s|" % valueIn
    field = print_field(valueIn)
    #print("a = |%s|" % field)

    val = interpret_value(field)
    if val != 0:
        p = (val - valueIn) / val
        if p > 0.01:
            raise ValueError('val=%s valueIn=%s' % (val, valueIn))

if __name__ == "__main__":
    unittest.main()
