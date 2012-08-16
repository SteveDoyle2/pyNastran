from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import random
import unittest

#from numpy import allclose

from pyNastran.bdf.fieldWriter import print_field
from pyNastran.bdf.bdfInterface.bdf_cardMethods import interpretValue


class TestFieldWriter(unittest.TestCase):

    def test_1(self):
        field = print_field(1e20)
        assert '   1.+20' == field, '|%s|' % (field)
        field = print_field(-.723476)
        assert '-.723476' == field, '|%s|' % (field)
        field = print_field(
            125000.)
        assert ' 125000.' == field, '|%s|' % (field)
        field = print_field(
            12500000.)
        assert '  1.25+7' == field, '|%s|' % (field)
        field = print_field(
            47.77267)
        assert '47.77267' == field, '|%s|' % (field)
        field = print_field(
            .001)
        assert '    .001' == field, '|%s|' % (field)
        field = print_field(
            .0000001)
        assert '.0000001' == field, '|%s|' % (field)
        field = print_field(
            -5.007e-3)
        assert '-5.007-3' == field, '|%s|' % (field)
        field = print_field(
            -0.0748662)
        assert '-.074866' == field, '|%s|' % (field)

    def a_test_2(self):
        for i in xrange(1000):
            a = random.uniform(-20, 20)
            a2 = 10 ** a
            self.compare(a2)
            self.compare(-a2)

    def test_ints(self):
        self.assertEquals(print_field(1), '       1', 'a')
        self.assertEquals(print_field(12345678), '12345678', 'b')
        self.assertEquals(print_field('12345678'), '12345678', 'c')
        #self.assertEquals(print_field('1       '),'       1','|%s|' %(printField('1       ')))

    def test_floats_greater_than_1(self):
        self.assertEquals(print_field(1.2), '     1.2', 'a')
        self.assertEquals(print_field(1.23456789), '1.234568', 'b')
        self.assertEquals(print_field(12.234568), '12.23457', 'c')
        self.assertEquals(print_field(123.23457), '123.2346', 'd')
        self.assertEquals(print_field(1234.23468), '1234.235', 'e')
        self.assertEquals(print_field(12345.238), '12345.24', 'f')
        self.assertEquals(print_field(123456.28), '123456.3', 'g')
        self.assertEquals(print_field(1234567.25), '1234567.',
                          print_field(1234567.25))  # 1.2346+6
        self.assertEquals(print_field(
            12345678.), '1.2346+7', '|%s|' % print_field(12345678.))
        self.assertEquals(print_field(
            123456789.), '1.2346+8', '|%s|' % print_field(12345678.))

    def test_floats_small(self):
        self.assertEquals(print_field(
            0.1), '      .1', 'A|%s|' % print_field(0.1))
        self.assertEquals(print_field(
            0.0001), '   .0001', 'B|%s|' % print_field(0.0001))
        self.assertEquals(print_field(
            0.00001), '  .00001', 'C|%s|' % print_field(0.00001))
        self.assertEquals(print_field(
            0.000001), ' .000001', 'D|%s|' % print_field(0.000001))
        self.assertEquals(print_field(0.0000001), '.0000001',
                          'E|%s|' % print_field(0.0000001))
        self.assertEquals(print_field(0.00000012), '   1.2-7',
                          'F|%s|' % print_field(0.00000012))
        self.assertEquals(print_field(0.000748519), '7.4852-4',
                          'G|%s|' % print_field(0.000748519))
        self.assertEquals(print_field(0.12345678), '.1234568',
                          'H|%s|' % print_field(0.12345678))
        self.assertEquals(print_field(0.00012349), '1.2349-4',
                          'I|%s|' % print_field(0.00012349))
        self.assertEquals(print_field(0.000012349), '1.2349-5',
                          'J|%s|' % print_field(0.000012349))
        self.assertEquals(print_field(
            1e-20), '   1.-20', 'K|%s|' % print_field(1e-20))

    def compare(self, valueIn):
        #print "a = |%s|" %(valueIn)
        field = print_field(valueIn)
        print("a = |%s|" % (field))

        val = interpretValue(field)
        if val != 0:
            p = (val - valueIn) / val
            if p > 0.01:
                raise ValueError('val=%s valueIn=%s' % (val, valueIn))
        #assert allclose(val,valueIn),

if __name__ == "__main__":
    unittest.main()
