from __future__ import absolute_import

import unittest
from pyNastran.f06.f06_formatting import writeFloats8p4F

class TestFormatting(unittest.TestCase):

    def test_write_floats_8p4F(self):
        val = 0.0
        expected = '  0.0   '
        self.check_float_8p4f(val, expected)

        val = 1e-50
        expected = '  0.0   '
        self.check_float_8p4f(val, expected)

        val = 1e50
        self.assertRaises(RuntimeError, writeFloats8p4F, [val])

        val = 89.83581
        expected = ' 89.8358'
        self.check_float_8p4f(val, expected)

        val = 89.83586
        expected = ' 89.8359'
        self.check_float_8p4f(val, expected)

        val = -89.83581
        expected = '-89.8358'
        self.check_float_8p4f(val, expected)

        val = -89.83586
        expected = '-89.8359'
        self.check_float_8p4f(val, expected)

        val = -101.23451
        expected = '-101.2345'
        self.assertRaises(RuntimeError, writeFloats8p4F, [val])

        val = 101.23451
        expected = '101.2345'
        self.check_float_8p4f(val, expected)

        val = 101.23458
        expected = '101.2346'
        self.check_float_8p4f(val, expected)

    def check_float_8p4f(self, val, expected):
        actual, isAllZero = writeFloats8p4F([val])
        actuali = actual[0]
        self.assertEqual(actuali, expected, msg='\nactual  =%r len(actual)=%i\nexpected=%r len(expected)=%i' % (actuali, len(actuali), expected, len(expected)))

if __name__ == '__main__':  # pragma: no cover
    unittest.main()