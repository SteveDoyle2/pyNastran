import unittest
from pyNastran.f06.f06_formatting import (
    write_floats_8p4f, write_floats_8p1e,
    write_floats_10e, write_floats_12e,
    write_imag_floats_13e)
from pyNastran.f06.f06_writer import (
    make_end, sorted_bulk_data_header, make_f06_header, make_stamp)

class TestF06Formatting(unittest.TestCase):

    def test_write_floats_8p4f(self):
        """testing write_floats_8p4f"""
        func = write_floats_8p4f
        val = 0.0
        expected = '  0.0   '
        self.check_floats(func, val, expected)

        val = 1e-50
        expected = '  0.0   '
        self.check_floats(func, val, expected)

        val = 1e50
        self.assertRaises(RuntimeError, write_floats_8p4f, [val])

        val = 89.83581
        expected = ' 89.8358'
        self.check_floats(func, val, expected)

        val = 89.83586
        expected = ' 89.8359'
        self.check_floats(func, val, expected)

        val = -89.83581
        expected = '-89.8358'
        self.check_floats(func, val, expected)

        val = -89.83586
        expected = '-89.8359'
        self.check_floats(func, val, expected)

        val = -101.23451
        expected = '-101.2345'
        self.assertRaises(RuntimeError, write_floats_8p4f, [val])

        val = 101.23451
        expected = '101.2345'
        self.check_floats(func, val, expected)

        val = 101.23458
        expected = '101.2346'
        self.check_floats(func, val, expected)

    def test_write_floats_8p1e(self):
        """testing write_floats_8p1e"""
        func = write_floats_8p1e

        val = 0.0
        expected = '       '
        self.check_floats(func, val, expected)

        val = 101.23458
        expected = ' 1.0E+02'
        self.check_floats(func, val, expected)

        val = -101.23458
        expected = '-1.0E+02'
        self.check_floats(func, val, expected)

    def test_write_floats_10e(self):
        """testing write_floats_10e"""
        func = write_floats_10e
        val = 0.0
        expected = ' 0.0'
        self.check_floats(func, val, expected)

        val = 1.0
        expected = ' 1.000E+00'
        self.check_floats(func, val, expected)

        val = -1.0
        expected = '-1.000E+00'
        self.check_floats(func, val, expected)

    def test_write_floats12e(self):
        """testing write_floats_12e"""
        func = write_floats_12e
        val = 0.0
        expected = ' 0.0'
        self.check_floats(func, val, expected)

        val = 1.0
        expected = ' 1.00000E+00'
        self.check_floats(func, val, expected)

        val = -1.0
        expected = '-1.00000E+00'
        self.check_floats(func, val, expected)

    def check_floats(self, func, val, expected):
        """helper method"""
        actual = func([val])
        actuali = actual[0]
        self.assertEqual(actuali, expected, msg='\nactual  =%r len(actual)=%i\nexpected=%r len(expected)=%i' % (
            actuali, len(actuali), expected, len(expected)))

    def test_write_imag_floats_13e(self):
        """testing write_imag_floats_13e"""
        func = write_imag_floats_13e

        val = 0.0 + 0.0j
        expected_real = ' 0.0'
        expected_imag = '   0.0'
        is_mag_phase = True
        self.check_imag_floats(func, val, expected_real, expected_imag, is_mag_phase)

        val = 1.0 + 0.0j
        expected_real = ' 1.000000E+00'
        expected_imag = '   0.0'
        is_mag_phase = True
        self.check_imag_floats(func, val, expected_real, expected_imag, is_mag_phase)

        val = 0.0 + 1.0j
        expected_real = ' 1.000000E+00'
        expected_imag = '90.0000      '
        is_mag_phase = True
        self.check_imag_floats(func, val, expected_real, expected_imag, is_mag_phase)

        val = 0.0 - 1.0j
        expected_real = ' 1.000000E+00'
        expected_imag = '270.0000     '
        is_mag_phase = True
        self.check_imag_floats(func, val, expected_real, expected_imag, is_mag_phase)
        #--------------------------------------------------------------------------
        val = 0.0 + 0.0j
        expected_real = ' 0.0'
        expected_imag = ' 0.0'
        is_mag_phase = False
        self.check_imag_floats(func, val, expected_real, expected_imag, is_mag_phase)

        val = 1.0 + 0.0j
        expected_real = ' 1.000000E+00'
        expected_imag = ' 0.0'
        is_mag_phase = False
        self.check_imag_floats(func, val, expected_real, expected_imag, is_mag_phase)

        val = 0.0 + 1.0j
        expected_real = ' 0.0'
        expected_imag = ' 1.000000E+00'
        is_mag_phase = False
        self.check_imag_floats(func, val, expected_real, expected_imag, is_mag_phase)

        val = 0.0 - 1.0j
        expected_real = ' 0.0'
        expected_imag = '-1.000000E+00'
        is_mag_phase = False
        self.check_imag_floats(func, val, expected_real, expected_imag, is_mag_phase)

    def check_imag_floats(self, func, val, expected_real, expected_imag, is_mag_phase):
        """helper method"""
        actual_real, actual_imag = func([val], is_mag_phase)
        self.assertEqual(actual_real, expected_real,
                         msg='\nreal:%s+%sj\nactual  =%r len(actual)=%i\nexpected=%r len(expected)=%i' % (
            val.real, val.imag, actual_real, len(actual_real), expected_real, len(expected_real)))
        self.assertEqual(actual_imag, expected_imag,
                         msg='\nimag %s+%sj:\nactual  =%r len(actual)=%i\nexpected=%r len(expected)=%i' % (
            val.real, val.imag, actual_imag, len(actual_imag), actual_imag, len(expected_imag)))

    def test_make_end(self):
        """miscellaneous F06 tester"""
        make_end(end_flag=True, options=None)
        options = ['SEMR']
        make_end(end_flag=True, options=options)
        sorted_bulk_data_header()
        make_f06_header()
        make_stamp(title=None)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
