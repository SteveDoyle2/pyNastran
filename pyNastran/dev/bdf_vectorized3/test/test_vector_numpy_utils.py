import unittest
import numpy as np
from functools import partial
from pyNastran.dev.bdf_vectorized3.cards.base_card import searchsorted_filter, searchsorted_filter_
from pyNastran.dev.bdf_vectorized3.bdf_interface.fast_float_print import print_float_8, compare

class TestNumpyExtensions(unittest.TestCase):
    def test_print_float_8(self):
        compare_print_float_8 = partial(compare, stop_on_error=True)
        value = -4.32693081e-05# ; field_old='-4.327-5' field_new='-.000043'
        s = compare_print_float_8(value)
        value = -3.755427e-06# ; field_old='-3.755-6' field_new='-.000004'
        s = compare_print_float_8(value)

        #value=-4.32693081e-05; field_old='-4.327-5' field_new='-.000043'
        #value=-3.755427e-06; field_old='-3.755-6' field_new='-.000004'
        value=-0.8732806594999999 #; field_old='-.873281' field_new='-8.733-1'
        s = compare_print_float_8(value)
        value=-0.6862676525499999 #; field_old='-.686268' field_new='-6.863-1'
        s = compare_print_float_8(value)
        value=-0.647050147 #; field_old=' -.64705' field_new='-6.471-1'
        s = compare_print_float_8(value)
        value=-0.9526646215 #; field_old='-.952665' field_new='-9.527-1'
        s = compare_print_float_8(value)
        value=-0.6291313619 #; field_old='-.629131' field_new='-6.291-1'
        s = compare_print_float_8(value)

        x = 2e-11
        s = compare_print_float_8(x)
        x = 12e-11
        s = compare_print_float_8(x)
        x = 12.1e-11
        s = compare_print_float_8(x)
        x = 12.12e-11
        s = compare_print_float_8(x)
        x = 12.123e-11
        s = compare_print_float_8(x)
        x = 12.1234e-11
        s = compare_print_float_8(x)
        x = 12e-12
        s = compare_print_float_8(x)
        x = 12.1e-12
        s = compare_print_float_8(x)
        x = 12.12e-12
        s = compare_print_float_8(x)
        x = 12.123e-12
        s = compare_print_float_8(x)
        x = 12.1234e-12
        s = compare_print_float_8(x)
        x = 1.0e-5
        s = compare_print_float_8(x)
        x = 1.0e-4
        s = compare_print_float_8(x)
        x = 1.0e-3
        s = compare_print_float_8(x)
        x = 1.0e-2
        s = compare_print_float_8(x)
        x = 1.1e-5
        s = compare_print_float_8(x)
        x = 1.1e-4
        s = compare_print_float_8(x)
        x = 1.1e-3
        s = compare_print_float_8(x)
        x = 1.1e-2
        s = compare_print_float_8(x)
        x = 1.12e-5
        s = compare_print_float_8(x)
        x = 1.12e-4
        s = compare_print_float_8(x)
        x = 1.12e-3
        s = compare_print_float_8(x)
        x = 1.12e-2
        s = compare_print_float_8(x)
        x = 1.123e-5
        s = compare_print_float_8(x)
        x = 1.123e-4
        s = compare_print_float_8(x)
        x = 1.123e-3
        s = compare_print_float_8(x)
        x = 1.123e-2
        s = compare_print_float_8(x)
        x = 1.e5
        s = compare_print_float_8(x)
        x = 1.e6
        s = compare_print_float_8(x)
        x = 1.e7
        s = compare_print_float_8(x)
        x = 1.e8
        s = compare_print_float_8(x)
        x = 1.2e8
        s = compare_print_float_8(x)
        x = 1.23e8
        s = compare_print_float_8(x)
        x = 1.234e8
        s = compare_print_float_8(x)
        x = 1.2345e8
        s = compare_print_float_8(x)
        x = 1.23456e8
        s = compare_print_float_8(x)
        #------------------------------------------------
        x = 12345678.
        s = compare_print_float_8(x)
        x = 1234567.8
        s = compare_print_float_8(x)
        x = 123456.78
        s = compare_print_float_8(x)
        x = 12345.678
        s = compare_print_float_8(x)
        x = 1234.5678
        s = compare_print_float_8(x)
        x = 123.45678
        s = compare_print_float_8(x)
        x = 12.345678
        s = compare_print_float_8(x)
        x = 1.2345678
        s = compare_print_float_8(x)
        x = .12345678
        s = compare_print_float_8(x)
        x = .012345678
        s = compare_print_float_8(x)
        x = .0012345678
        s = compare_print_float_8(x)
        x = .00012345678
        s = compare_print_float_8(x)
        x = .000012345678
        s = compare_print_float_8(x)
        x = .0000012345678
        s = compare_print_float_8(x)


    def test_searchsorted_filter1(self):
        all_ids = [1, 0]
        lookup_ids = [1]
        with self.assertRaises(RuntimeError):
            searchsorted_filter_(all_ids, lookup_ids)

        all_ids          = [0, 1, 2, 3, 4]
        all_ids_expected = [   1,       4]

        lookup_ids       = [   1,       4, 5]
        ilookup_expected = [   0,       1]
        ilookup, iall = searchsorted_filter_(all_ids, lookup_ids)
        assert np.array_equal(ilookup, ilookup_expected)
        assert np.array_equal(iall, all_ids_expected)

        all_ids          = [1, 2, 3, 4, 5]
        all_ids_expected = [   1,       4]

        lookup_ids       = [   2,       5, 6]
        ilookup_expected = [   0,       1]
        ilookup, iall = searchsorted_filter_(all_ids, lookup_ids)
        assert np.array_equal(ilookup, ilookup_expected)
        assert np.array_equal(iall, all_ids_expected)

    def test_searchsorted_filter2(self):
        all_ids = np.array([2, 20, 30, 300705], dtype='int32')
        v_all = all_ids
        #lookup_ids = [1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704]
        lookup_ids_ = [
            [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[1, 1, 1, 1, 2, 2, 1, 1, 1, 1],
            #[300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704], [300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704],
            #[300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704], [300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704],
            #[300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704], [300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704],
            #[300704, 300704, 300704, 300704, 300705, 300705, 300704, 300704, 300704, 300704]
        ]
        lookup_ids = np.ravel(lookup_ids_)
        all_ids_expected_value = np.where(all_ids == 2)[0]
        ilookup_expected = np.where(lookup_ids == 2)[0]
        all_ids_expected = np.ones(len(ilookup_expected), dtype='int32') * all_ids_expected_value

        ilookup, iall = searchsorted_filter(all_ids, lookup_ids)
        values_all = v_all[iall]
        assert np.array_equal(ilookup, ilookup_expected)
        assert np.array_equal(iall, all_ids_expected)



