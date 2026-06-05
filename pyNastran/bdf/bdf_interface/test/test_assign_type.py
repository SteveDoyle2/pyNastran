"""tests various assign_types.py functions"""
import unittest
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank,
    double, double_or_blank, double_from_str,
    integer_or_double, integer_double_or_blank,
    string, string_or_blank, double_or_string, double_string_or_blank,
    integer_or_string, integer_string_or_blank, integer_double_or_string,
    blank, parse_components, components_or_blank, integer_double_string_or_blank,
    _get_dtype, interpret_value, modal_components,
    parse_components_or_blank,
    modal_components_or_blank, string_choice_or_blank,
    exact_string_or_blank, filename_or_blank, check_string, loose_string,
)
from pyNastran.bdf.bdf_interface.assign_type_force import (
    force_integer, force_double, force_integer_or_blank,
    force_double_or_string, force_double_or_blank, lax_double_or_blank,
    parse_components as force_components
)

class TestAssignType(unittest.TestCase):
    """tests various assign_types.py functions"""

    def run_function_default(self, func, card, exact, default):
        """
        Helper function

        Parameters
        ----------
        func : function
           integer_or_blank
        card : list[varies]
            a series of values to add
        exacts : list[float]
            list of results
        default : list[float]
            list of default values
        """
        fieldname = 'f'
        assert len(card) == len(exact), 'len(card)=%s len(exact)=%s' % (len(card), len(exact))
        assert len(card) == len(default), 'len(card)=%s len(default)=%s' % (len(card), len(default))
        i = 0
        bdf_card = BDFCard(card)
        for i, exacti in enumerate(exact):
            defaulti = default[i]
            if exacti == SyntaxError:
                with self.assertRaises(exacti):
                    #msg = 'field=%r exact=%r default=%r' % (bdf_card.field(i), exacti, defaulti)
                    #print(msg)
                    func(bdf_card, i, fieldname, defaulti)
            else:
                value = func(bdf_card, i, fieldname, defaulti)
                self.assertEqual(value, exacti)
            i += 1

    def run_function(self, func, card, exact):
        """
        Helper function

        Parameters
        ----------
        func : function
           integer_or_blank
        card : list[varies]
            a series of values to add
        exacts : list[float]
            list of results
        """
        assert len(card) == len(exact), 'len(card)=%s len(exact)=%s' % (len(card), len(exact))
        i = 0
        fieldname = 'f'
        bdf_card = BDFCard(card)
        for i, exacti in enumerate(exact):
            if exacti == SyntaxError:
                with self.assertRaises(SyntaxError):
                    #msg = 'field=%r exact=%r' % (bdf_card.field(i), exacti)
                    #print msg
                    func(bdf_card, i, fieldname)
            else:

                value = func(bdf_card, i, fieldname)
                self.assertEqual(value, exacti)
            i += 1

    def test_blank(self):
        """
        value = integer(card, n, fieldname)
        """
        # integer
        with self.assertRaises(SyntaxError):
            blank(BDFCard([1]), 0, 'field')
        with self.assertRaises(SyntaxError):
            blank(BDFCard(['1']), 0, 'field')

        # float
        with self.assertRaises(SyntaxError):
            blank(BDFCard([1.]), 0, 'field')
        with self.assertRaises(SyntaxError):
            blank(BDFCard(['1.']), 0, 'field')

        # string
        with self.assertRaises(SyntaxError):
            blank(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            blank(BDFCard(['1b']), 0, 'field')

        # blank
        val = blank(BDFCard(['']), 0, 'field')
        self.assertEqual(val, None)
        val = blank(BDFCard([None]), 0, 'field')
        self.assertEqual(val, None)
        val = blank(BDFCard([None]), 0, 'field', 'default')
        self.assertEqual(val, 'default')

    def test_parse_components_or_blank(self):
        with self.assertRaises(SyntaxError):
            parse_components_or_blank(BDFCard([1.], has_none=False), 0, 'field')
        parse_components_or_blank(BDFCard([1], has_none=False), 0, 'field')

        # out of range
        # with self.assertRaises(SyntaxError):
        parse_components_or_blank(BDFCard([1.]), 1, 'field')

        # integer
        # with self.assertRaises(SyntaxError):
        parse_components_or_blank(BDFCard([1]), 0, 'field')
        # with self.assertRaises(SyntaxError):
        parse_components_or_blank(BDFCard(['1']), 0, 'field')
        with self.assertRaises(SyntaxError):
            # too large
            parse_components_or_blank(BDFCard(['7']), 0, 'field')

        # float
        with self.assertRaises(SyntaxError):
            self.assertEqual(1.e-9, parse_components_or_blank(BDFCard(['1-9']), 0, 'field'))
        with self.assertRaises(SyntaxError):
            self.assertEqual(1.e+9, parse_components_or_blank(BDFCard(['1+9']), 0, 'field'))

        # string
        with self.assertRaises(SyntaxError):
            parse_components_or_blank(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            parse_components_or_blank(BDFCard(['1b']), 0, 'field')

        # blank
        # with self.assertRaises(SyntaxError):
        parse_components_or_blank(BDFCard(['']), 0, 'field')
        # with self.assertRaises(SyntaxError):
        parse_components_or_blank(BDFCard([None]), 0, 'field')


        with self.assertRaises(SyntaxError):
            parse_components_or_blank(BDFCard(['.']), 0, 'field')
        parse_components_or_blank(BDFCard(['4']), 0, 'field')

        # card = [1.0, '2.0', '3.', 'C', None, '']
        # exact = [1.0, 2.0, 3.0, SyntaxError, SyntaxError, SyntaxError]
        # self.run_function(double, card, exact)

    def test_force_components(self):
        """tests the force_components function"""
        with self.assertRaises(SyntaxError):
            force_components(BDFCard(['1b']), 0, 'field')

        force_components(BDFCard([1], has_none=False), 0, 'field')

        # single ints
        val = force_components(BDFCard([0]), 0, 'field')
        self.assertEqual(val, '0')

        val = force_components(BDFCard([1]), 0, 'field')
        self.assertEqual(val, '1')

        # single strings
        val = force_components(BDFCard(['0']), 0, 'field')
        self.assertEqual(val, '0')

        val = force_components(BDFCard(['1']), 0, 'field')
        self.assertEqual(val, '1')

        # double ints
        val = force_components(BDFCard(['123']), 0, 'field')
        self.assertEqual(val, '123')

        val = force_components(BDFCard([123]), 0, 'field')
        self.assertEqual(val, '123')

        val = force_components(BDFCard([321]), 0, 'field')
        self.assertEqual(val, '123')

        # embedded whiteshape
        self.assertEqual(force_components(BDFCard(['12 3']), 0, 'field'), '123')

        # all numbers
        val = force_components(BDFCard(['123456']), 0, 'field')
        self.assertEqual(val, '123456')

        # invalid 0's defined with numbers
        with self.assertRaises(SyntaxError):
            force_components(BDFCard(['0123456']), 0, 'field')

        with self.assertRaises(SyntaxError):
            force_components(BDFCard(['01']), 0, 'field')

        # doubles
        with self.assertRaises(SyntaxError):
            force_components(BDFCard(['4524']), 0, 'field')

        # only 0 to 6
        with self.assertRaises(SyntaxError):
            force_components(BDFCard(['7']), 0, 'field')

        with self.assertRaises(SyntaxError):
            force_components(BDFCard([7]), 0, 'field')

        # dumb input
        with self.assertRaises(SyntaxError):
            force_components(BDFCard(['4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            force_components(BDFCard(['-4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            force_components(BDFCard(['asdf']), 0, 'field')

        with self.assertRaises(SyntaxError):
            force_components(BDFCard(['-1']), 0, 'field')

        # blank
        #force_components(BDFCard(['   ']), 0, 'field')
        #force_components(BDFCard([None]), 0, 'field')

    def test_force_double(self):
        """
        value = force_double(card, n, fieldname)
        """
        # out of range
        with self.assertRaises(SyntaxError):
            force_double(BDFCard([1.]), 1, 'field')

        # integer
        with self.assertRaises(SyntaxError):
            double(BDFCard([1]), 0, 'field')
        # with self.assertRaises(SyntaxError):
        self.assertEqual(force_double(BDFCard(['1']), 0, 'field'), 1.0)

        self.assertEqual(1.e-9, force_double(BDFCard(['1-9']), 0, 'field'))
        self.assertEqual(1.e+9, force_double(BDFCard(['1+9']), 0, 'field'))

        # float
        self.check_double(force_double)

        # string
        with self.assertRaises(SyntaxError):
            force_double(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            double(BDFCard(['1b']), 0, 'field')

        # blank
        with self.assertRaises(SyntaxError):
            force_double(BDFCard(['']), 0, 'field')
        with self.assertRaises(SyntaxError):
            force_double(BDFCard([None]), 0, 'field')

        with self.assertRaises(SyntaxError):
            force_double(BDFCard(['.']), 0, 'field')
        # with self.assertRaises(SyntaxError):
        self.assertEqual(force_double(BDFCard(['4']), 0, 'field'), 4.0)

        card = [1.0, '2.0', '3.', 'C', None, '']
        exact = [1.0, 2.0, 3.0, SyntaxError, SyntaxError, SyntaxError]
        self.run_function(force_double, card, exact)

    def test_force_integer(self):
        """
        value = force_integer(card, n, fieldname)
        """
        # integer
        self.check_integer(force_integer)

        # out of range
        with self.assertRaises(SyntaxError):
            integer(BDFCard([1.]), 1, 'field')

        # float
        self.assertEqual(force_integer(BDFCard([1.], has_none=False), 0, 'field'), 1)
        with self.assertRaises(SyntaxError):
            force_integer(BDFCard(['1-2']), 0, 'field')
        self.assertEqual(force_integer(BDFCard([1.]), 0, 'field'), 1)
        self.assertEqual(force_integer(BDFCard(['1.']), 0, 'field'), 1)
        self.assertEqual(force_integer(BDFCard([1.]), 0, 'field'), 1)

        # string
        with self.assertRaises(SyntaxError):
            force_integer(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            force_integer(BDFCard(['1b']), 0, 'field')

        # blank
        with self.assertRaises(SyntaxError):
            integer(BDFCard(['']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard([None]), 0, 'field')

        # with self.assertRaises(SyntaxError):
        self.assertEqual(force_integer(BDFCard(['7.0']), 0, 'field'), 7)

        with self.assertRaises(SyntaxError):
            force_integer(BDFCard(['1+2']), 0, 'field')

        card = [1, '2', '3.', 'C', None, '']
        exact = [1, 2, 3, SyntaxError, SyntaxError, SyntaxError]
        self.run_function(force_integer, card, exact)

    def test_force_integer_or_blank(self):
        """
        value = force_integer_or_blank(card, n, fieldname)
        """
        # integer
        self.check_integer(force_integer_or_blank)

        # float
        self.assertEqual(force_integer_or_blank(BDFCard([1.]), 0, 'field'), 1)
        self.assertEqual(force_integer_or_blank(BDFCard(['1.']), 0, 'field'), 1)
        #with self.assertRaises(SyntaxError):
            #integer(BDFCard(['cat']), 0, 'field')

        # string
        with self.assertRaises(SyntaxError):
            force_integer_or_blank(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            force_integer_or_blank(BDFCard(['1b']), 0, 'field')

        self.check_blank(force_integer_or_blank)

        with self.assertRaises(SyntaxError):
            force_integer_or_blank(BDFCard(['1+2']), 0, 'field')

        card = [1, '2', '3.', 'C', None, '']
        exact = [1, 2, 3, SyntaxError, None, None]
        default = [None, None, None, None, None, None]
        self.run_function_default(force_integer_or_blank, card, exact, default)

    def test_force_double_or_string(self):
        """tests the force_double_or_string function"""
        # out of range
        with self.assertRaises(SyntaxError):
            force_double_or_string(BDFCard([1.]), 1, 'field')
        self.assertEqual(1.e-9, force_double_or_string(BDFCard(['1-9']), 0, 'field'))
        self.assertEqual(1.e+9, force_double_or_string(BDFCard(['1+9']), 0, 'field'))

        self.check_double(force_double_or_string, is_strict=False)
        self.check_string(force_double_or_string, check_dash=False)

    def test_double(self):
        """
        value = double(card, n, fieldname)
        """
        # out of range
        with self.assertRaises(SyntaxError):
            double(BDFCard([1.]), 1, 'field')

        # integer
        with self.assertRaises(SyntaxError):
            double(BDFCard([1]), 0, 'field')
        with self.assertRaises(SyntaxError):
            double(BDFCard(['1']), 0, 'field')
        # signed integer strings
        with self.assertRaises(SyntaxError):
            double(BDFCard(['-400510']), 0, 'field')
        with self.assertRaises(SyntaxError):
            double(BDFCard(['+123']), 0, 'field')

        self.assertEqual(1.e-9, double(BDFCard(['1-9']), 0, 'field'))
        self.assertEqual(1.e+9, double(BDFCard(['1+9']), 0, 'field'))

        # D/E-notation without explicit sign on exponent (no decimal point)
        self.assertEqual(1.e3, double(BDFCard(['1D3']), 0, 'field'))
        self.assertEqual(1.e3, double(BDFCard(['1E3']), 0, 'field'))

        # float
        self.check_double(double)

        # string
        with self.assertRaises(SyntaxError):
            double(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            double(BDFCard(['1b']), 0, 'field')

        # blank
        with self.assertRaises(SyntaxError):
            double(BDFCard(['']), 0, 'field')
        with self.assertRaises(SyntaxError):
            double(BDFCard([None]), 0, 'field')


        with self.assertRaises(SyntaxError):
            double(BDFCard(['.']), 0, 'field')
        with self.assertRaises(SyntaxError):
            double(BDFCard(['4']), 0, 'field')

        card = [1.0, '2.0', '3.', 'C', None, '']
        exact = [1.0, 2.0, 3.0, SyntaxError, SyntaxError, SyntaxError]
        self.run_function(double, card, exact)

    def test_double_from_str(self):
        """tests the double_from_str function"""
        self.assertEqual(1.0, double_from_str('1.0'))
        self.assertEqual(-3.14, double_from_str('-3.14'))
        self.assertEqual(1.e-9, double_from_str('1-9'))
        self.assertEqual(1.e+9, double_from_str('1+9'))
        self.assertEqual(1.e3, double_from_str('1D3'))
        self.assertEqual(1.e3, double_from_str('1E3'))
        # integers must be rejected
        with self.assertRaises(SyntaxError):
            double_from_str('400510')
        with self.assertRaises(SyntaxError):
            double_from_str('-400510')
        with self.assertRaises(SyntaxError):
            double_from_str('+123')

    def test_integer(self):
        """
        value = integer(card, n, fieldname)
        """
        # integer
        self.check_integer(integer)

        # out of range
        with self.assertRaises(SyntaxError):
            integer(BDFCard([1.]), 1, 'field')

        # float
        with self.assertRaises(SyntaxError):
            integer(BDFCard([1.], has_none=False), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard(['1-2']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard([1.]), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard(['1.']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard([1.]), 0, 'field')

        # string
        with self.assertRaises(SyntaxError):
            integer(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard(['1b']), 0, 'field')

        # blank
        with self.assertRaises(SyntaxError):
            integer(BDFCard(['']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard([None]), 0, 'field')

        with self.assertRaises(SyntaxError):
            integer(BDFCard(['7.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            integer(BDFCard(['1+2']), 0, 'field')

        card = [1, '2', '3.', 'C', None, '']
        exact = [1, 2, SyntaxError, SyntaxError, SyntaxError, SyntaxError]
        self.run_function(integer, card, exact)

    def test_string(self):
        """
        value = string(card, n, fieldname)
        """
        # integer
        with self.assertRaises(SyntaxError):
            string(BDFCard([1]), 0, 'field')
        with self.assertRaises(SyntaxError):
            string(BDFCard(['1']), 0, 'field')

        # float
        with self.assertRaises(SyntaxError):
            string(BDFCard([1.]), 0, 'field')
        with self.assertRaises(SyntaxError):
            string(BDFCard(['1.']), 0, 'field')

        # string
        self.assertEqual('A', string(BDFCard(['a']), 0, 'field'))
        with self.assertRaises(SyntaxError):
            self.assertEqual('1B', string(BDFCard(['1b']), 0, 'field'))
        self.assertEqual('CAT', string(BDFCard([' cat ']), 0, 'field'))

        # blank
        with self.assertRaises(SyntaxError):
            string(BDFCard(['']), 0, 'field')
        with self.assertRaises(SyntaxError):
            string(BDFCard([None]), 0, 'field')

        with self.assertRaises(SyntaxError):
            string(BDFCard(['cat dog']), 0, 'field')

        with self.assertRaises(SyntaxError):
            string(BDFCard(['1+2']), 0, 'field')

        card = ['3.', None, '']
        exact = [SyntaxError, SyntaxError, SyntaxError]
        self.run_function(string, card, exact)

        self.check_string(string)

    def test_string_or_blank(self):
        with self.assertRaises(SyntaxError):
            string_or_blank(BDFCard(['-CAT']), 0, 'field')
        with self.assertRaises(SyntaxError):
            string_or_blank(BDFCard(['1+2']), 0, 'field')

        self.assertEqual('CAT', string_or_blank(BDFCard(['cat']), 0, 'field'))
        self.assertEqual(None, string_or_blank(BDFCard(['  ']), 0, 'field'))

        with self.assertRaises(SyntaxError):
            self.assertEqual(100, string_or_blank(BDFCard(['100']), 0, 'field'))

        with self.assertRaises(SyntaxError):
            string_or_blank(BDFCard(['1 2']), 0, 'field')
        with self.assertRaises(SyntaxError):
            string_or_blank(BDFCard(['c a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            string_or_blank(BDFCard(['1b']), 0, 'field')

        self.check_string(string_or_blank)
        self.check_blank(string_or_blank)

    #-------------------------------------------------------
    def test_integer_or_blank(self):
        """
        value = integer_or_blank(card, n, fieldname)
        """
        # integer
        self.check_integer(integer_or_blank)

        # float
        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard([1.]), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard(['1.']), 0, 'field')
        #with self.assertRaises(SyntaxError):
            #integer(BDFCard(['cat']), 0, 'field')

        # string
        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard(['1b']), 0, 'field')

        self.check_blank(integer_or_blank)

        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard(['1+2']), 0, 'field')

        card = [1, '2', '3.', 'C', None, '']
        exact = [1, 2, SyntaxError, SyntaxError, None, None]
        default = [None, None, None, None, None, None]
        self.run_function_default(integer_or_blank, card, exact, default)

    def check_integer(self, method):
        """common integer checks"""
        self.assertEqual(1, method(BDFCard([1]), 0, 'field'))
        self.assertEqual(2, method(BDFCard(['2']), 0, 'field'))
        self.assertEqual(-1, method(BDFCard(['-1']), 0, 'field'))

        #if check_space:
        with self.assertRaises(SyntaxError):
            method(BDFCard(['1 3']), 0, 'field')
        with self.assertRaises(SyntaxError):
            method(BDFCard(['-1 3']), 0, 'field')

    def check_double(self, method, is_strict: bool=True):
        """common double checks"""
        method(BDFCard([3.0]), 0, 'field')
        method(BDFCard(['4.0']), 0, 'field')
        method(BDFCard(['5.']), 0, 'field')

        self.assertEqual(1.0, method(BDFCard([1.]), 0, 'field'))
        self.assertEqual(1.0, method(BDFCard(['1.']), 0, 'field'))
        self.assertEqual(-9.31e-4, method(BDFCard(['-9.31-4']), 0, 'field'))

        # float
        val = method(BDFCard([1.]), 0, 'field')
        self.assertEqual(1., val)
        val = method(BDFCard(['1.']), 0, 'field')
        self.assertEqual(1., val)
        val = method(BDFCard(['1-3']), 0, 'field')
        self.assertEqual(1.e-3, val)

        self.assertEqual(1.e-9, method(BDFCard(['1-9']), 0, 'field'))
        self.assertEqual(1.e+9, method(BDFCard(['1+9']), 0, 'field'))

        self.assertEqual(1.e-9, method(BDFCard(['1.e-9']), 0, 'field'))
        self.assertEqual(1.e+9, method(BDFCard(['1.e+9']), 0, 'field'))
        self.assertEqual(1.e-9, method(BDFCard(['1.d-9']), 0, 'field'))
        self.assertEqual(1.e+9, method(BDFCard(['1.d+9']), 0, 'field'))

        self.assertEqual(1.e-9, method(BDFCard(['1.E-9']), 0, 'field'))
        self.assertEqual(1.e+9, method(BDFCard(['1.E+9']), 0, 'field'))
        self.assertEqual(1.e-9, method(BDFCard(['1.D-9']), 0, 'field'))
        self.assertEqual(1.e+9, method(BDFCard(['1.D+9']), 0, 'field'))

        # D-notation without explicit sign on exponent
        self.assertEqual(1.e3, method(BDFCard(['1.D3']), 0, 'field'))
        self.assertEqual(1.e3, method(BDFCard(['1.d3']), 0, 'field'))
        self.assertEqual(5.e2, method(BDFCard(['.5D3']), 0, 'field'))

        # E-notation without explicit sign on exponent
        self.assertEqual(1.e3, method(BDFCard(['1.E3']), 0, 'field'))

        # implicit exponent with leading dot
        self.assertEqual(5.e2, method(BDFCard(['.5+3']), 0, 'field'))
        self.assertEqual(-5.e2, method(BDFCard(['-.5+3']), 0, 'field'))
        self.assertEqual(5.e-4, method(BDFCard(['+.5-3']), 0, 'field'))

        #if check_space:
        if is_strict:
            with self.assertRaises(SyntaxError):
                method(BDFCard(['-9. 31-4']), 0, 'field')
        else:
            self.assertEqual(method(BDFCard(['-9. 31-4']), 0, 'field'), -9.31e-4)


    def check_string(self, method, check_dash=True):
        """common string checks"""
        self.assertEqual('A', method(BDFCard(['a']), 0, 'field'))
        self.assertEqual('B1', method(BDFCard(['b1']), 0, 'field'))
        self.assertEqual('C', method(BDFCard(['C']), 0, 'field'))
        self.assertEqual('FROG', method(BDFCard(['FROG']), 0, 'field'))

        if check_dash:
            self.assertEqual('VONE-MIS', method(BDFCard(['VONE-MIS']), 0, 'field'))
        self.assertEqual('VONE_MIS', method(BDFCard(['VONE_MIS']), 0, 'field'))

        with self.assertRaises(SyntaxError):
            method(BDFCard(['VON MISES']), 0, 'field')
        #with self.assertRaises(Syntt)

    def check_blank(self, method):
        """common blank checks"""
        assert method(BDFCard(['   ']), 0, 'field') is None
        assert method(BDFCard([None]), 0, 'field') is None
        #assert method(BDFCard(['.']), 0, 'field') == 0.
        self.assertEqual('a', method(BDFCard(['']), 0, 'field', 'a'))
        self.assertEqual('b', method(BDFCard([None]), 0, 'field', 'b'))

    def test_double_or_blank(self):
        """
        value = double_or_blank(card, n, fieldname, default=None)
        """
        # integer
        card = BDFCard([1])
        with self.assertRaises(SyntaxError):
            double_or_blank(card, 0, 'field')
        card = BDFCard(['2'])
        with self.assertRaises(SyntaxError):
            double_or_blank(card, 0, 'field')
        # signed integer strings must also be rejected
        with self.assertRaises(SyntaxError):
            double_or_blank(BDFCard(['-400510']), 0, 'field')
        with self.assertRaises(SyntaxError):
            double_or_blank(BDFCard(['+400510']), 0, 'field')

        self.check_double(double_or_blank)

        # string
        with self.assertRaises(SyntaxError):
            double_or_blank(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            double_or_blank(BDFCard(['1b']), 0, 'field')

        # blank
        self.check_blank(double_or_blank)
        assert double_or_blank(BDFCard(['.']), 0, 'field') == 0.

    def test_double_or_string(self):
        """tests the double_or_string function"""
        # out of range
        with self.assertRaises(SyntaxError):
            double_or_string(BDFCard([1.]), 1, 'field')
        self.assertEqual(1.e-9, double_or_string(BDFCard(['1-9']), 0, 'field'))
        self.assertEqual(1.e+9, double_or_string(BDFCard(['1+9']), 0, 'field'))

        self.check_double(double_or_string)
        self.check_string(double_or_string, check_dash=False)

    def test_double_string_or_blank(self):
        """tests the double_string_or_blank function"""
        # out of range
        #with self.assertRaises(SyntaxError):
        self.assertEqual(1., double_string_or_blank(BDFCard([1.]), 0, 'field'))
        self.assertEqual(1., double_string_or_blank(BDFCard(['1.']), 0, 'field'))
        self.assertEqual('CAT', double_string_or_blank(BDFCard(['CAT']), 0, 'field'))
        self.assertEqual('CAT', double_string_or_blank(BDFCard([' CAT ']), 0, 'field'))
        self.assertEqual(None, double_string_or_blank(BDFCard([None]), 0, 'field'))
        self.assertEqual(1.e-9, double_string_or_blank(BDFCard(['1-9']), 0, 'field'))
        self.assertEqual(1.e+9, double_string_or_blank(BDFCard(['1+9']), 0, 'field'))

        with self.assertRaises(SyntaxError):
            double_string_or_blank(BDFCard(['10']), 0, 'field')

        self.assertEqual(double_string_or_blank(BDFCard(['-1.0']), 0, 'field'), -1)
        self.assertEqual(double_string_or_blank(BDFCard(['cat']), 0, 'field'), 'CAT')
        self.assertEqual(double_string_or_blank(BDFCard(['  ']), 0, 'field'), None)

        self.check_double(double_string_or_blank)
        self.check_string(double_string_or_blank, check_dash=False)
        self.check_blank(double_string_or_blank)

    def test_integer_or_double(self):
        """tests the integer_or_double function"""
        integer_or_double(BDFCard([1.], has_none=False), 0, 'field')
        integer_or_double(BDFCard([1], has_none=False), 0, 'field')

        # out of range
        with self.assertRaises(SyntaxError):
            integer_or_double(BDFCard([1.]), 1, 'field')
        self.assertEqual(1.e-9, integer_or_double(BDFCard(['1-9']), 0, 'field'))
        self.assertEqual(1.e+9, integer_or_double(BDFCard(['1+9']), 0, 'field'))
        self.assertEqual(-1, integer_or_double(BDFCard(['-1']), 0, 'field'))

        with self.assertRaises(SyntaxError):
            integer_or_double(BDFCard(['cat']), 0, 'field')

        with self.assertRaises(SyntaxError):
            integer_or_double(BDFCard(['.']), 0, 'field')

        self.assertEqual(100., integer_or_double(BDFCard(['1+2']), 0, 'field'))

    def test_integer_or_string(self):
        """tests the integer_or_string function"""
        with self.assertRaises(SyntaxError):
            integer(BDFCard([1.], has_none=False), 0, 'field')
        integer(BDFCard([1], has_none=False), 0, 'field')

        # out of range
        with self.assertRaises(SyntaxError):
            integer_or_string(BDFCard([1.]), 1, 'field')

        self.assertEqual(1000, integer_or_string(BDFCard([1000]), 0, 'field'))
        self.assertEqual(-1, integer_or_string(BDFCard(['-1']), 0, 'field'))
        self.assertEqual(1000, integer_or_string(BDFCard(['1000']), 0, 'field'))
        self.assertEqual('CAT', integer_or_string(BDFCard(['cat']), 0, 'field'))
        self.assertEqual('CAT', integer_or_string(BDFCard([' cat ']), 0, 'field'))

        with self.assertRaises(SyntaxError):
            integer_or_string(BDFCard(['1b']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_or_string(BDFCard(['1+2']), 0, 'field')


    def test_integer_double_or_blank(self):
        """
        value = double_or_blank(card, n, fieldname, default=None)
        """
        with self.assertRaises(SyntaxError):
            integer(BDFCard([1.], has_none=False), 0, 'field')
        integer(BDFCard([1], has_none=False), 0, 'field')

        # integer
        self.check_integer(integer_double_or_blank)

        # float
        self.check_double(integer_double_or_blank)
        self.assertEqual(1.e-9, integer_double_or_blank(BDFCard(['1-9']), 0, 'field'))
        self.assertEqual(1.e+9, integer_double_or_blank(BDFCard(['1+9']), 0, 'field'))
        self.assertEqual(-1, integer_double_or_blank(BDFCard(['-1']), 0, 'field'))

        # error - string
        with self.assertRaises(SyntaxError):
            integer_double_or_blank(BDFCard(['C']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_double_or_blank(BDFCard(['1C']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_double_or_blank(BDFCard(['C1']), 0, 'field')

        # blank
        self.check_blank(integer_double_or_blank)
        self.assertEqual(1.e-9, double_or_blank(BDFCard(['1-9']), 0, 'field'))
        self.assertEqual(1.e+9, double_or_blank(BDFCard(['1+9']), 0, 'field'))

        self.assertEqual(1000, integer_double_or_blank(BDFCard(['1000']), 0, 'field'))

        #card    = [1,    2.0, '3.0', '4.', 'C',        None, None,          '', None, 'cat']
        #exact   = [1,    2.0,  3.0,   4.0, SyntaxError,None, 2.0,  SyntaxError, 1.0, SyntaxError]
        #default = [None, None, None, None, None,       None, 2.0,         None, 1.0, 1.0]
        #self.run_function_default(integer_double_or_blank, card, exact, default)

    def test_integer_string_or_blank(self):
        """tests the integer_string_or_blank function"""
        with self.assertRaises(SyntaxError):
            integer(BDFCard([1.], has_none=False), 0, 'field')
        integer(BDFCard([1], has_none=False), 0, 'field')

        # integer
        self.check_integer(integer_string_or_blank)

        # float
        #print type(integer_string_or_blank(BDFCard(['4.0']), 0, 'field'))

        with self.assertRaises(SyntaxError):
            integer_string_or_blank(BDFCard([3.0]), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_string_or_blank(BDFCard(['4.0']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_string_or_blank(BDFCard(['5.']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_string_or_blank(BDFCard(['1b']), 0, 'field')

        # string
        self.assertEqual('LOAD', integer_string_or_blank(BDFCard(['load']), 0, 'field'))
        self.assertEqual(1000, integer_string_or_blank(BDFCard([1000]), 0, 'field'))
        self.assertEqual(-1, integer_string_or_blank(BDFCard(['-1']), 0, 'field'))

        # blank
        self.check_blank(integer_string_or_blank)


    #def test_double_string_or_blank(self):

    def test_integer_double_or_string(self):
        """tests the integer_double_or_string function"""
        with self.assertRaises(SyntaxError):
            integer(BDFCard([1.], has_none=False), 0, 'field')
        integer(BDFCard([1], has_none=False), 0, 'field')

        # out of range
        with self.assertRaises(SyntaxError):
            integer_double_or_string(BDFCard([1.]), 1, 'field')
        with self.assertRaises(SyntaxError):
            integer_double_or_string(BDFCard(['1b']), 0, 'field')

        # integer
        self.check_integer(integer_double_or_string)

        # float
        self.check_double(integer_double_or_string)

        # string
        self.assertEqual('LOAD', integer_double_or_string(BDFCard(['load']), 0, 'field'))
        #self.assertEqual('MN-MM', integer_double_or_string(BDFCard(['MN-MM']), 0, 'field'))
        #self.assertEqual(-1, integer_double_or_string(BDFCard(['-1']), 0, 'field'))
        #self.assertEqual(1000, integer_double_or_string(BDFCard([1000]), 0, 'field'))

    def test_integer_double_string_or_blank(self):
        """tests the integer_double_string_or_blank function"""
        integer_double_string_or_blank(BDFCard([1.], has_none=False), 0, 'field')
        integer_double_string_or_blank(BDFCard([1], has_none=False), 0, 'field')

        # out of range
        self.assertEqual(None, integer_double_string_or_blank(BDFCard([1.]), 1, 'field'))
        #with self.assertRaises(SyntaxError):
            #print(integer_double_string_or_blank(BDFCard(['1b']), 0, 'field'))

        # integer
        self.check_integer(integer_double_string_or_blank)

        # float
        self.check_double(integer_double_string_or_blank)

        # string
        self.assertEqual('LOAD', integer_double_string_or_blank(BDFCard(['load']), 0, 'field'))
        self.assertEqual('MN-MM', integer_double_string_or_blank(BDFCard(['MN-MM']), 0, 'field'))
        #self.assertEqual(-1, integer_double_string_or_blank(BDFCard(['-1']), 0, 'field'))
        self.assertEqual(1000, integer_double_string_or_blank(BDFCard([1000]), 0, 'field'))
        self.assertEqual('CAT', integer_double_string_or_blank(BDFCard([100]), 1, 'field', 'CAT'))

        # int/float
        self.assertEqual(1000, integer_double_string_or_blank(BDFCard([1000]), 0, 'field'))
        self.assertEqual(1000.0, integer_double_string_or_blank(BDFCard([1000.0]), 0, 'field'))

        # blank
        #self.check_blank(integer_double_string_or_blank)

    def test_components(self):
        """tests the parse_components function"""
        with self.assertRaises(SyntaxError):
            integer_string_or_blank(BDFCard(['1b']), 0, 'field')

        parse_components(BDFCard([1], has_none=False), 0, 'field')

        # single ints
        val = parse_components(BDFCard([0]), 0, 'field')
        self.assertEqual(val, '0')

        val = parse_components(BDFCard([1]), 0, 'field')
        self.assertEqual(val, '1')

        # single strings
        val = parse_components(BDFCard(['0']), 0, 'field')
        self.assertEqual(val, '0')

        val = parse_components(BDFCard(['1']), 0, 'field')
        self.assertEqual(val, '1')

        # double ints
        val = parse_components(BDFCard(['123']), 0, 'field')
        self.assertEqual(val, '123')

        val = parse_components(BDFCard([123]), 0, 'field')
        self.assertEqual(val, '123')

        val = parse_components(BDFCard([321]), 0, 'field')
        self.assertEqual(val, '123')

        # embedded whiteshape
        with self.assertRaises(SyntaxError):
            parse_components(BDFCard(['12 3']), 0, 'field')

        # all numbers
        val = parse_components(BDFCard(['123456']), 0, 'field')
        self.assertEqual(val, '123456')

        # invalid 0's defined with numbers
        with self.assertRaises(SyntaxError):
            parse_components(BDFCard(['0123456']), 0, 'field')

        with self.assertRaises(SyntaxError):
            parse_components(BDFCard(['01']), 0, 'field')

        # doubles
        with self.assertRaises(SyntaxError):
            parse_components(BDFCard(['4524']), 0, 'field')

        # only 0 to 6
        with self.assertRaises(SyntaxError):
            parse_components(BDFCard(['7']), 0, 'field')

        with self.assertRaises(SyntaxError):
            parse_components(BDFCard([7]), 0, 'field')

        # dumb input
        with self.assertRaises(SyntaxError):
            parse_components(BDFCard(['4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            parse_components(BDFCard(['-4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            parse_components(BDFCard(['asdf']), 0, 'field')

        with self.assertRaises(SyntaxError):
            parse_components(BDFCard(['-1']), 0, 'field')

        # blank
        #parse_components(BDFCard(['   ']), 0, 'field')
        #parse_components(BDFCard([None]), 0, 'field')

    def test_components_or_blank_01(self):
        """tests the components_or_blank function"""
        # single ints
        val = components_or_blank(BDFCard([0]), 0, 'field')
        self.assertEqual(val, '0')

        val = components_or_blank(BDFCard([1]), 0, 'field')
        self.assertEqual(val, '1')

        # single strings
        val = components_or_blank(BDFCard(['0']), 0, 'field')
        self.assertEqual(val, '0')

        val = components_or_blank(BDFCard(['1']), 0, 'field')
        self.assertEqual(val, '1')

        # double ints
        val = components_or_blank(BDFCard(['123']), 0, 'field')
        self.assertEqual(val, '123')

        val = components_or_blank(BDFCard([123]), 0, 'field')
        self.assertEqual(val, '123')

        val = components_or_blank(BDFCard([321]), 0, 'field')
        self.assertEqual(val, '123')

        # embedded whiteshape
        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['12 3']), 0, 'field')

        # all numbers
        val = components_or_blank(BDFCard(['123456']), 0, 'field')
        self.assertEqual(val, '123456')

        # invalid 0's defined with numbers
        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['0123456']), 0, 'field')

        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['01']), 0, 'field')

        # doubles
        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['4524']), 0, 'field')

        # only 0 to 6
        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['7']), 0, 'field')

        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard([7]), 0, 'field')

        # dumb input
        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['-4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['asdf']), 0, 'field')

        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['-1']), 0, 'field')

        # blank
        #components_or_blank(BDFCard(['   ']), 0, 'field')
        #components_or_blank(BDFCard([None]), 0, 'field')

    def test_components_or_blank_02(self):
        # single ints
        val = components_or_blank(BDFCard([0]), 0, 'field', 'default')
        self.assertEqual(val, '0')

        val = components_or_blank(BDFCard([1]), 0, 'field', 'default')
        self.assertEqual(val, '1')

        # single strings
        val = components_or_blank(BDFCard(['0']), 0, 'field', 'default')
        self.assertEqual(val, '0')

        val = components_or_blank(BDFCard(['1']), 0, 'field', 'default')
        self.assertEqual(val, '1')

        # double ints
        val = components_or_blank(BDFCard(['123']), 0, 'field', 'default')
        self.assertEqual(val, '123')

        val = components_or_blank(BDFCard([123]), 0, 'field', 'default')
        self.assertEqual(val, '123')

        val = components_or_blank(BDFCard([321]), 0, 'field', 'default')
        self.assertEqual(val, '123')


        # embedded whiteshape
        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['12 3']), 0, 'field', 'default')

        # all numbers
        val = components_or_blank(BDFCard(['123456']), 0, 'field', 'default')
        self.assertEqual(val, '123456')

        # invalid 0's defined with numbers
        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['0123456']), 0, 'field', 'default')

        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['01']), 0, 'field', 'default')

        # doubles
        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['4524']), 0, 'field', 'default')

        # only 0 to 6
        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['7']), 0, 'field', 'default')

        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard([7]), 0, 'field', 'default')

        # dumb input
        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['4.0']), 0, 'field', 'default')

        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['-4.0']), 0, 'field', 'default')

        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['asdf']), 0, 'field', 'default')

        with self.assertRaises(SyntaxError):
            components_or_blank(BDFCard(['-1']), 0, 'field', 'default')

        # blank
        val = components_or_blank(BDFCard(['   ']), 0, 'field', 'default')
        self.assertEqual(val, 'default')
        val = components_or_blank(BDFCard([None]), 0, 'field', 'default')
        self.assertEqual(val, 'default')

    def test_bad(self):
        _get_dtype('1.000000000D+00')
        interpret_value('1.000000000D+00')

    def test_modal_components(self):
        """modal components"""
        card = BDFCard(['42'])

        with self.assertRaises(SyntaxError):
            modal_components(card, 0, 'field')

        self.assertEqual(modal_components(BDFCard(['-1']), 0, 'field'), -1)
        self.assertEqual(modal_components(BDFCard(['0']), 0, 'field'), 0)
        self.assertEqual(modal_components(BDFCard(['1']), 0, 'field'), 1)
        self.assertEqual(modal_components(BDFCard(['2']), 0, 'field'), 2)
        self.assertEqual(modal_components(BDFCard(['3']), 0, 'field'), 3)
        self.assertEqual(modal_components(BDFCard(['4']), 0, 'field'), 4)
        self.assertEqual(modal_components(BDFCard(['5']), 0, 'field'), 5)
        self.assertEqual(modal_components(BDFCard(['6']), 0, 'field'), 6)
        with self.assertRaises(SyntaxError):
            self.assertEqual(modal_components(BDFCard(['7']), 0, 'field'), 7)

    def test_interpret_value(self):
        """Test interpret_value covers all parsing branches:
        None, blank, int, float, string, D-notation, implicit exponent, special chars.
        """
        # None input
        self.assertIsNone(interpret_value(None))

        # blank / whitespace
        self.assertIsNone(interpret_value(''))
        self.assertIsNone(interpret_value('   '))
        self.assertIsNone(interpret_value('  *'))

        # integers
        self.assertEqual(interpret_value('1'), 1)
        self.assertEqual(interpret_value('-1'), -1)
        self.assertEqual(interpret_value('+42'), 42)
        self.assertEqual(interpret_value('  100  '), 100)
        self.assertEqual(interpret_value('0'), 0)

        # floats (standard notation)
        self.assertEqual(interpret_value('1.0'), 1.0)
        self.assertEqual(interpret_value('-3.14'), -3.14)
        self.assertEqual(interpret_value('1.5E+3'), 1500.0)
        self.assertEqual(interpret_value('1.5E-3'), 0.0015)
        self.assertEqual(interpret_value('.5'), 0.5)

        # D-notation (Nastran double precision)
        self.assertEqual(interpret_value('1.000000000D+00'), 1.0)
        self.assertEqual(interpret_value('2.5D+3'), 2500.0)
        self.assertEqual(interpret_value('-1.5D-2'), -0.015)
        self.assertEqual(interpret_value('1.D+0'), 1.0)

        # implicit exponent notation (Nastran shorthand: 1.5+3 means 1.5e3)
        self.assertEqual(interpret_value('1.5+3'), 1500.0)
        self.assertEqual(interpret_value('1.5-3'), 0.0015)
        self.assertAlmostEqual(interpret_value('-5.007-3') / -5.007e-3, 1.0, places=10)
        self.assertAlmostEqual(interpret_value('8.182-18') / 8.182e-18, 1.0, places=10)
        self.assertEqual(interpret_value('+3.0+2'), 300.0)

        # strings (start with alpha)
        self.assertEqual(interpret_value('GRID'), 'GRID')
        self.assertEqual(interpret_value('  mat1  '), 'MAT1')

        # special characters (=, (, *) returned as-is stripped
        self.assertEqual(interpret_value('='), '=')
        self.assertEqual(interpret_value('(1.0)'), '(1.0)')
        self.assertEqual(interpret_value('2*3.0'), '2*3.0')

        # already numeric input (int/float pass-through)
        self.assertEqual(interpret_value(42), 42)
        self.assertEqual(interpret_value(3.14), 3.14)

        # error case: can't find exponent sign
        with self.assertRaises(SyntaxError):
            interpret_value('1.5X3')

    def test_bdf_card_field(self):
        """Test BDFCard.field(i, default) for in-range, out-of-range, blank, and None fields."""
        card = BDFCard(['GRID', '1', '', '1.0', '2.0', '3.0'], has_none=False)

        # in-range, populated field
        self.assertEqual(card.field(0), 'GRID')
        self.assertEqual(card.field(1), '1')
        self.assertEqual(card.field(3), '1.0')

        # in-range, blank field -> returns default
        self.assertEqual(card.field(2), None)
        self.assertEqual(card.field(2, 'mydefault'), 'mydefault')

        # out-of-range -> returns default
        self.assertEqual(card.field(99), None)
        self.assertEqual(card.field(99, 42), 42)
        self.assertEqual(card.field(-100, 'fallback'), 'fallback')

        # None field (has_none=True path)
        card2 = BDFCard(['GRID', '1', None, '1.0'], has_none=True)
        self.assertEqual(card2.field(2), None)
        self.assertEqual(card2.field(2, 0), 0)

    def test_bdf_card_repr(self):
        """Test BDFCard.__repr__ returns a parseable list representation."""
        card = BDFCard(['GRID', '1', None, '1.0'], has_none=True)
        result = repr(card)
        # repr uses %r formatting of the internal card list
        self.assertIn('GRID', result)
        self.assertIn('1.0', result)
        # should be a valid Python list literal
        self.assertTrue(result.startswith('['))
        self.assertTrue(result.endswith(']'))

    def test_bdf_card_fields(self):
        """Test BDFCard.fields(i, j) for slicing subsets of card fields."""
        card = BDFCard(['GRID', '1', '', '10.0', '20.0', '30.0'], has_none=False)

        # fields(0) -> all fields
        all_fields = card.fields(0)
        self.assertEqual(len(all_fields), 6)
        self.assertEqual(all_fields[0], 'GRID')
        self.assertEqual(all_fields[5], '30.0')

        # fields(i, j) -> subset [i, j)
        subset = card.fields(3, 6)
        self.assertEqual(subset, ['10.0', '20.0', '30.0'])

        # fields(1, 3) -> indices 1 and 2
        subset2 = card.fields(1, 3)
        self.assertEqual(subset2[0], '1')
        # index 2 is blank ('') -> field() returns None
        self.assertEqual(subset2[1], None)

        # fields starting past end -> all None (default)
        subset3 = card.fields(10, 12)
        self.assertEqual(subset3, [None, None])

        # fields with j=None -> to end of card
        subset4 = card.fields(4)
        self.assertEqual(len(subset4), 2)
        self.assertEqual(subset4, ['20.0', '30.0'])


    def test_force_double_or_blank(self):
        """Test force_double_or_blank: coerces ints to floats with warning, parses Nastran floats."""
        import warnings

        # float passthrough
        self.assertEqual(force_double_or_blank(BDFCard([3.14]), 0, 'f'), 3.14)

        # integer coerced to float (with warning)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.assertEqual(force_double_or_blank(BDFCard([5]), 0, 'f'), 5.0)

        # string integer coerced to float (with warning)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.assertEqual(force_double_or_blank(BDFCard(['10']), 0, 'f'), 10.0)

        # Nastran float notation
        self.assertEqual(force_double_or_blank(BDFCard(['1.5']), 0, 'f'), 1.5)
        self.assertEqual(force_double_or_blank(BDFCard(['1.-3']), 0, 'f'), 1.0e-3)
        self.assertEqual(force_double_or_blank(BDFCard(['1.+3']), 0, 'f'), 1.0e3)

        # blank -> default
        self.assertIsNone(force_double_or_blank(BDFCard(['']), 0, 'f'))
        self.assertIsNone(force_double_or_blank(BDFCard([None]), 0, 'f'))
        self.assertEqual(force_double_or_blank(BDFCard(['']), 0, 'f', default=99.0), 99.0)

        # '.' -> 0.0
        self.assertEqual(force_double_or_blank(BDFCard(['.']), 0, 'f'), 0.0)

        # out of range -> default
        self.assertIsNone(force_double_or_blank(BDFCard([1.0]), 5, 'f'))
        self.assertEqual(force_double_or_blank(BDFCard([1.0]), 5, 'f', default=7.0), 7.0)

        # invalid string -> SyntaxError
        with self.assertRaises(SyntaxError):
            force_double_or_blank(BDFCard(['abc']), 0, 'f')

    def test_lax_double_or_blank(self):
        """Test lax_double_or_blank: wraps integer_double_or_blank, coerces int to float."""
        # float passthrough
        self.assertEqual(lax_double_or_blank(BDFCard([3.14]), 0, 'f'), 3.14)

        # integer coerced to float
        self.assertEqual(lax_double_or_blank(BDFCard([5]), 0, 'f'), 5.0)
        self.assertIsInstance(lax_double_or_blank(BDFCard([5]), 0, 'f'), float)

        # string float
        self.assertEqual(lax_double_or_blank(BDFCard(['1.5']), 0, 'f'), 1.5)

        # string integer -> int -> float
        self.assertEqual(lax_double_or_blank(BDFCard(['10']), 0, 'f'), 10.0)
        self.assertIsInstance(lax_double_or_blank(BDFCard(['10']), 0, 'f'), float)

        # blank -> default
        self.assertIsNone(lax_double_or_blank(BDFCard([None]), 0, 'f'))
        self.assertEqual(lax_double_or_blank(BDFCard([None]), 0, 'f', default=2.0), 2.0)

        # Nastran notation
        self.assertEqual(lax_double_or_blank(BDFCard(['1.-3']), 0, 'f'), 1e-3)

    def test_string_choice_or_blank(self):
        """Test string_choice_or_blank: validates string is in allowed set."""
        choices = ('YES', 'NO', 'MAYBE')

        # valid choice
        self.assertEqual(string_choice_or_blank(BDFCard(['yes']), 0, 'f', choices), 'YES')
        self.assertEqual(string_choice_or_blank(BDFCard(['NO']), 0, 'f', choices), 'NO')

        # blank -> default
        self.assertIsNone(string_choice_or_blank(BDFCard([None]), 0, 'f', choices))
        self.assertEqual(string_choice_or_blank(BDFCard([None]), 0, 'f', choices, default='YES'), 'YES')
        self.assertIsNone(string_choice_or_blank(BDFCard(['']), 0, 'f', choices))

        # invalid choice -> RuntimeError
        with self.assertRaises(RuntimeError):
            string_choice_or_blank(BDFCard(['INVALID']), 0, 'f', choices)

        # non-string input -> SyntaxError
        with self.assertRaises(SyntaxError):
            string_choice_or_blank(BDFCard([1.0], has_none=False), 0, 'f', choices)

    def test_modal_components_or_blank(self):
        """Test modal_components_or_blank: integer -1 to 6, or blank."""
        # valid values
        self.assertEqual(modal_components_or_blank(BDFCard(['-1']), 0, 'f', 0), -1)
        self.assertEqual(modal_components_or_blank(BDFCard(['0']), 0, 'f', 0), 0)
        self.assertEqual(modal_components_or_blank(BDFCard(['3']), 0, 'f', 0), 3)
        self.assertEqual(modal_components_or_blank(BDFCard(['6']), 0, 'f', 0), 6)

        # blank -> default
        self.assertEqual(modal_components_or_blank(BDFCard([None]), 0, 'f', 0), 0)
        self.assertEqual(modal_components_or_blank(BDFCard(['']), 0, 'f', 3), 3)

        # out of range -> SyntaxError
        with self.assertRaises(SyntaxError):
            modal_components_or_blank(BDFCard(['7']), 0, 'f', 0)
        with self.assertRaises(SyntaxError):
            modal_components_or_blank(BDFCard(['-2']), 0, 'f', 0)

    def test_exact_string_or_blank(self):
        """Test exact_string_or_blank: returns left-justified 8-char string or default."""
        # populated field -> 8-char left-justified
        result = exact_string_or_blank(BDFCard(['GRID'], has_none=False), 0, 'f')
        self.assertEqual(result, 'GRID    ')
        self.assertEqual(len(result), 8)

        # blank -> default
        self.assertIsNone(exact_string_or_blank(BDFCard([None]), 0, 'f'))
        self.assertEqual(exact_string_or_blank(BDFCard([None]), 0, 'f', default='DEF'), 'DEF')

        # short string -> padded to 8
        result = exact_string_or_blank(BDFCard(['AB'], has_none=False), 0, 'f')
        self.assertEqual(result, 'AB      ')
        self.assertEqual(len(result), 8)

    def test_filename_or_blank(self):
        """Test filename_or_blank: validates string as filename (no digits first, no +/- first)."""
        # valid filename
        self.assertEqual(filename_or_blank(BDFCard(['MYFILE'], has_none=False), 0, 'f'), 'MYFILE')
        self.assertEqual(filename_or_blank(BDFCard(['model'], has_none=False), 0, 'f'), 'MODEL')

        # blank -> default
        self.assertIsNone(filename_or_blank(BDFCard([None]), 0, 'f'))
        self.assertEqual(filename_or_blank(BDFCard([None]), 0, 'f', default='X'), 'X')

        # starts with digit -> SyntaxError
        with self.assertRaises(SyntaxError):
            filename_or_blank(BDFCard(['1FILE'], has_none=False), 0, 'f')

        # contains + -> SyntaxError
        with self.assertRaises(SyntaxError):
            filename_or_blank(BDFCard(['A+B'], has_none=False), 0, 'f')

        # starts with - -> SyntaxError
        with self.assertRaises(SyntaxError):
            filename_or_blank(BDFCard(['-FILE'], has_none=False), 0, 'f')

        # non-string -> SyntaxError
        with self.assertRaises(SyntaxError):
            filename_or_blank(BDFCard([1.0], has_none=False), 0, 'f')

        # space in string -> SyntaxError
        with self.assertRaises(SyntaxError):
            filename_or_blank(BDFCard(['MY FILE'], has_none=False), 0, 'f')

    def test_check_string(self):
        """Test check_string: validates string (no spaces, can't start with digit/./+/-)."""
        # valid strings
        self.assertEqual(check_string('GRID', 0, 'f'), 'GRID')
        self.assertEqual(check_string('mat1', 0, 'f'), 'MAT1')
        self.assertEqual(check_string('ABC', 0, 'f'), 'ABC')

        # space in string -> SyntaxError
        with self.assertRaises(SyntaxError):
            check_string('AB CD', 0, 'f')

        # starts with digit -> SyntaxError
        with self.assertRaises(SyntaxError):
            check_string('1ABC', 0, 'f')

        # starts with - -> SyntaxError
        with self.assertRaises(SyntaxError):
            check_string('-ABC', 0, 'f')

        # starts with + -> SyntaxError
        with self.assertRaises(SyntaxError):
            check_string('+ABC', 0, 'f')

        # contains . -> SyntaxError
        with self.assertRaises(SyntaxError):
            check_string('A.B', 0, 'f')

        # non-string input -> SyntaxError
        with self.assertRaises(SyntaxError):
            check_string(1.0, 0, 'f')

    def test_loose_string(self):
        """Test loose_string: lenient string check, only rejects leading digit."""
        # valid (starts with alpha) -> returns uppercased string
        result = loose_string(BDFCard(['LABEL'], has_none=False), 0, 'f')
        self.assertEqual(result, 'LABEL')

        result = loose_string(BDFCard(['myLabel'], has_none=False), 0, 'f')
        self.assertEqual(result, 'MYLABEL')

        # blank -> default
        self.assertIsNone(loose_string(BDFCard([None]), 0, 'f'))
        self.assertEqual(loose_string(BDFCard([None]), 0, 'f', default='X'), 'X')

        # starts with digit -> SyntaxError
        with self.assertRaises(SyntaxError):
            loose_string(BDFCard(['1LABEL'], has_none=False), 0, 'f')


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
