"""tests various assign_types.py functions"""
import unittest
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank,
    double, double_or_blank, integer_or_double, integer_double_or_blank,
    string, string_or_blank, double_or_string, double_string_or_blank,
    integer_or_string, integer_string_or_blank, integer_double_or_string,
    blank, parse_components, components_or_blank, integer_double_string_or_blank,
    _get_dtype, interpret_value, modal_components)


class TestAssignType(unittest.TestCase):
    """tests various assign_types.py functions"""

    def run_function_default(self, func, card, exact, default):
        """
        Helper function

        Parameters
        ----------
        func : function
           integer_or_blank
        card : List[varies]
            a series of values to add
        exacts : List[float]
            list of results
        default : List[float]
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
        card : List[varies]
            a series of values to add
        exacts : List[float]
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

        self.assertEqual(1.e-9, double(BDFCard(['1-9']), 0, 'field'))
        self.assertEqual(1.e+9, double(BDFCard(['1+9']), 0, 'field'))

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

    def check_double(self, method):
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

        #if check_space:
        with self.assertRaises(SyntaxError):
            method(BDFCard(['-9. 31-4']), 0, 'field')


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

        # blank
        #self.check_blank(integer_double_string_or_blank)

    def test_components(self):
        """tests the parse_components function"""
        with self.assertRaises(SyntaxError):
            integer_string_or_blank(BDFCard(['1b']), 0, 'field')

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
            val = parse_components(BDFCard(['12 3']), 0, 'field')

        # all numbers
        val = parse_components(BDFCard(['123456']), 0, 'field')
        self.assertEqual(val, '123456')

        # invalid 0's defined with numbers
        with self.assertRaises(SyntaxError):
            val = parse_components(BDFCard(['0123456']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = parse_components(BDFCard(['01']), 0, 'field')

        # doubles
        with self.assertRaises(SyntaxError):
            val = parse_components(BDFCard(['4524']), 0, 'field')

        # only 0 to 6
        with self.assertRaises(SyntaxError):
            val = parse_components(BDFCard(['7']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = parse_components(BDFCard([7]), 0, 'field')

        # dumb input
        with self.assertRaises(SyntaxError):
            val = parse_components(BDFCard(['4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = parse_components(BDFCard(['-4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = parse_components(BDFCard(['asdf']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = parse_components(BDFCard(['-1']), 0, 'field')

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
            val = components_or_blank(BDFCard(['12 3']), 0, 'field')

        # all numbers
        val = components_or_blank(BDFCard(['123456']), 0, 'field')
        self.assertEqual(val, '123456')

        # invalid 0's defined with numbers
        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['0123456']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['01']), 0, 'field')

        # doubles
        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['4524']), 0, 'field')

        # only 0 to 6
        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['7']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard([7]), 0, 'field')

        # dumb input
        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['-4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['asdf']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['-1']), 0, 'field')

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
            val = components_or_blank(BDFCard(['12 3']), 0, 'field', 'default')

        # all numbers
        val = components_or_blank(BDFCard(['123456']), 0, 'field', 'default')
        self.assertEqual(val, '123456')

        # invalid 0's defined with numbers
        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['0123456']), 0, 'field', 'default')

        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['01']), 0, 'field', 'default')

        # doubles
        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['4524']), 0, 'field', 'default')

        # only 0 to 6
        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['7']), 0, 'field', 'default')

        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard([7]), 0, 'field', 'default')

        # dumb input
        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['4.0']), 0, 'field', 'default')

        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['-4.0']), 0, 'field', 'default')

        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['asdf']), 0, 'field', 'default')

        with self.assertRaises(SyntaxError):
            val = components_or_blank(BDFCard(['-1']), 0, 'field', 'default')

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

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
