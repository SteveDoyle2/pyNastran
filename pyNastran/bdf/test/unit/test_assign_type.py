import unittest
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, integer_or_double, integer_double_or_blank,
    string, string_or_blank, double_or_string, double_string_or_blank,
    integer_or_string, integer_string_or_blank, integer_double_or_string,
    blank, components, components_or_blank, _get_dtype, interpret_value)

class ExtendedTestCase(unittest.TestCase):

    def assertRaisesWithMessage(self, msg, func, *args, **kwargs):
        try:
            func(*args, **kwargs)
            self.assertFail()
        except Exception as inst:
            self.assertEqual(inst.message, msg)

class Test(ExtendedTestCase):

    def run_function_default(self, f, card, exact, default):
        fieldname = 'f'
        assert len(card) == len(exact), 'len(card)=%s len(exact)=%s' % (len(card), len(exact))
        assert len(card) == len(default), 'len(card)=%s len(default)=%s' % (len(card), len(default))
        i = 0
        card = BDFCard(card)
        for i, exacti in enumerate(exact):
            defaulti = default[i]
            if exacti == SyntaxError:
                with self.assertRaises(exacti):
                    msg = 'field=%r exact=%r default=%r' % (card.field(i), exacti, defaulti)
                    #print msg
                    f(card, i, fieldname, defaulti)
            else:
                value = f(card, i, fieldname, defaulti)
                self.assertEqual(value, exacti)
            i += 1

    def run_function(self, f, card, exact):
        assert len(card) == len(exact), 'len(card)=%s len(exact)=%s' % (len(card), len(exact))
        i = 0
        fieldname = 'f'
        card = BDFCard(card)
        for i, exacti in enumerate(exact):
            if exacti == SyntaxError:
                with self.assertRaises(SyntaxError):
                    msg = 'field=%r exact=%r' % (card.field(i), exacti)
                    #print msg
                    f(card, i, fieldname)
            else:

                value = f(card, i, fieldname)
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
            integer(BDFCard([1.]), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard(['1.']), 0, 'field')

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
        self.assertEqual('1B', string(BDFCard(['1b']), 0, 'field'))

        # blank
        with self.assertRaises(SyntaxError):
            string(BDFCard(['']), 0, 'field')
        with self.assertRaises(SyntaxError):
            string(BDFCard([None]), 0, 'field')

        card = ['a', 'b1', '3.', 'C', None, '', 'frog']
        exact = ['A', 'B1', SyntaxError, 'C', SyntaxError, SyntaxError, 'FROG']
        self.run_function(string, card, exact)

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

        # string
        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard(['1b']), 0, 'field')

        # blank
        self.assertEqual('a', integer_or_blank(BDFCard(['']), 0, 'field', 'a'))
        self.assertEqual('b', integer_or_blank(BDFCard([None]), 0, 'field', 'b'))


        card = [1, '2', '3.', 'C', None, '']
        exact = [1, 2, SyntaxError, SyntaxError, None, None]
        default = [None, None, None, None, None, None]
        self.run_function_default(integer_or_blank, card, exact, default)


    def check_integer(self, method):
        # integer
        self.assertEqual(1, method(BDFCard([1]), 0, 'field'))
        self.assertEqual(2, method(BDFCard(['2']), 0, 'field'))

    def check_double(self, method):
        # float
        method(BDFCard([3.0]), 0, 'field')
        method(BDFCard(['4.0']), 0, 'field')
        method(BDFCard(['5.']), 0, 'field')

        self.assertEqual(1.0, double(BDFCard([1.]), 0, 'field'))
        self.assertEqual(1.0, double(BDFCard(['1.']), 0, 'field'))
        self.assertEqual(-9.31e-4, double(BDFCard(['-9.31-4']), 0, 'field'))

        # float
        val = method(BDFCard([1.]), 0, 'field')
        self.assertEqual(1., val)
        val = method(BDFCard(['1.']), 0, 'field')
        self.assertEqual(1., val)
        val = method(BDFCard(['1-3']), 0, 'field')
        self.assertEqual(1.e-3, val)

    def test_double_or_blank(self):
        """
        value = double_or_blank(card, n, fieldname, default=None)
        """
        # integer
        card = BDFCard([1])
        with self.assertRaises(SyntaxError):
            val = double_or_blank(card, 0, 'field')
        card = BDFCard(['2'])
        with self.assertRaises(SyntaxError):
            val = double_or_blank(card, 0, 'field')

        self.check_double(double_or_blank)

        # string
        with self.assertRaises(SyntaxError):
            double_or_blank(BDFCard(['a']), 0, 'field')
        with self.assertRaises(SyntaxError):
            double_or_blank(BDFCard(['1b']), 0, 'field')

        # blank
        double_or_blank(BDFCard(['   ']), 0, 'field')
        double_or_blank(BDFCard([None]), 0, 'field')

    def test_double_or_string(self):
        # out of range
        with self.assertRaises(SyntaxError):
            double_or_string(BDFCard([1.]), 1, 'field')

    def test_double_string_or_blank(self):
        # out of range
        #with self.assertRaises(SyntaxError):
        self.assertEqual(1., double_string_or_blank(BDFCard([1.]), 0, 'field'))
        self.assertEqual(1., double_string_or_blank(BDFCard(['1.']), 0, 'field'))
        self.assertEqual('CAT', double_string_or_blank(BDFCard(['CAT']), 0, 'field'))
        self.assertEqual(None, double_string_or_blank(BDFCard([None]), 0, 'field'))

    def test_integer_or_double(self):
        # out of range
        with self.assertRaises(SyntaxError):
            integer_or_double(BDFCard([1.]), 1, 'field')

    def test_integer_or_string(self):
        # out of range
        with self.assertRaises(SyntaxError):
            integer_or_string(BDFCard([1.]), 1, 'field')

    def test_integer_double_or_blank(self):
        """
        value = double_or_blank(card, n, fieldname, default=None)
        """
        # integer
        self.check_integer(integer_double_or_blank)

        # float
        self.check_double(integer_double_or_blank)

        # error - string
        with self.assertRaises(SyntaxError):
            integer_double_or_blank(BDFCard(['C']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_double_or_blank(BDFCard(['1C']), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_double_or_blank(BDFCard(['C1']), 0, 'field')

        # blank
        double_or_blank(BDFCard(['   ']), 0, 'field')
        double_or_blank(BDFCard([None]), 0, 'field')

        #card    = [1,    2.0, '3.0', '4.', 'C',        None, None,          '', None, 'cat']
        #exact   = [1,    2.0,  3.0,   4.0, SyntaxError,None, 2.0,  SyntaxError, 1.0, SyntaxError]
        #default = [None, None, None, None, None,       None, 2.0,         None, 1.0, 1.0]
        #self.run_function_default(integer_double_or_blank, card, exact, default)

    def test_integer_string_or_blank(self):
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

        # string
        self.assertEqual('LOAD', integer_string_or_blank(BDFCard(['load']), 0, 'field'))

        # blank
        integer_string_or_blank(BDFCard(['   ']), 0, 'field')
        integer_string_or_blank(BDFCard([None]), 0, 'field')


    #def test_double_string_or_blank(self):

    def test_integer_double_or_string(self):
        # out of range
        with self.assertRaises(SyntaxError):
            integer_or_double(BDFCard([1.]), 1, 'field')

        # integer
        self.check_integer(integer_double_or_string)

        # float
        self.check_double(integer_double_or_string)

        # string
        self.assertEqual('LOAD', integer_double_or_string(BDFCard(['load']), 0, 'field'))

    def test_components(self):
        # single ints
        val = components(BDFCard([0]), 0, 'field')
        self.assertEqual(val, '0')

        val = components(BDFCard([1]), 0, 'field')
        self.assertEqual(val, '1')

        # single strings
        val = components(BDFCard(['0']), 0, 'field')
        self.assertEqual(val, '0')

        val = components(BDFCard(['1']), 0, 'field')
        self.assertEqual(val, '1')

        # double ints
        val = components(BDFCard(['123']), 0, 'field')
        self.assertEqual(val, '123')

        val = components(BDFCard([123]), 0, 'field')
        self.assertEqual(val, '123')

        val = components(BDFCard([321]), 0, 'field')
        self.assertEqual(val, '123')

        # embedded whiteshape
        with self.assertRaises(SyntaxError):
            val = components(BDFCard(['12 3']), 0, 'field')

        # all numbers
        val = components(BDFCard(['123456']), 0, 'field')
        self.assertEqual(val, '123456')

        # invalid 0's defined with numbers
        with self.assertRaises(SyntaxError):
            val = components(BDFCard(['0123456']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = components(BDFCard(['01']), 0, 'field')

        # doubles
        with self.assertRaises(SyntaxError):
            val = components(BDFCard(['4524']), 0, 'field')

        # only 0 to 6
        with self.assertRaises(SyntaxError):
            val = components(BDFCard(['7']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = components(BDFCard([7]), 0, 'field')

        # dumb input
        with self.assertRaises(SyntaxError):
            val = components(BDFCard(['4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = components(BDFCard(['-4.0']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = components(BDFCard(['asdf']), 0, 'field')

        with self.assertRaises(SyntaxError):
            val = components(BDFCard(['-1']), 0, 'field')

        # blank
        #components(BDFCard(['   ']), 0, 'field')
        #components(BDFCard([None]), 0, 'field')

    def test_components_or_blank_01(self):
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
        val = _get_dtype('1.000000000D+00')
        #print "val = ", val
        val = interpret_value('1.000000000D+00')
        #print "val = ", val


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
