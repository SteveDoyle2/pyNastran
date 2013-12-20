import unittest
from pyNastran.bdf.bdfInterface.BDF_Card import BDFCard
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank,
    integer_or_double, integer_double_or_blank,
    string, string_or_blank,
    integer_or_string, integer_string_or_blank,
    double_or_string, double_string_or_blank,
    blank, components, getType, interpret_value)

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
            blank(BDFCard([1]  ), 0, 'field')
        with self.assertRaises(SyntaxError):
            blank(BDFCard(['1']), 0, 'field')

        # float
        with self.assertRaises(SyntaxError):
            blank(BDFCard([1.]  ), 0, 'field')
        with self.assertRaises(SyntaxError):
            blank(BDFCard(['1.']), 0, 'field')

        # string
        with self.assertRaises(SyntaxError):
            blank(BDFCard(['a'] ), 0, 'field')
        with self.assertRaises(SyntaxError):
            blank(BDFCard(['1b']), 0, 'field')

        # blank
        blank(BDFCard(['']   ), 0, 'field')
        blank(BDFCard([None] ), 0, 'field')

    def test_double(self):
        """
        value = double(card, n, fieldname)
        """
        # integer
        with self.assertRaises(SyntaxError):
            double(BDFCard([1]  ), 0, 'field')
        with self.assertRaises(SyntaxError):
            double(BDFCard(['1']), 0, 'field')

        # float
        double(BDFCard([1.]  ), 0, 'field')
        double(BDFCard(['1.']), 0, 'field')

        # string
        with self.assertRaises(SyntaxError):
            double(BDFCard(['a'] ), 0, 'field')
        with self.assertRaises(SyntaxError):
            double(BDFCard(['1b']), 0, 'field')

        # blank
        with self.assertRaises(SyntaxError):
            double(BDFCard([''] ), 0, 'field')
        with self.assertRaises(SyntaxError):
            double(BDFCard([None] ), 0, 'field')

        card =  [1.0, '2.0', '3.', 'C',        None,                '']
        exact = [1.0,  2.0,   3.0, SyntaxError,SyntaxError,SyntaxError]
        self.run_function(double, card, exact)

    def test_integer(self):
        """
        value = integer(card, n, fieldname)
        """
        # integer
        integer(BDFCard([1]  ), 0, 'field')
        integer(BDFCard(['1']), 0, 'field')

        # float
        with self.assertRaises(SyntaxError):
            integer(BDFCard([1.]  ), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard(['1.']), 0, 'field')

        # string
        with self.assertRaises(SyntaxError):
            integer(BDFCard(['a'] ), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard(['1b']), 0, 'field')

        # blank
        with self.assertRaises(SyntaxError):
            integer(BDFCard([''] ), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer(BDFCard([None] ), 0, 'field')

        card =  [1, '2', '3.',         'C',        None,                '']
        exact = [1,  2,   SyntaxError, SyntaxError,SyntaxError,SyntaxError]
        self.run_function(integer, card, exact)

    def test_string(self):
        """
        value = string(card, n, fieldname)
        """
        # integer
        with self.assertRaises(SyntaxError):
            string(BDFCard([1]  ), 0, 'field')
        with self.assertRaises(SyntaxError):
            string(BDFCard(['1']), 0, 'field')

        # float
        with self.assertRaises(SyntaxError):
            string(BDFCard([1.]  ), 0, 'field')
        with self.assertRaises(SyntaxError):
            string(BDFCard(['1.']), 0, 'field')

        # string
        string(BDFCard(['a'] ), 0, 'field')
        string(BDFCard(['1b']), 0, 'field')

        # blank
        with self.assertRaises(SyntaxError):
            string(BDFCard([''] ), 0, 'field')
        with self.assertRaises(SyntaxError):
            string(BDFCard([None] ), 0, 'field')

        card =  ['a', 'b1', '3.',        'C',       None,         '', 'frog']
        exact = ['a', 'b1', SyntaxError, 'C',SyntaxError,SyntaxError, 'frog']
        self.run_function(string, card, exact)

    #-------------------------------------------------------
    def test_integer_or_blank(self):
        """
        value = integer_or_blank(card, n, fieldname)
        """
        # integer
        self.assertEqual(1, integer_or_blank(BDFCard([1]   ), 0, 'field'))
        self.assertEqual(2, integer_or_blank(BDFCard(['2'] ), 0, 'field'))

        # float
        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard([1.]  ), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard(['1.']), 0, 'field')

        # string
        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard(['a'] ), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_or_blank(BDFCard(['1b']), 0, 'field')

        # blank
        self.assertEqual('a', integer_or_blank(BDFCard(['']   ), 0, 'field', 'a'))
        self.assertEqual('b', integer_or_blank(BDFCard([None] ), 0, 'field', 'b'))


        card    = [1,     '2',        '3.',         'C', None, '']
        exact   = [1,       2, SyntaxError, SyntaxError, None, None]
        default = [None, None,        None,        None, None, None]
        self.run_function_default(integer_or_blank, card, exact, default)

    def test_double_or_blank(self):
        """
        value = double_or_blank(card, n, fieldname, default=None)
        """
        # integer
        card = BDFCard([1])
        with self.assertRaises(SyntaxError):
            val = double_or_blank(card, 0, 'field')
        card = BDFCard(['2'] )
        with self.assertRaises(SyntaxError):
            val = double_or_blank(card, 0, 'field')

        # float
        val = double_or_blank(BDFCard([1.]  ), 0, 'field')
        self.assertEquals(1., val)
        val = double_or_blank(BDFCard(['1.']), 0, 'field')
        self.assertEquals(1., val)
        val = double_or_blank(BDFCard(['1-3']), 0, 'field')
        self.assertEquals(1.e-3, val)

        # string
        with self.assertRaises(SyntaxError):
            double_or_blank(BDFCard(['a'] ), 0, 'field')
        with self.assertRaises(SyntaxError):
            double_or_blank(BDFCard(['1b']), 0, 'field')

        # blank
        double_or_blank(BDFCard(['   ']), 0, 'field')
        double_or_blank(BDFCard([None]), 0, 'field')

    def test_integer_double_or_blank(self):
        """
        value = double_or_blank(card, n, fieldname, default=None)
        """
        # integer
        integer_double_or_blank(BDFCard([1]    ), 0, 'field')
        integer_double_or_blank(BDFCard(['2']  ), 0, 'field')

        # float
        integer_double_or_blank(BDFCard([3.0]  ), 0, 'field')
        integer_double_or_blank(BDFCard(['4.0']), 0, 'field')
        integer_double_or_blank(BDFCard(['5.'] ), 0, 'field')

        # error - string
        with self.assertRaises(SyntaxError):
            integer_double_or_blank(BDFCard(['C']  ), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_double_or_blank(BDFCard(['1C']  ), 0, 'field')
        with self.assertRaises(SyntaxError):
            integer_double_or_blank(BDFCard(['C1']  ), 0, 'field')

        #card    = [1,    2.0, '3.0', '4.', 'C',        None, None,          '', None, 'cat']
        #exact   = [1,    2.0,  3.0,   4.0, SyntaxError,None, 2.0,  SyntaxError, 1.0, SyntaxError]
        #default = [None, None, None, None, None,       None, 2.0,         None, 1.0, 1.0]
        #self.run_function_default(integer_double_or_blank, card, exact, default)

    def test_components(self):
        # single ints
        val = components(BDFCard([0]), 0, 'field')
        self.assertEquals(val, '0')

        val = components(BDFCard([1]), 0, 'field')
        self.assertEquals(val, '1')

        # single strings
        val = components(BDFCard(['0']), 0, 'field')
        self.assertEquals(val, '0')

        val = components(BDFCard(['1']), 0, 'field')
        self.assertEquals(val, '1')

        # double ints
        val = components(BDFCard(['123']), 0, 'field')
        self.assertEquals(val, '123')

        val = components(BDFCard([123]), 0, 'field')
        self.assertEquals(val, '123')

        val = components(BDFCard([321]), 0, 'field')
        self.assertEquals(val, '123')

        # embedded whiteshape
        with self.assertRaises(SyntaxError):
            val = components(BDFCard(['12 3']), 0, 'field')

        # all numbers
        val = components(BDFCard(['123456']), 0, 'field')
        self.assertEquals(val, '123456')

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

    def test_bad(self):
        val = getType('1.000000000D+00')
        print "val = ", val
        val = interpret_value('1.000000000D+00')
        print "val = ", val

if __name__ == '__main__':
    unittest.main()