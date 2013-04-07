import unittest
from pyNastran.bdf.bdfInterface.BDF_Card import BDFCard
from pyNastran.bdf.format import (integer, integer_or_blank,
                                  double, double_or_blank, integer_double_or_blank,
                                  blank)

class Test(unittest.TestCase):

    def test_double_or_blank(self):
        """
        value = double_or_blank(card, n, fieldname, default=None)
        """
        card    = [1.0, '1.0', '1.', 'C',        None, None,          '', None, 'cat']
        exact   = [1.0,  1.0,   1.0, SyntaxError,None, 2.0,  SyntaxError, 1.0, SyntaxError]
        default = [None, None, None, None,       None, 2.0,         None, 1.0, 1.0]
        card = BDFCard(card)
        
        assert len(card) == len(exact), 'len(card)=%s len(exact)=%s' % (len(card), len(exact))
        assert len(card) == len(default), 'len(card)=%s len(default)=%s' % (len(card), len(default))
        i = 0
        for i, exacti in enumerate(exact):
            defaulti = default[i]
            if exacti == SyntaxError:
                with self.assertRaises(exacti):
                    double_or_blank(card, i, i, defaulti)
            else:
                value = double_or_blank(card, i, i, defaulti)
                self.assertEqual(value, exacti)
            i += 1

    def test_double(self):
        """
        value = double(card, n, fieldname)
        """
        card =  [1.0, '1.0', '1.', 'C',        None,                '']
        exact = [1.0,  1.0,   1.0, SyntaxError,SyntaxError,SyntaxError]
        card = BDFCard(card)
        
        i = 0
        for i, exacti in enumerate(exact):
            if exacti == SyntaxError:
                try:
                    with self.assertRaises(exacti):
                        double(card, i, i)
                except:
                    msg  = 'cardi=%s\n' % (card.field(i))
                    msg += 'exacti=%s' % (exacti)
                    print msg
                    raise
            else:
                value = double_or_blank(card, i, i)
                self.assertEqual(value, exacti)
            i += 1

if __name__ == '__main__':
    unittest.main()