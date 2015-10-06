import unittest

import os
import StringIO
#import cStringIO
import pyNastran
from pyNastran.bdf.bdf import BDF

from numpy import zeros, array_equal

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')

comment_bad = 'this is a bad comment'
comment_good = 'this is a good comment\n'
class TestDEQATN(unittest.TestCase):
    """
    The cards tested are:
     * DEQATN
    """
    def test_deqatn_1(self):
        model = BDF(debug=None)
        #model.cards_to_read.add('DEQATN')
        #model.test_deqatn = True
        card = ["DEQATN", 1000, "MAXDIFF(t1,t2)=abs(t2-t1)/t1"]
        #with self.assertRaises(AttributeError): # TODO: fix this...
        model.add_card(card, "DEQATN", is_list=False)

        s = StringIO.StringIO()
        with self.assertRaises(AttributeError):
            # this is a result of the previous error
            model.write_bdf(s)
        s.getvalue()

    def test_deqatn_1b(self):
        model = BDF(debug=None)
        #model.cards_to_read.add('DEQATN')
        model.test_deqatn = True
        card = ["DEQATN", 1000, "MAXDIFF(t1,t2)=abs(t2-t1)/t1"]
        with self.assertRaises(ValueError):
            model.add_card(card, "DEQATN", is_list=False)

    def test_deqatn_2(self):
        model = BDF(debug=None)
        #model.cards_to_read.add('DEQATN')
        #model.test_deqatn = True
        card = ["DEQATN", 1000, "MAXDIFF(t1,t2)=abs(t2-t1)/t1"]
        model.add_card(card, "DEQATN", is_list=True)

        s = StringIO.StringIO()
        with self.assertRaises(AttributeError): # TODO: fix this...
            model.write_bdf(s)
        s.getvalue()

    def test_deqatn_2b(self):
        """the corrected version of 2"""
        model = BDF(debug=None)
        #model.cards_to_read.add('DEQATN')
        model.test_deqatn = True
        card = ["DEQATN", 1000, "MAXDIFF(t1,t2)=abs(t2-t1)/t1"]

        with self.assertRaises(ValueError):
            model.add_card(card, "DEQATN", is_list=True)

    def test_deqatn_3(self):
        model = BDF(debug=None)
        model.cards_to_read.add('DEQATN')
        model.test_deqatn = True
        card = ['DEQATN  1000',
                '        MAXDIFF(t1,t2)=abs(t2-t1)/t1']
        model.add_card(card, 'DEQATN', is_list=False)

        s = StringIO.StringIO()
        model.write_bdf(s, close=False)
        s.getvalue()
        #print(s.getvalue())
        s.close()

    def test_deqatn_4(self):
        """
        per nast/tpl/ptdmi1.dat
        """
        model = BDF(debug=None)
        #model.cards_to_read.add('DEQATN')
        model.test_deqatn = True
        card = [
            'deqatn  2       f(x,y,z)= 1.;',
            '        l=10.;',
            '        b= 4.;',
            '        h= 2.;',
            '        t1= 200.;',
            '        t2= 300.;',
            '        t=t1*(l-x)/l+t2*(x)/l',
         ]
        model.add_card(card, 'DEQATN', is_list=False)

        s = StringIO.StringIO()
        model.write_bdf(s, close=False)
        s.getvalue()
        s.close()

    def test_deqatn_5(self):
        """
        per nast/tpl/ptdmi1.dat
        """
        model = BDF(debug=None)
        model.cards_to_read.add('DEQATN')
        model.test_deqatn = True
        card = [
            'deqatn  2       f(x,y,z)= 1.;',
            '        l=10.;',
            '        b= 4.;',
            '        h= 2.;',
            '        t1= 200.;',
            '        t2= 300.;',
            '        t=t1*(l-x)/l+t2*(x)/l',
         ]
        model.add_card(card, 'DEQATN', is_list=False)

        f = open('junk.bdf', 'wb')
        model.write_bdf(f, close=False)
        #s.getvalue()
        f.close()
        os.remove('junk.bdf')

    def test_deqatn_6(self):
        func_str  = 'def f(x, y, z):\n'
        func_str += '    c = 3\n'
        func_str += '    return x + y + z + c\n'
        #func = exec(fnc_str)
        s = StringIO.StringIO()
        s.write(s)
        s.close()
        exec func_str
        f(1, 2, 3)
        #func = exec func_str
        assert f(1, 2, 3) == 9, func(1, 2, 3)

    def test_deqatn_7(self):
        """
        per nast/tpl/ptdmi1.dat
        """
        model = BDF(debug=None)
        model.cards_to_read.add('DEQATN')
        model.test_deqatn = True
        card = [
            'deqatn  2       f(x,y,z)= 1.;',
            '        L=1+2+3+',
            '        + 4/min(1,2);',
            '        b= 4.;',
            '        h= 2.;',
            '        t1= 200.;',
            '        t2= 300.;',
            '        t=t1*(L-x)/L+t2*x/L',
         ]
        model.add_card(card, 'DEQATN', is_list=False)

        s = StringIO.StringIO()
        model.write_bdf(s, close=False)
        s.getvalue()
        s.close()

    def test_deqatn_8(self):
        """
        per nast/tpl/ptdmi1.dat
        """
        model = BDF(debug=None)
        model.cards_to_read.add('DEQATN')
        model.test_deqatn = True
        card = [
            'deqatn  2       f(x,y,z)= 1.;',
            '        L=1+2+3+',
            '        + 4/min(1,2);',
            '        b= 4.;',
            '        h= 2.;',
            '        t1= 200.;',
            '        t2= 300.;',
            '        t=t1*(L-x)/L+t2*x/L;',
            '        +4'
        ]
        model.add_card(card, 'DEQATN', is_list=False)

        s = StringIO.StringIO()
        model.write_bdf(s, close=False)
        s.getvalue()
        s.close()

    def test_deqatn_8(self):
        """
        per nast/tpl/ptdmi1.dat
        """
        model = BDF(debug=None)
        model.cards_to_read.add('DEQATN')
        model.test_deqatn = True
        card = [
            'deqatn  2       f(x,y,z)= 1.;',
            '        L=1+2+3+',
            '        + 4/min(1,2);',
            '        b= 4.;',
            '        h= 2.;',
            '        t1= 200.;',
            '        t2= 300.;',
            '        t=t1*(L-x)/L+t2*x/L',
            '        +4'
        ]
        model.add_card(card, 'DEQATN', is_list=False)

        s = StringIO.StringIO()
        model.write_bdf(s, close=False)
        s.getvalue()
        s.close()

    def test_deqatn_8(self):
        """
        per nast/tpl/ptdmi1.dat
        """
        model = BDF(debug=None)
        model.cards_to_read.add('DEQATN')
        model.test_deqatn = True
        card = [
            'deqatn  2       f(x,y,z)= 1.;',
            '        L=x+y',
        ]
        model.add_card(card, 'DEQATN', is_list=False)

        s = StringIO.StringIO()
        model.write_bdf(s, close=False)
        s.getvalue()
        s.close()
        eq = model.dequations[2]
        x = zeros(10., dtype='float32')
        y = zeros(11., dtype='float32')
        z = zeros(12., dtype='float32')
        out = eq.func(x, y, z)
        out = eq.func(1.0, 2.0)
        print(out)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
