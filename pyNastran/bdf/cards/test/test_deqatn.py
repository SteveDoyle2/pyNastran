import os
import unittest
from io import StringIO

import numpy as np

#import pyNastran
from pyNastran.bdf.bdf import BDF, DEQATN
from pyNastran.bdf.cards.test.utils import save_load_deck


#root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')

#comment_bad = 'this is a bad comment'
#comment_good = 'this is a good comment\n'
class TestDEQATN(unittest.TestCase):
    """
    The cards tested are:
     * DEQATN
    """
    def _test_deqatn_1(self):
        """splitting a DEQATN doesnt work with is_list=False if it's a list"""
        model = BDF(debug=None)
        card = ["DEQATN", 1000, "MAXDIFF(t1,t2)=abs(t2-t1)/t1"]
        #with self.assertRaises(AttributeError): # TODO: fix this...
        model.add_card(card, "DEQATN", is_list=False)
        model.cross_reference()

        bdf_file = StringIO()
        with self.assertRaises(AttributeError):
            # this is a result of the previous error
            model.write_bdf(bdf_file)
        bdf_file.getvalue()

    def _test_deqatn_1b(self):
        """splitting a DEQATN doesnt work with is_list=False if it's a list"""
        model = BDF(debug=None)
        card = ["DEQATN", 1000, "MAXDIFF(t1,t2)=abs(t2-t1)/t1"]
        with self.assertRaises(ValueError):
            model.add_card(card, "DEQATN", is_list=False)

    def test_deqatn_1c(self):
        """
        works when is_list=True? and some magic flag is set..."

        def maxdiff(t1,t2):
            maxdiff = abs(t2-t1)/t1
            return maxdiff
        """
        model = BDF(debug=None)
        card = ["DEQATN      1000 MAXDIFF(t1,t2)=abs(t2-t1)/t1"]
        model.add_card(card, "DEQATN")
        model.cross_reference()
        #print(model.dequations[1000].func_str)

    def test_deqatn_2(self):
        """
        this works because is_list=True
        it doesn't work right now...
        """
        model = BDF(debug=None)
        card = ["DEQATN", 1000, "MAXDIFF(t1,t2)=abs(t2-t1)/t1"]
        with self.assertRaises(ValueError):  # this used to work....
            model.add_card(card, "DEQATN", is_list=True)
        #model.cross_reference()

        #bdf_file = StringIO()
        #with self.assertRaises(AttributeError): # TODO: fix this...
        #model.write_bdf(bdf_file)
        #bdf_file.getvalue()

    def test_deqatn_2b(self):
        """the corrected version of 2"""
        model = BDF(debug=None)
        card = ["DEQATN", 1000, "MAXDIFF(t1,t2)=abs(t2-t1)/t1"]

        with self.assertRaises(ValueError):
            model.add_card(card, "DEQATN", is_list=True)

    def test_deqatn_3(self):
        """
        Much simplier method of using add_card

        Creates the following equation:

        def maxdiff(t1,t2):
            maxdiff = abs(t2-t1)/t1
            return maxdiff
        """
        model = BDF(debug=None)
        card = ['DEQATN  1000',
                '        MAXDIFF(t1,t2)=abs(t2-t1)/t1']
        model.add_card(card, 'DEQATN', is_list=False)
        model.cross_reference()

        #print(model.dequations[1000].func_str)
        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.getvalue()
        #print(bdf_file.getvalue())
        bdf_file.close()

    def test_deqatn_4(self):
        """
        per nast/tpl/ptdmi1.dat

        def f(x, y, z):
            f =  1.
            l = 10.
            b = 4.
            h = 2.
            t1 = 200.
            t2 = 300.
            t = t1*(l-x)/l+t2*(x)/l
            return t
        """
        model = BDF(debug=None)
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
        model.cross_reference()
        #print(model.dequations[2].func_str)

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.getvalue()
        bdf_file.close()

    def test_deqatn_5(self):
        """
        per nast/tpl/ptdmi1.dat
        """
        model = BDF(debug=None)
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
        model.cross_reference()

        with open('junk.bdf', 'w') as bdf_file:
            model.write_bdf(bdf_file, close=False)
            #s.getvalue()
        os.remove('junk.bdf')

    def test_deqatn_6a(self):
        func_str = 'def f(x, y, z):\n'
        func_str += '    c = 3\n'
        func_str += '    return x + y + z + c\n'
        #func = exec(fnc_str)

        #s = StringIO()
        #s.write(s)
        #s.close()

        local_dict = {}
        exec(func_str, {}, local_dict)
        #print('locals() =', local_dict)
        f = local_dict['f']
        #out = f(1, 2, 3)
        #func = exec func_str
        assert f(1, 2, 3) == 9, func(1, 2, 3)

    def test_deqatn_6b(self):
        """this doesn't work because the commas are missing; see test_deqatn_6c"""
        model = BDF(debug=False)
        equation_id = 100
        eqs = [
            'f(x,y,z)=1.'  # the defaults
            'c = 3'
            'd = x + y + z + c'
        ]
        unused_deqatn = model.add_deqatn(equation_id, eqs, comment='')
        with self.assertRaises(SyntaxError):
            model.cross_reference()
        #print(deqatn.func_str)

    def test_deqatn_6c(self):
        """
        The add_deqatn version of test_deqatn_6a

        def f(x, y, z):
            f = 1.
            c = 3
            d = x + y + z + c
            return d
        """
        model = BDF(debug=False)
        equation_id = 100
        eqs = [
            'f(x,y,z)=1.',  # the defaults
            'c = 3',
            'd = x + y + z + c',
        ]
        unused_deqatn = model.add_deqatn(equation_id, eqs, comment='')
        model.cross_reference()
        #print(deqatn.func_str)

    def test_deqatn_6d(self):
        """
        The add_deqatn version of test_deqatn_6a with defaults

        def f(x, y, z, w=10.0):
            '''
            $deqatn
            DEQATN  100     f(x,y,z,w)=1.;
                    c = 3;
                    d = x + y + z + c
            '''
            f = 1.
            c = 3
            d = x + y + z + c
            return d
        """
        model = BDF(debug=False)
        equation_id = 100
        eqs = [
            'f(x,y,z,w)=1.',  # the defaults
            'c = 3',
            'd = x + y + z + c',
        ]
        unused_deqatn = model.add_deqatn(equation_id, eqs, comment='deqatn')
        default_values = {'w' : 10.}
        model.add_dtable(default_values, comment='dtable')
        model.cross_reference()
        #print(deqatn.func_str)

    def test_deqatn_6e(self):
        """
        tests comments
        """
        model = BDF(debug=False)
        card_lines = [
            'DEQATN  100     f(x,y,z,w)=1.;',
            'c = 3;',
            'd = x + y + z + c'
        ]
        with self.assertRaises(SyntaxError):
            model.add_card(card_lines, 'DEQATN', comment='deqatn', is_list=False)

        card_lines = [
            'DEQATN  100     f(x,y,z,w)=1.;',
            '        c = 3;',
            '        d = x + y + z + c'
        ]
        model.add_card(card_lines, 'DEQATN', comment='deqatn', is_list=False)
        deqatn = model.dequations[100]
        assert deqatn._comment == '$deqatn\n', '_comment=%r' % deqatn._comment

    def test_deqatn_7(self):
        """
        per nast/tpl/ptdmi1.dat

        def f(x, y, z):
            f = 1.
            l = 1+2+3++4/min(1,2)
            b = 4.
            h = 2.
            t1 = 200.
            t2 = 300.
            t = t1*(l-x)/l+t2*x/l
            return t
        """
        model = BDF(debug=None)
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
        model.cross_reference()
        #print(model.dequations[2].func_str)

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.getvalue()
        bdf_file.close()

    @unittest.skipUnless(2 < 1, 'skipping')
    def test_deqatn_8_skip(self):
        """
        based off nast/tpl/ptdmi1.dat

        What's going on with the last line?
        """
        model = BDF(debug=None)
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
        model.cross_reference()

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.getvalue()
        bdf_file.close()

    def test_deqatn_9(self):
        """
        based off nast/tpl/ptdmi1.dat
        """
        model = BDF(debug=None)
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
        model.cross_reference()

        s = StringIO()
        model.write_bdf(s, close=False)
        s.getvalue()
        s.close()

    def test_deqatn_9_tab(self):
        model = BDF(debug=None)
        card = [
            'deqatn  2       f(x,y,z)= 1.;',
            '\tL=1+2+3+',
            '\t+ 4/min(1,2);',
            '\tb= 4.;',
            '\th= 2.;',
            '\tt1= 200.;',
            '\tt2= 300.;',
            '\tt=t1*(L-x)/L+t2*x/L',
            '        +4'
        ]
        model.add_card(card, 'DEQATN', is_list=False)
        model.cross_reference()

        s = StringIO()
        model.write_bdf(s, close=False)
        s.getvalue()
        s.close()

    def test_deqatn_10(self):
        """
        based off nast/tpl/ptdmi1.dat

        def f(x, y, z):
            f = 1.
            l = x+y
            return l
        """
        model = BDF(debug=None)
        card = [
            'deqatn  2       f(x,y,z)= 1.;',
            '        L=x+y',
        ]
        model.add_card(card, 'DEQATN', is_list=False)
        model.cross_reference()

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.getvalue()
        bdf_file.close()

        eq = model.dequations[2]
        x = np.ones(10, dtype='float32')
        y = 2. * np.ones(10, dtype='float32')
        z = 3. * np.ones(10, dtype='float32')
        #print(model.dequations[2].func_str)
        #z = ones(10, dtype='float32')
        #out = eq.func(x, y, z)
        #out = eq.func(x, y, z)
        assert np.array_equal(eq.func(x, y, z), z)
        assert eq.func(1.0, 2.0, 1.0) == 3.0
        #print('out9 = %r' % out)


    def test_deqatn_11(self):
        """
        based off nast/tpl/ptdmi1.dat

        What's going on with this...
        """
        model = BDF(debug=None)
        deqatn_card = [
            'deqatn  2       f(x,y,z)= 1.;',
            '        L=x+y',
        ]
        model.add_card(deqatn_card, 'DEQATN', is_list=False)

        #dessub_desglb = 5
        dconstr_cards = [
            ['dconstr,5,10,',],
            ['dconstr,6,11,',],
        ]
        dresp_cards = [
            [ # card1
                'dresp2,10,respA,2',
                ',desvar,100,101,102',
            ],
            [ # card2
                'dresp2,11,respB,2',
                ',desvar,100,101,102',
            ],
            ## TODO: support this...
            #[
            #'dresp2,11,respB,F(A,B)=A+B**2*SIN(A*B)'
            #',desvar,100,101',
            #],
        ]
        #model.add_desvar(100, 'varA', 100.1, xlb=-1e20, xub=1e20,
                         #delx=None, ddval=None, comment='')
        desvar_cards = [
            ['desvar,100,varA,100.1',],
            ['desvar,101,varB,100.2',],
            ['desvar,102,varC,100.3',],
        ]

        for desvar in desvar_cards:
            model.add_card(desvar, 'DESVAR', is_list=False)
        for dconstr in dconstr_cards:
            model.add_card(dconstr, 'DCONSTR', is_list=False)
        for dresp in dresp_cards:
            model.add_card(dresp, 'DRESP2', is_list=False)
        model._verify_bdf()
        model.cross_reference()
        model._verify_bdf()
        save_load_deck(model)

    def test_deqatn_12(self):
        """
        made up

        def f1(a, b, c, d, r):
            f1 = a+b*c-(d**3+10.0)+sin(pi(1)*r)+a**2/(b-c)
            f = a+b-f1*d
            return f
        """
        model = BDF(debug=None)
        card_lines = [
            'DEQATN  41      F1(A,B,C,D,R) = A+B *C-(D**3 + 10.0) + sin(PI(1) * R)',
            '                + A**2 / (B - C); F = A + B - F1 * D',
        ]
        model.add_card(card_lines, 'DEQATN', is_list=False,
                       has_none=True)
        model.dequations[41].write_card()
        eqs = [
            'f(x,y)=1',
            'a = x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x',
            'b = a + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x'
        ]
        deqatn = model.add_deqatn(1000, eqs)
        deqatn.write_card()

        #print(deqatn)
        model.cross_reference()
        #print(model.dequations[41].eqs)
        #print(model.dequations[41].func_str)
        save_load_deck(model)

    def test_deqatn_13(self):
        """
        add_deqatn doesnt support semicolons (;) in the eqs
        You're defining them; break it up
        """
        eqs = [
            'f(x,y)=1',
            'a = x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x',
            'b = a + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x + x; c=42'
        ]
        model = BDF(debug=None)
        unused_deqatn = model.add_deqatn(1001, eqs)
        model.cross_reference()

    def test_deqatn_14(self):
        """
        def a(b, h):
            a = b * h
            return a
        """
        model = BDF(debug=None)
        deqatn_card = [
            'DEQATN  121     A(B,H) = B*H            $  ...equations for'
        ]
        model.add_card(deqatn_card, 'DEQATN', is_list=False)

    def test_deqatn_bad_1(self):
        """checks that a function name is not an argument"""
        model = BDF(debug=None)
        card = ["DEQATN      1000 A(a,b)=1"]
        model.add_card(card, "DEQATN")
        deqatn = model.dequations[1000]
        with self.assertRaises(RuntimeError):
            deqatn.cross_reference(model)
        #print(model.dequations[1000].func_str)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
