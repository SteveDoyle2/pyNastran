from __future__ import print_function
from itertools import count
import unittest
from pyNastran.bdf.bdf import BDF, BDFCard, SUPORT, SUPORT1, MPC

bdf = BDF(debug=False)
class TestConstraints(unittest.TestCase):
    def test_support_01(self):
        lines = ['SUPORT, 1,      126']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        con = SUPORT.add_card(card)
        con.write_card(size, 'dummy')
        con.raw_fields()

    def test_suport_02(self):
        card = ['SUPORT']
        cardi = BDFCard(card)
        with self.assertRaises(AssertionError):  # too short
            SUPORT.add_card(cardi)

        card = ['SUPORT', '2']
        cardi = BDFCard(card)
        con = SUPORT.add_card(cardi) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(len(con.IDs), 1)
        self.assertEqual(len(con.Cs), 1)

        card = ['SUPORT', '2', '432']
        cardi = BDFCard(card)
        con = SUPORT.add_card(cardi) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '234')
        self.assertEqual(len(con.IDs), 1)
        self.assertEqual(len(con.Cs), 1)

        card = ['SUPORT', '2', '432', 3]
        cardi = BDFCard(card)
        con = SUPORT.add_card(cardi) # default
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(len(con.IDs), 2)
        self.assertEqual(len(con.Cs), 2)

        card = ['SUPORT', '2', None, 3]
        cardi = BDFCard(card)
        con = SUPORT.add_card(cardi) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(con.Cs[1], '0')
        self.assertEqual(len(con.IDs), 2)
        self.assertEqual(len(con.Cs), 2)

        card = ['SUPORT', '2', '432', 3, '1325']
        cardi = BDFCard(card)
        con = SUPORT.add_card(cardi) # default
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[1], '1235')
        self.assertEqual(len(con.IDs), 2)
        self.assertEqual(len(con.Cs), 2)

    def test_support1_01(self):
        lines = ['SUPORT1, 1,      126']
        card = bdf.process_card(lines)
        cardi = BDFCard(card)

        size = 8
        con = SUPORT1.add_card(cardi)
        con.write_card(size, 'dummy')
        con.raw_fields()

    def test_suport1_02(self):
        card = ['SUPORT1', '1']
        card = BDFCard(card)
        with self.assertRaises(AssertionError):  # too short
            SUPORT1.add_card(card)

        card = ['SUPORT1', '1', '2']
        cardi = BDFCard(card)
        con = SUPORT1.add_card(cardi) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '0')

        card = ['SUPORT1', '1', '2', '432']
        cardi = BDFCard(card)
        con = SUPORT1.add_card(cardi) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '234')

        card = ['SUPORT1', '1', '2', '432', 3]
        cardi = BDFCard(card)
        con = SUPORT1.add_card(cardi) # default
        self.assertEqual(con.IDs[1], 3)

        card = ['SUPORT1', '1', '2', None, 3]
        cardi = BDFCard(card)
        con = SUPORT1.add_card(cardi) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(con.Cs[1], '0')

        card = ['SUPORT1', '1', '2', '432', 3, '1325']
        card = BDFCard(card)
        con = SUPORT1.add_card(card) # default
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[1], '1235')

    def test_mpc_01(self):
        card = ['MPC', 1, 1002, 1, 1., 1000, 1, -3.861003861]
        card = BDFCard(card)
        mpc = MPC.add_card(card)
        #print ' %r' % str(mpc)
        #print '%r' % mpc.write_card(size=8)
        #msg = mpc.write_card(size=8, double=False)
        self.assertEqual(
            'MPC            1    1002       1      1.    1000       1  -3.861\n',
            mpc.write_card(size=8))

    def test_mpc_02(self):
        model = BDF(debug=False)

        card = ['MPC            1    1002       4      1.    1000       4-.129394',
                '                    1000       5-7.152-3    1000       6-.013655']
        msg8 = ('MPC            1    1002       4      1.    1000       4-.129394\n'
                '                    1000       5-7.152-3    1000       6-.013655\n')
        msg16 = (
            'MPC*                   1            1002               4              1.\n'
            '*                   1000               4        -.129394\n'
            '*                                   1000               5        -.007152\n'
            '*                   1000               6        -.013655\n')


        card = model.add_card(card, 'MPC', is_list=False)
        assert card is not None
        mpc = MPC.add_card(card)
        msg_8_actual = mpc.write_card(size=8)
        msg_16_actual = mpc.write_card(size=16)
        self.check_card(msg8, msg_8_actual)
        self.check_card(msg16, msg_16_actual)

    def test_mpc_03(self):
        model = BDF(debug=False)

        card = [
            'MPC            1    1002       4      1.    1000       4-.129394',
            '                    1000       5-7.152-3    1000       6-.013655',
            '                    1004       2   123.3',]
        msg8 = (
            'MPC            1    1002       4      1.    1000       4-.129394\n'
            '                    1000       5-7.152-3    1000       6-.013655\n'
            '                    1004       2   123.3\n')
        msg16 = (
            'MPC*                   1            1002               4              1.\n'
            '*                   1000               4        -.129394\n'
            '*                                   1000               5        -.007152\n'
            '*                   1000               6        -.013655\n'
            '*                                   1004               2           123.3\n'
            '*\n')


        card = model.add_card(card, 'MPC', is_list=False)
        assert card is not None
        mpc = MPC.add_card(card)
        msg_8_actual = mpc.write_card(size=8)
        msg_16_actual = mpc.write_card(size=16)
        msg_16_double_actual = mpc.write_card(size=16, is_double=True)
        self.check_card(msg8, msg_8_actual)
        self.check_card(msg16, msg_16_actual)

    def check_card(self, msg_expected, msg_actual):
        if isinstance(msg_expected, tuple):
            msg_expected = msg_expected[0]
        msg_expected_lines = msg_expected.split('\n')
        msg_actual_lines = msg_actual.split('\n')
        for i, actual, expected  in zip(count(), msg_actual_lines, msg_expected_lines):
            msg = 'Error on line %i\n' % i
            msg += 'actual=\n%s\n' % actual
            msg += 'expected=\n%s\n\n' % expected
            msg += 'Actual Card =\n%s\n' % '\n'.join(msg_actual_lines)
            msg += '\nExpected Card =\n%s\n' % '\n'.join(msg_expected_lines)
            assert actual == expected, msg



if __name__ == '__main__':  # pragma: no cover
    unittest.main()
