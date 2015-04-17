import unittest
from pyNastran.bdf.bdf import BDF, BDFCard, SUPORT, SUPORT1, MPC

bdf = BDF(debug=False)
class TestConstraints(unittest.TestCase):
    def test_support_01(self):
        lines = ['SUPORT, 1,      126']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = SUPORT(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_suport_02(self):
        card = ['SUPORT']
        card = BDFCard(card)
        with self.assertRaises(AssertionError):  # too short
            SUPORT(card)

        card = ['SUPORT', '2']
        card = BDFCard(card)
        con = SUPORT(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(len(con.IDs), 1)
        self.assertEqual(len(con.Cs), 1)

        card = ['SUPORT', '2', '432']
        card = BDFCard(card)
        con = SUPORT(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '234')
        self.assertEqual(len(con.IDs), 1)
        self.assertEqual(len(con.Cs), 1)

        card = ['SUPORT', '2', '432', 3]
        card = BDFCard(card)
        con = SUPORT(card) # default
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(len(con.IDs), 2)
        self.assertEqual(len(con.Cs), 2)

        card = ['SUPORT', '2', None, 3]
        card = BDFCard(card)
        con = SUPORT(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(con.Cs[1], '0')
        self.assertEqual(len(con.IDs), 2)
        self.assertEqual(len(con.Cs), 2)

        card = ['SUPORT', '2', '432', 3, '1325']
        card = BDFCard(card)
        con = SUPORT(card) # default
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[1], '1235')
        self.assertEqual(len(con.IDs), 2)
        self.assertEqual(len(con.Cs), 2)

    def test_support1_01(self):
        lines = ['SUPORT1, 1,      126']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = SUPORT1(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_suport1_02(self):
        card = ['SUPORT1', '1']
        card = BDFCard(card)
        with self.assertRaises(AssertionError):  # too short
            SUPORT1(card)

        card = ['SUPORT1', '1', '2']
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '0')

        card = ['SUPORT1', '1', '2', '432']
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '234')

        card = ['SUPORT1', '1', '2', '432', 3]
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[1], 3)

        card = ['SUPORT1', '1', '2', None, 3]
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(con.Cs[1], '0')

        card = ['SUPORT1', '1', '2', '432', 3, '1325']
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[1], '1235')

    def test_mpc_01(self):
        card = ['MPC', 1, 1002, 1, 1., 1000, 1, -3.861003861]
        card = BDFCard(card)
        mpc = MPC(card)
        #print ' %r' % str(mpc)
        #print '%r' % mpc.write_card(size=8)
        #msg = mpc.write_card(size=8, double=False)
        self.assertEqual('MPC            1    1002       1      1.    1000       1  -3.861\n', mpc.write_card(size=8))

        model = BDF()

        card = ['MPC            1    1002       4      1.    1000       4-.129394',
                '                    1000       5-7.152-3    1000       6-.013655']
        msgA = ('MPC            1    1002       4      1.    1000       4-.129394\n'
                '                    1000       5-7.152-3    1000       6-.013655\n')

        card = model.add_card(card, 'MPC', is_list=False)
        mpc = MPC(card)
        #print('u%r' % msgA)
        #print('%r' % mpc.write_card(size=8))
        self.assertEqual(msgA, mpc.write_card(size=8))
        #print(mpc, type(mpc))

        card = []



if __name__ == '__main__':  # pragma: no cover
    unittest.main()
