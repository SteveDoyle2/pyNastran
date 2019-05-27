from itertools import count
import unittest
from pyNastran.bdf.bdf import BDF, BDFCard, SUPORT, SUPORT1, MPC
from pyNastran.bdf.cards.test.utils import save_load_deck

class TestConstraints(unittest.TestCase):
    def test_support_01(self):
        lines = ['SUPORT, 1,      126']
        model = BDF(debug=False)
        con = model.add_card(lines, 'SUPORT', comment='',
                             ifile=None, is_list=False, has_none=True)
        model.add_grid(1, [0., 0., 0.])

        con = model.suport[0]
        size = 8
        con.write_card(size, 'dummy')
        con.raw_fields()
        save_load_deck(model)

    def test_suport_02(self):
        card = ['SUPORT']
        cardi = BDFCard(card)
        with self.assertRaises(AssertionError):  # too short
            SUPORT.add_card(cardi)

        card = ['SUPORT', '2']
        cardi = BDFCard(card)
        con = SUPORT.add_card(cardi) # default
        self.assertEqual(con.nodes[0], 2)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(len(con.nodes), 1)
        self.assertEqual(len(con.Cs), 1)

        card = ['SUPORT', '2', '432']
        cardi = BDFCard(card)
        con = SUPORT.add_card(cardi) # default
        self.assertEqual(con.nodes[0], 2)
        self.assertEqual(con.Cs[0], '234')
        self.assertEqual(len(con.nodes), 1)
        self.assertEqual(len(con.Cs), 1)

        card = ['SUPORT', '2', '432', 3]
        cardi = BDFCard(card)
        con = SUPORT.add_card(cardi) # default
        self.assertEqual(con.nodes[1], 3)
        self.assertEqual(len(con.nodes), 2)
        self.assertEqual(len(con.Cs), 2)

        card = ['SUPORT', '2', None, 3]
        cardi = BDFCard(card)
        con = SUPORT.add_card(cardi) # default
        self.assertEqual(con.nodes[0], 2)
        self.assertEqual(con.nodes[1], 3)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(con.Cs[1], '0')
        self.assertEqual(len(con.nodes), 2)
        self.assertEqual(len(con.Cs), 2)

        card = ['SUPORT', '2', '432', 3, '1325']
        cardi = BDFCard(card)
        con = SUPORT.add_card(cardi) # default
        self.assertEqual(con.nodes[1], 3)
        self.assertEqual(con.Cs[1], '1235')
        self.assertEqual(len(con.nodes), 2)
        self.assertEqual(len(con.Cs), 2)

    def test_support1_01(self):
        lines = ['SUPORT1, 1,      126']
        model = BDF(debug=False)
        card = model._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        con = SUPORT1.add_card(cardi)
        con.write_card(size, 'dummy')
        con.raw_fields()
        save_load_deck(model)

    def test_suport1_02(self):
        card = ['SUPORT1', '1']
        card_obj = BDFCard(card)
        with self.assertRaises(AssertionError):  # too short
            SUPORT1.add_card(card_obj)

        card = ['SUPORT1', '1', '2']
        card_obj = BDFCard(card)
        con = SUPORT1.add_card(card_obj) # default
        self.assertEqual(con.nodes[0], 2)
        self.assertEqual(con.Cs[0], '0')

        card = ['SUPORT1', '1', '2', '432']
        card_obj = BDFCard(card)
        con = SUPORT1.add_card(card_obj) # default
        self.assertEqual(con.nodes[0], 2)
        self.assertEqual(con.Cs[0], '234')

        card = ['SUPORT1', '1', '2', '432', 3]
        card_obj = BDFCard(card)
        con = SUPORT1.add_card(card_obj) # default
        self.assertEqual(con.nodes[1], 3)

        card = ['SUPORT1', '1', '2', None, 3]
        card_obj = BDFCard(card)
        con = SUPORT1.add_card(card_obj) # default
        self.assertEqual(con.nodes[0], 2)
        self.assertEqual(con.nodes[1], 3)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(con.Cs[1], '0')

        card = ['SUPORT1', '1', '2', '432', 3, '1325']
        card = BDFCard(card)
        con = SUPORT1.add_card(card) # default
        self.assertEqual(con.nodes[1], 3)
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
        check_card(msg8, msg_8_actual)
        check_card(msg16, msg_16_actual)

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
        unused_msg_16_double_actual = mpc.write_card(size=16, is_double=True)
        check_card(msg8, msg_8_actual)
        check_card(msg16, msg_16_actual)

    def test_mpcadd(self):
        """tests MPCADD"""
        model = BDF(debug=False)
        mpc_id = 42
        sets = [1, 2]
        model.add_mpcadd(mpc_id, sets, comment='mpcadd')

        components = ['42', '3', '1']
        coefficients = [1000, 1, 101]
        node_ids = [2, 3, 4]
        mpc = model.add_mpc(1, node_ids, components, coefficients, comment='mpc')
        with self.assertRaises(AssertionError):
            mpc.validate()
        model.mpcs = {}

        coefficients = [1000., 1., 101.]
        mpc = model.add_mpc(1, node_ids, components, coefficients, comment='mpc')
        mpc.raw_fields()
        card = mpc.write_card(size=8)
        mpc.write_card(size=16, is_double=False)
        mpc.write_card(size=16, is_double=True)

        assert mpc.coefficients == coefficients
        assert mpc.nodes == node_ids

        model.pop_parse_errors()
        model.pop_xref_errors()

        model.mpcs = {}
        model.add_card(card.split('\n')[1:], 'MPC', is_list=False)
        with self.assertRaises(KeyError):
            model.cross_reference()

        model.add_mpc(2, node_ids, components, coefficients, comment='mpc')
        with self.assertRaises(KeyError):
            model.cross_reference()
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        model.add_grid(4, [0., 0., 0.])
        model.cross_reference()
        check_mpc_spc(model)
        save_load_deck(model)

    def test_spcadd(self):
        """tests SPCADD"""
        model = BDF(debug=False)
        spc_id = 42
        sets = [1, 2]
        model.add_spcadd(spc_id, sets, comment='spcadd')

        components = ['42', '3', '1']
        enforced = [1000, 1, 101]
        node_ids = [2, 3, 4]
        spc = model.add_spc(1, node_ids, components, enforced, comment='spc')
        with self.assertRaises(AssertionError):
            spc.validate()
        model.spcs = {}

        enforced = [1000., 1., 101.]
        node_ids = [2, 3, 4]
        spc = model.add_spc(1, node_ids, components, enforced, comment='spc')
        spc.validate()
        #assert spc.enforced == enforced
        assert spc.nodes == node_ids

        spc.raw_fields()
        card = spc.write_card(size=8)
        #print(card)
        spc.write_card(size=16, is_double=False)
        spc.write_card(size=16, is_double=True)
        model.pop_parse_errors()
        model.pop_xref_errors()

        model.spcs = {}
        model.add_card(card.split('\n')[1:], 'SPC', is_list=False)
        with self.assertRaises(KeyError):
            model.cross_reference()

        model.add_spc(2, node_ids, components, enforced, comment='spc')
        with self.assertRaises(KeyError):
            model.cross_reference()
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        model.add_grid(4, [0., 0., 0.])
        check_mpc_spc(model)
        save_load_deck(model)

    def test_spcoff(self):
        """tests SPCOFF/SPCOFF1"""
        model = BDF(debug=False)
        with self.assertRaises(KeyError):
            model.EmptyNodes([1, 2], msg='')
        #model.add_spcoff()
        card_lines = ['SPCOFF', 1]
        model.add_card(card_lines, 'SPCOFF', comment='spcoff', is_list=True, has_none=True)
        spcoff = model.spcoffs['SPCOFF'][0]
        spcoff.write_card(size=8)
        spcoff.write_card(size=16)
        spcoff.raw_fields()
        with self.assertRaises(KeyError):
            spcoff.cross_reference(model)
        model.add_grid(1, [0., 0., 0.])
        spcoff.cross_reference(model)

        card_lines = ['SPCOFF1', 24, 43]
        model.add_card(card_lines, 'SPCOFF1', comment='spcoff1', is_list=True, has_none=True)

        card_lines = ['SPCOFF1', 5, 50, 'THRU', 52]
        model.add_card(card_lines, 'SPCOFF1', comment='spcoff1', is_list=True, has_none=True)

        model.pop_parse_errors()
        model.pop_xref_errors()
        spcoff1 = model.spcoffs['SPCOFF1'][0]
        spcoff1.write_card(size=8)
        spcoff1.write_card(size=16)
        spcoff1.raw_fields()
        model.uncross_reference()

        with self.assertRaises(KeyError):
            spcoff1.cross_reference(model)
        model.add_grid(43, [0., 0., 0.])
        spcoff1.cross_reference(model)
        model.uncross_reference()

        model.add_grid(50, [0., 0., 0.])
        model.add_grid(51, [0., 0., 0.])
        model.add_grid(52, [0., 0., 0.])
        model.cross_reference()
        model.uncross_reference()
        model.validate()
        model.safe_cross_reference()

        save_load_deck(model, run_remove_unused=False, run_save_load_hdf5=False)

    def test_gmspc(self):
        """tests GMSPC"""
        model = BDF(debug=False, log=None, mode='msc')
        conid = 1
        component = 42
        entity = 'cat'
        entity_id = 100
        model.add_gmspc(conid, component, entity, entity_id, comment='gmspc')
        gmspc = model.spcs[1][0]
        card = gmspc.write_card(size=8)
        gmspc.raw_fields()
        gmspc.write_card(size=16)
        gmspc.write_card(size=16, is_double=True)
        model.add_card(card.split('\n')[1:], 'GMSPC', is_list=False)

        model.pop_parse_errors()
        model.validate()
        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        str(gmspc)
        check_mpc_spc(model)
        save_load_deck(model)

    def test_spcax(self):
        """tests SPCAX"""
        model = BDF(debug=False, log=None, mode='msc')
        conid = 1
        ringax = 42
        hid = 43
        component = 52
        enforced = 101.
        model.add_spcax(conid, ringax, hid, component, enforced, comment='spcax')
        spcax = model.spcs[1][0]
        card = spcax.write_card(size=8)
        spcax.raw_fields()
        spcax.write_card(size=16)
        spcax.write_card(size=16, is_double=True)
        model.add_card(card.split('\n')[1:], 'SPCAX', is_list=False)

        model.pop_parse_errors()
        model.validate()
        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        str(spcax)
        check_mpc_spc(model)
        save_load_deck(model)

def check_card(msg_expected, msg_actual):
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

def check_mpc_spc(model):
    """simple MPC/SPC checks"""
    mpc_ids = list(model.mpcs.keys()) + list(model.mpcadds.keys())
    spc_ids = list(model.spcs.keys()) + list(model.spcadds.keys())
    #node_ids = []
    #for nid in rbe.independent_nodes + rbe.dependent_nodes:
        #node_ids.append(nid)

    for mpc_id in mpc_ids:
        unused_mpcs1 = model.get_reduced_mpcs(mpc_id, consider_mpcadd=True, stop_on_failure=True)
        unused_mpcs2 = model.get_mpcs(mpc_id, consider_mpcadd=True, stop_on_failure=True)
    for spc_id in spc_ids:
        unused_spcs1 = model.get_reduced_spcs(spc_id, consider_spcadd=True, stop_on_failure=True)
        unused_spcs2 = model.get_spcs(spc_id, consider_nodes=False, stop_on_failure=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
