from itertools import count
import unittest
import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF, BDFCard # , SUPORT, SUPORT1, MPC
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck # , read_write_op2_geom

RUN_OP2 = False


class TestConstraints(unittest.TestCase):
    def test_support_01(self):
        model = BDF(debug=False)
        suport = model.suport

        lines = ['SUPORT, 1,      126']
        con = model.add_card(lines, 'SUPORT', comment='',
                             ifile=None, is_list=False, has_none=True)
        model.add_grid(1, [0., 0., 0.])

        size = 8
        suport.write(size, 'dummy')
        #con.raw_fields()
        save_load_deck(model)

    def test_suport_02(self):
        model = BDF(debug=False)
        suport = model.suport

        card = ['SUPORT']
        cardi = BDFCard(card)
        #with self.assertRaises(AssertionError):  # too short
            #SUPORT.add_card(cardi)

        card = ['SUPORT', '2']
        cardi = BDFCard(card)
        icon = suport.add_set_card(cardi) - 1 # default
        model.setup()
        self.assertEqual(suport.node_id[icon], 2)
        self.assertEqual(suport.component[icon], 0)
        self.assertEqual(len(suport.node_id), 1)
        self.assertEqual(len(suport.component), 1)

        card = ['SUPORT', '2', '432']
        cardi = BDFCard(card)
        icon = suport.add_set_card(cardi) - 1 # default
        model.setup()
        self.assertEqual(suport.node_id[icon], 2)
        self.assertEqual(suport.component[icon], 234)
        self.assertEqual(len(suport.node_id), 2)
        self.assertEqual(len(suport.component), 2)

        card = ['SUPORT', '3', '4325', 4]
        cardi = BDFCard(card)
        icon = suport.add_set_card(cardi) - 1 # default
        model.setup()
        self.assertEqual(suport.node_id[icon-1], 3)
        self.assertEqual(suport.node_id[icon], 4)
        self.assertEqual(suport.component[icon-1], 2345)
        self.assertEqual(suport.component[icon], 0)
        self.assertEqual(len(suport.node_id), 4)
        self.assertEqual(len(suport.component), 4)

        card = ['SUPORT', '2', None, 3]
        cardi = BDFCard(card)
        icon = suport.add_set_card(cardi) - 1 # default
        model.setup()
        self.assertEqual(suport.node_id[icon-1], 2)
        self.assertEqual(suport.node_id[icon], 3)
        self.assertEqual(suport.component[icon-1], 0)
        self.assertEqual(suport.component[icon], 0)
        self.assertEqual(len(suport.node_id), 6)
        self.assertEqual(len(suport.component), 6)

        card = ['SUPORT', '2', '432', 3, '1325']
        cardi = BDFCard(card)
        icon = suport.add_set_card(cardi) - 1 # default
        model.setup()
        self.assertEqual(suport.node_id[icon-1], 2)
        self.assertEqual(suport.node_id[icon], 3)
        self.assertEqual(suport.component[icon-1], 234)
        self.assertEqual(suport.component[icon], 1235)
        self.assertEqual(len(suport.node_id), 8)
        self.assertEqual(len(suport.component), 8)

    def test_suport1_01(self):
        model = BDF(debug=False)
        suport = model.suport

        lines = ['SUPORT1, 1,      126']
        card = model._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        con = suport.add_set1_card(cardi)
        suport.write(size, 'dummy')
        #con.raw_fields()
        save_load_deck(model)

    def test_suport1_02(self):
        model = BDF(debug=False)
        suport = model.suport

        card = ['SUPORT1', '1']
        card_obj = BDFCard(card)
        with self.assertRaises(AttributeError):  # removed
            suport.add_card(card_obj)
        #with self.assertRaises(AssertionError):  # too short
            #suport.add_set1_card(card_obj)

        card = ['SUPORT1', '1', '2']
        card_obj = BDFCard(card)
        con = suport.add_set1_card(card_obj) # default
        model.setup()
        #print(suport.write())
        suporti = suport.slice_card_by_id(1)
        self.assertEqual(suporti.node_id[0], 2)
        self.assertEqual(suporti.component[0], 0)

        card = ['SUPORT1', '1', '2', '432']
        card_obj = BDFCard(card)
        con = suport.add_set1_card(card_obj) # default
        model.setup()
        self.assertEqual(suport.node_id[1], 2)
        self.assertEqual(suport.component[1], 234)

        card = ['SUPORT1', '1', '2', '432', 3]
        card_obj = BDFCard(card)
        con = suport.add_set1_card(card_obj) # default
        model.setup()
        self.assertEqual(suport.node_id[3], 3)

        card = ['SUPORT1', '1', '2', None, 3]
        card_obj = BDFCard(card)
        con = suport.add_set1_card(card_obj) # default
        model.setup()
        self.assertEqual(suport.node_id[4], 2)
        self.assertEqual(suport.node_id[5], 3)
        self.assertEqual(suport.component[4], 0)
        self.assertEqual(suport.component[5], 0)

        card = ['SUPORT1', '1', '2', '432', 3, '1325']
        card = BDFCard(card)
        con = suport.add_set1_card(card) # default
        model.setup()
        self.assertEqual(suport.node_id[7], 3)
        self.assertEqual(suport.component[7], 1235)

    def test_mpc_01(self):
        card = ['MPC', 1, 1002, 1, 1., 1000, 1, -3.861003861]
        card = BDFCard(card)
        model = BDF(debug=False, mode='nx')
        mpc = model.mpc
        mpc_id = mpc.add_card(card)
        #print ' %r' % str(mpc)
        #print '%r' % mpc.write_card(size=8)
        #msg = mpc.write_card(size=8, double=False)
        model.setup(run_geom_check=True)
        msg_actual = mpc.write(size=8)
        self.assertEqual(
           #'MPC            1    1002       1      1.    1000       1  -3.861\n',
            'MPC            1    1002       1      1.    1000       1-3.86100\n',
            msg_actual)

    def test_mpc_02(self):
        model = BDF(debug=False)

        card = ['MPC            1    1002       4      1.    1000       4-.129394',
                '                    1000       5-7.152-3    1000       6-.013655']
        msg8 = ('MPC            1    1002       4      1.    1000       4-.129394\n'
               #'                    1000       5-7.152-3    1000       6-.013655\n'
                '                    1000       5-7.152-3    1000       6-.013655\n'
                )
        msg16 = (
            'MPC*                   1            1002               4              1.\n'
            '*                   1000               4        -.129394\n'
            '*                                   1000               5        -.007152\n'
            '*                   1000               6        -.013655\n')


        card = model.add_card(card, 'MPC', is_list=False)
        assert card is not None
        mpc = model.mpc
        mpc_id = mpc.add_card(card)
        model.setup(run_geom_check=True)
        msg_8_actual = mpc.write(size=8)
        msg_16_actual = mpc.write(size=16)
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
        mpc = model.mpc
        mpc_id = mpc.add_card(card)
        model.setup(run_geom_check=True)
        msg_8_actual = mpc.write(size=8)
        msg_16_actual = mpc.write(size=16)
        unused_msg_16_double_actual = mpc.write(size=16, is_double=True)
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
        mpc = model.mpc
        mpc_id = model.add_mpc(1, node_ids, components, coefficients, comment='mpc')
        #with self.assertRaises(AssertionError):
            #mpc.validate()
        #model.mpcs = {}

        coefficients = [1000., 1., 101.]
        mpc_id = model.add_mpc(1, node_ids, components, coefficients, comment='mpc')
        model.setup(run_geom_check=True)
        #mpc.raw_fields()
        card = mpc.write(size=8)
        mpc.write(size=16, is_double=False)
        mpc.write(size=16, is_double=True)

        #assert mpc.coefficients == coefficients
        #assert mpc.nodes == node_ids

        model.pop_parse_errors()

        if 0:
            #model.mpcs = {}
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

    def test_spc(self):
        model = BDF(debug=False, mode='nx')
        spc_id = 1
        nodes = [1, 2, 3]
        components = [123, 456, 123456]
        enforced = [0., 0., 0.]
        model.add_spc(spc_id, nodes, components, enforced, comment='spc')

        nodes = [10, 11, 12, 13]
        components = [4, 5, 3, 125]
        enforced = [1., 2., 3., 4.]
        model.add_spc(spc_id, nodes, components, enforced, comment='spc')
        model.setup(run_geom_check=False)
        save_load_deck(model, run_geom_check=False)

    def test_spcadd(self):
        """tests SPCADD"""
        model = BDF(debug=False)
        spc_id = 42
        sets = [1, 2]
        model.add_spcadd(spc_id, sets, comment='spcadd')

        components = ['42', '3', '1']
        #enforced = [1000, 1, 101]
        #node_ids = [2, 3, 4]
        #spc = model.add_spc(1, node_ids, components, enforced, comment='spc')
        #with self.assertRaises(AssertionError):
            #spc.validate()
        model.spc_cards.clear()

        enforced = [1000., 1., 101.]
        node_ids = [2, 3, 4]
        spc_id = model.add_spc(1, node_ids, components, enforced, comment='spc')
        spc = model.spc
        model.spc.validate()
        model.setup()
        #assert spc.enforced == enforced
        assert np.array_equal(spc.node_id, node_ids)

        #spc.raw_fields()
        card = spc.write(size=8)
        #print(card)
        spc.write(size=16, is_double=False)
        spc.write(size=16, is_double=True)
        model.pop_parse_errors()

        model.spc.clear()
        model.spcadd.clear()
        model.add_card(card.split('\n')[1:], 'SPC', is_list=False)
        model.setup(run_geom_check=True)
        #with self.assertRaises(KeyError):
            #model.cross_reference()

        model.add_spc(2, node_ids, components, enforced, comment='spc')
        #with self.assertRaises(KeyError):
            #model.cross_reference()
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        model.add_grid(4, [0., 0., 0.])
        model.setup()
        check_mpc_spc(model)
        save_load_deck(model)
        if RUN_OP2:
            read_write_op2_geom(
                model, run_op2_writer=True, run_op2_reader=True,
                nastran_format='msc')

    def test_spcoff(self):
        """tests SPCOFF/SPCOFF1"""
        model = BDF(debug=False)
        spcoff = model.spcoff
        #spcoff1 = model.spcoff1

        #with self.assertRaises(KeyError):
            #model.EmptyNodes([1, 2], msg='')
        #model.add_spcoff()

        # TOD: doesn't support default component=0
        #card_lines = ['SPCOFF', 1]
        #model.add_card(card_lines, 'SPCOFF', comment='spcoff', is_list=True, has_none=True)
        model.setup()
        spcoff.write(size=8)
        spcoff.write(size=16)
        #spcoff.raw_fields()
        #with self.assertRaises(KeyError):
            #spcoff.cross_reference(model)
        model.add_grid(1, [0., 0., 0.])
        model.setup()
        #spcoff.cross_reference(model)

        card_lines = ['SPCOFF1', 24, 43]
        model.add_card(card_lines, 'SPCOFF1', comment='spcoff1', is_list=True, has_none=True)

        card_lines = ['SPCOFF1', 5, 50, 'THRU', 52]
        model.add_card(card_lines, 'SPCOFF1', comment='spcoff1', is_list=True, has_none=True)
        model.setup()
        model.pop_parse_errors()

        spcoff.write(size=8)
        spcoff.write(size=16)
        #spcoff1.raw_fields()
        #model.uncross_reference()

        #with self.assertRaises(KeyError):
            #spcoff1.cross_reference(model)
        model.add_grid(43, [0., 0., 0.])
        #spcoff1.cross_reference(model)
        #model.uncross_reference()

        model.add_grid(50, [0., 0., 0.])
        model.add_grid(51, [0., 0., 0.])
        model.add_grid(52, [0., 0., 0.])
        model.setup()
        model.cross_reference()
        #model.uncross_reference()
        model.validate()
        #model.safe_cross_reference()

        save_load_deck(model, run_remove_unused=False, run_save_load_hdf5=False)

    def _test_gmspc(self):
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

    def _test_spcax(self):
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

def check_card(msg_expected: str, msg_actual: str):
    if isinstance(msg_expected, tuple):
        msg_expected = msg_expected[0]
    msg_expected_lines = msg_expected.rstrip().split('\n')
    msg_actual_lines = msg_actual.rstrip().split('\n')
    for i, actual, expected  in zip(count(), msg_actual_lines, msg_expected_lines):
        actual = actual.rstrip()
        expected = expected.rstrip()
        msg = 'Error on line %i\n' % i
        msg += 'actual=\n%s\n' % actual
        msg += 'expected=\n%s\n\n' % expected
        msg += 'Actual Card =\n%s\n' % '\n'.join(msg_actual_lines)
        msg += '\nExpected Card =\n%s\n' % '\n'.join(msg_expected_lines)
        assert actual == expected, msg

def check_mpc_spc(model: BDF):
    """simple MPC/SPC checks"""
    mpc_ids = []
    if hasattr(model, 'mpc'):
        mpc_ids = np.unique(np.hstack([mpc.mpc_id for mpc in [model.mpc, model.mpcadd]]))
        mpc_ids = mpc_ids.astype('int32')
    spc_ids = np.unique(np.hstack([spc.spc_id for spc in [model.spc, model.spc1, model.spcadd]]))
    spc_ids = spc_ids.astype('int32')

    #node_ids = []
    #for nid in rbe.independent_nodes + rbe.dependent_nodes:
        #node_ids.append(nid)

    for mpc_id in mpc_ids:
        unused_mpcs1 = model.get_reduced_mpcs(mpc_id, consider_mpcadd=True, stop_on_failure=True)
        unused_mpcs2 = model.get_mpcs(mpc_id, consider_mpcadd=True, stop_on_failure=True)
    for spc_id in spc_ids:
        spcs = model.spcadd.get_reduced_spcs(stop_on_failure=True)
        #unused_spcs1 = model.get_reduced_spcs(spc_id, consider_spcadd=True, stop_on_failure=True)
        #unused_spcs2 = model.get_spcs(spc_id, consider_nodes=False, stop_on_failure=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
