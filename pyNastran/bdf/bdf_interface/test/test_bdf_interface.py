import unittest
import numpy as np
from pyNastran.bdf.bdf import BDF, CaseControlDeck
from pyNastran.bdf.mesh_utils.forces_moments import get_load_arrays
from pyNastran.bdf.mesh_utils.mpc_dependency import get_mpc_node_ids_c1


class TestBDFInterface(unittest.TestCase):
    def test_get_cards_by_card_types(self):
        """tests the ``get_cards_by_card_types`` method"""
        model = BDF(debug=True, log=None, mode='msc')
        card_types = ['GRID']
        with self.assertRaises(RuntimeError):
            cards = model.get_cards_by_card_types(card_types, reset_type_to_slot_map=False,
                                                  stop_on_missing_card=True)

        cards = model.get_cards_by_card_types(card_types, reset_type_to_slot_map=False,
                                              stop_on_missing_card=False)
        assert cards == {'GRID' : []}, cards

        card_types = 'GRID'
        with self.assertRaises(TypeError):
            cards = model.get_cards_by_card_types(card_types, reset_type_to_slot_map=False,
                                                  stop_on_missing_card=False)


        card_types = ['GRID']
        cards = model.get_cards_by_card_types(card_types, reset_type_to_slot_map=True,
                                              stop_on_missing_card=False)
        assert cards == {'GRID' : []}, cards
        #-------------------------------------------------------
        model.add_grid(10, [0., 1., 2.])

        # TODO: not right...
        model.log.error('incorrect result...')
        #cards = model.get_cards_by_card_types(card_types, reset_type_to_slot_map=True,
                                              #stop_on_missing_card=False)
        #assert cards == {'GRID' : [10]}, model.nodes

    def test_get_spcs(self):
        """tests the ``get_spcs`` method"""
        spc_id = 3.9
        model = BDF(debug=True, log=None, mode='msc')
        with self.assertRaises(TypeError):
            nids, comps = model.get_spcs(spc_id, consider_nodes=False, stop_on_failure=True)
        with self.assertRaises(TypeError):
            nids, comps = model.get_spcs(spc_id, consider_nodes=False, stop_on_failure=False)

        spc_id = 42
        with self.assertRaises(KeyError):
            nids, comps = model.get_spcs(spc_id, consider_nodes=False, stop_on_failure=True)
        nids, comps = model.get_spcs(spc_id, consider_nodes=False, stop_on_failure=False)

        components = '53'
        nodes = [10]
        model.add_spc1(spc_id, components, nodes, comment='')

        model.add_grid(100, [0., 0., 0.], cp=0, cd=0, ps='6', seid=0, comment='')
        nids, comps = model.get_spcs(spc_id, consider_nodes=False, stop_on_failure=True)
        assert nids == [10], nids
        assert comps == ['53'], comps

        nids, comps = model.get_spcs(spc_id, consider_nodes=True, stop_on_failure=True)
        assert nids == [10, 100], nids
        assert comps == ['53', '6'], comps
        #nids, comps = model.get_spcs(spc_id, consider_nodes=False, stop_on_failure=True)
        #nids, comps = model.get_spcs(spc_id, consider_nodes=False, stop_on_failure=True)
        #nids, comps = model.get_spcs(spc_id, consider_nodes=False, stop_on_failure=True)


    def test_mpcs(self):
        """tests the ``get_MPCx_node_ids_c1`` method"""
        mpc_id = 3.9
        model = BDF(debug=True, log=None, mode='msc')

        # stop_on_failure stops if the data is formatted improperly
        with self.assertRaises(TypeError):
            get_mpc_node_ids_c1(model, mpc_id, consider_mpcadd=True, stop_on_failure=True)
        with self.assertRaises(TypeError):
            get_mpc_node_ids_c1(model, mpc_id, consider_mpcadd=True, stop_on_failure=False)

        mpc_id = 5
        with self.assertRaises(KeyError):
            get_mpc_node_ids_c1(model, mpc_id, consider_mpcadd=True, stop_on_failure=True)

        out = get_mpc_node_ids_c1(model, mpc_id, consider_mpcadd=True, stop_on_failure=False)
        independent_node_ids_c1, dependent_node_ids_c1 = out
        assert independent_node_ids_c1 == {}, independent_node_ids_c1
        assert dependent_node_ids_c1 == {}, dependent_node_ids_c1

        mpc_id = 10
        nodes = [11]
        components = ['4']
        coefficients = [1.]
        mpc = model.add_mpc(mpc_id, nodes, components, coefficients)
        mpc.validate()
        #print(mpc)
        out = get_mpc_node_ids_c1(model, mpc_id, consider_mpcadd=True, stop_on_failure=False)
        independent_node_ids_c1, dependent_node_ids_c1 = out
        #print('independent_node_ids_c1 =', independent_node_ids_c1)
        #print('dependent_node_ids_c1 =', dependent_node_ids_c1)
        assert independent_node_ids_c1 == {'4' : [11],}, independent_node_ids_c1
        assert dependent_node_ids_c1 == {}, dependent_node_ids_c1

    def test_loads(self):
        """tests the ``get_load_arrays`` and ``_reduce_dload_case`` methods"""
        model = BDF(debug=True, log=None, mode='msc')
        subcase_id = 10
        eid_map = None
        node_ids = []
        normals = None
        lines = []
        cc = CaseControlDeck(lines, log=None)
        #print(model.case_control_deck)
        model.case_control_deck = cc
        subcase = cc.create_new_subcase(subcase_id)
        out = get_load_arrays(model, subcase_id, eid_map, node_ids, normals,
                              nid_map=None)
        is_loads, is_temperatures, temperature_data, load_data = out
        assert is_loads is False, is_loads
        assert is_temperatures is False, is_temperatures
        assert temperature_data == (None, None), temperature_data
        assert load_data == (None, None, None), load_data

        key = 'LOAD'
        value = 10
        options = []
        param_type = '???'
        with self.assertRaises(TypeError):
            subcase.add(key, value, options, param_type)

        param_type = 'STRESS-type'
        subcase.add(key, value, options, param_type)

        with self.assertRaises(KeyError):
            get_load_arrays(model, subcase_id, eid_map, node_ids, normals, nid_map=None,
                            stop_on_failure=True)
        #is_loads, is_temperatures, temperature_data, load_data = out
        #assert is_loads is False, is_loads
        #assert is_temperatures is False, is_temperatures
        #assert temperature_data == (None, None), temperature_data
        #assert load_data == (None, None, None), load_data

        dload_case = []
        dloads, scale_factors = model._reduce_dload_case(
            dload_case, scale=1., unallowed_dload_ids=None, skip_scale_factor0=False, msg='')
        assert dloads == [], dloads
        assert scale_factors == [], scale_factors
        del dloads, scale_factors

        # ----------------------------------------------------------------------
        dload_id = 5
        scale = 1.
        scale_factors = [2]
        load_ids = [3]
        dload = model.add_dload(dload_id, scale, scale_factors, load_ids)
        dload_case = [dload]

        with self.assertRaises(KeyError):
            model._reduce_dload_case(
                dload_case, scale=1., unallowed_dload_ids=None, skip_scale_factor0=False, msg='')

        excite_id = 100
        tid = 2
        model.add_tload1(dload_id, excite_id, tid, delay=0, Type='LOAD', us0=0.0, vs0=0.0, comment='')
        with self.assertRaises(KeyError):
            model._reduce_dload_case(
                dload_case, scale=1., unallowed_dload_ids=None, skip_scale_factor0=False, msg='')
        #assert dloads == [], dloads
        #assert scale_factors == [], scale_factors

    def test_get_mklist(self):
        """tests the ``get_mklist`` method"""
        model = BDF(debug=True, log=None, mode='msc')
        mklist1 = model.get_mklist()
        assert np.array_equal(mklist1, []), mklist1

        machs = [0.8, 0.9, 0.95]
        reduced_freqs = [1., 2.]
        model.add_mkaero1(machs, reduced_freqs)
        mklist2 = model.get_mklist()
        mklist2_expected = [
            [0.80, 1.],
            [0.80, 2.],
            [0.90, 1.],
            [0.90, 2.],
            [0.95, 1.],
            [0.95, 2.],
        ]
        assert np.array_equal(mklist2, mklist2_expected), mklist2

if __name__ == '__main__':   # pragma: no cover
    unittest.main()
