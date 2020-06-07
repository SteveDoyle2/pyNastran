import os
import unittest
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck, _run_mass_properties, _run_loads, _run_hdf5


class TestSuperelements(unittest.TestCase):
    def test_superelements_1(self):
        """SEMPLN/SELOC/SEBULK test"""
        model = BDF(debug=False)
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 1.])
        model.add_grid(4, [1., 0., 0.])

        # we're going to mirror superelement 1 into superelement 2
        seid = 2
        p1 = 2
        p2 = 3
        p3 = 4
        nodes0 = [2, 3, 4]

        #  these come from seid=1 after the mirror and shift?
        nodes_seid = [10, 20, 30]
        model.add_sempln(seid, p1, p2, p3, comment='sempln')
        model.add_seloc(seid, nodes_seid, nodes0, comment='seloc')

        seid = 2
        superelement_type = 'MIRROR'
        rseid = 1
        model.add_sebulk(seid, superelement_type, rseid, comment='sebulk')

        #----------------------------------------------
        super1 = create_superelement(debug=False)
        model.superelement_models[1] = super1

        #----------------------------------------------

        model.validate()
        #model.cross_reference(
            #xref=True, xref_nodes=True, xref_elements=True, xref_nodes_with_elements=False,
            #xref_properties=True, xref_masses=True, xref_materials=True,
            #xref_loads=True, xref_constraints=True, xref_aero=True,
            #xref_sets=True, xref_optimization=True, word='')
        #model.uncross_reference()

        model.safe_cross_reference(
            xref=True, xref_nodes=True, xref_elements=True, xref_nodes_with_elements=False,
            xref_properties=True, xref_masses=True, xref_materials=True,
            xref_loads=True, xref_constraints=True, xref_aero=True,
            xref_sets=True, xref_optimization=True,
            create_superelement_geometry=True, debug=True, word='')
        os.remove('super_2.bdf')
        #save_load_deck(model, punch=True)

    def test_superelement_2(self):
        """
        tests the following cards:
         - SEBNDRY
         - SELOAD
         - SEEXCLD
        """
        model = BDF(debug=False)
        seid_a = 101
        seid_b = 102
        ids = [10, 20]
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(20, [0., 0., 1.])
        sebndry = model.add_sebndry(seid_a, seid_b, ids, comment='sebndry')

        lid_s0 = 1
        seid = 1
        lid_se = 1
        seload = model.add_seload(lid_s0, seid, lid_se, comment='seload')

        nodes = ids
        seexcld = model.add_seexcld(seid_a, seid_b, nodes, comment='seexcld')

        set_id = 42
        n = 43
        senqset = model.add_senqset(set_id, n, comment='senqset')
        sebndry.raw_fields()
        seload.raw_fields()
        seexcld.raw_fields()
        senqset.raw_fields()

        model.validate()

        save_load_deck(model, run_save_load_hdf5=False)

    def test_seexclude(self):
        model = BDF(debug=False)
        seid_a = 1
        seid_b = 2
        nodes = [10, 11, 12]

        super1 = create_superelement(debug=False)
        super2 = create_superelement(debug=False)
        model.superelement_models[seid_a] = super1
        model.superelement_models[seid_b] = super2

        card_fields = ['SEEXCLUD', seid_a, seid_b, 11, 30]
        seexcld = model.add_card_fields(card_fields, 'SEEXCLUD', comment='',
                                        has_none=True)

        seexcld = model.add_seexcld(seid_a, seid_b, nodes, comment='seexclud')
        seexcld.raw_fields()
        model.validate()

        nelements = len(model.elements) + len(model.masses)
        nnodes = len(model.nodes) + len(model.spoints) + len(model.epoints)
        _run_mass_properties(model, nnodes, nelements, run_mass_properties=True)
        _run_loads(model, nelements, run_loads=True)
        #_run_hdf5(model, model.log, run_save_load_hdf5=True)

        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()

        #save_load_deck(model, run_save_load_hdf5=False)

    def test_superelement_setree(self):
        """tests the SETREE"""
        model = BDF(debug=False)
        super1 = create_superelement(debug=True)
        super2 = create_superelement(debug=True)
        model.superelement_models[101] = super1
        model.superelement_models[102] = super2

        seid = 100
        seids = [101, 102]
        setree = model.add_setree(seid, seids, comment='setree')
        setree.raw_fields()
        model.validate()
        #TODO: enable this...fix error
        #save_load_deck(model, run_test_bdf=False, run_save_load_hdf5=False)

def create_superelement(debug=False):
    """creates a simple bar model"""
    super = BDF(debug=debug)
    super.add_grid(10, [0., 0., 0.])
    super.add_grid(20, [0., 0., 1.])
    super.add_grid(30, [1., 0., 0.])

    super.add_grid(11, [1., 1., 0.])
    super.add_grid(12, [2., 2., 0.])

    x = [0., 0., 1.]
    g0 = None
    super.add_cbar(100, 1000, [11, 12], x, g0)
    super.add_pbarl(1000, 2000, 'ROD', [1.,])
    super.add_mat1(2000, 3.0e7, None, 0.3)
    return super

if __name__ == '__main__':   # pragma: no cover
    unittest.main()
