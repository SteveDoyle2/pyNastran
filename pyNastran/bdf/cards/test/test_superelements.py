import unittest
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck


class TestSuperelements(unittest.TestCase):
    def test_superelements_1(self):
        """SEMPLN/SELOC/SEBULK test"""
        model = BDF()
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
        super1 = BDF()
        super1.add_grid(10, [0., 0., 0.])
        super1.add_grid(20, [0., 0., 1.])
        super1.add_grid(30, [1., 0., 0.])

        super1.add_grid(11, [1., 1., 0.])
        super1.add_grid(12, [2., 2., 0.])

        x = [0., 0., 1.]
        g0 = None
        super1.add_cbar(100, 1000, [11, 12], x, g0)
        super1.add_pbarl(1000, 2000, 'ROD', [1.,])
        super1.add_mat1(2000, 3.0e7, None, 0.3)
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
        #save_load_deck(model, punch=True)

if __name__ == '__main__':
    unittest.main()
