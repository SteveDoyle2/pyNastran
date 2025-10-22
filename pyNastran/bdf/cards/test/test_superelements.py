import os
from pathlib import Path
import unittest
import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck, _run_mass_properties, _run_loads, _run_hdf5

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / '..' / 'models'


class TestSuperelements(unittest.TestCase):

    def test_seconct(self):
        model = BDF(debug=None)
        model.sol = 103

        super_a = BDF(debug=None)
        super_a.add_grid(11, [1., 0., 0.])
        super_a.add_grid(12, [2., 0., 0.])
        super_a.add_grid(13, [3., 0., 0.])

        super_a.add_grid(101, [1., 0., 0.])
        super_a.add_grid(102, [2., 0., 0.])
        super_a.add_grid(103, [3., 0., 0.])
        seid_a = 10
        seid_b = 10
        model.superelement_models[('SUPER', 10, '')] = super_a

        tol = 0.1
        loc = 'NO'
        nodes_a = [11, 12, 13]
        nodes_b = [101, 102, 103]

        seconct = model.add_seconct(
            seid_a, seid_b,
            nodes_a, nodes_b,
            tol=tol, loc=loc, comment='seconct')
        seconct.raw_fields()
        model.cross_reference()
        save_load_deck(model)

    def test_seelt(self):
        model = BDF(debug=None)
        model.sol = 101

        super_a = BDF(debug=None)
        model.superelement_models[('SUPER', 10, '')] = super_a

        seid = 10
        eids = [4, 5]
        seelt = super_a.add_seelt(seid, eids, comment='seelt')
        super_a.add_grid(11, [1., 0., 0.])
        super_a.add_grid(12, [2., 0., 0.])
        super_a.add_grid(13, [3., 0., 0.])
        super_a.add_conrod(4, 100, [11, 12], A=1.0)
        super_a.add_conrod(5, 100, [12, 13], A=1.0)
        super_a.add_mat1(100, 3.0e7, None, 0.3)
        seelt.raw_fields()
        model.cross_reference()
        save_load_deck(model)

    def test_setree(self):
        model = BDF(debug=None)
        model.sol = 101

        seid_a = 10
        seid_b = seid_a + 1
        super_a = BDF(debug=None)
        super_b = BDF(debug=None)
        super_a.add_grid(11, [1., 0., 0.])
        super_b.add_grid(21, [1., 0., 0.])
        model.superelement_models[('SUPER', seid_a, '')] = super_a
        model.superelement_models[('SUPER', seid_b, '')] = super_b

        seid_c = 50
        seids = [seid_a, seid_b]
        setree = model.add_setree(seid_c, seids, comment='')
        setree.raw_fields()
        model.cross_reference()
        save_load_deck(model)

    def test_release(self):
        model = BDF(debug=None)
        comp = '3'
        psid = 42
        nodes = [11, 12, 13]

        label = 'CAT'
        comment = 'release'
        for seid in [10, 11]:
            release = model.add_release(seid, comp, nodes, comment=comment)
            csuper = model.add_csuper(
                seid, psid, nodes,
                comment=comment)
            csupext = model.add_csupext(seid, nodes, comment=comment)
            selabel = model.add_selabel(seid, label, comment=comment)
            comment = ''

        comps = ['A', 'B', 'C']
        sesup = model.add_sesup(nodes, comps, comment='sesup')

        model.add_grid(11, [1., 0., 0.])
        model.add_grid(12, [2., 0., 0.])
        model.add_grid(13, [3., 0., 0.])
        release.raw_fields()
        csuper.raw_fields()
        csupext.raw_fields()
        selabel.raw_fields()
        sesup.raw_fields()
        model.cross_reference()
        save_load_deck(model)

    def _test_superelements_pch(self):
        model = BDF(mode='nx')
        model.is_superelements = True
        bdf_filename = MODEL_PATH / 'bugs' / 'outboard_op4asmblk.pch'
        model.read_bdf(bdf_filename, punch=True)
        model.sol = 103
        save_load_deck(model)

    def test_superelements_1(self):
        """SEMPLN/SELOC/SEBULK test"""
        model = BDF(debug=False)
        model.sol = 103
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
        super1, super1_key = create_superelement(1, debug=False)
        model.superelement_models[super1_key] = super1

        #----------------------------------------------

        model.validate()
        #model.cross_reference(
            #xref=True, xref_nodes=True, xref_elements=True, xref_nodes_with_elements=False,
            #xref_properties=True, xref_masses=True, xref_materials=True,
            #xref_loads=True, xref_constraints=True, xref_aero=True,
            #xref_sets=True, xref_optimization=True, word='')
        #model.uncross_reference()

        assert len(model.superelement_models) == 1, model.superelement_models
        model.safe_cross_reference(
            xref=True, xref_nodes=True, xref_elements=True, xref_nodes_with_elements=False,
            xref_properties=True, xref_masses=True, xref_materials=True,
            xref_loads=True, xref_constraints=True, xref_aero=True,
            xref_sets=True, xref_optimization=True,
            create_superelement_geometry=True, debug=True, word='')
        model.uncross_reference()
        model.cross_reference()
        os.remove('super_2.bdf')
        save_load_deck(model)

    def test_superelement_2(self):
        """
        tests the following cards:
         - SEBNDRY
         - SELOAD
         - SEEXCLD
        """
        model = BDF(debug=None)
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
        model = BDF(debug=None)
        model.sol = 103
        seid_a = 1
        seid_b = 2
        nodes = [10, 11, 12]

        super1, super1_key = create_superelement(seid_a, debug=False)
        super2, super2_key = create_superelement(seid_b, debug=False)
        model.superelement_models[super1_key] = super1
        model.superelement_models[super2_key] = super2

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
        save_load_deck(model, run_save_load_hdf5=False)

    def test_superelement_setree(self):
        """tests the SETREE"""
        model = BDF(debug=False)
        model.sol = 101
        super1, super1_key = create_superelement(101, debug=True)
        super2, super2_key = create_superelement(102, debug=True)
        model.superelement_models[super1_key] = super1
        model.superelement_models[super2_key] = super2

        seid = 100
        seids = [101, 102]
        setree = model.add_setree(seid, seids, comment='setree')
        setree.raw_fields()
        model.validate()
        save_load_deck(model)

    def test_super_sets(self):
        model = BDF(debug=False, log=None, mode='msc')
        model.add_aset([1, 2, 3], '456', comment='aset')
        model.add_aset1([1, 2, 3], '123', comment='aset1')

        model.add_bset([3, 4], '456', comment='bset')
        model.add_bset1([3, 4], '123', comment='bset1')

        model.add_cset([5, 6], '456', comment='cset')
        model.add_cset1([5, 6], '123', comment='cset1')

        #model.add_omit([7, 8], '456', comment='aset')
        model.add_omit1([7, 8], '123', comment='omit1')

        model.add_sebset(103, [5, 6], '456', comment='sebset')
        model.add_sebset1(103, [5, 6], '123', comment='sebset1')
        model.add_secset(103, [5, 6], '456', comment='secset')
        model.add_secset1(103, [5, 6], '123', comment='secset1')

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        model.add_grid(4, [0., 0., 0.])
        model.add_grid(5, [0., 0., 0.])
        model.add_grid(6, [0., 0., 0.])
        model.add_grid(7, [0., 0., 0.])
        model.add_grid(8, [0., 0., 0.])
        save_load_deck(model)

def create_superelement(super_id: int,
                        debug=False) -> tuple[BDF, tuple[str, int, str]]:
    """creates a simple bar model"""
    super_model = BDF(debug=debug)
    super_model.add_grid(10, [0., 0., 0.])
    super_model.add_grid(20, [0., 0., 1.])
    super_model.add_grid(30, [1., 0., 0.])

    super_model.add_grid(11, [1., 1., 0.])
    super_model.add_grid(12, [2., 2., 0.])

    x = [0., 0., 1.]
    g0 = None
    super_model.add_cbar(100, 1000, [11, 12], x, g0)
    super_model.add_pbarl(1000, 2000, 'ROD', [1.,])
    super_model.add_mat1(2000, 3.0e7, None, 0.3)
    super_key = ('SUPER', super_id, '')
    return super_model, super_key


if __name__ == '__main__':   # pragma: no cover
    unittest.main()
