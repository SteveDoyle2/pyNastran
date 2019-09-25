import unittest
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck


class TestParametric(unittest.TestCase):
    def test_parametric(self):
        """tests xxx"""

        model = BDF(debug=False)
        idi = 10
        poly1 = 11
        poly2 = 12
        poly3 = 13
        cid = 0
        typei = 'cat'
        typeids = 15  #[-1, -1]
        pset = model.add_pset(idi, poly1, poly2, poly3, cid, typei, typeids,
                              comment='pset')
        #pval = model.add_pval(idi, poly1, poly2, poly3, cid, typei, typeids,
                              #comment='pval')

        edge_id = 100
        nids = [101, 102]
        geom_ids = [-1, -1]
        feedge = model.add_feedge(edge_id, nids, cid, geom_ids,
                                 geomin='POINT')
        #feface = model.add_feface()

        curve_id = 200
        group = 'MSCGRP0'
        data = ['some_string']
        gmcurv = model.add_gmcurv(curve_id, group, data, cid_in=0, cid_bc=0)
        #gmsurf = model.add_gmsurf()

        pset.raw_fields()
        #pval.raw_fields()
        feedge.raw_fields()
        gmcurv.raw_fields()

        #model.write_bdf
        save_load_deck(model, xref='standard', punch=True,
                      run_remove_unused=False,
                      run_convert=False,
                      run_renumber=False,
                      run_mirror=False,
                      run_save_load=True,
                      run_quality=False,
                      write_saves=True,
                      run_save_load_hdf5=True,
                      run_mass_properties=False,
                      run_loads=False,
                      run_test_bdf=False,
                      run_op2_writer=False,
                      run_op2_reader=False)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
