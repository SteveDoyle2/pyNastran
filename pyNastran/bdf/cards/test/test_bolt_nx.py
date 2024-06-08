import unittest
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck

class TestBolt(unittest.TestCase):
    def test_bolt_nx_1(self):
        model = BDF(debug=False)
        model.case_control_deck
        model.subcases

        bolt_id = 1
        element_type = 1
        eids = [11, 12, 13, 14, 15, 16]
        bolt = model.add_bolt_nx(bolt_id, element_type, eids=eids, csid=None, idir=None)
        str(bolt)

        bolt_id = 2
        element_type = 2
        eids = None
        nids = [101, 102, 103, 104, 105, 106]
        bolt = model.add_bolt_nx(bolt_id, element_type, nids=nids, csid=None, idir=1)
        str(bolt)

        bolt_id = 3
        element_type = 3
        eids = [11, 12, 13, 14, 15, 16]
        bolt = model.add_bolt_nx(bolt_id, element_type, eids=eids, csid=None, idir=None)
        str(bolt)

        #SID Bolt preload set identification number (Integer; No default)
        #S_NOi Sequence order number for the BOLTLD, BOLTFOR, and BOLTFRC IDs to be applied. (Integer; No default)
        #B_IDi SID of BOLTLD, BOLTFOR, or BOLTFRC bulk entries defining a boltpreload. (Integer; No default)
        #NINCi Number of increments in which to ramp up the bolt preload defined in

        sid = 4
        s_no = 5
        b_id = 6
        boltseq = model.add_boltseq_nx(sid, s_no, b_id, n_incs=1, comment='')

        sid = 5
        s_no = 5
        b_id = 6
        boltseq = model.add_boltseq_nx(sid, [s_no, s_no], [b_id, b_id], n_incs=1, comment='')

        load_value = 10.
        bolt_ids = [101, 102, 103, 104, 105, 106, 107]
        boltfor = model.add_boltfor_nx(sid, load_value, bolt_ids)

        for nid in nids:
            model.add_grid(nid, [0., 0., 0.])
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        pid = 10
        mid = 11
        x = [0., 0., 1.]
        g0 = None
        nidsi = [1, 2]
        for eid in eids:
            model.add_cbar(eid, pid, nidsi, x, g0, offt='GGG', pa=0, pb=0, wa=None, wb=None, comment='', validate=False)
        model.add_pbarl(pid, mid, 'BAR', [1., 2.], group='MSCBML0', nsm=0., comment='')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model2 = save_load_deck(model, xref='standard', punch=True, run_remove_unused=True, run_convert=True, run_renumber=False,
                                run_mirror=True, run_save_load=True, run_quality=True, write_saves=True, run_save_load_hdf5=True,
                                run_mass_properties=True, run_loads=True, run_test_bdf=True, run_op2_writer=True, run_op2_reader=True,
                                remove_disabled_cards=True, nastran_format='nx', op2_log_level='warning')
        boltfor = model2.boltfor[5]
        assert len(boltfor.bolt_ids) == 7, boltfor.bolt_ids
        #model.add_boltfor_nx()
        #model.add_boltd_nx()
