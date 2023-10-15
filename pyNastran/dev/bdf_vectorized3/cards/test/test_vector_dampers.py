import unittest
from io import StringIO
import numpy as np

from cpylog import get_logger
from pyNastran.dev.bdf_vectorized3.bdf import BDF, BDFCard, read_bdf#, get_logger2
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck
from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_shells import (
    make_dvprel_optimization,
    make_dvcrel_optimization,
    make_dvmrel_optimization,
)


class TestDampers(unittest.TestCase):
    def test_damper_01(self):
        """tests PDAMP"""
        lines = ['pdamp, 201, 1.e+5']
        log = get_logger(level='warning')
        model = BDF(log=log)
        card = model._process_card(lines)
        card = BDFCard(card)

        size = 8
        elem_ids = model.pdamp.add_card(card)
        model.setup()

        elem = model.pdamp
        elem.write(size, 'dummy')
        #elem.raw_fields()
        self.assertEqual(elem.property_id, 201)
        self.assertEqual(elem.b, 1e5)

        fields = ['pdamp', 201, 1.e+5, None, None, 202, 2.e+5]
        card_name = fields[0]
        model.add_card(fields, card_name, comment='', is_list=True,
                       has_none=True)
        assert len(model.pdamp) == 3, model.pdamp
        save_load_deck(model)

    def test_cdamp1_01(self):
        """tests a CDAMP1"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        lines = ['CDAMP1, 2001, 20, 1001, 1']
        card = model._process_card(lines)
        card = BDFCard(card)

        size = 8
        pid = 20
        b = 1.0

        elem_id = model.cdamp1.add_card(card)
        model.add_pdamp(pid, b)
        model.setup()

        elem = model.cdamp1
        self.assertEqual(elem.element_id, 2001)
        self.assertEqual(elem.property_id, 20)
        node_ids = elem.nodes[0]
        assert np.array_equal(node_ids, [1001, 0]), node_ids
        elem.write(size, 'dummy')
        #elem.raw_fields()

        pid = 1
        tbid = 2
        model.add_pdampt(pid, tbid, comment='pdampt')
        save_load_deck(model)

    def test_damper_02(self):
        """tests CDAMP1, CDAMP2, PDAMP, PDAMPT, GRID"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 1
        pid = 2
        nids = [3, 4]
        c1 = 1
        c2 = 1
        cpdamp1_id = model.add_cdamp1(eid, pid, nids, c1, c2, comment='cdamp1')

        b = 1.0e7
        pdamp_id = model.add_pdamp(pid, b, comment='pdamp')

        tbid = 10
        pdampt_id = model.add_pdampt(pid, tbid, comment='pdampt')

        eid = 5
        cdamp2_id = model.add_cdamp2(eid, b, nids, comment='cdamp2')
        #cdamp1.raw_fields()
        #pdamp.raw_fields()
        #pdampt.raw_fields()
        #cdamp2.raw_fields()
        model.setup()

        model.add_grid(3, [0., 0., 0.])
        model.add_grid(4, [0., 0., 0.])
        model.validate()
        model._verify_bdf(xref=False)
        model.cross_reference()
        model._verify_bdf(xref=True)

        model.cdamp1.write(size=8, is_double=False)
        model.pdamp.write(size=8, is_double=False)
        model.pdampt.write(size=8, is_double=False)
        model.cdamp2.write(size=8, is_double=False)

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.seek(0)
        unused_model2 = read_bdf(bdf_file, punch=True, debug=False)
        save_load_deck(model)

    def test_cdamp5(self):
        log = get_logger(level='warning')
        model = BDF(log=log)

        eid = 5
        pid = 5
        mid = 10
        s1 = 10
        s2 = 20
        bdamp = 3.0e3
        pdamp5_id = model.add_pdamp5(pid, mid, bdamp, comment='pdamp5')
        cdamp5_id = model.add_cdamp5(eid, pid, [s1, s2], comment='cdamp5')
        model.cdamp5.write(size=8, is_double=False)
        model.pdamp5.write(size=8, is_double=False)

    def test_damper_03(self):
        """tests the CDAMP4, PDAMP, CDAMP4, SPOINT"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 3
        pid = 2
        s1 = 3
        s2 = 4
        cdamp3_id = model.add_cdamp3(eid, pid, [s1, s2], comment='cdamp3')

        bdamp = 1.0e3
        pdamp_id = model.add_pdamp(pid, bdamp, comment='pdamp')
        spoints_id = model.add_spoint([3, 4, 5], comment='spoints')

        eid = 4
        bdamp = 2.0e3
        s1 = 3
        s2 = 4
        cdamp4_id = model.add_cdamp4(eid, bdamp, [s1, s2], comment='cdamp4')

        # two new CDAMP4s defined on one line
        bdamp2 = 3.0e3
        fields = ['CDAMP4', 104, bdamp, s1, s2, 105, bdamp2, s1, s2]
        model.add_card(fields, 'CDAMP4')

        eid = 6
        pid = 6
        ce = 1.
        cr = 1.
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [0., 0., 0.])
        nids = [10, 11]
        pvisc_id = model.add_pvisc(pid, ce, cr, comment='pvisc')
        cvisc_id = model.add_cvisc(eid, pid, nids, comment='cvisc')

        eid = 7
        pid = 7
        x = [1., 2., 3.]
        g0 = None
        cgap_id = model.add_cgap(eid, pid, nids, x, g0, cid=None, comment='cgap')
        pgap_id = model.add_pgap(pid, u0=0., f0=0., ka=1.e8, kb=None, mu1=0.,
                                 kt=None, mu2=None,
                                 tmax=0., mar=100., trmin=0.001,
                                 comment='pgap')

        eid = 8
        pid = 8
        k = [1.0]
        b = [2.0]
        ge = [0.01]
        cbush_id = model.add_cbush(eid, pid, nids, x, g0, cid=None, s=0.5,
                                   ocid=-1, si=None, comment='cbush')
        pbush_id = model.add_pbush(pid, k, b, ge, rcv=None, mass=None,
                                   comment='pbush')

        eid = 9
        pid = 9
        k = 1.
        c = 2.
        m = 3.
        sa = 4.
        se = 5.
        #optional_vars = None
        cbush1d_id = model.add_cbush1d(eid, pid, nids, cid=None, comment='cbush1d')
        pbush1d_id = model.add_pbush1d(pid, k=k, c=c, m=m, sa=sa, se=se,
                                       optional_vars=None, comment='pbush1d')

        card_lines = [
            'pbush1d, 204, 1.e+5, 1000., , , , , , +pb4',
            '+pb4, shocka, table, 1000., , 1., , 214, , +pb41',
            '+pb41, spring, table, 205',
        ]
        model.add_card_lines(card_lines, 'PBUSH1D', comment='', has_none=True)
        prop = model.Property(204)
        str(prop)

        mid = 10
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.validate()

        #cdamp3.raw_fields()
        #cdamp4.raw_fields()
        #cdamp5.raw_fields()
        #pdamp5.raw_fields()
        #cvisc.raw_fields()
        #pvisc.raw_fields()
        #cgap.raw_fields()
        #pgap.raw_fields()
        #cbush.raw_fields()
        #pbush.raw_fields()
        #cbush1d.raw_fields()
        #pbush1d.raw_fields()
        #spoints.raw_fields()

        model.cdamp3.write(size=8, is_double=False)
        model.cdamp4.write(size=8, is_double=False)
        model.cvisc.write(size=8, is_double=False)
        model.pvisc.write(size=8, is_double=False)
        model.cgap.write(size=8, is_double=False)
        model.pgap.write(size=8, is_double=False)
        model.cbush.write(size=8, is_double=False)
        model.pbush.write(size=8, is_double=False)
        model.cbush1d.write(size=8, is_double=False)
        model.pbush1d.write(size=8, is_double=False)

        run_opt = True
        run_opt_dvc = True
        if run_opt:
            params = [
                ('K1', 1.0), ('K2', 1.0), ('K3', 1.0), ('K4', 1.0), ('K5', 1.0), ('K6', 1.0),
                ('B1', 1.0), ('B2', 1.0), ('B3', 1.0), ('B4', 1.0), ('B5', 1.0), ('B6', 1.0),
                #('M1', 1.0), ('M2', 1.0), ('M3', 1.0), ('M4', 1.0), ('M5', 1.0), ('M6', 1.0),
            ]
            i = make_dvprel_optimization(model, params, 'PBUSH', pbush_id, i=1)

            params = [(5, 1.0)]
            i = make_dvprel_optimization(model, params, 'PGAP', pgap_id, i)

            params = [('K', 1.0), ('C', 1.0), ('M', 1.0)]
            i = make_dvprel_optimization(model, params, 'PBUSH1D', pbush1d_id, i)

            params = [('CE1', 1.0)]
            i = make_dvprel_optimization(model, params, 'PVISC', pvisc_id, i)

            params = [('B1', 1.0), (3, 1.0)]
            i = make_dvprel_optimization(model, params, 'PDAMP', pdamp_id, i)

            #-----------------------------------------
            if run_opt_dvc:
                params = []
                i = make_dvcrel_optimization(model, params, 'CVISC', cvisc_id, i)

                params = [('X1', 1.0), ('X2', 2.0), ('X3', 3.0),
                          ('S1', 1.0), ('S2', 2.0), ('S3', 3.0),
                          ('S', 1.0), ]
                i = make_dvcrel_optimization(model, params, 'CBUSH', cbush_id, i)

                params = []
                i = make_dvcrel_optimization(model, params, 'CBUSH1D', cbush_id, i)

                params = []
                i = make_dvcrel_optimization(model, params, 'CGAP', cgap_id, i)

        model.setup()
        model.spoint.write()

        model.cross_reference()
        run_opt_dvc2 = False
        if run_opt_dvc2:
            model.update_model_by_desvars()
        #assert 204 in model.properties, model.properties
        model.dvcrel1.write(size=8)
        model.dvprel1.write(size=8)
        save_load_deck(model, run_read_write=False)

    def test_pbusht(self):
        """tests CBUSH, PBUSH, PBUSHT"""
        model = BDF(debug=False, log=None, mode='msc')
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [0., 0., 0.])

        eid = 8
        pid = 8
        k = [1.0]
        b = [2.0]
        ge = [0.01]
        nids = [10, 11]
        x = [1., 0., 0.]
        g0 = None
        unused_cbush = model.add_cbush(eid, pid, nids, x, g0, cid=None, s=0.5,
                                       ocid=-1, si=None, comment='cbush')
        unused_pbush = model.add_pbush(pid, k, b, ge, rcv=None, mass=None,
                                       comment='pbush')

        k_tables = [2]
        b_tables = [2]
        ge_tables = [2]
        kn_tables = [2]
        pbusht_id = model.add_pbusht(pid, k_tables, b_tables, ge_tables, kn_tables,
                                  comment='')
        #pbusht.raw_fields()
        model.validate()
        model._verify_bdf()
        model.cross_reference()
        model._verify_bdf()
        save_load_deck(model, xref='standard', punch=True)


    def test_pdamp(self):
        """PDAMP"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid1 = 10
        eid2 = 20
        eid3 = 30
        eid4 = 40
        b1 = 1.0
        b2 = 2.0
        b3 = 3.0
        b4 = 4.0
        #nodes1 = [10, 20]
        #nodes2 = [20, 30]
        card_lines = ['PDAMP', eid1, b1, eid2, b2, eid3, b3, eid4, b4]
        model.add_card(card_lines, 'PDAMP', comment='', is_list=True, has_none=True)
        model.validate()
        model._verify_bdf()
        save_load_deck(model)

    def test_cgap_coord(self):
        """tests a CGAP/PGAP"""
        model = BDF(debug=False)
        eid = 1
        pid = 2
        nids = [10, 11]
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        x = [0., 0., 1.]
        g0 = None
        model.add_cgap(eid, pid, nids, x, g0, cid=0, comment='')
        model.add_pgap(pid, u0=0., f0=0., ka=1.e8, kb=None, mu1=0.,
                       kt=None, mu2=None, tmax=0., mar=100., trmin=0.001, comment='')
        save_load_deck(model)

    def test_cgap(self):
        """tests a CGAP/PGAP"""
        model = BDF(debug=False)
        eid = 1
        pid = 2
        nids = [10, 11]
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        x = [0., 0., 1.]
        g0 = None
        model.add_cgap(eid, pid, nids, x, g0, cid=None, comment='')
        model.add_pgap(pid, u0=0., f0=0., ka=1.e8, kb=None, mu1=0.,
                       kt=None, mu2=None, tmax=0., mar=100., trmin=0.001, comment='')

        eid += 1
        pid += 1
        k = 1.0
        b = 2.0
        ge = 0.0
        model.add_cbush(eid, pid, nids, x, g0, cid=None, s=0.5, ocid=-1, si=None, comment='')
        model.add_pbush(pid, k, b, ge, rcv=None, mass=None, alpha=0., tref=0.,
                        coincident_length=None, comment='')
        model.setup(run_geom_check=True)
        save_load_deck(model)

    def test_cvisc_01(self):
        """tests a CVISC/PVISC"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        size = 8
        eid = 2001
        nids = [4, 1]
        pid = 20
        ce = 1.0
        cr = 2.0

        model.add_grid(4, [0., 0., 0.])
        model.add_grid(1, [1., 0., 0.])
        model.add_cvisc(eid, pid, nids, comment='cvisc')
        model.add_pvisc(pid, ce, cr, comment='pvisc')
        model.setup()

        elem = model.cvisc
        self.assertEqual(elem.element_id, 2001)
        self.assertEqual(elem.property_id, 20)
        node_ids = elem.nodes[0]
        assert np.array_equal(node_ids, [4, 1]), node_ids
        elem.write(size, 'dummy')

        save_load_deck(model)


    def test_cfast_01(self):
        """tests a CFAST/PFASTE"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        size = 8
        eid = 2001
        nids = [4, 1]
        pid = 20
        ce = 1.0
        cr = 2.0

        fast_type = 'cat'
        ida = 10
        idb = 9

        model.add_grid(4, [0., 0., 0.])
        model.add_grid(1, [1., 0., 0.])

        #eid: int, pid: int, Type: str,
        #ida: int, idb: int
        model.add_cfast(eid, pid, fast_type, ida, idb, comment='cvisc')
        d = 42
        kt1 = 1e7
        kt2 = 2e7
        kt3 = 3e7
        model.add_pfast(pid, d, kt1, kt2, kt3, comment='pvisc')
        model.setup()

        elem = model.cfast
        self.assertEqual(elem.element_id, 2001)
        self.assertEqual(elem.property_id, 20)
        node_ids = elem.nodes[0]
        assert np.array_equal(node_ids, [0, 0]), node_ids
        elem.write(size, 'dummy')

        save_load_deck(model)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
