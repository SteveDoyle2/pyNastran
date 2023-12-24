"""tests dynamic cards and dynamic load cards"""
import unittest
from io import StringIO
import numpy as np

import pyNastran
from pyNastran.dev.bdf_vectorized3.bdf import BDF, read_bdf, CrossReferenceError
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck

#ROOT_PATH = pyNastran.__path__[0]
TEST_DYNAMIC = False

class TestDynamic(unittest.TestCase):
    """
    The cards tested are:
     * TSTEP
    """
    def test_tstep(self):
        """tests a TSTEP card"""
        model = BDF(debug=None)

        sid = 42
        n1 = n2 = 5
        dt1 = dt2 = 0.1
        no1 = no2 = 3
        card = ['TSTEP', sid,
                n1, dt1, no1, None, None, None, None, None,
                n2, dt2, no2]
        model.add_card(card, card[0], comment='tstep comment')
        model.validate()
        tstep = model.tsteps[42]
        tstep.raw_fields()
        tstep.write_card()
        tstep.write_card(size=16)

        sid = 43
        N = 5
        DT = 0.1
        NO = 3
        tstep2 = model.add_tstep(sid, N, DT, NO)
        tstep2.raw_fields()
        tstep2.write_card()
        tstep2.write_card(size=16)
        save_load_deck(model)

    def test_tstepnl(self):
        """tests a TSTEPNL card"""
        model = BDF(debug=None)
        card = ['TSTEPNL', 250, 100, .01, 1, 'ADAPT', 2, 10, 'PW',
                1.E-2, 1.E-3, 1.E-6, 2, 10, 2, .02, None,
                5, 5, 0, 0.75, 16.0, 0.1, 20.,]
        model.add_card(card, card[0], comment='tstepnl comment')
        model.validate()
        tstepnl = model.tstepnls[250]
        tstepnl.raw_fields()
        tstepnl.write_card()
        tstepnl.write_card(size=16)

        sid = 42
        ndt = 10
        dt = 3.
        no = 5
        tstepnl2 = model.add_tstepnl(sid, ndt, dt, no)
        tstepnl2.raw_fields()
        tstepnl2.write_card()
        tstepnl2.write_card(size=16)
        save_load_deck(model)

    def test_delay(self):
        """tests a two field DELAY card"""
        model = BDF(debug=False)
        node1, c1, t1 = 100, 3, 0.3
        node2, c2, t2 = 101, 4, 0.4
        sid = 42
        card_lines = ['DELAY', sid, node1, c1, t1, node2, c2, t2]
        model.add_card(card_lines, card_lines[0], comment='', is_list=True,
                       has_none=True)
        model.add_grid(100, [0., 0., 0.])
        model.add_grid(101, [0., 0., 0.])
        model.validate()
        model.cross_reference()
        #print(model.delays[42])
        save_load_deck(model)

    def test_dphase(self):
        """tests a two field DPHASE card"""
        model = BDF(debug=False)
        node1, c1, t1 = 100, 3, 0.3
        node2, c2, t2 = 101, 4, 0.4
        sid = 42
        card_lines = ['DPHASE', sid, node1, c1, t1, node2, c2, t2]
        model.add_card(card_lines, card_lines[0], comment='', is_list=True,
                       has_none=True)
        model.add_grid(100, [0., 0., 0.])
        model.add_grid(101, [0., 0., 0.])
        model.validate()
        model.cross_reference()
        #print(model.dphases[42])
        save_load_deck(model)

    def test_freq(self):
        """tests FREQ, FREQ1, FREQ2, FREQ4"""
        model = BDF(debug=False)
        sid = 101
        freqs = 0.1
        freq = model.add_freq(sid, freqs, comment='freq')
        #print(freq)

        freqs = [2.0, 3.0]
        freq = model.add_freq(sid, freqs, comment='freq')
        #print(freq)

        f1 = 0.
        df = 2.0
        freq1 = model.add_freq1(sid, f1, df, ndf=5, comment='freq1')
        assert len(freq1.freqs) == 6, 'freqs=%s' % freq1.freqs
        #print(freq1)

        f1 = 1.
        f2 = 8.0
        freq2 = model.add_freq2(sid, f1, f2, nf=6, comment='freq2')
        assert len(freq2.freqs) == 7, 'freqs=%s' % freq2.freqs
        assert np.allclose(freq2.freqs.max(), f2), freq2.freqs
        #print(freq2)

        freq4 = model.add_freq4(sid, f1, f2, fspread=0.1, nfm=3, comment='freq4')
        #print(model.frequencies[sid])
        #print(freq4)

        fractions = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        freq5 = model.add_freq5(sid, fractions, f1=0., f2=100., comment='freq5')

        fractions = np.linspace(0., 1.)
        unused_freq5b = model.add_freq5(sid, fractions, f1=0., f2=100., comment='freq5')
        model.validate()

        freq.raw_fields()
        freq.write_card()
        freq.write_card(size=16)

        freq1.raw_fields()
        freq1.write_card()
        freq1.write_card(size=16)

        freq2.raw_fields()
        freq2.write_card()
        freq2.write_card(size=16)

        freq4.raw_fields()
        freq4.write_card()
        freq4.write_card(size=16)

        freq5.raw_fields()
        freq5.write_card()
        freq5.write_card(size=16)

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        unused_out = bdf_file.getvalue()
        bdf_file.seek(0)

        model2 = read_bdf(bdf_file, punch=True, debug=False)
        #model2.uncross_reference()
        #model2.safe_cross_reference()
        #model2.uncross_reference()
        save_load_deck(model)

    def test_tload(self):
        """tests DLOAD, TLOAD1, TLOAD2, TABLED2 cards"""
        model = BDF(debug=False)
        tload1 = model.tload1
        tload2 = model.tload2
        dload = model.dload
        delay = model.delay

        model.set_error_storage(nparse_errors=0, stop_on_parsing_error=True,
                                nxref_errors=0, stop_on_xref_error=True)

        sid = 2
        excite_id = 20
        #delay = 0
        tid = 42
        model.add_tload1(sid, excite_id, tid, delay=0, load_type='LOAD',
                         us0=0.0, vs0=0.0, comment='tload1')
        model.add_tload1(sid, excite_id, tid, delay=1., load_type='DISP',
                         us0=0.0, vs0=0.0, comment='')
        model.add_tload1(sid, excite_id, tid, delay=2, load_type='VELO',
                         us0=0.0, vs0=0.0, comment='')
        model.add_tload1(sid, excite_id, tid, delay=0, load_type='ACC',
                         us0=0.0, vs0=0.0, comment='')

        nid = 100
        model.add_grid(nid, [0., 0., 0.])

        darea_id = excite_id
        component = 4
        scale = 1.
        model.add_darea(darea_id, nid, component, scale, comment='')

        sid = 3
        excite_id = 30
        model.add_tload2(sid, excite_id, delay=0, load_type='LOAD',
                         T1=0., T2=None, frequency=0., phase=0.,
                         c=0., b=0., us0=0., vs0=0., comment='tload2')
        model.add_tload2(sid, excite_id, delay=1., load_type='D',
                         T1=0., T2=None, frequency=0., phase=0.,
                         c=0., b=0., us0=0., vs0=0., comment='')
        model.add_tload2(sid, excite_id, delay=2, load_type='V',
                         T1=0., T2=None, frequency=0., phase=0.,
                         c=0., b=0., us0=0., vs0=0., comment='')
        model.add_tload2(sid, excite_id, delay=0, load_type='A',
                         T1=0., T2=1., frequency=0., phase=0.,
                         c=0., b=0., us0=0., vs0=0., comment='')

        darea_id = excite_id
        component = 4
        scale = 1.
        model.add_darea(darea_id, nid, component, scale, comment='')

        delay_id = 2
        nodes = 100
        components = 2
        delays = 1.5

        sid = 1
        scale = 1.0
        scale_factors = 1.
        load_ids = 2
        dloadi = model.add_dload(sid, scale, scale_factors, load_ids,
                                 comment='dload')

        delayi = model.add_delay(delay_id, nodes, components, delays)
        x1 = 0.1
        x = np.linspace(0., 1.)
        y = np.sin(x)
        tabled2 = model.add_tabled2(tid, x1, x, y, comment='tabled2')

        model.pop_parse_errors()

        #delay.validate()
        #delay.raw_fields()
        delay.write()
        delay.write(size=16)

        #tload1.validate()
        #tload1.raw_fields()
        tload1.write()
        tload1.write(size=16)

        #tload2.validate()
        #tload2.raw_fields()
        tload2.write()
        tload2.write(size=16)

        dload.validate()
        #dload.raw_fields()
        dload.write()
        dload.write(size=16)

        if TEST_DYNAMIC:
            #tabled2.validate()
            #tabled2.raw_fields()
            tabled2.write()
            tabled2.write(size=16)

        model.validate()
        model.cross_reference()
        #model.pop_xref_errors()

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        unused_out = bdf_file.getvalue()
        bdf_file.seek(0)
        unused_outs = model.get_bdf_stats(return_type='list')
        unused_outs = model.get_bdf_stats(return_type='string')

        if TEST_DYNAMIC:
            time = 0.5
            out1 = tload1.get_load_at_time(time, scale=1.)
            out2 = tload2.get_load_at_time(time, scale=1.)
            #print(out1)
            assert len(out1) == 1, out1
            assert len(out2) == 1, out2
            #print(out1)
            #print(out2)

            time = [0.5, 0.9]
            out1 = tload1.get_load_at_time(time, scale=1.)
            out2 = tload2.get_load_at_time(time, scale=1.)
            assert len(out1) == 2, out1
            assert len(out2) == 2, out2
            #print(out1)
            #print(out2)

        model2 = read_bdf(bdf_file, punch=True, debug=False)
        #model2.uncross_reference()
        #model2.safe_cross_reference()
        #model2.uncross_reference()
        save_load_deck(model, run_renumber=False, run_convert=False)

    def test_rload1(self):
        """tests DLOAD, RLOAD1, TABLED2 cards"""
        model = BDF(debug=False)
        rload1 = model.rload1
        dload = model.dload
        delay = model.delay
        #model.case_control_deck = CaseControlDeck(['DLOAD=2', 'BEGIN BULK'])
        sid = 2
        excite_id = 20
        delayi = 0
        tid = 42
        rload1i = model.add_rload1(sid, excite_id, delay=0, dphase=0, tc=0,
                                   td=0, load_type='LOAD', comment='rload1')
        rload1i = model.add_rload1(sid, excite_id, delay=1., dphase=0, tc=0,
                                   td=0, load_type='DISP', comment='rload1')
        rload1i = model.add_rload1(sid, excite_id, delay=2, dphase=0, tc=0,
                                   td=0, load_type='VELO', comment='rload1')
        rload1i = model.add_rload1(sid, excite_id, delay=0, dphase=0, tc=0,
                                   td=0, load_type='ACC', comment='rload1')

        excite_id = 20
        nid = 21
        c = 1
        scale = 1.0
        model.add_darea(excite_id, nid, c, scale, comment='darea')
        model.add_grid(nid, [0., 0., 0.])


        sid = 1
        scale = 1.0
        scale_factors = 1.
        load_ids = 2
        dloadi = model.add_dload(sid, scale, scale_factors, load_ids,
                                 comment='dload')

        delay_id = 2
        nodes = 100
        components = 2
        delays = 1.5
        delayi = model.add_delay(delay_id, nodes, components, delays)
        model.add_grid(nodes, [0., 0., 0.])

        x1 = 0.1
        x = np.linspace(0., 1.)
        y = np.sin(x)
        tabled2i = model.add_tabled2(tid, x1, x, y, comment='tabled2')

        model.setup()
        model.pop_parse_errors()

        rload1.validate()
        #rload1.raw_fields()
        rload1.write()
        rload1.write(size=16)

        dload.validate()
        #dload.raw_fields()
        dload.write()
        dload.write(size=16)

        #delay.validate()
        #delay.raw_fields()
        delay.write()
        delay.write(size=16)
        if TEST_DYNAMIC:
            tabled2.validate()
            tabled2.raw_fields()
            tabled2.write_card()
            tabled2.write_card(size=16)

        model.validate()
        model.cross_reference()
        #model.pop_xref_errors()
        #print(model.dareas)

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        unused_out = bdf_file.getvalue()
        bdf_file.seek(0)
        unused_outs = model.get_bdf_stats(return_type='list')
        unused_outs = model.get_bdf_stats(return_type='string')

        if TEST_DYNAMIC:
            freq = 0.5
            out1 = rload1.get_load_at_freq(freq, scale=1.)
            #out2 = rload2.get_load_at_time(freq, scale=1.)
            #print(out1)
            #print(out2)
            assert len(out1) == 1, out1
            #assert len(out2) == 1, out2

            freq = [0.5, 0.9]
            out1 = rload1.get_load_at_freq(freq, scale=1.)
            #out2 = rload2.get_load_at_freq(freq, scale=1.)
            #print(out1)
            #print(out2)
            assert len(out1) == 2, out1
            #assert len(out2) == 2, out2

        model2 = read_bdf(bdf_file, punch=True, debug=False)
        #model2.uncross_reference()
        #model2.safe_cross_reference()
        #model2.uncross_reference()
        #print(out)
        #print(outs)
        save_load_deck(model, run_renumber=False, run_convert=False)

    def test_rload2(self):
        """tests DLOAD, RLOAD2, TABLED2 cards"""
        model = BDF(debug=False)
        rload2 = model.rload2
        dload = model.dload
        delay = model.delay
        #model.case_control_deck = CaseControlDeck(['DLOAD=2', 'BEGIN BULK'])

        sid = 3
        excite_id = 30
        rload2i = model.add_rload2(sid, excite_id, delay=0, dphase=0, tb=0,
                                   tphi=0, load_type='LOAD', comment='rload2')
        rload2i = model.add_rload2(sid, excite_id, delay=1., dphase=0, tb=0,
                                   tphi=0, load_type='D', comment='rload2')
        rload2i = model.add_rload2(sid, excite_id, delay=2, dphase=0, tb=0,
                                   tphi=0, load_type='V', comment='rload2')
        rload2i = model.add_rload2(sid, excite_id, delay=0, dphase=0, tb=0,
                                   tphi=0, load_type='A', comment='rload2')

        #excite_id = 20
        nid = 21
        c = 1
        scale = 1.0
        #model.add_darea(excite_id, nid, c, scale, comment='darea')
        model.add_grid(nid, [0., 0., 0.])

        excite_id = 30
        model.add_darea(excite_id, nid, c, scale, comment='darea')

        sid = 1
        scale = 1.0
        scale_factors = 1.
        load_ids = 3
        dloadi = model.add_dload(sid, scale, scale_factors, load_ids,
                                 comment='dload')

        delay_id = 2
        nodes = 100
        components = 2
        delays = 1.5
        delayi = model.add_delay(delay_id, nodes, components, delays)
        model.add_grid(nodes, [0., 0., 0.])

        model.pop_parse_errors()

        rload2.validate()
        #rload2.raw_fields()
        rload2.write()
        rload2.write(size=16)

        dload.validate()
        #dload.raw_fields()
        dload.write()
        dload.write(size=16)

        #delay.validate()
        #delay.raw_fields()
        delay.write()
        delay.write(size=16)

        model.validate()
        model.cross_reference()
        #model.pop_xref_errors()
        #print(model.dareas)

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        unused_out = bdf_file.getvalue()
        bdf_file.seek(0)
        unused_outs = model.get_bdf_stats(return_type='list')
        unused_outs = model.get_bdf_stats(return_type='string')

        if TEST_DYNAMIC:
            freq = 0.5
            out1 = rload2.get_load_at_freq(freq, scale=1.)
            #out2 = rload2.get_load_at_time(freq, scale=1.)
            #print(out1)
            #print(out2)
            assert len(out1) == 1, out1
            #assert len(out2) == 1, out2
            assert np.allclose(out1, np.zeros(1, dtype='complex128'))

            freq = [0.5, 0.9]
            out1 = rload2.get_load_at_freq(freq, scale=1.)
            #out2 = rload2.get_load_at_freq(freq, scale=1.)
            #print(out1)
            #print(out2)
            assert len(out1) == 2, out1
            assert np.allclose(out1, np.zeros(2, dtype='complex128'))
            #assert len(out2) == 2, out2

        model2 = read_bdf(bdf_file, punch=True, debug=False)
        #model2.uncross_reference()
        #model2.safe_cross_reference()
        #model2.uncross_reference()
        #print(out)
        #print(outs)
        save_load_deck(model, run_renumber=False, run_convert=False)

    def test_ascre(self):
        """tests ASCRE, DELAY, DPHASE, TABLED2"""
        model = BDF(debug=False)
        sid = 1
        excite_id = 2
        rho = 1.0
        b = 2.0
        acsrce_id1 = model.add_acsrce(sid, excite_id, rho, b, delay=0, dphase=0, power=0,
                                      comment='acsrce')
        #acsrce.raw_fields()
        delay = 3
        dphase = 4
        power = 5

        sid = 3
        excite_id = 4
        rho = 1.0
        b = 2.0
        unused_acsrce2 = model.add_acsrce(sid, excite_id, rho, b, delay=delay,
                                          dphase=dphase, power=power)

        nodes = 4
        components = 5
        delays = 6.0
        delay = model.add_delay(sid, nodes, components, delays, comment='')

        nodes = 4
        components = 6
        phase_leads = 2.0
        delay = model.add_dphase(sid, nodes, components, phase_leads)

        tid = power
        x1 = 1.
        x = np.linspace(0., 1.) + 10.
        y = np.sin(x) + 2.
        model.add_tabled2(tid, x1, x, y, comment='tabled2')
        model.add_grid(4, [0., 0., 0.])

        model.validate()
        model.pop_parse_errors()
        model.cross_reference()
        #model.pop_xref_errors()

        save_load_deck(model, run_convert=False)

    def test_nlparm(self):
        """tests NLPARM"""
        model = BDF(debug=False)
        nlparm_id = 42
        model.add_nlparm(nlparm_id, comment='nlparm')
        save_load_deck(model)

    def test_nlpci(self):
        """tests NLPCI"""
        model = BDF(debug=False)
        nlpci_id = 42
        nlpci = model.add_nlpci(nlpci_id, Type='CRIS', minalr=0.25, maxalr=4.,
                                scale=0., desiter=12, mxinc=20,
                                comment='nlpci')
        nlpci.raw_fields()
        #print(nlpci)
        save_load_deck(model)

    def _test_rotorg(self):
        """tests ROTORD"""
        model = BDF(debug=False)
        rotorg = model.rotorg
        card_lines = ['ROTORG', '10', '1', '2']
        card_lines = ['ROTORG', '10', '10', ' THRU', '20']
        model.add_card(card_lines, 'ROTORG', comment='rotorg')
        model.setup()
        print(rotorg.write())

    def _test_rotord(self):
        """tests ROTORD"""
        model = BDF(debug=False)
        sid = 42
        rstart = 10.0
        rstep = 11.0
        numstep = 10
        rids = []
        rsets = [31]
        rspeeds = [None]
        rcords = []
        w3s = []
        w4s = []
        rforces = []
        brgsets = []
        rotord = model.add_rotord(
            sid, rstart, rstep, numstep,
            rids, rsets, rspeeds, rcords, w3s, w4s, rforces, brgsets,
            refsys='ROT', cmout=0.0, runit='RPM',
            funit='RPM', zstein='NO', orbeps=1.e-6,
            roprt=0, sync=1, etype=1, eorder=1.0,
            threshold=0.02, maxiter=10, comment='rotord')
        rotord.validate()
        save_load_deck(model)

    def _test_loadcyn(self):
        """tests LOADCYN"""
        model = BDF(debug=False, log=None, mode='msc')
        sid = 42
        scale = 4.
        segment_id = 10
        scales = [1.]
        load_ids = [3]
        loadcyn = model.add_loadcyn(sid, scale, segment_id, scales, load_ids,
                                    segment_type=None, comment='loadcyn')
        loadcyn.validate()
        model.pop_parse_errors()
        card = loadcyn.write_card(size=8)
        loadcyn.write_card(size=16, is_double=False)
        loadcyn.write_card(size=16, is_double=True)
        loadcyn.raw_fields()
        str(loadcyn)
        #print(model.loads)
        model.loads = {}
        model.add_card(card.split('\n')[1:], 'LOADCYN', comment='', is_list=False, has_none=True)

        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(model, run_convert=False)

    def test_deform(self):
        """tests DEFORM"""
        model = BDF(debug=False, log=None, mode='msc')
        tload1 = model.tload1

        sid = 42
        eid = 10
        deformation = 32.
        deform = model.deform
        model.add_deform(sid, eid, deformation, comment='deform')
        deform.validate()
        model.pop_parse_errors()
        card = deform.write(size=8)
        deform.write(size=16, is_double=False)
        deform.write(size=16, is_double=True)
        with self.assertRaises(ValueError):
            model.setup()
        #deform.raw_fields()
        str(deform)
        #model.loads = {}
        #model.add_card(card.split('\n')[1:], 'DEFORM', comment='', is_list=False, has_none=True)
        #model.pop_parse_errors()

        #with self.assertRaises(CrossReferenceError):
            #model.cross_reference()
        #with self.assertRaises(CrossReferenceError):
            #model.pop_xref_errors()
        #model.uncross_reference()
        #model.reset_errors()
        #model.safe_cross_reference()

        delta = 0.1
        eid1 = 11
        eid2 = 12
        eid3 = 13
        fields = ['DEFORM', sid, eid1, delta, eid2, delta, eid3, delta]
        model.add_card(fields, 'DEFORM')

        eid = 10
        nids = [2, 3]
        mid = 100
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.add_conrod(eid, mid, nids, A=0.0, j=0.0, c=0.0, nsm=0.0, comment='')
        model.add_conrod(eid1, mid, nids, A=0.0, j=0.0, c=0.0, nsm=0.0, comment='')
        model.add_conrod(eid2, mid, nids, A=0.0, j=0.0, c=0.0, nsm=0.0, comment='')
        model.add_conrod(eid3, mid, nids, A=0.0, j=0.0, c=0.0, nsm=0.0, comment='')
        model.setup()
        save_load_deck(model)

    def test_tload1(self):
        model = BDF(debug=False, log=None, mode='msc')
        tload1 = model.tload1
        excite_id = 11
        node = 1
        mag = 1.
        xyz = [1., 0., 0.]
        model.add_grid(1, [0., 0., 0.])
        model.add_force(excite_id, node, mag, xyz,
                        cid=0, comment='')

        tabled2_id = 100
        x1 = 0.
        x = t = np.linspace(0., 1., num=11)
        y = np.sin(t)
        model.add_tabled2(tabled2_id, x1, x, y, extrap=0, comment='')

        tload_id = 12
        tload1a = model.add_tload1(tload_id, excite_id, tabled2_id, delay=0,
                                  load_type='LOAD', us0=0.0, vs0=0.0, comment='')

        tload_id = 13
        tload1b = model.add_tload1(tload_id, excite_id, tabled2_id, delay=.1,
                                  load_type='LOAD', us0=0.0, vs0=0.0, comment='')
        model.setup()
        #model.cross_reference()

        if TEST_DYNAMIC:
            time = np.linspace(0., 1., num=11)
            loads = tload1a.get_load_at_time(time, scale=1.)
            loads_expected = y
            assert np.allclose(loads, loads_expected)

            loads = tload1b.get_load_at_time(time, scale=1.)
            loads_expected = np.hstack([0., y[:-1]])
            assert np.allclose(loads, loads_expected)

        print(tload1.write())
        save_load_deck(model, run_renumber=False)

    def test_tload2(self):
        model = BDF(debug=False, log=None, mode='msc')
        tload2 = model.tload2
        excite_id = 11
        node = 1
        mag = 1.
        xyz = [1., 0., 0.]
        model.add_grid(1, [0., 0., 0.])
        model.add_force(excite_id, node, mag, xyz,
                        cid=0, comment='')
        tload_id = 12
        model.add_tload2(
            tload_id, excite_id, delay=0, load_type='LOAD', T1=0., T2=None,
            frequency=0., phase=0., c=0., b=0., us0=0., vs0=0., comment='')
        num = 11
        if TEST_DYNAMIC:
            time = np.linspace(0., 1., num=num, dtype='float64')
            loads = tload2.get_load_at_time(time, scale=1.)
            loads_expected = np.ones(num)
            assert np.array_equal(loads, loads_expected)
            assert np.array_equal(loads, [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])

        tload_id = 13
        model.add_tload2(
            tload_id, excite_id, delay=0, load_type='LOAD', T1=0., T2=0.5,
            frequency=0., phase=0., c=0., b=0., us0=0., vs0=0., comment='')
        if TEST_DYNAMIC:
            loads = tload2.get_load_at_time(time, scale=1.)
            assert np.allclose(loads, [1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0.])

        tload_id = 14
        model.add_tload2(
            tload_id, excite_id, delay=0, load_type='LOAD', T1=0.15, T2=0.5,
            frequency=0., phase=0., c=0., b=0., us0=0., vs0=0., comment='')
        if TEST_DYNAMIC:
            loads = tload2.get_load_at_time(time, scale=1.)
            assert np.allclose(loads, [0., 0., 1., 1., 1., 1., 0., 0., 0., 0., 0.])

        tload_id = 15
        model.add_tload2(
            tload_id, excite_id, delay=0, load_type='LOAD', T1=0.15, T2=0.5,
            frequency=2., phase=0., c=0., b=0., us0=0., vs0=0., comment='')
        if TEST_DYNAMIC:
            loads = tload2.get_load_at_time(time, scale=1.)
            val = -0.8090169943749472
            val2 = 0.30901699437494723
            assert np.allclose(loads, [0., 0., val, val, val2, 1., 0., 0., 0., 0., 0.])

        tload_id = 16
        model.add_tload2(
            tload_id, excite_id, delay=0, load_type='LOAD', T1=0.15, T2=0.5,
            frequency=10., phase=45., c=0., b=0., us0=0., vs0=0., comment='')
        if TEST_DYNAMIC:
            loads = tload2.get_load_at_time(time, scale=1.)
            val = 0.5 ** 0.5 # sqrt(1/2)
            assert np.allclose(loads, [0., 0., val, val, val, val, 0., 0., 0., 0., 0.])

        tload_id = 17
        model.add_tload2(
            tload_id, excite_id, delay=0, load_type='LOAD', T1=0.15, T2=0.5,
            frequency=10., phase=45., c=3., b=5., us0=0., vs0=0., comment='')
        if TEST_DYNAMIC:
            loads = tload2.get_load_at_time(time, scale=1.)
            #print(loads)
            loads_expected = np.array([0.0, 0.0, 0.0004122984191627603, 0.004226260954591647, 0.02404021442968808, 0.09903227290017788, 0.0, 0.0, 0.0, 0.0, 0.0], dtype='float64')
            assert np.allclose(loads, loads_expected)

        #------------------------------------------------
        # phase int
        # delay int
        tload_id = 26
        delay_id = 15
        dphase_id = 45
        model.add_tload2(
            tload_id, excite_id, delay=delay_id, load_type='LOAD', T1=0.15, T2=0.5,
            frequency=10., phase=dphase_id, c=0., b=0., us0=0., vs0=0., comment='')

        delay = 0.15
        phase = 45.
        nodes = [1]
        components = [1]
        delays = [delay]
        model.add_delay(delay_id, nodes, components, delays, comment='')
        nodes = [1]
        components = [1]
        phase_leads = [phase]
        model.add_dphase(dphase_id, nodes, components, phase_leads, comment='')

        nodes = [1]
        components = [2]
        delays = [delay]
        model.add_delay(delay_id, nodes, components, delays, comment='')
        nodes = [1]
        components = [2]
        phase_leads = [phase]
        model.add_dphase(dphase_id, nodes, components, phase_leads, comment='')

        nodes = [1]
        components = [3]
        delays = [delay]
        model.add_delay(delay_id, nodes, components, delays, comment='')
        nodes = [1]
        components = [3]
        phase_leads = [phase]
        model.add_dphase(dphase_id, nodes, components, phase_leads, comment='')

        if TEST_DYNAMIC:
            model.cross_reference()
            loads = tload2.get_load_at_time(time, scale=1.)
            val = 0.5 ** 0.5 # sqrt(1/2)
            with self.assertRaises(AssertionError):
                assert np.allclose(loads, [0., 0., val, val, val, val, 0., 0., 0., 0., 0.])

        model.cross_reference()
        print(tload2.write())
        #save_load_deck(model)

    def test_rforce(self):
        """tests RFORCE"""
        model = BDF(debug=False, log=None, mode='msc')
        rforce = model.rforce
        #model._nxref_errors = 0
        sid = 42
        nid = 2
        cid = 1
        scale = 2.
        r123 = [0., 1., 2.]
        model.add_rforce(sid, nid, scale, r123, cid=cid,
                         method=1, racc=0., main_bulk=0, idrf=0, comment='rforce')
        #rforce.validate()
        model.setup()
        card = rforce.write(size=8)
        rforce.write(size=16, is_double=False)
        rforce.write(size=16, is_double=True)
        #rforce.raw_fields()
        str(rforce)
        #model.loads = {}
        #model.add_card(card.split('\n')[1:], 'RFORCE', comment='', is_list=False, has_none=True)
        #model.pop_parse_errors()

        #with self.assertRaises(CrossReferenceError):
            #model.cross_reference()
        #with self.assertRaises(CrossReferenceError):
            #model.pop_xref_errors()
        #model.uncross_reference()
        #model.reset_errors()
        #with self.assertRaises(KeyError):
            #model.safe_cross_reference()
        #model.reset_errors()

        model.add_grid(2, [0., 0., 0.])
        model.add_cord2r(cid, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.], rid=0, comment='')
        model.cross_reference()
        save_load_deck(model, run_convert=False)

    def test_rforce1(self):
        """tests RFORCE1"""
        model = BDF(debug=False, log=None, mode='msc')
        rforce1 = model.rforce1
        sid = 42
        nid = 2
        scale = 2.
        #r123 = None
        group_id = -4
        cid = 1
        rforce1_id1 = model.add_rforce1(sid, nid, scale, group_id, cid=cid, r123=None,
                                        racc=0., main_bulk=0, method=2, comment='rforce1')
        rforce1.validate()
        rforce1_id2 = model.add_rforce1(sid, nid, scale, group_id, cid=0, r123=[1., 2., 3.],
                                        racc=0., main_bulk=0, method=2, comment='rforce1')
        #rforce1b.validate()
        model.pop_parse_errors()
        model.setup()

        card1 = rforce1.slice_card_by_index([rforce1_id1]).write(size=8)
        card2 = rforce1.slice_card_by_index([rforce1_id2]).write(size=8)
        rforce1.write(size=16, is_double=False)
        rforce1.write(size=16, is_double=True)
        #rforce1.raw_fields()
        str(rforce1)
        #model.loads = {}
        card1i = card1.split('\n')
        model.add_card(card1i, 'RFORCE1', comment='', is_list=False, has_none=True)
        model.pop_parse_errors()

        #with self.assertRaises(CrossReferenceError):
            #model.cross_reference()
        #with self.assertRaises(CrossReferenceError):
            #model.pop_xref_errors()
        #model.uncross_reference()
        #model.reset_errors()
        #with self.assertRaises(KeyError):
            #model.safe_cross_reference()
        #model.reset_errors()

        model.add_grid(2, [0., 0., 0.])
        model.add_cord2r(cid, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.], rid=0, comment='')
        #model.cross_reference()
        model.setup(run_geom_check=True)
        save_load_deck(model, run_convert=False)

    def test_transfer_function(self):
        model = BDF(debug=False, log=None, mode='msc')
        tf = model.tf

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [2., 0., 0.])
        tf_id = 10
        nid0 = 1
        component = 3
        b0 = 1.
        b1 = 2.
        b2 = 3.
        nids = [2, 3]
        components = [4, 5]
        a = [
            (10., 11., 12.),
            (13., 14., 15.),
        ]
        model.add_tf(tf_id, nid0, component, b0, b1, b2, nids, components, a)
        model.add_tf(tf_id+1, nid0, component, b0, b1, b2, nids, components, a)
        model.add_tf(tf_id+2, nid0, component, b0, b1, b2, nids, components, a)
        model.setup()

        tf.write(size=8)
        tf.write(size=16)
        save_load_deck(model)

    def test_dynamic1(self):
        """
        xref test for:
         - DLOAD -> DAREA -> NID

        DLOAD take priority
        useful for dynamic nodal forces/disp/vel/acc
        """
        msg = """
SOL 108
CEND
SUBCASE 1
    DLOAD = 33
    DISP(PLOT) = ALL
BEGIN BULK
$DLOAD SID S    S1   L1   S2  L2
DLOAD, 33, 1.0, 1.0, 35, 1.0, 36
$RLOAD1 SID EXCITEID DELAY DPHASE TC   TD   TYPE
RLOAD1, 35, 29,      0.2,  5.0,   40,  0.0, 0
RLOAD1, 36, 29,      31,   32,    4.0, 41,  0
$DAREA SID GRID COMP SCALE
DAREA, 29, 30,  1,   5.2
$DELAY SID GRID COMP LAG
DELAY, 31, 30,  1,   0.2
$DPHASE SID GRID COMP ANGLE
DPHASE, 32, 30,  1,   5.0
$TABLED1 TID XAXIS YAXIS
$ x1 y1 x2 y2 x3 y3 x4 y4
TABLED1, 40, LINEAR, LINEAR
,0.0, 4.0, 2.0, 8.0, 6.0, 8.0, ENDT
TABLED1, 41, LINEAR, LINEAR
,0.0, 0.5, 0.6, 0.4, 0.8, 0.7, ENDT
GRID,30
"""
        model = BDF(debug=False)
        bdf_file = StringIO()
        bdf_file.write(msg)
        bdf_file.seek(0)
        model.read_bdf(bdf_file)
        #In the example:
        # * The DLOAD case control command selects the loading reference
        #   by the DLOAD bulk entry having SID = 33 as the dynamic
        #   loading for the analysis.
        # * The DLOAD bulk entry combines the dynamic loads defined by
        #   two RLOAD1 entries having SIDs of 35 and 36. Neither dynamic
        #   load is scaled using the DLOAD entry.
        # * Both RLOAD1 entries reference the same DAREA entry. Thus,
        #   both dynamic loads are applied to the same degree-of-freedom.
        #   In this example, it is a single degree-of-freedom, Component 1
        #   of Grid 30. Both dynamic loads are scaled 5.2 times by the
        #   DAREA entry.
        # * Because the dynamic loads are applied at only one
        #   degree-of-freedom, the time delay and phase angle can be
        #   defined directly on the RLOAD1 entries. This is the case
        #   for the RLOAD1 entry having SID = 35. However, for
        #   demonstration purposes, the RLOAD1 entry having SID = 36
        #   references DELAY and DPHASE bulk entries. Both approaches
        #   define a delay of 0.2 and a phase angle of 5.0 for the
        #   corresponding dynamic load.
        # * C(f) for the RLOAD1 entry having SID = 35 is defined by the
        #   TABLED1 entry having TID = 40. (See Figure 6-6.) D(f) for
        #    this same RLOAD1 entry is defined as zero.
        # * C(f) for the RLOAD1 entry having SID = 36 is a constant
        #   value of 4.0. D(f) for this same RLOAD entry is defined by
        #   the TABLED1 entry having TID = 41.

    def test_dynamic2(self):
        """
        xref test for:
         - LOADSET -> LSEQ   -> FORCE, PLOAD
         - DLOAD   -> RLOAD1 -> TABLED1

        LOADSET take priority
        useful for generalized dynamic forces/disp/vel/acc
        """
        msg = """
SOL 108
CEND
SUBCASE 1
    LOADSET = 27
    DLOAD = 25
    DISP(PLOT) = ALL
BEGIN BULK
$LSEQ   SID EXCITEID LID
LSEQ,   27, 28,      26
$RLOAD1 SID EXCITEID DELAY DPHASE TC TD
RLOAD1, 25, 28,      0.0,  10.0,  29
$FORCE SID GRID CID F    N1 N2 N3
FORCE, 26, 425, ,   2.5, 1.0
$PLOAD SID PRES  GRID1 GRID2 GRID3 GRID4
PLOAD, 26, 50.0, 63,   64,   88,   91
$TABLED1 TID XAXIS YAXIS
$ x1 y1 x2 y2 x3 y3 x4 y4
TABLED1, 29, LINEAR, LINEAR
,0.0, 0.5, 0.6, 0.4, 0.8, 0.7, ENDT
"""
        model = BDF(debug=False)
        bdf_file = StringIO()
        bdf_file.write(msg)
        bdf_file.seek(0)
        model.read_bdf(bdf_file)
        #In the example:
        # * The LOADSET request in case control selects the LSEQ entry
        #   having SID = 27.
        # * The LSEQ entry references the static loads having SID = 26.
        #   These loads include the FORCE and PLOAD entries. The FORCE
        #   and PLOAD entries provide the spatial distribution of the
        #   dynamic loading.
        # * The DLOAD request in case control selects the RLOAD1 entry
        #   having SID = 25.
        # * The RLOAD1 entry references a TABLED1 entry having TID = 29.
        #   This TABLED1 entry defines C(f) for the RLOAD1 entry. Because
        #   the TD field on the RLOAD1 entry is undefined, D(f) defaults
        #   to zero.
        # * The EXCITEID fields of the LSEQ and RLOAD1 entries are both
        #   28, thereby linking the temporal and spatial distributions of
        #   the dynamic loading. Thus, the dynamic load defined by the
        #   RLOAD1 entry is:
        #   o Scaled by 2.5 and applied as a force to Component 1 of Grid 425.
        #   o Scaled by 50.0 and applied as a pressure to the quadrilateral
        #     element face defined by Grids 63, 64, 88, and 91.

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
