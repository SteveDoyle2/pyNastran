"""tests dynamic cards and dynamic load cards"""
import unittest
from io import StringIO
import numpy as np

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf, CrossReferenceError
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused

#ROOT_PATH = pyNastran.__path__[0]

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
        model2.uncross_reference()
        model2.safe_cross_reference()
        model2.uncross_reference()
        save_load_deck(model)

    def test_tload(self):
        """tests DLOAD, TLOAD1, TLOAD2, TABLED2 cards"""
        model = BDF(debug=False)
        model.set_error_storage(nparse_errors=0, stop_on_parsing_error=True,
                                nxref_errors=0, stop_on_xref_error=True)

        sid = 2
        excite_id = 20
        delay = 0
        tid = 42
        tload1 = model.add_tload1(sid, excite_id, tid, delay=0, Type='LOAD',
                                  us0=0.0, vs0=0.0, comment='tload1')
        tload1 = model.add_tload1(sid, excite_id, tid, delay=1., Type='DISP',
                                  us0=0.0, vs0=0.0, comment='')
        tload1 = model.add_tload1(sid, excite_id, tid, delay=2, Type='VELO',
                                  us0=0.0, vs0=0.0, comment='')
        tload1 = model.add_tload1(sid, excite_id, tid, delay=0, Type='ACC',
                                  us0=0.0, vs0=0.0, comment='')

        nid = 100
        model.add_grid(nid, [0., 0., 0.])

        darea_id = excite_id
        component = 4
        scale = 1.
        model.add_darea(darea_id, nid, component, scale, comment='')

        sid = 3
        excite_id = 30
        tload2 = model.add_tload2(sid, excite_id, delay=0, Type='LOAD',
                                  T1=0., T2=None, frequency=0., phase=0.,
                                  c=0., b=0., us0=0., vs0=0., comment='tload2')
        tload2 = model.add_tload2(sid, excite_id, delay=1., Type='D',
                                  T1=0., T2=None, frequency=0., phase=0.,
                                  c=0., b=0., us0=0., vs0=0., comment='')
        tload2 = model.add_tload2(sid, excite_id, delay=2, Type='V',
                                  T1=0., T2=None, frequency=0., phase=0.,
                                  c=0., b=0., us0=0., vs0=0., comment='')
        tload2 = model.add_tload2(sid, excite_id, delay=0, Type='A',
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
        delay = model.add_delay(delay_id, nodes, components, delays)

        sid = 1
        scale = 1.0
        scale_factors = 1.
        load_ids = 2
        dload = model.add_dload(sid, scale, scale_factors, load_ids,
                                comment='dload')

        x1 = 0.1
        x = np.linspace(0., 1.)
        y = np.sin(x)
        tabled2 = model.add_tabled2(tid, x1, x, y, comment='tabled2')

        model.pop_parse_errors()

        delay.validate()
        delay.raw_fields()
        delay.write_card()
        delay.write_card(size=16)

        tload1.validate()
        tload1.raw_fields()
        tload1.write_card()
        tload1.write_card(size=16)

        tload2.validate()
        tload2.raw_fields()
        tload2.write_card()
        tload2.write_card(size=16)

        dload.validate()
        dload.raw_fields()
        dload.write_card()
        dload.write_card(size=16)

        tabled2.validate()
        tabled2.raw_fields()
        tabled2.write_card()
        tabled2.write_card(size=16)

        model.validate()
        model.cross_reference()
        model.pop_xref_errors()

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        unused_out = bdf_file.getvalue()
        bdf_file.seek(0)
        unused_outs = model.get_bdf_stats(return_type='list')
        unused_outs = model.get_bdf_stats(return_type='string')

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
        model2.uncross_reference()
        model2.safe_cross_reference()
        model2.uncross_reference()
        save_load_deck(model, run_renumber=False, run_convert=False)

    def test_rload1(self):
        """tests DLOAD, RLOAD1, TABLED2 cards"""
        model = BDF(debug=False)
        #model.case_control_deck = CaseControlDeck(['DLOAD=2', 'BEGIN BULK'])
        sid = 2
        excite_id = 20
        delay = 0
        tid = 42
        rload1 = model.add_rload1(sid, excite_id, delay=0, dphase=0, tc=0,
                                  td=0, Type='LOAD', comment='rload1')
        rload1 = model.add_rload1(sid, excite_id, delay=1., dphase=0, tc=0,
                                  td=0, Type='DISP', comment='rload1')
        rload1 = model.add_rload1(sid, excite_id, delay=2, dphase=0, tc=0,
                                  td=0, Type='VELO', comment='rload1')
        rload1 = model.add_rload1(sid, excite_id, delay=0, dphase=0, tc=0,
                                  td=0, Type='ACC', comment='rload1')

        excite_id = 20
        nid = 21
        c = 1
        scale = 1.0
        model.add_darea(excite_id, nid, c, scale, comment='darea')
        model.add_grid(nid, [0., 0., 0.])

        delay_id = 2
        nodes = 100
        components = 2
        delays = 1.5
        delay = model.add_delay(delay_id, nodes, components, delays)

        sid = 1
        scale = 1.0
        scale_factors = 1.
        load_ids = 2
        dload = model.add_dload(sid, scale, scale_factors, load_ids,
                                comment='dload')

        x1 = 0.1
        x = np.linspace(0., 1.)
        y = np.sin(x)
        tabled2 = model.add_tabled2(tid, x1, x, y, comment='tabled2')

        model.pop_parse_errors()

        delay.validate()
        delay.raw_fields()
        delay.write_card()
        delay.write_card(size=16)

        rload1.validate()
        rload1.raw_fields()
        rload1.write_card()
        rload1.write_card(size=16)

        dload.validate()
        dload.raw_fields()
        dload.write_card()
        dload.write_card(size=16)

        tabled2.validate()
        tabled2.raw_fields()
        tabled2.write_card()
        tabled2.write_card(size=16)

        model.validate()
        model.cross_reference()
        model.pop_xref_errors()
        #print(model.dareas)

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        unused_out = bdf_file.getvalue()
        bdf_file.seek(0)
        unused_outs = model.get_bdf_stats(return_type='list')
        unused_outs = model.get_bdf_stats(return_type='string')

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
        model2.uncross_reference()
        model2.safe_cross_reference()
        model2.uncross_reference()
        #print(out)
        #print(outs)
        save_load_deck(model, run_renumber=False, run_convert=False)

    def test_rload2(self):
        """tests DLOAD, RLOAD2, TABLED2 cards"""
        model = BDF(debug=False)
        #model.case_control_deck = CaseControlDeck(['DLOAD=2', 'BEGIN BULK'])

        sid = 3
        excite_id = 30
        rload2 = model.add_rload2(sid, excite_id, delay=0, dphase=0, tb=0,
                                  tp=0, Type='LOAD', comment='rload2')
        rload2 = model.add_rload2(sid, excite_id, delay=1., dphase=0, tb=0,
                                  tp=0, Type='D', comment='rload2')
        rload2 = model.add_rload2(sid, excite_id, delay=2, dphase=0, tb=0,
                                  tp=0, Type='V', comment='rload2')
        rload2 = model.add_rload2(sid, excite_id, delay=0, dphase=0, tb=0,
                                  tp=0, Type='A', comment='rload2')

        #excite_id = 20
        nid = 21
        c = 1
        scale = 1.0
        #model.add_darea(excite_id, nid, c, scale, comment='darea')
        model.add_grid(nid, [0., 0., 0.])

        excite_id = 30
        model.add_darea(excite_id, nid, c, scale, comment='darea')

        delay_id = 2
        nodes = 100
        components = 2
        delays = 1.5
        delay = model.add_delay(delay_id, nodes, components, delays)

        sid = 1
        scale = 1.0
        scale_factors = 1.
        load_ids = 3
        dload = model.add_dload(sid, scale, scale_factors, load_ids,
                                comment='dload')

        model.pop_parse_errors()

        delay.validate()
        delay.raw_fields()
        delay.write_card()
        delay.write_card(size=16)

        rload2.validate()
        rload2.raw_fields()
        rload2.write_card()
        rload2.write_card(size=16)

        dload.validate()
        dload.raw_fields()
        dload.write_card()
        dload.write_card(size=16)

        model.validate()
        model.cross_reference()
        model.pop_xref_errors()
        #print(model.dareas)

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        unused_out = bdf_file.getvalue()
        bdf_file.seek(0)
        unused_outs = model.get_bdf_stats(return_type='list')
        unused_outs = model.get_bdf_stats(return_type='string')

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
        model2.uncross_reference()
        model2.safe_cross_reference()
        model2.uncross_reference()
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
        acsrce = model.add_acsrce(sid, excite_id, rho, b, delay=0, dphase=0, power=0,
                                  comment='acsrce')
        acsrce.raw_fields()
        sid = 3
        excite_id = 4
        rho = 1.0
        b = 2.0
        delay = 3
        dphase = 4
        power = 5
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
        model.pop_xref_errors()

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

    #def test_rotord(self):
        #"""tests ROTORD"""
        #model = BDF(debug=False)
        #sid = 42
        #rstart = 10.0
        #rstep = 11.0
        #numstep = 10
        #rids = []
        #rsets = [31]
        #rspeeds = [None]
        #rcords = []
        #w3s = []
        #w4s = []
        #rforces = []
        #brgsets = []
        #rotord = model.add_rotord(
            #sid, rstart, rstep, numstep,
            #rids, rsets, rspeeds, rcords, w3s, w4s, rforces, brgsets,
            #refsys='ROT', cmout=0.0, runit='RPM',
            #funit='RPM', zstein='NO', orbeps=1.e-6,
            #roprt=0, sync=1, etype=1, eorder=1.0,
            #threshold=0.02, maxiter=10, comment='rotord')
        #rotord.validate()
        #save_load_deck(model)

    def test_loadcyn(self):
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
        sid = 42
        eid = 10
        deformation = 32.
        deform = model.add_deform(sid, eid, deformation, comment='deform')
        deform.validate()
        model.pop_parse_errors()
        card = deform.write_card(size=8)
        deform.write_card(size=16, is_double=False)
        deform.write_card(size=16, is_double=True)
        deform.raw_fields()
        str(deform)
        model.loads = {}
        model.add_card(card.split('\n')[1:], 'DEFORM', comment='', is_list=False, has_none=True)
        model.pop_parse_errors()

        with self.assertRaises(CrossReferenceError):
            model.cross_reference()
        with self.assertRaises(CrossReferenceError):
            model.pop_xref_errors()
        model.uncross_reference()
        model.reset_errors()
        model.safe_cross_reference()

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
        model.cross_reference()
        save_load_deck(model)

    def test_tload1(self):
        model = BDF(debug=True, log=None, mode='msc')
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
                                  Type='LOAD', us0=0.0, vs0=0.0, comment='')

        tload_id = 13
        tload1b = model.add_tload1(tload_id, excite_id, tabled2_id, delay=.1,
                                  Type='LOAD', us0=0.0, vs0=0.0, comment='')
        model.cross_reference()

        time = np.linspace(0., 1., num=11)
        loads = tload1a.get_load_at_time(time, scale=1.)
        loads_expected = y
        assert np.allclose(loads, loads_expected)

        loads = tload1b.get_load_at_time(time, scale=1.)
        loads_expected = np.hstack([0., y[:-1]])
        assert np.allclose(loads, loads_expected)

        save_load_deck(model, run_renumber=False)

    def test_tload2(self):
        model = BDF(debug=True, log=None, mode='msc')
        excite_id = 11
        node = 1
        mag = 1.
        xyz = [1., 0., 0.]
        model.add_grid(1, [0., 0., 0.])
        model.add_force(excite_id, node, mag, xyz,
                        cid=0, comment='')
        tload_id = 12
        tload2 = model.add_tload2(
            tload_id, excite_id, delay=0, Type='LOAD', T1=0., T2=None,
            frequency=0., phase=0., c=0., b=0., us0=0., vs0=0., comment='')
        num = 11
        time = np.linspace(0., 1., num=num, dtype='float64')
        loads = tload2.get_load_at_time(time, scale=1.)
        loads_expected = np.ones(num)
        assert np.array_equal(loads, loads_expected)
        assert np.array_equal(loads, [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])

        tload_id = 13
        tload2 = model.add_tload2(
            tload_id, excite_id, delay=0, Type='LOAD', T1=0., T2=0.5,
            frequency=0., phase=0., c=0., b=0., us0=0., vs0=0., comment='')
        loads = tload2.get_load_at_time(time, scale=1.)
        assert np.allclose(loads, [1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0.])

        tload_id = 14
        tload2 = model.add_tload2(
            tload_id, excite_id, delay=0, Type='LOAD', T1=0.15, T2=0.5,
            frequency=0., phase=0., c=0., b=0., us0=0., vs0=0., comment='')
        loads = tload2.get_load_at_time(time, scale=1.)
        assert np.allclose(loads, [0., 0., 1., 1., 1., 1., 0., 0., 0., 0., 0.])

        tload_id = 15
        tload2 = model.add_tload2(
            tload_id, excite_id, delay=0, Type='LOAD', T1=0.15, T2=0.5,
            frequency=2., phase=0., c=0., b=0., us0=0., vs0=0., comment='')
        loads = tload2.get_load_at_time(time, scale=1.)
        val = -0.8090169943749472
        val2 = 0.30901699437494723
        assert np.allclose(loads, [0., 0., val, val, val2, 1., 0., 0., 0., 0., 0.])

        tload_id = 16
        tload2 = model.add_tload2(
            tload_id, excite_id, delay=0, Type='LOAD', T1=0.15, T2=0.5,
            frequency=10., phase=45., c=0., b=0., us0=0., vs0=0., comment='')
        loads = tload2.get_load_at_time(time, scale=1.)
        val = 0.5 ** 0.5 # sqrt(1/2)
        assert np.allclose(loads, [0., 0., val, val, val, val, 0., 0., 0., 0., 0.])

        tload_id = 17
        tload2 = model.add_tload2(
            tload_id, excite_id, delay=0, Type='LOAD', T1=0.15, T2=0.5,
            frequency=10., phase=45., c=3., b=5., us0=0., vs0=0., comment='')
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
        tload2 = model.add_tload2(
            tload_id, excite_id, delay=delay_id, Type='LOAD', T1=0.15, T2=0.5,
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

        model.cross_reference()
        loads = tload2.get_load_at_time(time, scale=1.)
        val = 0.5 ** 0.5 # sqrt(1/2)
        with self.assertRaises(AssertionError):
            assert np.allclose(loads, [0., 0., val, val, val, val, 0., 0., 0., 0., 0.])

        model.cross_reference()
        #save_load_deck(model)

    def test_rforce(self):
        """tests RFORCE"""
        model = BDF(debug=False, log=None, mode='msc')
        #model._nxref_errors = 0
        sid = 42
        nid = 2
        cid = 1
        scale = 2.
        r123 = [0., 1., 2.]
        rforce = model.add_rforce(sid, nid, scale, r123, cid=cid,
                                  method=1, racc=0., mb=0, idrf=0, comment='rforce')
        rforce.validate()
        card = rforce.write_card(size=8)
        rforce.write_card(size=16, is_double=False)
        rforce.write_card(size=16, is_double=True)
        rforce.raw_fields()
        str(rforce)
        model.loads = {}
        model.add_card(card.split('\n')[1:], 'RFORCE', comment='', is_list=False, has_none=True)
        model.pop_parse_errors()

        with self.assertRaises(CrossReferenceError):
            model.cross_reference()
        with self.assertRaises(CrossReferenceError):
            model.pop_xref_errors()
        model.uncross_reference()
        model.reset_errors()
        with self.assertRaises(KeyError):
            model.safe_cross_reference()
        model.reset_errors()

        model.add_grid(2, [0., 0., 0.])
        model.add_cord2r(cid, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.], rid=0, comment='')
        model.cross_reference()
        save_load_deck(model, run_convert=False)

    def test_rforce1(self):
        """tests RFORCE1"""
        model = BDF(debug=False, log=None, mode='msc')
        sid = 42
        nid = 2
        scale = 2.
        #r123 = None
        group_id = -4
        cid = 1
        rforce1 = model.add_rforce1(sid, nid, scale, group_id, cid=cid, r123=None,
                                    racc=0., mb=0, method=2, comment='rforce1')
        rforce1.validate()
        rforce1b = model.add_rforce1(sid, nid, scale, group_id, cid=0, r123=[1., 2., 3.],
                                     racc=0., mb=0, method=2, comment='rforce1')
        rforce1b.validate()
        model.pop_parse_errors()
        card = rforce1.write_card(size=8)
        rforce1.write_card(size=16, is_double=False)
        rforce1.write_card(size=16, is_double=True)
        rforce1.raw_fields()
        str(rforce1)
        model.loads = {}
        model.add_card(card.split('\n')[1:], 'RFORCE1', comment='', is_list=False, has_none=True)
        model.pop_parse_errors()

        with self.assertRaises(CrossReferenceError):
            model.cross_reference()
        with self.assertRaises(CrossReferenceError):
            model.pop_xref_errors()
        model.uncross_reference()
        model.reset_errors()
        with self.assertRaises(KeyError):
            model.safe_cross_reference()
        model.reset_errors()

        model.add_grid(2, [0., 0., 0.])
        model.add_cord2r(cid, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.], rid=0, comment='')
        model.cross_reference()
        save_load_deck(model, run_convert=False)

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

    def _test_dynamic2(self):
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

    def test_nolin(self) -> None:
        """
        Tests the NOLIN1/NOLIN2/NOLIN3/NOLIN4 cards.

        Tests:
         - parsing the NX QRG examples (blank components -> 0)
         - add_nolin1..4, grouping by SID in model.nolins
         - validate (component ranges; 10-16 velocity on CJ/CK only)
         - NOLIN1 TABLED1 cross-reference / uncross-reference
         - write (8/16) and write/read round trip (save_load_deck)
         - remove_unused keeps the NOLIN-referenced points
        Tolerances: exact (ids/components); floats are round-tripped
        through the 8-character writer, compared with assertAlmostEqual
        """
        model = BDF(debug=None)
        for nid in [2, 3, 4, 101, 102]:
            model.add_grid(nid, [float(nid), 0., 0.])
        model.add_tabled1(6, [0., 1.], [0., 2.])

        # NX QRG examples
        model.add_card(['NOLIN1', 21, 3, 4, 2.1, 3, 10, 6], 'NOLIN1', comment='nolin1')
        model.add_card(['NOLIN2', 14, 2, 1, 2.9, 2, 1, 2], 'NOLIN2')
        model.add_card(['NOLIN3', 4, 102, None, -6.1, 2, 15, -3.5], 'NOLIN3')
        model.add_card(['NOLIN4', 2, 4, 6, 2.0, 101, None, 16.3], 'NOLIN4')

        nolin2 = model.nolins[14][0]
        self.assertEqual(nolin2.type, 'NOLIN2')
        self.assertEqual(nolin2.ck, 0)  # blank
        nolin3 = model.nolins[4][0]
        self.assertEqual(nolin3.ci, 0)
        self.assertEqual(nolin3.cj, 15)
        self.assertAlmostEqual(nolin3.a, -3.5)
        nolin4 = model.nolins[2][0]
        self.assertEqual(nolin4.cj, 0)
        self.assertEqual(nolin4.node_ids, [4, 101])

        # same SID shares a list (NONLINEAR=21)
        nolin1b = model.add_nolin1(21, 2, 1, 1.0, 2, 11, 6)
        model.add_nolin2(21, 2, 2, 0.5, 3, 12, 4, 13)
        model.add_nolin3(21, 3, 3, 1.5, 4, 3, 2.0)
        model.add_nolin4(21, 3, 3, 1.5, 4, 3, 2.0)
        self.assertEqual(len(model.nolins[21]), 5)
        self.assertEqual(
            [card.type for card in model.nolins[21]],
            ['NOLIN1', 'NOLIN1', 'NOLIN2', 'NOLIN3', 'NOLIN4'])

        model.validate()
        for nolins in model.nolins.values():
            for nolin in nolins:
                nolin.raw_fields()
                nolin.write_card(size=8)
                nolin.write_card(size=16)

        # velocity components are not allowed on the load point (CI)
        bad = model.nolins[21][0].__class__(1, 2, 11, 1.0, 2, 1, 6)
        with self.assertRaises(AssertionError):
            bad.validate()

        model.cross_reference()
        self.assertEqual(nolin1b.tid_ref.tid, 6)
        self.assertEqual(nolin1b.Tid(), 6)
        model.uncross_reference()
        self.assertIsNone(nolin1b.tid_ref)
        model.safe_cross_reference()

        model2 = save_load_deck(model, run_remove_unused=False)
        self.assertEqual(sorted(model2.nolins), [2, 4, 14, 21])
        self.assertEqual(len(model2.nolins[21]), 5)
        nolin3b = model2.nolins[4][0]
        self.assertEqual(nolin3b.cj, 15)
        self.assertAlmostEqual(nolin3b.s, -6.1)

        remove_unused(model2)
        self.assertEqual(len(model2.nodes), 5)

    def test_rcross(self) -> None:
        """
        Tests the RCROSS card.

        Tests:
         - parsing the two NX QRG examples (incl. blank CURID)
         - add_rcross, CURID=0 (OP2 blank) -> None
         - validate rejects unknown RTYPE and both-blank RTYPEs
         - write (8/16) and write/read round trip (save_load_deck)
         - remove_unused keeps the RCROSS grids (DISP/VELO/PRESS and a
           blank RTYPE2 that inherits VELO) and removes an unused grid;
           STRESS ids are treated as elements, not grids
        Tolerances: exact (ids/strings)
        """
        model = BDF(debug=None)
        for nid in [50, 51, 52, 101, 150, 999]:
            model.add_grid(nid, [float(nid), 0., 0.])
        model.add_card(['RCROSS', 20, 'DISP', 50, 2, 'STRESS', 150, 8, 4], 'RCROSS',
                       comment='rcross')
        model.add_card(['RCROSS', 200, 'PRESS', 101, 1, 'DISP', 50, 2], 'RCROSS')
        rcross_op2 = model.add_rcross(20, 'VELO', 51, 3, None, 52, 1, 0)
        self.assertIsNone(rcross_op2.curid)
        self.assertEqual(len(model.rcrosses[20]), 2)

        rcross1 = model.rcrosses[20][0]
        self.assertEqual(rcross1.rtype2, 'STRESS')
        self.assertEqual(rcross1.curid, 4)
        rcross2 = model.rcrosses[200][0]
        self.assertIsNone(rcross2.curid)

        model.validate()
        for rcrosses in model.rcrosses.values():
            for rcross in rcrosses:
                rcross.write_card(size=8)
                rcross.write_card(size=16)

        bad_rtype = model.rcrosses[20][0].__class__(1, 'DISPL', 1, 1, 'DISP', 2, 1)
        with self.assertRaises(AssertionError):
            bad_rtype.validate()
        both_blank = model.rcrosses[20][0].__class__(1, None, 1, 1, None, 2, 1)
        with self.assertRaises(AssertionError):
            both_blank.validate()

        model2 = save_load_deck(model, run_remove_unused=False)
        self.assertEqual(sorted(model2.rcrosses), [20, 200])
        self.assertEqual(len(model2.rcrosses[20]), 2)
        self.assertIsNone(model2.rcrosses[20][1].rtype2)
        self.assertEqual(model2.rcrosses[20][1].comp1, 3)
        self.assertEqual(model2.rcrosses[20][1].response_types_ids(),
                         [('VELO', 51), ('VELO', 52)])

        # grid 150 is only referenced as a STRESS element id -> removed
        remove_unused(model2)
        self.assertEqual(sorted(model2.nodes), [50, 51, 52, 101])

    def test_elist(self) -> None:
        """
        Tests the ELIST card.

        Tests:
         - parsing the NX QRG example (signed ids + THRU)
         - THRU expansion keeps the sign; xref keeps only the THRU ids
           that exist, but explicitly listed ids must exist
         - validate rejects eid=0, leading/trailing THRU, and THRU
           ranges with mixed signs
         - the writer never puts THRU in the first/last data field of
           a row (NX Remark 3), and the ids survive a write/read
         - write/read round trip (save_load_deck)
        Tolerances: exact (integer ids)
        """
        model = BDF(debug=None)
        _add_quads(model, [51, 62, 68, 70, 102, 122])

        model.add_card(['ELIST', 3, 51, -62, 68, 'THRU', 102, 122], 'ELIST',
                       comment='elist')
        elist = model.elists[3]
        self.assertEqual(elist.eids, [51, -62, 68, 'THRU', 102, 122])
        signed = elist.signed_element_ids
        self.assertEqual(signed[:3], [51, -62, 68])
        self.assertEqual(len(signed), 3 + (102 - 68) + 1)
        self.assertEqual(signed[-1], 122)
        model.validate()

        model.cross_reference()
        self.assertEqual(elist.element_ids, [51, 62, 68, 70, 102, 122])
        model.uncross_reference()
        model.safe_cross_reference()

        # negative THRU range keeps the sign
        elist_neg = model.add_elist(4, [-68, 'THRU', -70])
        self.assertEqual(elist_neg.signed_element_ids, [-68, -69, -70])

        # explicitly listed missing element -> error; THRU gap -> ok
        model.add_elist(5, [51, 999])
        model.uncross_reference()
        with self.assertRaises(KeyError):
            model.elists[5].cross_reference(model)
        del model.elists[5]

        cls = elist.__class__
        for bad_eids in ([0], ['THRU', 2], [1, 'THRU'], [1, 'THRU', -3], [5, 'THRU', 3]):
            with self.assertRaises(AssertionError):
                cls(10, bad_eids).validate()

        # THRU must not land in the first/last data field of a row
        eids = [1, 'THRU', 3, 5, 'THRU', 7, 9, 'THRU', 11, 20, 21, 30, 'THRU', 33]
        elist_layout = cls(11, eids)
        elist_layout.validate()
        msg = elist_layout.write_card(size=8)
        for irow, line in enumerate(msg.rstrip().split('\n')):
            row = [line[i:i+8].strip() for i in range(0, len(line), 8)]
            ifirst = 2 if irow == 0 else 1
            self.assertNotEqual(row[ifirst], 'THRU', msg)
            if len(row) == 9:
                self.assertNotEqual(row[8], 'THRU', msg)
        model.add_card(msg.rstrip().split('\n'), 'ELIST', is_list=False)
        self.assertEqual(sorted(model.elists[11].signed_element_ids),
                         sorted(elist_layout.signed_element_ids))
        del model.elists[11]  # layout check only; its elements don't exist

        for elisti in model.elists.values():
            elisti.write_card(size=8)
            elisti.write_card(size=16)

        model.cross_reference()
        model2 = save_load_deck(model, run_remove_unused=False)
        self.assertEqual(model2.elists[3].eids, [51, -62, 68, 'THRU', 102, 122])
        self.assertEqual(model2.elists[4].signed_element_ids, [-68, -69, -70])

    def test_mfluid(self) -> None:
        """
        Tests the MFLUID card.

        Tests:
         - parsing the NX QRG example (with continuation: RMAX)
         - defaults (CID=0, ZFS=None, PLANE='N', RMAX=1e10, FMEXACT=1e15)
         - add_mfluid, two fluid volumes sharing one SID
         - validate rejects ELIST1+ELIST2=0 and bad PLANE values
         - CORD2R and ELIST cross-reference / uncross-reference
         - missing ELIST raises on cross_reference
         - write (8/16) and write/read round trip (save_load_deck)
        Tolerances: exact (ids/strings); floats via assertAlmostEqual
        (8-character field writer)
        """
        model = BDF(debug=None)
        model.add_cord2r(2, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.])
        _add_quads(model, [1, 2, 3, 4])
        model.add_elist(3, [1, -2])
        model.add_elist(4, [3])
        model.add_elist(5, [4])
        card = ['MFLUID', 3, 2, 15.73, 1006., 3, 4, 'S', 'N',
                100.]
        model.add_card(card, 'MFLUID', comment='mfluid')
        mfluid1 = model.mfluids[3][0]
        self.assertEqual(mfluid1.cid, 2)
        self.assertAlmostEqual(mfluid1.zfs, 15.73)
        self.assertAlmostEqual(mfluid1.rho, 1006.)
        self.assertEqual((mfluid1.elist1, mfluid1.elist2), (3, 4))
        self.assertEqual((mfluid1.plane1, mfluid1.plane2), ('S', 'N'))
        self.assertAlmostEqual(mfluid1.rmax, 100.)
        self.assertAlmostEqual(mfluid1.fmexact, 1.0e15)

        mfluid2 = model.add_mfluid(3, 1.2, elist1=5)
        self.assertEqual(mfluid2.cid, 0)
        self.assertIsNone(mfluid2.zfs)
        self.assertEqual((mfluid2.plane1, mfluid2.plane2), ('N', 'N'))
        self.assertAlmostEqual(mfluid2.rmax, 1.0e10)
        self.assertEqual(len(model.mfluids[3]), 2)
        self.assertNotIn('1.+10', mfluid2.write_card(size=8))  # defaults are blank

        model.validate()
        for mfluid in model.mfluids[3]:
            mfluid.write_card(size=8)
            mfluid.write_card(size=16)

        no_elist = mfluid2.__class__(1, 1.0)
        with self.assertRaises(AssertionError):
            no_elist.validate()
        bad_plane = mfluid2.__class__(1, 1.0, elist1=1, plane1='X')
        with self.assertRaises(AssertionError):
            bad_plane.validate()

        model.cross_reference()
        self.assertEqual(mfluid1.cid_ref.cid, 2)
        self.assertEqual(mfluid1.Cid(), 2)
        self.assertIs(mfluid1.elist1_ref, model.elists[3])
        self.assertIs(mfluid1.elist2_ref, model.elists[4])
        self.assertEqual((mfluid1.Elist1(), mfluid1.Elist2()), (3, 4))
        self.assertIsNone(mfluid2.elist2_ref)  # ELIST2=0 is unused
        model.uncross_reference()
        self.assertIsNone(mfluid1.cid_ref)
        self.assertIsNone(mfluid1.elist1_ref)
        self.assertEqual((mfluid1.elist1, mfluid1.elist2), (3, 4))
        model.safe_cross_reference()

        missing_elist = mfluid2.__class__(9, 1.0, elist1=99)
        with self.assertRaises(KeyError):
            missing_elist.cross_reference(model)

        model2 = save_load_deck(model, run_remove_unused=False)
        self.assertEqual(len(model2.mfluids[3]), 2)
        mfluid1b = model2.mfluids[3][0]
        self.assertAlmostEqual(mfluid1b.zfs, 15.73)
        self.assertAlmostEqual(mfluid1b.rmax, 100.)
        self.assertEqual(mfluid1b.plane1, 'S')

        # the MFLUID coordinate system is "used"
        remove_unused(model2)
        self.assertIn(2, model2.coords)
        self.assertEqual(sorted(model2.elists), [3, 4, 5])


def _add_quads(model: BDF, eids: list[int]) -> None:
    """adds one unit CQUAD4 per element id (with a PSHELL/MAT1)"""
    model.add_mat1(1, 3.0e7, None, 0.3)
    model.add_pshell(1, mid1=1, t=0.1)
    for eid in eids:
        nid = 4 * eid
        model.add_grid(nid, [0., 0., float(eid)])
        model.add_grid(nid + 1, [1., 0., float(eid)])
        model.add_grid(nid + 2, [1., 1., float(eid)])
        model.add_grid(nid + 3, [0., 1., float(eid)])
        model.add_cquad4(eid, 1, [nid, nid + 1, nid + 2, nid + 3])


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
