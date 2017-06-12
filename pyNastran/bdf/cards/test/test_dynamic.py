from __future__ import print_function
import unittest
from six.moves import StringIO
import numpy as np

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')

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

    def test_delay(self):
        """tests a two field DELAY card"""
        model = BDF(debug=False)
        node1, c1, t1 = 100, 3, 0.3
        node2, c2, t2 = 101, 4, 0.4
        sid = 42
        card_lines = ['DELAY', sid, node1, c1, t1, node2, c2, t2]
        model.add_card(card_lines, card_lines[0], comment='', is_list=True,
                       has_none=True)
        model.add_grid(100)
        model.add_grid(101)
        model.validate()
        model.cross_reference()
        #print(model.delays[42])

    def test_dphase(self):
        """tests a two field DPHASE card"""
        model = BDF(debug=False)
        node1, c1, t1 = 100, 3, 0.3
        node2, c2, t2 = 101, 4, 0.4
        sid = 42
        card_lines = ['DPHASE', sid, node1, c1, t1, node2, c2, t2]
        model.add_card(card_lines, card_lines[0], comment='', is_list=True,
                       has_none=True)
        model.add_grid(100)
        model.add_grid(101)
        model.validate()
        model.cross_reference()
        #print(model.dphases[42])

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

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        out = bdf_file.getvalue()
        bdf_file.seek(0)

        model2 = read_bdf(bdf_file, punch=True, debug=False)
        model2.uncross_reference()
        model2.safe_cross_reference()
        model2.uncross_reference()

    def test_tload(self):
        """tests DLOAD, TLOAD1, TLOAD2, TABLED2 cards"""
        model = BDF(debug=False)
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
        out = bdf_file.getvalue()
        bdf_file.seek(0)
        outs = model.get_bdf_stats(return_type='list')
        outs = model.get_bdf_stats(return_type='string')

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
        #print(out)
        #print(outs)

    def test_rload(self):
        """tests DLOAD, RLOAD1, RLOAD2, TABLED2 cards"""
        model = BDF(debug=False)
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

        rload2.validate()
        rload2.raw_fields()
        rload2.write_card()
        rload2.write_card(size=16)

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
        out = bdf_file.getvalue()
        bdf_file.seek(0)
        outs = model.get_bdf_stats(return_type='list')
        outs = model.get_bdf_stats(return_type='string')

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


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
