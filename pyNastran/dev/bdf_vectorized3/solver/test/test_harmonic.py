import os
from pathlib import Path
import unittest
import numpy as np
from cpylog import SimpleLogger
import pyNastran
from pyNastran.dev.bdf_vectorized3.solver.solver import Solver, partition_vector2
from pyNastran.dev.bdf_vectorized3.bdf import BDF, read_bdf

#from pyNastran.bdf.cards.params import PARAM
#from pyNastran.f06.errors import FatalError
from pyNastran.bdf.case_control_deck import CaseControlDeck, Subcase

PKG_PATH = Path(pyNastran.__path__[0])
TEST_DIR = Path(__file__).parent / 'examples'


class TestHarmonic(unittest.TestCase):
    def test_spring_mass_damper_direct(self):
        """spring-mass-damper problem"""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'freq_direct_response.bdf'
        mid = 3
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [2., 0., 0.])

        eid = 10
        pid = 20
        nids = [1, 2]
        #x = [1., 0., 0.]
        #g0 = None
        #model.add_cbar(eid, pid, nids, x, g0, offt='GGG',
                       #pa=0, pb=0, wa=None, wb=None, comment='', validate=False)
        #mid = 30
        #Type = 'ROD'
        #dim = [0.5]
        #model.add_pbarl(pid, mid, Type, dim, group='MSCBML0', nsm=0., comment='')

        k = 1e6
        model.add_celas2(eid, k, nids, c1=1, c2=1, ge=0., s=0., comment='')

        eid += 1
        nids = [2, 3]
        model.add_celas2(eid, k, nids, c1=1, c2=1, ge=0., s=0., comment='')

        #E = 3.0e7
        #G = None
        #nu = 0.3
        #model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0,
                       #St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')

        eid += 1
        nid = 2
        mass = 2.0
        model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')

        eid += 1
        nid = 3
        model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')

        sid = 103
        nmodes = 2
        model.add_eigrl(sid, v1=None, v2=None, nd=nmodes, msglvl=0,
                        maxset=None, shfscl=None, norm=None, options=None, values=None,comment='')

        darea_id = 10
        component = 1
        scale = 1.0
        model.add_darea(darea_id, nid, component, scale, comment='')

        sid = 3
        excite_id = darea_id
        model.add_rload1(sid, excite_id, delay=0, dphase=0, tc=1., td=0, load_type='LOAD', comment='')

        spc_id = 3
        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')

        freq_id = 100
        f1 = 1.
        df = 1.
        ndf = 2000
        model.add_freq1(freq_id, f1, df, ndf=ndf, comment='')

        freqs = [1000.5, 0.5]
        model.add_freq(freq_id, freqs)

        f1 = 20.
        f2 = 40.
        model.add_freq2(freq_id, f1, f2, nf=6, comment='')

        f1 = 10.
        f2 = 2000.
        model.add_freq3(freq_id, f1, f2, Type='LINEAR', nef=10, cluster=1.0, comment='')

        f1 = 500.
        model.add_freq3(freq_id, f1, f2=None, Type='LINEAR', nef=10, cluster=1.0, comment='')

        model.add_freq4(freq_id, f1=1., f2=2000., fspread=0.1, nfm=3, comment='')

        fractions = [0.6, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2]
        model.add_freq5(freq_id, fractions, f1=20., f2=200., comment='')

        setup_direct_frequency_case_control(model)
        solver = Solver(model)
        #with self.assertRaises(RuntimeError):
        solver.run()

    def test_spring_mass_damper_freq_response(self):
        """spring-mass-damper problem"""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'freq_response.bdf'
        # mid = 3
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [2., 0., 0.])

        eid = 10
        # pid = 20
        nids = [1, 2]
        #x = [1., 0., 0.]
        #g0 = None
        #model.add_cbar(eid, pid, nids, x, g0, offt='GGG',
                       #pa=0, pb=0, wa=None, wb=None, comment='', validate=False)
        #mid = 30
        #Type = 'ROD'
        #dim = [0.5]
        #model.add_pbarl(pid, mid, Type, dim, group='MSCBML0', nsm=0., comment='')

        k = 1e6
        model.add_celas2(eid, k, nids, c1=1, c2=1, ge=0., s=0., comment='')

        eid += 1
        nids = [2, 3]
        model.add_celas2(eid, k, nids, c1=1, c2=1, ge=0., s=0., comment='')

        #E = 3.0e7
        #G = None
        #nu = 0.3
        #model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0,
                       #St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')

        eid += 1
        nid = 2
        mass = 3.0
        model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')

        eid += 1
        nid = 3
        model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')

        sid = 103
        nmodes = 2
        model.add_eigrl(sid, v1=None, v2=None, nd=nmodes, msglvl=0,
                        maxset=None, shfscl=None, norm=None, options=None, values=None, comment='')

        darea_id = 10
        component = 1
        scale = 1.0
        model.add_darea(darea_id, nid, component, scale, comment='')

        sid = 3
        excite_id = darea_id
        model.add_rload1(sid, excite_id, delay=0, dphase=0, tc=1., td=0, load_type='LOAD', comment='')

        spc_id = 3
        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')

        freq_id = 100
        f1 = 1.
        df = 1.
        ndf = 2000
        model.add_freq1(freq_id, f1, df, ndf=ndf, comment='')

        freqs = [1000.5, 0.5]
        model.add_freq(freq_id, freqs)

        f1 = 20.
        f2 = 40.
        model.add_freq2(freq_id, f1, f2, nf=6, comment='')

        f1 = 10.
        f2 = 2000.
        model.add_freq3(freq_id, f1, f2, Type='LINEAR', nef=10, cluster=1.0, comment='')

        f1 = 500.
        model.add_freq3(freq_id, f1, f2=None, Type='LINEAR', nef=10, cluster=1.0, comment='')

        model.add_freq4(freq_id, f1=1., f2=2000., fspread=0.1, nfm=3, comment='')

        fractions = [0.6, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2]
        model.add_freq5(freq_id, fractions, f1=20., f2=200., comment='')

        setup_frequency_response_case_control(model)
        solver = Solver(model)
        #with self.assertRaises(RuntimeError):
        solver.run()


def setup_direct_frequency_case_control(model: BDF, extra_case_lines=None):
    lines = [
        'STRESS(PLOT,PRINT) = ALL',
        'STRAIN(PLOT,PRINT) = ALL',
        'FORCE(PLOT,PRINT) = ALL',
        'DISP(PLOT,PRINT) = ALL',
        #'GPFORCE(PLOT,PRINT) = ALL',
        'SPCFORCE(PLOT,PRINT) = ALL',
        'MPCFORCE(PLOT,PRINT) = ALL',
        #'OLOAD(PLOT,PRINT) = ALL',
        #'ESE(PLOT,PRINT) = ALL',
        'SUBCASE 1',
        '  LOAD = 2',
        '  SPC = 3',
        '  METHOD = 103',
        '  FREQUENCY = 100',
    ]
    if extra_case_lines is not None:
        lines += extra_case_lines
    cc = CaseControlDeck(lines, log=model.log)
    model.sol = 108
    model.case_control_deck = cc


def setup_frequency_response_case_control(model: BDF, extra_case_lines=None):
    lines = [
        'STRESS(PLOT,PRINT) = ALL',
        'STRAIN(PLOT,PRINT) = ALL',
        'FORCE(PLOT,PRINT) = ALL',
        'DISP(PLOT,PRINT) = ALL',
        'SET 2 = 2,3',
        'SET 3 = 2',
        'VELO(PLOT,PRINT) = 2',
        'ACCEL(PLOT,PRINT) = 3',
        #'GPFORCE(PLOT,PRINT) = ALL',
        'SPCFORCE(PLOT,PRINT) = ALL',
        # 'MPCFORCE(PLOT,PRINT) = ALL',
        'OLOAD(PLOT,PRINT) = ALL',
        'FORCE(PLOT,PRINT) = ALL',
        #'ESE(PLOT,PRINT) = ALL',
        'SUBCASE 1',
        '  LOAD = 2',
        '  SPC = 3',
        '  METHOD = 103',
        '  FREQUENCY = 100',
    ]
    if extra_case_lines is not None:
        lines += extra_case_lines
    cc = CaseControlDeck(lines, log=model.log)
    model.sol = 111
    model.case_control_deck = cc


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
