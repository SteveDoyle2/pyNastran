import os
from pathlib import Path
import unittest

import numpy as np
from cpylog import SimpleLogger
import pyNastran
from pyNastran.dev.bdf_vectorized3.solver.solver import Solver, partition_vector2
from pyNastran.dev.bdf_vectorized3.bdf import BDF, read_bdf

from pyNastran.bdf.cards.params import PARAM
from pyNastran.f06.errors import FatalError
from pyNastran.bdf.case_control_deck import CaseControlDeck, Subcase
from pyNastran.dev.bdf_vectorized3.solver.test.setup_cc import (
    setup_static_case_control,
    setup_modes_case_control,
    setup_buckling_case_control)

PKG_PATH = Path(pyNastran.__path__[0])
# TEST_DIR = PKG_PATH / 'dev' / 'solver'
TEST_DIR = Path(__file__).parent / 'examples'
MODEL_PATH = PKG_PATH / '..' / 'models'


class TestStaticBar(unittest.TestCase):
    """tests the CBARs"""
    def _test_cbar_rbe3_ex(self):
        """Tests a CBAR/PBAR"""
        log = SimpleLogger(level='debug', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'cbar_rbe3.bdf'
        model.add_grid(100, [3., 3., 0.])
        model.add_grid(101, [3., 3., 4.])
        L = 1.0

        nids = [100, 101]
        eid = 1
        pid = 2
        mid = 3
        E = 8.0e6
        G = 3.2e6
        nu = None
        model.add_mat1(mid, E, G, nu, rho=0.0, alpha=0.0, tref=0.0,
            ge=0.0, St=0.0, Sc=0.0, Ss=0.0, mcsid=0)

        x = [0., 1., 0.]
        g0 = None
        model.add_cbar(eid, pid, nids, x, g0,
                       offt='GGG', pa=0, pb=0, wa=None, wb=None, comment='')

        A = 1.
        k_axial = A * E / L
        model.add_pbar(pid, mid,
                       area=A, i1=1., i2=1., i12=1., j=1.,
                       nsm=0.,
                       c1=0., c2=0., d1=0., d2=0.,
                       e1=0., e2=0., f1=0., f2=0.,
                       k1=1.e8, k2=1.e8, comment='')
        load_id = 2
        nid = 101
        mag = 1.
        fxyz = np.array([0., 0., 1.])
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        spc_id = 3
        components = 123456
        nodes = 100
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

        A = np.array([
            [ 15,   0,  0,  0, -30, 0],
            [  0,  15,  0, 30,   0, 0],
            [  0,   0, 20,  0,   0, 0],
            [  0,  30,  0, 80,   0, 0],
            [-30,   0,  0,  0,  80, 0],
            [  0,   0,  0,  0,   0, 8],
        ])

        B = np.array([
            [-15,   0,   0,   0,-30,  0],
            [  0, -15,   0, -30,  0,  0],
            [  0,   0, -20,   0,  0,  0],
            [  0, -30,   0,  40,  0,  0],
            [ 30,   0,   0,   0, 40,  0],
            [  0,   0,   0,   0,  0, -8],
        ])

        C = np.array([
            [-15,   0,   0,   0, 30,  0],
            [  0, -15,   0, -30,  0,  0],
            [  0,   0, -20,   0,  0,  0],
            [  0,  30,   0,  40,  0,  0],
            [-30,   0,   0,   0, 40,  0],
            [  0,   0,   0,   0,  0, -8],
        ])

        D = np.array([
            [ 15,   0,  0,  0,  30, 0],
            [  0,  15,  0,-30,   0, 0],
            [  0,   0, 20,  0,   0, 0],
            [  0, -30,  0, 80,   0, 0],
            [ 30,   0,  0,  0,  80, 0],
            [  0,   0,  0,  0,   0, 8],
        ])
        Kgg_expected = np.vstack([
            np.hstack([A, B,]),
            np.hstack([C, D]),
        ])
        Kgg = solver.Kgg / 1e5
        np.set_printoptions(linewidth=100)
        A_actual = Kgg[:6, :6] ###
        B_actual = Kgg[:6, 6:]
        C_actual = Kgg[6:, :6] ###
        D_actual = Kgg[6:, 6:]

        # F = kx
        fmag = 1.0
        dx = fmag / k_axial
        nnodes_6 = len(solver.xg)
        nnodes = nnodes_6 // 6
        xg = solver.xg.reshape(nnodes, 6)
        assert np.allclose(dx, xg[1, 2]), f'dx={dx} xgi={xg[1,2]} xg={xg}'

    def test_cbar(self):
        """Tests a CBAR/PBAR"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'cbar.bdf'
        k_axial = build_static_cbar(model)

        load_id = 2
        nid = 2
        mag = 1.
        fxyz = np.array([1., 0., 0.])
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        solver = Solver(model)
        solver.run()

        # F = kx
        fmag = 1.0
        dx = fmag / k_axial
        xg = solver.xg
        assert dx == xg[6], f'dx={dx} xg={xg}'

    def test_cbar_grav(self):
        """Tests a CBAR/PBAR"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'cbar_grav.bdf'
        k_axial = build_static_cbar(model)

        load_id = 2
        nid = 2
        mag = 1.
        gxyz = np.array([1., 0., 0.])
        model.add_grav(load_id, scale=mag, N=gxyz, cid=0)

        solver = Solver(model)
        solver.run()
        read_bdf(solver._bdf_filename)

        # F = kx
        #fmag = 1.0
        #dx = fmag / k_axial
        #xg = solver.xg
        #assert dx == xg[6], f'dx={dx} xg={xg}'

    def test_cbar_pload1(self):
        """Tests a CBAR/PBAR"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'cbar_pload1.bdf'
        k_axial = build_static_cbar(model)

        load_id = 2
        nid = 2
        mag = 1.
        fxyz = np.array([1., 0., 0.])
        eid = 1

        # constant
        model.add_pload1(
            load_id, eid, x1=0., x2=1.0, p1=1., p2=1.,
            load_type='FZ', scale='FR',)
        model.add_pload1(
            load_id, eid, x1=0., x2=1.0, p1=1., p2=1.,
            load_type='FY', scale='FR',)
        model.add_pload1(
            load_id, eid, x1=0., x2=1.0, p1=1., p2=1.,
            load_type='FX', scale='FR',)

        # ramp
        model.add_pload1(
            load_id, eid, x1=0., x2=1.0, p1=2., p2=1.,
            load_type='FZ', scale='FR',)
        model.add_pload1(
            load_id, eid, x1=0., x2=1.0, p1=2., p2=1.,
            load_type='FY', scale='FR',)
        model.add_pload1(
            load_id, eid, x1=0., x2=1.0, p1=2., p2=1.,
            load_type='FX', scale='FR',)

        # point
        model.add_pload1(
            load_id, eid, x1=0., x2=1.0, p1=0.4, p2=0.4,
            load_type='FZ', scale='FR',)
        model.add_pload1(
            load_id, eid, x1=0., x2=1.0, p1=0.4, p2=0.4,
            load_type='FY', scale='FR',)
        model.add_pload1(
            load_id, eid, x1=0., x2=1.0, p1=0.4, p2=0.4,
            load_type='FX', scale='FR',)

        solver = Solver(model)
        solver.run()
        read_bdf(solver._bdf_filename)

        # F = kx
        #fmag = 1.0
        #dx = 1. # fmag / k_axial
        #xg = solver.xg
        #assert dx == xg[6], f'dx={dx} xg={xg}'

    def test_cbar2(self):
        """Tests a CBAR/PBAR"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'cbar.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0.5, 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        L = 1.0
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)

        eid = 1
        x = [0., 1., 0.]
        g0 = None
        nids = [1, 2]
        model.add_cbar(eid, pid, nids, x, g0,
                       offt='GGG', pa=0, pb=0, wa=None, wb=None, comment='')

        eid = 2
        nids = [2, 3]
        model.add_cbar(eid, pid, nids, x, g0,
                       offt='GGG', pa=0, pb=0, wa=None, wb=None, comment='')

        A = 1.
        k_axial = A * E / L
        model.add_pbar(pid, mid,
                       area=A, i1=1., i2=1., i12=1., j=1.,
                       nsm=0.,
                       c1=0., c2=0., d1=0., d2=0.,
                       e1=0., e2=0., f1=0., f2=0.,
                       k1=1.e8, k2=1.e8, comment='')
        load_id = 2
        spc_id = 3
        nid = 3
        mag = 1.
        fxyz = np.array([1., 0., 0.])
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

        # F = kx
        fmag = 1.0
        dx = fmag / k_axial
        xg = solver.xg
        assert dx == xg[6*2], f'dx={dx} xg={xg}'

    def test_cbeam(self):
        """Tests a CBEAM/PBEAM"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'cbeam.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        L = 1.0
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)

        x = [0., 1., 0.]
        g0 = None
        model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                        pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment='')
        #beam_type = 'BAR'
        xxb = [0.]
        #dims = [[1., 2.]]
        #model.add_pbeaml(pid, mid, beam_type, xxb, dims,
                         #so=None, nsm=None, group='MSCBML0', comment='')
        A = 1.0
        so = ['YES']
        area = [1.0]
        i1 = [1.]
        i2 = [2.]
        i12 = [0.]
        j = [2.]
        model.add_pbeam(pid, mid, xxb, so, area, i1, i2, i12, j, nsm=None,
                        c1=None, c2=None, d1=None, d2=None, e1=None, e2=None, f1=None, f2=None,
                        k1=1., k2=1., s1=0., s2=0., nsia=0., nsib=None, cwa=0., cwb=None,
                        m1a=0., m2a=0., m1b=None, m2b=None, n1a=0., n2a=0., n1b=None, n2b=None,
                        comment='')
        load_id = 2
        spc_id = 3
        nid = 2
        mag = 1.
        fxyz = np.array([1., 0., 0.])
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

        # F = kx
        k_axial = A * E / L
        fmag = 1.0
        dx = fmag / k_axial
        xg = solver.xg
        assert dx == xg[6], f'dx={dx} xg={xg}'


    def test_cbeam2(self):
        """Tests a CBEAM/PBEAM"""
        log = SimpleLogger(level='warning')
        model = BDF(debug=True, log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'cbeam.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [.5, 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        L = 1.0
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)

        eid = 1
        nids = [1, 2]
        x = [0., 1., 0.]
        g0 = None
        model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                        pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment='')

        eid = 2
        nids = [2, 3]
        x = [0., 1., 0.]
        g0 = None
        model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                        pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment='')


        #beam_type = 'BAR'
        xxb = [0.]
        #dims = [[1., 2.]]
        #model.add_pbeaml(pid, mid, beam_type, xxb, dims,
                         #so=None, nsm=None, group='MSCBML0', comment='')
        A = 1.0
        so = ['YES']
        area = [1.0]
        i1 = [1.]
        i2 = [2.]
        i12 = [0.]
        j = [2.]
        model.add_pbeam(pid, mid, xxb, so, area, i1, i2, i12, j, nsm=None,
                        c1=None, c2=None, d1=None, d2=None, e1=None, e2=None, f1=None, f2=None,
                        k1=1., k2=1., s1=0., s2=0., nsia=0., nsib=None, cwa=0., cwb=None,
                        m1a=0., m2a=0., m1b=None, m2b=None, n1a=0., n2a=0., n1b=None, n2b=None,
                        comment='')
        load_id = 2
        spc_id = 3
        nid = 3
        mag = 1.
        fxyz = np.array([1., 0., 0.])
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

        # F = kx
        k_axial = A * E / L
        fmag = 1.0
        dx = fmag / k_axial

        nnode = 3
        xg = solver.xg.reshape(nnode, 6)
        dx_actual = solver.xg[6*2]
        assert dx == xg[2, 0], f'dx={dx} dx_actual={dx_actual} xg:\n{xg}'


def build_static_cbar(model: BDF):
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    L = 1.0

    nids = [1, 2]
    eid = 1
    pid = 2
    mid = 3
    E = 3.0e7
    G = None
    nu = 0.3
    model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                   Sc=0.0, Ss=0.0, mcsid=0)

    x = [0., 1., 0.]
    g0 = None
    model.add_cbar(eid, pid, nids, x, g0,
                   offt='GGG', pa=0, pb=0, wa=None, wb=None, comment='')

    A = 1.
    k_axial = A * E / L
    model.add_pbar(pid, mid,
                   area=A, i1=1., i2=1., i12=1., j=1.,
                   nsm=0.,
                   c1=0., c2=0., d1=0., d2=0.,
                   e1=0., e2=0., f1=0., f2=0.,
                   k1=1.e8, k2=1.e8, comment='')
    spc_id = 3

    components = 123456
    nodes = 1
    model.add_spc1(spc_id, components, nodes, comment='')
    setup_static_case_control(model)
    return k_axial


if __name__ == '__main__':
    import unittest
    unittest.main()
