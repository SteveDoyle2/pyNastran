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

PKG_PATH = Path(pyNastran.__path__[0])
# TEST_DIR = PKG_PATH / 'dev' / 'solver'
TEST_DIR = Path(__file__).parent / 'examples'
MODEL_PATH = PKG_PATH / '..' / 'models'


def setup_static_case_control(model: BDF, extra_case_lines=None):
    lines = [
        'STRESS(PLOT,PRINT) = ALL',
        'STRAIN(PLOT,PRINT) = ALL',
        'FORCE(PLOT,PRINT) = ALL',
        'DISP(PLOT,PRINT) = ALL',
        'GPFORCE(PLOT,PRINT) = ALL',
        'SPCFORCE(PLOT,PRINT) = ALL',
        'MPCFORCE(PLOT,PRINT) = ALL',
        'OLOAD(PLOT,PRINT) = ALL',
        'ESE(PLOT,PRINT) = ALL',
        'SUBCASE 1',
        '  LOAD = 2',
        '  SPC = 3',
    ]
    if extra_case_lines is not None:
        lines += extra_case_lines
    cc = CaseControlDeck(lines, log=model.log)
    model.sol = 101
    model.case_control_deck = cc


class TestSolverTools(unittest.TestCase):
    def test_partition_vector(self):
        xg = np.array([0., 1., 2., 3., 4.], dtype='float64')
        gset = np.ones(5, dtype='bool')
        sset = np.array([True, False, True, False, False])
        aset = gset & ~sset
        xa, xs = partition_vector2(xg, [['a', aset], ['s', sset]])
        assert len(xa) == 3 and np.allclose(xa, [1., 3., 4.]), xa
        assert len(xs) == 2 and np.allclose(xs, [0., 2.]), xs


class TestStaticSpring(unittest.TestCase):
    def test_celas2_conm2(self):
        """Tests a CELAS2/CONM2"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'celas1.bdf'
        model.add_grid(1, [0., 0., 0.], cd=1)
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0.5, 1., 0.], cd=3)
        model.add_grid(4, [0.5, 0., 1.])

        origin = [0., 0., 0.]
        zaxis = [0., 0., 1]
        xzplane = [1., 1., 0.]
        model.add_cord2r(1, origin, zaxis, xzplane, rid=0, setup=True, comment='')

        origin = [0.5, 1., 0.]
        zaxis = [0., 0., 1]
        xzplane = [0.5, 2., 0.]
        model.add_cord2r(3, origin, zaxis, xzplane, rid=0, setup=True, comment='')

        model.add_conm2(1, 1, mass=2.0, cid=0, X=None, I=None, comment='')
        model.add_conm2(2, 2, mass=3.0, cid=0, X=None, I=None, comment='')
        model.add_conm2(3, 3, mass=3.0, cid=0, X=None, I=None, comment='')
        model.add_conm2(4, 4, mass=5.0, cid=0, X=None, I=None, comment='')

        nids = [1, 2]
        eid = 1
        k = 1000.
        model.add_celas2(eid, k, nids, c1=1, c2=2, ge=0., s=0., comment='')

        load_id = 2
        spc_id = 3
        # model.add_sload(load_id, 2, 20.)
        fxyz = np.array([1., 0., 0.])
        mag = 20.
        model.add_force(load_id, 2, mag, fxyz, cid=0, comment='')

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)

        solver = Solver(model)
        model.sol = 101
        solver.run()

    def _test_celas2_cd(self):
        """Tests a CELAS2"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')

        model.bdf_filename = TEST_DIR / 'celas2.bdf'
        model.add_grid(1, [0., 0., 0.])
        cd = 100
        model.add_grid(2, [0., 0., 0.], cd=cd)
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [0., 1., 0.]
        model.add_cord2r(cd, origin, zaxis, xzplane, rid=0, setup=True, comment='')
        nids = [1, 2]
        eid = 1
        k = 1000.
        model.add_celas2(eid, k, nids, c1=1, c2=2, ge=0., s=0., comment='')

        load_id = 2
        spc_id = 3
        # model.add_sload(load_id, 2, 20.)
        fxyz = np.array([1., 0., 0.])
        mag = 20.
        model.add_force(load_id, 2, mag, fxyz, cid=0, comment='')

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)

        solver = Solver(model)
        solver.run()

        # F = k * d
        d = mag / k
        Fg = solver.Fg
        Kgg = solver.Kgg
        Kaa2 = solver.Kaa.toarray()
        Fg_expected = [ 0., 0.0, 0., 0., 0., 0.,
                        20., 0., 0., 0., 0., 0.]
        assert np.allclose(Fg, Fg_expected), f'Force Error:\n Fg2={Fg}\n Fg_expected={Fg_expected}'
        #assert np.allclose(Kgg, Kgg)
        #assert np.allclose(Kaa, Kaa)
        assert np.allclose(solver.xa_[0], d)

    def test_celas3(self):
        """Tests a CELAS3/PELAS"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'celas3.bdf'
        #model.add_grid(1, [0., 0., 0.])
        #model.add_grid(2, [0., 0., 0.])
        model.add_spoint([1, 2])
        nids = [1, 2]
        eid = 1
        pid = 2
        model.add_celas3(eid, pid, nids, comment='')
        k = 1000.

        load_id = 2
        spc_id = 3
        model.add_pelas(pid, k, ge=0., s=0., comment='')
        # model.add_sload(load_id, 2, 20.)
        #fxyz = np.array([1., 0., 0.])
        mag = 20.
        # model.add_force(load_id, 2, mag, fxyz, cid=0, comment='')
        nids = 2
        mags = mag
        model.add_sload(load_id, nids, mags, comment='')

        components = 0
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)

        solver = Solver(model)
        #model.sol = 103
        #solver.run()

        model.sol = 101
        solver.run()

        # F = k * d
        # d = F / k
        # 20=1000.*d
        d = mag / k
        assert np.allclose(solver.xa_[0], d), solver.xa_
        os.remove(solver.f06_filename)
        os.remove(solver.op2_filename)

    def test_celas4_cd(self):
        """Tests a CELAS4"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'celas4.bdf'
        model.add_spoint([1, 2])
        #origin = [0., 0., 0.]
        #zaxis = [0., 0., 1.]
        #xzplane = [0., 1., 0.]
        #model.add_cord2r(cd, origin, zaxis, xzplane, rid=0, setup=True, comment='')
        nids = [1, 2]
        eid = 1
        k = 1000.
        model.add_celas4(eid, k, nids, comment='')

        load_id = 2
        spc_id = 3
        # model.add_sload(load_id, 2, 20.)
        #fxyz = np.array([1., 0., 0.])
        mag = 20.
        model.add_sload(load_id, 2, mag)

        components = 0
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)

        solver = Solver(model)
        solver.run()

        # F = k * d
        d = mag / k
        assert np.allclose(solver.xa_[0], d)

    def test_ese_all_springs(self):
        """Tests ESE for all spring types in one model.

        4 springs in series (one of each type), all with k=1000:
          SPOINT 1 --[CELAS4]-- SPOINT 2 --[CELAS3]-- SPOINT 3
                   --[CELAS2]-- SPOINT 4 --[CELAS1]-- SPOINT 5

        SPCs at endpoints (1, 5), force at the middle junction (3).
        MPC ties SPOINT 2 to SPOINT 4 (making it one continuous chain).

        Equivalent: 4 springs in series, k_eq = k/4 = 250.
        But with both ends fixed and force at mid-chain junction (node 3):
          Left side: 2 springs in series from node 1 to node 3: k_left = k/2
          Right side: 2 springs in series from node 3 to node 5: k_right = k/2
          Combined at node 3: k_total = k_left + k_right = k
          u3 = F / k_total = F / k

        Per-element energies depend on individual dx values.
        Total SE = 0.5 * F * u3
        """
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'test_ese_all_springs.bdf'
        model.add_spoint([1, 2, 3, 4, 5])

        k = 1000.
        # CELAS4: SPOINT 1 -- SPOINT 2 (k on card, no property)
        model.add_celas4(1, k, [1, 2])
        # CELAS3: SPOINT 2 -- SPOINT 3 (uses PELAS)
        model.add_pelas(10, k)
        model.add_celas3(2, 10, [2, 3])
        # CELAS2: SPOINT 3 -- SPOINT 4 (k on card, GRID-like with comp=0)
        model.add_celas2(3, k, [3, 4], c1=0, c2=0)
        # CELAS1: SPOINT 4 -- SPOINT 5 (uses PELAS)
        model.add_pelas(20, k)
        model.add_celas1(4, 20, [4, 5], c1=0, c2=0)

        F = 2.0
        model.add_sload(2, 3, F)
        model.add_spc1(3, 0, [1, 5])

        lines = [
            'DISP(PLOT,PRINT) = ALL',
            'SPCFORCE(PLOT,PRINT) = ALL',
            'ESE(PLOT,PRINT) = ALL',
            'SUBCASE 1',
            '  LOAD = 2',
            '  SPC = 3',
        ]
        cc = CaseControlDeck(lines, log=model.log)
        model.sol = 101
        model.case_control_deck = cc
        solver = Solver(model)
        solver.run()

        # Analytical solution:
        # Chain: 1--[k]--2--[k]--3--[k]--4--[k]--5
        # SPC at 1 and 5, force F at 3.
        # Left pair (springs 1,2): k_left = k/2, carries force F_left
        # Right pair (springs 3,4): k_right = k/2, carries force F_right
        # Equilibrium at node 3: F = k_left*u3 + k_right*u3 = k*u3
        # u3 = F/k
        u3 = F / k
        # Left pair: u2 = u3/2 (symmetric sub-chain with equal springs)
        u2 = u3 / 2
        # Right pair: u4 = u3/2
        u4 = u3 / 2

        # Verify displacements (SPOINT indices: 1->0, 2->1, 3->2, 4->3, 5->4)
        assert np.allclose(solver.xg[2], u3, rtol=1e-10), \
            f"u3={solver.xg[2]}, expected={u3}"

        # Per-element strain energies
        # Spring 1 (1-2): dx = u2 - 0 = u3/2
        # Spring 2 (2-3): dx = u3 - u2 = u3/2
        # Spring 3 (3-4): dx = u4 - u3 = -(u3/2)
        # Spring 4 (4-5): dx = 0 - u4 = -(u3/2)
        se_each = 0.5 * k * (u3 / 2) ** 2
        total_expected = 4 * se_each
        # Also check: total = 0.5 * F * u3
        assert np.allclose(total_expected, 0.5 * F * u3, rtol=1e-10)

        ese_results = solver.op2.op2_results.strain_energy
        all_energies = []
        for attr in ['celas4_strain_energy', 'celas3_strain_energy',
                     'celas2_strain_energy', 'celas1_strain_energy']:
            ese_dict = getattr(ese_results, attr)
            if 1 in ese_dict:
                ese = ese_dict[1]
                energies = ese.data[0, :-1, 0]
                all_energies.extend(energies.tolist())
                for e in energies:
                    assert np.allclose(e, se_each, rtol=1e-10), \
                        f"{attr}: SE={e}, expected={se_each}"

        assert len(all_energies) == 4, f"Expected 4 elements, got {len(all_energies)}"
        assert np.allclose(sum(all_energies), total_expected, rtol=1e-10), \
            f"Total SE={sum(all_energies)}, expected={total_expected}"


    def test_ese_all_springs_cantilever(self):
        """4 springs in series (one of each type), cantilevered.

        SPOINT 1 --[CELAS4,k=1000]-- 2 --[CELAS3,k=500]-- 3
                 --[CELAS2,k=2000]-- 4 --[CELAS1,k=800]-- 5
        SPC at node 1 only. Force F=10 at node 5.
        k_series = 1/(1/k1 + 1/k2 + 1/k3 + 1/k4)
        u5 = F / k_series
        Total SE = 0.5 * F * u5
        Per-element: SE_i = 0.5 * k_i * dx_i^2, where dx_i = F / k_i
        """
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'test_ese_cantilever.bdf'
        model.add_spoint([1, 2, 3, 4, 5])

        k1, k2, k3, k4 = 1000., 500., 2000., 800.
        model.add_celas4(1, k1, [1, 2])
        model.add_pelas(10, k2)
        model.add_celas3(2, 10, [2, 3])
        model.add_celas2(3, k3, [3, 4], c1=0, c2=0)
        model.add_pelas(20, k4)
        model.add_celas1(4, 20, [4, 5], c1=0, c2=0)

        F = 10.0
        model.add_sload(2, 5, F)
        model.add_spc1(3, 0, [1])

        lines = [
            'DISP(PLOT,PRINT) = ALL',
            'SPCFORCE(PLOT,PRINT) = ALL',
            'ESE(PLOT,PRINT) = ALL',
            'SUBCASE 1',
            '  LOAD = 2',
            '  SPC = 3',
        ]
        cc = CaseControlDeck(lines, log=model.log)
        model.sol = 101
        model.case_control_deck = cc
        solver = Solver(model)
        solver.run()

        # All springs carry the same force F (series chain, one free end)
        ks = [k1, k2, k3, k4]
        k_series = 1.0 / sum(1.0 / ki for ki in ks)
        u5_expected = F / k_series

        # Intermediate displacements: u_i = sum(F/k_j for j=1..i)
        u2 = F / k1
        u3 = u2 + F / k2
        u4 = u3 + F / k3
        u5 = u4 + F / k4
        assert np.allclose(u5, u5_expected, rtol=1e-10)

        # Check solver displacements (SPOINT order: 1->0, 2->1, ..., 5->4)
        assert np.allclose(solver.xg[1], u2, rtol=1e-10), \
            f"u2={solver.xg[1]}, expected={u2}"
        assert np.allclose(solver.xg[4], u5, rtol=1e-10), \
            f"u5={solver.xg[4]}, expected={u5}"

        # Per-element SE: each spring stretches by dx_i = F/k_i
        se_expected = [0.5 * ki * (F / ki)**2 for ki in ks]  # = 0.5*F^2/k_i
        total_expected = sum(se_expected)
        assert np.allclose(total_expected, 0.5 * F * u5, rtol=1e-10)

        # Verify from op2 ESE tables
        ese_results = solver.op2.op2_results.strain_energy
        ese_attrs = [
            ('celas4_strain_energy', se_expected[0]),
            ('celas3_strain_energy', se_expected[1]),
            ('celas2_strain_energy', se_expected[2]),
            ('celas1_strain_energy', se_expected[3]),
        ]
        total_from_op2 = 0.0
        for attr, se_ref in ese_attrs:
            ese_dict = getattr(ese_results, attr)
            assert 1 in ese_dict, f"{attr} not found in ESE results"
            ese = ese_dict[1]
            energy = float(ese.data[0, 0, 0])  # single element per type
            assert np.allclose(energy, se_ref, rtol=1e-6), \
                f"{attr}: SE={energy}, expected={se_ref}"
            total_from_op2 += energy

        assert np.allclose(total_from_op2, total_expected, rtol=1e-6), \
            f"Total SE={total_from_op2}, expected={total_expected}"


class TestStaticRod(unittest.TestCase):
    """tests the rods"""
    def test_crod_conrod_ctube_axial(self):
        """Tests a CROD/PROD"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'crod_axial.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(
            mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0,
            St=0.0, Sc=0.0, Ss=0.0, mcsid=0)
        model.add_crod(eid, pid, nids)
        model.add_prod(pid, mid, A=1.0, j=0., c=0., nsm=0.1)

        model.add_conrod(eid, mid, nids, A=1.0)

        model.add_ctube(eid, pid+1, nids)
        model.add_ptube(pid+1, mid, OD1=1.0, t=0.1, nsm=0.1)

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
        #with self.assertRaises(RuntimeError):
        solver.run()

    def test_crod_axial(self):
        """Tests a CROD/PROD"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'crod_axial.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_crod(eid, pid, nids)
        model.add_prod(pid, mid, A=1.0, j=0., c=0., nsm=0.)

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
        #with self.assertRaises(RuntimeError):
        solver.run()

    def test_crod_torsion(self):
        """Tests a CROD/PROD"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'crod_torsion.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_crod(eid, pid, nids)
        model.add_prod(pid, mid, A=0.0, j=2., c=0., nsm=1.)

        load_id = 2
        spc_id = 3
        nid = 2
        mag = 1.
        fxyz = np.array([1., 0., 0.])
        model.add_moment(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

    def _test_crod_spcd(self):
        """Tests a CROD/PROD with an SPCD and no free DOFs"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'crod_spcd.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 42.
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_crod(eid, pid, nids)
        model.add_prod(pid, mid, A=1.0, j=0., c=0., nsm=0.)

        load_id = 2
        spc_id = 3
        nid = 2
        mag = 1.
        fxyz = np.array([0., 0., 0.])
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        nodes = 1
        components = 123456
        #model.add_spc1(spc_id, components, nodes, comment='')
        model.add_spc(spc_id, nodes, components, 0., comment='')

        #components = 123456
        #nodes = 1
        #enforced = 0.1
        #model.add_spc(spc_id, nodes, components, enforced, comment='')
        nodes = 2
        components = 1
        enforced = 0.1
        model.add_spcd(load_id, nodes, components, enforced, comment='')

        components = '1'
        model.add_spcd(9999999, nodes, components, enforced, comment='')

        components = 123456
        #model.add_spc1(spc_id, components, nodes, comment='')
        model.add_spc(spc_id, nodes, components, 0., comment='')
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()
        # F = kx
        dx = enforced
        k = E
        F = k * dx
        assert len(solver.xa_) == 0, solver.xa_
        assert len(solver.Fa_) == 0, solver.Fa_
        assert np.allclose(solver.xg[6], dx), solver.xa_
        assert np.allclose(solver.Fg[6], F), solver.xa_
        #assert np.allclose(solver.Fa_[0], 0.), solver.Fa_

    def test_crod_simple(self):
        """Tests a CROD/PROD"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'crod.bdf'
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
        model.add_mat1(mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_crod(eid, pid, nids)
        A = 1.0
        J = 2.0
        model.add_prod(pid, mid, A=A, j=J, c=0., nsm=0.)

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

        G = model.mat1.G[0]

        ka_expected = A * E / L
        kt_expected = G * J / L
        Kgg_expected = np.zeros((12, 12))
        # axial
        Kgg_expected[0, 0] = Kgg_expected[6, 6] = ka_expected
        Kgg_expected[0, 6] = Kgg_expected[6, 0] = -ka_expected

        # torsion
        Kgg_expected[3, 3] = Kgg_expected[9, 9] = kt_expected
        Kgg_expected[3, 9] = Kgg_expected[9, 3] = -kt_expected

        dof = [0, 3, 6, 9]
        Kgg_expectedi = Kgg_expected[dof, :][:, dof]
        Kgg = solver.Kgg[dof, :][:, dof].toarray()
        assert np.allclose(Kgg, Kgg_expectedi)


    def _test_crod_rotate(self):
        """Tests a CROD/PROD"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        # log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'crod.bdf'
        thetad = 30.
        theta = np.radians(thetad)
        model.add_grid(1, [0., 0., 0.])
        x = c = np.cos(theta)
        y = s = np.sin(theta)
        model.add_grid(2, [x, y, 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        #E = 1.
        #G = 1.
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_crod(eid, pid, nids)
        A = 1.0
        J = 2.0
        L = 1.0

        # ka=1
        # kt=2
        model.add_prod(pid, mid, A=A, j=J, c=0., nsm=0.)

        load_id = 2
        spc_id = 3
        nid = 2
        mag = 1.
        fxyz = np.array([x, y, 0.])
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()
        #---------------------------
        # results
        G = model.mat1.G[0]

        ka_expected = A * E / L
        kt_expected = G * J / L
        Kgg_expected = np.zeros((12, 12))
        # axial
        Kgg_expected[0, 0] = Kgg_expected[6, 6] = ka_expected
        Kgg_expected[0, 6] = Kgg_expected[6, 0] = -ka_expected

        # torsion
        Kgg_expected[3, 3] = Kgg_expected[9, 9] = kt_expected
        Kgg_expected[3, 9] = Kgg_expected[9, 3] = -kt_expected
        T3 = np.array([
            [c, -s, 0.],
            [s, c, 0.],
            [0., 0., 1.]
        ])
        z = np.zeros((3, 3))
        T = np.block([
            [T3, z, z, z],
            [z, T3, z, z],
            [z, z, T3, z],
            [z, z, z, T3],
        ])
        Kgg_expectedi = T.T @ Kgg_expected @ T
        Kgg_expectedi = np.round(Kgg_expectedi)
        #dof = [0, 3, 6, 9]
        #Kgg_expectedi = Kgg_expected[dof, :][:, dof]
        Kgg = np.round(solver.Kgg) # [dof, :][:, dof]
        assert np.allclose(Kgg, Kgg_expectedi)

    def test_crod_aset(self):
        """
        Tests a CROD/PROD using an ASET

        same answer as ``test_crod``
        """
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'crod_aset.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu,
                       rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_crod(eid, pid, nids)
        model.add_prod(pid, mid, A=1.0, j=2., c=0., nsm=0.)

        load_id = 2
        spc_id = 3
        nid = 2
        mag = 1.
        fxyz = np.array([1., 0., 0.])
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')

        aset_ids = 2
        aset_components = '456'
        model.add_aset(aset_ids, aset_components)

        aset_ids = [2, 2]
        aset_components = ['456', '123']
        model.add_aset1(aset_ids, aset_components, comment='')

        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_crod_mpc(self):
        """Tests a CROD/PROD with MPC tying node 2 DOF 1 to node 3 DOF 3.

        MPC: 1.0*u2_1 + (-1.0)*u3_3 = 0  =>  u2_x = u3_z
        Rod from node 1 to 2 (axial along X), force at node 3 in Z.
        Expected: u2_x = u3_z = F / (EA/L) = 1.0 / 3e7
        """
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'crod_mpc.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_crod(eid, pid, nids)
        model.add_prod(pid, mid, A=1.0, j=2., c=0., nsm=0.)
        mpc_id = 10
        nodes = [2, 3]
        components = [1, 3]
        coefficients = [1., -1.]
        model.add_mpc(mpc_id, nodes, components, coefficients)

        load_id = 2
        spc_id = 3
        nid = 3
        mag = 1.
        fxyz = np.array([0., 0., 1.])
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model, extra_case_lines=['MPC=10'])
        solver = Solver(model)
        solver.run()

        # Verify MPC result
        # Rod: K = EA/L = 3e7 * 1.0 / 1.0 = 3e7
        # u2_x = F/K = 1.0 / 3e7
        # MPC: u2_x = u3_z
        A = 1.0
        L = 1.0
        F = 1.0
        expected_disp = F / (E * A / L)
        # DOF indices: node 2 DOF 1 = index 6, node 3 DOF 3 = index 14
        u2_x = solver.xg[6]
        u3_z = solver.xg[14]
        assert np.allclose(u2_x, expected_disp, rtol=1e-10), \
            f"u2_x={u2_x}, expected={expected_disp}"
        assert np.allclose(u3_z, expected_disp, rtol=1e-10), \
            f"u3_z={u3_z}, expected={expected_disp}"
        assert np.allclose(u2_x, u3_z, rtol=1e-10), \
            f"MPC violated: u2_x={u2_x} != u3_z={u3_z}"

    def test_celas_mpc_series(self):
        """Two springs in series connected by an MPC.

        Node 1 (fixed) --[k1=1000]-- Node2 == Node3 --[k2=1000]-- Node4 (fixed)
        MPC: u2_1 = u3_1
        Force F=1.0 at node 2.
        Expected: u_mid = F / (k1 + k2) = 0.0005
        """
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'celas_mpc_series.bdf'
        model.add_spoint([1, 2, 3, 4])

        k1 = 1000.
        k2 = 1000.
        model.add_celas4(1, k1, [1, 2])
        model.add_celas4(2, k2, [3, 4])

        # MPC: 1.0*u2 + (-1.0)*u3 = 0  =>  u2 = u3
        model.add_mpc(10, [2, 3], [0, 0], [1., -1.])

        F = 1.0
        model.add_sload(2, 2, F)

        model.add_spc1(3, 0, [1, 4])

        lines = [
            'DISP(PLOT,PRINT) = ALL',
            'SPCFORCE(PLOT,PRINT) = ALL',
            'SUBCASE 1',
            '  LOAD = 2',
            '  SPC = 3',
            'MPC=10',
        ]
        cc = CaseControlDeck(lines, log=model.log)
        model.sol = 101
        model.case_control_deck = cc

        solver = Solver(model)
        solver.run()

        # u_mid = F / (k1 + k2) for force at the junction with both ends fixed
        expected = F / (k1 + k2)
        # SPOINTs: DOF indices are 0-based by SPOINT ID order
        # SPOINT 1 -> idx 0, SPOINT 2 -> idx 1, SPOINT 3 -> idx 2, SPOINT 4 -> idx 3
        u2 = solver.xg[1]
        u3 = solver.xg[2]
        assert np.allclose(u2, expected, rtol=1e-10), \
            f"u2={u2}, expected={expected}"
        assert np.allclose(u3, expected, rtol=1e-10), \
            f"u3={u3}, expected={expected}"
        assert np.allclose(u2, u3, rtol=1e-10), \
            f"MPC violated: u2={u2} != u3={u3}"

    def test_ctube(self):
        """Tests a CTUBE/PTUBE"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'ctube.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_ctube(eid, pid, nids)
        OD1 = 1.0
        t = 0.1
        model.add_ptube(pid, mid, OD1, t=t, nsm=0., OD2=None, comment='')

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

        ID = OD1 - 2 * t
        # F = k * x
        # x = F / k
        L = 1.0
        A = np.pi * (OD1 ** 2 - ID ** 2) / 4
        kaxial = A * E / L
        F = 1.0
        dx = F / kaxial
        assert np.allclose(solver.xg[6], dx), f'dx={dx} xg[6]={solver.xg[6]}'
        assert np.allclose(solver.Fg[6], F), solver.Fg

    def test_conrod(self):
        """Tests a CONROD"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'conrod.bdf'
        L = 1.
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [L, 0., 0.])

        nids = [1, 2]
        eid = 1
        mid = 3
        E = 200.
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        A = 1.0
        J = 2.0
        model.add_conrod(eid, mid, nids, A=A, j=J, c=0.0, nsm=0.0)

        load_id = 2
        spc_id = 3
        nid = 2
        mag_axial = 10000.
        fxyz = np.array([1., 0., 0.])
        model.add_force(load_id, nid, mag_axial, fxyz, cid=0)

        mag_torsion = 20000.
        fxyz = np.array([1., 0., 0.])
        model.add_moment(load_id, nid, mag_torsion, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')

        #nodes = 2
        #model.add_spc1(spc_id, 23456, nodes, comment='')
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

        # TODO: why is this a torsional result???
        G = E / (2 * (1 + nu))
        kaxial = A * E / L
        ktorsion = G * J / L
        # F = k * d
        daxial = mag_axial / kaxial
        dtorsion = mag_torsion / ktorsion
        #print(solver.xa_)
        assert np.allclose(solver.xa_[0], daxial), f'daxial={daxial} kaxial={kaxial} xa_={solver.xa_}'
        assert np.allclose(solver.xa_[1], dtorsion), f'dtorsion={dtorsion} ktorsion={ktorsion} xa_={solver.xa_}'

        #L = 1.0
        #A = np.pi * (OD1 ** 2 - ID ** 2) / 4
        #kaxial = A * E / L
        #F = 1.0
        #dx = F / kaxial
        assert np.allclose(solver.xg[6], daxial), f'dx={daxial} xg[6]={solver.xg[6]}'
        assert np.allclose(solver.xg[9], dtorsion), f'dx={dtorsion} xg[6]={solver.xg[9]}'
        assert np.allclose(solver.Fg[6], mag_axial), f'F={mag_axial} Fg[6]={solver.Fg[6]}'
        assert np.allclose(solver.Fg[9], mag_torsion), f'F={mag_torsion} Fg[9]={solver.Fg[9]}'


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
        model.add_mat1(mid, E, G, nu, rho=0.0, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
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
        assert dx == solver.xg[6*2], f'dx={dx} xg={solver.xg}'


#Class TestStaticAero(unittest.TestCase):
#    def test_bwb_saero_modes(self):
#        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero.bdf'
#        model = read_bdf(bdf_filename)
#        model.add_eigrl(31, nd=15)
#        model.setup()
#
#        model.sol = 103
#        subcases = model.subcases
#        subcase: Subcase = subcases[1]
#        subcase.add_integer_type('METHOD', 31)
#        solver = Solver(model)
#        solver.run()

#class TestStatic2(unittest.TestCase):
#    def test_static(self):
#        bdf_filename = MODEL_PATH / 'sol_101_elements' / #'buckling_solid_shell_bar.bdf'
#        model = read_bdf(bdf_filename)
#        solver = Solver(model)
#        solver.run()

class TestStaticSolid(unittest.TestCase):
    def test_ctetra_10(self):
        bdf_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.bdf'
        model = read_bdf(bdf_filename)
        solver = Solver(model)
        solver.run()

    def test_ctetra4(self):
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'ctetra4.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 0., 1.])

        mid = 3
        pid = 4
        model.add_ctetra(10, pid, [1, 2, 3, 4])
        model.add_psolid(pid, mid)

        spc_id = 3
        load_id = 2
        model.add_spc1(spc_id, '123456', [1, 2, 3], comment='')
        model.add_force(load_id, 4, 1.0, [0., 0., 1.], cid=0, comment='')

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        model.add_param('GRDPNT', 0)

        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_hexa8(self):
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'chexa8.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        model.add_grid(5, [0., 0., 1.])
        model.add_grid(6, [1., 0., 1.])
        model.add_grid(7, [1., 1., 1.])
        model.add_grid(8, [0., 1., 1.])

        mid = 3
        pid = 4
        model.add_chexa(10, pid, [1, 2, 3, 4, 5, 6, 7, 8])
        model.add_psolid(pid, mid)

        spc_id = 3
        load_id = 2
        model.add_spc1(spc_id, '123456', [1, 2, 3], comment='')
        model.add_force(load_id, 4, 1.0, [0., 0., 1.], cid=0, comment='')

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        model.add_param('GRDPNT', 0)

        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_penta6(self):
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cpenta6.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])

        model.add_grid(4, [0., 0., 1.])
        model.add_grid(5, [1., 0., 1.])
        model.add_grid(6, [1., 1., 1.])

        mid = 3
        pid = 4
        model.add_cpenta(10, pid, [1, 2, 3, 4, 5, 6])
        model.add_psolid(pid, mid)

        spc_id = 3
        load_id = 2
        model.add_spc1(spc_id, '123456', [1, 2, 3], comment='')
        model.add_force(load_id, 4, 1.0, [0., 0., 1.], cid=0, comment='')

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        model.add_param('GRDPNT', 0)
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()


class TestStaticShear(unittest.TestCase):
    def test_cshear1(self):
        """Tests CSHEAR."""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cshear1.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 0., 2.])
        model.add_grid(4, [0., 0., 2.])
        thickness = 0.3
        mid = 3
        nids = [1, 2, 3, 4]
        model.add_cshear(1, 1, nids)
        model.add_cshear(2, 1, nids)
        model.add_pshear(1, mid, thickness)

        model.add_conrod(10, mid, [1, 2], A=1.0)
        model.add_conrod(11, mid, [2, 3], A=1.0)
        model.add_conrod(12, mid, [3, 4], A=1.0)
        model.add_conrod(13, mid, [4, 1], A=1.0)

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        spc_id = 3
        load_id = 2
        model.add_spc1(spc_id, '123456', [1, 2], comment='')
        model.add_force(load_id, 3, 1.0, [0., 0., 1.], cid=0, comment='')
        model.add_force(load_id, 4, 1.0, [0., 0., 1.], cid=0, comment='')

        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()


class TestStaticShell(unittest.TestCase):
    """tests the shells"""
    def test_cquad4_bad_normal(self):
        """test that the code crashes with a bad normal"""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad4_bad_normal.bdf'
        mid = 3
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        model.add_grid(2, [1., 0., 2.])
        model.add_grid(4, [0., 0., 2.])
        nids = [1, 2, 3, 4]
        model.add_cquad4(5, 5, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_pcomp(5, [mid], [0.1], thetas=None,
                        souts=None, nsm=0., sb=0., ft=None, tref=0., ge=0., lam=None, z0=None, comment='')

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        spc_id = 3
        load_id = 2
        components = '123456'
        nodes = [1, 2]
        model.add_spc1(spc_id, components, nodes, comment='')
        model.add_force(load_id, 3, 1.0, [0., 0., 1.], cid=0, comment='')
        model.add_force(load_id, 4, 1.0, [0., 0., 1.], cid=0, comment='')

        setup_static_case_control(model)
        solver = Solver(model)

        with self.assertRaises(RuntimeError):
            # invalid normal vector
            solver.run()
        #os.remove(model.bdf_filename)

    def test_cquad4_bad_jacobian(self):
        """
        Tests that the code crashes with a degenerate CQUAD4
        (crossed/bowtie element with negative Jacobian).
        """
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad4_bad_jacobian.bdf'
        mid = 3
        # Bowtie element: nodes 2 and 4 are swapped to create negative Jacobian
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0., 1., 0.])
        model.add_grid(4, [1., 1., 0.])
        nids = [1, 2, 3, 4]
        model.add_cquad4(5, 5, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_pcomp(5, [mid], [0.1], thetas=None,
                        souts=None, nsm=0., sb=0., ft=None, tref=0., ge=0., lam=None, z0=None, comment='')

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        spc_id = 3
        load_id = 2
        components = '123456'
        nodes = [1, 2]
        model.add_spc1(spc_id, components, nodes, comment='')
        model.add_force(load_id, 3, 1.0, [0., 0., 1.], cid=0, comment='')
        model.add_force(load_id, 4, 1.0, [0., 0., 1.], cid=0, comment='')

        setup_static_case_control(model)
        solver = Solver(model)
        with self.assertRaises(RuntimeError):
            solver.run()
        #os.remove(model.bdf_filename)

    def test_ctria3_pshell_mat1(self):
        """Tests CQUAD4/PSHELL/MAT1 with CROD, CONROD, and CONM2."""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'ctria3_pshell_mat1.bdf'
        # 6---5---4
        # |       |
        # |       |
        # 1---2---3
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        thickness = 0.3
        mid = 3
        E = 1.0E7
        G = None
        nu = 0.3

        # --- CQUAD4 elements ---
        model.add_ctria3(1, 1, [1, 2, 3], theta_mcid=0.0, zoffset=0.,
                         tflag=0, T1=None, T2=None, T3=None)
        model.add_ctria3(2, 1, [1, 2, 3], theta_mcid=0.0, zoffset=0.,
                         tflag=0, T1=None, T2=None, T3=None)
        model.add_pshell(1, mid1=mid, t=thickness, mid2=mid, twelveIt3=1.0,
                         mid3=mid, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None)

        # --- CROD element (node 5 to node 6, stiffening bar) ---
        model.add_crod(10, 10, [1, 2])
        model.add_prod(10, mid, A=0.5, j=0.01)

        # --- CONROD element (node 4 to node 6, diagonal brace) ---
        model.add_conrod(20, mid, [2, 3], A=0.3, j=0.005)

        # --- CONM2 (point mass at node 5) ---
        model.add_conm2(30, 3, mass=2.0, cid=0, X=None, I=None)

        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0)

        # BC/loads: fix nodes 1, 2, 3; load nodes 4, 5, 6
        spc_id = 3
        load_id = 2
        model.add_spc1(spc_id, '123456', [1, 2])
        model.add_force(load_id, 1, 1.0, [0., 0., 1.], cid=0)
        model.add_force(load_id, 2, 1.0, [0., 0., 1.], cid=0)
        model.add_force(load_id, 3, 1.0, [0., 0., 1.], cid=0)
        
        model.add_pload4(load_id, eids=1, pressures=2.0)

        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_cquad4_pshell_mat1(self):
        """Tests CQUAD4/PSHELL/MAT1 with CROD, CONROD, and CONM2."""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad4_pshell_mat1.bdf'
        # 6---5---4
        # |       |
        # |       |
        # 1---2---3
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [2., 0., 0.])
        model.add_grid(4, [2., 0., 2.])
        model.add_grid(5, [1., 0., 2.])
        model.add_grid(6, [0., 0., 2.])
        thickness = 0.3
        mid = 3
        E = 1.0E7
        G = None
        nu = 0.3

        # --- CQUAD4 elements ---
        model.add_cquad4(1, 1, [1, 2, 5, 6], theta_mcid=0.0, zoffset=0.,
                         tflag=0, T1=None, T2=None, T3=None, T4=None)
        model.add_cquad4(2, 1, [2, 3, 4, 5], theta_mcid=0.0, zoffset=0.,
                         tflag=0, T1=None, T2=None, T3=None, T4=None)
        model.add_pshell(1, mid1=mid, t=thickness, mid2=mid, twelveIt3=1.0,
                         mid3=mid, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None)

        # --- CROD element (node 5 to node 6, stiffening bar) ---
        model.add_crod(10, 10, [5, 6])
        model.add_prod(10, mid, A=0.5, j=0.01)

        # --- CONROD element (node 4 to node 6, diagonal brace) ---
        model.add_conrod(20, mid, [4, 6], A=0.3, j=0.005)

        # --- CONM2 (point mass at node 5) ---
        model.add_conm2(30, 5, mass=2.0, cid=0, X=None, I=None)

        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0)

        # BC/loads: fix nodes 1, 2, 3; load nodes 4, 5, 6
        spc_id = 3
        load_id = 2
        model.add_spc1(spc_id, '123456', [1, 2, 3])
        model.add_force(load_id, 4, 1.0, [0., 0., 1.], cid=0)
        model.add_force(load_id, 5, 1.0, [0., 0., 1.], cid=0)
        model.add_force(load_id, 6, 1.0, [0., 0., 1.], cid=0)
        
        model.add_pload4(load_id, eids=1, pressures=2.0)

        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_cquad4_macneal(self):
        """Tests CQUAD4 with PARAM,MYQUAD,MACN (MacNeal formulation)."""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad4_macneal.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 0., 2.])
        model.add_grid(4, [0., 0., 2.])
        thickness = 0.3
        mid = 3
        nids = [1, 2, 3, 4]
        model.add_cquad4(1, 1, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_pshell(1, mid1=mid, t=thickness, mid2=mid, twelveIt3=1.0,
                         mid3=mid, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None, comment='')
        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        spc_id = 3
        load_id = 2
        model.add_spc1(spc_id, '123456', [1, 2], comment='')
        model.add_force(load_id, 3, 1.0, [0., 0., 1.], cid=0, comment='')
        model.add_force(load_id, 4, 1.0, [0., 0., 1.], cid=0, comment='')

        # Set PARAM,MYQUAD,MACN
        model.params['MYQUAD'] = PARAM('MYQUAD', ['MACN'])

        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_cquad4_allman_drilling(self):
        """Tests CQUAD4 with PARAM,MYQDRIL,ALLMAN."""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad4_allman.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 0., 2.])
        model.add_grid(4, [0., 0., 2.])
        mid = 3
        nids = [1, 2, 3, 4]
        model.add_cquad4(1, 1, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None)
        model.add_pshell(1, mid1=mid, t=0.3, mid2=mid, twelveIt3=1.0,
                         mid3=mid, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None)
        model.add_mat1(mid, 1.0E7, None, 0.3, rho=1.0)
        model.add_spc1(3, '123456', [1, 2])
        model.add_force(2, 3, 1.0, [0., 0., 1.], cid=0)
        model.add_force(2, 4, 1.0, [0., 0., 1.], cid=0)
        model.params['MYQDRIL'] = PARAM('MYQDRIL', ['ALLMAN'])
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_cquad4_hb_drilling(self):
        """Tests CQUAD4 with PARAM,MYQDRIL,HB (Hughes-Brezzi)."""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad4_hb.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 0., 2.])
        model.add_grid(4, [0., 0., 2.])
        mid = 3
        nids = [1, 2, 3, 4]
        model.add_cquad4(1, 1, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None)
        model.add_pshell(1, mid1=mid, t=0.3, mid2=mid, twelveIt3=1.0,
                         mid3=mid, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None)
        model.add_mat1(mid, 1.0E7, None, 0.3, rho=1.0)
        model.add_spc1(3, '123456', [1, 2])
        model.add_force(2, 3, 1.0, [0., 0., 1.], cid=0)
        model.add_force(2, 4, 1.0, [0., 0., 1.], cid=0)
        model.params['MYQDRIL'] = PARAM('MYQDRIL', ['HB'])
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_cquad4_pcomp_mat8(self):
        """Tests CQUAD4/PCOMP/MAT8 (orthotropic composite) with all quad types."""
        for quad_type in ['MITC4', 'MACN', 'MACN2']:
            model = BDF(debug=None, log=None, mode='msc')
            model.bdf_filename = TEST_DIR / f'cquad4_pcomp_mat8_{quad_type}.bdf'
            model.add_grid(1, [0., 0., 0.])
            model.add_grid(2, [1., 0., 0.])
            model.add_grid(3, [1., 0., 2.])
            model.add_grid(4, [0., 0., 2.])
            mid = 3
            nids = [1, 2, 3, 4]
            model.add_cquad4(1, 1, nids, theta_mcid=0.0, zoffset=0.,
                             tflag=0, T1=None, T2=None, T3=None, T4=None)
            # 4-ply quasi-isotropic: [0/45/-45/90]
            model.add_pcomp(1, [mid, mid, mid, mid],
                            [0.075, 0.075, 0.075, 0.075],
                            thetas=[0., 45., -45., 90.])
            model.add_mat8(mid, e11=1.4e7, e22=1.0e6, nu12=0.3,
                           g12=5.0e5, g1z=5.0e5, g2z=3.0e5, rho=1.0)
            model.add_spc1(3, '123456', [1, 2])
            model.add_force(2, 3, 1.0, [0., 0., 1.], cid=0)
            model.add_force(2, 4, 1.0, [0., 0., 1.], cid=0)
            model.params['MYQUAD'] = PARAM('MYQUAD', [quad_type])
            setup_static_case_control(model)
            solver = Solver(model)
            solver.run()

    def test_cquad4_pshell_mat8(self):
        """Tests CQUAD4/PSHELL/MAT8 (orthotropic) with all quad types."""
        for quad_type in ['MITC4', 'MACN', 'MACN2']:
            model = BDF(debug=None, log=None, mode='msc')
            model.bdf_filename = TEST_DIR / f'cquad4_mat8_{quad_type}.bdf'
            model.add_grid(1, [0., 0., 0.])
            model.add_grid(2, [1., 0., 0.])
            model.add_grid(3, [1., 0., 2.])
            model.add_grid(4, [0., 0., 2.])
            mid = 3
            nids = [1, 2, 3, 4]
            model.add_cquad4(1, 1, nids, theta_mcid=0.0, zoffset=0.,
                             tflag=0, T1=None, T2=None, T3=None, T4=None)
            model.add_pshell(1, mid1=mid, t=0.3, mid2=mid, twelveIt3=1.0,
                             mid3=mid, tst=0.833333, nsm=0.0, z1=None,
                             z2=None, mid4=None)
            # Orthotropic: carbon/epoxy-like
            model.add_mat8(mid, e11=1.4e7, e22=1.0e6, nu12=0.3,
                           g12=5.0e5, g1z=5.0e5, g2z=3.0e5, rho=1.0)
            model.add_spc1(3, '123456', [1, 2])
            model.add_force(2, 3, 1.0, [0., 0., 1.], cid=0)
            model.add_force(2, 4, 1.0, [0., 0., 1.], cid=0)
            model.params['MYQUAD'] = PARAM('MYQUAD', [quad_type])
            setup_static_case_control(model)
            solver = Solver(model)
            solver.run()

    def _test_cquad8_pshell_mat1(self):
        """Tests a CQUAD8/PSHELL/MAT1"""
        # 4--7--3
        # |     |
        # 8     6
        # |     |
        # 1--5--2
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad8_pshell_mat1.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 0., 2.])
        model.add_grid(4, [0., 0., 2.])

        model.add_grid(5, [0.5, 0., 0.])
        model.add_grid(6, [1., 0., 1.])
        model.add_grid(7, [0.5, 0., 2.])
        model.add_grid(8, [0., 0., 1.])

        mid = 3
        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        model.add_cquad8(1, 1, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_cquad8(2, 2, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_cquad8(3, 3, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_cquad8(4, 4, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_cquad8(5, 5, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')

        model.add_pshell(1, mid1=mid, t=0.3, mid2=None, twelveIt3=1.0,
                         mid3=None, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None, comment='')
        model.add_pshell(2, mid1=None, t=0.3, mid2=mid, twelveIt3=1.0,
                         mid3=None, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None, comment='')
        model.add_pshell(3, mid1=None, t=0.3, mid2=None, twelveIt3=1.0,
                         mid3=mid, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None, comment='')
        model.add_pshell(4, mid1=None, t=0.3, mid2=None, twelveIt3=1.0,
                         mid3=None, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=mid, comment='')
        model.add_pcomp(5, [mid], [0.1], thetas=None,
                        souts=None, nsm=0., sb=0., ft=None, tref=0., ge=0., lam=None, z0=None, comment='')

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        spc_id = 3
        load_id = 2
        components = '123456'
        nodes = [1, 2, 5]
        model.add_spc1(spc_id, components, nodes, comment='')
        model.add_force(load_id, 3, 1.0, [0., 0., 1.], cid=0, comment='')
        model.add_force(load_id, 4, 1.0, [0., 0., 1.], cid=0, comment='')

        setup_static_case_control(model)
        solver = Solver(model)
        with self.assertRaises(RuntimeError):
            solver.run()
        #os.remove(model.bdf_filename)
        #os.remove(solver.f06_filename)
        #os.remove(solver.op2_filename)


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

if __name__ == '__main__':  # pragma: no cover
    np.seterr(all='raise')
    unittest.main()
