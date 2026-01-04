import os
import pathlib
import unittest
import numpy as np
from cpylog import SimpleLogger
import pyNastran
from pyNastran.dev.bdf_vectorized3.solver.solver import Solver, BDF, partition_vector2
from pyNastran.dev.solver.solver import Solver as SolverOld, BDF as BDFold
#from .solver import Solver, BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck

PKG_PATH = pathlib.Path(pyNastran.__path__[0])
TEST_DIR = PKG_PATH / 'dev' / 'solver'


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


def setup_modal_case_control(model: BDF, extra_case_lines=None):
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
        #'  FREQUENCY = 100',
    ]
    if extra_case_lines is not None:
        lines += extra_case_lines
    cc = CaseControlDeck(lines, log=model.log)
    model.sol = 103
    model.case_control_deck = cc

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
        fxyz = [1., 0., 0.]
        mag = 20.
        model.add_force(load_id, 2, mag, fxyz, cid=0, comment='')

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)

        solver = Solver(model)
        model.sol = 101
        solver.run()

    def test_celas1(self):
        """Tests a CELAS1/PELAS"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'celas1.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        model.add_celas1(eid, pid, nids, c1=1, c2=1, comment='')
        k = 1000.

        load_id = 2
        spc_id = 3
        model.add_pelas(pid, k, ge=0., s=0., comment='')
        # model.add_sload(load_id, 2, 20.)
        fxyz = [1., 0., 0.]
        mag = 20.
        model.add_force(load_id, 2, mag, fxyz, cid=0, comment='')

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)

        solver = Solver(model)
        model.sol = 101
        solver.run()

        eid = 2
        nid = 2
        mass = 1.
        model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')
        setup_modal_case_control(model)
        sid = 103
        nmodes = 2
        model.add_eigrl(sid, v1=None, v2=None, nd=nmodes, msglvl=0,
                        maxset=None, shfscl=None, norm=None, options=None, values=None, comment='')
        model.sol = 103
        solver.run()

        # F = k * d
        # d = F / k
        # 20=1000.*d
        d = mag / k
        assert np.allclose(solver.xa_[0], d)
        os.remove(solver.f06_filename)
        os.remove(solver.op2_filename)

    def test_celas2_cd(self):
        """Tests a CELAS2"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model_old = BDFold(log=log, mode='msc')
        modelv = BDF(log=log, mode='msc')
        models = [model_old, modelv]
        for model in models:
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
            fxyz = [1., 0., 0.]
            mag = 20.
            model.add_force(load_id, 2, mag, fxyz, cid=0, comment='')

            components = 123456
            nodes = 1
            model.add_spc1(spc_id, components, nodes, comment='')
            setup_static_case_control(model)

        solver_old = SolverOld(model_old)
        solver_old.run()

        solverv = Solver(modelv)
        solverv.run()

        # F = k * d
        d = mag / k
        Fg1 = solver_old.Fg
        Fg2 = solverv.Fg

        Kgg1 = solver_old.Kgg
        Kgg2 = solverv.Kgg

        Kaa1 = solver_old.Kaa.toarray()
        Kaa2 = solverv.Kaa.toarray()
        Fg_expected = [ 0., 0.0, 0., 0., 0., 0.,
                        20., 0., 0., 0., 0., 0.]
        assert np.allclose(Fg2, Fg_expected), f'Force Error:\n Fg2={Fg2}\n Fg_expected={Fg_expected}'
        assert np.allclose(Kgg1, Kgg2)
        assert np.allclose(Kaa1, Kaa2)
        assert np.allclose(solver_old.xa_[0], d)
        assert np.allclose(solverv.xa_[0], d)

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
        #fxyz = [1., 0., 0.]
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
        #fxyz = [1., 0., 0.]
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


class TestStaticRod(unittest.TestCase):
    """tests the rods"""
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
        fxyz = [1., 0., 0.]
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
        fxyz = [1., 0., 0.]
        model.add_moment(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_crod_spcd(self):
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
        fxyz = [0., 0., 0.]
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
        fxyz = [1., 0., 0.]
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
        Kgg = solver.Kgg[dof, :][:, dof]
        assert np.allclose(Kgg, Kgg_expectedi)


    def test_crod_rotate(self):
        """Tests a CROD/PROD"""
        log = SimpleLogger(level='warning', encoding='utf-8')
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
        fxyz = [x, y, 0.]
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
        model.add_mat1(mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_crod(eid, pid, nids)
        model.add_prod(pid, mid, A=1.0, j=2., c=0., nsm=0.)

        load_id = 2
        spc_id = 3
        nid = 2
        mag = 1.
        fxyz = [1., 0., 0.]
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
        """Tests a CROD/PROD"""
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
        fxyz = [0., 0., 1.]
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_static_case_control(model, extra_case_lines=['MPC=10'])
        solver = Solver(model)
        solver.run()

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
        fxyz = [1., 0., 0.]
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
        fxyz = [1., 0., 0.]
        model.add_force(load_id, nid, mag_axial, fxyz, cid=0)

        mag_torsion = 20000.
        fxyz = [1., 0., 0.]
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


def mesh_line_by_points(xyz1: np.ndarray, xyz2: np.ndarray, nelement: int,
                        nid0: int=1, eid0: int=1,
                        idtype='int32',
                        fdtype='float64') -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    assert nelement >= 1
    dxyz = xyz2 - xyz1
    length = np.linalg.norm(dxyz)
    assert length > 0, length
    nnodes = nelement + 1
    nids = np.arange(0, nnodes)
    eids = np.arange(eid0, eid0+nelement)
    assert len(nids) == nnodes, nnodes
    assert len(eids) == nelement, nelement

    line_nodes = nid0 + np.column_stack([nids[:-1], nids[1:]])
    xyzs = (nids / length)[:, np.newaxis] * dxyz
    nids += nid0
    return nids, xyzs, eids, line_nodes


class TestModalBar(unittest.TestCase):
    """tests the CBARs"""
    def test_cbar_modes(self):
        """Tests a CBAR/PBAR"""
        log = SimpleLogger(level='debug', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'cbar_modes.bdf'
        xyz1 = np.array([0., 0., 0.])
        xyz2 = np.array([1., 0., 0.])
        nelement = 100
        nmodes = 16

        nid0 = 1
        eid0 = 1
        pid = 10
        mid = 100
        node_id, xyz, element_id, line_nodes = mesh_line_by_points(xyz1, xyz2, nelement, nid0=nid0, eid0=eid0)

        nnode = len(node_id)
        zero = np.zeros(nnode, dtype='int32')
        cp = zero
        cd = zero
        ps = zero
        seid = zero
        model.grid._save(node_id, cp, cd, xyz, ps, seid, comment=None)
        model.grid._xyz_cid0 = xyz

        one = np.ones(nelement, dtype='int32')
        property_id = pid * one

        # y-axis
        x = np.zeros((nelement, 3), dtype='float64')
        x[:, 1] = 1.
        model.cbar._save(element_id, property_id, line_nodes, x=x)

        #bar_type = 'BAR'
        #dim = [1., 2.]
        #model.add_pbarl(pid, mid, bar_type, dim, group='MSCBML0', nsm=0., comment='')

        A = 1.0
        b = 3.0
        h = 3.0
        i1 = 1 / 12 * b * h ** 3
        i2 = 1 / 12 * h * b ** 3
        i12 = 0.
        j = i1 + i2
        model.add_pbar(pid, mid, area=A, i1=i1, i2=i2, i12=i12, j=j, nsm=0.,
                       c1=0., c2=0., d1=0., d2=0., e1=0., e2=0., f1=0., f2=0.,
                       k1=1.e8, k2=1.e8, comment='')

        E = 3.0e7
        G = None
        nu = 0.3
        rho = 0.1
        model.add_mat1(mid, E, G, nu, rho=rho,
                       alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)

        model.add_eigrl(103, nd=nmodes)
        setup_modal_case_control(model)
        model.setup()
        assert len(model.cbar)
        solver = Solver(model)
        solver.run()
        x = 1


class TestStaticBar(unittest.TestCase):
    """tests the CBARs"""
    def test_cbar_rbe3_ex(self):
        """Tests a CBAR/PBAR"""
        log = SimpleLogger(level='debug', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'cbar.bdf'
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
        fxyz = [0., 0., 1.]
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
        log = SimpleLogger(level='debug', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'cbar.bdf'
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
        nid = 2
        mag = 1.
        fxyz = [1., 0., 0.]
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
        assert dx == solver.xg[6], f'dx={dx} xg={xg}'

    def test_cbar2(self):
        """Tests a CBAR/PBAR"""
        log = SimpleLogger(level='debug', encoding='utf-8')
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
        fxyz = [1., 0., 0.]
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
        assert dx == solver.xg[6*2], f'dx={dx} xg={xg}'

    def test_cbeam(self):
        """Tests a CBEAM/PBEAM"""
        log = SimpleLogger(level='debug', encoding='utf-8')
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
        fxyz = [1., 0., 0.]
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
        assert dx == solver.xg[6], f'dx={dx} xg={xg}'


    def test_cbeam2(self):
        """Tests a CBEAM/PBEAM"""
        model = BDF(debug=True, log=None, mode='msc')
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
        fxyz = [1., 0., 0.]
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


class TestHarmonic(unittest.TestCase):
    def test_spring_mass_damper_direct(self):
        """spring-mass-damper problem"""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad4_bad_normal.bdf'
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

        setup_direct_frequency_case_control(model)
        solver = Solver(model)
        #with self.assertRaises(RuntimeError):
        solver.run()

    def test_spring_mass_damper_freq_response(self):
        """spring-mass-damper problem"""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad4_bad_normal.bdf'
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
        Tests that the code crashes with a really terrible CQUAD4

        The Jacobian is defined between [-1, 1]
        """
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad4_bad_jacobian.bdf'
        mid = 3
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0.5, 100., 0.])
        model.add_grid(3, [1., 0., 0.])
        model.add_grid(4, [0.5, 1., 0.])
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

    def _test_cquad4_pshell_mat1(self):
        """Tests a CQUAD4/PSHELL/MAT1"""
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cquad4_pshell_mat1.bdf'
        model.add_grid(1, [0., 0., 0., ])
        model.add_grid(2, [1., 0., 0., ])
        model.add_grid(3, [1., 0., 2., ])
        model.add_grid(4, [0., 0., 2., ])
        area = 2.0
        nsm_area = 0.0
        thickness = 0.3
        rho = 0.1
        mass = area * (thickness * rho + nsm_area)
        mid = 3
        nids = [1, 2, 3, 4]
        model.add_cquad4(1, 1, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_cquad4(2, 2, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_cquad4(3, 3, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_cquad4(4, 4, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_cquad4(5, 5, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')

        model.add_pshell(1, mid1=mid, t=thickness, mid2=None, twelveIt3=1.0,
                         mid3=None, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None, comment='')
        model.add_pshell(2, mid1=None, t=thickness, mid2=mid, twelveIt3=1.0,
                         mid3=None, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None, comment='')
        model.add_pshell(3, mid1=None, t=thickness, mid2=None, twelveIt3=1.0,
                         mid3=mid, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=None, comment='')
        model.add_pshell(4, mid1=None, t=thickness, mid2=None, twelveIt3=1.0,
                         mid3=None, tst=0.833333, nsm=0.0, z1=None, z2=None,
                         mid4=mid, comment='')
        model.add_pcomp(5, [mid], [thickness], thetas=None,
                        souts=None, nsm=0., sb=0., ft=None, tref=0., ge=0., lam=None, z0=None, comment='')

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
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
        #with self.assertRaises(RuntimeError):
        solver.run()
        #print('total_mass =', mass * 5)
        #os.remove(model.bdf_filename)
        #os.remove(solver.f06_filename)
        #os.remove(solver.op2_filename)

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


if __name__ == '__main__':   # pragma: no cover
    unittest.main()
