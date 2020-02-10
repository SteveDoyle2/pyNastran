import unittest
import numpy as np
from pyNastran.dev.solver.solver import Solver, BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck


class TestSolverSpring(unittest.TestCase):
    def test_celas1(self):
        """Tests a CELAS1/PELAS"""
        model = BDF(debug=True, log=None, mode='msc')
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
        setup_case_control(model)

        solver = Solver(model)
        model.sol = 103
        solver.run()

        model.sol = 101
        solver.run()

        # F = k * d
        # d = F / k
        # 20=1000.*d
        d = mag / k
        assert np.allclose(solver.xa_[0], d)

    def test_celas2_cd(self):
        """Tests a CELAS2"""
        model = BDF(debug=True, log=None, mode='msc')
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
        setup_case_control(model)

        solver = Solver(model)
        solver.run()

        # F = k * d
        d = mag / k
        assert np.allclose(solver.xa_[0], d)


class TestSolverRod(unittest.TestCase):
    """tests the rods"""
    def test_crod_axial(self):
        """Tests a CROD/PROD"""
        model = BDF(debug=True, log=None, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0, St=0.0,
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
        setup_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_crod_torsion(self):
        """Tests a CROD/PROD"""
        model = BDF(debug=True, log=None, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_crod(eid, pid, nids)
        model.add_prod(pid, mid, A=0.0, j=2., c=0., nsm=0.)

        load_id = 2
        spc_id = 3
        nid = 2
        mag = 1.
        fxyz = [1., 0., 0.]
        model.add_moment(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_crod(self):
        """Tests a CROD/PROD"""
        model = BDF(debug=True, log=None, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0, St=0.0,
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
        setup_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_crod_aset(self):
        """
        Tests a CROD/PROD using an ASET

        same answer as ``test_crod``
        """
        model = BDF(debug=True, log=None, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0, St=0.0,
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

        setup_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_crod_mpc(self):
        """Tests a CROD/PROD"""
        model = BDF(debug=True, log=None, mode='msc')
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
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0, St=0.0,
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
        setup_case_control(model, extra_case_lines=['MPC=10'])
        solver = Solver(model)
        solver.run()

    def test_ctube(self):
        """Tests a CTUBE/PTUBE"""
        model = BDF(debug=True, log=None, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_ctube(eid, pid, nids)
        OD1 = 1.0
        model.add_ptube(pid, mid, OD1, t=0.1, nsm=0., OD2=None, comment='')

        load_id = 2
        spc_id = 3
        nid = 2
        mag = 1.
        fxyz = [1., 0., 0.]
        model.add_force(load_id, nid, mag, fxyz, cid=0)

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_case_control(model)
        solver = Solver(model)
        solver.run()

        # F = k * x
        # x = F / k

    def test_conrod(self):
        """Tests a CONROD"""
        model = BDF(debug=True, log=None, mode='msc')
        L = 1.
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [L, 0., 0.])

        nids = [1, 2]
        eid = 1
        mid = 3
        E = 200.
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0, St=0.0,
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
        setup_case_control(model)
        solver = Solver(model)
        solver.run()

        # TODO: why is this a torsional result???
        G = E / (2 * (1 + nu))
        kaxial = A * E / L
        ktorsion = G * J / L
        # F = k * d
        daxial = mag_axial / kaxial
        dtorsion = mag_torsion / ktorsion
        assert np.allclose(solver.xa_[0], daxial), f'daxial={daxial} kaxial={kaxial} xa_={solver.xa_}'
        assert np.allclose(solver.xa_[3], dtorsion), f'dtorsion={dtorsion} ktorsion={ktorsion} xa_={solver.xa_}'

    def test_cbar(self):
        """Tests a CBAR/PBAR"""
        model = BDF(debug=True, log=None, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)

        x = [0., 1., 0.]
        g0 = None
        model.add_cbar(eid, pid, nids, x, g0,
                       offt='GGG', pa=0, pb=0, wa=None, wb=None, comment='')
        model.add_pbar(pid, mid,
                       A=1., i1=1., i2=1., i12=1., j=1.,
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
        setup_case_control(model)
        solver = Solver(model)
        solver.run()

def setup_case_control(model, extra_case_lines=None):
    lines = [
        'STRESS(PLOT,PRINT) = ALL',
        'STRAIN(PLOT,PRINT) = ALL',
        'FORCE(PLOT,PRINT) = ALL',
        'DISP(PLOT,PRINT) = ALL',
        'GPFORCE(PLOT,PRINT) = ALL',
        'SPCFORCE(PLOT,PRINT) = ALL',
        'MPCFORCE(PLOT,PRINT) = ALL',
        'OLOAD(PLOT,PRINT) = ALL',
        'SUBCASE 1',
        '  LOAD = 2',
        '  SPC = 3',
    ]
    if extra_case_lines is not None:
        lines += extra_case_lines
    cc = CaseControlDeck(lines, log=model.log)
    model.sol = 101
    model.case_control_deck = cc

if __name__ == '__main__':
    unittest.main()
