import unittest
from pyNastran.dev.solver.solver import Solver, BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck


class TestSolverSprings(unittest.TestCase):
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
        model.add_force(load_id, 2, 20., fxyz, cid=0, comment='')

        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')
        setup_case_control(model)

        solver = Solver(model)
        solver.run()

    def test_celas2(self):
        """Tests a CELAS2"""
        model = BDF(debug=True, log=None, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        nids = [1, 2]
        eid = 1
        k = 1000.
        model.add_celas2(eid, k, nids, c1=1, c2=1, ge=0., s=0., comment='')

        load_id = 2
        spc_id = 3
        # model.add_sload(load_id, 2, 20.)
        fxyz = [1., 0., 0.]
        model.add_force(load_id, 2, 20., fxyz, cid=0, comment='')

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

    def test_conrod(self):
        """Tests a CONROD"""
        model = BDF(debug=True, log=None, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        nids = [1, 2]
        eid = 1
        #pid = 2
        mid = 3
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        model.add_conrod(eid, mid, nids, A=1.0, j=2.0, c=0.0, nsm=0.0)

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

def setup_case_control(model):
    lines = [
        'STRESS(PLOT,PRINT) = ALL',
        'STRAIN(PLOT,PRINT) = ALL',
        'FORCE(PLOT,PRINT) = ALL',
        'DISP(PLOT,PRINT) = ALL',
        'GPFORCE(PLOT,PRINT) = ALL',
        'SPCFORCE(PLOT,PRINT) = ALL',
        'SUBCASE 1',
        '  LOAD = 2',
        '  SPC = 3',
    ]
    cc = CaseControlDeck(lines, log=model.log)
    model.sol = 101
    model.case_control_deck = cc

if __name__ == '__main__':
    unittest.main()
