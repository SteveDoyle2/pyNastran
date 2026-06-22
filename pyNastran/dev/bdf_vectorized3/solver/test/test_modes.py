import os
from pathlib import Path
import unittest
import numpy as np
from cpylog import SimpleLogger
import pyNastran
from pyNastran.dev.bdf_vectorized3.solver.solver import Solver, partition_vector2
from pyNastran.dev.bdf_vectorized3.bdf import BDF, read_bdf

#from pyNastran.bdf.cards.params import PARAM
from pyNastran.f06.errors import FatalError
from pyNastran.bdf.case_control_deck import CaseControlDeck, Subcase
from pyNastran.dev.bdf_vectorized3.solver.test.test_vector_solver_springs import setup_static_case_control, MODEL_PATH

PKG_PATH = Path(pyNastran.__path__[0])
TEST_DIR = Path(__file__).parent / 'examples'


class TestModes(unittest.TestCase):

    def test_celas1_modes(self):
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
        #model.add_sload(load_id, 2, 20.)
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

        eid = 2
        nid = 2
        mass = 1.
        model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')
        os.remove(solver.f06_filename)
        os.remove(solver.op2_filename)

        model.bdf_filename = TEST_DIR / 'celas1_modes.bdf'
        setup_modal_case_control(model)
        sid = 103
        nmodes = 2
        model.add_eigrl(sid, v1=None, v2=None, nd=nmodes, msglvl=0,
                        maxset=None, shfscl=None, norm=None, options=None, values=None,comment='')
        model.sol = 103
        solver.run()

        # F = k * d
        # d = F / k
        # 20=1000.*d
        d = mag / k
        assert np.allclose(solver.xa_[0], d)
        os.remove(solver.f06_filename)
        os.remove(solver.op2_filename)

    def test_rod_modes(self):
        """Tests a CONROD"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(log=log, mode='msc')
        model.bdf_filename = TEST_DIR / 'rods_modes.bdf'
        L = 1.
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [L, 0., 0.])

        nids = [1, 2]
        eid = 10
        pid = 20
        mid = 3
        E = 200.
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1, alpha=0.0, tref=0.0, ge=0.0, St=0.0,
                       Sc=0.0, Ss=0.0, mcsid=0)
        A = 1.0
        J = 2.0
        model.add_conrod(eid, mid, nids, A=A, j=J, c=0.0, nsm=0.0)

        OD1 = 1.0
        t = 0.1
        model.add_ctube(eid+1, pid+1, nids)
        model.add_ptube(pid+1, mid, OD1, t=t, nsm=0., OD2=None, comment='')

        # --- CROD element (node 5 to node 6, stiffening bar) ---
        model.add_crod(eid+2, pid+2, nids)
        model.add_prod(pid+2, mid, A=0.5, j=0.01)

        spc_id = 3
        components = 123456
        nodes = 1
        model.add_spc1(spc_id, components, nodes, comment='')

        #nodes = 2
        #model.add_spc1(spc_id, 23456, nodes, comment='')
        setup_modal_case_control(model, nmodes=6)
        solver = Solver(model)
        solver.run()
        os.remove(solver.f06_filename)
        os.remove(solver.op2_filename)

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
        node_id, xyz, element_id, line_nodes = mesh_line_by_points(
            xyz1, xyz2, nelement, nid0=nid0, eid0=eid0)

        nnode = len(node_id)
        zero = np.zeros(nnode, dtype='int32')
        cp = zero
        cd = zero
        ps = zero
        seid = zero
        model.grid._save(node_id, cp, cd, xyz, ps, seid, comment=None)
        model.grid._xyz_cid0 = xyz

        property_id = np.full(nelement, pid, dtype='int32')

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

        setup_modal_case_control(model, nmodes=nmodes)
        model.setup()
        assert len(model.cbar)
        solver = Solver(model)
        with self.assertRaises(FatalError):
            # arpack convergence error
            solver.run()

    def test_ctetra10(self):
        bdf_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.bdf'
        model = read_bdf(bdf_filename)
        model.bdf_filename = TEST_DIR / 'solid_bending_modes.bdf'
        setup_modal_case_control(model, nmodes=6)
        solver = Solver(model)
        solver.run()

    def test_ctetra10_cb(self):
        bdf_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.bdf'
        model = read_bdf(bdf_filename)
        model.bdf_filename = TEST_DIR / 'ctetra10_cb.bdf'

        model.add_suport([10], ['123'])
        setup_modal_case_control(model, nmodes=6)

        solver = Solver(model)
        model.sol = 31
        solver.run()

    def test_ctetra4(self):
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'ctetra4_modes.bdf'
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 0., 1.])

        mid = 3
        pid = 4
        model.add_ctetra(10, pid, [1, 2, 3, 4])
        model.add_psolid(pid, mid)

        spc_id = 3
        model.add_spc1(spc_id, '123456', [1, 2, 3], comment='')

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        model.add_param('GRDPNT', 0)

        setup_modal_case_control(model)
        solver = Solver(model)
        solver.run()
        os.remove(solver.f06_filename)
        os.remove(solver.op2_filename)

    def test_hexa8(self):
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'chexa8_modes.bdf'
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
        model.add_spc1(spc_id, '123456', [1, 2, 3], comment='')

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        model.add_param('GRDPNT', 0)

        setup_modal_case_control(model)
        solver = Solver(model)
        solver.run()

    def test_penta6(self):
        model = BDF(debug=None, log=None, mode='msc')
        model.bdf_filename = TEST_DIR / 'cpenta6_modes.bdf'
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

        E = 1.0E7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.0, alpha=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        model.add_param('GRDPNT', 0)
        setup_modal_case_control(model)
        solver = Solver(model)
        solver.run()
        os.remove(solver.f06_filename)
        os.remove(solver.op2_filename)


def setup_modal_case_control(model: BDF, extra_case_lines=None,
                             nmodes: int=2):
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
    ]
    if extra_case_lines is not None:
        lines += extra_case_lines
    cc = CaseControlDeck(lines, log=model.log)
    model.sol = 103
    model.case_control_deck = cc

    sid = 103
    model.add_eigrl(sid, v1=None, v2=None, nd=nmodes, msglvl=0,
                    maxset=None, shfscl=None, norm=None, options=None, values=None,comment='')

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



if __name__ == '__main__':  # pragma: no cover
    unittest.main()
