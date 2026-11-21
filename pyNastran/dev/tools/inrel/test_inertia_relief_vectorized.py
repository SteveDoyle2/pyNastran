"""Tests for the bdf_vectorized3 inertia relief adapter and the 6x6 core."""

import unittest

import numpy as np
from cpylog import SimpleLogger

from pyNastran.bdf.bdf import BDF as BDF_old
from pyNastran.bdf.mesh_utils.loads import sum_forces_moments
from pyNastran.bdf.mesh_utils.mass_properties import mass_properties
from pyNastran.dev.bdf_vectorized3.bdf import BDF
from pyNastran.dev.bdf_vectorized3.solver.matrices.mass_matrix import build_mgg_lumped
from pyNastran.dev.tools.inrel.inertia_relief2 import (
    get_mass_blocks,
    inertia_relief,
    inertia_relief_6x6,
    lumped_blocks_to_point_masses,
)
from pyNastran.dev.tools.inrel.inertia_relief_vectorized import (
    add_inertia_relief_cards,
    get_grid_loads,
    inertia_relief_from_model,
    inertia_relief_vectorized,
)


def _add_cards(model: BDF | BDF_old) -> None:
    """Free-free model: 4 grids, offset/rotated CONM2s with products of inertia, a CONROD."""
    model.add_cord2r(10, [1.0, 2.0, 3.0], [1.0, 2.0, 4.0], [2.0, 3.0, 3.0])
    model.add_grid(1, [0.0, 0.0, 0.0])
    model.add_grid(2, [10.0, 1.0, 0.0])
    model.add_grid(3, [3.0, 8.0, 2.0])
    model.add_grid(4, [5.0, 3.0, 7.0])
    model.add_conm2(11, 1, 2.0, X=[0.5, -0.3, 0.2], I=[3.0, 0.4, 2.0, -0.3, 0.2, 1.5])
    model.add_conm2(12, 2, 1.5, cid=10, X=[0.1, 0.7, -0.4], I=[1.0, 0.3, 2.5, 0.1, -0.2, 0.8])
    model.add_conm2(13, 3, 3.0)
    model.add_conm2(14, 4, 0.5, cid=-1, X=[5.2, 2.6, 7.1], I=[0.2, 0.0, 0.3, 0.0, 0.0, 0.4])
    model.add_conrod(21, 100, [1, 3], A=0.5)
    model.add_mat1(100, 1.0e7, None, 0.3, rho=0.1)

    # loads: two FORCEs on grid 2 (one in CID 10) must accumulate, and
    # LOAD 100 scales the sets
    model.add_force(1, 2, 100.0, [0.0, 0.0, 1.0])
    model.add_force(1, 2, 50.0, [1.0, 0.0, 0.0], cid=10)
    model.add_force(1, 4, 30.0, [0.0, 1.0, 1.0])
    model.add_moment(2, 3, 20.0, [1.0, 1.0, 0.0])
    model.add_load(100, 2.0, [1.5, 0.5], [1, 2])
    model.add_grav(200, 9.81, [0.0, 0.0, -1.0])


def _add_offset_conm2_cards(model: BDF) -> None:
    """Three grids, each carrying only one offset CONM2 (CID=0, rotated CID, CID=-1)."""
    model.add_cord2r(10, [1.0, 2.0, 3.0], [1.0, 2.0, 4.0], [2.0, 3.0, 3.0])
    model.add_grid(1, [1.0, 2.0, 3.0])
    model.add_grid(2, [-2.0, 1.0, 4.0])
    model.add_grid(3, [3.0, -1.0, 0.5])
    model.add_conm2(1, 1, 2.0, X=[0.5, -0.3, 0.2], I=[3.0, 0.4, 2.0, -0.3, 0.2, 1.5])
    model.add_conm2(2, 2, 1.5, cid=10, X=[0.1, 0.7, -0.4], I=[1.0, 0.3, 2.5, 0.1, -0.2, 0.8])
    model.add_conm2(3, 3, 0.5, cid=-1, X=[3.4, -0.6, 1.1], I=[0.2, 0.05, 0.3, 0.0, 0.1, 0.4])


def _resultant(xyz: np.ndarray, force: np.ndarray, moment: np.ndarray) -> np.ndarray:
    """Net force and moment about the origin."""
    return np.hstack([force.sum(axis=0), moment.sum(axis=0) + np.cross(xyz, force).sum(axis=0)])


class TestInertiaReliefVectorized(unittest.TestCase):
    def test_inertia_relief_vectorized(self):
        """
        Builds the same free-free model in bdf_vectorized3 and pyNastran.bdf.

        Tests (tolerances):
        - lumped point masses reproduce the classic mass and cg, and the QRG
          inertia tensor (rtol=1e-10); this checks CONM2 offsets, the CID=10
          rotation, CID=-1 absolute location, and the product-of-inertia
          sign in build_mgg_lumped
        - force and moment equilibrium of the relieved loads (atol=1e-9)
        - grid-level offsets: inertia_relief_6x6 matches inertia_relief when
          the masses are moved to the grids (atol=1e-12)
        """
        log = SimpleLogger(level="warning")
        model = BDF(log=log)
        _add_cards(model)
        model.setup()

        model_old = BDF_old(log=log)
        _add_cards(model_old)
        model_old.cross_reference()
        mass_expected, cg_expected, unused_inertia = mass_properties(model_old)

        # Expected inertia tensor about the origin, built from the QRG CONM2
        # tensor (CONM2.Inertia(), rotated to basic) plus the CONROD lumped
        # half-mass at each end (mass_properties treats a CONROD differently
        # than the lumped matrix, so its inertia isn't used).  CONM2 Iij follow
        # QRG CONM2 remark 5 (I21 = int(rho x1 x2 dV); the matrix term is -I21).
        def point_tensor(mass: float, xyz: np.ndarray) -> np.ndarray:
            return mass * (xyz @ xyz * np.eye(3) - np.outer(xyz, xyz))

        tensor_expected = np.zeros((3, 3))
        for elem in model_old.masses.values():
            tensor_expected += elem.Inertia() + point_tensor(elem.mass, elem.Centroid())
        rod = model_old.elements[21]
        rod_mass = rod.Mass()
        for nid in rod.node_ids:
            tensor_expected += point_tensor(rod_mass / 2, model_old.nodes[nid].get_position())

        mbb, grid_id, unused_total_mass = build_mgg_lumped(model=model)
        mass_blocks, coupling_ratio = get_mass_blocks(mbb)
        assert coupling_ratio == 0.0, coupling_ratio
        mass, offset, inertia = lumped_blocks_to_point_masses(mass_blocks)
        xyz = model.grid.xyz_cid0() + offset

        # rigid-body mass properties about the origin
        mass_total = mass.sum()
        cg = (mass[:, np.newaxis] * xyz).sum(axis=0) / mass_total
        tensor = np.zeros((3, 3))
        for massi, xyzi, (ixx, iyy, izz, ixy, ixz, iyz) in zip(mass, xyz, inertia):
            tensor += point_tensor(massi, xyzi) + np.array(
                [
                    [ixx, -ixy, -ixz],
                    [-ixy, iyy, -iyz],
                    [-ixz, -iyz, izz],
                ]
            )
        assert np.allclose(mass_total, mass_expected, rtol=1e-10), (mass_total, mass_expected)
        assert np.allclose(cg, cg_expected, rtol=1e-10), (cg, cg_expected)
        assert np.allclose(tensor, tensor_expected, rtol=1e-10), (tensor, tensor_expected)

        # equilibrium of the relieved loads
        rng = np.random.default_rng(3)
        ngrid = len(grid_id)
        force = rng.normal(scale=100.0, size=(ngrid, 3))
        moment = rng.normal(scale=100.0, size=(ngrid, 3))
        grid_id2, dforce, dmoment = inertia_relief_vectorized(model, force, moment)
        assert np.array_equal(grid_id, grid_id2)
        resultant = _resultant(model.grid.xyz_cid0(), force + dforce, moment + dmoment)
        assert np.allclose(resultant, 0.0, atol=1e-9), resultant

        # with the offsets folded in, the 6x6 path is the point-mass path
        xyz_grid = model.grid.xyz_cid0()
        blocks = np.zeros((ngrid, 6, 6))
        blocks[:, [0, 1, 2], [0, 1, 2]] = mass[:, np.newaxis]
        ixx, iyy, izz, ixy, ixz, iyz = inertia.T
        blocks[:, 3:, 3:] = np.stack(
            [
                np.column_stack([ixx, -ixy, -ixz]),
                np.column_stack([-ixy, iyy, -iyz]),
                np.column_stack([-ixz, -iyz, izz]),
            ],
            axis=1,
        )
        dforce1, dmoment1 = inertia_relief_6x6(blocks, xyz_grid, force, moment)
        dforce2, dmoment2 = inertia_relief(mass, xyz_grid, inertia, force, moment)
        assert np.allclose(dforce1, dforce2, atol=1e-12)
        assert np.allclose(dmoment1, dmoment2, atol=1e-12)

    def test_inertia_relief_from_model(self):
        """
        inertia_relief_from_model on a load set id.

        Tests (tolerances):
        - the assembled applied load resultant about the origin matches the
          classic pyNastran.bdf sum_forces_moments for FORCE sets (one in a
          rotated CID, two on one grid), a MOMENT, and a scaled LOAD
          combination (rtol=1e-10)
        - the relieved loads are in force and moment equilibrium (atol=1e-9)
        - GRAV: applied force = total mass * g, and the relief cancels it
          exactly (atol=1e-9)
        """
        log = SimpleLogger(level="warning")
        model = BDF(log=log)
        _add_cards(model)
        model.setup()
        xyz = model.grid.xyz_cid0()

        model_old = BDF_old(log=log)
        _add_cards(model_old)
        model_old.cross_reference()

        p0 = np.zeros(3)
        for load_id in (1, 2, 100):
            grid_id, force, moment, dforce, dmoment = inertia_relief_from_model(model, load_id)
            assert np.array_equal(grid_id, model.grid.node_id)
            force_expected, moment_expected = sum_forces_moments(model_old, p0, load_id)
            resultant = _resultant(xyz, force, moment)
            assert np.allclose(resultant[:3], force_expected, rtol=1e-10), (
                load_id,
                resultant[:3],
                force_expected,
            )
            assert np.allclose(resultant[3:], moment_expected, rtol=1e-10), (
                load_id,
                resultant[3:],
                moment_expected,
            )
            assert np.abs(resultant).max() > 0.0, load_id

            balanced = _resultant(xyz, force + dforce, moment + dmoment)
            assert np.allclose(balanced, 0.0, atol=1e-9), (load_id, balanced)

        # GRAV: a uniform acceleration is fully relieved
        grid_id, force, moment, dforce, dmoment = inertia_relief_from_model(model, 200)
        mass_total = mass_properties(model_old)[0]
        assert np.allclose(force.sum(axis=0), [0.0, 0.0, -9.81 * mass_total], atol=1e-9)
        assert np.allclose(force + dforce, 0.0, atol=1e-9), force + dforce

    def test_add_inertia_relief_cards(self):
        """
        GRAV + RFORCE + LOAD cards balance a load set.

        Tests (tolerances):
        - for FORCE sets (one in a rotated CID), a MOMENT, and a LOAD
          combination, the assembled balanced LOAD has a zero force/moment
          resultant (atol=1e-9); the model has offset and rotated CONM2s with
          rotary inertia, which RFORCE ignores and the solve accounts for
        - a LOAD load_id is expanded, so the new LOAD never references a LOAD
        - with point masses only (no offsets/rotary inertia), the GRAV+RFORCE
          grid loads equal the inertia_relief_vectorized relief (atol=1e-9)
        """
        log = SimpleLogger(level="warning")
        for load_id in (1, 2, 100):
            model = BDF(log=log)
            _add_cards(model)
            model.setup()
            unused_g, unused_alpha, residual = add_inertia_relief_cards(
                model, load_id, 1000, 1001, 1002
            )
            assert np.allclose(residual, 0.0, atol=1e-9), (load_id, residual)
            irow = np.where(model.load.load_id == 1000)[0][0]
            i0, i1 = model.load.iload[irow]
            assert not np.isin(model.load.load_ids[i0:i1], model.load.load_id).any()

            force, moment = get_grid_loads(model, 1000)
            balanced = _resultant(model.grid.xyz_cid0(), force, moment)
            assert np.allclose(balanced, 0.0, atol=1e-9), (load_id, balanced)

        # point masses: RFORCE/GRAV physics == the exact relief
        model = BDF(log=log)
        rng = np.random.default_rng(4)
        for nid, xyz in enumerate(rng.normal(scale=5.0, size=(6, 3)), start=1):
            model.add_grid(nid, xyz)
            model.add_conm2(nid, nid, rng.uniform(0.5, 2.0))
            model.add_force(1, nid, 1.0, rng.normal(size=3))
            model.add_moment(1, nid, 1.0, rng.normal(size=3))
        model.setup()
        grid_id, force, moment, dforce, dmoment = inertia_relief_from_model(model, 1)
        add_inertia_relief_cards(model, 1, 1000, 1001, 1002)
        model.add_load(1003, 1.0, [1.0, 1.0], [1001, 1002])
        model.setup()
        force_cards, moment_cards = get_grid_loads(model, 1003)
        assert np.allclose(force_cards, dforce, atol=1e-9), force_cards - dforce
        assert np.allclose(moment_cards, dmoment, atol=1e-9), moment_cards - dmoment

    def test_rforce_load(self):
        """
        Vectorized RFORCE assembly vs the QRG formula.

        F_i = -m_i [w x (w x (r_i - r_a)) + a x (r_i - r_a)],
        w = 2 pi A R, a = 2 pi RACC R (R not normalized).

        Tests (atol=1e-10): grid forces for A and RACC both non-zero, a
        non-unit R, G = a grid, and LOAD scaling of the RFORCE; no grid moments.
        """
        log = SimpleLogger(level="warning")
        model = BDF(log=log)
        xyz = np.array([[0.0, 0.0, 0.0], [2.0, 1.0, 0.5], [-1.0, 3.0, 2.0], [4.0, -2.0, 1.0]])
        mass = np.array([1.0, 2.0, 0.5, 3.0])
        for nid, (xyzi, massi) in enumerate(zip(xyz, mass), start=1):
            model.add_grid(nid, xyzi)
            model.add_conm2(nid, nid, massi)
        r123 = np.array([0.3, -1.2, 2.0])
        scale_a, racc = 0.7, 1.3
        model.add_rforce(5, 2, scale_a, r123, method=2, racc=racc)
        model.add_load(6, 2.0, [1.5], [5])
        model.setup()

        omega = 2 * np.pi * scale_a * r123
        alpha = 2 * np.pi * racc * r123
        dxyz = xyz - xyz[1]
        accel = np.cross(omega, np.cross(omega, dxyz)) + np.cross(alpha, dxyz)
        force_expected = -mass[:, np.newaxis] * accel
        for load_id, factor in ((5, 1.0), (6, 3.0)):
            force, moment = get_grid_loads(model, load_id)
            assert np.allclose(force, factor * force_expected, atol=1e-10), load_id
            assert np.allclose(moment, 0.0, atol=1e-10), load_id

    def test_rforce_method1_matches_nastran(self):
        """
        RFORCE METHOD=1 angular acceleration vs NX 2412 OLOAD reference data.

        The reference is NX OLOAD for RFORCE (A=0, G=0, RACC=1/(2 pi), R = x,
        y, z) on the _add_cards model plus three CBARs; offset, rotated, and
        CID=-1 CONM2s with rotary inertia.  NX double-counts the CONM2 offset
        (see rforce_alpha_operator).  Tolerance: atol=1e-4 (the OLOAD values
        are printed to 5 decimals).
        """
        # nastran_oload[grid, dof, alpha_direction]
        nastran_oload = np.array(
            [
                [
                    [0.0, -0.8, -1.2],
                    [0.8, -0.0, -2.0],
                    [1.2, 2.0, -0.0],
                    [-3.45104, -0.12044, 0.04696],
                    [-0.12044, -3.00618, -0.00817],
                    [0.04696, -0.00817, -2.67966],
                ],
                [
                    [-0.0, 1.2, 4.20452],
                    [-1.2, 0.0, -23.80189],
                    [-4.20452, 23.80189, -0.0],
                    [-4.04924, 8.65992, -5.38103],
                    [-0.46175, 4.09887, -1.21316],
                    [0.61897, -0.61316, 3.51677],
                ],
                [
                    [-0.0, -8.18354, 32.73417],
                    [8.18354, 0.0, -12.27532],
                    [-32.73417, 12.27532, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                ],
                [
                    [0.0, -6.17196, 2.20227],
                    [6.17196, 0.0, -4.53712],
                    [-2.20227, 4.53712, -0.0],
                    [-0.084, -1.06306, 0.26576],
                    [0.23694, -1.18941, 0.11847],
                    [0.71576, -1.43153, -0.45764],
                ],
            ]
        )
        log = SimpleLogger(level="warning")
        model = BDF(log=log)
        _add_cards(model)
        model.add_pbar(200, 100, area=1.0, i1=1.0, i2=1.0, j=1.0)
        for eid, (n1, n2) in enumerate([(1, 2), (2, 3), (3, 4)], start=31):
            model.add_cbar(eid, 200, [n1, n2], [0.0, 0.0, 1.0], None)
        for j in range(3):
            r123 = np.zeros(3)
            r123[j] = 1.0
            model.add_rforce(3001 + j, 0, 0.0, r123, method=1, racc=1.0 / (2 * np.pi))
        model.setup()
        for j in range(3):
            force, moment = get_grid_loads(model, 3001 + j)
            loads = np.hstack([force, moment])
            assert np.allclose(loads, nastran_oload[:, :, j], atol=1e-4), (
                j,
                loads - nastran_oload[:, :, j],
            )

    def test_offset_conm2(self):
        """
        Offset CONM2s alone on their grids: GRAV, RFORCE METHOD=1, and balancing.

        Each grid carries only one CONM2 (offset in CID=0, a rotated CID, and
        CID=-1), so nothing else at the grid hides the offset terms.

        Tests (tolerances):
        - GRAV with a non-unit N (|a| = A*|N|, QRG GRAV remark 1):
          F = m g, M = m e x g at the grid (atol=1e-10 vs hand formula)
        - RFORCE METHOD=1 RACC, unit alpha along x, y, z:
          F = -m alpha x (r + 2e), M = -[m e x (alpha x (r + e)) +
          (I0 + m((e.e) I3 - e e^T)) alpha]; NX counts the offset twice
          (atol=1e-10 vs hand formula)
        - both match NX 2412 OLOAD reference data for this model
          (atol=1e-5; OLOAD is float32)
        - add_inertia_relief_cards balances an applied load on this model
          (atol=1e-9)
        """
        # NX 2412 OLOAD, [load, grid, dof]: GRAV 4001, then RFORCE 3001-3003
        nastran_oload = np.array(
            [
                [
                    [5.886, -9.81, 19.62, -3.924, -8.6328, -3.1392],
                    [4.4145, -7.3575, 14.715, 5.381061, 4.477246, 0.6243046],
                    [1.4715, -2.4525, 4.905, 3.4335, -1.0791, -1.5696],
                ],
                [
                    [0.0, 6.8, -2.8, -3.52, 1.8, 3.1],
                    [0.0, 4.8, -3.197056, -1.938528, -0.6063961, -1.824335],
                    [0.0, 0.85, 0.1, -0.67, 0.01, 0.34],
                ],
                [
                    [-6.8, 0.0, 4.0, -0.8, -5.36, -1.84],
                    [-4.8, 0.0, -4.272792, -1.667056, -1.342792, 2.644579],
                    [-0.85, 0.0, 1.9, 0.81, -1.57, 0.44],
                ],
                [
                    [2.8, -4.0, 0.0, 0.5, 0.76, -2.66],
                    [3.197056, 4.272792, 0.0, 1.921249, -1.349533, -4.42132],
                    [-0.1, -1.9, 0.0, 1.14, 0.04, -1.12],
                ],
            ]
        )
        log = SimpleLogger(level="warning")
        model = BDF(log=log)
        _add_offset_conm2_cards(model)
        grav_n = np.array([0.3, -0.5, 1.0])
        model.add_grav(4001, 9.81, grav_n)
        for j in range(3):
            r123 = np.zeros(3)
            r123[j] = 1.0
            model.add_rforce(3001 + j, 0, 0.0, r123, method=1, racc=1.0 / (2 * np.pi))
        model.setup()
        xyz = model.grid.xyz_cid0()
        beta = model.coord.xyz_to_global_transform[10]  # rows: CID 10 axes

        # hand formulas from the card values
        cards = [  # mass, X, [I11, I21, I22, I31, I32, I33], cid
            (2.0, [0.5, -0.3, 0.2], [3.0, 0.4, 2.0, -0.3, 0.2, 1.5], 0),
            (1.5, [0.1, 0.7, -0.4], [1.0, 0.3, 2.5, 0.1, -0.2, 0.8], 10),
            (0.5, [3.4, -0.6, 1.1], [0.2, 0.05, 0.3, 0.0, 0.1, 0.4], -1),
        ]
        accel = 9.81 * grav_n
        expected = np.zeros((4, 3, 6))
        for igrid, (mass, x123, (i11, i21, i22, i31, i32, i33), cid) in enumerate(cards):
            r = xyz[igrid]
            tensor = np.array([[i11, -i21, -i31], [-i21, i22, -i32], [-i31, -i32, i33]])
            if cid == 0:
                e = np.array(x123)
            elif cid == -1:
                e = np.array(x123) - r  # X1-X3 are the cg location in basic
            else:
                e = beta.T @ np.array(x123)
                tensor = beta.T @ tensor @ beta
            expected[0, igrid, :3] = mass * accel
            expected[0, igrid, 3:] = mass * np.cross(e, accel)
            parallel = mass * (e @ e * np.eye(3) - np.outer(e, e))
            for j in range(3):
                alpha = np.zeros(3)
                alpha[j] = 1.0
                expected[1 + j, igrid, :3] = -mass * np.cross(alpha, r + 2 * e)
                expected[1 + j, igrid, 3:] = -(
                    mass * np.cross(e, np.cross(alpha, r + e)) + (tensor + parallel) @ alpha
                )

        for k, load_id in enumerate([4001, 3001, 3002, 3003]):
            force, moment = get_grid_loads(model, load_id)
            loads = np.hstack([force, moment])
            assert np.allclose(loads, expected[k], atol=1e-10), (load_id, loads - expected[k])
            assert np.allclose(loads, nastran_oload[k], atol=1e-5), (
                load_id,
                loads - nastran_oload[k],
            )

        # balance an applied load on the offset-only model
        rng = np.random.default_rng(9)
        for nid in (1, 2, 3):
            model.add_force(1, nid, 1.0, rng.normal(scale=10.0, size=3))
            model.add_moment(1, nid, 1.0, rng.normal(scale=10.0, size=3))
        model.setup()
        unused_g, unused_alpha, residual = add_inertia_relief_cards(model, 1, 1000, 1001, 1002)
        assert np.allclose(residual, 0.0, atol=1e-9), residual
        force, moment = get_grid_loads(model, 1000)
        balanced = _resultant(xyz, force, moment)
        assert np.allclose(balanced, 0.0, atol=1e-9), balanced

    def test_cards_classic_sum_forces_moments(self):
        """
        The GRAV/RFORCE balancing cards also balance in pyNastran.bdf.

        The cards from the vectorized add_inertia_relief_cards are added to a
        classic BDF of the same deck (offset/rotated CONM2s with rotary
        inertia, CBARs, a CONROD).

        Tests (tolerances):
        - classic sum_forces_moments of the balanced LOAD is ~0
          (atol=1e-8 relative to the applied load)
        - classic and vectorized give the same applied, GRAV, and RFORCE
          resultants (atol=1e-9)
        """
        log = SimpleLogger(level="warning")

        def add_bars(model: BDF | BDF_old) -> None:
            model.add_pbar(200, 100, 1.0, 1.0, 1.0, 0.0, 1.0)
            for eid, (n1, n2) in enumerate([(1, 2), (2, 3), (3, 4)], start=31):
                model.add_cbar(eid, 200, [n1, n2], [0.0, 0.0, 1.0], None)

        model = BDF(log=log)
        _add_cards(model)
        add_bars(model)
        model.setup()
        add_inertia_relief_cards(model, 100, 1000, 1001, 1002)
        grav = model.grav.slice_card_by_id(1001)
        rforce = model.rforce.slice_card_by_id(1002)

        model_old = BDF_old(log=log)
        _add_cards(model_old)
        add_bars(model_old)
        model_old.add_grav(1001, grav.scale[0], grav.N[0])
        model_old.add_rforce(
            1002, 0, 0.0, rforce.r[0], method=int(rforce.method[0]), racc=rforce.racc[0]
        )
        # LOAD 100 = 2.0 * (1.5 * set 1 + 0.5 * set 2), copied into LOAD 1000
        model_old.add_load(1000, 1.0, [3.0, 1.0, 1.0, 1.0], [1, 2, 1001, 1002])
        model_old.cross_reference()

        p0 = np.zeros(3)
        xyz = model.grid.xyz_cid0()
        applied = None
        for load_id in (100, 1001, 1002, 1000):
            force_old, moment_old = sum_forces_moments(
                model_old, p0, load_id, include_grav=True, include_rforce=True
            )
            force, moment = get_grid_loads(model, load_id)
            resultant = _resultant(xyz, force, moment)
            assert np.allclose(np.hstack([force_old, moment_old]), resultant, atol=1e-9), (
                load_id,
                np.hstack([force_old, moment_old]) - resultant,
            )
            if load_id == 100:
                applied = np.abs(resultant).max()
        assert np.allclose(np.hstack([force_old, moment_old]), 0.0, atol=1e-8 * applied), (
            force_old,
            moment_old,
        )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
