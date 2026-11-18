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
    inertia_relief,
    inertia_relief_6x6,
    lumped_blocks_to_point_masses,
)
from pyNastran.dev.tools.inrel.inertia_relief_vectorized import (
    get_mass_blocks,
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


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
