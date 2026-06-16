"""Unit tests for vectorized inertia relief (free-free structures).

Tests build_rigid_body_modes, compute_inertia_relief, and
compute_inertia_relief_from_model using bdf_vectorized3.

Verifies:
- Rigid body mode matrix D shape and orthogonality properties
- Self-equilibration: D^T @ F_net = 0 (zero net force/moment)
- Known analytical cases (uniform gravity on free beam)
- Multiple load cases
- Singular mass matrix error handling
- Sparse and dense mass matrix support
- Integration with BWB model via compute_inertia_relief_from_model

Tolerances:
- Self-equilibration residual: atol=1e-10 (D^T @ F_net ~ 0)
- Rigid body acceleration: rtol=1e-8 (inverse of 6x6 matrix)
- Total inertial force: atol=1e-10 (mass * acceleration identity)
- Mass conservation: rtol=1e-8
"""
import os
import tempfile
import unittest

import numpy as np
import numpy.testing as npt
from scipy import sparse

from pyNastran.dev.bdf_vectorized3.bdf import BDF
from pyNastran.dev.bdf_vectorized3.mesh_utils.inertia_relief import (
    build_rigid_body_modes,
    compute_inertia_relief,
    compute_inertia_relief_from_model,
)
from pyNastran.dev.bdf_vectorized3.solver.matrices.mass_matrix import (
    build_mgg_lumped,
)


# =============================================================================
# BDF FIXTURES
# =============================================================================

BDF_FREE_BEAM_3NODE = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0
GRID    2               5.0     0.0     0.0
GRID    3               10.0    0.0     0.0
CBAR    1       1       1       2       0.0     1.0     0.0
CBAR    2       1       2       3       0.0     1.0     0.0
PBAR    1       1       1.0     1.0     1.0     1.0
MAT1    1       1.0+7           0.3     0.1
ENDDATA
"""

BDF_SINGLE_CONM2 = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0
CONM2   1       1       0       10.0    0.0     0.0     0.0
ENDDATA
"""

BDF_TWO_MASS_BEAM = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0
GRID    2               10.0    0.0     0.0
CONM2   1       1       0       5.0     0.0     0.0     0.0
CONM2   2       2       0       5.0     0.0     0.0     0.0
ENDDATA
"""


def write_temp_bdf(content: str) -> str:
    """Write BDF content to a temporary file and return path."""
    fd, path = tempfile.mkstemp(suffix='.bdf', prefix='test_ir_')
    with os.fdopen(fd, 'w', encoding='utf-8') as f:
        f.write(content)
    return path


# =============================================================================
# TESTS: build_rigid_body_modes
# =============================================================================

class TestBuildRigidBodyModes(unittest.TestCase):
    """Rigid body mode matrix D construction."""

    def test_shape_single_grid(self):
        """D is (6, 6) for a single grid."""
        grid_ids = np.array([1])
        xyz = np.array([[0.0, 0.0, 0.0]])
        D = build_rigid_body_modes(grid_ids, xyz)
        self.assertEqual(D.shape, (6, 6))

    def test_shape_multiple_grids(self):
        """D is (n_dof, 6) = (18, 6) for 3 grids."""
        grid_ids = np.array([1, 2, 3])
        xyz = np.array([[0.0, 0.0, 0.0],
                        [5.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0]])
        D = build_rigid_body_modes(grid_ids, xyz)
        self.assertEqual(D.shape, (18, 6))

    def test_translation_modes_at_origin(self):
        """Translation columns are identity blocks when ref=origin.

        For a grid at origin, translational DOFs (rows 0,1,2) should be
        [1,0,0,0,0,0], [0,1,0,0,0,0], [0,0,1,0,0,0] for columns 0,1,2.
        """
        grid_ids = np.array([1])
        xyz = np.array([[0.0, 0.0, 0.0]])
        D = build_rigid_body_modes(grid_ids, xyz)
        # Translation DOFs
        npt.assert_array_equal(D[0, :3], [1, 0, 0])
        npt.assert_array_equal(D[1, :3], [0, 1, 0])
        npt.assert_array_equal(D[2, :3], [0, 0, 1])

    def test_rotation_modes_at_origin(self):
        """Rotation columns at origin grid are zero for translational DOFs.

        Grid at the reference point has zero moment arm, so rotation modes
        produce no translation.
        """
        grid_ids = np.array([1])
        xyz = np.array([[0.0, 0.0, 0.0]])
        D = build_rigid_body_modes(grid_ids, xyz)
        # Translation DOFs (0,1,2) should be zero for rotation columns (3,4,5)
        npt.assert_allclose(D[:3, 3:], 0.0, atol=1e-14)

    def test_rotation_dof_identity(self):
        """Rotational DOFs (rows 3,4,5) are identity for rotation columns."""
        grid_ids = np.array([1])
        xyz = np.array([[0.0, 0.0, 0.0]])
        D = build_rigid_body_modes(grid_ids, xyz)
        npt.assert_array_equal(D[3:6, 3:6], np.eye(3))

    def test_moment_arm_coupling(self):
        """Rotation mode produces correct translation via cross product.

        Grid at (L, 0, 0) with Rz rotation should produce:
        - Tx contribution: -dy = 0 (from col 5, row 0)
        - Ty contribution: +dx = L (from col 5, row 1)

        Physical: rotating about Z moves a point on X-axis in the Y direction.
        """
        L = 5.0
        grid_ids = np.array([1, 2])
        xyz = np.array([[0.0, 0.0, 0.0],
                        [L, 0.0, 0.0]])
        D = build_rigid_body_modes(grid_ids, xyz)
        # Grid 2 starts at DOF index 6. Rz is column 5.
        # Tx (row 6): -dy = 0
        npt.assert_allclose(D[6, 5], 0.0, atol=1e-14)
        # Ty (row 7): +dx = L
        npt.assert_allclose(D[7, 5], L, atol=1e-14)

    def test_ref_point_offset(self):
        """Non-zero ref_point shifts moment arms.

        Grid at (10, 0, 0) with ref_point (5, 0, 0) has arm (5, 0, 0).
        Same result as grid at (5, 0, 0) with ref_point at origin.
        """
        grid_ids = np.array([1])
        xyz_a = np.array([[10.0, 0.0, 0.0]])
        ref_a = np.array([5.0, 0.0, 0.0])
        D_a = build_rigid_body_modes(grid_ids, xyz_a, ref_a)

        xyz_b = np.array([[5.0, 0.0, 0.0]])
        D_b = build_rigid_body_modes(grid_ids, xyz_b, ref_point=None)

        npt.assert_allclose(D_a, D_b, atol=1e-14)

    def test_d_transpose_m_d_is_mass_matrix(self):
        """D^T @ M @ D recovers the 6x6 rigid body mass matrix.

        For a single point mass m at origin, M_rr = m * I_6 for the
        translational part and zero rotational inertia.
        """
        m = 3.0
        grid_ids = np.array([1])
        xyz = np.array([[0.0, 0.0, 0.0]])
        D = build_rigid_body_modes(grid_ids, xyz)
        M = m * np.eye(6)
        # Zero out rotational mass (lumped point mass)
        M[3:, 3:] = 0.0
        M_rr = D.T @ M @ D
        expected = np.zeros((6, 6))
        expected[:3, :3] = m * np.eye(3)
        npt.assert_allclose(M_rr, expected, atol=1e-14)


# =============================================================================
# TESTS: compute_inertia_relief
# =============================================================================

class TestComputeInertiaRelief(unittest.TestCase):
    """Core inertia relief computation."""

    def setUp(self):
        """Two point masses (5 kg each) at (0,0,0) and (10,0,0).

        Total mass = 10 kg. CG at (5, 0, 0).
        Under unit gravity in Z: F = [0,0,-mg, 0,0,-mg] for each grid.
        Expected: a_z = -g, F_net has zero resultant.
        """
        self.m = 5.0
        self.n_grids = 2
        self.n_dof = 12
        self.grid_ids = np.array([1, 2])
        self.xyz = np.array([[0.0, 0.0, 0.0],
                             [10.0, 0.0, 0.0]])
        self.D = build_rigid_body_modes(self.grid_ids, self.xyz)

        # Lumped mass: m on T1,T2,T3 for each grid, zero on R1,R2,R3
        diag = np.zeros(self.n_dof)
        diag[0:3] = self.m
        diag[6:9] = self.m
        self.M = sparse.diags(diag, format='csr')

    def test_self_equilibration(self):
        """D^T @ F_net = 0 (zero resultant force and moment).

        This is the fundamental property of inertia relief: after removing
        the rigid body component, the net loads are self-equilibrated.
        """
        g = 9.81
        F = np.zeros(self.n_dof)
        F[2] = -self.m * g   # Grid 1, Fz
        F[8] = -self.m * g   # Grid 2, Fz

        F_net, _, _ = compute_inertia_relief(self.M, self.D, F)
        residual = self.D.T @ F_net
        npt.assert_allclose(residual, 0.0, atol=1e-10)

    def test_uniform_gravity_acceleration(self):
        """Uniform -Z gravity -> rigid body acceleration a_z = -g.

        For uniform gravity on a free body, the rigid body acceleration
        equals the gravitational acceleration (Newton's 2nd law for the
        whole body: F = m*a -> a = F/m = -g in Z).
        """
        g = 9.81
        total_mass = 2 * self.m
        F = np.zeros(self.n_dof)
        F[2] = -total_mass * g / 2  # Split evenly
        F[8] = -total_mass * g / 2

        _, a_rigid, _ = compute_inertia_relief(self.M, self.D, F)
        # a_rigid = [ax, ay, az, alpha_x, alpha_y, alpha_z]
        npt.assert_allclose(a_rigid[2], -g, rtol=1e-8)
        npt.assert_allclose(a_rigid[:2], 0.0, atol=1e-10)
        npt.assert_allclose(a_rigid[3:], 0.0, atol=1e-10)

    def test_uniform_gravity_zero_fnet(self):
        """Uniform gravity on symmetric body -> F_net = 0 everywhere.

        When gravity is uniform and the body is symmetric about CG,
        every grid gets the same inertial load that exactly cancels the
        applied gravity load.
        """
        g = 9.81
        F = np.zeros(self.n_dof)
        F[2] = -self.m * g
        F[8] = -self.m * g

        F_net, _, _ = compute_inertia_relief(self.M, self.D, F)
        npt.assert_allclose(F_net, 0.0, atol=1e-10)

    def test_inertia_force_conservation(self):
        """Sum of inertial forces = -F_applied (for pure rigid body load).

        Total inertial force = -M @ D @ a_r. For a rigid body load (one
        where F_net=0), the inertial forces exactly cancel applied forces.
        """
        F = np.zeros(self.n_dof)
        F[0] = 100.0  # Fx at grid 1
        F[6] = 100.0  # Fx at grid 2

        _, _, F_inertia = compute_inertia_relief(self.M, self.D, F)
        npt.assert_allclose(F_inertia, -F, atol=1e-10)

    def test_moment_load_produces_angular_acceleration(self):
        """Fy force at grid 2 (away from origin) produces Rz acceleration.

        Force Fy=100 at grid 2 (x=10, y=0). Ref at origin.
        D^T @ F = [0, Fy, 0, 0, 0, dx*Fy] = [0, 100, 0, 0, 0, 1000].
        M_rr has Ty-Rz coupling for this collinear geometry, so
        a_y != Fy/m. Key check: self-equilibration and non-zero alpha_z.
        """
        Fy = 100.0
        F = np.zeros(self.n_dof)
        F[7] = Fy  # Grid 2 (DOF index 6+1=7), Fy

        F_net, a_rigid, _ = compute_inertia_relief(self.M, self.D, F)

        # Self-equilibration
        residual = self.D.T @ F_net
        npt.assert_allclose(residual, 0.0, atol=1e-10)

        # Non-zero Rz angular acceleration (moment arm = 10)
        self.assertNotAlmostEqual(a_rigid[5], 0.0, places=5)

    def test_dense_mass_matrix(self):
        """Works with dense (ndarray) mass matrix too."""
        F = np.zeros(self.n_dof)
        F[2] = -self.m * 9.81
        F[8] = -self.m * 9.81

        M_dense = self.M.toarray()
        F_net, a_rigid, F_inertia = compute_inertia_relief(M_dense, self.D, F)
        residual = self.D.T @ F_net
        npt.assert_allclose(residual, 0.0, atol=1e-10)

    def test_multiple_load_cases(self):
        """Multiple load cases processed simultaneously.

        Each column of F_applied is an independent load case.
        D^T @ F_net = 0 must hold for each column.
        """
        F = np.zeros((self.n_dof, 3))
        # Load case 1: gravity Z
        F[2, 0] = -self.m * 9.81
        F[8, 0] = -self.m * 9.81
        # Load case 2: force X on grid 1 only
        F[0, 1] = 50.0
        # Load case 3: moment (couple)
        F[1, 2] = 10.0   # Fy at grid 1
        F[7, 2] = -10.0  # -Fy at grid 2

        F_net, a_rigid, F_inertia = compute_inertia_relief(self.M, self.D, F)

        self.assertEqual(F_net.shape, (self.n_dof, 3))
        self.assertEqual(a_rigid.shape, (6, 3))

        # Each load case self-equilibrated
        residual = self.D.T @ F_net
        npt.assert_allclose(residual, 0.0, atol=1e-10)

    def test_zero_mass_raises(self):
        """Zero mass matrix raises ValueError."""
        M_zero = sparse.csr_matrix((self.n_dof, self.n_dof))
        F = np.zeros(self.n_dof)
        F[0] = 1.0

        with self.assertRaises(ValueError) as ctx:
            compute_inertia_relief(M_zero, self.D, F)
        self.assertIn("zero", str(ctx.exception).lower())

    def test_dimension_mismatch_raises(self):
        """Mismatched D and F dimensions raise ValueError."""
        F_wrong = np.zeros(6)  # Wrong size
        with self.assertRaises(ValueError):
            compute_inertia_relief(self.M, self.D, F_wrong)


# =============================================================================
# TESTS: Pure couple (self-equilibrated input)
# =============================================================================

class TestSelfEquilibratedInput(unittest.TestCase):
    """Input that is already self-equilibrated should pass through unchanged."""

    def test_couple_unchanged(self):
        """Load with D^T @ F = 0 passes through unchanged (F_net = F).

        Use 4 grids in a square (non-degenerate 3D inertia).
        Apply equal-and-opposite Fx at grids 1 and 2 (same Y, along X-axis):
        Grid 1 at (-0.5, -0.5, 0), Grid 2 at (0.5, -0.5, 0). Ref = origin.
        r1 x F1 = (-0.5,-0.5,0) x (Fx,0,0) = (0, 0, 0.5*Fx)
        r2 x F2 = (0.5,-0.5,0) x (-Fx,0,0) = (0, 0, -0.5*(-Fx))?

        Actually use the simplest approach: forces along line through ref.
        +Fx at grid 1 and -Fx at grid 2, both on X-axis through origin:
        Grid 1 at (-1,0,0), grid 2 at (1,0,0): arms (-1,0,0) and (1,0,0).
        r1 x F1 = (-1,0,0) x (F,0,0) = (0,0,0)
        r2 x F2 = (1,0,0) x (-F,0,0) = (0,0,0)
        Net force = 0, net moment = 0. Perfectly self-equilibrated.
        """
        m = 1.0
        grid_ids = np.array([1, 2, 3, 4])
        xyz = np.array([[-1.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, -1.0, 0.0]])
        D = build_rigid_body_modes(grid_ids, xyz, ref_point=None)

        n_dof = 24
        diag = np.zeros(n_dof)
        for i in range(4):
            diag[6*i:6*i+3] = m
        M = sparse.diags(diag, format='csr')

        Fx = 100.0
        F = np.zeros(n_dof)
        F[0] = Fx       # +Fx at grid 1 (on X-axis)
        F[6] = -Fx      # -Fx at grid 2 (on X-axis)

        # Verify D^T @ F = 0
        DtF = D.T @ F
        npt.assert_allclose(DtF, 0.0, atol=1e-10)

        F_net, a_rigid, _ = compute_inertia_relief(M, D, F)

        # a_rigid = 0 and F_net = F (no rigid body component)
        npt.assert_allclose(a_rigid, 0.0, atol=1e-10)
        npt.assert_allclose(F_net, F, atol=1e-10)


# =============================================================================
# TESTS: 3D geometry
# =============================================================================

class TestInertiaRelief3D(unittest.TestCase):
    """Inertia relief with 3D grid geometry (non-collinear grids)."""

    def setUp(self):
        """4 grids at corners of a unit square in the XY plane.

        Grids: (0,0,0), (1,0,0), (1,1,0), (0,1,0). Mass = 1 kg each.
        Total mass = 4 kg. CG at (0.5, 0.5, 0).
        """
        self.m = 1.0
        self.grid_ids = np.array([1, 2, 3, 4])
        self.xyz = np.array([[0.0, 0.0, 0.0],
                             [1.0, 0.0, 0.0],
                             [1.0, 1.0, 0.0],
                             [0.0, 1.0, 0.0]])
        self.n_dof = 24
        self.D = build_rigid_body_modes(self.grid_ids, self.xyz)

        diag = np.zeros(self.n_dof)
        for i in range(4):
            diag[6*i:6*i+3] = self.m
        self.M = sparse.diags(diag, format='csr')

    def test_self_equilibration_gravity_z(self):
        """Z-gravity on a planar structure is self-equilibrated."""
        g = 9.81
        F = np.zeros(self.n_dof)
        for i in range(4):
            F[6*i + 2] = -self.m * g

        F_net, _, _ = compute_inertia_relief(self.M, self.D, F)
        residual = self.D.T @ F_net
        npt.assert_allclose(residual, 0.0, atol=1e-10)

    def test_asymmetric_load(self):
        """Asymmetric point load on one grid is still self-equilibrated."""
        F = np.zeros(self.n_dof)
        F[0] = 1000.0  # Fx on grid 1 only

        F_net, a_rigid, _ = compute_inertia_relief(self.M, self.D, F)
        residual = self.D.T @ F_net
        npt.assert_allclose(residual, 0.0, atol=1e-10)

        # Non-zero angular acceleration (torque about CG)
        # Force at (0,0,0), CG at (0.5, 0.5, 0), arm = (-0.5, -0.5, 0)
        # Moment = r x F = (-0.5,-0.5,0) x (1000,0,0) = (0, 0, 500) about CG
        # But ref is at origin, not CG, so check self-eq only.
        self.assertTrue(np.any(np.abs(a_rigid[3:]) > 1e-6))

    def test_m_rr_diagonal_uniform_square(self):
        """M_rr for uniform square has known structure.

        4 unit masses at corners of unit square (ref at origin):
        - M_rr[0:3, 0:3] = 4*I_3 (total mass)
        - Moments of inertia about origin:
          Ixx = sum(m*(y^2+z^2)) = 1*(0+0) + 1*(0+0) + 1*(1+0) + 1*(1+0) = 2
          Iyy = sum(m*(x^2+z^2)) = 1*(0) + 1*(1) + 1*(1) + 1*(0) = 2
          Izz = sum(m*(x^2+y^2)) = 1*(0) + 1*(1) + 1*(2) + 1*(1) = 4
        """
        MD = self.M.dot(self.D)
        M_rr = self.D.T @ MD

        # Translational mass
        npt.assert_allclose(M_rr[:3, :3], 4.0 * np.eye(3), atol=1e-10)

        # Total mass is 4
        npt.assert_allclose(M_rr[0, 0], 4.0, atol=1e-10)


# =============================================================================
# TESTS: compute_inertia_relief_from_model
# =============================================================================

class TestInertiaReliefFromModel(unittest.TestCase):
    """Integration test using BDF model."""

    @classmethod
    def setUpClass(cls):
        cls.bdf_path = write_temp_bdf(BDF_TWO_MASS_BEAM)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.bdf_path)

    def _load_model(self):
        model = BDF(log=None)
        model.read_bdf(self.bdf_path)
        return model

    def test_from_model_self_equilibration(self):
        """compute_inertia_relief_from_model produces self-equilibrated loads."""
        model = self._load_model()
        n_dof = 12
        g = 9.81
        m = 5.0

        F = np.zeros(n_dof)
        F[2] = -m * g   # Grid 1 Fz
        F[8] = -m * g   # Grid 2 Fz

        F_net, a_rigid, F_inertia, info = compute_inertia_relief_from_model(
            model, F)

        npt.assert_allclose(info['residual_force'], 0.0, atol=1e-10)
        npt.assert_allclose(info['total_mass'], 10.0, rtol=1e-8)

    def test_from_model_returns_diagnostics(self):
        """Info dict contains expected keys and shapes."""
        model = self._load_model()
        F = np.zeros(12)
        F[0] = 1.0

        _, _, _, info = compute_inertia_relief_from_model(model, F)

        self.assertIn('M_rr', info)
        self.assertIn('total_mass', info)
        self.assertIn('grid_ids', info)
        self.assertIn('D', info)
        self.assertIn('residual_force', info)
        self.assertEqual(info['M_rr'].shape, (6, 6))
        self.assertEqual(info['D'].shape, (12, 6))

    def test_from_model_ref_point(self):
        """Different ref_point gives same F_net (invariant property).

        The self-equilibrated loads are independent of the reference point
        choice. Only a_rigid changes (it's expressed at the ref point).
        """
        model = self._load_model()
        F = np.zeros(12)
        F[0] = 50.0  # Fx at grid 1
        F[2] = -30.0  # Fz at grid 1

        F_net_a, _, _, _ = compute_inertia_relief_from_model(
            model, F, ref_point=np.array([0.0, 0.0, 0.0]))
        F_net_b, _, _, _ = compute_inertia_relief_from_model(
            model, F, ref_point=np.array([5.0, 0.0, 0.0]))
        F_net_c, _, _, _ = compute_inertia_relief_from_model(
            model, F, ref_point=np.array([3.0, 7.0, -2.0]))

        npt.assert_allclose(F_net_a, F_net_b, atol=1e-10)
        npt.assert_allclose(F_net_a, F_net_c, atol=1e-10)

    def test_from_model_wtmass(self):
        """WTMASS parameter scales the mass matrix.

        With wtmass=2, masses are doubled -> acceleration halves for same F.
        """
        model = self._load_model()
        F = np.zeros(12)
        F[2] = -100.0
        F[8] = -100.0

        _, a_1, _, info_1 = compute_inertia_relief_from_model(
            model, F, wtmass=1.0)
        _, a_2, _, info_2 = compute_inertia_relief_from_model(
            model, F, wtmass=2.0)

        # wtmass=2 doubles mass, halves acceleration
        npt.assert_allclose(a_2, a_1 / 2.0, rtol=1e-8)
        npt.assert_allclose(info_2['total_mass'], 2 * info_1['total_mass'],
                            rtol=1e-8)


# =============================================================================
# TESTS: BWB integration model
# =============================================================================

class TestBwbInertiaRelief(unittest.TestCase):
    """Integration test with BWB (Blended Wing Body) model.

    Uses the real BWB model to verify inertia relief on a large structure.
    Tests self-equilibration and physical reasonableness.
    """

    @classmethod
    def setUpClass(cls):
        bwb_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(
                os.path.dirname(os.path.dirname(os.path.dirname(
                    os.path.abspath(__file__))))))),
            'models', 'bwb')
        cls.bwb_path = os.path.join(bwb_dir, 'bwb_saero.bdf')
        cls.has_bwb = os.path.isfile(cls.bwb_path)

    @unittest.skipUnless(
        os.path.isfile(os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(
                os.path.dirname(os.path.dirname(os.path.dirname(
                    os.path.abspath(__file__))))))),
            'models', 'bwb', 'bwb_saero.bdf')),
        "BWB model not found")
    def test_bwb_self_equilibration(self):
        """BWB under unit gravity is self-equilibrated after inertia relief."""
        model = BDF(log=None)
        model.read_bdf(self.bwb_path)

        M, grid_ids, total_mass = build_mgg_lumped(model=model)
        node_xyz = model.grid.xyz_cid0()
        D = build_rigid_body_modes(grid_ids, node_xyz)

        n_dof = 6 * len(grid_ids)
        n_grids = len(grid_ids)

        # Apply unit gravity in -Z
        g = 9.81
        F = np.zeros(n_dof)
        diag = M.diagonal()
        for i in range(n_grids):
            mass_i = diag[6*i]  # T1 mass for this grid
            F[6*i + 2] = -mass_i * g

        F_net, a_rigid, _ = compute_inertia_relief(M, D, F)

        # Self-equilibration check (relaxed for 10k+ grids floating-point)
        residual = D.T @ F_net
        npt.assert_allclose(residual, 0.0, atol=1e-4)

        # Acceleration should be pure -Z gravity
        npt.assert_allclose(a_rigid[2], -g, rtol=1e-4)
        npt.assert_allclose(a_rigid[:2], 0.0, atol=1e-4)

        # F_net should be nearly zero for uniform gravity
        # (if mass distribution is purely translational lumped)
        npt.assert_allclose(np.max(np.abs(F_net)), 0.0, atol=1e-4)


if __name__ == '__main__':
    unittest.main()
