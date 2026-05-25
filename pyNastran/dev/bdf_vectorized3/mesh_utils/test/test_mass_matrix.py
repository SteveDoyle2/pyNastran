"""Unit tests for vectorized lumped mass matrix (MGG) assembly.

Tests build_mgg_lumped and get_grid_point_weight using bdf_vectorized3.
Verifies:
- Lumped mass distribution from structural elements (CBAR, CQUAD4)
- CONM2 6x6 mass matrix (with offset and inertia)
- CMASS2 scalar mass elements
- Grid point weight computation (rigid body mass)
- Total mass conservation
- Sparse matrix structure

Tolerances:
- Lumped mass values: atol=1e-10 (floating-point arithmetic)
- Total mass: rtol=1e-8 (sum over many elements)
- CONM2 inertia coupling: atol=1e-10 (exact formula)
- Grid point weight: atol=1e-8 (D^T M D multiplication)
"""
import os
import tempfile
import unittest

import numpy as np
import numpy.testing as npt
from scipy import sparse

from pyNastran.dev.bdf_vectorized3.bdf import read_bdf
from pyNastran.dev.bdf_vectorized3.mesh_utils.mass_matrix import (
    build_mgg_lumped,
    get_grid_point_weight,
)


# =============================================================================
# BDF FIXTURES
# =============================================================================

BDF_TWO_NODE_CBAR = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0
GRID    2               10.0    0.0     0.0
CBAR    1       1       1       2       0.0     1.0     0.0
PBAR    1       1       1.0     1.0     1.0     1.0
MAT1    1       1.0+7           0.3     0.01
ENDDATA
"""

BDF_CONM2_SIMPLE = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0
GRID    2               10.0    0.0     0.0
CONM2   1       1       0       5.0     0.0     0.0     0.0
ENDDATA
"""

BDF_CONM2_OFFSET = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0
GRID    2               5.0     0.0     0.0
CONM2   1       1       0       2.0     3.0     4.0     0.0     +C1
+C1     10.0    0.0     20.0    0.0     0.0     30.0
ENDDATA
"""

BDF_CMASS2 = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0
GRID    2               5.0     0.0     0.0
CMASS2  1       7.5     1       3
ENDDATA
"""

BDF_CQUAD4 = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0
GRID    2               2.0     0.0     0.0
GRID    3               2.0     3.0     0.0
GRID    4               0.0     3.0     0.0
CQUAD4  1       1       1       2       3       4
PSHELL  1       1       0.1
MAT1    1       1.0+7           0.3     2.0
ENDDATA
"""

BDF_MIXED = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0
GRID    2               5.0     0.0     0.0
GRID    3               2.5     3.0     0.0
CBAR    1       1       1       2       0.0     1.0     0.0
PBAR    1       1       0.5     1.0     1.0     1.0
MAT1    1       1.0+7           0.3     0.1
CONM2   1       3       0       10.0    0.0     0.0     0.0
ENDDATA
"""


def write_temp_bdf(content: str) -> str:
    """Write BDF content to a temporary file and return path."""
    fd, path = tempfile.mkstemp(suffix='.bdf', prefix='test_mgg_')
    with os.fdopen(fd, 'w', encoding='utf-8') as f:
        f.write(content)
    return path


# =============================================================================
# TESTS
# =============================================================================

class TestBuildMggLumpedCbar(unittest.TestCase):
    """Lumped mass from CBAR elements."""

    @classmethod
    def setUpClass(cls):
        cls.bdf_path = write_temp_bdf(BDF_TWO_NODE_CBAR)
        cls.M, cls.grid_ids, cls.total_mass = build_mgg_lumped(cls.bdf_path)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.bdf_path)

    def test_matrix_shape(self):
        """MGG shape is (12, 12) for 2 grids x 6 DOFs."""
        self.assertEqual(self.M.shape, (12, 12))

    def test_grid_ids_sorted(self):
        """Grid IDs returned in sorted order."""
        npt.assert_array_equal(self.grid_ids, [1, 2])

    def test_sparse_type(self):
        """Result is CSR sparse matrix."""
        self.assertTrue(sparse.issparse(self.M))
        self.assertEqual(self.M.format, 'csr')

    def test_total_mass(self):
        """Total mass = rho * A * L = 0.01 * 1.0 * 10.0 = 0.1.

        CBAR length=10, area=1.0, density=0.01.
        """
        expected_mass = 0.01 * 1.0 * 10.0
        npt.assert_allclose(self.total_mass, expected_mass, rtol=1e-8)

    def test_lumped_distribution(self):
        """Mass evenly split: 0.1 / 2 = 0.05 per node on T1, T2, T3."""
        mass_per_node = 0.1 / 2
        diag = self.M.diagonal()
        # Node 1: DOFs 0,1,2
        npt.assert_allclose(diag[0], mass_per_node, atol=1e-10)
        npt.assert_allclose(diag[1], mass_per_node, atol=1e-10)
        npt.assert_allclose(diag[2], mass_per_node, atol=1e-10)
        # Node 2: DOFs 6,7,8
        npt.assert_allclose(diag[6], mass_per_node, atol=1e-10)
        npt.assert_allclose(diag[7], mass_per_node, atol=1e-10)
        npt.assert_allclose(diag[8], mass_per_node, atol=1e-10)

    def test_rotational_dofs_zero(self):
        """Rotational DOFs get zero mass from lumped bar."""
        diag = self.M.diagonal()
        npt.assert_allclose(diag[3:6], 0.0, atol=1e-14)
        npt.assert_allclose(diag[9:12], 0.0, atol=1e-14)

    def test_diagonal_dominant(self):
        """Lumped mass is purely diagonal (no off-diagonal from bars)."""
        M_dense = self.M.toarray()
        off_diag = M_dense - np.diag(np.diag(M_dense))
        npt.assert_allclose(off_diag, 0.0, atol=1e-14)


class TestBuildMggLumpedCquad4(unittest.TestCase):
    """Lumped mass from a CQUAD4 shell element."""

    @classmethod
    def setUpClass(cls):
        cls.bdf_path = write_temp_bdf(BDF_CQUAD4)
        cls.M, cls.grid_ids, cls.total_mass = build_mgg_lumped(cls.bdf_path)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.bdf_path)

    def test_total_mass(self):
        """Element mass = area * t * rho = 6.0 * 0.1 * 2.0 = 1.2."""
        npt.assert_allclose(self.total_mass, 1.2, rtol=1e-6)

    def test_lumped_per_node(self):
        """1.2 / 4 nodes = 0.3 per node on T1, T2, T3."""
        mass_per_node = 1.2 / 4.0
        diag = self.M.diagonal()
        for node_idx in range(4):
            base = 6 * node_idx
            for dof in range(3):
                npt.assert_allclose(
                    diag[base + dof], mass_per_node, rtol=1e-6,
                    err_msg=f"Node {node_idx+1} DOF {dof+1}")

    def test_shape(self):
        """4 grids -> (24, 24) matrix."""
        self.assertEqual(self.M.shape, (24, 24))


class TestConm2MassMatrix(unittest.TestCase):
    """Tests for CONM2 mass assembly."""

    def test_simple_no_offset(self):
        """CONM2 mass=5.0, no offset: diagonal [5,5,5,0,0,0]."""
        bdf_path = write_temp_bdf(BDF_CONM2_SIMPLE)
        try:
            M, grid_ids, total_mass = build_mgg_lumped(bdf_path)
        finally:
            os.unlink(bdf_path)

        npt.assert_allclose(total_mass, 5.0, atol=1e-10)
        diag = M.diagonal()
        npt.assert_allclose(diag[0:3], 5.0, atol=1e-10)
        npt.assert_allclose(diag[3:6], 0.0, atol=1e-10)
        npt.assert_allclose(diag[6:12], 0.0, atol=1e-14)

    def test_offset_coupling(self):
        """CONM2 with offset=[3,4,0] produces off-diagonal coupling.

        mass=2, offset=[3, 4, 0]:
        M[0,5] = -m*x2 = -8, M[1,5] = m*x1 = 6
        M[2,3] = m*x2 = 8,   M[2,4] = -m*x1 = -6
        """
        bdf_path = write_temp_bdf(BDF_CONM2_OFFSET)
        try:
            M, _, total_mass = build_mgg_lumped(bdf_path)
        finally:
            os.unlink(bdf_path)

        npt.assert_allclose(total_mass, 2.0, atol=1e-10)
        M_dense = M.toarray()
        m = 2.0
        x1, x2, x3 = 3.0, 4.0, 0.0

        npt.assert_allclose(M_dense[0, 0], m, atol=1e-10)
        npt.assert_allclose(M_dense[1, 1], m, atol=1e-10)
        npt.assert_allclose(M_dense[2, 2], m, atol=1e-10)

        npt.assert_allclose(M_dense[0, 4], m * x3, atol=1e-10)
        npt.assert_allclose(M_dense[0, 5], -m * x2, atol=1e-10)
        npt.assert_allclose(M_dense[1, 3], -m * x3, atol=1e-10)
        npt.assert_allclose(M_dense[1, 5], m * x1, atol=1e-10)
        npt.assert_allclose(M_dense[2, 3], m * x2, atol=1e-10)
        npt.assert_allclose(M_dense[2, 4], -m * x1, atol=1e-10)

    def test_offset_inertia(self):
        """CONM2 rotational inertia includes parallel axis theorem.

        mass=2, offset=[3,4,0], I11=10, I22=20, I33=30
        I_grid[3,3] = I11 + m*(x2^2+x3^2) = 10 + 2*16 = 42
        I_grid[4,4] = I22 + m*(x1^2+x3^2) = 20 + 2*9 = 38
        I_grid[5,5] = I33 + m*(x1^2+x2^2) = 30 + 2*25 = 80
        """
        bdf_path = write_temp_bdf(BDF_CONM2_OFFSET)
        try:
            M, _, _ = build_mgg_lumped(bdf_path)
        finally:
            os.unlink(bdf_path)

        M_dense = M.toarray()
        npt.assert_allclose(M_dense[3, 3], 42.0, atol=1e-10)
        npt.assert_allclose(M_dense[4, 4], 38.0, atol=1e-10)
        npt.assert_allclose(M_dense[5, 5], 80.0, atol=1e-10)

    def test_conm2_symmetric(self):
        """CONM2 6x6 mass block is symmetric."""
        bdf_path = write_temp_bdf(BDF_CONM2_OFFSET)
        try:
            M, _, _ = build_mgg_lumped(bdf_path)
        finally:
            os.unlink(bdf_path)

        M_dense = M.toarray()
        M66 = M_dense[:6, :6]
        npt.assert_allclose(M66, M66.T, atol=1e-10)


class TestCmass(unittest.TestCase):
    """Tests for CMASS2 scalar mass."""

    def test_cmass2_single_dof(self):
        """CMASS2 mass=7.5 on grid 1, DOF 3 (T3).

        Only diagonal[2] should be 7.5.
        total_mass = (0 + 0 + 7.5) / 3 = 2.5
        """
        bdf_path = write_temp_bdf(BDF_CMASS2)
        try:
            M, grid_ids, total_mass = build_mgg_lumped(bdf_path)
        finally:
            os.unlink(bdf_path)

        diag = M.diagonal()
        npt.assert_allclose(diag[2], 7.5, atol=1e-10)
        npt.assert_allclose(diag[0], 0.0, atol=1e-14)
        npt.assert_allclose(diag[1], 0.0, atol=1e-14)
        npt.assert_allclose(total_mass, 2.5, atol=1e-10)


class TestMixedModel(unittest.TestCase):
    """Test with both structural elements and concentrated mass."""

    @classmethod
    def setUpClass(cls):
        cls.bdf_path = write_temp_bdf(BDF_MIXED)
        cls.M, cls.grid_ids, cls.total_mass = build_mgg_lumped(cls.bdf_path)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.bdf_path)

    def test_total_mass(self):
        """Total = CBAR(0.25) + CONM2(10.0) = 10.25.

        CBAR: L=5, A=0.5, rho=0.1 -> mass = 0.25
        """
        npt.assert_allclose(self.total_mass, 10.25, rtol=1e-6)

    def test_conm2_at_grid3(self):
        """CONM2 mass=10 at grid 3 (index 2, base DOF=12)."""
        diag = self.M.diagonal()
        npt.assert_allclose(diag[12], 10.0, atol=1e-10)
        npt.assert_allclose(diag[13], 10.0, atol=1e-10)
        npt.assert_allclose(diag[14], 10.0, atol=1e-10)

    def test_bar_nodes_get_lumped(self):
        """Nodes 1,2 get CBAR lumped: 0.25/2 = 0.125 each."""
        diag = self.M.diagonal()
        mass_per_node = 0.25 / 2.0
        npt.assert_allclose(diag[0:3], mass_per_node, rtol=1e-6)
        npt.assert_allclose(diag[6:9], mass_per_node, rtol=1e-6)


class TestWtmass(unittest.TestCase):
    """Test WTMASS parameter scaling."""

    def test_wtmass_scales_mass(self):
        """WTMASS=0.5 halves all mass entries."""
        bdf_path = write_temp_bdf(BDF_CONM2_SIMPLE)
        try:
            M1, _, mass1 = build_mgg_lumped(bdf_path, wtmass=1.0)
            M2, _, mass2 = build_mgg_lumped(bdf_path, wtmass=0.5)
        finally:
            os.unlink(bdf_path)

        npt.assert_allclose(mass2, mass1 * 0.5, atol=1e-10)
        npt.assert_allclose(M2.diagonal(), M1.diagonal() * 0.5, atol=1e-10)


class TestGetGridPointWeight(unittest.TestCase):
    """Tests for get_grid_point_weight (D^T M D)."""

    def test_single_point_mass_at_origin(self):
        """m=3 at origin: M0 = diag([3,3,3,0,0,0])."""
        grid_ids = np.array([1])
        node_xyz = np.array([[0.0, 0.0, 0.0]])
        M = sparse.diags([3., 3., 3., 0., 0., 0.], format='csr')

        M0 = get_grid_point_weight(M, grid_ids, node_xyz)
        expected = np.diag([3., 3., 3., 0., 0., 0.])
        npt.assert_allclose(M0, expected, atol=1e-10)

    def test_point_mass_offset_from_ref(self):
        """m=4 at (5,0,0), ref at origin: Iyy=Izz=100."""
        grid_ids = np.array([1])
        node_xyz = np.array([[5.0, 0.0, 0.0]])
        M = sparse.diags([4., 4., 4., 0., 0., 0.], format='csr')

        M0 = get_grid_point_weight(M, grid_ids, node_xyz, ref_point=np.zeros(3))
        npt.assert_allclose(M0[0, 0], 4.0, atol=1e-10)
        npt.assert_allclose(M0[3, 3], 0.0, atol=1e-10)
        npt.assert_allclose(M0[4, 4], 100.0, atol=1e-10)
        npt.assert_allclose(M0[5, 5], 100.0, atol=1e-10)

    def test_two_masses_symmetric(self):
        """Two m=1 at (+/-3, 0, 0): Iyy=Izz=18, symmetric."""
        grid_ids = np.array([1, 2])
        node_xyz = np.array([[3.0, 0.0, 0.0], [-3.0, 0.0, 0.0]])
        diag_vals = np.array([1., 1., 1., 0., 0., 0., 1., 1., 1., 0., 0., 0.])
        M = sparse.diags(diag_vals, format='csr')

        M0 = get_grid_point_weight(M, grid_ids, node_xyz)
        npt.assert_allclose(M0, M0.T, atol=1e-12)
        npt.assert_allclose(M0[0, 0], 2.0, atol=1e-10)
        npt.assert_allclose(M0[3, 3], 0.0, atol=1e-10)
        npt.assert_allclose(M0[4, 4], 18.0, atol=1e-10)
        npt.assert_allclose(M0[5, 5], 18.0, atol=1e-10)

    def test_ref_point_shifts_inertia(self):
        """m=1 at (2,0,0): ref at origin Iyy=4, ref at mass Iyy=0."""
        grid_ids = np.array([1])
        node_xyz = np.array([[2.0, 0.0, 0.0]])
        M = sparse.diags([1., 1., 1., 0., 0., 0.], format='csr')

        M0_origin = get_grid_point_weight(M, grid_ids, node_xyz,
                                          ref_point=np.array([0., 0., 0.]))
        M0_at_mass = get_grid_point_weight(M, grid_ids, node_xyz,
                                           ref_point=np.array([2., 0., 0.]))

        npt.assert_allclose(M0_origin[4, 4], 4.0, atol=1e-10)
        npt.assert_allclose(M0_at_mass[4, 4], 0.0, atol=1e-10)


class TestCoordinateTransforms(unittest.TestCase):
    """Tests for CID rotation and Mbb-to-Mgg transform."""

    def test_output_frame_basic_is_default(self):
        """output_frame='basic' is the default and matches previous behavior."""
        bdf_path = write_temp_bdf(BDF_CONM2_SIMPLE)
        try:
            M_basic, _, _ = build_mgg_lumped(bdf_path, output_frame='basic')
            M_default, _, _ = build_mgg_lumped(bdf_path)
        finally:
            os.unlink(bdf_path)

        npt.assert_allclose(M_basic.toarray(), M_default.toarray(), atol=1e-14)

    def test_global_frame_cd0_equals_basic(self):
        """When all grids have CD=0, Mgg = Mbb (no rotation applied)."""
        bdf_path = write_temp_bdf(BDF_CONM2_SIMPLE)
        try:
            M_basic, _, mass_b = build_mgg_lumped(bdf_path, output_frame='basic')
            M_global, _, mass_g = build_mgg_lumped(bdf_path, output_frame='global')
        finally:
            os.unlink(bdf_path)

        npt.assert_allclose(M_global.toarray(), M_basic.toarray(), atol=1e-14)
        npt.assert_allclose(mass_g, mass_b, atol=1e-14)

    def test_total_mass_invariant_under_cd_rotation(self):
        """Total mass doesn't change when applying Mbb->Mgg transform.

        CD rotation is orthogonal: trace(T^T M T) for diagonal blocks
        preserves the sum m*(R^T@I@R)[i,i] = m (since R^T R = I for m*I3).
        """
        # Use a model with non-trivial content but all CD=0 for comparison
        bdf_path = write_temp_bdf(BDF_MIXED)
        try:
            _, _, mass_b = build_mgg_lumped(bdf_path, output_frame='basic')
            _, _, mass_g = build_mgg_lumped(bdf_path, output_frame='global')
        finally:
            os.unlink(bdf_path)

        npt.assert_allclose(mass_g, mass_b, rtol=1e-10)

    def test_conm2_cid0_no_rotation(self):
        """CONM2 with CID=0: offset used directly in basic frame.

        mass=2, offset=[3,4,0], CID=0 -> same as before.
        """
        bdf_path = write_temp_bdf(BDF_CONM2_OFFSET)
        try:
            M, _, _ = build_mgg_lumped(bdf_path)
        finally:
            os.unlink(bdf_path)

        M_dense = M.toarray()
        # Coupling M[0,5] = -m*x2 = -8
        npt.assert_allclose(M_dense[0, 5], -8.0, atol=1e-10)
        # Coupling M[1,5] = m*x1 = 6
        npt.assert_allclose(M_dense[1, 5], 6.0, atol=1e-10)

    def test_conm2_cid_rotation(self):
        """CONM2 with non-zero CID rotates offset to basic.

        CID=1 rotates 90 deg about Z: local X -> basic Y, local Y -> basic -X.
        Offset in CID=1: [1, 0, 0] -> basic: [0, 1, 0]
        mass=3, offset_basic=[0, 1, 0]

        Expected coupling in basic:
            M[0,5] = -m*x2_basic = -3*1 = -3
            M[1,5] = m*x1_basic = 3*0 = 0
            M[2,3] = m*x2_basic = 3*1 = 3
        """
        bdf_content = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0
GRID    2               5.0     0.0     0.0
CORD2R  1       0       0.0     0.0     0.0     0.0     0.0     1.0     +CR1
+CR1    0.0     1.0     0.0
CONM2   1       1       1       3.0     1.0     0.0     0.0
ENDDATA
"""
        bdf_path = write_temp_bdf(bdf_content)
        try:
            M, _, total_mass = build_mgg_lumped(bdf_path, output_frame='basic')
        finally:
            os.unlink(bdf_path)

        npt.assert_allclose(total_mass, 3.0, atol=1e-10)
        M_dense = M.toarray()

        # offset in basic = [0, 1, 0] after 90-deg Z rotation
        # M[0,5] = -m*y_basic = -3*1 = -3
        npt.assert_allclose(M_dense[0, 5], -3.0, atol=1e-10)
        # M[1,5] = m*x_basic = 3*0 = 0
        npt.assert_allclose(M_dense[1, 5], 0.0, atol=1e-10)
        # M[2,3] = m*y_basic = 3*1 = 3
        npt.assert_allclose(M_dense[2, 3], 3.0, atol=1e-10)

    def test_mbb_to_mgg_with_cd(self):
        """Grid with CD != 0: Mgg rotates DOFs into CD frame.

        Grid 1 at origin with CD=1 (90-deg rotation about Z).
        CONM2 mass=1, no offset, no inertia -> Mbb = diag([1,1,1,0,0,0]).

        In Mgg: T_1 = R^T for the 3x3 translational block.
        Since m*I3 is rotation-invariant, Mgg diagonal should still be
        [1,1,1,0,0,0] — isotropic mass is unchanged by rotation.
        """
        bdf_content = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0     1
GRID    2               5.0     0.0     0.0
CORD2R  1       0       0.0     0.0     0.0     0.0     0.0     1.0     +CR1
+CR1    0.0     1.0     0.0
CONM2   1       1       0       1.0     0.0     0.0     0.0
ENDDATA
"""
        bdf_path = write_temp_bdf(bdf_content)
        try:
            M_bb, _, _ = build_mgg_lumped(bdf_path, output_frame='basic')
            M_gg, _, _ = build_mgg_lumped(bdf_path, output_frame='global')
        finally:
            os.unlink(bdf_path)

        # Isotropic mass: diagonal unchanged by rotation
        diag_bb = M_bb.diagonal()
        diag_gg = M_gg.diagonal()
        npt.assert_allclose(diag_gg[:3], diag_bb[:3], atol=1e-10)

    def test_mgg_with_anisotropic_conm2(self):
        """CONM2 with inertia: Mgg rotates the 6x6 block.

        Grid 1 with CD=1 (90 deg Z rotation).
        CONM2: mass=1, offset=[0,0,0], I11=10, I22=20, I33=30.
        In basic frame: Mbb[3,3]=10, Mbb[4,4]=20, Mbb[5,5]=30.

        CD=1: local X = basic Y, local Y = -basic X, local Z = basic Z.
        Mgg = R @ Mbb @ R^T, so:
            Mgg[3,3] = Ixx_local = Iyy_basic = 20
            Mgg[4,4] = Iyy_local = Ixx_basic = 10
            Mgg[5,5] = Izz_local = Izz_basic = 30
        """
        bdf_content = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0     1
GRID    2               5.0     0.0     0.0
CORD2R  1       0       0.0     0.0     0.0     0.0     0.0     1.0     +CR1
+CR1    0.0     1.0     0.0
CONM2   1       1       0       1.0     0.0     0.0     0.0     +C1
+C1     10.0    0.0     20.0    0.0     0.0     30.0
ENDDATA
"""
        bdf_path = write_temp_bdf(bdf_content)
        try:
            M_bb, _, _ = build_mgg_lumped(bdf_path, output_frame='basic')
            M_gg, _, _ = build_mgg_lumped(bdf_path, output_frame='global')
        finally:
            os.unlink(bdf_path)

        # Basic frame: I11=10, I22=20, I33=30
        npt.assert_allclose(M_bb.toarray()[3, 3], 10.0, atol=1e-10)
        npt.assert_allclose(M_bb.toarray()[4, 4], 20.0, atol=1e-10)
        npt.assert_allclose(M_bb.toarray()[5, 5], 30.0, atol=1e-10)

        # Global frame: rotated by CD=1 (90 deg Z)
        # R rotates basic->local: local_x=basic_y, local_y=-basic_x
        # Mgg = T^T @ Mbb @ T, where T = R^T (basic-to-local)
        # For inertia: I_local = R^T @ I_basic @ R
        # With 90-deg Z: Ixx_local = Iyy_basic = 20
        #                Iyy_local = Ixx_basic = 10
        #                Izz_local = Izz_basic = 30
        npt.assert_allclose(M_gg.toarray()[3, 3], 20.0, atol=1e-10)
        npt.assert_allclose(M_gg.toarray()[4, 4], 10.0, atol=1e-10)
        npt.assert_allclose(M_gg.toarray()[5, 5], 30.0, atol=1e-10)

    def test_invalid_output_frame(self):
        """Invalid output_frame raises ValueError."""
        bdf_path = write_temp_bdf(BDF_CONM2_SIMPLE)
        try:
            with self.assertRaises(ValueError):
                build_mgg_lumped(bdf_path, output_frame='invalid')
        finally:
            os.unlink(bdf_path)


class TestBwbModel(unittest.TestCase):
    """Integration test on the BWB model."""

    BWB_PATH = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))),
        'models', 'bwb', 'bwb_saero.bdf')

    @classmethod
    def setUpClass(cls):
        if not os.path.exists(cls.BWB_PATH):
            raise unittest.SkipTest("BWB model not found")
        cls.M, cls.grid_ids, cls.total_mass = build_mgg_lumped(cls.BWB_PATH)

    def test_matrix_shape(self):
        """MGG shape = (6*n_grids, 6*n_grids)."""
        n_grids = len(self.grid_ids)
        self.assertEqual(self.M.shape, (6 * n_grids, 6 * n_grids))

    def test_total_mass_positive(self):
        """Total mass > 0."""
        self.assertGreater(self.total_mass, 0.0)

    def test_total_mass_reasonable(self):
        """BWB mass should be > 1 (thousands of kg expected)."""
        self.assertGreater(self.total_mass, 1.0)

    def test_diagonal_non_negative(self):
        """All diagonal entries >= 0."""
        diag = self.M.diagonal()
        self.assertTrue(np.all(diag >= -1e-14),
                        f"Negative diagonal: min={diag.min()}")

    def test_sparse_efficiency(self):
        """Matrix density < 1% (lumped mass is highly sparse)."""
        n_dof = self.M.shape[0]
        density = self.M.nnz / (n_dof * n_dof)
        self.assertLess(density, 0.01)

    def test_grid_point_weight(self):
        """Grid point weight is symmetric with correct total mass."""
        model = read_bdf(self.BWB_PATH, xref=True, debug=False)
        node_xyz = model.grid.xyz_cid0()

        M0 = get_grid_point_weight(self.M, self.grid_ids, node_xyz)

        npt.assert_allclose(M0, M0.T, atol=1e-6)
        npt.assert_allclose(M0[0, 0], self.total_mass, rtol=1e-6)
        npt.assert_allclose(M0[1, 1], self.total_mass, rtol=1e-6)
        npt.assert_allclose(M0[2, 2], self.total_mass, rtol=1e-6)
        self.assertGreater(M0[3, 3], 0.0)
        self.assertGreater(M0[4, 4], 0.0)
        self.assertGreater(M0[5, 5], 0.0)


class TestTransformVectorized(unittest.TestCase):
    """Tests that the vectorized Mbb->Mgg transform handles multiple CIDs."""

    def test_multiple_cids(self):
        """Multiple grids with different CDs produce correct rotated mass.

        Grid 1: CD=0, CONM2 mass=1, I11=5, I22=10, I33=15
        Grid 2: CD=1 (90-deg Z), CONM2 mass=2, no inertia
        Grid 3: CD=2 (90-deg Y), CONM2 mass=3, no inertia

        Isotropic translational mass unchanged by rotation.
        """
        bdf_content = """\
SOL 101
CEND
BEGIN BULK
GRID    1               0.0     0.0     0.0     0
GRID    2               5.0     0.0     0.0     1
GRID    3               0.0     5.0     0.0     2
CORD2R  1       0       0.0     0.0     0.0     0.0     0.0     1.0     +CR1
+CR1    0.0     1.0     0.0
CORD2R  2       0       0.0     0.0     0.0     0.0     1.0     0.0     +CR2
+CR2    0.0     0.0     1.0
CONM2   1       1       0       1.0     0.0     0.0     0.0     +C1
+C1     5.0     0.0     10.0    0.0     0.0     15.0
CONM2   2       2       0       2.0     0.0     0.0     0.0
CONM2   3       3       0       3.0     0.0     0.0     0.0
ENDDATA
"""
        bdf_path = write_temp_bdf(bdf_content)
        try:
            M_bb, _, mass_b = build_mgg_lumped(bdf_path, output_frame='basic')
            M_gg, _, mass_g = build_mgg_lumped(bdf_path, output_frame='global')
        finally:
            os.unlink(bdf_path)

        # Total mass invariant
        npt.assert_allclose(mass_g, mass_b, rtol=1e-10)
        npt.assert_allclose(mass_b, 6.0, atol=1e-10)

        # Translational diag unchanged (isotropic)
        diag_bb = M_bb.diagonal()
        diag_gg = M_gg.diagonal()
        for i in range(3):
            base = 6 * i
            npt.assert_allclose(diag_gg[base:base+3], diag_bb[base:base+3],
                                atol=1e-10)

        # Grid 1 (CD=0): inertia unchanged
        npt.assert_allclose(diag_gg[3], 5.0, atol=1e-10)
        npt.assert_allclose(diag_gg[4], 10.0, atol=1e-10)
        npt.assert_allclose(diag_gg[5], 15.0, atol=1e-10)


class TestPerformanceLargeModel(unittest.TestCase):
    """Performance regression test with synthetic large model.

    Verifies correctness and that assembly completes in reasonable time
    for 10000+ element models. Not a strict timing assertion — just
    ensures the vectorized paths produce correct results at scale.

    Tolerances: total mass rtol=1e-8, per-node mass atol=1e-10.
    Interactions: CQUAD4 lumped mass distribution, searchsorted indexing.
    """

    @classmethod
    def setUpClass(cls):
        nx, ny = 50, 50
        lines = ['SOL 101', 'CEND', 'BEGIN BULK']
        nid = 1
        for j in range(ny):
            for i in range(nx):
                lines.append(
                    f'GRID    {nid:8d}        '
                    f'{float(i):8.1f}{float(j):8.1f}     0.0')
                nid += 1
        eid = 1
        for j in range(ny - 1):
            for i in range(nx - 1):
                n1 = j * nx + i + 1
                n2 = n1 + 1
                n3 = n2 + nx
                n4 = n1 + nx
                lines.append(
                    f'CQUAD4  {eid:8d}       1'
                    f'{n1:8d}{n2:8d}{n3:8d}{n4:8d}')
                eid += 1
        lines.append('PSHELL  1       1       0.1')
        lines.append('MAT1    1       1.0+7           0.3     2.7')
        lines.append('ENDDATA')

        cls.bdf_path = write_temp_bdf('\n'.join(lines))
        cls.n_grids = nx * ny
        cls.n_elem = (nx - 1) * (ny - 1)
        cls.expected_mass = 1.0 * 1.0 * 0.1 * 2.7 * cls.n_elem

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.bdf_path)

    def test_total_mass(self):
        """Total mass = n_elem * area * thickness * density."""
        M, grid_ids, total_mass = build_mgg_lumped(self.bdf_path)
        npt.assert_allclose(total_mass, self.expected_mass, rtol=1e-8)

    def test_shape(self):
        """Matrix shape matches 6 DOFs per grid."""
        M, grid_ids, _ = build_mgg_lumped(self.bdf_path)
        self.assertEqual(M.shape, (6 * self.n_grids, 6 * self.n_grids))
        self.assertEqual(len(grid_ids), self.n_grids)

    def test_diagonal_non_negative(self):
        """All diagonal entries non-negative."""
        M, _, _ = build_mgg_lumped(self.bdf_path)
        self.assertTrue(np.all(M.diagonal() >= -1e-14))


if __name__ == '__main__':
    unittest.main()
