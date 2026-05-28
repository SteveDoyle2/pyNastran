"""Tests for spline_methods.py — aerodynamic-structural interpolation.

Tests the core spline math with raw numpy arrays (no BDF/dataclass coupling).
Validates fundamental properties:
- Partition of unity (rows sum to 1 => rigid body motion preserved)
- Linear field reproduction (polynomial augmentation guarantees)
- Correct matrix shapes
- Coincident point handling
- Singular-matrix fallback path (lstsq)

Tolerances:
- Partition of unity: atol=1e-12
- Linear field reproduction: atol=1e-10
- Coincident point: atol=1e-12
"""

import numpy as np
import numpy.testing as npt
import pytest

from pyNastran.dev.bdf_vectorized3.cards.aero.spline_methods import (
    ips_spline,
    tps_spline,
    fps_spline,
    ris_spline,
    beam_spline,
    rbf_spline,
    nearest_spline,
    clip_to_convex_hull,
    validate_spline,
)


# =============================================================================
# FIXTURES
# =============================================================================


def _rect_grid(
    nx: int, ny: int, spacing: float = 1.0, offset: tuple = (0.0, 0.0, 0.0)
) -> np.ndarray:
    """Generate rectangular grid of (N, 3) points."""
    pts = []
    for iy in range(ny):
        for ix in range(nx):
            pts.append([offset[0] + ix * spacing, offset[1] + iy * spacing, offset[2]])
    return np.array(pts)


@pytest.fixture
def struct_pts():
    """4x3 = 12 structural points, spacing=1.0, at origin."""
    return _rect_grid(4, 3, spacing=1.0)


@pytest.fixture
def aero_pts():
    """5x4 = 20 aero points, spacing=0.7, offset (0.25, 0.25)."""
    return _rect_grid(5, 4, spacing=0.7, offset=(0.25, 0.25, 0.0))


# =============================================================================
# IPS (Infinite Plate Spline)
# =============================================================================


class TestIpsSpline:
    """IPS: 2D projected thin plate spline (SPLINE1 METHOD=IPS)."""

    def test_shape(self, struct_pts, aero_pts):
        """G is (M_aero, N_struct)."""
        G = ips_spline(struct_pts, aero_pts)
        assert G.shape == (20, 12)

    def test_partition_of_unity(self, struct_pts, aero_pts):
        """Rows sum to 1 (constant field preserved). atol=1e-12."""
        G = ips_spline(struct_pts, aero_pts)
        npt.assert_allclose(G.sum(axis=1), 1.0, atol=1e-12)

    def test_linear_field_x(self, struct_pts, aero_pts):
        """Linear field u=x reproduced exactly. atol=1e-10."""
        G = ips_spline(struct_pts, aero_pts)
        u_struct = struct_pts[:, 0]
        u_aero = G @ u_struct
        npt.assert_allclose(u_aero, aero_pts[:, 0], atol=1e-10)

    def test_linear_field_y(self, struct_pts, aero_pts):
        """Linear field u=y reproduced exactly. atol=1e-10."""
        G = ips_spline(struct_pts, aero_pts)
        u_struct = struct_pts[:, 1]
        u_aero = G @ u_struct
        npt.assert_allclose(u_aero, aero_pts[:, 1], atol=1e-10)

    def test_coincident_point(self, struct_pts):
        """Aero point at struct location gives unit weight. atol=1e-12."""
        aero_at_5 = struct_pts[5:6, :]
        G = ips_spline(struct_pts, aero_at_5)
        expected = np.zeros(12)
        expected[5] = 1.0
        npt.assert_allclose(G[0], expected, atol=1e-12)


# =============================================================================
# TPS (Thin Plate Spline — 3D)
# =============================================================================


class TestTpsSpline:
    """TPS: 3D distance thin plate spline (SPLINE1 METHOD=TPS)."""

    def test_shape(self, struct_pts, aero_pts):
        G = tps_spline(struct_pts, aero_pts)
        assert G.shape == (20, 12)

    def test_partition_of_unity(self, struct_pts, aero_pts):
        """atol=1e-12."""
        G = tps_spline(struct_pts, aero_pts)
        npt.assert_allclose(G.sum(axis=1), 1.0, atol=1e-12)

    def test_linear_field(self, struct_pts, aero_pts):
        """u = 2x + 3y reproduced exactly. atol=1e-10."""
        G = tps_spline(struct_pts, aero_pts)
        u_struct = 2.0 * struct_pts[:, 0] + 3.0 * struct_pts[:, 1]
        u_aero = G @ u_struct
        u_expected = 2.0 * aero_pts[:, 0] + 3.0 * aero_pts[:, 1]
        npt.assert_allclose(u_aero, u_expected, atol=1e-10)

    def test_coplanar_matches_ips(self, struct_pts, aero_pts):
        """For z=0 coplanar points, TPS ≈ IPS. atol=1e-8."""
        G_ips = ips_spline(struct_pts, aero_pts)
        G_tps = tps_spline(struct_pts, aero_pts)
        npt.assert_allclose(G_tps, G_ips, atol=1e-8)


# =============================================================================
# RBF (Radial Basis Function)
# =============================================================================


class TestRbfSpline:
    """RBF spline with various kernels."""

    def test_thin_plate_matches_ips(self, struct_pts, aero_pts):
        """RBF thin_plate kernel ≈ IPS for coplanar z=0. atol=1e-8."""
        G_ips = ips_spline(struct_pts, aero_pts)
        G_rbf = rbf_spline(struct_pts, aero_pts, rbf="thin_plate")
        npt.assert_allclose(G_rbf, G_ips, atol=1e-8)

    def test_multiquadric_shape(self, struct_pts, aero_pts):
        """Multiquadric produces finite matrix of correct shape."""
        G = rbf_spline(struct_pts, aero_pts, rbf="multiquadric")
        assert G.shape == (20, 12)
        assert np.all(np.isfinite(G))

    def test_partition_of_unity(self, struct_pts, aero_pts):
        """All RBF kernels preserve constant fields. atol=1e-10."""
        for kernel in ("thin_plate", "multiquadric", "linear"):
            G = rbf_spline(struct_pts, aero_pts, rbf=kernel)
            npt.assert_allclose(G.sum(axis=1), 1.0, atol=1e-10, err_msg=f"kernel={kernel}")

    def test_singular_fallback(self):
        """Duplicate struct points trigger lstsq fallback without crashing."""
        struct = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        aero = np.array([[0.5, 0.5, 0.0]])
        G = rbf_spline(struct, aero, rbf="thin_plate")
        assert G.shape == (1, 4)
        assert np.all(np.isfinite(G))


# =============================================================================
# RIS (Radial Interpolation Spline — Wendland compact support)
# =============================================================================


class TestRisSpline:
    """RIS: Wendland compact-support RBF (SPLINE5 METHOD=RIS)."""

    def test_shape_wf0(self, struct_pts, aero_pts):
        G = ris_spline(struct_pts, aero_pts, ftype="WF0", rcore=3.0)
        assert G.shape == (20, 12)

    def test_shape_wf2(self, struct_pts, aero_pts):
        G = ris_spline(struct_pts, aero_pts, ftype="WF2", rcore=3.0)
        assert G.shape == (20, 12)

    def test_partition_of_unity_wf0(self, struct_pts, aero_pts):
        """atol=1e-10."""
        G = ris_spline(struct_pts, aero_pts, ftype="WF0", rcore=5.0)
        npt.assert_allclose(G.sum(axis=1), 1.0, atol=1e-10)

    def test_partition_of_unity_wf2(self, struct_pts, aero_pts):
        """atol=1e-10."""
        G = ris_spline(struct_pts, aero_pts, ftype="WF2", rcore=5.0)
        npt.assert_allclose(G.sum(axis=1), 1.0, atol=1e-10)

    def test_linear_field_wf2(self, struct_pts, aero_pts):
        """WF2 reproduces linear fields. atol=1e-8."""
        G = ris_spline(struct_pts, aero_pts, ftype="WF2", rcore=5.0)
        u_struct = struct_pts[:, 0]
        u_aero = G @ u_struct
        npt.assert_allclose(u_aero, aero_pts[:, 0], atol=1e-8)

    def test_finite_values(self, struct_pts, aero_pts):
        """RIS with various rcore values always produces finite output."""
        for rcore in (1.5, 3.0, 5.0):
            G = ris_spline(struct_pts, aero_pts, ftype="WF0", rcore=rcore)
            assert np.all(np.isfinite(G)), f"rcore={rcore} produced non-finite"

    def test_invalid_ftype_raises(self, struct_pts, aero_pts):
        with pytest.raises((ValueError, KeyError)):
            ris_spline(struct_pts, aero_pts, ftype="INVALID", rcore=3.0)

    def test_auto_rcore(self, struct_pts, aero_pts):
        """rcore=None triggers automatic estimation."""
        G = ris_spline(struct_pts, aero_pts, ftype="WF2", rcore=None)
        assert G.shape == (20, 12)
        assert np.all(np.isfinite(G))


# =============================================================================
# NEAREST-NEIGHBOR
# =============================================================================


class TestNearestSpline:
    """Nearest-neighbor / inverse-distance weighting (SPLINE3-like)."""

    def test_partition_of_unity(self, struct_pts, aero_pts):
        """Normalized weights sum to 1. atol=1e-14."""
        G = nearest_spline(struct_pts, aero_pts, n_nearest=4, weighting="inverse_distance")
        npt.assert_allclose(G.sum(axis=1), 1.0, atol=1e-14)

    def test_n_nearest_1(self, struct_pts, aero_pts):
        """n_nearest=1: exactly one non-zero per row."""
        G = nearest_spline(struct_pts, aero_pts, n_nearest=1, weighting="inverse_distance")
        for i in range(G.shape[0]):
            nonzero = np.count_nonzero(G[i])
            assert nonzero == 1, f"Row {i} has {nonzero} nonzeros, expected 1"

    def test_coincident_gives_unit(self, struct_pts):
        """Aero point at struct location gives weight 1 on that node."""
        aero = struct_pts[3:4, :]
        G = nearest_spline(struct_pts, aero, n_nearest=4, weighting="inverse_distance")
        assert abs(G[0, 3] - 1.0) < 1e-12


# =============================================================================
# BEAM SPLINE
# =============================================================================


class TestBeamSpline:
    """Beam spline (SPLINE2) along a beam axis."""

    def test_shape(self):
        """Beam spline between line of struct and offset aero points."""
        struct = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], dtype=float)
        aero = np.array([[0.5, 0.5, 0], [1.5, 0.5, 0], [2.5, 0.5, 0]], dtype=float)
        G = beam_spline(struct, aero, beam_axis="x")
        assert G.shape == (3, 4)

    def test_partition_of_unity(self):
        """atol=1e-12."""
        struct = np.array([[0, 0, 0], [2, 0, 0], [4, 0, 0], [6, 0, 0]], dtype=float)
        aero = np.array([[1, 1, 0], [3, 1, 0], [5, 1, 0]], dtype=float)
        G = beam_spline(struct, aero, beam_axis="x")
        npt.assert_allclose(G.sum(axis=1), 1.0, atol=1e-12)

    def test_linear_along_axis(self):
        """Linear field along beam axis reproduced. atol=1e-10."""
        struct = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], dtype=float)
        aero = np.array([[0.5, 0.5, 0], [1.5, 0.5, 0], [2.5, 0.5, 0]], dtype=float)
        G = beam_spline(struct, aero, beam_axis="x")
        u_struct = struct[:, 0]  # u = x
        u_aero = G @ u_struct
        npt.assert_allclose(u_aero, aero[:, 0], atol=1e-10)


# =============================================================================
# CLIP TO CONVEX HULL
# =============================================================================


class TestClipToConvexHull:
    """Clip spline weights for aero points outside structural convex hull."""

    def test_interior_preserved(self, struct_pts):
        """Interior aero points keep original weights."""
        # Aero points fully inside struct hull [0,3]x[0,2]
        aero_inside = _rect_grid(3, 2, spacing=0.5, offset=(0.5, 0.5, 0.0))
        G = ips_spline(struct_pts, aero_inside)
        G_clipped = clip_to_convex_hull(G, struct_pts, aero_inside)
        npt.assert_allclose(G_clipped, G, atol=1e-14)

    def test_exterior_zeroed(self, struct_pts):
        """Aero point far outside hull gets zeroed row."""
        aero_outside = np.array([[100.0, 100.0, 0.0]])
        G = ips_spline(struct_pts, aero_outside)
        G_clipped = clip_to_convex_hull(G, struct_pts, aero_outside)
        npt.assert_allclose(G_clipped[0], 0.0, atol=1e-14)

    def test_returns_copy(self, struct_pts, aero_pts):
        """Original G is not modified."""
        G = ips_spline(struct_pts, aero_pts)
        G_orig = G.copy()
        _ = clip_to_convex_hull(G, struct_pts, aero_pts)
        npt.assert_array_equal(G, G_orig)


# =============================================================================
# VALIDATE SPLINE
# =============================================================================


class TestValidateSpline:
    """validate_spline diagnostics."""

    def test_good_spline(self, struct_pts, aero_pts):
        """Well-conditioned spline passes validation."""
        G = ips_spline(struct_pts, aero_pts)
        result = validate_spline(G)
        assert result["partition_of_unity_ok"]
        assert result["partition_of_unity_error"] < 1e-10

    def test_bad_spline(self):
        """Random matrix fails partition-of-unity check."""
        G_bad = np.random.default_rng(42).standard_normal((5, 3))
        result = validate_spline(G_bad)
        assert not result["partition_of_unity_ok"]


# =============================================================================
# FPS (Finite Plate Spline)
# =============================================================================


class TestFpsSpline:
    """FPS: BFS bicubic Hermite finite-element spline."""

    def test_shape(self, struct_pts, aero_pts):
        G = fps_spline(struct_pts, aero_pts)
        assert G.shape == (20, 12)

    def test_partition_of_unity(self, struct_pts, aero_pts):
        """atol=1e-10."""
        G = fps_spline(struct_pts, aero_pts)
        npt.assert_allclose(G.sum(axis=1), 1.0, atol=1e-10)
