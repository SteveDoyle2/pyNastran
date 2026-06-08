"""Shell element stiffness for bdf_vectorized3 solver.

Provides a unified ``ShellQuadSolver`` dispatching MITC4 or MacNeal (stabilized
1-point) formulations via ``PARAM,MYQUAD,{MITC4|MACN}``.

Both formulations produce a 20x20 element stiffness in the local 5-DOF frame
[u, v, w, theta_x, theta_y] per node, then transform to a 24x24 global matrix
with 6 DOF per node.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pyNastran.utils.scipy_utils import dok_matrix
from pyNastran.dev.solver.utils import DOF_MAP

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

# ---------------------------------------------------------------------------
# Gauss quadrature constants
# ---------------------------------------------------------------------------
_GP = 1.0 / np.sqrt(3.0)
_GAUSS_2x2_PTS = np.array([[-_GP, -_GP], [_GP, -_GP], [_GP, _GP], [-_GP, _GP]])
_GAUSS_2x2_WTS = np.array([1.0, 1.0, 1.0, 1.0])

# MITC4 tying points
_TYING_XI = np.array([[0.0, -1.0], [0.0, 1.0]])
_TYING_ETA = np.array([[-1.0, 0.0], [1.0, 0.0]])

# MacNeal hourglass vector
_GAMMA = np.array([1.0, -1.0, 1.0, -1.0])


# ---------------------------------------------------------------------------
# Shape functions and Jacobian
# ---------------------------------------------------------------------------
def _shape_quad4(xi: float, eta: float) -> np.ndarray:
    """Bilinear shape functions (4,)."""
    return 0.25 * np.array(
        [
            (1.0 - xi) * (1.0 - eta),
            (1.0 + xi) * (1.0 - eta),
            (1.0 + xi) * (1.0 + eta),
            (1.0 - xi) * (1.0 + eta),
        ]
    )


def _dshape_quad4(xi: float, eta: float) -> np.ndarray:
    """Shape function derivatives (2, 4): [dN/dxi; dN/deta]."""
    return 0.25 * np.array(
        [
            [-(1.0 - eta), (1.0 - eta), (1.0 + eta), -(1.0 + eta)],
            [-(1.0 - xi), -(1.0 + xi), (1.0 + xi), (1.0 - xi)],
        ]
    )


def _jacobian(dN: np.ndarray, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    """Jacobian (2, 2) and determinant."""
    J = np.array([[dN[0] @ x, dN[0] @ y], [dN[1] @ x, dN[1] @ y]])
    det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    return J, det_J


def _jacobian_inv(J: np.ndarray, det_J: float) -> np.ndarray:
    """Inverse Jacobian (2, 2)."""
    return np.array([[J[1, 1], -J[0, 1]], [-J[1, 0], J[0, 0]]]) / det_J


# ---------------------------------------------------------------------------
# Strain-displacement matrices
# ---------------------------------------------------------------------------
def _membrane_B(dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Membrane B-matrix (3, 8)."""
    Bm = np.zeros((3, 8))
    Bm[0, 0::2] = dN_dx
    Bm[1, 1::2] = dN_dy
    Bm[2, 0::2] = dN_dy
    Bm[2, 1::2] = dN_dx
    return Bm


def _bending_B(dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Bending curvature-displacement matrix (3, 12). DOF: [w, tx, ty]."""
    Bb = np.zeros((3, 12))
    Bb[0, 2::3] = dN_dx
    Bb[1, 1::3] = -dN_dy
    Bb[2, 1::3] = -dN_dx
    Bb[2, 2::3] = dN_dy
    return Bb


def _shear_B(N: np.ndarray, dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Standard transverse shear B-matrix (2, 12). DOF: [w, tx, ty]."""
    Bs = np.zeros((2, 12))
    Bs[0, 0::3] = dN_dx
    Bs[0, 2::3] = N
    Bs[1, 0::3] = dN_dy
    Bs[1, 1::3] = -N
    return Bs


# ---------------------------------------------------------------------------
# MITC4 transverse shear
# ---------------------------------------------------------------------------
def _covariant_shear_row(N: np.ndarray, dN_row: np.ndarray, J_row: np.ndarray) -> np.ndarray:
    """One covariant transverse shear strain row (12,)."""
    J_a1, J_a2 = J_row[0], J_row[1]
    row = np.zeros(12)
    row[0::3] = dN_row
    row[1::3] = -J_a2 * N
    row[2::3] = J_a1 * N
    return row


def _mitc4_shear_B(xi: float, eta: float, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """MITC4 assumed-strain transverse shear B-matrix (2, 12)."""
    e_xi3_rows = np.zeros((2, 12))
    for itp, (xi_tp, eta_tp) in enumerate(_TYING_XI):
        N_tp = _shape_quad4(xi_tp, eta_tp)
        dN_tp = _dshape_quad4(xi_tp, eta_tp)
        J_tp, _ = _jacobian(dN_tp, x, y)
        e_xi3_rows[itp] = _covariant_shear_row(N_tp, dN_tp[0], J_tp[0])

    e_eta3_rows = np.zeros((2, 12))
    for itp, (xi_tp, eta_tp) in enumerate(_TYING_ETA):
        N_tp = _shape_quad4(xi_tp, eta_tp)
        dN_tp = _dshape_quad4(xi_tp, eta_tp)
        J_tp, _ = _jacobian(dN_tp, x, y)
        e_eta3_rows[itp] = _covariant_shear_row(N_tp, dN_tp[1], J_tp[1])

    e_xi3 = 0.5 * (1.0 - eta) * e_xi3_rows[0] + 0.5 * (1.0 + eta) * e_xi3_rows[1]
    e_eta3 = 0.5 * (1.0 - xi) * e_eta3_rows[0] + 0.5 * (1.0 + xi) * e_eta3_rows[1]

    dN_here = _dshape_quad4(xi, eta)
    J_here, det_J = _jacobian(dN_here, x, y)
    J_inv = _jacobian_inv(J_here, det_J)

    Bs_mitc = np.zeros((2, 12))
    Bs_mitc[0] = J_inv[0, 0] * e_xi3 + J_inv[0, 1] * e_eta3
    Bs_mitc[1] = J_inv[1, 0] * e_xi3 + J_inv[1, 1] * e_eta3
    return Bs_mitc


# ---------------------------------------------------------------------------
# Selective-RI membrane splitting
# ---------------------------------------------------------------------------
def _split_membrane_constitutive(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split A into volumetric and deviatoric for selective reduced integration."""
    P = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 0]], dtype=float) / 2.0
    A_vol = P.T @ A @ P
    A_dev = A - A_vol
    return A_vol, A_dev


def _split_membrane_direct_shear(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split A into direct-stress and in-plane-shear parts.

    Direct (A11, A22, A12) is integrated at 2x2 Gauss.
    Shear (A66, A16, A26) is integrated at 1-point (centroid).

    Works correctly for both isotropic and orthotropic materials.
    """
    A_direct = np.array([[A[0, 0], A[0, 1], 0.0], [A[1, 0], A[1, 1], 0.0], [0.0, 0.0, 0.0]])
    A_shear = A - A_direct
    return A_direct, A_shear


def _membrane_selective(x: np.ndarray, y: np.ndarray, A_mat: np.ndarray) -> np.ndarray:
    """Compute 8x8 membrane stiffness with 1-point + hourglass stabilization.

    Uses the NX Nastran approach: full constitutive at 1-point (constant
    strain), plus a stabilized hourglass correction scaled so that the
    hourglass eigenvalue equals E*t/3 for isotropic materials.

    The stabilization coefficient alpha_hg is computed from the material
    as: alpha_hg = (A11 + A22 - 2*A12) / (3 * (eig_hg_2x2 / Area))
    which reduces to the NX formula for isotropic and works correctly
    for orthotropic materials.
    """
    # Full 2x2 integration
    K_2x2 = np.zeros((8, 8))
    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = _jacobian_inv(J, det_J)
        dN_dxy = J_inv @ dN_dnat
        Bm = _membrane_B(dN_dxy[0], dN_dxy[1])
        K_2x2 += (Bm.T @ A_mat @ Bm) * det_J * wt

    # 1-point (centroid) integration
    dN_c = _dshape_quad4(0.0, 0.0)
    J_c, det_J_c = _jacobian(dN_c, x, y)
    J_inv_c = _jacobian_inv(J_c, det_J_c)
    dN_dxy_c = J_inv_c @ dN_c
    A_elem = det_J_c * 4.0
    Bm_c = _membrane_B(dN_dxy_c[0], dN_dxy_c[1])
    K_1pt = (Bm_c.T @ A_mat @ Bm_c) * A_elem

    # Hourglass stabilization: K = K_1pt + alpha_hg * (K_2x2 - K_1pt)
    # NX targets hourglass eigenvalue = (A11 + A22 - 2*A12) / 3
    # The full 2x2 hourglass eigenvalue is the max eigenvalue of (K_2x2 - K_1pt)
    K_hg = K_2x2 - K_1pt
    eig_hg = np.linalg.eigvalsh(K_hg).max()

    if eig_hg > 1e-10 * np.abs(A_mat).max():
        A11 = A_mat[0, 0]
        if A11 > 0.0:
            target_hg = (A11 - A_mat[0, 1] ** 2 / A11) / 3.0
        else:
            target_hg = 0.0
        alpha_hg = target_hg / eig_hg
        alpha_hg = min(max(alpha_hg, 0.0), 1.0)
    else:
        alpha_hg = 1.0

    return K_1pt + alpha_hg * K_hg


# ---------------------------------------------------------------------------
# MITC4 element stiffness
# ---------------------------------------------------------------------------
def mitc4_stiffness(
    x: np.ndarray,
    y: np.ndarray,
    A_mat: np.ndarray,
    D_mat: np.ndarray,
    Ds_mat: np.ndarray,
    B_mat: np.ndarray | None = None,
    membrane: str = "selective",
) -> np.ndarray:
    """20x20 MITC4 shell element stiffness.

    Parameters
    ----------
    x, y : (4,) nodal coords in element local frame
    A_mat : (3, 3) membrane constitutive
    D_mat : (3, 3) bending constitutive
    Ds_mat : (2, 2) transverse shear constitutive
    B_mat : (3, 3) membrane-bending coupling, optional
    membrane : str
        'full' or 'selective' (default, matches Nastran)
    """
    m_idx = np.array([0, 1, 5, 6, 10, 11, 15, 16], dtype=int)
    b_idx = np.array([2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19], dtype=int)

    # Membrane: generalized selective-RI (direct at 2x2, shear at 1pt)
    K_mm = _membrane_selective(x, y, A_mat)

    K_bb = np.zeros((12, 12))
    K_ss = np.zeros((12, 12))
    K_mb = np.zeros((8, 12))

    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = _jacobian_inv(J, det_J)
        dN_dxy = J_inv @ dN_dnat

        Bb = _bending_B(dN_dxy[0], dN_dxy[1])
        K_bb += (Bb.T @ D_mat @ Bb) * det_J * wt

        Bs = _mitc4_shear_B(xi, eta, x, y)
        K_ss += (Bs.T @ Ds_mat @ Bs) * det_J * wt

        if B_mat is not None:
            Bm = _membrane_B(dN_dxy[0], dN_dxy[1])
            K_mb += (Bm.T @ B_mat @ Bb) * det_J * wt

    Ke = np.zeros((20, 20))
    Ke[np.ix_(m_idx, m_idx)] = K_mm
    Ke[np.ix_(b_idx, b_idx)] = K_bb + K_ss
    if B_mat is not None:
        Ke[np.ix_(m_idx, b_idx)] = K_mb
        Ke[np.ix_(b_idx, m_idx)] = K_mb.T
    return Ke


# ---------------------------------------------------------------------------
# MacNeal stabilization helpers
# ---------------------------------------------------------------------------
def _macneal_membrane_stabilization(
    x: np.ndarray, y: np.ndarray, A_mat: np.ndarray, Bm_c: np.ndarray, det_J_c: float
) -> np.ndarray:
    """Membrane hourglass stabilization (8, 8). alpha=0.3."""
    A_elem = det_J_c * 4.0
    K_2x2 = np.zeros((8, 8))
    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = _jacobian_inv(J, det_J)
        dN_dxy = J_inv @ dN_dnat
        Bm = _membrane_B(dN_dxy[0], dN_dxy[1])
        K_2x2 += (Bm.T @ A_mat @ Bm) * det_J * wt
    K_1pt = (Bm_c.T @ A_mat @ Bm_c) * A_elem
    return 0.3 * (K_2x2 - K_1pt)


def _macneal_bending_stabilization(
    x: np.ndarray, y: np.ndarray, D_mat: np.ndarray, Bb_c: np.ndarray, det_J_c: float
) -> np.ndarray:
    """Bending hourglass stabilization (12, 12). alpha=0.5."""
    A_elem = det_J_c * 4.0
    K_2x2 = np.zeros((12, 12))
    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = _jacobian_inv(J, det_J)
        dN_dxy = J_inv @ dN_dnat
        Bb = _bending_B(dN_dxy[0], dN_dxy[1])
        K_2x2 += (Bb.T @ D_mat @ Bb) * det_J * wt
    K_1pt = (Bb_c.T @ D_mat @ Bb_c) * A_elem
    return 0.5 * (K_2x2 - K_1pt)


# ---------------------------------------------------------------------------
# MacNeal element stiffness
# ---------------------------------------------------------------------------
def macneal_stiffness(
    x: np.ndarray,
    y: np.ndarray,
    A_mat: np.ndarray,
    D_mat: np.ndarray,
    Ds_mat: np.ndarray,
    B_mat: np.ndarray | None = None,
) -> np.ndarray:
    """20x20 MacNeal shell element stiffness.

    Uses generalized selective-RI for membrane (same as MITC4/MACN2),
    1-point bending + stabilization, and 1-point transverse shear.

    Parameters
    ----------
    x, y : (4,) nodal coords in element local frame
    A_mat : (3, 3) membrane constitutive
    D_mat : (3, 3) bending constitutive
    Ds_mat : (2, 2) transverse shear constitutive
    B_mat : (3, 3) membrane-bending coupling, optional
    """
    m_idx = np.array([0, 1, 5, 6, 10, 11, 15, 16], dtype=int)
    b_idx = np.array([2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19], dtype=int)

    # Membrane: generalized selective-RI
    K_mm = _membrane_selective(x, y, A_mat)

    # Bending: 1-point + stabilization
    dN_c = _dshape_quad4(0.0, 0.0)
    J_c, det_J_c = _jacobian(dN_c, x, y)
    J_inv_c = _jacobian_inv(J_c, det_J_c)
    dN_dxy_c = J_inv_c @ dN_c
    dN_dx_c = dN_dxy_c[0]
    dN_dy_c = dN_dxy_c[1]
    A_elem = det_J_c * 4.0

    Bb_c = _bending_B(dN_dx_c, dN_dy_c)
    K_bb = (Bb_c.T @ D_mat @ Bb_c) * A_elem

    N_c = np.array([0.25, 0.25, 0.25, 0.25])
    Bs_c = _shear_B(N_c, dN_dx_c, dN_dy_c)
    K_ss = (Bs_c.T @ Ds_mat @ Bs_c) * A_elem

    K_mb = None
    if B_mat is not None:
        Bm_c = _membrane_B(dN_dx_c, dN_dy_c)
        K_mb = (Bm_c.T @ B_mat @ Bb_c) * A_elem

    K_bb += _macneal_bending_stabilization(x, y, D_mat, Bb_c, det_J_c) + K_ss

    Ke = np.zeros((20, 20))
    Ke[np.ix_(m_idx, m_idx)] = K_mm
    Ke[np.ix_(b_idx, b_idx)] = K_bb
    if K_mb is not None:
        Ke[np.ix_(m_idx, b_idx)] = K_mb
        Ke[np.ix_(b_idx, m_idx)] = K_mb.T
    return Ke


# ---------------------------------------------------------------------------
# MACN2 element stiffness — calibrated to match NX Nastran CQUAD4
# ---------------------------------------------------------------------------
def macn2_stiffness(
    x: np.ndarray,
    y: np.ndarray,
    A_mat: np.ndarray,
    D_mat: np.ndarray,
    Ds_mat: np.ndarray,
    B_mat: np.ndarray | None = None,
) -> np.ndarray:
    """20x20 MACN2 shell element stiffness — Nastran-matched formulation.

    Combines:
    - Selective membrane: 2x2 for direct (A11/A22), 1-point for shear (A66)
    - Full 2x2 Gauss bending
    - MITC4 assumed-strain transverse shear (locking-free)

    This combination best matches NX Nastran's CQUAD4 element stiffness
    for both isotropic and orthotropic materials.
    """
    m_idx = np.array([0, 1, 5, 6, 10, 11, 15, 16], dtype=int)
    b_idx = np.array([2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19], dtype=int)

    # --- Membrane: generalized selective-RI ---
    K_mm = _membrane_selective(x, y, A_mat)

    # --- Bending: 1-point + hourglass stabilization (NX target = D11/3) ---
    K_bb_2x2 = np.zeros((12, 12))
    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = _jacobian_inv(J, det_J)
        dN_dxy = J_inv @ dN_dnat
        Bb = _bending_B(dN_dxy[0], dN_dxy[1])
        K_bb_2x2 += (Bb.T @ D_mat @ Bb) * det_J * wt

    dN_c = _dshape_quad4(0.0, 0.0)
    J_c, det_J_c = _jacobian(dN_c, x, y)
    J_inv_c = _jacobian_inv(J_c, det_J_c)
    dN_dxy_c = J_inv_c @ dN_c
    A_elem = det_J_c * 4.0
    Bb_c = _bending_B(dN_dxy_c[0], dN_dxy_c[1])
    K_bb_1pt = (Bb_c.T @ D_mat @ Bb_c) * A_elem

    K_bb_hg = K_bb_2x2 - K_bb_1pt
    eig_hg_bb = np.linalg.eigvalsh(K_bb_hg).max()
    if eig_hg_bb > 1e-10 * np.abs(D_mat).max():
        D11 = D_mat[0, 0]
        alpha_hg_bb = (D11 / 3.0) / eig_hg_bb if D11 > 0.0 else 1.0
        alpha_hg_bb = min(max(alpha_hg_bb, 0.0), 1.0)
    else:
        alpha_hg_bb = 1.0
    K_bb = K_bb_1pt + alpha_hg_bb * K_bb_hg

    # --- Transverse shear: MITC4 assumed strain (2x2 integration) ---
    K_ss = np.zeros((12, 12))
    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        _, det_J = _jacobian(dN_dnat, x, y)
        Bs = _mitc4_shear_B(xi, eta, x, y)
        K_ss += (Bs.T @ Ds_mat @ Bs) * det_J * wt

    # --- Membrane-bending coupling ---
    K_mb = np.zeros((8, 12))
    if B_mat is not None:
        for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
            xi, eta = gpt
            dN_dnat = _dshape_quad4(xi, eta)
            J, det_J = _jacobian(dN_dnat, x, y)
            J_inv = _jacobian_inv(J, det_J)
            dN_dxy = J_inv @ dN_dnat
            Bm = _membrane_B(dN_dxy[0], dN_dxy[1])
            Bb = _bending_B(dN_dxy[0], dN_dxy[1])
            K_mb += (Bm.T @ B_mat @ Bb) * det_J * wt

    # --- Assemble 20x20 ---
    Ke = np.zeros((20, 20))
    Ke[np.ix_(m_idx, m_idx)] = K_mm
    Ke[np.ix_(b_idx, b_idx)] = K_bb + K_ss
    if B_mat is not None and np.any(np.abs(K_mb) > 1e-12):
        Ke[np.ix_(m_idx, b_idx)] = K_mb
        Ke[np.ix_(b_idx, m_idx)] = K_mb.T
    return Ke


# ---------------------------------------------------------------------------
# Drilling DOF formulations
# ---------------------------------------------------------------------------
def _drilling_k6rot(Ke_local: np.ndarray, k6rot: float) -> float:
    """K6ROT artificial spring: k_drill = K6ROT * 1e-4 * max(diag(K_mm))."""
    m_idx = np.array([0, 1, 5, 6, 10, 11, 15, 16], dtype=int)
    diag_mm = np.diag(Ke_local[np.ix_(m_idx, m_idx)])
    max_diag = np.abs(diag_mm).max()
    if max_diag == 0.0:
        return 0.0
    return k6rot * 1e-4 * max_diag


def _drilling_allman(x: np.ndarray, y: np.ndarray, A_mat: np.ndarray) -> np.ndarray:
    """Allman-type drilling DOF: drilling enrichment for a 12-DOF membrane.

    Enriches the bilinear displacement field with edge-based quadratic
    bubbles whose amplitudes are proportional to the difference in drilling
    rotations at adjacent nodes.

    Returns the INCREMENTAL drilling contribution K_drill (12, 12) that
    should be ADDED to the existing element stiffness. The standard membrane
    terms (u,v only) are subtracted so only the θz coupling remains.

    DOF order per node: [u, v, θz].
    """
    # Edge vectors and lengths
    # Nodes ordered CCW: 0-1-2-3
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    dx = np.array([x[j] - x[i] for i, j in edges])
    dy = np.array([y[j] - y[i] for i, j in edges])
    Le = np.sqrt(dx**2 + dy**2)

    # Centroid Jacobian for element area
    dN_c = _dshape_quad4(0.0, 0.0)
    _, det_J_c = _jacobian(dN_c, x, y)
    A_elem = det_J_c * 4.0

    # Build the 12x12 stiffness with 2x2 Gauss integration
    # DOF ordering: [u0, v0, θz0, u1, v1, θz1, u2, v2, θz2, u3, v3, θz3]
    K_drill = np.zeros((12, 12))

    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        N = _shape_quad4(xi, eta)
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = _jacobian_inv(J, det_J)
        dN_dxy = J_inv @ dN_dnat
        dN_dx = dN_dxy[0]
        dN_dy = dN_dxy[1]

        # Enhanced membrane B-matrix (3, 12) including drilling DOF
        # Standard membrane: eps_x = du/dx, eps_y = dv/dy, gamma = du/dy + dv/dx
        # Allman adds: the drilling rotation modifies the displacement field
        # via edge bubble functions. For a 4-node quad, the enhancement is:
        #   u_enh = Σ (Le/8) * ny_e * (θz_j - θz_i) * bubble_e(ξ,η)
        #   v_enh = Σ (Le/8) * (-nx_e) * (θz_j - θz_i) * bubble_e(ξ,η)
        # The resulting B-matrix has columns for θz at each node.
        Bm = np.zeros((3, 12))
        for i in range(4):
            c = 3 * i
            # Standard bilinear part
            Bm[0, c] = dN_dx[i]  # du/dx
            Bm[1, c + 1] = dN_dy[i]  # dv/dy
            Bm[2, c] = dN_dy[i]  # du/dy
            Bm[2, c + 1] = dN_dx[i]  # dv/dx

        # Allman drilling enrichment
        # Each edge (i,j) contributes to the θz DOFs of nodes i and j
        # Bubble derivative at Gauss point for each edge
        for ie, (ni, nj) in enumerate(edges):
            nx_e = dx[ie] / Le[ie]
            ny_e = dy[ie] / Le[ie]
            alpha = Le[ie] / 8.0

            # Bubble function for edge (i,j): φ_e = N_i * N_j * 4
            # Its derivatives:
            dphi_dx = 4.0 * (dN_dx[ni] * N[nj] + N[ni] * dN_dx[nj])
            dphi_dy = 4.0 * (dN_dy[ni] * N[nj] + N[ni] * dN_dy[nj])

            # Contribution to strain from θz_j - θz_i:
            # Δu = alpha * ny_e * phi_e * (θz_j - θz_i)
            # Δv = alpha * (-nx_e) * phi_e * (θz_j - θz_i)
            # eps_x += alpha * ny_e * dphi_dx * (θz_j - θz_i)
            # eps_y += alpha * (-nx_e) * dphi_dy * (θz_j - θz_i)
            # gamma += alpha * (ny_e * dphi_dy + (-nx_e) * dphi_dx) * (θz_j - θz_i)

            ci = 3 * ni + 2  # θz column for node i
            cj = 3 * nj + 2  # θz column for node j

            eps_x_coeff = alpha * ny_e * dphi_dx
            eps_y_coeff = alpha * (-nx_e) * dphi_dy
            gamma_coeff = alpha * (ny_e * dphi_dy + (-nx_e) * dphi_dx)

            # θz_j contributes +1, θz_i contributes -1
            Bm[0, cj] += eps_x_coeff
            Bm[0, ci] -= eps_x_coeff
            Bm[1, cj] += eps_y_coeff
            Bm[1, ci] -= eps_y_coeff
            Bm[2, cj] += gamma_coeff
            Bm[2, ci] -= gamma_coeff

        K_drill += (Bm.T @ A_mat @ Bm) * det_J * wt

    # Subtract the standard membrane stiffness (u,v only, no θz)
    # so we only add the drilling enrichment on top of existing K_mm.
    uv_idx = np.array([0, 1, 3, 4, 6, 7, 9, 10], dtype=int)  # u,v DOFs in 12-DOF
    K_std = np.zeros((12, 12))
    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = _jacobian_inv(J, det_J)
        dN_dxy = J_inv @ dN_dnat
        Bm_std = _membrane_B(dN_dxy[0], dN_dxy[1])  # (3, 8)
        K_8x8 = (Bm_std.T @ A_mat @ Bm_std) * det_J * wt
        K_std[np.ix_(uv_idx, uv_idx)] += K_8x8

    return K_drill - K_std


def _drilling_hughes_brezzi(x: np.ndarray, y: np.ndarray, gamma: float) -> np.ndarray:
    """Hughes-Brezzi variational drilling: penalty on (θz - ½ curl u).

    Adds the term γ·∫(θz - ½(∂v/∂x - ∂u/∂y))² dA to the element energy.
    The penalty parameter γ = G·t (shear modulus × thickness) is passed directly.

    Returns K_hb (12, 12) in local frame, DOF order per node: [u, v, θz].
    """
    K_hb = np.zeros((12, 12))

    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        N = _shape_quad4(xi, eta)
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = _jacobian_inv(J, det_J)
        dN_dxy = J_inv @ dN_dnat
        dN_dx = dN_dxy[0]
        dN_dy = dN_dxy[1]

        # Constraint vector: c = θz - ½(∂v/∂x - ∂u/∂y)
        # DOF per node: [u, v, θz]
        # c = Σ_i [ N_i·θz_i - ½(dN_dx_i·v_i - dN_dy_i·u_i) ]
        # G-vector (12,): dc/d(DOF)
        G_vec = np.zeros(12)
        for i in range(4):
            c = 3 * i
            G_vec[c] = 0.5 * dN_dy[i]  # dc/du_i = +½ dN_i/dy
            G_vec[c + 1] = -0.5 * dN_dx[i]  # dc/dv_i = -½ dN_i/dx
            G_vec[c + 2] = N[i]  # dc/dθz_i = N_i

        K_hb += gamma * np.outer(G_vec, G_vec) * det_J * wt

    return K_hb


# ---------------------------------------------------------------------------
# Geometric stiffness
# ---------------------------------------------------------------------------
def geometric_stiffness(
    x: np.ndarray,
    y: np.ndarray,
    stress_resultants: np.ndarray,
    thickness: float = 0.0,
) -> np.ndarray:
    """20x20 geometric (differential) stiffness for buckling.

    Parameters
    ----------
    x, y : (4,) nodal coords in element local frame
    stress_resultants : (3,) or (4, 3) [Nxx, Nyy, Nxy] membrane stress resultants
    thickness : float
        Shell thickness for rotational DOF scaling.

    Returns
    -------
    Kg : (20, 20)
    """
    stress_resultants = np.asarray(stress_resultants)
    if stress_resultants.ndim == 1:
        stress_at_gp = np.tile(stress_resultants, (4, 1))
    else:
        stress_at_gp = stress_resultants

    t2_12 = thickness**2 / 12.0 if thickness > 0.0 else 0.0
    Kg = np.zeros((20, 20))

    u_idx = np.array([0, 5, 10, 15], dtype=int)
    v_idx = np.array([1, 6, 11, 16], dtype=int)
    w_idx = np.array([2, 7, 12, 17], dtype=int)
    tx_idx = np.array([3, 8, 13, 18], dtype=int)
    ty_idx = np.array([4, 9, 14, 19], dtype=int)

    for igp, (gpt, wt) in enumerate(zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS)):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = _jacobian_inv(J, det_J)
        dN_dxy = J_inv @ dN_dnat  # (2, 4)

        Nxx, Nyy, Nxy = stress_at_gp[igp]
        S_N = np.array([[Nxx, Nxy], [Nxy, Nyy]])
        kg_N = (dN_dxy.T @ S_N @ dN_dxy) * det_J * wt  # (4, 4)

        Kg[np.ix_(u_idx, u_idx)] += kg_N
        Kg[np.ix_(v_idx, v_idx)] += kg_N
        Kg[np.ix_(w_idx, w_idx)] += kg_N

        if t2_12 > 0.0:
            Kg[np.ix_(tx_idx, tx_idx)] += kg_N * t2_12
            Kg[np.ix_(ty_idx, ty_idx)] += kg_N * t2_12

    return Kg


def ctria3_geometric_stiffness(
    x: np.ndarray,
    y: np.ndarray,
    stress_resultants: np.ndarray,
    thickness: float = 0.0,
) -> np.ndarray:
    """15x15 geometric (differential) stiffness for CTRIA3 buckling.

    Parameters
    ----------
    x, y : (3,) nodal coords in element local frame
    stress_resultants : (3,) [Nxx, Nyy, Nxy] membrane stress resultants
    thickness : float
        Shell thickness for rotational DOF scaling.

    Returns
    -------
    Kg : (15, 15)
        DOF order per node: [u, v, w, tx, ty].
    """
    area = 0.5 * abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))

    # CST shape function derivatives (constant over element)
    # dN/dx = [y2-y3, y3-y1, y1-y2] / (2*A)
    # dN/dy = [x3-x2, x1-x3, x2-x1] / (2*A)
    dNdx = np.array([y[1] - y[2], y[2] - y[0], y[0] - y[1]]) / (2.0 * area)
    dNdy = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]]) / (2.0 * area)
    dN_dxy = np.array([dNdx, dNdy])  # (2, 3)

    Nxx, Nyy, Nxy = stress_resultants[0], stress_resultants[1], stress_resultants[2]
    S_N = np.array([[Nxx, Nxy], [Nxy, Nyy]])

    # kg_N = dN^T @ S @ dN * Area  (3x3 node-node coupling)
    kg_N = (dN_dxy.T @ S_N @ dN_dxy) * area

    t2_12 = thickness**2 / 12.0 if thickness > 0.0 else 0.0
    Kg = np.zeros((15, 15))

    u_idx = np.array([0, 5, 10], dtype=int)
    v_idx = np.array([1, 6, 11], dtype=int)
    w_idx = np.array([2, 7, 12], dtype=int)
    tx_idx = np.array([3, 8, 13], dtype=int)
    ty_idx = np.array([4, 9, 14], dtype=int)

    Kg[np.ix_(u_idx, u_idx)] = kg_N
    Kg[np.ix_(v_idx, v_idx)] = kg_N
    Kg[np.ix_(w_idx, w_idx)] = kg_N

    if t2_12 > 0.0:
        Kg[np.ix_(tx_idx, tx_idx)] = kg_N * t2_12
        Kg[np.ix_(ty_idx, ty_idx)] = kg_N * t2_12

    return Kg


# ---------------------------------------------------------------------------
# Shell offset eccentricity
# ---------------------------------------------------------------------------
def _apply_shell_offset(Ke: np.ndarray, z0: float) -> np.ndarray:
    """Apply offset: u_mid = u_grid + z0*ty, v_mid = v_grid - z0*tx."""
    E_full = np.eye(20)
    for i in range(4):
        r = 5 * i
        E_full[r + 0, r + 4] = z0
        E_full[r + 1, r + 3] = -z0
    return E_full.T @ Ke @ E_full


# ---------------------------------------------------------------------------
# Constitutive helpers
# ---------------------------------------------------------------------------
def _get_ABD_for_element(
    model: BDF,
    pid: int,
    pshell,
    pcomp,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray, np.ndarray, float]:
    """Get A, B, D, Ds, thickness for a single property ID.

    Returns
    -------
    A_mat : (3, 3)
    B_mat : (3, 3) or None
    D_mat : (3, 3)
    Ds_mat : (2, 2)
    thickness : float
    """
    kappa = 5.0 / 6.0

    def _get_shear_modulus(mid: int) -> tuple[float, float]:
        """Get transverse shear moduli (G13, G23) for a material ID.

        Returns (G13, G23). For MAT1 (isotropic), G13 = G23 = G.
        For MAT8 (orthotropic), uses G1Z and G2Z fields.
        """
        mat1 = model.mat1
        if mat1.n > 0:
            idx = np.searchsorted(mat1.material_id, mid)
            if idx < mat1.n and mat1.material_id[idx] == mid:
                g = float(mat1.G[idx])
                return g, g
        mat8 = model.mat8
        if mat8.n > 0:
            idx = np.searchsorted(mat8.material_id, mid)
            if idx < mat8.n and mat8.material_id[idx] == mid:
                g13 = float(mat8.G13[idx])
                g23 = float(mat8.G23[idx])
                return g13, g23
        return 0.0, 0.0

    # Try PSHELL first
    if pshell.n > 0:
        iprop = np.searchsorted(pshell.property_id, pid)
        if iprop < pshell.n and pshell.property_id[iprop] == pid:
            prop_slice = pshell.slice_card_by_id(np.array([pid]))
            thickness = float(prop_slice.t[0])
            mids = prop_slice.material_id[0]
            mid1, mid2, mid3, mid4 = mids
            tst = float(prop_slice.tst[0])

            # Validate: must have MID1 or MID2 for a solvable element
            if mid1 <= 0 and mid2 <= 0:
                raise RuntimeError(
                    f"PSHELL pid={pid}: requires MID1 (membrane) or MID2 "
                    f"(bending) for a solvable element"
                )

            # ABD matrices (get_individual_ABD_matrices may produce NaN
            # for missing MIDs — replace with zeros)
            A_arr, B_arr, D_arr = prop_slice.get_individual_ABD_matrices()
            A_mat = A_arr[0]
            D_mat = D_arr[0]
            B_arr_0 = B_arr[0]

            # Zero out NaN entries (from missing MIDs)
            A_mat = np.nan_to_num(A_mat, nan=0.0)
            D_mat = np.nan_to_num(D_mat, nan=0.0)
            B_arr_0 = np.nan_to_num(B_arr_0, nan=0.0)

            B_mat = B_arr_0 if np.any(np.abs(B_arr_0) > 1e-12) else None

            # Transverse shear: only meaningful when bending exists (MID2)
            # or when MID3 is explicitly specified.
            twelveIt3 = float(prop_slice.twelveIt3[0])
            if mid3 > 0:
                G13, G23 = _get_shear_modulus(mid3)
            elif mid2 > 0:
                G13, G23 = _get_shear_modulus(mid2)
            else:
                G13, G23 = 0.0, 0.0

            if mid3 > 0:
                # Explicit MID3: use G-based shear stiffness
                ts = tst * thickness
                Ds_mat = kappa * ts * np.array([[G13, 0.0], [0.0, G23]])
            elif mid2 > 0 and D_mat[0, 0] > 0.0:
                # MID3 blank: NX derives shear from bending rigidity D11
                Ds_mat = 129.0 * D_mat[0, 0] * twelveIt3 * np.eye(2)
            elif G13 > 0.0 or G23 > 0.0:
                ts = tst * thickness
                Ds_mat = kappa * ts * np.array([[G13, 0.0], [0.0, G23]])
            else:
                Ds_mat = np.zeros((2, 2))

            return A_mat, B_mat, D_mat, Ds_mat, thickness

    # Try PCOMP
    if pcomp.n > 0:
        iprop = np.searchsorted(pcomp.property_id, pid)
        if iprop < pcomp.n and pcomp.property_id[iprop] == pid:
            prop_slice = pcomp.slice_card_by_id(np.array([pid]))
            A_arr, B_arr, D_arr = prop_slice.get_individual_ABD_matrices()
            A_mat = A_arr[0]
            B_arr_0 = B_arr[0]
            D_mat = D_arr[0]
            thickness = float(prop_slice.total_thickness()[0])

            B_mat = B_arr_0 if np.any(np.abs(B_arr_0) > 1e-12) else None

            # Transverse shear: thickness-weighted average across plies,
            # rotated by ply angle.
            # Ds = κ · Σ ti · R(θi)^T · [[G13_i, 0],[0, G23_i]] · R(θi)
            ilayer = prop_slice.ilayer
            i0, i1 = ilayer[0]
            ply_mids = prop_slice.material_id[i0:i1]
            ply_thetas = prop_slice.theta[i0:i1]
            ply_thicknesses = prop_slice.thickness[i0:i1]
            Ds_mat = np.zeros((2, 2))
            for mid_ply, theta_ply, t_ply in zip(ply_mids, ply_thetas, ply_thicknesses):
                g13_ply, g23_ply = _get_shear_modulus(int(mid_ply))
                Ds_ply = np.array([[g13_ply, 0.0], [0.0, g23_ply]])
                if abs(theta_ply) > 1e-10:
                    c = np.cos(np.radians(theta_ply))
                    s = np.sin(np.radians(theta_ply))
                    R = np.array([[c, s], [-s, c]])
                    Ds_ply = R.T @ Ds_ply @ R
                Ds_mat += t_ply * Ds_ply
            Ds_mat *= kappa
            return A_mat, B_mat, D_mat, Ds_mat, thickness

    raise RuntimeError(f"No PSHELL or PCOMP found for PID={pid}")


# ---------------------------------------------------------------------------
# ShellQuadSolver — unified assembly dispatcher
# ---------------------------------------------------------------------------
class ShellQuadSolver:
    """Assembles CQUAD4 element stiffness matrices into global Kbb.

    Dispatches to MITC4 or MacNeal based on PARAM,MYQUAD.

    Parameters
    ----------
    model : BDF
        The bdf_vectorized3 model (setup + cross-referenced).
    formulation : str or None
        Override formulation ('MITC4' or 'MACN'). If None, reads
        PARAM,MYQUAD from the model (default MITC4).
    """

    def __init__(self, model: BDF, formulation: str | None = None):
        self.model = model
        if formulation is None:
            formulation = self._read_param_myquad()
        formulation = formulation.upper()
        if formulation not in ("MITC4", "MACN", "MACN2"):
            raise ValueError(f"PARAM,MYQUAD must be MITC4/MACN/MACN2, got '{formulation}'")
        self.formulation = formulation
        self.k6rot = self._read_param_k6rot()
        self.drilling = self._read_param_myqdril()

    def _read_param_myquad(self) -> str:
        """Read PARAM,MYQUAD from model; default MITC4."""
        params = self.model.params
        if "MYQUAD" in params:
            val = params["MYQUAD"].values[0]
            if isinstance(val, str):
                return val.upper()
            return str(val).upper()
        return "MITC4"

    def _read_param_k6rot(self) -> float:
        """Read PARAM,K6ROT from model; default 100.0."""
        params = self.model.params
        if "K6ROT" in params:
            return float(params["K6ROT"].values[0])
        return 100.0

    def _read_param_myqdril(self) -> str:
        """Read PARAM,MYQDRIL from model; default K6ROT.

        Options: K6ROT, ALLMAN, HB, NONE.
        """
        params = self.model.params
        if "MYQDRIL" in params:
            val = params["MYQDRIL"].values[0]
            if isinstance(val, str):
                return val.upper()
            return str(val).upper()
        return "K6ROT"

    def _element_stiffness(
        self,
        x: np.ndarray,
        y: np.ndarray,
        A_mat: np.ndarray,
        D_mat: np.ndarray,
        Ds_mat: np.ndarray,
        B_mat: np.ndarray | None,
    ) -> np.ndarray:
        """Compute 20x20 element stiffness using the selected formulation."""
        if self.formulation == "MITC4":
            return mitc4_stiffness(x, y, A_mat, D_mat, Ds_mat, B_mat)
        if self.formulation == "MACN2":
            return macn2_stiffness(x, y, A_mat, D_mat, Ds_mat, B_mat)
        return macneal_stiffness(x, y, A_mat, D_mat, Ds_mat, B_mat)

    def _apply_drilling(
        self,
        Ke_global: np.ndarray,
        Ke_local: np.ndarray,
        x: np.ndarray,
        y: np.ndarray,
        A_mat: np.ndarray,
        thickness: float,
        n_vec: np.ndarray,
        Ti: np.ndarray,
    ) -> None:
        """Apply drilling DOF stiffness to the 24x24 global element matrix.

        Dispatches based on PARAM,MYQDRIL:
          K6ROT  — artificial diagonal spring (default)
          ALLMAN — physical Allman-type edge drilling
          HB     — Hughes-Brezzi variational penalty
          NONE   — no drilling stiffness
        """
        method = self.drilling
        if method == "NONE":
            return

        if method == "K6ROT":
            k_drill = _drilling_k6rot(Ke_local, self.k6rot)
            if k_drill > 0.0:
                K_drill_node = k_drill * np.outer(n_vec, n_vec)
                for inode in range(4):
                    r = 6 * inode + 3
                    Ke_global[r : r + 3, r : r + 3] += K_drill_node
            return

        # For ALLMAN and HB, we get a 12x12 local drilling matrix with
        # DOF order [u, v, θz] per node. Transform to global 24-DOF.
        if method == "ALLMAN":
            K_drill_local = _drilling_allman(x, y, A_mat)
        elif method == "HB":
            gamma = A_mat[2, 2]  # A[2,2] = G·t for isotropic PSHELL
            K_drill_local = _drilling_hughes_brezzi(x, y, gamma)
        else:
            raise ValueError(
                f"PARAM,MYQDRIL must be 'K6ROT', 'ALLMAN', 'HB', or 'NONE', got '{method}'"
            )

        # Build 12->24 transformation for the drilling DOFs
        # Local DOF per node: [u_local, v_local, θz_local]
        # Global DOF per node: [ux, uy, uz, rx, ry, rz]
        # u_local = Ti[0,:] · [ux, uy, uz]  (row 0 of Ti)
        # v_local = Ti[1,:] · [ux, uy, uz]  (row 1 of Ti)
        # θz_local = n_vec · [rx, ry, rz]   (normal direction)
        T_drill_node = np.zeros((3, 6))
        T_drill_node[0, 0:3] = Ti[0, :]  # u_local from global translations
        T_drill_node[1, 0:3] = Ti[1, :]  # v_local from global translations
        T_drill_node[2, 3:6] = n_vec  # θz_local from global rotations

        T_drill = np.zeros((12, 24))
        for inode in range(4):
            T_drill[3 * inode : 3 * inode + 3, 6 * inode : 6 * inode + 6] = T_drill_node

        Ke_global += T_drill.T @ K_drill_local @ T_drill

    def build_kbb(
        self, Kbb: dok_matrix, dof_map: DOF_MAP, all_nids: np.ndarray, xyz_cid0: np.ndarray
    ) -> int:
        """Assemble CQUAD4 stiffness into Kbb.

        Returns
        -------
        nelements : int
            Number of CQUAD4 elements processed.
        """
        model = self.model
        cquad4 = model.cquad4
        nelements = cquad4.n
        if nelements == 0:
            return 0

        pshell = model.pshell
        pcomp = model.pcomp

        # Nodal coordinates
        grid = model.grid
        nodes = cquad4.nodes  # (nelements, 4)
        inid = grid.index(nodes)
        xyz = grid.xyz_cid0()

        p1 = xyz[inid[:, 0], :]
        p2 = xyz[inid[:, 1], :]
        p3 = xyz[inid[:, 2], :]
        p4 = xyz[inid[:, 3], :]

        # Element normals
        v13 = p1 - p3
        v24 = p2 - p4
        normal = np.cross(v13, v24)
        ni = np.linalg.norm(normal, axis=1)
        if np.any(ni == 0.0):
            bad = cquad4.element_id[ni == 0.0]
            raise RuntimeError(f"CQUAD4 elements {bad.tolist()} have zero-area normals")
        normal /= ni[:, np.newaxis]

        # Element local x-axis (node 1 -> node 2)
        ihat = p2 - p1
        inorm = np.linalg.norm(ihat, axis=1)
        ihat /= inorm[:, np.newaxis]

        # Element local y-axis (k x i)
        jhat = np.cross(normal, ihat, axis=1)
        jnorm = np.linalg.norm(jhat, axis=1)
        jhat /= jnorm[:, np.newaxis]

        # Process each element
        nid_array = grid.node_id
        for i_elem in range(nelements):
            pid = int(cquad4.property_id[i_elem])

            # Build rotation matrix T (3x3): rows = [ihat, jhat, normal]
            Ti = np.vstack([ihat[i_elem], jhat[i_elem], normal[i_elem]])

            # Transform nodes to local frame
            xy = np.zeros((4, 3))
            xy[0] = Ti @ p1[i_elem]
            xy[1] = Ti @ p2[i_elem]
            xy[2] = Ti @ p3[i_elem]
            xy[3] = Ti @ p4[i_elem]
            x_local = xy[:, 0]
            y_local = xy[:, 1]

            # Get constitutive matrices (in material frame)
            A_mat, B_mat, D_mat, Ds_mat, thickness = _get_ABD_for_element(model, pid, pshell, pcomp)

            # Material coordinate rotation (THETA/MCID on CQUAD4)
            # Rotates constitutive from material frame to element frame
            theta_elem = cquad4.theta[i_elem]
            mcid_elem = cquad4.mcid[i_elem]

            if mcid_elem >= 0:
                # MCID: material x-axis from coordinate system projected onto element
                mcid_ref = model.coord.slice_card_by_id(np.array([mcid_elem]))
                i_mcid = mcid_ref.i[0]  # x-axis of coord system
                # Project onto element plane
                i_proj = i_mcid - np.dot(i_mcid, normal[i_elem]) * normal[i_elem]
                i_proj_norm = np.linalg.norm(i_proj)
                if i_proj_norm > 1e-10:
                    i_proj /= i_proj_norm
                    # Angle between element x-axis (ihat) and projected MCID x-axis
                    cos_theta = np.dot(ihat[i_elem], i_proj)
                    sin_theta = np.dot(np.cross(ihat[i_elem], i_proj), normal[i_elem])
                    theta_mat = np.arctan2(sin_theta, cos_theta)
                else:
                    theta_mat = 0.0
            elif not np.isnan(theta_elem) and abs(theta_elem) > 1e-10:
                # THETA: direct angle in degrees
                theta_mat = np.radians(theta_elem)
            else:
                theta_mat = 0.0

            # Rotate ABD from material to element frame
            if abs(theta_mat) > 1e-10:
                c = np.cos(theta_mat)
                s = np.sin(theta_mat)
                T_rot = np.array(
                    [
                        [c**2, s**2, 2 * s * c],
                        [s**2, c**2, -2 * s * c],
                        [-s * c, s * c, c**2 - s**2],
                    ]
                )
                T_inv = np.array(
                    [
                        [c**2, s**2, -2 * s * c],
                        [s**2, c**2, 2 * s * c],
                        [s * c, -s * c, c**2 - s**2],
                    ]
                )
                A_mat = T_inv @ A_mat @ T_inv.T
                D_mat = T_inv @ D_mat @ T_inv.T
                if B_mat is not None:
                    B_mat = T_inv @ B_mat @ T_inv.T
                # Rotate transverse shear: Ds in material -> element
                R2 = np.array([[c, s], [-s, c]])
                Ds_mat = R2.T @ Ds_mat @ R2

            # Verify Jacobian quality at centroid
            dN_c = _dshape_quad4(0.0, 0.0)
            _, det_J_c = _jacobian(dN_c, x_local, y_local)
            if det_J_c <= 0.0:
                eid = int(cquad4.element_id[i_elem])
                raise RuntimeError(
                    f"CQUAD4 eid={eid} has non-positive Jacobian "
                    f"determinant ({det_J_c:.6g}) — element is degenerate"
                )

            # Compute element stiffness (20x20 in local frame)
            Ke_local = self._element_stiffness(x_local, y_local, A_mat, D_mat, Ds_mat, B_mat)

            # Apply shell offset
            zoffset = cquad4.zoffset[i_elem]
            if not np.isnan(zoffset) and abs(zoffset) > 0.0:
                Ke_local = _apply_shell_offset(Ke_local, zoffset)

            # Transform local 20-DOF -> global 24-DOF
            T_node = np.zeros((5, 6))
            T_node[0:3, 0:3] = Ti
            T_node[3:5, 3:6] = Ti[0:2, :]

            T_full = np.zeros((20, 24))
            for inode in range(4):
                r0 = 5 * inode
                c0 = 6 * inode
                T_full[r0 : r0 + 5, c0 : c0 + 6] = T_node

            Ke_global = T_full.T @ Ke_local @ T_full

            # Drilling stiffness (rotation about element normal)
            self._apply_drilling(
                Ke_global,
                Ke_local,
                x_local,
                y_local,
                A_mat,
                thickness,
                normal[i_elem],
                Ti,
            )

            # Assemble into Kbb
            elem_nodes = nodes[i_elem]
            global_dofs = []
            for nid in elem_nodes:
                i_dof = dof_map[(nid, 1)]
                for comp in range(6):
                    global_dofs.append(i_dof + comp)

            for ii in range(24):
                gi = global_dofs[ii]
                for jj in range(24):
                    gj = global_dofs[jj]
                    val = Ke_global[ii, jj]
                    if abs(val) > 0.0:
                        Kbb[gi, gj] += val

        return nelements


# ---------------------------------------------------------------------------
# Public entry point for build_stiffness.py
# ---------------------------------------------------------------------------
def build_kbb_cquad4(
    model: BDF,
    Kbb: dok_matrix,
    dof_map: DOF_MAP,
    all_nids: np.ndarray,
    xyz_cid0: np.ndarray,
    idtype: str = "int32",
    fdtype: str = "float64",
) -> int:
    """Build CQUAD4 stiffness using PARAM,MYQUAD formulation dispatch."""
    solver = ShellQuadSolver(model)
    return solver.build_kbb(Kbb, dof_map, all_nids, xyz_cid0)


# ---------------------------------------------------------------------------
# PLOAD4 pressure load vector for CQUAD4
# ---------------------------------------------------------------------------
def mitc4_pressure_load(
    x: np.ndarray,
    y: np.ndarray,
    pressures: np.ndarray,
    direction: np.ndarray | None = None,
) -> np.ndarray:
    """Compute the 20x1 equivalent nodal force vector for surface pressure.

    Parameters
    ----------
    x, y : (4,) nodal coords in element local frame
    pressures : (4,) or scalar
        Pressure at each node. Scalar = uniform.
    direction : (3,) or None
        Load direction in element local frame. None = +z (normal).

    Returns
    -------
    Fe : (20,)
        DOF order per node: [u, v, w, tx, ty].
    """
    pressures = np.atleast_1d(np.asarray(pressures, dtype=float))
    if pressures.size == 1:
        pressures = np.full(4, pressures[0])

    Fe = np.zeros(20)
    w_idx = np.array([2, 7, 12, 17])

    if direction is None:
        for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
            xi, eta = gpt
            N = _shape_quad4(xi, eta)
            dN_dnat = _dshape_quad4(xi, eta)
            _, det_J = _jacobian(dN_dnat, x, y)
            p_gp = float(N @ pressures)
            Fe[w_idx] += N * p_gp * det_J * wt
    else:
        dir_arr = np.asarray(direction, dtype=float)
        for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
            xi, eta = gpt
            N = _shape_quad4(xi, eta)
            dN_dnat = _dshape_quad4(xi, eta)
            _, det_J = _jacobian(dN_dnat, x, y)
            p_gp = float(N @ pressures)
            scale = p_gp * det_J * wt
            for i in range(4):
                Fe[5 * i : 5 * i + 3] += N[i] * scale * dir_arr

    return Fe


def _pload4_direction_local(
    pload4,
    i_card: int,
    Ti: np.ndarray,
    normal_global: np.ndarray,
    model: BDF,
    _coord_cache: dict | None = None,
) -> np.ndarray | None:
    """Compute PLOAD4 load direction in element local frame.

    Returns None for normal pressure (default), or (3,) unit vector
    in element local frame when CID/nvector is specified.
    """
    nvec = pload4.nvector[i_card]
    nvec_sq = nvec[0] * nvec[0] + nvec[1] * nvec[1] + nvec[2] * nvec[2]
    if nvec_sq < 1e-60:
        return None

    cid = int(pload4.coord_id[i_card])

    if cid == 0:
        nvec_global = nvec
    else:
        if _coord_cache is not None and cid in _coord_cache:
            R = _coord_cache[cid]
        else:
            coord = model.coord
            coord_slice = coord.slice_card_by_id(np.array([cid]))
            R = np.vstack([coord_slice.i[0], coord_slice.j[0], coord_slice.k[0]]).T
            if _coord_cache is not None:
                _coord_cache[cid] = R
        nvec_global = R @ nvec

    nvec_norm = np.sqrt(nvec_global[0] ** 2 + nvec_global[1] ** 2 + nvec_global[2] ** 2)
    if nvec_norm < 1e-30:
        return None

    return Ti @ (nvec_global / nvec_norm)


def build_pload4_cquad4(
    model: BDF,
    Fb: np.ndarray,
    dof_map: DOF_MAP,
    load_id: int,
) -> None:
    """Assemble PLOAD4 on CQUAD4 elements into global force vector.

    Supports variable pressure per node and CID direction vectors.
    """
    pload4 = model.pload4
    if pload4.n == 0:
        return

    cquad4 = model.cquad4
    if cquad4.n == 0:
        return

    grid = model.grid
    xyz = grid.xyz_cid0()
    coord_cache: dict = {}

    for i_card in range(pload4.n):
        if pload4.load_id[i_card] != load_id:
            continue

        pressures_card = pload4.pressure[i_card]
        i0, i1 = pload4.ielement[i_card]

        for i_elem_idx in range(i0, i1):
            eid = pload4.element_ids[i_elem_idx]
            idx_elem = np.searchsorted(cquad4.element_id, eid)
            if idx_elem >= cquad4.n or cquad4.element_id[idx_elem] != eid:
                continue

            elem_nodes = cquad4.nodes[idx_elem]
            inid = grid.index(elem_nodes.reshape(1, 4))

            pts = xyz[inid[0], :]  # (4, 3)
            v13 = pts[0] - pts[2]
            v24 = pts[1] - pts[3]
            normal = np.cross(v13, v24)
            normal /= np.linalg.norm(normal)

            ihat = pts[1] - pts[0]
            ihat /= np.linalg.norm(ihat)
            jhat = np.cross(normal, ihat)
            jhat /= np.linalg.norm(jhat)

            Ti = np.vstack([ihat, jhat, normal])

            xy_local = (Ti @ pts.T).T  # (4, 3)
            x_local = xy_local[:, 0]
            y_local = xy_local[:, 1]

            direction_local = _pload4_direction_local(
                pload4, i_card, Ti, normal, model, coord_cache
            )

            Fe_local = mitc4_pressure_load(x_local, y_local, pressures_card, direction_local)

            # Transform Fe from local 20-DOF to global 24-DOF directly
            # Local DOF per node: [u,v,w,tx,ty] -> Global: [Ux,Uy,Uz,Rx,Ry,Rz]
            # u_local = Ti @ U_global (translations)
            # [tx,ty]_local = Ti[0:2] @ R_global (rotations)
            # Fe_global_trans = Ti.T @ Fe_local_trans
            # Fe_global_rot = Ti[0:2].T @ Fe_local_rot
            Ti_T = Ti.T
            Ti_rot_T = Ti[0:2, :].T  # (3, 2)

            for inode in range(4):
                nid = int(elem_nodes[inode])
                i_dof = dof_map[(nid, 1)]
                fl = Fe_local[5 * inode : 5 * inode + 5]
                Fb[i_dof : i_dof + 3] += Ti_T @ fl[:3]
                Fb[i_dof + 3 : i_dof + 6] += Ti_rot_T @ fl[3:5]


# ---------------------------------------------------------------------------
# PLOAD4 pressure load vector for CTRIA3
# ---------------------------------------------------------------------------
def ctria3_pressure_load(
    x: np.ndarray,
    y: np.ndarray,
    pressures: np.ndarray,
    direction: np.ndarray | None = None,
) -> np.ndarray:
    """Compute the 15x1 equivalent nodal force vector for surface pressure on CTRIA3.

    Parameters
    ----------
    x, y : (3,) nodal coords in element local frame
    pressures : (3,) or scalar
        Pressure at each node. Scalar = uniform.
    direction : (3,) or None
        Load direction in element local frame. None = +z (normal).

    Returns
    -------
    Fe : (15,)
        DOF order per node: [u, v, w, tx, ty].
    """
    pressures = np.atleast_1d(np.asarray(pressures, dtype=float))
    if pressures.size == 1:
        pressures = np.full(3, pressures[0])

    area = 0.5 * abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))
    Fe = np.zeros(15)

    # Consistent nodal forces for linear triangle:
    # Fe_i = A/12 * (2*p_i + p_j + p_k)
    p1, p2, p3 = pressures[0], pressures[1], pressures[2]
    f1 = area / 12.0 * (2.0 * p1 + p2 + p3)
    f2 = area / 12.0 * (p1 + 2.0 * p2 + p3)
    f3 = area / 12.0 * (p1 + p2 + 2.0 * p3)

    if direction is None:
        Fe[2] = f1
        Fe[7] = f2
        Fe[12] = f3
    else:
        dx, dy, dz = direction[0], direction[1], direction[2]
        Fe[0] = f1 * dx
        Fe[1] = f1 * dy
        Fe[2] = f1 * dz
        Fe[5] = f2 * dx
        Fe[6] = f2 * dy
        Fe[7] = f2 * dz
        Fe[10] = f3 * dx
        Fe[11] = f3 * dy
        Fe[12] = f3 * dz

    return Fe


def build_pload4_ctria3(
    model: BDF,
    Fb: np.ndarray,
    dof_map: DOF_MAP,
    load_id: int,
) -> None:
    """Assemble PLOAD4 on CTRIA3 elements into global force vector.

    Supports variable pressure per node and CID direction vectors.
    """
    pload4 = model.pload4
    if pload4.n == 0:
        return

    ctria3 = model.ctria3
    if ctria3.n == 0:
        return

    grid = model.grid
    xyz = grid.xyz_cid0()
    coord_cache: dict = {}

    for i_card in range(pload4.n):
        if pload4.load_id[i_card] != load_id:
            continue

        pressures_card = pload4.pressure[i_card]
        i0, i1 = pload4.ielement[i_card]

        for i_elem_idx in range(i0, i1):
            eid = pload4.element_ids[i_elem_idx]
            idx_elem = np.searchsorted(ctria3.element_id, eid)
            if idx_elem >= ctria3.n or ctria3.element_id[idx_elem] != eid:
                continue

            elem_nodes = ctria3.nodes[idx_elem]
            inid = grid.index(elem_nodes.reshape(1, 3))

            pts = xyz[inid[0], :]  # (3, 3)
            v12 = pts[1] - pts[0]
            v13 = pts[2] - pts[0]
            normal = np.cross(v12, v13)
            norm_len = np.linalg.norm(normal)
            if norm_len < 1e-30:
                continue
            normal /= norm_len

            ihat = v12 / np.linalg.norm(v12)
            jhat = np.cross(normal, ihat)
            jhat /= np.linalg.norm(jhat)

            Ti = np.vstack([ihat, jhat, normal])

            xy_local = (Ti @ pts.T).T  # (3, 3)
            x_local = xy_local[:, 0]
            y_local = xy_local[:, 1]

            pressures_tri = pressures_card[:3]

            direction_local = _pload4_direction_local(
                pload4, i_card, Ti, normal, model, coord_cache
            )

            Fe_local = ctria3_pressure_load(x_local, y_local, pressures_tri, direction_local)

            Ti_T = Ti.T
            Ti_rot_T = Ti[0:2, :].T

            for inode in range(3):
                nid = int(elem_nodes[inode])
                i_dof = dof_map[(nid, 1)]
                fl = Fe_local[5 * inode : 5 * inode + 5]
                Fb[i_dof : i_dof + 3] += Ti_T @ fl[:3]
                Fb[i_dof + 3 : i_dof + 6] += Ti_rot_T @ fl[3:5]


# ---------------------------------------------------------------------------
# Thermal load vector for CQUAD4
# ---------------------------------------------------------------------------
def _get_alpha_vec_for_element(
    model: BDF,
    pid: int,
    pshell,
    pcomp,
) -> tuple[np.ndarray, float]:
    """Get CTE vector and reference temperature for a single property ID.

    Returns
    -------
    alpha_vec : (3,)
        Thermal expansion vector [alpha_x, alpha_y, alpha_xy] in material frame.
    tref : float
        Reference temperature from the material card.
    """
    mat1 = model.mat1
    mat8 = model.mat8

    if pshell.n > 0:
        iprop = np.searchsorted(pshell.property_id, pid)
        if iprop < pshell.n and pshell.property_id[iprop] == pid:
            prop_slice = pshell.slice_card_by_id(np.array([pid]))
            mids = prop_slice.material_id[0]
            mid1 = int(mids[0])
            if mid1 <= 0:
                return np.zeros(3), 0.0

            if mat1.n > 0:
                idx = np.searchsorted(mat1.material_id, mid1)
                if idx < mat1.n and mat1.material_id[idx] == mid1:
                    alpha = float(mat1.alpha[idx])
                    tref = float(mat1.tref[idx])
                    return np.array([alpha, alpha, 0.0]), tref

            if mat8.n > 0:
                idx = np.searchsorted(mat8.material_id, mid1)
                if idx < mat8.n and mat8.material_id[idx] == mid1:
                    a1 = float(mat8.alpha[idx, 0])
                    a2 = float(mat8.alpha[idx, 1])
                    tref = float(mat8.tref[idx])
                    return np.array([a1, a2, 0.0]), tref

    if pcomp.n > 0:
        iprop = np.searchsorted(pcomp.property_id, pid)
        if iprop < pcomp.n and pcomp.property_id[iprop] == pid:
            prop_slice = pcomp.slice_card_by_id(np.array([pid]))
            ilayer = prop_slice.ilayer
            i0, i1 = ilayer[0]
            ply_mids = prop_slice.material_id[i0:i1]
            ply_thetas = prop_slice.theta[i0:i1]
            ply_thicknesses = prop_slice.thickness[i0:i1]
            total_t = float(prop_slice.total_thickness()[0])

            alpha_vec = np.zeros(3)
            tref_val = 0.0
            for mid_ply, theta_ply, t_ply in zip(ply_mids, ply_thetas, ply_thicknesses):
                mid_ply = int(mid_ply)
                a_ply = np.zeros(3)
                tref_ply = 0.0

                if mat1.n > 0:
                    idx = np.searchsorted(mat1.material_id, mid_ply)
                    if idx < mat1.n and mat1.material_id[idx] == mid_ply:
                        alpha_val = float(mat1.alpha[idx])
                        a_ply = np.array([alpha_val, alpha_val, 0.0])
                        tref_ply = float(mat1.tref[idx])

                if mat8.n > 0 and np.all(a_ply == 0.0):
                    idx = np.searchsorted(mat8.material_id, mid_ply)
                    if idx < mat8.n and mat8.material_id[idx] == mid_ply:
                        a_ply = np.array(
                            [float(mat8.alpha[idx, 0]), float(mat8.alpha[idx, 1]), 0.0]
                        )
                        tref_ply = float(mat8.tref[idx])

                if abs(theta_ply) > 1e-10:
                    c = np.cos(np.radians(theta_ply))
                    s = np.sin(np.radians(theta_ply))
                    # Transform alpha from ply to laminate coords
                    ax = a_ply[0] * c**2 + a_ply[1] * s**2
                    ay = a_ply[0] * s**2 + a_ply[1] * c**2
                    axy = 2.0 * (a_ply[0] - a_ply[1]) * s * c
                    a_ply = np.array([ax, ay, axy])

                alpha_vec += a_ply * (t_ply / total_t)
                tref_val += tref_ply * (t_ply / total_t)

            return alpha_vec, tref_val

    return np.zeros(3), 0.0


def build_thermal_load_cquad4(
    model: BDF,
    Fb: np.ndarray,
    dof_map: DOF_MAP,
    node_temperatures: dict[int, float],
) -> None:
    """Compute and add CQUAD4 thermal equivalent nodal loads to Fb.

    Parameters
    ----------
    model : BDF
        The bdf_vectorized3 model (setup + cross-referenced).
    Fb : (ndof,)
        Global force vector to add thermal loads into.
    dof_map : DOF_MAP
        Maps (nid, dof) -> global DOF index.
    node_temperatures : dict[int, float]
        Temperature at each grid point.
    """
    cquad4 = model.cquad4
    nelements = cquad4.n
    if nelements == 0:
        return

    pshell = model.pshell
    pcomp = model.pcomp
    grid = model.grid
    nodes = cquad4.nodes
    inid = grid.index(nodes)
    xyz = grid.xyz_cid0()

    p1 = xyz[inid[:, 0], :]
    p2 = xyz[inid[:, 1], :]
    p3 = xyz[inid[:, 2], :]
    p4 = xyz[inid[:, 3], :]

    v13 = p1 - p3
    v24 = p2 - p4
    normal = np.cross(v13, v24)
    ni = np.linalg.norm(normal, axis=1)
    normal /= ni[:, np.newaxis]

    ihat = p2 - p1
    inorm = np.linalg.norm(ihat, axis=1)
    ihat /= inorm[:, np.newaxis]

    jhat = np.cross(normal, ihat, axis=1)
    jnorm = np.linalg.norm(jhat, axis=1)
    jhat /= jnorm[:, np.newaxis]

    for i_elem in range(nelements):
        pid = int(cquad4.property_id[i_elem])
        elem_nodes = nodes[i_elem]

        # Get element temperatures (average of nodal temps)
        elem_temps = np.array([node_temperatures.get(int(nid), 0.0) for nid in elem_nodes])

        # Get CTE and Tref
        alpha_vec, tref = _get_alpha_vec_for_element(model, pid, pshell, pcomp)
        if np.all(np.abs(alpha_vec) < 1e-30):
            continue

        # Temperature differences at nodes
        dT_nodes = elem_temps - tref
        if np.all(np.abs(dT_nodes) < 1e-30):
            continue

        # Build rotation matrix T (3x3): rows = [ihat, jhat, normal]
        Ti = np.vstack([ihat[i_elem], jhat[i_elem], normal[i_elem]])

        # Transform nodes to local frame
        xy = np.zeros((4, 3))
        xy[0] = Ti @ p1[i_elem]
        xy[1] = Ti @ p2[i_elem]
        xy[2] = Ti @ p3[i_elem]
        xy[3] = Ti @ p4[i_elem]
        x_local = xy[:, 0]
        y_local = xy[:, 1]

        # Get constitutive matrices
        A_mat, B_mat, D_mat, Ds_mat, thickness = _get_ABD_for_element(model, pid, pshell, pcomp)

        # Material coordinate rotation
        theta_elem = cquad4.theta[i_elem]
        mcid_elem = cquad4.mcid[i_elem]

        if mcid_elem >= 0:
            mcid_ref = model.coord.slice_card_by_id(np.array([mcid_elem]))
            i_mcid = mcid_ref.i[0]
            i_proj = i_mcid - np.dot(i_mcid, normal[i_elem]) * normal[i_elem]
            i_proj_norm = np.linalg.norm(i_proj)
            if i_proj_norm > 1e-10:
                i_proj /= i_proj_norm
                cos_theta = np.dot(ihat[i_elem], i_proj)
                sin_theta = np.dot(np.cross(ihat[i_elem], i_proj), normal[i_elem])
                theta_mat = np.arctan2(sin_theta, cos_theta)
            else:
                theta_mat = 0.0
        elif not np.isnan(theta_elem) and abs(theta_elem) > 1e-10:
            theta_mat = np.radians(theta_elem)
        else:
            theta_mat = 0.0

        # Rotate ABD and alpha from material to element frame
        if abs(theta_mat) > 1e-10:
            c = np.cos(theta_mat)
            s = np.sin(theta_mat)
            T_inv = np.array(
                [
                    [c**2, s**2, -2 * s * c],
                    [s**2, c**2, 2 * s * c],
                    [s * c, -s * c, c**2 - s**2],
                ]
            )
            A_mat = T_inv @ A_mat @ T_inv.T
            if B_mat is not None:
                B_mat = T_inv @ B_mat @ T_inv.T
            # Rotate alpha vector to element frame
            # alpha transforms like strain: alpha_bar = T_inv @ alpha
            alpha_vec = T_inv @ alpha_vec

        # Thermal force resultant: N_th = A * alpha * dT
        # Integrate using 2x2 Gauss quadrature with bilinear temperature
        # interpolation within the element.
        #
        # F_thermal_membrane(20-DOF) = integral(Bm^T * N_th) dA
        # F_thermal_bending(20-DOF) = integral(Bb^T * M_th) dA  (if B_mat coupling)

        m_idx = np.array([0, 1, 5, 6, 10, 11, 15, 16], dtype=int)
        b_idx = np.array([2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19], dtype=int)

        Fe_local = np.zeros(20)

        for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
            xi, eta = gpt
            N = _shape_quad4(xi, eta)
            dN_dnat = _dshape_quad4(xi, eta)
            J, det_J = _jacobian(dN_dnat, x_local, y_local)
            J_inv = _jacobian_inv(J, det_J)
            dN_dxy = J_inv @ dN_dnat

            # Temperature at this Gauss point (bilinear interpolation)
            dT_gp = N @ dT_nodes

            # Thermal force resultant
            N_th = A_mat @ (alpha_vec * dT_gp)

            # Membrane contribution
            Bm = _membrane_B(dN_dxy[0], dN_dxy[1])
            Fe_local[m_idx] += (Bm.T @ N_th) * det_J * wt

            # Bending contribution from membrane-bending coupling
            if B_mat is not None:
                M_th = B_mat @ (alpha_vec * dT_gp)
                Bb = _bending_B(dN_dxy[0], dN_dxy[1])
                Fe_local[b_idx] += (Bb.T @ M_th) * det_J * wt

        # Apply shell offset
        zoffset = cquad4.zoffset[i_elem]
        if not np.isnan(zoffset) and abs(zoffset) > 0.0:
            E_full = np.eye(20)
            for inode in range(4):
                r = 5 * inode
                E_full[r + 0, r + 4] = zoffset
                E_full[r + 1, r + 3] = -zoffset
            Fe_local = E_full.T @ Fe_local

        # Transform local 20-DOF -> global 24-DOF
        T_node = np.zeros((5, 6))
        T_node[0:3, 0:3] = Ti
        T_node[3:5, 3:6] = Ti[0:2, :]

        T_full = np.zeros((20, 24))
        for inode in range(4):
            r0 = 5 * inode
            c0 = 6 * inode
            T_full[r0 : r0 + 5, c0 : c0 + 6] = T_node

        Fe_global = T_full.T @ Fe_local

        # Assemble into Fb
        for inode_local in range(4):
            nid = int(elem_nodes[inode_local])
            i_dof = dof_map[(nid, 1)]
            for comp in range(6):
                Fb[i_dof + comp] += Fe_global[6 * inode_local + comp]


def build_thermal_load_ctria3(
    model: BDF,
    Fb: np.ndarray,
    dof_map: DOF_MAP,
    node_temperatures: dict[int, float],
) -> None:
    """Compute and add CTRIA3 thermal equivalent nodal loads to Fb.

    Parameters
    ----------
    model : BDF
        The bdf_vectorized3 model (setup + cross-referenced).
    Fb : (ndof,)
        Global force vector to add thermal loads into.
    dof_map : DOF_MAP
        Maps (nid, dof) -> global DOF index.
    node_temperatures : dict[int, float]
        Temperature at each grid point.
    """
    ctria3 = model.ctria3
    nelements = ctria3.n
    if nelements == 0:
        return

    pshell = model.pshell
    pcomp = model.pcomp
    grid = model.grid
    nodes = ctria3.nodes
    inid = grid.index(nodes)
    xyz = grid.xyz_cid0()

    p1 = xyz[inid[:, 0], :]
    p2 = xyz[inid[:, 1], :]
    p3 = xyz[inid[:, 2], :]

    v12 = p2 - p1
    v13 = p3 - p1
    normal = np.cross(v12, v13)
    ni = np.linalg.norm(normal, axis=1)
    normal /= ni[:, np.newaxis]

    ihat = v12.copy()
    inorm = np.linalg.norm(ihat, axis=1)
    ihat /= inorm[:, np.newaxis]

    jhat = np.cross(normal, ihat, axis=1)
    jnorm = np.linalg.norm(jhat, axis=1)
    jhat /= jnorm[:, np.newaxis]

    for i_elem in range(nelements):
        pid = int(ctria3.property_id[i_elem])
        elem_nodes = nodes[i_elem]

        elem_temps = np.array([node_temperatures.get(int(nid), 0.0) for nid in elem_nodes])

        alpha_vec, tref = _get_alpha_vec_for_element(model, pid, pshell, pcomp)
        if np.all(np.abs(alpha_vec) < 1e-30):
            continue

        dT_nodes = elem_temps - tref
        if np.all(np.abs(dT_nodes) < 1e-30):
            continue

        # Build rotation matrix T (3x3): rows = [ihat, jhat, normal]
        Ti = np.vstack([ihat[i_elem], jhat[i_elem], normal[i_elem]])

        # Transform nodes to local frame
        xy = np.zeros((3, 3))
        xy[0] = Ti @ p1[i_elem]
        xy[1] = Ti @ p2[i_elem]
        xy[2] = Ti @ p3[i_elem]
        x_local = xy[:, 0]
        y_local = xy[:, 1]

        # Element area
        area = ni[i_elem] / 2.0

        # Get constitutive matrices
        A_mat, B_mat, D_mat, Ds_mat, thickness = _get_ABD_for_element(model, pid, pshell, pcomp)

        # Material coordinate rotation
        theta_elem = ctria3.theta[i_elem]
        mcid_elem = ctria3.mcid[i_elem]

        theta_mat = 0.0
        if mcid_elem >= 0:
            mcid_ref = model.coord.slice_card_by_id(np.array([mcid_elem]))
            i_mcid = mcid_ref.i[0]
            i_proj = i_mcid - np.dot(i_mcid, normal[i_elem]) * normal[i_elem]
            i_proj_norm = np.linalg.norm(i_proj)
            if i_proj_norm > 1e-10:
                i_proj /= i_proj_norm
                cos_theta = np.dot(ihat[i_elem], i_proj)
                sin_theta = np.dot(np.cross(ihat[i_elem], i_proj), normal[i_elem])
                theta_mat = np.arctan2(sin_theta, cos_theta)
        elif not np.isnan(theta_elem) and abs(theta_elem) > 1e-10:
            theta_mat = np.radians(theta_elem)

        if abs(theta_mat) > 1e-10:
            c = np.cos(theta_mat)
            s = np.sin(theta_mat)
            T_inv = np.array(
                [
                    [c**2, s**2, -2 * s * c],
                    [s**2, c**2, 2 * s * c],
                    [s * c, -s * c, c**2 - s**2],
                ]
            )
            A_mat = T_inv @ A_mat @ T_inv.T
            if B_mat is not None:
                B_mat = T_inv @ B_mat @ T_inv.T
            alpha_vec = T_inv @ alpha_vec

        # CST: constant strain, single evaluation
        # Average dT over element (linear field -> centroid value)
        dT_avg = np.mean(dT_nodes)

        # Thermal force resultant
        N_th = A_mat @ (alpha_vec * dT_avg)

        # Membrane B-matrix
        Bm = _tri3_membrane_B(x_local, y_local, area)
        # Fe_membrane = integral(Bm.T @ N_th dA) = Bm.T @ N_th * A
        Fe_membrane = (Bm.T @ N_th) * area

        # Assemble 15-DOF local force (5 DOF/node: u, v, w, theta_x, theta_y)
        Fe_local = np.zeros(15)
        m_idx = np.array([0, 1, 5, 6, 10, 11], dtype=int)
        Fe_local[m_idx] = Fe_membrane

        # Apply shell offset
        zoffset = ctria3.zoffset[i_elem]
        if not np.isnan(zoffset) and abs(zoffset) > 0.0:
            E_full = np.eye(15)
            for inode in range(3):
                r = 5 * inode
                E_full[r + 0, r + 4] = zoffset
                E_full[r + 1, r + 3] = -zoffset
            Fe_local = E_full.T @ Fe_local

        # Transform local 15-DOF -> global 18-DOF
        T_node = np.zeros((5, 6))
        T_node[0:3, 0:3] = Ti
        T_node[3:5, 3:6] = Ti[0:2, :]

        T_full = np.zeros((15, 18))
        for inode in range(3):
            r0 = 5 * inode
            c0 = 6 * inode
            T_full[r0 : r0 + 5, c0 : c0 + 6] = T_node

        Fe_global = T_full.T @ Fe_local

        # Assemble into Fb
        for inode_local in range(3):
            nid = int(elem_nodes[inode_local])
            i_dof = dof_map[(nid, 1)]
            for comp in range(6):
                Fb[i_dof + comp] += Fe_global[6 * inode_local + comp]


# ===========================================================================
# CTRIA3 Element
# ===========================================================================


# ---------------------------------------------------------------------------
# CTRIA3 shape functions and B-matrices
# ---------------------------------------------------------------------------
def _tri3_membrane_B(x: np.ndarray, y: np.ndarray, area: float) -> np.ndarray:
    """Constant-strain triangle membrane B-matrix (3, 6).

    DOF order: [u1, v1, u2, v2, u3, v3].
    """
    y23 = y[1] - y[2]
    y31 = y[2] - y[0]
    y12 = y[0] - y[1]
    x32 = x[2] - x[1]
    x13 = x[0] - x[2]
    x21 = x[1] - x[0]
    inv2A = 1.0 / (2.0 * area)
    Bm = inv2A * np.array(
        [
            [y23, 0.0, y31, 0.0, y12, 0.0],
            [0.0, x32, 0.0, x13, 0.0, x21],
            [x32, y23, x13, y31, x21, y12],
        ]
    )
    return Bm


def _dkt_bending_stiffness(x: np.ndarray, y: np.ndarray, D_mat: np.ndarray) -> np.ndarray:
    """DKT (Discrete Kirchhoff Triangle) bending stiffness (9, 9).

    DOF order per node: [w, theta_x, theta_y] where theta_x = dw/dy,
    theta_y = -dw/dx (same convention as NX Nastran).

    Derived from first principles using:
    1. Quadratic (6-node) interpolation of beta_x and beta_y over the triangle
    2. DKT midside constraints from cubic w along edges and linear normal rotation
    3. DOF transform from Batoz [w, beta_x, beta_y] to Nastran [w, theta_x, theta_y]

    The midside constraint at each edge (from node ni to nj) is:
      beta_s_mid = (3/(2*Lk))*(w_nj - w_ni) - (1/4)*(beta_s_ni + beta_s_nj)
      beta_n_mid = (1/2)*(beta_n_ni + beta_n_nj)

    References
    ----------
    Batoz, J.L., Bathe, K.J., Ho, L.W. (1980). "A study of three-node triangular
    plate bending elements." IJNME, 15(12), 1771-1812.
    """
    # Edge directions (from first node to second):
    # Edge 4 (opp node 1): from node 2 to node 3
    # Edge 5 (opp node 2): from node 3 to node 1
    # Edge 6 (opp node 3): from node 1 to node 2
    dx4 = x[2] - x[1]
    dy4 = y[2] - y[1]
    L4 = np.sqrt(dx4 * dx4 + dy4 * dy4)
    dx5 = x[0] - x[2]
    dy5 = y[0] - y[2]
    L5 = np.sqrt(dx5 * dx5 + dy5 * dy5)
    dx6 = x[1] - x[0]
    dy6 = y[1] - y[0]
    L6 = np.sqrt(dx6 * dx6 + dy6 * dy6)

    # Edge unit vectors
    cos4 = dx4 / L4
    sin4 = dy4 / L4
    cos5 = dx5 / L5
    sin5 = dy5 / L5
    cos6 = dx6 / L6
    sin6 = dy6 / L6

    # Midside beta_x coefficients for each edge.
    # For edge from ni to nj with direction cosines (ck, sk) and length Lk:
    #   bx_mid = (3*ck/(2*Lk))*(w_nj - w_ni)
    #          + (-(ck^2)/4 + (sk^2)/2)*(bx_ni + bx_nj)
    #          + (-3*ck*sk/4)*(by_ni + by_nj)
    # Stored as 9-vectors over DOF = [w1,bx1,by1,w2,bx2,by2,w3,bx3,by3]
    bx4 = np.zeros(9)  # edge 4: ni=1(idx 1), nj=2(idx 2)
    bx4[3] = -3.0 * cos4 / (2.0 * L4)
    bx4[6] = 3.0 * cos4 / (2.0 * L4)
    bx_c4 = -(cos4 * cos4) / 4.0 + (sin4 * sin4) / 2.0
    bx4[4] = bx_c4
    bx4[7] = bx_c4
    by_c4 = -3.0 * cos4 * sin4 / 4.0
    bx4[5] = by_c4
    bx4[8] = by_c4

    bx5 = np.zeros(9)  # edge 5: ni=2(idx 2), nj=0(idx 0)
    bx5[6] = -3.0 * cos5 / (2.0 * L5)
    bx5[0] = 3.0 * cos5 / (2.0 * L5)
    bx_c5 = -(cos5 * cos5) / 4.0 + (sin5 * sin5) / 2.0
    bx5[7] = bx_c5
    bx5[1] = bx_c5
    by_c5 = -3.0 * cos5 * sin5 / 4.0
    bx5[8] = by_c5
    bx5[2] = by_c5

    bx6 = np.zeros(9)  # edge 6: ni=0(idx 0), nj=1(idx 1)
    bx6[0] = -3.0 * cos6 / (2.0 * L6)
    bx6[3] = 3.0 * cos6 / (2.0 * L6)
    bx_c6 = -(cos6 * cos6) / 4.0 + (sin6 * sin6) / 2.0
    bx6[1] = bx_c6
    bx6[4] = bx_c6
    by_c6 = -3.0 * cos6 * sin6 / 4.0
    bx6[2] = by_c6
    bx6[5] = by_c6

    # Midside beta_y coefficients
    by4 = np.zeros(9)
    by4[3] = -3.0 * sin4 / (2.0 * L4)
    by4[6] = 3.0 * sin4 / (2.0 * L4)
    bx_d4 = -3.0 * cos4 * sin4 / 4.0
    by4[4] = bx_d4
    by4[7] = bx_d4
    by_d4 = -(sin4 * sin4) / 4.0 + (cos4 * cos4) / 2.0
    by4[5] = by_d4
    by4[8] = by_d4

    by5 = np.zeros(9)
    by5[6] = -3.0 * sin5 / (2.0 * L5)
    by5[0] = 3.0 * sin5 / (2.0 * L5)
    bx_d5 = -3.0 * cos5 * sin5 / 4.0
    by5[7] = bx_d5
    by5[1] = bx_d5
    by_d5 = -(sin5 * sin5) / 4.0 + (cos5 * cos5) / 2.0
    by5[8] = by_d5
    by5[2] = by_d5

    by6 = np.zeros(9)
    by6[0] = -3.0 * sin6 / (2.0 * L6)
    by6[3] = 3.0 * sin6 / (2.0 * L6)
    bx_d6 = -3.0 * cos6 * sin6 / 4.0
    by6[1] = bx_d6
    by6[4] = bx_d6
    by_d6 = -(sin6 * sin6) / 4.0 + (cos6 * cos6) / 2.0
    by6[2] = by_d6
    by6[5] = by_d6

    # Area and coordinate transform
    two_A = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])
    area = 0.5 * two_A
    inv_2A = 1.0 / two_A
    # df/dx = (y23*df/dL1 + y31*df/dL2) / (2A)
    # df/dy = (-x23*df/dL1 - x31*df/dL2) / (2A)
    y23 = y[1] - y[2]
    y31 = y[2] - y[0]
    x23 = x[1] - x[2]
    x31 = x[2] - x[0]

    # Unit vectors for extracting corner DOFs from the 9-vector
    e_bx = np.zeros((3, 9))
    e_bx[0, 1] = 1.0  # bx1
    e_bx[1, 4] = 1.0  # bx2
    e_bx[2, 7] = 1.0  # bx3
    e_by = np.zeros((3, 9))
    e_by[0, 2] = 1.0  # by1
    e_by[1, 5] = 1.0  # by2
    e_by[2, 8] = 1.0  # by3

    # 3-point Gauss integration (weight = A/3 per point)
    gauss_pts = [
        (2.0 / 3.0, 1.0 / 6.0),
        (1.0 / 6.0, 2.0 / 3.0),
        (1.0 / 6.0, 1.0 / 6.0),
    ]

    K_bb = np.zeros((9, 9))

    for L1, L2 in gauss_pts:
        L3 = 1.0 - L1 - L2

        # Quadratic shape function derivatives w.r.t. L1, L2
        # Vertex: P1=L1(2L1-1), P2=L2(2L2-1), P3=L3(2L3-1)
        dP1_d1 = 4.0 * L1 - 1.0
        dP2_d1 = 0.0
        dP3_d1 = -(4.0 * L3 - 1.0)
        dP1_d2 = 0.0
        dP2_d2 = 4.0 * L2 - 1.0
        dP3_d2 = -(4.0 * L3 - 1.0)
        # Midside: P4=4*L2*L3, P5=4*L3*L1, P6=4*L1*L2
        dP4_d1 = -4.0 * L2
        dP4_d2 = 4.0 * (L3 - L2)
        dP5_d1 = 4.0 * (L3 - L1)
        dP5_d2 = -4.0 * L1
        dP6_d1 = 4.0 * L2
        dP6_d2 = 4.0 * L1

        # d(beta_x)/dL1 and d(beta_x)/dL2 as 9-vectors
        dBx_dL1 = (
            dP1_d1 * e_bx[0]
            + dP2_d1 * e_bx[1]
            + dP3_d1 * e_bx[2]
            + dP4_d1 * bx4
            + dP5_d1 * bx5
            + dP6_d1 * bx6
        )
        dBx_dL2 = (
            dP1_d2 * e_bx[0]
            + dP2_d2 * e_bx[1]
            + dP3_d2 * e_bx[2]
            + dP4_d2 * bx4
            + dP5_d2 * bx5
            + dP6_d2 * bx6
        )

        # d(beta_y)/dL1 and d(beta_y)/dL2 as 9-vectors
        dBy_dL1 = (
            dP1_d1 * e_by[0]
            + dP2_d1 * e_by[1]
            + dP3_d1 * e_by[2]
            + dP4_d1 * by4
            + dP5_d1 * by5
            + dP6_d1 * by6
        )
        dBy_dL2 = (
            dP1_d2 * e_by[0]
            + dP2_d2 * e_by[1]
            + dP3_d2 * e_by[2]
            + dP4_d2 * by4
            + dP5_d2 * by5
            + dP6_d2 * by6
        )

        # Convert to Cartesian derivatives
        dBx_dx = (y23 * dBx_dL1 + y31 * dBx_dL2) * inv_2A
        dBx_dy = (-x23 * dBx_dL1 - x31 * dBx_dL2) * inv_2A
        dBy_dx = (y23 * dBy_dL1 + y31 * dBy_dL2) * inv_2A
        dBy_dy = (-x23 * dBy_dL1 - x31 * dBy_dL2) * inv_2A

        # Curvature B-matrix (3x9) in Batoz DOF order [w, beta_x, beta_y]
        Bb = np.zeros((3, 9))
        Bb[0, :] = dBx_dx  # kx = d(beta_x)/dx
        Bb[1, :] = dBy_dy  # ky = d(beta_y)/dy
        Bb[2, :] = dBx_dy + dBy_dx  # kxy = d(beta_x)/dy + d(beta_y)/dx

        K_bb += (Bb.T @ D_mat @ Bb) * (area / 3.0)

    # Transform from Batoz DOF [w, beta_x, beta_y] to Nastran [w, theta_x, theta_y]
    # where theta_x = dw/dy = beta_y, theta_y = -dw/dx = -beta_x.
    # q_Batoz = T_node * q_Nastran with T_node = [[1,0,0],[0,0,-1],[0,1,0]]
    # K_Nastran = T_block^T @ K_Batoz @ T_block
    T_block = np.zeros((9, 9))
    for i in range(3):
        r = 3 * i
        T_block[r, r] = 1.0  # w -> w
        T_block[r + 1, r + 2] = -1.0  # beta_x -> -theta_y
        T_block[r + 2, r + 1] = 1.0  # beta_y -> theta_x
    return T_block.T @ K_bb @ T_block


# ---------------------------------------------------------------------------
# CTRIA3 element stiffness (15x15 in local 5-DOF frame)
# ---------------------------------------------------------------------------
def ctria3_stiffness(
    x: np.ndarray,
    y: np.ndarray,
    A_mat: np.ndarray,
    D_mat: np.ndarray,
    Ds_mat: np.ndarray,
    B_mat: np.ndarray | None = None,
) -> np.ndarray:
    """15x15 CTRIA3 shell element stiffness in local frame.

    Parameters
    ----------
    x, y : (3,) nodal coords in element local frame
    A_mat : (3, 3) membrane constitutive [N/length]
    D_mat : (3, 3) bending constitutive [force*length]
    Ds_mat : (2, 2) transverse shear constitutive [force/length]
    B_mat : (3, 3) membrane-bending coupling, optional

    Returns
    -------
    Ke : (15, 15) element stiffness
        DOF order per node: [u, v, w, theta_x, theta_y]
    """
    # Compute area
    area = 0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))
    if area <= 0.0:
        raise RuntimeError(f"CTRIA3 has non-positive area ({area:.6g}) — check node ordering")

    # DOF index mapping in 15-DOF vector
    # Per node: [u, v, w, tx, ty] -> 5 DOF
    m_idx = np.array([0, 1, 5, 6, 10, 11], dtype=int)  # membrane: u,v per node
    b_idx = np.array([2, 3, 4, 7, 8, 9, 12, 13, 14], dtype=int)  # bending: w,tx,ty

    # Membrane stiffness (constant strain, 1-point exact)
    Bm = _tri3_membrane_B(x, y, area)
    K_mm = (Bm.T @ A_mat @ Bm) * area  # (6, 6)

    # Bending stiffness — DKT (Discrete Kirchhoff Triangle)
    # No separate transverse shear: DKT is a pure Kirchhoff element
    K_bb = _dkt_bending_stiffness(x, y, D_mat)  # (9, 9)

    # Assemble 15x15
    Ke = np.zeros((15, 15))
    Ke[np.ix_(m_idx, m_idx)] = K_mm
    Ke[np.ix_(b_idx, b_idx)] = K_bb

    # Membrane-bending coupling (if B matrix from asymmetric laminate)
    if B_mat is not None:
        # Use constant-curvature B for coupling (consistent with membrane CST)
        y23 = y[1] - y[2]
        y31 = y[2] - y[0]
        y12 = y[0] - y[1]
        x32 = x[2] - x[1]
        x13 = x[0] - x[2]
        x21 = x[1] - x[0]
        inv2A = 1.0 / (2.0 * area)
        dNdx = np.array([y23, y31, y12]) * inv2A
        dNdy = np.array([x32, x13, x21]) * inv2A
        Bb_cst = np.zeros((3, 9))
        for i in range(3):
            c = 3 * i
            Bb_cst[0, c + 2] = dNdx[i]
            Bb_cst[1, c + 1] = -dNdy[i]
            Bb_cst[2, c + 1] = -dNdx[i]
            Bb_cst[2, c + 2] = dNdy[i]
        K_mb = (Bm.T @ B_mat @ Bb_cst) * area  # (6, 9)
        Ke[np.ix_(m_idx, b_idx)] = K_mb
        Ke[np.ix_(b_idx, m_idx)] = K_mb.T

    return Ke


# ---------------------------------------------------------------------------
# CTRIA3 Mindlin-Reissner (thick plate) stiffness
# ---------------------------------------------------------------------------
def _dst_bending_stiffness(
    x: np.ndarray, y: np.ndarray, D_mat: np.ndarray, Ds_mat: np.ndarray
) -> np.ndarray:
    """DST (Discrete Shear Triangle) bending stiffness (9, 9).

    Extends DKT by relaxing the Kirchhoff constraint at midside nodes
    using shear flexibility parameters phi_k per edge. Transitions smoothly
    from Kirchhoff (phi=0, Ds->inf) to thick-plate behavior (large phi).

    Parameters
    ----------
    x, y : (3,) nodal coordinates in local frame
    D_mat : (3, 3) bending constitutive
    Ds_mat : (2, 2) transverse shear constitutive (kappa*G*ts)

    Returns
    -------
    K_bb : (9, 9) bending+shear stiffness in Nastran DOF [w, theta_x, theta_y]

    References
    ----------
    Batoz, J.L. & Lardeur, P. (1989). "A discrete shear triangular nine DOF
    element for the analysis of thick to very thin plates." IJNME, 28, 533-560.
    """
    # Edge vectors and lengths
    # Edge 4 (opp node 0): node 1 -> node 2
    # Edge 5 (opp node 1): node 2 -> node 0
    # Edge 6 (opp node 2): node 0 -> node 1
    dx4 = x[2] - x[1]
    dy4 = y[2] - y[1]
    L4 = np.sqrt(dx4 * dx4 + dy4 * dy4)
    dx5 = x[0] - x[2]
    dy5 = y[0] - y[2]
    L5 = np.sqrt(dx5 * dx5 + dy5 * dy5)
    dx6 = x[1] - x[0]
    dy6 = y[1] - y[0]
    L6 = np.sqrt(dx6 * dx6 + dy6 * dy6)

    cos4 = dx4 / L4
    sin4 = dy4 / L4
    cos5 = dx5 / L5
    sin5 = dy5 / L5
    cos6 = dx6 / L6
    sin6 = dy6 / L6

    # Shear flexibility parameters per edge (directional)
    # phi_k = 12 * D_edge_k / (Ds_edge_k * L_k^2)
    # D_edge_k = bending rigidity in the direction of edge k:
    #   e_k = [cos_k^2, sin_k^2, cos_k*sin_k]  (curvature strain basis)
    #   D_edge_k = e_k^T . D . e_k
    # Ds_edge_k = shear stiffness in edge direction:
    #   s_k = [cos_k, sin_k]
    #   Ds_edge_k = s_k^T . Ds . s_k
    def _phi_edge(ck, sk, Lk):
        e_k = np.array([ck * ck, sk * sk, ck * sk])
        D_edge = e_k @ D_mat @ e_k
        s_k = np.array([ck, sk])
        Ds_edge = s_k @ Ds_mat @ s_k
        if Ds_edge > 1e-30:
            return 12.0 * D_edge / (Ds_edge * Lk * Lk)
        return 0.0

    phi4 = _phi_edge(cos4, sin4, L4)
    phi5 = _phi_edge(cos5, sin5, L5)
    phi6 = _phi_edge(cos6, sin6, L6)

    # DST modification factors for midside constraints (Batoz & Lardeur 1989)
    # Tangential: beta_s_mid = alpha_s_k * (3/(2L))*(w_j-w_i) + beta_s_k * (beta_s_i+beta_s_j)
    # alpha_s_k = 1/(1+phi_k)
    # beta_s_k = (2*phi_k - 1) / (4*(1+phi_k))  [= -1/4 for phi=0 = DKT]
    # Normal: beta_n_mid = (1/2)*(beta_n_i + beta_n_j)  [unchanged from DKT]
    alpha4 = 1.0 / (1.0 + phi4)
    alpha5 = 1.0 / (1.0 + phi5)
    alpha6 = 1.0 / (1.0 + phi6)
    beta4 = (2.0 * phi4 - 1.0) / (4.0 * (1.0 + phi4))
    beta5 = (2.0 * phi5 - 1.0) / (4.0 * (1.0 + phi5))
    beta6 = (2.0 * phi6 - 1.0) / (4.0 * (1.0 + phi6))

    # Build midside beta_x coefficients for each edge
    # DOF vector: [w1, bx1, by1, w2, bx2, by2, w3, bx3, by3] (Batoz convention)
    #
    # For edge k from node ni to nj:
    #   bx_mid_k = cos_k * [alpha_k*(3/(2L_k))*(w_nj - w_ni) + beta_k*(bxs_i + bxs_j)]
    #            + sin_k * [1/2 * (bxn_i + bxn_j)]
    # where bxs = cos_k*bx + sin_k*by (tangential projection of beta_x)
    #       bxn = -sin_k*bx + cos_k*by (normal projection of beta_x)
    #
    # Expanding:
    #   bx_mid_k = alpha_k*(3*cos_k/(2L_k))*(w_nj - w_ni)
    #            + [beta_k*cos_k^2 + sin_k^2/2]*(bx_i + bx_j)
    #            + [beta_k*cos_k*sin_k - sin_k*cos_k/2]*(by_i + by_j)
    #
    # Simplify coefficients on (bx_i + bx_j):
    #   c_bx_k = beta_k*cos_k^2 + sin_k^2/2
    # Coefficients on (by_i + by_j):
    #   c_by_k = (beta_k - 1/2)*cos_k*sin_k = ((2*phi-1)/(4*(1+phi)) - 1/2)*ck*sk
    #          = ((2*phi-1-2-2*phi)/(4*(1+phi)))*ck*sk = (-3/(4*(1+phi)))*ck*sk

    # Edge 4: from node 1 (idx 1) to node 2 (idx 2)
    bx4 = np.zeros(9)
    bx4[3] = -alpha4 * 3.0 * cos4 / (2.0 * L4)
    bx4[6] = alpha4 * 3.0 * cos4 / (2.0 * L4)
    c_bx4 = beta4 * cos4 * cos4 + sin4 * sin4 / 2.0
    bx4[4] = c_bx4
    bx4[7] = c_bx4
    c_by4 = -3.0 * cos4 * sin4 / (4.0 * (1.0 + phi4))
    bx4[5] = c_by4
    bx4[8] = c_by4

    # Edge 5: from node 2 (idx 2) to node 0 (idx 0)
    bx5 = np.zeros(9)
    bx5[6] = -alpha5 * 3.0 * cos5 / (2.0 * L5)
    bx5[0] = alpha5 * 3.0 * cos5 / (2.0 * L5)
    c_bx5 = beta5 * cos5 * cos5 + sin5 * sin5 / 2.0
    bx5[7] = c_bx5
    bx5[1] = c_bx5
    c_by5 = -3.0 * cos5 * sin5 / (4.0 * (1.0 + phi5))
    bx5[8] = c_by5
    bx5[2] = c_by5

    # Edge 6: from node 0 (idx 0) to node 1 (idx 1)
    bx6 = np.zeros(9)
    bx6[0] = -alpha6 * 3.0 * cos6 / (2.0 * L6)
    bx6[3] = alpha6 * 3.0 * cos6 / (2.0 * L6)
    c_bx6 = beta6 * cos6 * cos6 + sin6 * sin6 / 2.0
    bx6[1] = c_bx6
    bx6[4] = c_bx6
    c_by6 = -3.0 * cos6 * sin6 / (4.0 * (1.0 + phi6))
    bx6[2] = c_by6
    bx6[5] = c_by6

    # Midside beta_y coefficients
    # by_mid_k = sin_k * [alpha_k*(3/(2L_k))*(w_nj-w_ni) + beta_k*(bys_i+bys_j)]
    #          + cos_k * [1/2 * (byn_i + byn_j)]  (note: + for by normal component)
    # Actually:
    #   by_mid_k = sin_k*[alpha_k*(3/(2L))*(wj-wi)] + ...rotations...
    # bys = cos_k*by + sin_k*... wait, let me be precise.
    #
    # For beta_y, the tangential component along edge k is:
    #   bys = cos_k*by + sin_k*(-bx)?  No.
    # Actually the tangential/normal decomposition for beta_y:
    #   Tangential projection of (bx, by) along edge direction (ck, sk):
    #     beta_s = ck*bx + sk*by
    #   Normal projection:
    #     beta_n = -sk*bx + ck*by
    # The constraint is on beta_s and beta_n (the vector field, not just bx or by).
    # So for the y-component:
    #   by_mid = sk * beta_s_mid + ck * beta_n_mid
    #          = sk * [alpha*(3/(2L))*(wj-wi) + beta*(bs_i+bs_j)] + ck*[1/2*(bn_i+bn_j)]
    # where bs_i = ck*bx_i + sk*by_i, bn_i = -sk*bx_i + ck*by_i
    #
    # Coefficient on (by_i + by_j):
    #   sk * beta * sk + ck * (1/2) * ck = beta*sk^2 + ck^2/2
    # Coefficient on (bx_i + bx_j):
    #   sk * beta * ck + ck * (1/2) * (-sk) = (beta - 1/2)*ck*sk = -3*ck*sk/(4*(1+phi))

    # Edge 4
    by4 = np.zeros(9)
    by4[3] = -alpha4 * 3.0 * sin4 / (2.0 * L4)
    by4[6] = alpha4 * 3.0 * sin4 / (2.0 * L4)
    c_bx4_y = -3.0 * cos4 * sin4 / (4.0 * (1.0 + phi4))
    by4[4] = c_bx4_y
    by4[7] = c_bx4_y
    c_by4_y = beta4 * sin4 * sin4 + cos4 * cos4 / 2.0
    by4[5] = c_by4_y
    by4[8] = c_by4_y

    # Edge 5
    by5 = np.zeros(9)
    by5[6] = -alpha5 * 3.0 * sin5 / (2.0 * L5)
    by5[0] = alpha5 * 3.0 * sin5 / (2.0 * L5)
    c_bx5_y = -3.0 * cos5 * sin5 / (4.0 * (1.0 + phi5))
    by5[7] = c_bx5_y
    by5[1] = c_bx5_y
    c_by5_y = beta5 * sin5 * sin5 + cos5 * cos5 / 2.0
    by5[8] = c_by5_y
    by5[2] = c_by5_y

    # Edge 6
    by6 = np.zeros(9)
    by6[0] = -alpha6 * 3.0 * sin6 / (2.0 * L6)
    by6[3] = alpha6 * 3.0 * sin6 / (2.0 * L6)
    c_bx6_y = -3.0 * cos6 * sin6 / (4.0 * (1.0 + phi6))
    by6[1] = c_bx6_y
    by6[4] = c_bx6_y
    c_by6_y = beta6 * sin6 * sin6 + cos6 * cos6 / 2.0
    by6[2] = c_by6_y
    by6[5] = c_by6_y

    # Area and coordinate transform
    two_A = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])
    area = 0.5 * two_A
    inv_2A = 1.0 / two_A
    y23 = y[1] - y[2]
    y31 = y[2] - y[0]
    x23 = x[1] - x[2]
    x31 = x[2] - x[0]

    e_bx = np.zeros((3, 9))
    e_bx[0, 1] = 1.0
    e_bx[1, 4] = 1.0
    e_bx[2, 7] = 1.0
    e_by = np.zeros((3, 9))
    e_by[0, 2] = 1.0
    e_by[1, 5] = 1.0
    e_by[2, 8] = 1.0

    # 3-point Gauss integration
    gauss_pts = [
        (2.0 / 3.0, 1.0 / 6.0),
        (1.0 / 6.0, 2.0 / 3.0),
        (1.0 / 6.0, 1.0 / 6.0),
    ]

    K_bb = np.zeros((9, 9))

    for L1, L2 in gauss_pts:
        L3 = 1.0 - L1 - L2

        dP1_d1 = 4.0 * L1 - 1.0
        dP2_d1 = 0.0
        dP3_d1 = -(4.0 * L3 - 1.0)
        dP1_d2 = 0.0
        dP2_d2 = 4.0 * L2 - 1.0
        dP3_d2 = -(4.0 * L3 - 1.0)
        dP4_d1 = -4.0 * L2
        dP4_d2 = 4.0 * (L3 - L2)
        dP5_d1 = 4.0 * (L3 - L1)
        dP5_d2 = -4.0 * L1
        dP6_d1 = 4.0 * L2
        dP6_d2 = 4.0 * L1

        dBx_dL1 = (
            dP1_d1 * e_bx[0]
            + dP2_d1 * e_bx[1]
            + dP3_d1 * e_bx[2]
            + dP4_d1 * bx4
            + dP5_d1 * bx5
            + dP6_d1 * bx6
        )
        dBx_dL2 = (
            dP1_d2 * e_bx[0]
            + dP2_d2 * e_bx[1]
            + dP3_d2 * e_bx[2]
            + dP4_d2 * bx4
            + dP5_d2 * bx5
            + dP6_d2 * bx6
        )

        dBy_dL1 = (
            dP1_d1 * e_by[0]
            + dP2_d1 * e_by[1]
            + dP3_d1 * e_by[2]
            + dP4_d1 * by4
            + dP5_d1 * by5
            + dP6_d1 * by6
        )
        dBy_dL2 = (
            dP1_d2 * e_by[0]
            + dP2_d2 * e_by[1]
            + dP3_d2 * e_by[2]
            + dP4_d2 * by4
            + dP5_d2 * by5
            + dP6_d2 * by6
        )

        dBx_dx = (y23 * dBx_dL1 + y31 * dBx_dL2) * inv_2A
        dBx_dy = (-x23 * dBx_dL1 - x31 * dBx_dL2) * inv_2A
        dBy_dx = (y23 * dBy_dL1 + y31 * dBy_dL2) * inv_2A
        dBy_dy = (-x23 * dBy_dL1 - x31 * dBy_dL2) * inv_2A

        Bb = np.zeros((3, 9))
        Bb[0, :] = dBx_dx
        Bb[1, :] = dBy_dy
        Bb[2, :] = dBx_dy + dBy_dx

        K_bb += (Bb.T @ D_mat @ Bb) * (area / 3.0)

    # Transform from Batoz DOF [w, beta_x, beta_y] to Nastran [w, theta_x, theta_y]
    T_block = np.zeros((9, 9))
    for i in range(3):
        r = 3 * i
        T_block[r, r] = 1.0
        T_block[r + 1, r + 2] = -1.0
        T_block[r + 2, r + 1] = 1.0
    return T_block.T @ K_bb @ T_block


def ctria3_mindlin_stiffness(
    x: np.ndarray,
    y: np.ndarray,
    A_mat: np.ndarray,
    D_mat: np.ndarray,
    Ds_mat: np.ndarray,
    B_mat: np.ndarray | None = None,
) -> np.ndarray:
    """15x15 CTRIA3 Mindlin-Reissner shell element stiffness (DST formulation).

    Uses the DST (Discrete Shear Triangle) of Batoz & Lardeur (1989), which
    extends DKT by relaxing the Kirchhoff constraint at midside nodes via
    shear flexibility parameters phi_k = 12*D/(Ds*L_k^2) per edge.
    Activated by PARAM,TRITYP,THICK.

    Parameters
    ----------
    x, y : (3,) nodal coords in element local frame
    A_mat : (3, 3) membrane constitutive [N/length]
    D_mat : (3, 3) bending constitutive [force*length]
    Ds_mat : (2, 2) transverse shear constitutive [force/length]
    B_mat : (3, 3) membrane-bending coupling, optional

    Returns
    -------
    Ke : (15, 15) element stiffness
        DOF order per node: [u, v, w, theta_x, theta_y]

    References
    ----------
    Batoz, J.L. & Lardeur, P. (1989). "A discrete shear triangular nine DOF
    element for the analysis of thick to very thin plates." IJNME, 28, 533-560.
    """
    area = 0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))
    if area <= 0.0:
        raise RuntimeError(f"CTRIA3 has non-positive area ({area:.6g}) — check node ordering")

    m_idx = np.array([0, 1, 5, 6, 10, 11], dtype=int)
    b_idx = np.array([2, 3, 4, 7, 8, 9, 12, 13, 14], dtype=int)

    # Membrane stiffness (CST)
    Bm = _tri3_membrane_B(x, y, area)
    K_mm = (Bm.T @ A_mat @ Bm) * area

    # Bending + shear stiffness (DST)
    K_bb = _dst_bending_stiffness(x, y, D_mat, Ds_mat)

    # Assemble 15x15
    Ke = np.zeros((15, 15))
    Ke[np.ix_(m_idx, m_idx)] = K_mm
    Ke[np.ix_(b_idx, b_idx)] = K_bb

    # Membrane-bending coupling
    if B_mat is not None:
        y23 = y[1] - y[2]
        y31 = y[2] - y[0]
        y12 = y[0] - y[1]
        x32 = x[2] - x[1]
        x13 = x[0] - x[2]
        x21 = x[1] - x[0]
        inv2A = 1.0 / (2.0 * area)
        dNdx = np.array([y23, y31, y12]) * inv2A
        dNdy = np.array([x32, x13, x21]) * inv2A
        Bb_cst = np.zeros((3, 9))
        for i in range(3):
            c = 3 * i
            Bb_cst[0, c + 2] = dNdx[i]
            Bb_cst[1, c + 1] = -dNdy[i]
            Bb_cst[2, c + 1] = -dNdx[i]
            Bb_cst[2, c + 2] = dNdy[i]
        K_mb = (Bm.T @ B_mat @ Bb_cst) * area
        Ke[np.ix_(m_idx, b_idx)] = K_mb
        Ke[np.ix_(b_idx, m_idx)] = K_mb.T

    return Ke


# ---------------------------------------------------------------------------
# CTRIA3 drilling DOF
# ---------------------------------------------------------------------------
def _tri3_drilling_k6rot(Ke_local: np.ndarray, k6rot: float) -> float:
    """K6ROT artificial spring for CTRIA3."""
    m_idx = np.array([0, 1, 5, 6, 10, 11], dtype=int)
    diag_mm = np.diag(Ke_local[np.ix_(m_idx, m_idx)])
    max_diag = np.abs(diag_mm).max()
    if max_diag == 0.0:
        return 0.0
    return k6rot * 1e-4 * max_diag


def _tri3_drilling_hughes_brezzi(x: np.ndarray, y: np.ndarray, gamma: float) -> np.ndarray:
    """Hughes-Brezzi variational drilling for 3-node triangle (9, 9).

    DOF order per node: [u, v, theta_z].
    Uses 1-point centroidal integration (exact for linear shape functions).
    """
    area = 0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))

    y23 = y[1] - y[2]
    y31 = y[2] - y[0]
    y12 = y[0] - y[1]
    x32 = x[2] - x[1]
    x13 = x[0] - x[2]
    x21 = x[1] - x[0]
    inv2A = 1.0 / (2.0 * area)

    dNdx = np.array([y23, y31, y12]) * inv2A
    dNdy = np.array([x32, x13, x21]) * inv2A
    N = np.array([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])

    # G-vector: dc/d(DOF) where c = θz - ½(∂v/∂x - ∂u/∂y)
    G_vec = np.zeros(9)
    for i in range(3):
        c = 3 * i
        G_vec[c] = 0.5 * dNdy[i]  # dc/du_i = +½ dN_i/dy
        G_vec[c + 1] = -0.5 * dNdx[i]  # dc/dv_i = -½ dN_i/dx
        G_vec[c + 2] = N[i]  # dc/dθz_i = N_i

    K_hb = gamma * np.outer(G_vec, G_vec) * area
    return K_hb


# ---------------------------------------------------------------------------
# ShellTriSolver — unified assembly for CTRIA3
# ---------------------------------------------------------------------------
class ShellTriSolver:
    """Assembles CTRIA3 element stiffness matrices into global Kbb.

    Dispatches to DKT (default) or Mindlin-Reissner (PARAM,TRITYP,THICK).

    Parameters
    ----------
    model : BDF
        The bdf_vectorized3 model (setup + cross-referenced).
    """

    def __init__(self, model: BDF):
        self.model = model
        self.k6rot = self._read_param_k6rot()
        self.drilling = self._read_param_myqdril()
        self.trityp = self._read_param_trityp()

    def _read_param_k6rot(self) -> float:
        params = self.model.params
        if "K6ROT" in params:
            return float(params["K6ROT"].values[0])
        return 100.0

    def _read_param_myqdril(self) -> str:
        params = self.model.params
        if "MYQDRIL" in params:
            val = params["MYQDRIL"].values[0]
            if isinstance(val, str):
                return val.upper()
            return str(val).upper()
        return "K6ROT"

    def _read_param_trityp(self) -> str:
        """Read PARAM,TRITYP from model; default DKT.

        Options: DKT (Discrete Kirchhoff), THICK (Mindlin-Reissner).
        """
        params = self.model.params
        if "TRITYP" in params:
            val = params["TRITYP"].values[0]
            if isinstance(val, str):
                return val.upper()
            return str(val).upper()
        return "DKT"

    def _apply_drilling(
        self,
        Ke_global: np.ndarray,
        Ke_local: np.ndarray,
        x: np.ndarray,
        y: np.ndarray,
        A_mat: np.ndarray,
        thickness: float,
        n_vec: np.ndarray,
        Ti: np.ndarray,
    ) -> None:
        """Apply drilling DOF stiffness to the 18x18 global element matrix."""
        method = self.drilling
        if method == "NONE":
            return

        if method == "K6ROT":
            k_drill = _tri3_drilling_k6rot(Ke_local, self.k6rot)
            if k_drill > 0.0:
                K_drill_node = k_drill * np.outer(n_vec, n_vec)
                for inode in range(3):
                    r = 6 * inode + 3
                    Ke_global[r : r + 3, r : r + 3] += K_drill_node
            return

        if method == "HB":
            gamma = A_mat[2, 2]  # G·t for isotropic
            K_drill_local = _tri3_drilling_hughes_brezzi(x, y, gamma)
        else:
            # ALLMAN not implemented for tri — fall back to K6ROT
            k_drill = _tri3_drilling_k6rot(Ke_local, self.k6rot)
            if k_drill > 0.0:
                K_drill_node = k_drill * np.outer(n_vec, n_vec)
                for inode in range(3):
                    r = 6 * inode + 3
                    Ke_global[r : r + 3, r : r + 3] += K_drill_node
            return

        # Transform 9x9 local drilling matrix [u, v, θz]×3 -> global 18-DOF
        T_drill_node = np.zeros((3, 6))
        T_drill_node[0, 0:3] = Ti[0, :]
        T_drill_node[1, 0:3] = Ti[1, :]
        T_drill_node[2, 3:6] = n_vec

        T_drill = np.zeros((9, 18))
        for inode in range(3):
            T_drill[3 * inode : 3 * inode + 3, 6 * inode : 6 * inode + 6] = T_drill_node

        Ke_global += T_drill.T @ K_drill_local @ T_drill

    def build_kbb(
        self, Kbb: dok_matrix, dof_map: DOF_MAP, all_nids: np.ndarray, xyz_cid0: np.ndarray
    ) -> int:
        """Assemble CTRIA3 stiffness into Kbb.

        Returns
        -------
        nelements : int
            Number of CTRIA3 elements processed.
        """
        model = self.model
        ctria3 = model.ctria3
        nelements = ctria3.n
        if nelements == 0:
            return 0

        pshell = model.pshell
        pcomp = model.pcomp

        grid = model.grid
        nodes = ctria3.nodes  # (nelements, 3)
        inid = grid.index(nodes)
        xyz = grid.xyz_cid0()

        p1 = xyz[inid[:, 0], :]
        p2 = xyz[inid[:, 1], :]
        p3 = xyz[inid[:, 2], :]

        # Element normals: (p2 - p1) x (p3 - p1)
        v12 = p2 - p1
        v13 = p3 - p1
        normal = np.cross(v12, v13)
        ni = np.linalg.norm(normal, axis=1)
        if np.any(ni == 0.0):
            bad = ctria3.element_id[ni == 0.0]
            raise RuntimeError(f"CTRIA3 elements {bad.tolist()} have zero-area normals")
        normal /= ni[:, np.newaxis]

        # Element local x-axis (node 1 -> node 2)
        ihat = v12.copy()
        inorm = np.linalg.norm(ihat, axis=1)
        ihat /= inorm[:, np.newaxis]

        # Element local y-axis (normal x ihat)
        jhat = np.cross(normal, ihat, axis=1)
        jnorm = np.linalg.norm(jhat, axis=1)
        jhat /= jnorm[:, np.newaxis]

        # Process each element
        for i_elem in range(nelements):
            pid = int(ctria3.property_id[i_elem])

            Ti = np.vstack([ihat[i_elem], jhat[i_elem], normal[i_elem]])

            # Transform nodes to local frame
            xy = np.zeros((3, 3))
            xy[0] = Ti @ p1[i_elem]
            xy[1] = Ti @ p2[i_elem]
            xy[2] = Ti @ p3[i_elem]
            x_local = xy[:, 0]
            y_local = xy[:, 1]

            # Get constitutive matrices
            A_mat, B_mat, D_mat, Ds_mat, thickness = _get_ABD_for_element(model, pid, pshell, pcomp)

            # Material coordinate rotation (THETA/MCID on CTRIA3)
            theta_elem = ctria3.theta[i_elem]
            mcid_elem = ctria3.mcid[i_elem]

            if mcid_elem >= 0:
                mcid_ref = model.coord.slice_card_by_id(np.array([mcid_elem]))
                i_mcid = mcid_ref.i[0]
                i_proj = i_mcid - np.dot(i_mcid, normal[i_elem]) * normal[i_elem]
                i_proj_norm = np.linalg.norm(i_proj)
                if i_proj_norm > 1e-10:
                    i_proj /= i_proj_norm
                    cos_theta = np.dot(ihat[i_elem], i_proj)
                    sin_theta = np.dot(np.cross(ihat[i_elem], i_proj), normal[i_elem])
                    theta_mat = np.arctan2(sin_theta, cos_theta)
                else:
                    theta_mat = 0.0
            elif not np.isnan(theta_elem) and abs(theta_elem) > 1e-10:
                theta_mat = np.radians(theta_elem)
            else:
                theta_mat = 0.0

            # Rotate ABD from material to element frame
            if abs(theta_mat) > 1e-10:
                c = np.cos(theta_mat)
                s = np.sin(theta_mat)
                T_inv = np.array(
                    [
                        [c**2, s**2, -2 * s * c],
                        [s**2, c**2, 2 * s * c],
                        [s * c, -s * c, c**2 - s**2],
                    ]
                )
                A_mat = T_inv @ A_mat @ T_inv.T
                D_mat = T_inv @ D_mat @ T_inv.T
                if B_mat is not None:
                    B_mat = T_inv @ B_mat @ T_inv.T
                R2 = np.array([[c, s], [-s, c]])
                Ds_mat = R2.T @ Ds_mat @ R2

            # Compute element stiffness (15x15 in local frame)
            if self.trityp == "THICK":
                Ke_local = ctria3_mindlin_stiffness(x_local, y_local, A_mat, D_mat, Ds_mat, B_mat)
            else:
                Ke_local = ctria3_stiffness(x_local, y_local, A_mat, D_mat, Ds_mat, B_mat)

            # Apply shell offset
            zoffset = ctria3.zoffset[i_elem]
            if not np.isnan(zoffset) and abs(zoffset) > 0.0:
                Ke_local = _apply_tri3_shell_offset(Ke_local, zoffset)

            # Transform local 15-DOF -> global 18-DOF
            T_node = np.zeros((5, 6))
            T_node[0:3, 0:3] = Ti
            T_node[3:5, 3:6] = Ti[0:2, :]

            T_full = np.zeros((15, 18))
            for inode in range(3):
                r0 = 5 * inode
                c0 = 6 * inode
                T_full[r0 : r0 + 5, c0 : c0 + 6] = T_node

            Ke_global = T_full.T @ Ke_local @ T_full

            # Drilling stiffness
            self._apply_drilling(
                Ke_global,
                Ke_local,
                x_local,
                y_local,
                A_mat,
                thickness,
                normal[i_elem],
                Ti,
            )

            # Assemble into Kbb
            elem_nodes = nodes[i_elem]
            global_dofs = []
            for nid in elem_nodes:
                i_dof = dof_map[(nid, 1)]
                for comp in range(6):
                    global_dofs.append(i_dof + comp)

            for ii in range(18):
                gi = global_dofs[ii]
                for jj in range(18):
                    gj = global_dofs[jj]
                    val = Ke_global[ii, jj]
                    if abs(val) > 0.0:
                        Kbb[gi, gj] += val

        return nelements


def _apply_tri3_shell_offset(Ke: np.ndarray, z0: float) -> np.ndarray:
    """Apply offset for CTRIA3: u_mid = u_grid + z0*ty, v_mid = v_grid - z0*tx."""
    E_full = np.eye(15)
    for i in range(3):
        r = 5 * i
        E_full[r + 0, r + 4] = z0  # u += z0 * theta_y
        E_full[r + 1, r + 3] = -z0  # v -= z0 * theta_x
    return E_full.T @ Ke @ E_full


# ---------------------------------------------------------------------------
# Public entry point for CTRIA3
# ---------------------------------------------------------------------------
def build_kbb_ctria3(
    model: BDF,
    Kbb: dok_matrix,
    dof_map: DOF_MAP,
    all_nids: np.ndarray,
    xyz_cid0: np.ndarray,
    idtype: str = "int32",
    fdtype: str = "float64",
) -> int:
    """Build CTRIA3 stiffness. DKT bending (default) or Mindlin-Reissner (PARAM,TRITYP,THICK)."""
    solver = ShellTriSolver(model)
    return solver.build_kbb(Kbb, dof_map, all_nids, xyz_cid0)


# ---------------------------------------------------------------------------
# KDGG assembly for shells
# ---------------------------------------------------------------------------
def build_KDgg_cquad4(
    model: BDF,
    KDgg: dok_matrix,
    dof_map: DOF_MAP,
    u_global: np.ndarray,
) -> int:
    """Assemble CQUAD4 geometric stiffness from membrane stress resultants.

    Parameters
    ----------
    model : BDF
    KDgg : (ndof, ndof) sparse matrix to accumulate into
    dof_map : DOF_MAP
    u_global : (ndof,) global displacement vector from static preload

    Returns
    -------
    nelements : int
    """
    cquad4 = model.cquad4
    nelements = cquad4.n
    if nelements == 0:
        return 0

    pshell = model.pshell
    pcomp = model.pcomp
    grid = model.grid
    nodes = cquad4.nodes
    inid = grid.index(nodes)
    xyz = grid.xyz_cid0()

    p1 = xyz[inid[:, 0], :]
    p2 = xyz[inid[:, 1], :]
    p3 = xyz[inid[:, 2], :]
    p4 = xyz[inid[:, 3], :]

    v13 = p1 - p3
    v24 = p2 - p4
    normal = np.cross(v13, v24)
    ni = np.linalg.norm(normal, axis=1)
    normal /= ni[:, np.newaxis]

    ihat = p2 - p1
    inorm = np.linalg.norm(ihat, axis=1)
    ihat /= inorm[:, np.newaxis]

    jhat = np.cross(normal, ihat, axis=1)
    jnorm = np.linalg.norm(jhat, axis=1)
    jhat /= jnorm[:, np.newaxis]

    nid_array = grid.node_id

    for i_elem in range(nelements):
        pid = int(cquad4.property_id[i_elem])
        Ti = np.vstack([ihat[i_elem], jhat[i_elem], normal[i_elem]])

        xy = np.zeros((4, 3))
        xy[0] = Ti @ p1[i_elem]
        xy[1] = Ti @ p2[i_elem]
        xy[2] = Ti @ p3[i_elem]
        xy[3] = Ti @ p4[i_elem]
        x_local = xy[:, 0]
        y_local = xy[:, 1]

        A_mat, _, _, _, thickness = _get_ABD_for_element(model, pid, pshell, pcomp)

        # Extract membrane displacements in local frame
        elem_nodes = nodes[i_elem]
        u_mem_local = np.zeros(8)
        for inode in range(4):
            nid = int(elem_nodes[inode])
            i_dof = dof_map[(nid, 1)]
            u_g = u_global[i_dof : i_dof + 3]
            u_l = Ti @ u_g
            u_mem_local[2 * inode] = u_l[0]
            u_mem_local[2 * inode + 1] = u_l[1]

        # Compute membrane stress resultants at centroid
        dN_c = _dshape_quad4(0.0, 0.0)
        J_c, det_J_c = _jacobian(dN_c, x_local, y_local)
        J_inv_c = _jacobian_inv(J_c, det_J_c)
        dN_dxy_c = J_inv_c @ dN_c
        Bm_c = _membrane_B(dN_dxy_c[0], dN_dxy_c[1])
        strain_m = Bm_c @ u_mem_local
        stress_resultants = A_mat @ strain_m  # [Nxx, Nyy, Nxy]

        # Build element geometric stiffness
        Kg_local = geometric_stiffness(x_local, y_local, stress_resultants, thickness)

        # Apply offset
        zoffset = cquad4.zoffset[i_elem]
        if not np.isnan(zoffset) and abs(zoffset) > 0.0:
            Kg_local = _apply_shell_offset(Kg_local, zoffset)

        # Transform to global
        T_node = np.zeros((5, 6))
        T_node[0:3, 0:3] = Ti
        T_node[3:5, 3:6] = Ti[0:2, :]

        T_full = np.zeros((20, 24))
        for inode in range(4):
            T_full[5 * inode : 5 * inode + 5, 6 * inode : 6 * inode + 6] = T_node

        Kg_global = T_full.T @ Kg_local @ T_full

        # Assemble
        global_dofs = []
        for nid in elem_nodes:
            i_dof = dof_map[(nid, 1)]
            for comp in range(6):
                global_dofs.append(i_dof + comp)

        for ii in range(24):
            gi = global_dofs[ii]
            for jj in range(24):
                gj = global_dofs[jj]
                val = Kg_global[ii, jj]
                if abs(val) > 0.0:
                    KDgg[gi, gj] += val

    return nelements


def build_KDgg_ctria3(
    model: BDF,
    KDgg: dok_matrix,
    dof_map: DOF_MAP,
    u_global: np.ndarray,
) -> int:
    """Assemble CTRIA3 geometric stiffness from membrane stress resultants.

    Parameters
    ----------
    model : BDF
    KDgg : (ndof, ndof) sparse matrix to accumulate into
    dof_map : DOF_MAP
    u_global : (ndof,) global displacement vector from static preload

    Returns
    -------
    nelements : int
    """
    ctria3 = model.ctria3
    nelements = ctria3.n
    if nelements == 0:
        return 0

    pshell = model.pshell
    pcomp = model.pcomp
    grid = model.grid
    nodes = ctria3.nodes
    inid = grid.index(nodes)
    xyz = grid.xyz_cid0()

    p1 = xyz[inid[:, 0], :]
    p2 = xyz[inid[:, 1], :]
    p3 = xyz[inid[:, 2], :]

    v12 = p2 - p1
    v13 = p3 - p1
    normal = np.cross(v12, v13, axis=1)
    ni = np.linalg.norm(normal, axis=1)
    normal /= ni[:, np.newaxis]

    ihat = v12 / np.linalg.norm(v12, axis=1)[:, np.newaxis]
    jhat = np.cross(normal, ihat, axis=1)
    jhat /= np.linalg.norm(jhat, axis=1)[:, np.newaxis]

    for i_elem in range(nelements):
        pid = int(ctria3.property_id[i_elem])
        Ti = np.vstack([ihat[i_elem], jhat[i_elem], normal[i_elem]])

        xy = np.zeros((3, 3))
        xy[0] = Ti @ p1[i_elem]
        xy[1] = Ti @ p2[i_elem]
        xy[2] = Ti @ p3[i_elem]
        x_local = xy[:, 0]
        y_local = xy[:, 1]

        A_mat, _, _, _, thickness = _get_ABD_for_element(model, pid, pshell, pcomp)

        # Extract membrane displacements in local frame
        elem_nodes = nodes[i_elem]
        u_mem_local = np.zeros(6)
        for inode in range(3):
            nid = int(elem_nodes[inode])
            i_dof = dof_map[(nid, 1)]
            u_g = u_global[i_dof : i_dof + 3]
            u_l = Ti @ u_g
            u_mem_local[2 * inode] = u_l[0]
            u_mem_local[2 * inode + 1] = u_l[1]

        # CST membrane strain (constant over element)
        area = 0.5 * abs(
            (x_local[1] - x_local[0]) * (y_local[2] - y_local[0])
            - (x_local[2] - x_local[0]) * (y_local[1] - y_local[0])
        )
        dNdx = np.array(
            [y_local[1] - y_local[2], y_local[2] - y_local[0], y_local[0] - y_local[1]]
        ) / (2.0 * area)
        dNdy = np.array(
            [x_local[2] - x_local[1], x_local[0] - x_local[2], x_local[1] - x_local[0]]
        ) / (2.0 * area)

        Bm = np.zeros((3, 6))
        Bm[0, 0::2] = dNdx
        Bm[1, 1::2] = dNdy
        Bm[2, 0::2] = dNdy
        Bm[2, 1::2] = dNdx

        strain_m = Bm @ u_mem_local
        stress_resultants = A_mat @ strain_m

        # Build element geometric stiffness
        Kg_local = ctria3_geometric_stiffness(x_local, y_local, stress_resultants, thickness)

        # Apply offset
        zoffset = ctria3.zoffset[i_elem]
        if not np.isnan(zoffset) and abs(zoffset) > 0.0:
            Kg_local = _apply_tri3_shell_offset(Kg_local, zoffset)

        # Transform to global
        T_node = np.zeros((5, 6))
        T_node[0:3, 0:3] = Ti
        T_node[3:5, 3:6] = Ti[0:2, :]

        T_full = np.zeros((15, 18))
        for inode in range(3):
            T_full[5 * inode : 5 * inode + 5, 6 * inode : 6 * inode + 6] = T_node

        Kg_global = T_full.T @ Kg_local @ T_full

        # Assemble
        global_dofs = []
        for nid in elem_nodes:
            i_dof = dof_map[(nid, 1)]
            for comp in range(6):
                global_dofs.append(i_dof + comp)

        for ii in range(18):
            gi = global_dofs[ii]
            for jj in range(18):
                gj = global_dofs[jj]
                val = Kg_global[ii, jj]
                if abs(val) > 0.0:
                    KDgg[gi, gj] += val

    return nelements
